using MAT
# using Plots
using Dates
using Random
using Dierckx
using Cambrian
using SeisNoise
using Statistics
using NumericalIntegration
using CartesianGeneticProgramming
include("model_template_utils.jl")
# include("shorefor_functions.jl")
include("create_shorefor_ind.jl")
include("evaluation_metrics.jl")
include("dataloader_monthly.jl")


cfg_filename = "cfg/shorefor_dxdt_MO-1.yaml"
cfg = get_config(cfg_filename; columns=50, n_in=9, n_out=1, n_population=200, d_fitness=5)
# cfg = get_config(ARGS[1])

ind_shorefor = get_fullmodel_shorefor_dxdt_v2(cfg)
pure_shorefor_output_trace =  get_output_trace(ind_shorefor)

function intersect_lists(shoreline_dnums, wave_dnums)
    shoreline_dnums = vec(shoreline_dnums)
    wave_dnums = vec(wave_dnums)
    i_shorelines = findall(x -> x in wave_dnums, shoreline_dnums)
    i_waves = findall(x -> x in shoreline_dnums, wave_dnums)
    return i_shorelines, i_waves
end

function inputs_from_dataset(dataset)
    Tp, Hsb, Dir, dates_waves, loc_shore, t_shore, E, Sla, rivdis, omega, P,  t, phi = dataset
    dt_waves = dates_waves[2] - dates_waves[1]
    phi = ceil(phi / dt_waves)
    
    inputs = [omega, P, phi, 2phi, Dir, Hsb, Tp, Sla, rivdis]
    return inputs
end

function make_shorefor_individual(cfg::NamedTuple)
    get_fullmodel_shorefor_dxdt_v2(cfg)
end

function evaluate_dataset(ind::CGPInd, dataset)
    CartesianGeneticProgramming.reset!(ind)
    Tp, Hsb, Dir, dates_waves, loc_shore, t_shore, E, Sla, rivdis, omega, P,  t, phi = dataset
    dt_waves = dates_waves[2] - dates_waves[1]
    phi = ceil(phi / dt_waves)
    local inputs
    
    inputs = [omega, P, phi, 2phi, Dir, Hsb, Tp, Sla, rivdis]

    res = process(ind, inputs)
    dxdts = res[1]

    if dxdts == nothing || length(dxdts) != length(omega)
        return -Inf
    end

    i_pred, i_target = intersect_lists(dates_waves, t_shore)

    vx = collect(1.0:length(dxdts))
    X = NumericalIntegration.cumul_integrate(vx, dxdts) .* dt_waves
    B = hcat(ones(length(dxdts)), t, X)

    Â = B[i_pred, :]\loc_shore[i_target]
    xs = B * Â
    
    eval_pred = xs[i_pred]
    eval_target = loc_shore[i_target]
    score = mielke(eval_target, eval_pred)

    if isequal(score, NaN)
        return -Inf
    else
        return score
    end
end

dataset_grandpopo = load_monthly_dataset_GRANDPOPO(;calibration=true, remove_means=false, normalize=true)
dataset_narrabeen = load_monthly_dataset_NARRABEEN(;calibration=true, remove_means=false, normalize=true)
dataset_duck = load_monthly_dataset_DUCK(;calibration=true, remove_means=false, normalize=true)
dataset_torreypines = load_monthly_dataset_TORREYPINES(;calibration=true, remove_means=false, normalize=true)
dataset_trucvert = load_monthly_dataset_TRUCVERT(;calibration=true, remove_means=false, normalize=true)

function evaluate_individual(ind)
    score_grandpop = evaluate_dataset(ind, dataset_grandpopo)
    score_narrabeen = evaluate_dataset(ind, dataset_narrabeen)
    score_duck = evaluate_dataset(ind, dataset_duck)
    score_torreypines = evaluate_dataset(ind, dataset_torreypines)
    score_trucvert = evaluate_dataset(ind, dataset_trucvert)
    scores = [score_grandpop, score_narrabeen, score_duck, score_torreypines, score_trucvert]
    if -Inf in scores
        scores .= -Inf
    end
    return scores 
end

init_seed = cfg.seed
n_runs = 5 # tryparse(Int, ARGS[2])
start_seed = (init_seed-1) * n_runs
end_seed = start_seed + n_runs - 1
println("seeds: $(collect(start_seed:end_seed))")

for s in collect(start_seed:end_seed)
    Random.seed!(s)
    
    job_name = "cgpshorefor-nsgaii-5obj-run$s"
    cfg = get_config(cfg_filename; n_population=200, id=job_name, n_in=9, d_fitness=5, n_offsprings=200, n_gen=10)  # ARGS[1]
    
    fit(i::CGPInd) = evaluate_individual(i)
    test_inputs = inputs_from_dataset(dataset_grandpopo)
    CartesianGeneticProgramming.mutate(i::CGPInd) = goldman_mutate(cfg, i, test_inputs, evaluate_individual)
    logfile = joinpath(cfg.output_dir, "logs", string(cfg.id, ".csv"))
    evo_dir = joinpath(cfg.output_dir, "gens", job_name)
    if isdir(evo_dir)
        gen_nums = tryparse.(Int64, readdir(evo_dir))
        max_gen = maximum(gen_nums)
        if max_gen < cfg.n_gen
            e = NSGA2Evolution(cfg, fit; logfile=logfile, population_dir="$evo_dir/$max_gen")
            run!(e)
        end
    else
        e = NSGA2Evolution(cfg, fit; logfile=logfile, init_function=make_shorefor_individual)
        run!(e)
    end
end