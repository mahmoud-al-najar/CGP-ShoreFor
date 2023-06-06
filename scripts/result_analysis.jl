using MAT
using Plots
using Plots.PlotMeasures
using Dates
using Random
using Dierckx
using Cambrian
using SeisNoise
using NaNStatistics
using NumericalIntegration
using CartesianGeneticProgramming
include("model_template_utils.jl")
include("create_shorefor_ind.jl")
include("evaluation_metrics.jl")
include("dataloader_monthly.jl")
include("../graphing/graphs.jl")

plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing, grid=true, guidefontsize=12, labelfontsize=12, legendfontsize=12, tickfontsize=12, dpi=200, foreground_color_legend=nothing, background_color_legend=nothing)#, yguidefontrotation=-90)

cfg = get_config("cfg/shorefor_dxdt_MO-1.yaml"; columns=50, n_in=9, n_out=1, n_population=200, d_fitness=5)

ind_shorefor = get_fullmodel_shorefor_dxdt_v2(cfg)

function intersect_lists(shoreline_dnums, wave_dnums)
    shoreline_dnums = vec(shoreline_dnums)
    wave_dnums = vec(wave_dnums)
    i_shorelines = findall(x -> x in wave_dnums, shoreline_dnums)
    i_waves = findall(x -> x in shoreline_dnums, wave_dnums)
    return i_shorelines, i_waves
end

function evaluate_dataset_calib_forecast(ind::CGPInd, dataset_calib, dataset_forecast; eval_start_index=1, eval_end_index=nothing)
    CartesianGeneticProgramming.reset!(ind)
    Tp, Hsb, Dir, dates_waves, loc_shore, t_shore, E, Sla, rivdis, omega, P,  t, phi = dataset_calib
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
    if eval_end_index === nothing
        eval_pred = xs[i_pred][eval_start_index:end]
        eval_target = loc_shore[i_target][eval_start_index:end]
    else
        eval_pred = xs[i_pred][eval_start_index:eval_end_index]
        eval_target = loc_shore[i_target][eval_start_index:eval_end_index]
    end
    
    calib_score = mielke(eval_target, eval_pred)

    CartesianGeneticProgramming.reset!(ind)
    Tp, Hsb, Dir, dates_waves, loc_shore, t_shore, E, Sla, rivdis, omega, P,  t, phi = dataset_forecast
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

    xs = B * Â
    eval_start_index = length(dataset_calib[1])
    if eval_end_index === nothing
        eval_pred = xs[i_pred][eval_start_index:end]
        eval_target = loc_shore[i_target][eval_start_index:end]
    else
        eval_pred = xs[i_pred][eval_start_index:eval_end_index]
        eval_target = loc_shore[i_target][eval_start_index:eval_end_index]
    end
    
    forec_score = mielke(eval_target, eval_pred)
    
    if isequal(calib_score, NaN) || isequal(forec_score, NaN)
        return -Inf, -Inf
    else
        return calib_score, forec_score
    end
end

sites = ["Grand Popo", "Narrabeen", "Duck", "Torrey Pines", "Truc Vert"]

dataset_grandpopo_calib = load_monthly_dataset_GRANDPOPO(;calibration=true, remove_means=false, normalize=true)
dataset_narrabeen_calib = load_monthly_dataset_NARRABEEN(;calibration=true, remove_means=false, normalize=true)
dataset_duck_calib = load_monthly_dataset_DUCK(;calibration=true, remove_means=false, normalize=true)
dataset_torreypines_calib = load_monthly_dataset_TORREYPINES(;calibration=true, remove_means=false, normalize=true)
dataset_trucvert_calib = load_monthly_dataset_TRUCVERT(;calibration=true, remove_means=false, normalize=true)
datasets_calib = [
    dataset_grandpopo_calib, 
    dataset_narrabeen_calib, 
    dataset_duck_calib, 
    dataset_torreypines_calib, 
    dataset_trucvert_calib
]

dataset_grandpopo_full = load_monthly_dataset_GRANDPOPO(;calibration=false, remove_means=false, normalize=true)
dataset_narrabeen_full = load_monthly_dataset_NARRABEEN(;calibration=false, remove_means=false, normalize=true)
dataset_duck_full = load_monthly_dataset_DUCK(;calibration=false, remove_means=false, normalize=true)
dataset_torreypines_full = load_monthly_dataset_TORREYPINES(;calibration=false, remove_means=false, normalize=true)
dataset_trucvert_full = load_monthly_dataset_TRUCVERT(;calibration=false, remove_means=false, normalize=true)
datasets_full = [
    dataset_grandpopo_full, 
    dataset_narrabeen_full, 
    dataset_duck_full, 
    dataset_torreypines_full, 
    dataset_trucvert_full
]

function evaluate_calibration_forecast_all_individuals(cfg, inds)
    cal_scores = Array{Float64}(undef, length(inds), cfg.d_fitness)
    all_scores = Array{Float64}(undef, length(inds), cfg.d_fitness)
    for i in eachindex(inds)
        ind = inds[i]
        if i % 2000 == 0
            println("i: $i")
        end
        try
            score_grandpop_calib, score_grandpop_forec = evaluate_dataset_calib_forecast(ind, dataset_grandpopo_calib, dataset_grandpopo_full)
            score_narrabeen_calib, score_narrabeen_forec = evaluate_dataset_calib_forecast(ind, dataset_narrabeen_calib, dataset_narrabeen_full)
            score_duck_calib, score_duck_forec = evaluate_dataset_calib_forecast(ind, dataset_duck_calib, dataset_duck_full)
            score_torreypines_calib, score_torreypines_forec = evaluate_dataset_calib_forecast(ind, dataset_torreypines_calib, dataset_torreypines_full)
            score_trucvert_calib, score_trucvert_forec = evaluate_dataset_calib_forecast(ind, dataset_trucvert_calib, dataset_trucvert_full)
            
            cal_scores[i, :] .= [score_grandpop_calib, score_narrabeen_calib, score_duck_calib, score_torreypines_calib, score_trucvert_calib]
            all_scores[i, :] .= [score_grandpop_forec, score_narrabeen_forec, score_duck_forec, score_torreypines_forec, score_trucvert_forec]
        catch e
            println(i, "  ERROR::  ", e)
            score_grandpop_calib = score_narrabeen_calib = score_duck_calib = score_torreypines_calib = score_trucvert_calib = -Inf    
            score_grandpop_forec = score_narrabeen_forec = score_duck_forec = score_torreypines_forec = score_trucvert_forec = -Inf    
            cal_scores[i, :] .= [score_grandpop_calib, score_narrabeen_calib, score_duck_calib, score_torreypines_calib, score_trucvert_calib]
            all_scores[i, :] .= [score_grandpop_forec, score_narrabeen_forec, score_duck_forec, score_torreypines_forec, score_trucvert_forec]
        end
        # if ind.fitness != all_scores[i, :]
        #     println("NOT EQUAL =========== $i")
        # end
        ind.fitness .= all_scores[i, :]
    end
    return cal_scores, all_scores
end

function sortrows(M, by=zeros(0))
    # from: https://discourse.julialang.org/t/sort-matrix-based-on-the-elements-of-a-specific-column/23475/5
    if by == zeros(0)
        order = copy(M)
    else
        order = copy(by)
    end
    if size(order,2) > 1
        order = Float64.(order.-minimum(order, dims = 1))
        order = (order./maximum(order,dims=1))*(10).^(size(order,2):-1:1)
    end
    order = sortperm(order[:,1])
    return M[order,:], order
end


lowF = T -> 1.5 * T
highF = T -> 0.75 * T
passband = (S,lowf,highf) -> movmean(S - movmean(S, lowf), highf)

function _spectral_eval(x, x̂)
    window_lengths = collect(1:length(x))
    ma_x = Array{Float64}(undef, length(window_lengths), length(x))
    ma_x̂ = Array{Float64}(undef, length(window_lengths), length(x))
    mielke_scores = Array{Float64}(undef, length(window_lengths))
    for l in window_lengths
        ma_x[l, :] = passband(x, lowF(l), highF(l)) ./ l
        ma_x̂[l, :] = passband(x̂, lowF(l), highF(l)) ./ l
        mielke_scores[l] = mielke(ma_x[l, :], ma_x̂[l, :])
    end
    return ma_x, ma_x̂, window_lengths, mielke_scores
end

function singel_eval_plot(t, t_calib, ts, X_shorefor, X_inds, sitename; left_margin=12mm, right_margin=4mm, figsize=(1000, 180), bottom_margin=4mm, legend=false, linewidth=1, ind_labels=["CGP"])
    colors=palette(:rainbow, size(X_inds)[1] + 1)
    x_vals = num2date.(t)

    x_m_ticks = map(m -> m[1:3], Dates.monthname.(num2date.(t)))
    x_y_ticks = map(m -> "’" * string(m)[end-1:end], Dates.year.(num2date.(t)))
    i_start_month = findfirst(x -> x in ["Jan"], x_m_ticks)
    i_last_month = findlast(x -> x in ["Jan"], x_m_ticks)
    dm = 12

    if length(x_vals[i_start_month:dm:i_last_month]) > 8
        dm = 24
    end
        
    xticks = (x_vals[i_start_month:dm:i_last_month], x_y_ticks[i_start_month:dm:i_last_month])

    p_temp = Plots.plot(num2date.(t), ts, title=LaTeXString(sitename), ylabel="Shoreline", label="Target", color=:black, linewidth=linewidth, linealpha=0.9)
    
    Plots.plot!(p_temp, num2date.(t), X_shorefor, ylabel="Shoreline", label="ShoreFor", color=:red, linewidth=linewidth, linealpha=0.9)
    for i in 1:size(X_inds)[1]
        Plots.plot!(p_temp, num2date.(t), X_inds[i, :], ylabel="Shoreline", label=ind_labels[i], color=colors[i], linewidth=linewidth, linealpha=0.9)
    end
    yticks = (collect(0:0.25:1), rpad.(string.(collect(0:0.25:1)), 4, "0"))

    Plots.vspan!(p_temp, [num2date.(t[length(t_calib)]), num2date.(t[end])], linecolor = :grey, fillcolor = :grey, alpha=0.5, yticks=false)
    # Plots.plot!(p_temp, [num2date.(t[length(t_calib)]), num2date.(t[end])], yticks=false, seriestype=:vspan, alpha=0.5, linecolor = :grey, fillcolor = :grey)
    # Plots.vspan!(p_temp, [num2date.(t[length(t_calib)]), num2date.(t[end])], [0.0, 1.0], xticks=false, seriestype=:vspan, alpha=0.5, linecolor = :grey, fillcolor = :grey)
    
    return Plots.plot(p_temp, size=figsize, left_margin=left_margin, right_margin=right_margin, bottom_margin=bottom_margin, xticks=xticks, xminorticks=4, legend=legend, legend_columns=-1, yticks=yticks)
end

function evaluate_multiple_runs_dataset_v4_single_plots(inds, datasets_calib, datasets_fullseries, site_names; pdf_path="", timeseries_name="", plot_size=(1080, 720), linewidth=1, ind_labels)
    # colors=palette(:rainbow, length(inds) + 1)
    colors=palette(:rainbow, length(inds) + 1)

    all_inds_cal_scores = Array{Float64, 2}(undef, length(inds), length(site_names))
    all_inds_for_scores = Array{Float64, 2}(undef, length(inds), length(site_names))

    for di in eachindex(site_names)
        ts_X = datasets_fullseries[di][5]
        ts_X_shorefor = Array{Float64, 1}(undef, length(datasets_fullseries[di][5]))
        ts_X_inds = Array{Float64, 2}(undef, length(inds), length(datasets_fullseries[di][5]))

        site_name = site_names[di]
        # inds = inds_per_site[di]
        Tp, Hsb, Dir, dates_waves, loc_shore, t_shore, E, Sla, rivdis, omega, P,  t, phi = datasets_calib[di]
        calibration_shores_length = length(loc_shore)
        dt_waves = dates_waves[2] - dates_waves[1]
        phi = ceil(phi / dt_waves)
        local inputs = [omega, P, phi, 2phi, Dir, Hsb, Tp, Sla, rivdis]
        # inputs = [omega, P, phi, 2phi, Dir, Hsb, Tp]

        x_vals = num2date.(t_shore)
        x_ym = Dates.format.(num2date.(t_shore), "yyyy-mm")
        x_m_ticks = map(m -> m[1:3], Dates.monthname.(num2date.(t_shore)))
        x_y_ticks = map(m -> "’" * string(m)[end-1:end], Dates.year.(num2date.(t_shore)))
        i_start_month = findfirst(x -> x in ["Jan"], x_m_ticks)
        i_last_month = findlast(x -> x in ["Jan"], x_m_ticks)
        dm = 12
        if length(x_vals[i_start_month:dm:i_last_month]) > 5
            dm = 24
        end
            
        xticks = (x_vals[i_start_month:dm:i_last_month], x_y_ticks[i_start_month:dm:i_last_month])
        
        p1 = Plots.plot(x_vals, loc_shore, xticks=xticks, label="Target", color=:black, legend=:outerbottomright, linewidth=linewidth*2, linealpha=0.7, ylabel=L"X", xlabel=LaTeXString(L"$Years$"), title=site_name, xminorticks=4, bottom_margin=8mm)
        
        p_spec = Plots.plot()

        ### ShoreFor results
        CartesianGeneticProgramming.reset!(ind_shorefor)
        res = process(ind_shorefor, inputs)
        dxdts = res[1]
        i_pred, i_target = intersect_lists(dates_waves, t_shore)
        vx = collect(1.0:length(dxdts))
        X = NumericalIntegration.cumul_integrate(vx, dxdts) .* dt_waves
        B = hcat(ones(length(dxdts)), t, X)
        Â = B[i_pred, :]\loc_shore[i_target]
        xs = B * Â
        eval_pred = xs[i_pred]
        eval_target = loc_shore[i_target]
        score = mielke(eval_target, eval_pred; digits=2)
        ma_x, ma_x̂, window_lengths, shorefor_mielke_scores = _spectral_eval(eval_target, eval_pred)
        # Plots.plot!(p_spec, window_lengths, shorefor_mielke_scores, label="ShoreFor", linecolor=:red, ylabel=L"Mielke", xlabel=LaTeXString(L"$Timescale [months]$"), legend = false, linewidth=linewidth, linealpha=0.9, bottom_margin=8mm)
        ### ECJ:
        # Plots.plot!(p1, x_vals, xs[i_pred], xticks=xticks, label="ShoreFor, λ=$score", linecolor=:red, linewidth=linewidth, linealpha=0.9)
        Plots.plot!(p1, x_vals, xs[i_pred], xticks=xticks, label="ShoreFor", linecolor=:red, linewidth=linewidth, linealpha=0.9)
        ### ICLR:
        # Plots.plot!(p1, x_vals, xs[i_pred], xticks=xticks, label="ShoreFor, λ=$score", linecolor=:red, linewidth=linewidth*2, linealpha=0.9)
        ts_X_shorefor[1:calibration_shores_length] .= xs[i_pred]
        ###

        ### p_spec zero line
        zero_line = zeros(length(shorefor_mielke_scores))
        Plots.plot!(p_spec, window_lengths, zero_line, label=false, linecolor=:grey, linestyle=:dot, linewidth=linewidth*.7, linealpha=0.8)
        ###

        scores = Array{Float64}(undef, length(inds))
        Â_s = []
        for i in eachindex(inds)
            ind = inds[i]
            CartesianGeneticProgramming.reset!(ind)
            
            res = process(ind, inputs)
            dxdts = res[1]

            if dxdts == nothing || length(dxdts) != length(omega)
                println("dxdts == nothing || length(dxdts) != length(omega)  ::  $(objectid(ind))")
            else
                i_pred, i_target = intersect_lists(dates_waves, t_shore)

                vx = collect(1.0:length(dxdts))
                X = NumericalIntegration.cumul_integrate(vx, dxdts) .* dt_waves
                B = hcat(ones(length(dxdts)), t, X)

                Â = B[i_pred, :]\loc_shore[i_target]
                push!(Â_s, Â)
                xs = B * Â
                
                eval_pred = xs[i_pred]
                eval_target = loc_shore[i_target]
                score = mielke(eval_target, eval_pred; digits=2)
                # score = corr(eval_target, eval_pred; digits=2) ^2
                scores[i] = score

                
                ma_x, ma_x̂, window_lengths, mielke_scores = _spectral_eval(eval_target, eval_pred)
                ### ECJ::
                Plots.plot!(p_spec, window_lengths, mielke_scores .- shorefor_mielke_scores, label="$(ind_labels[i])", linecolor=colors[i], ylabel="", xlabel="", legend = false, linewidth=linewidth, linealpha=0.9, bottom_margin=8mm)
                # Plots.plot!(p1, x_vals, xs[i_pred], xticks=xticks, label="$(ind_labels[i]), λ=$score", linecolor=colors[i], linewidth=linewidth, linealpha=0.9)
                Plots.plot!(p1, x_vals, xs[i_pred], xticks=xticks, label="$(ind_labels[i])", linecolor=colors[i], linewidth=linewidth, linealpha=0.9)
                ### ICLR::
                ## relative
                # Plots.plot!(p_spec, window_lengths, mielke_scores .- shorefor_mielke_scores, label="$(ind_labels[i])", linecolor=:blue, ylabel="", xlabel="", legend = false, linewidth=linewidth*2, linealpha=0.9, bottom_margin=8mm)
                ## raw
                # Plots.plot!(p_spec, window_lengths, shorefor_mielke_scores, label="ShoreFor", linecolor=:red, ylabel="λ", xlabel="Months", legend = false, linewidth=linewidth*2, linealpha=0.9, bottom_margin=8mm)
                # Plots.plot!(p_spec, window_lengths, mielke_scores, label="$(ind_labels[i])", linecolor=:blue, ylabel="λ", xlabel="Months", legend = false, linewidth=linewidth*2, linealpha=0.9, bottom_margin=8mm)
                # Plots.plot!(p1, x_vals, xs[i_pred], xticks=xticks, label="$(ind_labels[i]), λ=$score", linecolor=:blue, linewidth=linewidth*2, linealpha=0.9)
                ts_X_inds[i, 1:calibration_shores_length] .= xs[i_pred]
                # Plots.plot!(p2, x_vals, eval_pred, xticks=xticks, label="ind$i, $score", linecolor=colors[i], linewidth=linewidth, linealpha=0.9)
            end
        end
        # push!(cal_plots, p1)
        # push!(cal_spec_plots, p_spec)
        score_text = "Model, Mielke\tCALIBRATION\n"
        for i in eachindex(scores)
            score_text *= "$(ind_labels[i]), $(scores[i])\n"
        end
        println(score_text)
        all_inds_cal_scores[:, di] .= scores

        # annotate!(p2, text(score_text, :red, :right, 3))
        plot!(p1, size=(680, 180), dpi=200)
        plot!(p_spec, size=(380, 180), dpi=200)
        p3 = Plots.plot(p1, p_spec, layout=Plots. grid(1,2, widths=(0.72, 0.28)), size=(1000,200), left_margin=6mm)
        savefig(p3, "$pdf_path/$timeseries_name-calibration-d$di.pdf")

        # ####################
        site_name = site_names[di]
        # inds = inds_per_site[di]
        Tp, Hsb, Dir, dates_waves, loc_shore, t_shore, E, Sla, rivdis, omega, P,  t, phi = datasets_fullseries[di]
        dt_waves = dates_waves[2] - dates_waves[1]
        phi = ceil(phi / dt_waves)
        local inputs = [omega, P, phi, 2phi, Dir, Hsb, Tp, Sla, rivdis]

        x_vals = num2date.(t_shore[calibration_shores_length:end])
        x_ym = Dates.format.(num2date.(t_shore[calibration_shores_length:end]), "yyyy-mm")
        x_m_ticks = map(m -> m[1:3], Dates.monthname.(num2date.(t_shore[calibration_shores_length:end])))
        x_y_ticks = map(m -> "’" * string(m)[end-1:end], Dates.year.(num2date.(t_shore[calibration_shores_length:end])))
        i_start_month = findfirst(x -> x in ["Jan"], x_m_ticks)
        i_last_month = findlast(x -> x in ["Jan"], x_m_ticks)
        dm = 12
            
        xticks = (x_vals[i_start_month:dm:i_last_month], x_y_ticks[i_start_month:dm:i_last_month])

        p2 = Plots.plot(x_vals, loc_shore[calibration_shores_length:end], xticks=xticks, label="Target", color=:black, legend=:outerbottomright, linewidth=linewidth*2, linealpha=0.7, ylabel=L"X", xlabel=LaTeXString(L"$Years$"), title=site_name, xminorticks=4, bottom_margin=8mm)
        p_spec = Plots.plot()

        ### ShoreFor results
        CartesianGeneticProgramming.reset!(ind_shorefor)
        res = process(ind_shorefor, inputs)
        dxdts = res[1]
        i_pred, i_target = intersect_lists(dates_waves, t_shore)
        vx = collect(1.0:length(dxdts))
        X = NumericalIntegration.cumul_integrate(vx, dxdts) .* dt_waves
        B = hcat(ones(length(dxdts)), t, X)
        Â = B[i_pred, :]\loc_shore[i_target]
        xs = B * Â
        eval_pred = xs[i_pred][calibration_shores_length:end]
        eval_target = loc_shore[i_target][calibration_shores_length:end]
        score = mielke(eval_target, eval_pred; digits=2)
        ma_x, ma_x̂, window_lengths, shorefor_mielke_scores = _spectral_eval(eval_target, eval_pred)        
        # Plots.plot!(p_spec, window_lengths, shorefor_mielke_scores, label="ShoreFor", linecolor=:red, ylabel=L"Mielke", xlabel=LaTeXString(L"$Timescale [months]$"), legend = false, linewidth=linewidth, linealpha=0.9, bottom_margin=8mm)
        ### ECJ::
        # Plots.plot!(p2, x_vals, eval_pred, xticks=xticks, label="ShoreFor, λ=$score", linecolor=:red, linewidth=linewidth, linealpha=0.9)
        Plots.plot!(p2, x_vals, eval_pred, xticks=xticks, label="ShoreFor", linecolor=:red, linewidth=linewidth, linealpha=0.9)
        ### ICLR::
        # Plots.plot!(p2, x_vals, eval_pred, xticks=xticks, label="ShoreFor, λ=$score", linecolor=:red, linewidth=linewidth*2, linealpha=0.9)
        ts_X_shorefor[calibration_shores_length:end] .= eval_pred
        ###
        
        ### p_spec zero line
        zero_line = zeros(length(shorefor_mielke_scores))
        Plots.plot!(p_spec, window_lengths, zero_line, label=false, linecolor=:grey, linestyle=:dot, linewidth=linewidth*.7, linealpha=0.8)
        ###

        scores = Array{Float64}(undef, length(inds))
        for i in eachindex(inds)
            
            ind = inds[i]
            CartesianGeneticProgramming.reset!(ind)
            
            res = process(ind, inputs)
            dxdts = res[1]
            
            if dxdts == nothing || length(dxdts) != length(omega)
                println("dxdts == nothing || length(dxdts) != length(omega)  ::  $(objectid(ind))")
            else
                i_pred, i_target = intersect_lists(dates_waves, t_shore)

                vx = collect(1.0:length(dxdts))
                X = NumericalIntegration.cumul_integrate(vx, dxdts) .* dt_waves
                B = hcat(ones(length(dxdts)), t, X)

                Â = Â_s[i] # B[i_pred, :]\loc_shore[i_target]
                # println("old: $(Â_s[i]) -- new: $Â")
                xs = B * Â
                
                eval_pred = xs[i_pred][calibration_shores_length:end]
                eval_target = loc_shore[i_target][calibration_shores_length:end]
                score = mielke(eval_target, eval_pred; digits=2)
                scores[i] = score
                
                ma_x, ma_x̂, window_lengths, mielke_scores = _spectral_eval(eval_target, eval_pred)
                ### ECJ::
                Plots.plot!(p_spec, window_lengths, mielke_scores .- shorefor_mielke_scores, label="$(ind_labels[i])", linecolor=colors[i], ylabel="", xlabel="", legend = false, linewidth=linewidth, linealpha=0.9, bottom_margin=8mm)
                # Plots.plot!(p2, x_vals, eval_pred, xticks=xticks, label="$(ind_labels[i]), λ=$score", linecolor=colors[i], linewidth=linewidth, linealpha=0.9)
                Plots.plot!(p2, x_vals, eval_pred, xticks=xticks, label="$(ind_labels[i])", linecolor=colors[i], linewidth=linewidth, linealpha=0.9)
                ### ICLR::
                ## relative
                # Plots.plot!(p_spec, window_lengths, mielke_scores .- shorefor_mielke_scores, label="$(ind_labels[i])", linecolor=:blue, ylabel="", xlabel="", legend = false, linewidth=linewidth*2, linealpha=0.9, bottom_margin=8mm)
                ## raw
                # Plots.plot!(p_spec, window_lengths, shorefor_mielke_scores, label="ShoreFor", linecolor=:red, ylabel="λ", xlabel="Months", legend = false, linewidth=linewidth*2, linealpha=0.9, bottom_margin=8mm)
                # Plots.plot!(p_spec, window_lengths, mielke_scores, label="$(ind_labels[i])", linecolor=:blue, ylabel="λ", xlabel="Months", legend = false, linewidth=linewidth*2, linealpha=0.9, bottom_margin=8mm)
                # Plots.plot!(p2, x_vals, eval_pred, xticks=xticks, label="$(ind_labels[i]), λ=$score", linecolor=:blue, linewidth=linewidth*2, linealpha=0.9)
                ts_X_inds[i, calibration_shores_length:end] .= eval_pred
            end
            
        end
        
        score_text = "Model, Mielke\tFORECAST-d$di\n"
        for i in eachindex(scores)
            score_text *= "$(ind_labels[i]), $(scores[i])\n"
        end
        println(score_text)
        all_inds_for_scores[:, di] .= scores

        # annotate!(p2, text(score_text, :red, :right, 3))
        plot!(p2, size=(680, 180), dpi=200)
        plot!(p_spec, size=(380, 180), dpi=200)
        p4 = Plots.plot(p2,p_spec, layout=Plots.grid(1,2, widths=(0.72, 0.28)), size=(1000,200), left_margin=6mm)
        savefig(p4, "$pdf_path/$timeseries_name-forecast-d$di.pdf")
        p5 = Plots.plot(p3,p4, layout=Plots.grid(2,1), size=(1000,200*2))
        savefig(p5, "$pdf_path/$timeseries_name-BOTH-d$di.pdf")

        if di == 5
            legend = :outerbottom
            figsize = (1400, 180)
            
            sp = singel_eval_plot(datasets_fullseries[di][4], datasets_calib[di][4], datasets_fullseries[di][5], ts_X_shorefor, ts_X_inds, site_name, legend=legend, figsize=figsize, ind_labels=ind_labels, linewidth=3)
            savefig(sp, "$pdf_path/$timeseries_name-SP-d$di-CBAR.pdf")

            legend = false
            figsize = (1000, 140)
        else
            legend = false
            figsize = (1000, 140)
        end
        # sp = singel_eval_plot(datasets_fullseries[di][4], datasets_calib[di][4], datasets_fullseries[di][5], ts_X_shorefor, ts_X_inds, site_name, legend=legend, figsize=figsize, ind_labels=ind_labels)
        # savefig(sp, "$pdf_path/$timeseries_name-SP-d$di.pdf")
    end

    # heatmap_plots(all_inds_cal_scores, all_inds_for_scores; pdf_path=outdir_nsgaii, plot_name="heatmaps_selected_$timeseries_name", png=true, ind_labels=ind_labels)
    return all_inds_cal_scores, all_inds_for_scores
end


function uniqueidx(x::AbstractArray{T}) where T
    uniqueset = Set{T}()
    ex = eachindex(x)
    idxs = Vector{eltype(ex)}()
    for i in ex
        xi = x[i]
        if !(xi in uniqueset)
            push!(idxs, i)
            push!(uniqueset, xi)
        end
    end
    idxs
end

function remove_duplicates(inds::Array{CGPInd}, fits)
    unique_indices = uniqueidx(fits)
    return inds[unique_indices], unique_indices
end

function remove_negative_fitness(inds::Array{CGPInd})
    # fits[map(x -> !(0 in (x .>= 0)), eachslice(fits, dims=1)), :]
    idx = map(x -> !(0 in (x.fitness .>= 0)), inds)
    return inds[idx], idx
end

#### load regular population inds
all_inds = Array{CGPInd}(undef, 0)
evolution_fits = Array{Float64}(undef, 20_000, cfg.d_fitness)
for run in readdir("data/models/gens")
    n_run = split(run, '-')[end]
    
    run_search = 500000
    gen_path = "data/models/gens/cgpshorefor-5sites-nsgaii-$n_run/$run_search"
    while !isdir(gen_path)
        run_search -= 25000
        gen_path = "data/models/gens/cgpshorefor-5sites-nsgaii-$n_run/$run_search"
    end
    println("$n_run: $gen_path")
    
    inds = Array{CGPInd}(undef, cfg.n_population)
    fitness = Array{Float64}(undef, cfg.n_population, cfg.d_fitness)

    i = 1
    if ispath(gen_path)
        for dna_file in readdir(gen_path)
            path_dna = "$gen_path/$dna_file"
            dna = JSON.parsefile(path_dna)
            ind_fit = dna["fitness"]
            chromo = convert(Array{Float64,1}, dna["chromosome"])
            ind = CGPInd(cfg, chromo)
            ind.fitness .= replace!(ind_fit, nothing => -Inf)
            inds[i] = ind
            fitness[i, :] .= ind_fit
            i += 1
            evolution_fits[i, :] .= ind_fit
        end
        push!(all_inds, inds...)
    else
        println("Doesn't exist: $gen_path")
    end
end

function heatmap_plots(cal_scores, for_scores; pdf_path="", plot_name="", plot_size=(950,300), right_margin=2.5mm, left_margin=0mm, png=false, ind_labels=[])
    if length(ind_labels) == 0
        ind_labels = map(x->"Ind-$x", collect(1:size(for_scores)[1]))
    end
    
    cal_scores = replace(cal_scores, -Inf=>NaN)
    for_scores = replace(for_scores, -Inf=>NaN)
    vmin = nanmin(nanminimum(cal_scores), nanminimum(for_scores))
    vmax = nanmax(nanmaximum(cal_scores), nanmaximum(for_scores))
    
    println("vmin: $vmin, vmax: $vmax")
    colors=palette(:rainbow, size(for_scores)[1])
    hm_for = heatmap(["GP", "NB", "DU", "TP", "TV"], ind_labels, for_scores, c=:rainbow, title="Forecast scores", clims=(vmin, vmax))
    hm_cal = heatmap(["GP", "NB", "DU", "TP", "TV"], ind_labels, cal_scores, c=:rainbow, title="Calibration scores", clims=(vmin, vmax))
    hm_plot = Plots.plot(hm_cal, hm_for, layout=Plots. grid(1,2, widths=(0.5, 0.5)), size=plot_size, right_margin=right_margin, left_margin=left_margin)
    if png
        savefig(hm_plot, "$pdf_path/$plot_name.png")
    else
        savefig(hm_plot, "$pdf_path/$plot_name.pdf")
    end
end

function print_minmax(cal_scores, for_scores)
    cal_scores = replace(cal_scores, -Inf=>NaN)
    for_scores = replace(for_scores, -Inf=>NaN)
    vmin = nanmin(nanminimum(cal_scores), nanminimum(for_scores))
    println("nanminimum(cal_scores), nanminimum(for_scores): $(nanminimum(cal_scores)) $(nanminimum(for_scores))")
    vmax = nanmax(nanmaximum(cal_scores), nanmaximum(for_scores))
    println("nanmaximum(cal_scores), nanmaximum(for_scores): $(nanmaximum(cal_scores)) $(nanmaximum(for_scores))")
    println("vmin: $vmin, vmax: $vmax")
end

function matrix_to_rows_array(mat)
    rows = Array{Array{Float64}}(undef, size(mat)[1])
    for i in 1:size(mat)[1]
        rows[i] = mat[i, :]
    end
    rows
end

outdir_nsgaii = "outputs/plots/EMS_plots"

cal_scores, all_scores = evaluate_calibration_forecast_all_individuals(cfg, all_inds)

heatmap_plots(cal_scores, all_scores; pdf_path=outdir_nsgaii, plot_name="heatmaps_raw", png=true)

println("size(all_inds): $(size(all_inds))")
println("remove_negative_fitness, remove_duplicates")
# all_inds, unique_idx = remove_duplicates(all_inds, matrix_to_rows_array(all_scores))
# all_inds, neg_idx = remove_negative_fitness(all_inds)
println("size(all_inds): $(size(all_inds))")

# cal_scores = evaluate_calibration_all_individuals(cfg, all_inds)
# all_scores = evaluate_forecast_all_individuals(cfg, all_inds)
cal_scores, all_scores = evaluate_calibration_forecast_all_individuals(cfg, all_inds)
heatmap_plots(cal_scores, all_scores; pdf_path=outdir_nsgaii, plot_name="heatmaps_reduced", png=true)
# all_scores = evaluate_calibration_all_individuals(cfg, all_inds)
mean_scores = Array{Float64}(undef, size(all_scores)[1])
min_scores = Array{Float64}(undef, size(all_scores)[1])
max_scores = Array{Float64}(undef, size(all_scores)[1])
median_scores = Array{Float64}(undef, size(all_scores)[1])

#### calibration sort
for i in 1:size(all_scores)[1]
    mean_scores[i] = nanmean(cal_scores[i, :])
    min_scores[i] = nanminimum(cal_scores[i, :])
    max_scores[i] = nanmaximum(cal_scores[i, :])
    median_scores[i] = nanmedian(cal_scores[i, :])
end
meansorted_inds, meanorder = sortrows(all_inds[:], mean_scores[:]); meansorted_inds = reverse(meansorted_inds[:]); reverse!(meanorder);
minsorted_inds, minorder = sortrows(all_inds[:], min_scores[:]); minsorted_inds = reverse(minsorted_inds[:]); reverse!(minorder);
maxsorted_inds, maxorder = sortrows(all_inds[:], max_scores[:]); maxsorted_inds = reverse(maxsorted_inds[:]); reverse!(maxorder);
mediansorted_inds, medianorder = sortrows(all_inds[:], median_scores[:]); mediansorted_inds = reverse(mediansorted_inds[:]); reverse!(medianorder);

heatmap_plots(cal_scores[meanorder, :], all_scores[meanorder, :]; pdf_path=outdir_nsgaii, plot_name="heatmaps_reduced_calibrationmeansort", png=true)

#### forecast sort
for i in 1:size(all_scores)[1]
    mean_scores[i] = nanmean(all_scores[i, :])
    min_scores[i] = nanminimum(all_scores[i, :])
    max_scores[i] = nanmaximum(all_scores[i, :])
    median_scores[i] = nanmedian(all_scores[i, :])
end
meansorted_inds, meanorder = sortrows(all_inds[:], mean_scores[:]); meansorted_inds = reverse(meansorted_inds[:]); reverse!(meanorder);
minsorted_inds, minorder = sortrows(all_inds[:], min_scores[:]); minsorted_inds = reverse(minsorted_inds[:]); reverse!(minorder);
maxsorted_inds, maxorder = sortrows(all_inds[:], max_scores[:]); maxsorted_inds = reverse(maxsorted_inds[:]); reverse!(maxorder);
mediansorted_inds, medianorder = sortrows(all_inds[:], median_scores[:]); mediansorted_inds = reverse(mediansorted_inds[:]); reverse!(medianorder);

heatmap_plots(cal_scores[meanorder, :], all_scores[meanorder, :]; pdf_path=outdir_nsgaii, plot_name="heatmaps_reduced_forecastmeansort", png=true)

prefix = "nf-" # "nf-"
algorithm = "NSGA-ii"  # "NS"

### Experts
n_experts = 2
best_inds = Array{CGPInd}(undef, n_experts, cfg.d_fitness)
site_acronyms = ["GP", "NB", "DU", "TP", "TV"]
site_cal_scores = []
site_for_scores = []
site_ind_labels = []
for i in collect(1:cfg.d_fitness)
    plot_inds = Array{CGPInd}(undef, 0)
    labels = []
    
    fitsorted_inds, fitorder = sortrows(all_inds[:], all_scores[:, i]); fitsorted_inds = fitsorted_inds[:];
    reverse!(fitsorted_inds); reverse!(fitorder);
    
    best_inds[:, i] .= fitsorted_inds[1:n_experts]
    j = 1
    for ind in best_inds[:, i]
        push!(plot_inds, ind)
        # push!(labels, "$(site_acronyms[i])-expert-$j")
        push!(labels, "CGP-expert-$(site_acronyms[i])-$j")
        j += 1
    end
    j = 1
    for ind in meansorted_inds[1:n_experts]
        push!(plot_inds, ind)
        push!(labels, "CGP-generalist-$j")
        j += 1
    end

    # evaluate_multiple_runs_dataset_v4_single_plots(plot_inds, datasets_calib, datasets_full, sites; pdf_path=outdir_nsgaii, timeseries_name="$(prefix)results-d$i-$(algorithm)", plot_size=(3280, 1280), linewidth=1, ind_labels=labels)
    cal_scores, for_scores = evaluate_multiple_runs_dataset_v4_single_plots(plot_inds, datasets_calib, datasets_full, sites; pdf_path=outdir_nsgaii, timeseries_name="$(prefix)results-d$i-$(algorithm)", plot_size=(3780, 1280), linewidth=1, ind_labels=labels)
    push!(site_cal_scores, cal_scores)
    push!(site_for_scores, for_scores)
    push!(site_ind_labels, labels)
end

# hm_expert_cal_perfs = Array{Float64, 2}(undef, length(site_ind_labels)+2, length(site_cal_scores))
# hm_expert_for_perfs = Array{Float64, 2}(undef, length(site_ind_labels)+2, length(site_for_scores))
# hm_expert_labels = Array{Any}(undef, length(site_ind_labels)+2)
# additional_ind_labels = ["ShoreFor", "Generalist-1"]
# for ind_i in 1:2
#     for s_i in 1:lastindex(site_cal_scores)
#         hm_expert_cal_perfs[ind_i, s_i] = site_cal_scores[ind_i][1, s_i]
#         hm_expert_for_perfs[ind_i, s_i] = site_for_scores[ind_i][1, s_i]
#     end
#     hm_expert_labels[ind_i] = additional_ind_labels[ind_i]
# end

# for ind_i in 1:lastindex(site_ind_labels)
#     for s_i in 1:lastindex(site_cal_scores)
#         hm_expert_cal_perfs[ind_i+2, s_i] = site_cal_scores[ind_i][1, s_i]
#         hm_expert_for_perfs[ind_i+2, s_i] = site_for_scores[ind_i][1, s_i]
#     end
#     hm_expert_labels[ind_i+2] = site_ind_labels[ind_i][1]
# end
# reverse!(hm_expert_cal_perfs, dims=1)
# reverse!(hm_expert_for_perfs, dims=1)
# reverse!(hm_expert_labels)
# heatmap_plots(hm_expert_cal_perfs, hm_expert_for_perfs; pdf_path=outdir_nsgaii, plot_name="heatmaps_topexperts", png=true, ind_labels=hm_expert_labels)

### Graphs
# Experts
for i in 1:cfg.d_fitness
    for j in 1:n_experts
        chromo_draw(best_inds[j, i], "$outdir_nsgaii/$(prefix)expert-d$i-graph-$(algorithm)-ind$j")
    end
end
# General models
for i in 1:n_experts
    chromo_draw(meansorted_inds[i], "$outdir_nsgaii/$(prefix)generalmodel-graph-$(algorithm)-ind$i")
end
