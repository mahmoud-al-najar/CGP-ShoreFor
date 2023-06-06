using CartesianGeneticProgramming
using Random
using JSON

function genes_to_chromo(cfg::NamedTuple, genes, outputs)
    genes, outputs = deepcopy(genes), deepcopy(outputs)
    R = cfg.rows
    C = cfg.columns
    P = cfg.n_parameters

    maxs = get_maxs(cfg)
    # genes to chromo
    x_chromo = (genes[:, :, 1] ./ maxs)[:]
    y_chromo = (genes[:, :, 2] ./ maxs)[:]
    f_chromo = (genes[:, :, 3] ./ length(cfg.functions))[:]
    p_chromo = (genes[:, :, 4:3+P] ./cfg.param_max)[:]
    # outputs
    o_chromo = outputs ./ (R * C + cfg.n_in)
    # order
    chromosome = 0 * ones(cfg.rows * cfg.columns * (3 + cfg.n_parameters) + cfg.n_out)
    chromosome[1:C] .= x_chromo
    chromosome[C+1:2C] .= y_chromo
    chromosome[2C+1:3C] .= f_chromo
    chromosome[3C+1:R*C*(3+P)] = p_chromo
    chromosome[(R*C*(3+P)+1):end] .= o_chromo
    chromosome
end

function get_ind_from_template(cfg, model_template, inputs, outputs)
    # global buffer, ids, fs, xs, ys, ps, indices
    rand_ind = CGPInd(cfg)

    ids = []
    fs = []
    xs = []
    ys = []
    ps = []
    for r in model_template
        id = r.id
        f = r.f
        x = r.x
        y = r.y
        p = haskey(r, :p) ? r.p : []
        
        push!(ids, id)
        push!(fs, f)
        push!(xs, x)
        push!(ys, y)
        push!(ps, p)
    end

    buffer = Array{String}(undef, cfg.columns)
    buffer .= "-"
    buffer = vcat(inputs, buffer)
    
    lnrng = LinRange(cfg.n_in+1,cfg.columns, length(ids))
    indices = unique(Int.(ceil.(lnrng)))
    # println(indices)
    @assert length(indices) == length(ids)
    buffer[indices] .= ids
    
    xs = map((x1) -> findfirst(x2 -> x2==x1, buffer), xs)
    ys = map((x1) -> findfirst(x2 -> x2==x1, buffer), ys)
    fs = map((x1) -> findfirst(x2 -> x2==x1, string.(cfg.functions)), fs)
    
    # fill empty y's with random values
    if cfg.recur == 0.0
        map((i) -> ys[i] = rand(1:i), findall(x -> isequal(x, nothing), ys))
    else
        ys[ys.==nothing] .= 1.0
    end

    rand_ind.genes[:, indices .- cfg.n_in, 1] .= xs'
    rand_ind.genes[:, indices .- cfg.n_in, 2] .= ys'
    rand_ind.genes[:, indices .- cfg.n_in, 3] .= fs'

    outputs = map((x1) -> findfirst(x2 -> x2==x1, buffer), outputs)

    chromo = genes_to_chromo(cfg, rand_ind.genes, outputs)

    ind = CGPInd(cfg, chromo)
    
    # @assert sum(ind.genes[:,indices .- cfg.n_in,:] .!= rand_ind.genes[:,indices .- cfg.n_in,:]) == 0
    return ind
end

function get_template_from_ind(cfg, ind, inputs)
    template = []
    ot = get_output_trace(ind)
    nodes = ind.nodes[ot]
    output_nodes = ind.nodes[ind.outputs]
    outputs_ids = []
    i = 1
    for node in nodes # start template from inputs 
        f = node.f
        x = node.x
        y = node.y
        id = string(f, "-", x, "-", y, "--oti:", i)
        y_id = ""
        # println(node)

        if node in output_nodes
            push!(outputs_ids, id)
        end

        if string(f) != "f_null"
            ar = CGPFunctions.arity[string(f)]
            # println(ot)
            node_x = ind.nodes[x]
            x_id = string(node_x.f, "-", node_x.x, "-", node_x.y)
            if ar > 1
                if y <= cfg.n_in
                    y_id = inputs[y]
                else
                    node_y = ind.nodes[y]
                    y_id = string(node_y.f, "-", node_y.x, "-", node_y.y)
                end
            end
        else
            ar = 0
            x_id = inputs[ot[i]]
        end
        row = (id=id, f=string(f), x=x_id, y=y_id)
        push!(template, row)
        i += 1
    end

    return template, outputs_ids
end

function get_template_from_dna(cfg, path_dna, inputs)
    dna = JSON.parsefile(path_dna)
    chromo = convert(Array{Float64,1}, dna["chromosome"])
    ind = CGPInd(cfg, chromo);
    return get_template_from_ind(cfg, ind, inputs)
end

function print_dna(cfg, path_dna, inputs)
    dna = JSON.parsefile(path_dna)
    chromo = convert(Array{Float64,1}, dna["chromosome"])
    fs = dna["fs"]
    xs = dna["xs"]
    ys = dna["ys"]
    @assert length(fs) == length(xs) == length(ys)

    buffer = []
    for i in 1:length(inputs)
        s = "$i\tinput:$(inputs[i])"
        push!(buffer, s)
    end

    for i in 1:length(fs)
        func = cfg.functions[Int(fs[i])]
        x = xs[i]
        if CGPFunctions.arity[string(func)] > 1
            y = ys[i]
            push!(buffer, "$(i + cfg.n_in)\t$func\t$x\t$y")
        else
            push!(buffer, "$(i + cfg.n_in)\t$func\t$x")
        end
    end

    ind = CGPInd(cfg, chromo)
    buffer = buffer[get_output_trace(ind)]
    for a in reverse(buffer)
        println(a)
    end
    println("outputs: $(ind.outputs)")
end
