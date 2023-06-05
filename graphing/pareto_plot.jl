using Plots
using Plots.PlotMeasures
using LaTeXStrings
plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing, grid=true, guidefontsize=14, tickfontsize=12, dpi=200, left_margin=10mm) #labelfontsize=14, 

function evolution_to_fits(e::NSEvolution)
    fits = Array{Float64}(undef, e.config.n_population, e.config.d_fitness)
    for i in eachindex(e.population)
        fits[i, :] .= e.population[i].fitness
    end
    fits
end

function evolution_archive_fits_and_novelty(e::NSEvolution)
    fits = Array{Float64}(undef, length(e.inds_archive), e.config.d_fitness)
    archive = reverse(e.inds_archive)
    for i in eachindex(archive)
        fits[i, :] .= archive[i].fitness
    end
    fits, reverse(e.novelty_archive)
end

function pcp_plot(e::NSEvolution; plot_size=(720, 480), pdf_path="")
    fits = evolution_to_fits(e)
    p1 = Plots.plot()
    # fits: (n_population, d_fitness)
    
    colors=palette(:rainbow, size(fits)[1])
    # colors = palette(:rainbow, unique(e.novelty_archive))
    # novelty_colors = 

    xs = collect(1:size(fits)[2])
    stds = []
    medians = []
    in_archive = Array{Int}(undef, size(fits)[1])
    for i in collect(1:size(fits)[1])
        push!(medians, median(fits[i, :]))
        push!(stds, std(fits[i, :]))
        in_archive[i] = fits[i, :] in map(x->x.fitness, e.inds_archive)
        Plots.plot!(p1, xs, fits[i, :], label="ind$i", linecolor=colors[i], linealpha=0.9)
        # Plots.plot!(p1, xs, fits[i, :], label="ind$i", linealpha=0.9)
    end
    
    p2 = scatter(medians, stds, xlabel="median(fitness)", ylabel="std(fitness)", yflip=true, marker_z=in_archive, xlim=(-0.5, 1.0), ylim=(0.0, 1.0))

    p3 = heatmap(["GP", "NB", "DU", "TP", "TV"], map(x->"Ind-$x", collect(1:size(fits)[1])), fits, c=colors)
    
    p4 = Plots.plot(p1, p2, p3, layout=@layout([C B A]), plot_title="Gen-$(e.gen) population", dpi=200, left_margin=5mm, bottom_margin=5mm, legend=false)  # , size=plot_size
    
    if pdf_path != ""
        println("saving")
        savefig(p4, pdf_path)
    else
        return p4
    end
end

function archive_plot(e::NSEvolution; plot_size=(720, 480), pdf_path="")
    fits, novelty = evolution_archive_fits_and_novelty(e)
    p1 = Plots.plot()
    # fits: (n_population, d_fitness)
    if size(fits)[1] > 1
        colors=palette(:rainbow, size(fits)[1])
    else
        colors = [:red]
    end
    xs = collect(1:size(fits)[2])
    stds = []
    medians = []
    for i in collect(1:size(fits)[1])
        push!(medians, median(fits[i, :]))
        push!(stds, std(fits[i, :]))
        # Plots.plot!(p1, xs, fits[i, :], label="ind$i", linecolor=colors[i], linealpha=0.9)
        Plots.plot!(p1, xs, fits[i, :], label="ind$i", linealpha=0.9, line_z=novelty[i])
    end

    # nov = e.novelty_archive
    # nov[nov.==-Inf] .= 0.0
    p2 = scatter(medians, stds, xlabel="median(fitness)", ylabel="std(fitness)", yflip=true, marker_z=novelty, colorbar_title="Novelty", xlim=(-0.5, 1.0), ylim=(0.0, 1.0))
    # println(e.novelty_archive)
    p3 = heatmap(["GP", "NB", "DU", "TP", "TV"], map(x->"Ind-$x", collect(1:size(fits)[1])), fits, c=colors, colorbar_title="Mielke")
    
    p4 = Plots.plot(p1, p2, p3, layout=@layout([C B A]), plot_title="Gen-$(e.gen) archive", dpi=200, left_margin=5mm, bottom_margin=5mm, legend=false, size=plot_size)
    
    if pdf_path != ""
        println("saving")
        savefig(p4, pdf_path)
    else
        return p4
    end
end

function population_and_archive(e::NSEvolution; plot_size=(1380, 320), pdf_path="", dpi=300)
    p1 = pcp_plot(e, plot_size=plot_size)
    p2 = archive_plot(e, plot_size=plot_size)
    
    p3 = Plots.plot(p1, p2, layout=@layout([E ; D]), plot_title="Gen-$(e.gen)", dpi=dpi, legend=false, size=plot_size)  # , left_margin=5mm, bottom_margin=10mm
    
    if pdf_path != ""
        println("saving")
        savefig(p3, pdf_path)
    else
        return p3
    end
end
