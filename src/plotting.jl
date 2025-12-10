seic_color(a) =
    a.status == :S ? :blue :
    a.status == :E ? :orange :
    a.status == :I ? :red :
    a.status == :C ? :purple : :black

function plot_dashboard(model)
    fig = Figure(size=(1600, 1200), fontsize=12)

    weeks = 1:length(model.weekly_prevalence)
    years = weeks ./ 52

    # Row 1: Population and Prevalence
    ax1 = Axis(fig[1, 1], xlabel="Years", ylabel="Population", title="Host Population")
    lines!(ax1, years, model.weekly_population, color=:steelblue, linewidth=2)

    ax2 = Axis(fig[1, 2], xlabel="Years", ylabel="Prevalence (%)", title="Disease Prevalence")
    lines!(ax2, years, model.weekly_prevalence .* 100, color=:firebrick, linewidth=2)

    # Row 2: SEIC compartments and Environmental reservoir
    ax3 = Axis(fig[2, 1], xlabel="Years", ylabel="Count", title="Disease Compartments")
    lines!(ax3, years, model.weekly_S, label="S", color=:green)
    lines!(ax3, years, model.weekly_E, label="E", color=:orange)
    lines!(ax3, years, model.weekly_I, label="I", color=:red)
    lines!(ax3, years, model.weekly_C, label="C", color=:purple)
    axislegend(ax3, position=:rt)

    ax4 = Axis(fig[2, 2], xlabel="Years", ylabel="Total Prions", title="Environmental Reservoir")
    lines!(ax4, years, model.weekly_V_total, color=:darkgreen, linewidth=2)

    # Row 3: Cumulative infections by source
    ax5 = Axis(fig[3, 1], xlabel="Years", ylabel="Cumulative", title="Cumulative Infections by Source")
    lines!(ax5, years, cumsum(model.weekly_direct), label="Direct", color=:coral, linewidth=2)
    lines!(ax5, years, cumsum(model.weekly_indirect), label="Indirect", color=:teal, linewidth=2)
    axislegend(ax5, position=:lt)

    # Rolling transmission rate
    window = 13
    if length(model.weekly_direct) >= window
        roll_d = [mean(model.weekly_direct[max(1, i - window + 1):i]) for i in 1:length(model.weekly_direct)]
        roll_i = [mean(model.weekly_indirect[max(1, i - window + 1):i]) for i in 1:length(model.weekly_indirect)]

        ax6 = Axis(fig[3, 2], xlabel="Years", ylabel="Weekly (13-wk avg)", title="Infection Rate by Source")
        lines!(ax6, years, roll_d, label="Direct", color=:coral, linewidth=2)
        lines!(ax6, years, roll_i, label="Indirect", color=:teal, linewidth=2)
        axislegend(ax6, position=:rt)
    end

    # Row 4: Heatmaps
    ax7 = Axis(fig[4, 1], xlabel="X", ylabel="Y", title="Current Prion Distribution", aspect=DataAspect())
    hm1 = heatmap!(ax7, model.V, colormap=:YlOrRd)
    Colorbar(fig[4, 1][1, 2], hm1)

    # ax8 = Axis(fig[4, 2], xlabel="X", ylabel="Y", title="Cumulative Prion Deposition (log)", aspect=DataAspect())
    # hm2 = heatmap!(ax8, log1p.(model.V_cumulative), colormap=:viridis)
    # Colorbar(fig[4, 2][1, 2], hm2)

    ax9 = Axis(fig[4, 2], xlabel="X", ylabel="Y", title="Population Density", aspect=DataAspect())
    hm3 = heatmap!(ax9, model.population_grid, colormap=:lajolla)
    Colorbar(fig[4, 2][1, 2], hm3)

    return fig
end

function plot_transmission_analysis(model)
    temporal = temporal_analysis(model; window_weeks=52)
    direct_map, indirect_map = spatial_infection_maps(model)

    fig = Figure(size=(1600, 800), fontsize=12)

    # Temporal: indirect fraction over time
    years = [w.year for w in temporal]
    indirect_frac = [w.indirect_fraction for w in temporal]
    valid_idx = findall(!isnan, indirect_frac)

    ax1 = Axis(fig[1, 1], xlabel="Year", ylabel="Indirect Fraction",
        title="Proportion from Indirect Transmission Over Time")

    if !isempty(valid_idx)
        scatter!(ax1, years[valid_idx], indirect_frac[valid_idx], color=:steelblue, markersize=10)

        if length(valid_idx) > 2
            x = Float64.(years[valid_idx])
            y = indirect_frac[valid_idx]
            slope = sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^ 2)
            intercept = mean(y) - slope * mean(x)
            x_line = range(minimum(x), maximum(x), length=100)
            lines!(ax1, x_line, intercept .+ slope .* x_line, color=:red, linestyle=:dash,
                label="Trend (slope=$(round(slope, digits=3))/yr)")
            axislegend(ax1, position=:rb)
        end
    end
    hlines!(ax1, [0.5], color=:gray, linestyle=:dot)
    ylims!(ax1, 0, 1)

    # Stacked annual infections
    n_direct = [w.n_direct for w in temporal]
    n_indirect = [w.n_indirect for w in temporal]

    ax2 = Axis(fig[1, 2], xlabel="Year", ylabel="Infections", title="Annual Infections by Source")
    band!(ax2, years, zeros(length(years)), Float64.(n_direct), color=(:coral, 0.7), label="Direct")
    band!(ax2, years, Float64.(n_direct), Float64.(n_direct .+ n_indirect), color=(:teal, 0.7), label="Indirect")
    axislegend(ax2, position=:rt)

    # Spatial heatmaps
    ax3 = Axis(fig[2, 1], xlabel="X", ylabel="Y", title="Direct Infections", aspect=DataAspect())
    hm3 = heatmap!(ax3, direct_map, colormap=:Reds)
    Colorbar(fig[2, 1][1, 2], hm3)

    ax4 = Axis(fig[2, 2], xlabel="X", ylabel="Y", title="Indirect Infections", aspect=DataAspect())
    hm4 = heatmap!(ax4, indirect_map, colormap=:Blues)
    Colorbar(fig[2, 2][1, 2], hm4)

    # Ratio map
    total_map = direct_map .+ indirect_map
    ratio_map = similar(total_map, Float64)
    for i in eachindex(total_map)
        ratio_map[i] = total_map[i] > 0 ? indirect_map[i] / total_map[i] : NaN
    end

    ax5 = Axis(fig[2, 3], xlabel="X", ylabel="Y", title="Indirect Fraction by Location", aspect=DataAspect())
    hm5 = heatmap!(ax5, ratio_map, colormap=:RdBu, colorrange=(0, 1))
    Colorbar(fig[2, 3][1, 2], hm5)

    return fig
end
