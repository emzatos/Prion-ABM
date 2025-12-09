cd(@__DIR__) #src
using Agents, Random
using Agents.DataFrames, Agents.Graphs
using DrWatson: @dict
using CairoMakie
using Distributions: truncated, Normal, mean
using GLMakie


@agent struct DeerAgent(GridAgent{2})
    status::Symbol              # :S, :E, :I, :C
    weeks_in_state::Int         # weeks since entering current state
    state_duration::Int         # sampled duration for current state
    home_center::Tuple{Int,Int} # home range center for movement fidelity
    infection_source::Symbol    # :none, :direct, :indirect, :initial
    infection_tick::Int         # when infected (-1 if never)
    lifespan::Int 
end


function sample_duration(rng, min_weeks::Int, max_weeks::Int)
    μ = (min_weeks + max_weeks) / 2
    σ = (max_weeks - min_weeks) / 4
    d = truncated(Normal(μ, σ), min_weeks, max_weeks)
    return round(Int, rand(rng, d))
end

function model_initiation(;
    nx=100,
    ny=100,
    n_deer=5000,
    n_infected=10,
    exposure_min=15,
    exposure_max=35,
    infectious_min=25,
    infectious_max=44,
    clinical_min=1,
    clinical_max=36,
    prion_survival=0.9967,    # ~4 year half-life
    shedding_rate=1.0,
    carcass_load=100.0,
    # Transmission parameters
    bd=5.2e-6,               # direct transmission coefficient
    bi=5.5e-8,               # indirect transmission coefficient
    k=0.1,                    # aggregation parameter (0.01, 1, or 10000)
    ε=0.0001,                    # 0=density-dependent, 1=frequency-dependent
    transmission_radius=2,
    # Demographics
    weekly_survival=0.991,
    clinical_mortality_mult=2.0,
    annual_birth_rate=0.55,
    offspring_per_birth=1.7,
    carrying_capacity=20000,
    seed=0,
)
    rng = Xoshiro(seed)

    properties = @dict(
        exposure_min,
        exposure_max,
        infectious_min,
        infectious_max,
        clinical_min,
        clinical_max,
        prion_survival,
        shedding_rate,
        carcass_load,
        bd,
        bi,
        k,
        ε,
        transmission_radius,
        weekly_survival,
        clinical_mortality_mult,
        annual_birth_rate,
        offspring_per_birth,
        carrying_capacity,
        V = zeros(Float64, nx, ny),
        V_cumulative = zeros(Float64, nx, ny),
        tick = 0,
        total_direct = 0,
        total_indirect = 0,
        prev_total_direct = 0,
        prev_total_indirect = 0,
        events = NamedTuple{(:tick, :source, :agent_id, :pos),Tuple{Int,Symbol,Int,Tuple{Int,Int}}}[],
        weekly_direct = Int[],
        weekly_indirect = Int[],
        weekly_prevalence = Float64[],
        weekly_population = Int[],
        weekly_V_total = Float64[],
        weekly_S = Int[],
        weekly_E = Int[],
        weekly_I = Int[],
        weekly_C = Int[],
    )

    space = GridSpace((nx, ny); periodic=false)
    model = StandardABM(DeerAgent, space; agent_step!, model_step!, properties, rng)

    # Initialize deer population
    for i in 1:n_deer
        home = (rand(abmrng(model), 1:nx), rand(abmrng(model), 1:ny))
        if i <= n_infected
            status = :I
            duration = sample_duration(abmrng(model), infectious_min, infectious_max)
            source = :initial
            inf_tick = 0
        else
            status = :S
            duration = 0
            source = :none
            inf_tick = -1
        end

        add_agent!(home, DeerAgent, model;
            status,
            weeks_in_state=0,
            state_duration=duration,
            home_center=home,
            infection_source=source,
            infection_tick=inf_tick,
            lifespan=0
        )
    end

    return model
end


function agent_step!(agent, model)
    move_deer!(agent, model)
    transmit!(agent, model)
    progress_disease!(agent, model)
    shed_prions!(agent, model)
end

function model_step!(model)
    model.tick += 1

    # Prion decay
    model.V .*= model.prion_survival

    # Survival check (background + disease mortality)
    to_remove = Int[]
    n_C = 0
    n_E = 0
    n_S = 0
    n_I = 0
    for agent in allagents(model)
        if agent.status == :C
            n_C += 1
        end
        if agent.status == :E
            n_E += 1
        end
        if agent.status == :S
            n_S += 1
        end
        if agent.status == :I
            n_I += 1
        end
        if !check_survival!(agent, model)
            push!(to_remove, agent.id)
        end
    end
    for id in to_remove
        remove_agent!(id, model)
    end

    # Reproduction
    reproduction!(model)

    # Record weekly statistics
    push!(model.weekly_direct, model.total_direct - model.prev_total_direct)
    push!(model.weekly_indirect, model.total_indirect - model.prev_total_indirect)

    model.prev_total_direct = model.total_direct
    model.prev_total_indirect = model.total_indirect

    n = nagents(model)
    push!(model.weekly_population, n)

    if n > 0
        push!(model.weekly_prevalence, (n_E + n_I + n_C) / n)
        push!(model.weekly_S, n_S)
        push!(model.weekly_E, n_E)
        push!(model.weekly_I, n_I)
        push!(model.weekly_C, n_C)
    else
        push!(model.weekly_prevalence, 0.0)
        push!(model.weekly_S, 0)
        push!(model.weekly_E, 0)
        push!(model.weekly_I, 0)
        push!(model.weekly_C, 0)
    end

    push!(model.weekly_V_total, sum(model.V))
end



function move_deer!(agent, model)
    randomwalk!(agent, model)
end


function compute_local_lambda(agent, model)
    nx, ny = size(model.V)
    pos = agent.pos

    # Count local infectious deer (Z) and local population (N)
    Z_local = 0
    N_local = 0

    ids = nearby_ids(agent, model, model.transmission_radius)
    Z_local = count(id -> model[id].status in (:I, :C), ids)

    # Sum local environmental prion load
    V_local = model.V[pos...]
    # for dx in -model.V_radius:model.V_radius
    #     for dy in -model.V_radius:model.V_radius
    #         x, y = pos[1] + dx, pos[2] + dy
    #         if 1 <= x <= nx && 1 <= y <= ny
    #             V_local += model.V[x, y]
    #         end
    #     end
    # end

    N_local = max(N_local, 1)  # avoid division by zero

    k = model.k
    ε = model.ε
    denom = (1 - ε) + ε * N_local

    # Separate λ for source attribution
    λ_direct = (Z_local > 0 && model.bd > 0) ?
               k * log(1 + (model.bd * Z_local) / k) / denom : 0.0
    λ_indirect = (V_local > 0 && model.bi > 0) ?
                 k * log(1 + (model.bi * V_local) / k) / denom : 0.0

    return λ_direct, λ_indirect
end

function transmit!(agent, model)
    agent.status == :S || return

    λ_direct, λ_indirect = compute_local_lambda(agent, model)

    p_direct = 1 - exp(-λ_direct)
    p_indirect = 1 - exp(-λ_indirect)

    direct_infection = rand(abmrng(model)) < p_direct
    indirect_infection = rand(abmrng(model)) < p_indirect

    if !direct_infection && !indirect_infection
        return  # no infection
    end

    if direct_infection && !indirect_infection
        source = :direct
    elseif !direct_infection && indirect_infection
        source = :indirect
    else
        # Both occurred - choose based on relative probabilities
        total_p = p_direct + p_indirect
        if rand(abmrng(model)) < p_direct / total_p
            source = :direct
        else
            source = :indirect
        end
    end

    model.total_direct += (source == :direct)
    model.total_indirect += (source == :indirect)

    # Transition to Exposed
    agent.status = :E
    agent.weeks_in_state = 0
    agent.state_duration = sample_duration(abmrng(model), model.exposure_min, model.exposure_max)
    agent.infection_source = source
    agent.infection_tick = model.tick

    # Record event
    push!(model.events, (tick=model.tick, source=source, agent_id=agent.id, pos=agent.pos))
end



function progress_disease!(agent, model)
    agent.status == :S && return

    agent.weeks_in_state += 1

    if agent.weeks_in_state >= agent.state_duration
        if agent.status == :E
            # Transition to Infectious
            agent.status = :I
            agent.weeks_in_state = 0
            agent.state_duration = sample_duration(abmrng(model), model.infectious_min, model.infectious_max)

        elseif agent.status == :I
            # Transition to Clinical
            agent.status = :C
            agent.weeks_in_state = 0
            agent.state_duration = sample_duration(abmrng(model), model.clinical_min, model.clinical_max)

        elseif agent.status == :C
            # End of clinical = death (handled in check_survival!)
            agent.state_duration = 0  # signal for removal
        end
    end
end



function shed_prions!(agent, model)
    if agent.status in (:I, :C)
        x, y = agent.pos
        model.V[x, y] += model.shedding_rate
        model.V_cumulative[x, y] += model.shedding_rate
    end
end

function deposit_carcass!(agent, model)
    x, y = agent.pos
    model.V[x, y] += model.carcass_load
    model.V_cumulative[x, y] += model.carcass_load
end



function check_survival!(agent, model)
    # Disease death at end of clinical phase
    if agent.status == :C && agent.weeks_in_state >= agent.state_duration
        deposit_carcass!(agent, model)
        return false
    end

    # Background mortality (elevated during clinical)
    survival = model.weekly_survival
    if agent.status == :C
        survival /= model.clinical_mortality_mult
    end

    if rand(abmrng(model)) >= survival
        if agent.status in (:I, :C)
            deposit_carcass!(agent, model)
        end
        return false
    end

    return true
end

function reproduction!(model)
    nx, ny = size(model.V)
    N = nagents(model)
    K = model.carrying_capacity

    # Density-dependent probability of reproducing (from paper)
    pr_reproduce = model.annual_birth_rate / (1 + (N / K)^2)
    weekly_pr = pr_reproduce / 52
    expected_offspring = weekly_pr * model.offspring_per_birth

    new_agents = []

    for agent in allagents(model)
        if rand(abmrng(model)) < expected_offspring
            # Offspring inherits home range with small variation
            home_offset = (rand(abmrng(model), -3:3), rand(abmrng(model), -3:3))
            new_home = (
                clamp(agent.home_center[1] + home_offset[1], 1, nx),
                clamp(agent.home_center[2] + home_offset[2], 1, ny),
            )

            push!(new_agents, (pos=new_home, home_center=new_home))
        end
    end

    for na in new_agents
        add_agent!(na.pos, DeerAgent, model;
            status=:S,
            weeks_in_state=0,
            state_duration=0,
            home_center=na.home_center,
            infection_source=:none,
            infection_tick=-1,
        )
    end
end



function compute_metrics(model)
    total = model.total_direct + model.total_indirect
    n_final = nagents(model)

    final_prev = n_final > 0 ? count(a -> a.status != :S, allagents(model)) / n_final : NaN
    direct_frac = total > 0 ? model.total_direct / total : NaN
    peak_prev = isempty(model.weekly_prevalence) ? 0.0 : maximum(model.weekly_prevalence)
    init_pop = isempty(model.weekly_population) ? n_final : first(model.weekly_population)
    pop_frac = init_pop > 0 ? n_final / init_pop : 0.0

    return (
        total_direct=model.total_direct,
        total_indirect=model.total_indirect,
        direct_fraction=direct_frac,
        final_prevalence=final_prev,
        peak_prevalence=peak_prev,
        final_population=n_final,
        population_fraction=pop_frac,
        weeks_run=model.tick,
    )
end

function temporal_analysis(model; window_weeks::Int=52)
    events = model.events
    max_tick = model.tick
    windows = NamedTuple[]

    for start_tick in 1:window_weeks:max_tick
        end_tick = min(start_tick + window_weeks - 1, max_tick)
        window_events = filter(e -> start_tick <= e.tick <= end_tick, events)

        n_direct = count(e -> e.source == :direct, window_events)
        n_indirect = count(e -> e.source == :indirect, window_events)
        total = n_direct + n_indirect

        push!(windows, (
            start_week=start_tick,
            end_week=end_tick,
            year=(start_tick - 1) ÷ 52 + 1,
            n_direct=n_direct,
            n_indirect=n_indirect,
            indirect_fraction=total > 0 ? n_indirect / total : NaN,
        ))
    end
    return windows
end

function spatial_infection_maps(model)
    dims = size(model.V)
    direct_map = zeros(Int, dims)
    indirect_map = zeros(Int, dims)

    for e in model.events
        if e.source == :direct
            direct_map[e.pos...] += 1
        elseif e.source == :indirect
            indirect_map[e.pos...] += 1
        end
    end
    return direct_map, indirect_map
end



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

    ax8 = Axis(fig[4, 2], xlabel="X", ylabel="Y", title="Cumulative Prion Deposition (log)", aspect=DataAspect())
    hm2 = heatmap!(ax8, log1p.(model.V_cumulative), colormap=:viridis)
    Colorbar(fig[4, 2][1, 2], hm2)

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


function runner()
    susceptible(x) = count(i == :S for i in x)
    exposed(x) = count(i == :E for i in x)
    infected(x) = count(i == :I for i in x)
    clinical(x) = count(i == :C for i in x)

    CairoMakie.activate!()

    println("Initializing model...")
    model = model_initiation(;
        nx=100,
        ny=100,
        n_deer=5000,
        n_infected=10,
        # Prion dynamics
        prion_survival=0.9967,    # ~4 year half-life
        shedding_rate=1.0,
        carcass_load=100.0,
        # Transmission parameters
        bd=9.2e-4,               # direct transmission coefficient
        bi=5.5e-5,               # indirect transmission coefficient
        k=10000,                    # aggregation parameter (0.01, 1, or 10000)
        ε=0.0001,                    # 0=density-dependent, 1=frequency-dependent
        transmission_radius=2,
    )

    # println("Running simulation for 10 years (520 weeks)...")
    # for week in 1:520
    #     step!(model, 1)
    #     if week % 52 == 0
    #         m = compute_metrics(model)
    #         println("  Year $(week÷52): Pop=$(m.final_population), " *
    #                 "Prev=$(round(m.peak_prevalence*100, digits=1))%, " *
    #                 "Direct=$(m.total_direct), Indirect=$(m.total_indirect)")
    #     end
    # end

    # # Final metrics
    # metrics = compute_metrics(model)
    # println("\nFinal Results:")
    # println("  Population: $(metrics.final_population) ($(round(metrics.population_fraction*100, digits=1))% remaining)")
    # println("  Peak prevalence: $(round(metrics.peak_prevalence*100, digits=1))%")
    # println("  Total infections: $(metrics.total_direct + metrics.total_indirect)")
    # println("  Direct: $(metrics.total_direct), Indirect: $(metrics.total_indirect)")
    # if !isnan(metrics.direct_fraction)
    #     mode = metrics.direct_fraction > 0.5 ? "DIRECT" : "INDIRECT"
    #     println("  Dominant mode: $mode ($(round(metrics.direct_fraction*100, digits=1))% direct)")
    # end

    # # Temporal analysis
    # println("\nTemporal transmission analysis (by year):")
    # temporal = temporal_analysis(model; window_weeks=52)
    # for w in temporal
    #     if !isnan(w.indirect_fraction) && (w.n_direct + w.n_indirect) > 0
    #         println("  Year $(w.year): Direct=$(w.n_direct), Indirect=$(w.n_indirect), " *
    #                 "Indirect fraction=$(round(w.indirect_fraction*100, digits=0))%")
    #     end
    # end

    # # Generate plots
    # println("\nGenerating dashboard...")
    # fig1 = plot_dashboard(model)
    # save("cwd_dashboard.png", fig1)

    # println("Generating transmission analysis...")
    # fig2 = plot_transmission_analysis(model)
    # save("cwd_transmission_analysis.png", fig2)

    # Video generation (optional - uncomment to use)

    # println("Generating video...")
    model_video = model_initiation(       nx=100,
        ny=100,
        n_deer=5000,
        n_infected=10,
        # Prion dynamics
        prion_survival=0.9967,    # ~4 year half-life
        shedding_rate=1.0,
        carcass_load=100.0,
        # Transmission parameters
        #bd=9.2e-4,               # direct transmission coefficient
        #bi=5.5e-5,   
        bd=0.017,
        bi = 0.0003,            # indirect transmission coefficient
        k=10000,                    # aggregation parameter (0.01, 1, or 10000)
        ε=0.0001,                    # 0=density-dependent, 1=frequency-dependent
        transmission_radius=2)
    # abmvideo("cwd_evolution.mp4", model_video;
    #     agent_step! = agent_step!,
    #     model_step! = model_step!,
    #     agent_color = seic_color,
    #     agent_size = 4,
    #     frames = 2600,
    #     framerate = 100,
    #     title = "CWD Spatial ABM",
    #     heatarray = :V,
    #     heatkwargs = (colormap = :inferno, colorrange = (0, 50)),
    # )


    # Interactive exploration (requires GLMakie - uncomment to use)

    GLMakie.activate!()
    model_explore = model_initiation(seed = 456)
    params = Dict(
        :bd => 1e-4:1e-4:2e-3,
        :bi => 1e-7:1e-7:2e-3,
        :prion_survival => 0.9:0.01:0.9967,
    )

    fig, abmobs = abmexploration(model_video;
        agent_color = seic_color,
        agent_size = 4,
        heatarray = :V,
        heatkwargs = (colormap = :inferno, colorrange = (0, 50)),
        params = params
    )
    fig


    println("\nDone!")
end

#     using GLMakie
#     GLMakie.activate!()
#     model_explore = model_initiation(seed = 456 + randn(Int))
#     params = Dict(
#         :bd => 1e-4:1e-4:2e-3,
#         :bi => 1e-7:1e-7:2e-3,
#         :prion_survival => 0.9:0.01:0.9967,
#     )

#     healthy(a) = a.status == :S
#     unhealthy(a) = a.status in (:E, :I, :C)

#     adata = [(healthy, count), (unhealthy, count)]

#     weeklydirect(m) = mean(m.total_direct)
#     weeklyindirect(m) = mean(m.total_indirect)

#     mdata = [weeklydirect, weeklyindirect]

#     fig, abmobs = abmexploration(model_explore;
#         agent_color = seic_color,
#         agent_size = 4,
#         heatarray = :V,
#         heatkwargs = (colormap = :inferno, colorrange = (0, 50)),
#         params = params,
#         adata, alabels = ["Healthy", "Unhealthy"],
#         mdata, mlabels = ["Avg Weekly Direct Infections", "Avg Weekly Indirect Infections"],
#     )
#     fig
# GLMakie.closeall()

runner()