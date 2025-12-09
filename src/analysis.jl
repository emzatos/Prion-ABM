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
