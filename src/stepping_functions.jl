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

    AgentT = NamedTuple{(:pos, :home_center), Tuple{Tuple{Int,Int}, Tuple{Int,Int}}}
    new_agents = AgentT[]

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
