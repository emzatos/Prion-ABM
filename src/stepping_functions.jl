 #=
 
 Step function for agent. Describes what every agent does every timestep
 The steps are:
 Moving: Choose a new location for the next timestep up to 1 grid cell away. Update the population grids after the move
 Transmission: Determine if the current healthy deer becomes infected at this timestep
 Progression: For infected deer, move them through the stage compartments described in the paper (e.g., E1 -> E2). 
   If the deer has reached the max length of their compartment length, move them to the next stage (e.g., I -> C)
 Shed prions: For infected deer, add to the prion reservoir 

 =#

function agent_step!(agent, model)
    old_pos = agent.pos
    was_infectious = agent.status in (:I, :C)
    move_deer!(agent, model)

    if was_infectious && agent.pos != old_pos
        model.infectious_grid[old_pos...] -= 1
        model.infectious_grid[agent.pos...] += 1
    end

    model.population_grid[old_pos...] -= 1
    model.population_grid[agent.pos...] += 1
    agent.age += 1
    transmit!(agent, model)
    progress_disease!(agent, model)
    shed_prions!(agent, model)
end

#= 

Step function for the global model. Key steps:
- Track counts of population at each stage
- Determine which agents will die (from disease or natural causes) and remove them
- Reproduction: determine how many new deer to add
- Record statistics 
=# 

function model_step!(model)
    model.tick += 1

    # Prion decay
    model.V .*= model.prion_survival

    # Survival check (background + disease mortality)
    to_remove = Int[]

    for agent in allagents(model)
        if !check_survival!(agent, model)
            push!(to_remove, agent.id)
        end
    end

    for id in to_remove
        dying_agent = model[id]
        if dying_agent.status in (:I, :C)
            model.infectious_grid[dying_agent.pos...] -= 1
        end
        model.population_grid[dying_agent.pos...] -= 1
        remove_agent!(id, model)
    end

    n_total = 0
    n_S = 0; n_E = 0; n_I = 0; n_C = 0
    
    # Sex specific counts
    m_total = 0; m_inf = 0
    f_total = 0; f_inf = 0
    
    # Age specific counts
    fawn_total = 0; fawn_inf = 0
    yrl_total  = 0; yrl_inf  = 0
    ad_total   = 0; ad_inf   = 0

    for agent in allagents(model)
        n_total += 1

        if agent.status == :S; n_S += 1
        elseif agent.status == :E; n_E += 1
        elseif agent.status == :I; n_I += 1
        elseif agent.status == :C; n_C += 1
        end

        is_inf = agent.status in (:E, :I, :C)

        # Sex Counts
        if agent.sex == :male
            m_total += 1
            if is_inf; m_inf += 1; end
        else
            f_total += 1
            if is_inf; f_inf += 1; end
        end

        # Age Counts (in weeks)
        if agent.age < 52
            fawn_total += 1
            if is_inf; fawn_inf += 1; end
        elseif agent.age < 104
            yrl_total += 1
            if is_inf; yrl_inf += 1; end
        else
            ad_total += 1
            if is_inf; ad_inf += 1; end
        end
    end

    push!(model.weekly_population, n_total)
    
    prev = n_total > 0 ? (n_E + n_I + n_C) / n_total : 0.0
    push!(model.weekly_prevalence, prev)
    
    push!(model.weekly_S, n_S)
    push!(model.weekly_E, n_E)
    push!(model.weekly_I, n_I)
    push!(model.weekly_C, n_C)

    # Demographic Histories
    push!(model.prev_male, m_total > 0 ? m_inf / m_total : 0.0)
    push!(model.prev_female, f_total > 0 ? f_inf / f_total : 0.0)
    
    push!(model.prev_fawn, fawn_total > 0 ? fawn_inf / fawn_total : 0.0)
    push!(model.prev_yearling, yrl_total > 0 ? yrl_inf / yrl_total : 0.0)
    push!(model.prev_adult, ad_total > 0 ? ad_inf / ad_total : 0.0)


    # Reproduction
    reproduction!(model)

    # Record weekly statistics
    push!(model.weekly_direct, model.total_direct - model.prev_total_direct)
    push!(model.weekly_indirect, model.total_indirect - model.prev_total_indirect)

    model.prev_total_direct = model.total_direct
    model.prev_total_indirect = model.total_indirect
    push!(model.weekly_V_total, sum(model.V))
end



function move_deer!(agent, model)
    randomwalk!(agent, model)
end

# Compute force of infection based on local neighbors
function compute_local_lambda(agent, model)
    nx, ny = size(model.V)
    pos = agent.pos
    x, y = pos
    infectious_grid = model.infectious_grid
    population_grid = model.population_grid
    r = model.transmission_radius


    # Count local infectious deer (Z) and local population (N)
    Z_local = 0
    N_local = 0
    
    for dy in -r:r
        for dx in -r:r
            cx, cy = x + dx, y + dy
            if 1 <= cx <= nx && 1 <= cy <= ny
                Z_local += infectious_grid[cx, cy]
                N_local += population_grid[cx, cy]
            end
        end
    end

    if agent.status in (:I, :C)
        Z_local = max(0, Z_local - 1)
    end

    V_local = model.V[pos...]
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
            model.infectious_grid[agent.pos...] += 1
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

    survival = 0.0
    
    if agent.age < 8
        survival = model.surv_fawn_early
    elseif agent.age < model.age_yearling
        survival = model.surv_fawn_late
    elseif agent.age < model.age_adult
        survival = model.surv_yearling
    else
        survival = (agent.sex == :male) ? model.surv_adult_male : model.surv_adult_female
    end

    if agent.status == :C
        survival /= model.clinical_mortality_mult
    end

    if agent.age >= model.max_age
        survival = 0.0
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

    density_factor = 1.0 / (1.0 + (N / K)^2)

    AgentT = NamedTuple{(:pos, :home_center), Tuple{Tuple{Int,Int}, Tuple{Int,Int}}}
    new_agents = AgentT[]

    for agent in allagents(model)
        if agent.sex == :female
            max_rate = 0.0
            if agent.age < model.age_yearling
                max_rate = model.birth_rate_fawn
            elseif agent.age < model.age_adult
                max_rate = model.birth_rate_yearling
            else
                max_rate = model.birth_rate_adult
            end
            
            prob = (max_rate * density_factor) / 52.0
            
            if rand(abmrng(model)) < prob
                offset = (rand(abmrng(model), -3:3), rand(abmrng(model), -3:3))
                new_home = (
                    clamp(agent.home_center[1] + offset[1], 1, nx),
                    clamp(agent.home_center[2] + offset[2], 1, ny)
                )
                push!(new_agents, (pos=new_home, home_center=new_home))
            end
        end
    end

    for na in new_agents
        new_sex = rand(abmrng(model)) < 0.5 ? :male : :female
        model.population_grid[na.pos...] += 1
        add_agent!(na.pos, DeerAgent, model;
            status=:S,
            weeks_in_state=0,
            state_duration=0,
            home_center=na.home_center,
            infection_source=:none,
            infection_tick=-1,
            sex=new_sex,
            age=0
        )
    end
end
