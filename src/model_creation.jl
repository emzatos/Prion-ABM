# Parameters the model will keep track of during the simulation with default values set
@kwdef mutable struct ModelConstants
    ## SEIC Stage durations
    # Min/Max durations for each of the SEIC stages. We sample a normal distribution between the min/max values to find how long each agent should be in a stage
    exposure_min::Int = 15
    exposure_max::Int = 35
    infectious_min::Int = 25
    infectious_max::Int = 44
    clinical_min::Int = 1
    clinical_max::Int = 36

    # Prion related information
    prion_survival::Float64 = 0.9967    # Prion load decay rate each time step, equal to a half-life of 48 months
    shedding_rate::Float64 = 1.0        # How many prions dropped by infected deer per timestep
    carcass_load::Float64 = 100.0       # How many prions dropped at death

    # Transmission related information, used when calculating force of infection
    bd::Float64 = 9.2e-4               # direct transmission coefficient
    bi::Float64 = 5.5e-5               # indirect transmission coefficient
    k::Float64 = 0.1                    # aggregation factor from paper
    ε::Float64 = 0.0001                 # frequency/density dependence
    transmission_radius::Int = 2        # determines how far away an infected deer can be while still causing direct infections (not used for indirect)

    # Demographic information
    weekly_survival::Float64 = 0.991        # non-disease related death rate
    clinical_mortality_mult::Float64 = 2.0  # death rate multiplier for sick deer (they have a higher chance of spontaneous death)
    annual_birth_rate::Float64 = 0.55       # used to determine probability a deer has offspring at a timestep
    offspring_per_birth::Float64 = 1.7      # used to calculate expected value of new deer per timestep
    carrying_capacity::Int = 20000          # carrying capacity of space, used in reproduction probability 

    # Statistical measures
    tick::Int = 0                   # Keep track of timestep when recording events
    total_direct::Int = 0
    total_indirect::Int = 0
    prev_total_direct::Int = 0
    prev_total_indirect::Int = 0
    
    V::Matrix{Float64} = zeros(100, 100)                # Matrix storing current prion load 
    V_cumulative::Matrix{Float64} = zeros(100, 100)     # Matrix storing cumulative prion load 
    infectious_grid::Matrix{Int} = zeros(Int, 100, 100) # Matrix storing number of infected or clinical deer in each cell
    population_grid::Matrix{Int} = zeros(Int, 100,100)  # Matrix storing number of total deer in a cell
    
    events::Vector{NamedTuple{(:tick, :source, :agent_id, :pos), Tuple{Int, Symbol, Int, Tuple{Int, Int}}}} = [] # Used to track when and how infections occurred
    
    weekly_direct::Vector{Int} = []
    weekly_indirect::Vector{Int} = []
    weekly_prevalence::Vector{Float64} = []
    weekly_population::Vector{Int} = []
    weekly_V_total::Vector{Float64} = []
    weekly_S::Vector{Int} = []
    weekly_E::Vector{Int} = []
    weekly_I::Vector{Int} = []
    weekly_C::Vector{Int} = []
end

function model_initiation(;
    nx=100,
    ny=100,
    n_deer=5000,
    n_infected=10,
    seed=0,
    kw_args...
)
    rng = Xoshiro(seed)

    properties = ModelConstants(;
        V = zeros(Float64, nx, ny),
        V_cumulative = zeros(Float64, nx, ny),
        infectious_grid = zeros(Int, nx, ny),
        population_grid = zeros(Int, nx, ny),
        kw_args...
    )

    space = GridSpace((nx, ny); periodic=false)
    model = StandardABM(DeerAgent, space; agent_step!, model_step!, properties, rng)

    # Initialize deer population
    for i in 1:n_deer
        home = (rand(abmrng(model), 1:nx), rand(abmrng(model), 1:ny))
        if i <= n_infected
            status = :I
            duration = sample_duration(abmrng(model), properties.infectious_min, properties.infectious_max)
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
        )
    end

    for agent in allagents(model)
        if agent.status in (:I, :C)
            model.infectious_grid[agent.pos...] += 1
        end
        model.population_grid[agent.pos...] += 1
    end


    return model
end
