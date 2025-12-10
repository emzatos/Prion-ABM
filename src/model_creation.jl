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
    surv_adult_male::Float64 = 0.9947
    surv_adult_female::Float64 = 0.9985
    surv_yearling::Float64 = 0.9937
    surv_fawn_early::Float64 = 0.9920     # First 8 weeks
    surv_fawn_late::Float64 = 0.9984     # Weeks 9-52

    birth_rate_adult::Float64 = 1.80
    birth_rate_yearling::Float64 = 1.25
    birth_rate_fawn::Float64 = 0.40

    # Age Thresholds (in weeks)
    age_yearling::Int = 52  # 1 year
    age_adult::Int = 104 # 2 years
    max_age::Int = 624 # 12 years (~12 years max age)
    
    clinical_mortality_mult::Float64 = 2.0  # death rate multiplier for sick deer (they have a higher chance of spontaneous death)
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
    
    prev_male::Vector{Float64} = []
    prev_female::Vector{Float64} = []

    prev_fawn::Vector{Float64} = []     # 0-1 yr
    prev_yearling::Vector{Float64} = [] # 1-2 yr
    prev_adult::Vector{Float64} = []    # 2+ yr

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
        r = rand(rng)
        
        # Default values
        sex = :female
        age = 0
        
        if r < 0.12 # Adult Male (12%)
            sex = :male
            age = rand(rng, 104:624)
        elseif r < 0.12 + 0.33 # Adult Female (33%)
            sex = :female
            age = rand(rng, 104:624)
        elseif r < 0.45 + 0.09 # Yearling Male (9%)
            sex = :male
            age = rand(rng, 52:103)
        elseif r < 0.54 + 0.10 # Yearling Female (10%)
            sex = :female
            age = rand(rng, 52:103)
        else # Fawn (36%) - Split 50/50 M/F
            sex = rand(rng) < 0.5 ? :male : :female
            age = rand(rng, 0:51)
        end

        home = (rand(abmrng(model), 1:nx), rand(abmrng(model), 1:ny))
        if i <= n_infected
            status = :I
            duration = sample_duration(abmrng(model), properties.infectious_min, properties.infectious_max)
            source = :initial
            inf_tick = 0
            model.infectious_grid[home...] += 1
        else
            status = :S
            duration = 0
            source = :none
            inf_tick = -1
        end
        model.population_grid[home...] += 1

        add_agent!(home, DeerAgent, model;
            status,
            weeks_in_state=0,
            state_duration=duration,
            home_center=home,
            infection_source=source,
            infection_tick=inf_tick,
            sex=sex,
            age=age
        )
    end
    return model
end
