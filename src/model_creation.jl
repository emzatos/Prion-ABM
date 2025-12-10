@kwdef mutable struct ModelConstants
    exposure_min::Int = 15
    exposure_max::Int = 35
    infectious_min::Int = 25
    infectious_max::Int = 44
    clinical_min::Int = 1
    clinical_max::Int = 36
    prion_survival::Float64 = 0.9967
    shedding_rate::Float64 = 1.0
    carcass_load::Float64 = 100.0
    bd::Float64 = 5.2e-6
    bi::Float64 = 5.5e-8
    k::Float64 = 0.1
    ε::Float64 = 0.0001
    transmission_radius::Int = 2
    weekly_survival::Float64 = 0.991
    clinical_mortality_mult::Float64 = 2.0
    annual_birth_rate::Float64 = 0.55
    offspring_per_birth::Float64 = 1.7
    carrying_capacity::Int = 20000
    tick::Int = 0
    total_direct::Int = 0
    total_indirect::Int = 0
    prev_total_direct::Int = 0
    prev_total_indirect::Int = 0
    
    V::Matrix{Float64} = zeros(100, 100)
    V_cumulative::Matrix{Float64} = zeros(100, 100)
    
    events::Vector{NamedTuple{(:tick, :source, :agent_id, :pos), Tuple{Int, Symbol, Int, Tuple{Int, Int}}}} = []
    
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

    return model
end
