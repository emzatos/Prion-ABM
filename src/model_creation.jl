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
        )
    end

    return model
end
