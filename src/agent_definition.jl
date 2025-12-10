# Definition of agent and its parameters
@agent struct DeerAgent(GridAgent{2})
    status::Symbol              # :S, :E, :I, :C
    weeks_in_state::Int         # weeks since entering current state
    state_duration::Int         # sampled duration for current state
    home_center::Tuple{Int,Int} # home range center for movement fidelity
    infection_source::Symbol    # :none, :direct, :indirect, :initial
    infection_tick::Int         # when infected (-1 if never)
    sex::Symbol                 # M or F
    age::Int
end

# Function used to determine how long an agent will stay in each of the SEIC stages
function sample_duration(rng, min_weeks::Int, max_weeks::Int)
    μ = (min_weeks + max_weeks) / 2
    σ = (max_weeks - min_weeks) / 4
    d = truncated(Normal(μ, σ), min_weeks, max_weeks)
    return round(Int, rand(rng, d))
end
