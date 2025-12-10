module PrionABM

using Agents, Random
using Agents.DataFrames, Agents.Graphs
using DrWatson: @dict
using CairoMakie
using Distributions: truncated, Normal, mean
using GLMakie

# Include the separated source files
include("agent_definition.jl")
include("model_creation.jl")
include("stepping_functions.jl")
include("analysis.jl")
include("plotting.jl")

# Export the core functions that a user might interact with
export model_initiation, agent_step!, model_step!
export compute_metrics, plot_dashboard, plot_transmission_analysis, plot_spatial_scenario
export DeerAgent, step!
export seic_color


end # module PrionABM
