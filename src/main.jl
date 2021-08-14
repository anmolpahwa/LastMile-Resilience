using JuMP
using Juniper
using Ipopt
using Statistics
using Plots
using JLD
using DataFrames


# Modeling disruption
double_logistic(x, α₁, α₂, μ₁, μ₂, θ₁, θ₂) = α₁/(1 + exp(-(x - μ₁)/θ₁)) - α₂/(1 + exp(-(x - μ₂)/θ₂))
Nₜ(t) = (1 + double_logistic(t, 0.685, 0.486, 49.373, 89.512, 8.447, 7.885)) * 30000
φₜ(t) = (1 + double_logistic(t, 0.1274, 0.1274, 49.373, 89.512, 8.447, 7.885)) * 0.887
f̅ₜ(f, t) = Int(floor((1 - double_logistic(t, 0.0, 0.0, 49.373, 89.512, 8.447, 7.885)) * f))
    

include("optimize.jl")       # Optimizer
include("simulate.jl")       # Simulator
include("analyze.jl")        # Analyzer


"""
    experiment()

Workflow for running last-mile resilience analysis.
"""
function experiment()
    outsource_channel = :x
    parameters = (Nᶠ = 0, pᵤ = 0.0, f̅′ = 0, re_optimize = [:operational])
    sim = joinpath(@__DIR__, "Results\\sim_data.jld")
    simulate(outsource_channel, parameters)
    analyze(sim, true)
end
experiment()