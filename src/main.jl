using JuMP
using Juniper
using Ipopt
using Statistics
using Plots
using JLD
using DataFrames
using CSV
using CSVFiles

double_logistic(x, α₁, α₂, μ₁, μ₂, θ₁, θ₂) = α₁/(1 + exp(-(x - μ₁)/θ₁)) - α₂/(1 + exp(-(x - μ₂)/θ₂))
Nₜ(t) = (1 + double_logistic(t, 0.685, 0.486, 49.373, 89.512, 8.447, 7.885)) * 30000
φₜ(t) = (1 + double_logistic(t, 0.1274, 0.1274, 49.373, 89.512, 8.447, 7.885)) * 0.887
v, d = findmax([Nₜ(t) for t in 1:500])
include("optimize.jl")       # Optimizer
include("simulate.jl")       # Simulator
include("analyze.jl")        # Analyzer


function experiment()
    outsource_channel = :x
    Nᶠ = 0
    pᵤ = 0.0
    f̅′ = 0
    re_optimize = :operational
 
    parameters = (Nᶠ = Nᶠ, pᵤ = pᵤ, f̅′ = f̅′, re_optimize = re_optimize)
    sim_name = "sim_data - $outsource_channel"
    simulate(outsource_channel, parameters, sim_name) 
    sim = joinpath(@__DIR__, "Results\\$sim_name.jld")
    analyze(sim, true)
end


function batch_analyze(dir)
    df = DataFrame(file = String[], start = Int64[], recovery = Int64[], 
                   rob = Float64[], red = Float64[], res = Float64[], rap = Float64[],
                   TD = Float64[], APD₁ = Float64[], APD₂ = Float64[], APD₃ = Float64[],
                   DL = Float64[], IL = Float64[], TL = Float64[],
                   DLp = Float64[], ILp = Float64[], TLp = Float64[])
    for jld in readdir(dir)
        println("# ─── $jld ────────────────────────────────────────────────────────────────────────")
        sim = joinpath(dir, jld)
        tˢ, tʳ, rob, red, res, rap, TD, APD₁, APD₂, APD₃, DL, IL, TL, DLp, ILp, TLp = analyze(sim, false)
        push!(df[!, :file], jld)
        push!(df[!, :start], tˢ)
        push!(df[!, :recovery], tʳ)
        push!(df[!, :rob], rob)
        push!(df[!, :red], red)
        push!(df[!, :res], res)
        push!(df[!, :rap], rap)
        push!(df[!, :TD], TD)
        push!(df[!, :APD₁], APD₁)
        push!(df[!, :APD₂], APD₂)
        push!(df[!, :APD₃], APD₃)
        push!(df[!, :DL], DL)
        push!(df[!, :IL], IL)
        push!(df[!, :TL], TL)
        push!(df[!, :DLp], DLp)
        push!(df[!, :ILp], ILp)
        push!(df[!, :TLp], TLp)
    end
    CSV.write("df - sensitivity analysis - θ₂.csv", df)
end


function figures()
    # figure 1. distribution capacity in base case  ────────────────────────────────────────────────────────────────────────────────
    tˢ = 14
    tᵉ = 118
    x  = load(joinpath(@__DIR__, "Results/#1. Base Case/#1. No outsourcing/sim_data - x - operational.jld"))
    cs = load(joinpath(@__DIR__, "Results/#1. Base Case/#2. Crowdsourcing/sim_data - cs - operational.jld"))
    cp = load(joinpath(@__DIR__, "Results/#1. Base Case/#3. Collection-points/sim_data - cp - operational.jld"))
    mh = load(joinpath(@__DIR__, "Results/#1. Base Case/#4. Micro-hubs/sim_data - mh - operational.jld"))
    fig = plot(xlab="day", ylab="distribution capacity (customers served/day)",
               xlims=(0, tᵉ+14), fontfamily="times", 
               legend=(0.625, 0.25), fg_legend=:transparent)
    plot!(x["N̄"] , label="w/o outsourcing")
    plot!(cs["N̄"], label="w/ crowdsourced fleet")
    plot!(cp["N̄"], label="w/ collection-points")
    plot!(mh["N̄"], label="w/ LSP's micro-hubs")
    vline!([tˢ], linestyle=:dash, color=:black, label="")
    vline!([tᵉ], linestyle=:dash, color=:black, label="")
    annotate!(tˢ-2, 0.77 * 5e4, text("Disruption start day", 8, family="times", rotation=90))
    annotate!(tᵉ+2, 0.77 * 5e4, text("Disruption end day"  , 8, family="times", rotation=90))
    display(fig)

    df = DataFrame(x = x["N̄"][1:132], cs = cs["N̄"][1:132], cp = cp["N̄"][1:132], mh = mh["N̄"][1:132])
    CSV.write("df - distribution capacity.csv", df)
    
    ## figure 2. distribution costs in base case  ────────────────────────────────────────────────────────────────────────────────
    tˢ = 14
    tᵉ = 118
    x  = load(joinpath(@__DIR__, "Results/#1. Base Case/#1. No outsourcing/sim_data - x - operational.jld"))
    cs = load(joinpath(@__DIR__, "Results/#1. Base Case/#2. Crowdsourcing/sim_data - cs - operational.jld"))
    cp = load(joinpath(@__DIR__, "Results/#1. Base Case/#3. Collection-points/sim_data - cp - operational.jld"))
    mh = load(joinpath(@__DIR__, "Results/#1. Base Case/#4. Micro-hubs/sim_data - mh - operational.jld"))
    fig = plot(xlab="day", ylab="distribution cost per package (\$)",
               xlims=(0, tᵉ+14), fontfamily="times",
               legend=(0.625, 0.9), fg_legend=:transparent)
    plot!(x["TC⃰"] .* x["N"]  , label="w/o outsourcing")
    plot!(cs["TC⃰"] .* cs["N"], label="w/ crowdsourced fleet")
    plot!(cp["TC⃰"] .* cp["N"], label="w/ collection-points")
    plot!(mh["TC⃰"] .* mh["N"], label="w/ LSP's micro-hubs")
    vline!([tˢ], linestyle=:dash, color=:black, label="")
    vline!([tᵉ], linestyle=:dash, color=:black, label="")
    annotate!(tˢ-2, 0.77 * 1.4e5, text("Disruption start day", 8, family="times", rotation=90))
    annotate!(tᵉ+2, 0.77 * 1.4e5, text("Disruption end day"  , 8, family="times", rotation=90))
    display(fig)

    df = DataFrame(x = (x["TC⃰"] .* x["N"])[1:132], cs = (cs["TC⃰"] .* cs["N"])[1:132], cp = (cp["TC⃰"] .* cp["N"])[1:132], mh = (mh["TC⃰"] .* mh["N"])[1:132])
    CSV.write("df - distribution cost.csv", df)
end

#experiment()

directory = joinpath(@__DIR__, "Results/#3. Sensitivity Analysis/#6. θ₂/")
batch_analyze(directory)

#figures()

LMO(:minTC, (N=Nₜ(-Inf), φ=φₜ(-Inf)), false)
LMO(:maxN, (ρₓ=5.44, φ=φₜ(d), f̅=152, Ψᵐʰ=1, Nᵐʰ=10, f̅′=220, pᵤ = 1.0), false)

#LMO(:minTC, (N=Nₜ(d), φ=φₜ(d), ρₓ=6.447, f̅=98, Ψᶜˢ=1, f̅′=565, pᵤ=1.0), false)
#LMO(:minTC, (N=Nₜ(d), φ=φₜ(d), ρₓ=6.447, f̅=98, Ψᶜᵖ=1, Nᶜᵖ=200, f̅′=v, pᵤ=0.85), false)
#LMO(:minTC, (N=Nₜ(d), φ=φₜ(d), ρₓ=6.447, f̅=98, Ψᵐʰ=1, Nᵐʰ=10, f̅′=220, pᵤ = 1.0), false)

#LMO(:minTC, (N=3.6e4, φ=0.887, ρₓ=6.447, f̅=98, Ψᶜˢ=1, f̅′=565, pᵤ=1.0), false)
#LMO(:minTC, (N=3.6e4, φ=0.887, ρₓ=6.447, f̅=98, Ψᶜᵖ=1, Nᶜᵖ=200, f̅′=v, pᵤ=0.85), false)
#LMO(:minTC, (N=3.6e4, φ=0.887, ρₓ=6.447, f̅=100, Ψᶜᵖ=1, Nᶜᵖ=200, f̅′=v, pᵤ=0.27), false)
#LMO(:minTC, (N=3.6e4, φ=0.887, ρₓ=6.447, f̅=117, Ψᶜᵖ=1, Nᶜᵖ=200, f̅′=v, pᵤ=0.0), false)
#LMO(:minTC, (N=3.6e4, φ=0.887, ρₓ=6.447, f̅=98, Ψᵐʰ=1, Nᵐʰ=10, f̅′=220, pᵤ = 1.0), false)
#LMO(:minTC, (N=3.6e4, φ=0.887, ρₓ=6.447, f̅=100, Ψᵐʰ=1, Nᵐʰ=10, f̅′=220, pᵤ = 1.0), false)
#LMO(:minTC, (N=3.6e4, φ=0.887, ρₓ=6.447, f̅=117, Ψᵐʰ=1, Nᵐʰ=10, f̅′=220, pᵤ = 1.0), false)

#= ────────────────────────────────────────────────────────────────────────────────
# NOTE:
Double Logistic Model
        yₜ = Qₜ/Qₒ - 1
        yₜ = α₁/(1 + exp(-(t - μ₁)/θ₁)) - α₂/(1 + exp(-(t - μ₂)/θ₂))

   Double Logistic Model Non-Linear Least-Squares Regression:
        min z = Σ ε²
            z = Σₜ (α₁/(1 + exp(-(t - μ₁)/θ₁)) - α₂/(1 + exp(-(t - μ₂)/θ₂)) - yₜ)²

   Interpretation of parameters:
       1. α₁: Growth factor
       2. α₂: Decay factor
       3. μ₁: Half-life of growth
       4. μ₂: Half-life of decay
       5. θ₁: Inverse of average growth rate
       6. θ₂: Inverse of average decay rateō
──────────────────────────────────────────────────────────────────────────────── =#
#= ────────────────────────────────────────────────────────────────────────────────
# NOTE:
    Scenarios
    - X: No outsourcing

    - CS: Crowd-sourcing
        Nᶠ: 0
        pᵤ: 1.0 (0.352)
        f̅′: 565 delivery vehicles
            
    - CP:
        Nᶠ: 200
        pᵤ: 0.85
        f̅′: max(Nₜ)

    - MH:
        Nᶠ: 10
        pᵤ: 1.0 (0.825)
        f̅′: 220 cargo bikes
──────────────────────────────────────────────────────────────────────────────── =#