"""
    R₄(Y)

Return R₄ metrics on Y.
"""
function R₄(Y, converged)
    ε = 1e-4

    tˢ, tᵉ= 1, 1
    t = 1
    Nₒ = Nₜ(-Inf)
    while true
        N = Nₜ(t)
        if N/Nₒ >= 1.01
            tˢ = t
            break
        end
        t += 1
    end

    t = tˢ + 1
    Nₒ = Nₜ(Inf)
    while true
        N = Nₜ(t)
        if Nₜ(t-1) > Nₜ(t) && N/Nₒ <= 1.01
            tᵉ = t - 1
            break
        end
        t += 1
    end
    
    tᵖ = argmin(Y[(tˢ-1):end]) + (tˢ - 2)
    tʳ = findfirst(x -> (Y[end] - ε <= x <= Y[end] + ε), Y[tᵖ:end]) + (tᵖ - 1)   
    
    rob = converged ? Y[tᵖ] : 0.0
    red = converged ? atan(((tᵖ - (tˢ - 1))/(tᵉ - tˢ + 1))/(Y[tˢ-1] - Y[tᵖ]))/(π/2) : NaN
    res = converged ? (Y[tʳ] - Y[tᵖ])/(1 - Y[tᵖ]) : 0.0
    rap = converged ? atan((Y[tʳ] - Y[tᵖ])/((tʳ - tᵖ)/(tᵉ - tˢ + 1)))/(π/2) : NaN

    return round.((rob, red, res, rap), digits=3)
end

    
"""
    analyze(sim, saveplot=true)

Analyzes simulated data, creating response plots and rendering R₄ metric.
"""
function analyze(sim, saveplot=true)
    function abbreviate(x)
        if length(x) <= 3 return "$x."
        else
            i = 3
            while true
                if i > length(x) return "$(x[1] * x[2] * x[3])." end 
                if x[i] ∉ ['a', 'e', 'i', 'o', 'u'] return "$(x[1] * x[2] * x[i])." end
                i += 1
            end
        end
    end

    tˢ, tᵉ= 1, 1
    t = 1
    Nₒ = Nₜ(-Inf)
    while true
        N = Nₜ(t)
        if N/Nₒ >= 1.01
            tˢ = t
            break
        end
        t += 1
    end

    t = tˢ + 1
    Nₒ = Nₜ(Inf)
    while true
        N = Nₜ(t)
        if Nₜ(t-1) > Nₜ(t) && N/Nₒ <= 1.01
            tᵉ = t - 1
            break
        end
        t += 1
    end

    ε = 1e-4

    data = load(sim)
    N̄, Ñ, TC⃰, re_optimize = data["N̄"], data["Ñ"], data["TC⃰"], data["re_optimize"]
    Label = [abbreviate(String(re_opt)) for re_opt in re_optimize]

    R = size(Ñ)[1]
    C = size(Ñ)[2]
    df = DataFrame(Metric          = fill("", 2R), 
                    ReOpt           = fill("", 2R), 
                    Robustness      = zeros(2R), 
                    Redundancy      = zeros(2R), 
                    Resourcefulness = zeros(2R), 
                    Rapidity        = zeros(2R))
    
    # Demand
    fig = plot(xlab="day", ylab="demand", title="demand", titlefontsize=11)
    plot!([1:0.1:(tᵉ+14)], [Nₜ(t) for t ∈ 1:0.1:(tᵉ+14)], label="")
    vline!([tˢ], linestyle=:dash, color=:black, label="")
    vline!([tᵉ], linestyle=:dash, color=:black, label="")
    display(fig)
    if saveplot savefig(fig, joinpath(@__DIR__, "Results\\demand.png")) end
    
    # Congestion
    fig = plot(xlab="day", ylab="congestion factor", title="congestion factor", titlefontsize=11)
    plot!([1:0.1:(tᵉ+14)], [φₜ(t) for t ∈ 1:0.1:(tᵉ+14)], label="")
    vline!([tˢ], linestyle=:dash, color=:black, label="")
    vline!([tᵉ], linestyle=:dash, color=:black, label="")
    display(fig)
    if saveplot savefig(fig, joinpath(@__DIR__, "Results\\congestion_factor.png")) end
    
    # Driver availability
    fig = plot(xlab="day", ylab="driver availability %", title="driver availability", titlefontsize=11)
    plot!([1:0.1:(tᵉ+14)], [f̅ₜ(100, t) for t ∈ 1:0.1:(tᵉ+14)], label="")
    vline!([tˢ], linestyle=:dash, color=:black, label="")
    vline!([tᵉ], linestyle=:dash, color=:black, label="")
    display(fig)
    if saveplot savefig(fig, joinpath(@__DIR__, "Results\\driver_availabilty.png")) end
    
    # Level of Service
    Z = zeros(size(Ñ))
    for r ∈ 1:R
        for c ∈ 1:C
            c < tˢ ? t = -Inf : t = c 
            Ñ[r, c] == 0 ? Z[r, c] = 1.0 : Z[r, c] = N̄[r, c]/(Ñ[r, c] + N̄[r, c])  
        end
    end
    println()

    for r ∈ 1:R
        ix = C - findfirst(x -> !isnan(x), reverse(TC⃰[r, :])) + 1
        bool = N̄[r, ix] >= Nₜ(Inf)
        Y = filter(x -> !isnan(x), Z[r, :])
        rob, red, res, rap     = R₄(Y, bool) 
        df[r,:Metric]          = "LOS"
        df[r,:ReOpt]           = Label[r]
        df[r,:Robustness]      = rob
        df[r,:Redundancy]      = red
        df[r,:Resourcefulness] = res
        df[r,:Rapidity]        = rap
    end
    y = transpose(Z)
    x = [t for t ∈ 1:size(y)[1], _ ∈ 1:size(y)[2]]
    fig = plot(xlab="day", ylab="LOS", title="Level of Service (LOS)", titlefontsize=11)
    plot!(x, y, label="")#, label=hcat(Label), legend=:bottomleft)
    vline!([tˢ], linestyle=:dash, color=:black, label="")
    vline!([tᵉ], linestyle=:dash, color=:black, label="")
    display(fig)
    if saveplot savefig(fig, joinpath(@__DIR__, "Results\\LOS.png")) end
    
    # Level of efficiency
    Z = zeros(size(Ñ))
    for r ∈ 1:R 
        for c ∈ 1:C 
            Z[r, c] = TC⃰[r,1]/TC⃰[r,c] 
        end 
    end
    
    for r ∈ 1:R
        ix = C - findfirst(x -> !isnan(x), reverse(TC⃰[r, :])) + 1
        bool = bool = N̄[r, ix] >= Nₜ(Inf)
        Y = filter(x -> !isnan(x), Z[r, :])
        rob, red, res, rap         = R₄(Y, bool) 
        df[r + R,:Metric]          = "LOE"
        df[r + R,:ReOpt]           = Label[r]
        df[r + R,:Robustness]      = rob
        df[r + R,:Redundancy]      = red
        df[r + R,:Resourcefulness] = res
        df[r + R,:Rapidity]        = rap
    end
    y = transpose(Z)
    x = [t for t ∈ 1:size(y)[1], _ ∈ 1:size(y)[2]]
    fig = plot(xlab="day", ylab="LOE", title="Level of Efficiency (LOE)", titlefontsize=11)
    plot!(x, y, label="")#, label=hcat(Label), legend=:bottomleft)
    vline!([tˢ], linestyle=:dash, color=:black, label="")
    vline!([tᵉ], linestyle=:dash, color=:black, label="")
    display(fig)
    if saveplot savefig(fig, joinpath(@__DIR__, "Results\\LOE.png")) end

    display(df)
end