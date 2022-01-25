"""
    analyze(sim, saveplot=true)

Analyzes simulated data, creating response plots and rendering R₄ metric.
"""
function analyze(sim, showplots=true)
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

    ε = 1e-4

    tˢ= 1                                                                                                       # Disruption start day
    t = 2
    Nₒ = Nₜ(-Inf)
    while true
        N = Nₜ(t)
        if N/Nₒ >= 1.01
            tˢ = t
            break
        end
        t += 1
    end
    
    tᵉ= 1                                                                                                       # Disruption end day
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
    

    data = load(sim)
    N̄, Ñ, N, TC⃰, re_optimize = data["N̄"], data["Ñ"], data["N"], data["TC⃰"], data["re_optimize"]
    Label = abbreviate(String(re_optimize))
    T = length(Ñ)

    N = zeros(T)
    for t ∈ 1:T N[t] = (t ≤ tˢ ? Nₜ(-Inf) : (Ñ[t] == 0 ? Nₜ(t) + Ñ[t-1] : N̄[t])) end

    # ─── LEVEL OF SERVICE ───────────────────────────────────────────────────────────
    Z = zeros(T)
    for t ∈ 1:T Z[t] = Ñ[t] == 0 ? 1.0 : N̄[t]/(N̄[t] + Ñ[t]) end
    LOS = copy(Z)
    converged = Ñ[T] == 0
    Y  = copy(Z)
    tᵖ = argmin(Y[(tˢ-1):end]) + (tˢ - 2)                                                                       # Peak peri-disruption day
    tʳ = findfirst(x -> (Y[end] - ε <= x <= Y[end] + ε), Y[tᵖ:end]) + (tᵖ - 1)                                  # Recovery day
    rob = converged ? round(Y[tᵖ], digits=3) : 0.0                                                              # Robustness
    red = converged ? round(atan(((tᵖ - (tˢ - 1))/(tᵉ - tˢ + 1))/(Y[tˢ-1] - Y[tᵖ]))/(π/2), digits=3) : NaN      # Redundancy
    res = converged ? round((Y[tʳ] - Y[tᵖ])/(1 - Y[tᵖ]), digits=3) : 0.0                                        # Resourcefulness
    rap = converged ? round(atan((Y[tʳ] - Y[tᵖ])/((tʳ - tᵖ)/(tᵉ - tˢ + 1)))/(π/2), digits=3) : NaN              # Rapidity
    TD   = sum(Ñ[tˢ:tʳ])
    APD₁ = sum(Ñ[tˢ:tʳ])/(tʳ - tˢ)
    APD₂ = sum(Ñ[tˢ:tʳ])/(sum(Nₜ.(tˢ:tʳ)))
    APD₃ = tʳ > tˢ ? sum([Nₜ(t) - N[t] ≥ 0.0 ? Nₜ(t) - N[t] : 0.0 for t ∈ tˢ:tʳ])/(tʳ - tˢ) : NaN
    println("Level of Service - $Label")
    println("   Disruption start day   : $tˢ")
    println("   Peak day               : $tᵖ")
    println("   Recovery day           : $tʳ")
    println("   Time to recovery       : $(tʳ - tˢ)")
    println()
    println("   Robustness             = $rob")
    println("   Redundancy             = $red")
    println("   Resourcefulness        = $res")
    println("   Rapidity               = $rap")
    println()
    println("   Total delay (pkg-days) = $TD")
    println("   Avg. pkg delay (pkg)   = $APD₁")
    println("   Avg. pkg delay (days)  = $APD₂")
    println("   Avg. pkgs delayed      = $APD₃")

    # ─── LOSS ───────────────────────────────────────────────────────────────────────
    Z = zeros(T)
    for t ∈ 1:T  Z[t] = t ≤ tˢ ? TC⃰[1] * Nₜ(-Inf) : TC⃰[t] * N[t] end
    C = copy(Z)
    Y = copy(Z) .- TC⃰[1] * Nₜ(-Inf)
    tᵖ = argmax(Y[(tˢ-1):end]) + (tˢ - 2)
    tʳ = findfirst(x -> (Y[end] - ε <= x <= Y[end] + ε), Y[tᵖ:end]) + (tᵖ - 1)   
    DL = round(sum(Y[1:tʳ]), digits=1)                                                                          # Direct loss
    IL = round(5 * sum(Ñ[1:tʳ]), digits=3)                                                                      # Indirect loss
    TL = round((sum(Y[1:tʳ]) + 5 * sum(Ñ[1:tʳ])), digits=3)                                                     # Total loss    
    println("\n\nNet Loss")
    println("   Disruption start day   : $tˢ")
    println("   Peak day               : $tᵖ")
    println("   Recovery day           : $tʳ")
    println("   Time to recovery       : $(tʳ - tˢ)")
    println()
    println("   Direct Loss (\$)      = $DL")
    println("   Indirect Loss (\$)    = $IL")
    println("   Total Loss (\$)       = $TL")
    println("\n\n")

    if showplots
        # ─── DEMAND ──────────────────────────────────────────────────────
        fig = plot(yrange=(0,1), xlab="day", ylab="demand", title="demand", titlefontsize=11)
        plot!([1:0.1:(tᵉ+14)], [Nₜ(t)/Nₜ(-Inf) - 1 for t ∈ 1:0.1:(tᵉ+14)], label="")
        vline!([tˢ], linestyle=:dash, color=:black, label="")
        vline!([tᵉ], linestyle=:dash, color=:black, label="")
        display(fig)
    
        # ─── CONGESTION ──────────────────────────────────────────────────
        fig = plot(xlab="day", ylab="congestion factor", title="congestion factor", titlefontsize=11)
        plot!([1:0.1:(tᵉ+14)], [φₜ(t) for t ∈ 1:0.1:(tᵉ+14)], label="")
        vline!([tˢ], linestyle=:dash, color=:black, label="")
        vline!([tᵉ], linestyle=:dash, color=:black, label="")
        display(fig)

        # ─── DISTRIBUTION CAPACITY ──────────────────────────────────────
        y = N̄
        x = 1:T
        fig = plot(xlab="day", ylab="Distribution capacity", title="distribution capacity", titlefontsize=11)
        plot!(x, y, label="")
        vline!([tˢ], linestyle=:dash, color=:black, label="")
        vline!([tᵉ], linestyle=:dash, color=:black, label="")
        display(fig)

        # ─── LEVEL OF SERVICE ───────────────────────────────────────────────────────────
        y = LOS
        x = 1:T
        fig = plot(yrange=(0,1), xlab="day", ylab="LOS", title="Level of Service (LOS)", titlefontsize=11)
        plot!(x, y, label="")
        vline!([tˢ], linestyle=:dash, color=:black, label="")
        vline!([tᵉ], linestyle=:dash, color=:black, label="")
        display(fig)

        # ─── DISTRIBUTION COST ───────────────────────────────────────────────────────────
        y = C
        x = 1:T
        fig = plot(xlab="day", ylab="TC (\$)", title="Total Cost (\$)", titlefontsize=11)
        plot!(x, y, label="")
        vline!([tˢ], linestyle=:dash, color=:black, label="")
        vline!([tᵉ], linestyle=:dash, color=:black, label="")
        display(fig)
    end

    return tˢ, tʳ, rob, red, res, rap, TD, APD₁, APD₂, APD₃, DL, IL, TL
end