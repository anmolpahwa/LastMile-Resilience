"""
    simulate(outsource_channel::Symbol, parameters::NamedTuple)

Simulates last-mile response and stores simulation as a jld file.

# Arguments
- outsource_channel: One of `:cs`, `:mh` or `:cp`.
- parameters: A named tuple of simulation parameters.
    - `pᵤ` (Max outsource demand share)                 : Required.
    - `Nᶠ` (Number of facilities)                       : Required.
    - `f̅′` (Available outsource fleet)                  : Required.
    - `re_optimize` (Post-distruption re-optimization)  : Required. Subset of `[:operational, :tactical, :strategic]`.
"""
function simulate(outsource_channel::Symbol, parameters::NamedTuple)    
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
    
    ε = 1e-4             # Cut-off tolerance
    Tᴸ = tᵉ + 7 * 2      # Lower Cut-off time
    Tᴴ = Tᴸ              # Higher Cut-off time

    t = tᵉ
    Nₒ = Nₜ(Inf)
    while true
        N = Nₜ(t)
        if N/Nₒ <= 1 + ε
            Tᴴ = t
            break
        end
        t += 1
    end

    Nᶠ = parameters[:Nᶠ]
    pᵤ = parameters[:pᵤ]
    f̅′ = parameters[:f̅′]
    
    l   = length(parameters[:re_optimize])
    N̄   = zeros((l, tᵉ))
    Ñ   = zeros((l, tᵉ))
    TC⃰ = zeros((l, tᵉ)) 
    Ψᶜˢ = 0 
    Ψᵐʰ = 0
    Ψᶜᵖ = 0
    Nᵐʰ = 0
    Nᶜᵖ = 0
    if outsource_channel == :cs Ψᶜˢ = 1 end
    if outsource_channel == :mh Ψᵐʰ = 1 end
    if outsource_channel == :cp Ψᶜᵖ = 1 end
    if outsource_channel == :mh Nᵐʰ = Nᶠ end
    if outsource_channel == :cp Nᶜᵖ = Nᶠ end
    println()

    
    # ─────────────────────────────────────────────────────────────────
    ## Pre-disruption phase
    
    println("\n-----PRE - DISRUPTION-----")
    TS, _, ρₓ⃰, f₁⃰, _, TC₁⃰ = LMO(:minTC, 
                                  (N = Nₜ(-Inf), φ = φₜ(-Inf)), 
                                   true)
    f⃰ = Int(f₁⃰)
    println("   Minimizing total cost per package: $TS")
    println("       Total cost per package: $(round(TC₁⃰, digits=2))\n")
    
    TS, N̄₁, _, _, _, T̄C̄₁   = LMO(:maxN, 
                                (ρₓ = ρₓ⃰, φ = φₜ(-Inf), 
                                 fₒ = f⃰, f̅ = f̅ₜ(f⃰, -Inf), f̅′ = 0), 
                                 true)
    println("   Maximizing distribution capacity without outsourcing: $TS")
    println("       Distribution capacity: $(round(N̄₁, digits=2))\n")
    
    TS′, N̄₁′, _, _, _, T̄C̄₁′ = LMO(:maxN, 
                                 (ρₓ = ρₓ⃰, φ = φₜ(-Inf),
                                  Ψᶜˢ = Ψᶜˢ, Ψᵐʰ = Ψᵐʰ, Ψᶜᵖ = Ψᶜᵖ,
                                  Nᵐʰ = Nᵐʰ, Nᶜᵖ = Nᶜᵖ,
                                  fₒ = f⃰, f̅ = f̅ₜ(f⃰, -Inf), f̅′ = f̅ₜ(f̅′, -Inf),
                                  pᵤ = pᵤ), 
                                  true)
    println("   Maximizing distribution capacity with outsourcing: $TS′")
    println("       Distribution capacity: $(round(N̄₁′, digits=2))\n")

    for t ∈ 1:(tˢ-1)
        println("   \nDay-$t")
        println("   Demand: $(Nₜ(-Inf))")
        println("   Total cost per package: $(round(TC₁⃰, digits=2))")
        N̄[:,t]   .= N̄₁
        Ñ[:,t]   .= 0.0
        TC⃰[:,t] .= TC₁⃰
    end
    

    
    # ─────────────────────────────────────────────────────────────────
    ## Peri-disruption phase
    println("\n\n\n-----PERI - DISRUPTION-----")
    N₂ᵘ = 0.0
    for t ∈ tˢ:tᵉ
        println("   \nDay-$t")
        # Distribution capacity without outsourcing
        TS, N̄₂, _, _, _, T̄C̄₂ = LMO(:maxN, 
                                  (ρₓ = ρₓ⃰, φ = φₜ(t), 
                                   fₒ = f⃰, f̅ = f̅ₜ(f⃰, t), f̅′ = 0), 
                                   true)
        println("       Maximizing distribution capacity without outsourcing: $TS")
        println("           Distribution capacity: $(round(N̄₂, digits=2))\n")

        # Distribution capacity with outsourcing
        TS′, N̄₂′, _, _, _, T̄C̄₂′ = LMO(:maxN, 
                                     (ρₓ = ρₓ⃰, φ = φₜ(t),
                                      Ψᶜˢ = Ψᶜˢ, Ψᵐʰ = Ψᵐʰ, Ψᶜᵖ = Ψᶜᵖ,
                                      Nᵐʰ = Nᵐʰ, Nᶜᵖ = Nᶜᵖ,
                                      fₒ = f⃰, f̅ = f̅ₜ(f⃰, t), f̅′ = f̅ₜ(f̅′, t),
                                      pᵤ = pᵤ), 
                                      true)
        println("       Maximizing distribution capacity with outsourcing: $TS′")
        println("           Distribution capacity: $(round(N̄₂′, digits=2))\n")

        N₂ = Nₜ(t) + N₂ᵘ
        println("       Total demand: $(round(N₂, digits=2))")
        if N₂ >= N̄₂′
            println("       Total demand >= Combined distribution capacity")
            TC   = T̄C̄₂ * N₂/N̄₂
            TC′  = T̄C̄₂′ * N₂/N̄₂′
            Nᵘ   = N₂ - N̄₂
            Nᵘ′  = N₂ - N̄₂′
            TC₂⃰ = TC′
            N₂ᵘ  = Nᵘ′
            bool = (Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ == 1) && (abs(TC′- TC) > 1e-5) ? true : false 
            println("           Outsourced: $bool\n")
            N̄[:,t] .= N̄₂′
        elseif N₂ >= N̄₂
            println("       Total demand ∈ (Combined, E-retailer's] distribution capacity")
            TC = T̄C̄₂ * N₂/N̄₂
            TS′, _, _, _, _, TC′ = LMO(:minTC, 
                                      (N = N₂, ρₓ = ρₓ⃰, φ = φₜ(t),
                                       Ψᶜˢ = Ψᶜˢ, Ψᵐʰ = Ψᵐʰ, Ψᶜᵖ = Ψᶜᵖ,
                                       Nᵐʰ = Nᵐʰ, Nᶜᵖ = Nᶜᵖ,
                                       fₒ = f⃰, f̅ = f̅ₜ(f⃰, t), f̅′ = f̅ₜ(f̅′, t),
                                       pₗ = floor(1 - N̄₂/N₂, digits=1), pᵤ = pᵤ),
                                       true)
            Nᵘ   = N₂ - N̄₂
            Nᵘ′  = 0.0
            TC₂⃰ = TC′
            N₂ᵘ  = Nᵘ′
            bool = (Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ == 1) && (abs(TC′- TC) > 1e-5) ? true : false
            println("           Minimizing total cost per package: $TS′")
            println("           Outsourced: $bool\n")
            N̄[:,t] .= N̄₂′
        else
            println("       Total demand < E-retailer's distribution capacity")
            TS, _, _, _, _, TC = LMO(:minTC, 
                                    (N = N₂, ρₓ = ρₓ⃰, φ = φₜ(t), 
                                     fₒ = f⃰, f̅ = f̅ₜ(f⃰, t)), 
                                     true)                        
            Nᵘ = 0.0
            TC₂⃰ = TC
            N₂ᵘ = Nᵘ
            println("       Minimizing total cost per package: $TS")
            println("       Outsourced: false\n")
            N̄[:,t] .= N̄₂
        end
        println("       Unmet demand: $(round(N₂ᵘ, digits=2))")
        println("       Total cost per package: $(round(TC₂⃰, digits=2))\n")
        Ñ[:,t]   .= N₂ᵘ
        TC⃰[:,t] .= TC₂⃰
    end

    
    # ────────────────────────────────────────────────────────────────────────────────
    ## Post-disruption
    ρₓ⃰₁ = ρₓ⃰
    for re_optimize ∈ parameters[:re_optimize]
        r = findfirst(x -> x == re_optimize, parameters[:re_optimize])::Int64
        println("\n\n\n-----POST - DISRUPTION-----")
        println("   Re-optimizing operations: $re_optimize")
        if re_optimize == :operational
            ρₓ⃰ = ρₓ⃰₁
            TS, = LMO(:minTC, 
                     (N = Nₜ(Inf), ρₓ = ρₓ⃰, φ = φₜ(Inf), 
                      fₒ = f₁⃰, f̅ = f̅ₜ(f₁⃰, Inf)),
                      true)
            f⃰ = Int(f₁⃰)
            println("       Minimizing total cost per package: $TS\n")
        end
        if re_optimize == :tactical
            ρₓ⃰ = ρₓ⃰₁
            TS, _, _, f₃⃰, = LMO(:minTC, 
                                (N = Nₜ(Inf), ρₓ = ρₓ⃰, φ = φₜ(Inf)), 
                                true)
            f⃰ = Int(f₃⃰)
            println("       Minimizing total cost per package: $TS\n")
        end
        if re_optimize == :strategic
            TS, _, ρₓ⃰, f₃⃰, = LMO(:minTC, 
                                  (N = Nₜ(Inf), φ = φₜ(Inf)), 
                                  true)
            f⃰ = Int(f₃⃰)
            println("       Minimizing total cost per package: $TS\n")
        end
        
        N₃ᵘ = N₂ᵘ
        t = tᵉ + 1
        while true
            if t > length(N̄[r, :])   N̄  = [N̄   fill(NaN, l)] end
            if t > length(Ñ[r, :])   Ñ  = [Ñ   fill(NaN, l)] end
            if t > length(TC⃰[r, :]) TC⃰= [TC⃰ fill(NaN, l)] end
            
            println("   \nDay-$t")
            # Distribution capacity without outsourcing
            TS, N̄₃, _, _, _, T̄C̄₃ = LMO(:maxN, 
                                      (ρₓ = ρₓ⃰, φ = φₜ(t), 
                                       f̅ = f⃰, fₒ = f̅ₜ(f⃰, t), f̅′ = 0), 
                                       true)
            println("       Maximizing distribution capacity without outsourcing: $TS")
            println("           Distribution capacity: $(round(N̄₃, digits=2))\n")

            # Distribution capacity with outsourcing
            TS′, N̄₃′, _, _, _, T̄C̄₃′ = LMO(:maxN, 
                                         (ρₓ = ρₓ⃰, φ = φₜ(t),
                                          Ψᶜˢ = Ψᶜˢ, Ψᵐʰ = Ψᵐʰ, Ψᶜᵖ = Ψᶜᵖ,
                                          Nᵐʰ = Nᵐʰ, Nᶜᵖ = Nᶜᵖ,
                                          fₒ = f⃰, f̅ = f̅ₜ(f⃰, t), f̅′ = f̅ₜ(f̅′, t),
                                          pᵤ = pᵤ), 
                                          true)
            println("       Maximizing distribution capacity with outsourcing: $TS′")
            println("           Distribution capacity: $(round(N̄₃′, digits=2))\n")

            N₃ = Nₜ(t) + N₃ᵘ
            println("       Total demand: $(round(N₃, digits=2))")
            if N₃ >= N̄₃′
                println("       Total demand >= Combined distribution capacity")
                TC   = T̄C̄₃ * N₃/N̄₃
                TC′  = T̄C̄₃′ * N₃/N̄₃′
                Nᵘ   = N₃ - N̄₃
                Nᵘ′  = N₃ - N̄₃′
                TC₃⃰ = TC′
                N₃ᵘ  = Nᵘ′
                bool = (Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ == 1) && (abs(TC′- TC) > 1e-5) ? true : false
                println("           Outsourced: $bool\n")
                N̄[r, t]= N̄₃′
            elseif N₃ >= N̄₃
                println("       Total demand ∈ (Combined, E-retailer's] distribution capacity")
                TC = T̄C̄₃ * N₃/N̄₃
                TS′, _, _, _, _, TC′ = LMO(:minTC, 
                                          (N = N₃, ρₓ = ρₓ⃰, φ = φₜ(t),
                                           Ψᶜˢ = Ψᶜˢ, Ψᵐʰ = Ψᵐʰ, Ψᶜᵖ = Ψᶜᵖ,
                                           Nᵐʰ = Nᵐʰ, Nᶜᵖ = Nᶜᵖ,
                                           fₒ = f⃰, f̅ = f̅ₜ(f⃰, t), f̅′ = f̅ₜ(f̅′, t),
                                           pₗ = floor(1 - N̄₃/N₃, digits=1), pᵤ = pᵤ),
                                           true)
                Nᵘ   = N₃ - N̄₃
                Nᵘ′  = 0.0
                TC₃⃰ = TC′
                N₃ᵘ  = Nᵘ′
                bool = (Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ == 1) && (abs(TC′- TC) > 1e-5) ? true : false
                println("           Minimizing total cost per package: $TS′")
                println("           Outsourced: $bool\n")
                N̄[r, t]= N̄₃′
            else
                println("       Total demand < E-retailer's distribution capacity")
                TS, _, _, _, _, TC   = LMO(:minTC, 
                                          (N = N₃, ρₓ = ρₓ⃰, φ = φₜ(t), 
                                           fₒ = f⃰, f̅ = f̅ₜ(f⃰, t)), 
                                           true)
                TC₃⃰ = TC
                N₃ᵘ = 0.0
                println("           Minimizing total cost per package: $TS")
                println("           Outsourced: false\n")
                N̄[r, t]= N̄₃
            end
            println("       Unmet demand: $(round(N₃ᵘ, digits=2))")
            println("       Total cost per package: $(round(TC₃⃰, digits=2))")
            Ñ[r, t]  = N₃ᵘ
            TC⃰[r, t]= TC₃⃰

            
            if t >= Tᴸ
                convergence = max(N̄₃, N̄₃′) >= Nₜ(Inf)
                if convergence
                    if (TC⃰[r, t-1] - ε <= TC⃰[r, t] <= TC⃰[r, t-1] + ε) break end
                else break
                end
            end
            t += 1
        end
    end
    # ─────────────────────────────────────────────────────────────────
    filename = joinpath(@__DIR__, "Results\\sim_data.jld")
    save(filename, "N̄", N̄, "Ñ", Ñ, "TC⃰", TC⃰, "re_optimize", parameters[:re_optimize])
    return N̄, Ñ, TC⃰
end