"""
    LMO(objective::Symbol, given::NamedTuple, silent::Bool=true)

Last-Mile Optimization for a generic multi-echelon distribution structure.

# Decision variables
- Number of stops to be made in a delivery tour
- Number of delivery tours to be made by a delveiry vehicle
- Size of fleet
- Facility location
- Outsourcing channel usage

# Arguments
- `objective::Symbol`                       : The objective must either be to minimize total cost `:minTC`, or to maximize distribution capacity `:maxN`.
- `given::NamedTuple`    
    ### To minimize total cost:
    - `N` (Market size)                     : Required.
    - `ρₓ` (Facility location)              : Optional. If given, facility location is treated as a parameter else as a decision variable.
    - `φ` (Congestion factor)               : Optional. Defaulted at `φ = 1.0`.
    - `Ψᶜˢ` (Crowd-shipping binary)         : Optional. Defaulted at `Ψᶜˢ = 0`. To outsource via crowd-shipping `Ψᶜˢ = 1` must be given.
    - `Ψᵐʰ` (Micro-hub binary)              : Optional. Defaulted at `Ψᵐʰ = 0`. To outsource via micro-hub `Ψᵐʰ = 1` must be given.
    - `Ψᶜᵖ` (Collection-point binary)       : Optional. Defaulted at `Ψᶜᵖ = 0`. To outsource via collection-point `Ψᶜᵖ = 1` must be given.
    - `Nᵐʰ` (Number of micro-hubs)          : Optional. Defaulted at `Nᵐʰ = Ψᵐʰ`.
    - `Nᶜᵖ` (Number of collection-points)   : Optional. Defaulted at `Nᶜᵖ = Ψᶜᵖ`.
    - `fₒ` (Available trucks)               : Optional. Defaulted at decision variable `f`.
    - `f̅` (Available truck drivers)         : Optional. Defaulted at `N`.
    - `f̅′` (Available outsource drivers)    : Optional. Defaulted at `f̅′ = Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ`.
    - `pₗ` (Min outsource demand share)     : Optional. Defaulted at `pₗ = 0`.
    - `pᵤ` (Max outsource demand share)     : Optional. Defaulted at `pᵤ = Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ`.
    ### To maximize distribution capacity:
    - `ρₓ` (Facility location)              : Optional. If given, facility location is treated as a parameter else as a decision variable.
    - `φ` (Congestion factor)               : Optional. Defaulted at `φ = 1.0`.
    - `Ψᶜˢ` (Crowd-shipping binary)         : Optional. Defaulted at `Ψᶜˢ = 0`. When outsourcing via crowd-shipping `Ψᶜˢ = 1` must be given.
    - `Ψᵐʰ` (Micro-hub binary)              : Optional. Defaulted at `Ψᵐʰ = 0`. When outsourcing via micro-hub `Ψᵐʰ = 1` must be given.
    - `Ψᶜᵖ` (Collection-point binary)       : Optional. Defaulted at `Ψᶜᵖ = 0`. When outsourcing via collection-point `Ψᶜᵖ = 1` must be given.
    - `Nᵐʰ` (Number of micro-hubs)          : Optional. Defaulted at `Nᵐʰ = Ψᵐʰ`.
    - `Nᶜᵖ` (Number of collection-points)   : Optional. Defaulted at `Nᶜᵖ = Ψᶜᵖ`.
    - `fₒ` (Available trucks)               : Required.
    - `f̅` (Available truck drivers)         : Required.
    - `f̅′` (Available outsource drivers)    : Required.
    - `pₗ` (Min outsource demand share)     : Optional. Defaulted at `pₗ = 0`.
    - `pᵤ` (Max outsource demand share)     : Optional. Defaulted at `pᵤ = Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ`.
- `silent::Bool=true` if `false` prints output of the optimization model.

"""
function LMO(objective::Symbol, given::NamedTuple, silent::Bool)
    ## Solver
    optimizer = Juniper.Optimizer
    nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    model = Model(optimizer_with_attributes(optimizer,"nl_solver" => nl_solver))
    set_optimizer_attributes(model)
    set_silent(model)
    register(model, :ceil, 1, ceil; autodiff = true)
    ϵ = 1e-4
    # ────────────────────────────────────────────────────────────────────────────────
    ## Parameters
    # General parameters
    :φ ∈ keys(given) ? φ = given[:φ] : φ = 1.0                                          # Congestion factor
    :Ψᶜˢ ∈ keys(given) ? Ψᶜˢ = given[:Ψᶜˢ] : Ψᶜˢ = 0.0                                  # Crowd-shipping binary
    :Ψᵐʰ ∈ keys(given) ? Ψᵐʰ = given[:Ψᵐʰ] : Ψᵐʰ = 0.0                                  # Micro-hub binary
    :Ψᶜᵖ ∈ keys(given) ? Ψᶜᵖ = given[:Ψᶜᵖ] : Ψᶜᵖ = 0.0                                  # Collection-point binary
    :Nᵐʰ ∈ keys(given) ? Nᵐʰ = given[:Nᵐʰ] : Nᵐʰ = Ψᵐʰ                                  # Number of micro-hubs
    :Nᶜᵖ ∈ keys(given) ? Nᶜᵖ = given[:Nᶜᵖ] : Nᶜᵖ = Ψᶜᵖ                                  # Number of collection-points
    :f̅ ∈ keys(given) ? f̅ = given[:f̅] : f̅ = given[:N]                                    # Available fleet size
    :f̅′ ∈ keys(given) ? f̅′= given[:f̅′] : f̅′ = given[:N] * (Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ)             # Available outsource fleet
    :pₗ ∈ keys(given) ? pₗ = given[:pₗ] : pₗ = 0.0                                         # Min outsource demand share
    :pᵤ ∈ keys(given) ? pᵤ = given[:pᵤ] : pᵤ = Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ                          # Max outsource demand share
    if :fₒ ∈ keys(given) f̅ = min(f̅, given[:fₒ]) end
    pₗ = min(pₗ, pᵤ)
    # Market parameters
    A = 475                                                                             # Service region size (square miles)
    θ = 3                                                                               # Number of customers servied per vehicle stop
    W = 9                                                                               # Working hours in a day
    d = 330                                                                             # Number of working days
    Y = 10                                                                              # Time-horizon (Planned years of operation)
    r = 0.03                                                                            # Discount rate
    # E-retailer fleet parameters
    vᵢ = 20                                                                             # Vehicle free-flow speed inside the service region (mph)
    vₒ = 55                                                                             # Vehicle free-flow speed outside the service region (mph)
    VC = 360                                                                            # Vehicle capacity (number of customers)
    PC = 80000                                                                          # Vehicle purchase cost
    τᶜ = 1.0/60                                                                         # Service time at customer (hours)
    τᶠ = 0.3/60                                                                         # Service time at facility (hours)
    rᶠ = 0.1                                                                            # Rate of fuel consumption (gallon/mile)
    rᵉ = 1                                                                              # Rate of emission (grams/mile)
    # E-retailer cost parameters (class 5 diesel truck)
    πᵈ = 35                                                                             # Hourly driver cost
    πᵐ = 0.2                                                                            # Maintenance cost per mile
    πᶠ = 3.86                                                                           # Fuel cost per gallon
    πᵉ = 0.471                                                                          # Emission cost per gram
    # Outsource fleet parameters (cs:pickup-truck, mh:cargo-bike/diesel-van, cp:personal car)
    vᵢ′ = sum([    24,      9,   24] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ]) + 1                            # Vehicle free-flow speed inside the service region (mph)
    vₒ′ = sum([    59,      9,   59] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ]) + 1                            # Vehicle free-flow speed outside the service region (mph)
    VC′ = sum([    30,     30,    1] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                                # Vehicle capacity (number of customers)
    PC′ = sum([     0,   9500,    0] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                                # Vehcile purchase cost
    τᶜ′ = sum([0.5/60, 0.5/60, 1/60] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                                # Service time at customer (hours)
    τᶠ′ = sum([0.5/60, 0.3/60,    0] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                                # Service time at facility (hours)
    rᶠ′ = sum([  0.05,   0.29, 0.03] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                                # Rate of fuel consumption (gallon/mile)
    rᵉ′ = sum([     1,      0,    1] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                                # Rate of emission (grams/mile)
    # Outsource cost parameters (cs:pickup-truck, mh:cargo-bike/diesel-van, cp:personal car)
    πᵈ′ = sum([    35,     35,      0] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                              # Hourly driver cost
    πᵐ′ = sum([   0.0,   0.02,      0] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                              # Maintenance cost per mile
    πᶠ′ = sum([   0.0,   0.12,      0] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                              # Fuel cost per gallon
    πᵉ′ = sum([0.0405,    0.0, 0.0273] .* [Ψᶜˢ, Ψᵐʰ, Ψᶜᵖ])                              # Emission cost per gram
    # ────────────────────────────────────────────────────────────────────────────────
    ## Variables
    objective == :minTC ? λ = 1 : λ = 0
    # Facility location
    @variable(model, ϕ, Bin, start = 0)                                                                                                     # Binary variable for ρ and t; ϕ = 0 if depot is inside the service region else 1                                                                                                # Binary variable for k; α₄ = 1 if S >= 4 else 0
    :ρₓ ∈ keys(given) ? @variable(model, ρₓ == given[:ρₓ]) : @variable(model, 0.001 <= ρₓ <= 50)                                            # Facility location                   
    # Outsourcing
    pₗ == pᵤ ? @variable(model, p == λ * pₗ + (1 - λ) * pᵤ) : @variable(model, pₗ <= p <= pᵤ, start = λ * pₗ + (1 - λ) * pᵤ)                   # Share of customers served via crowdshipping service
    γ′ = @NLexpression(model, ceil(p))                                                                                                      # Binary variable (simplified to a parameter) for use of outsourcing channel       
    # E-retailer variables
    @variable(model, 0 <= f <= f̅, Int, start = f̅)                                                                                           # Fleet size
    @variable(model, (1 - λ) * Ψᶜˢ <= m <= W * vᵢ/sqrt(A), Int)                                                                             # Delivery tours per vehicle 
    objective == :minTC ? Cᶜ = @NLexpression(model, given[:N] * (1 - p)/(m * f + ϵ)) : @variable(model, 0 <= Cᶜ <= VC, start = VC)          # Customers served directly per delivery tour
    Cᵐʰ = @NLexpression(model, Nᵐʰ * γ′/(m * f + ϵ))                                                                                        # Micro-hubs per delivery tour
    Cᶜᵖ = @NLexpression(model, Nᶜᵖ * γ′/(m * f + ϵ))                                                                                        # Collecton-points per delivery tour
    @variable(model, 0 <= Sᶜ <= VC/θ, Int)                                                                                                  # Customer stops per delivery tour
    Sᵐʰ = @NLexpression(model, ceil(Cᵐʰ))                                                                                                   # Micro-hub stops per delivery tour
    Sᶜᵖ = @NLexpression(model, ceil(Cᶜᵖ))                                                                                                   # Collection-point stops per delivery tour
    γ  = @NLexpression(model, ceil(1 - Ψᶜˢ * p))                                                                                            # Binary variable (simplified to a parameter) for use of e-retailer's fleet
    # Outsource variables
    m̄′ = (min(ceil(W * vᵢ′/sqrt(A)), 1) * Ψᶜˢ + ceil(W * vᵢ′/sqrt(A) * 2/3 * sqrt(Nᵐʰ)) * Ψᵐʰ) + Ψᶜᵖ                                                                            
    if Ψᶜˢ == 1
        @variable(model, 0 <= f′ <= f̅′, Int, start = (1 - λ) * f̅′)                                                                          # Fleet size
        @variable(model, 0 <= m′ <= m̄′, Int)                                                                                                # Delivery tours per vehicle 
        @variable(model, 0 <= Cᶜ′<= VC′, start = VC′)                                                                                       # Customers served per delivery tour
    elseif Ψᵐʰ == 1
        @variable(model, Nᵐʰ <= f′ <= f̅′, Int, start = (1 - λ) * f̅′ + λ * Nᵐʰ)                                                              # Fleet size
        @variable(model, 0 <= m′ <= m̄′, Int)                                                                                                # Delivery tours per vehicle 
        @variable(model, 0 <= Cᶜ′<= VC′, start = VC′)                                                                                       # Customers served per delivery tour
    elseif Ψᶜᵖ == 1
        @variable(model, 0 <= f′ <= f̅′, start = (1 - λ) * f̅′)                                                                               # Fleet size
        m′ = @NLexpression(model, 1)                                                                                                        # Delivery tours per vehicle 
        Cᶜ′= @NLexpression(model, 1)                                                                                                        # Customers served per delivery tour
    else
        f′ = @NLexpression(model, 0)                                                                                                        # Fleet size
        m′ = @NLexpression(model, 0)                                                                                                        # Delivery tours per vehicle 
        Cᶜ′= @NLexpression(model, 0)                                                                                                        # Customers served per delivery tour
    end
    Sᶜ′ = @NLexpression(model, ceil(Cᶜ′/θ))                                                                                                 # Customer stops per delivery tour
    # ────────────────────────────────────────────────────────────────────────────────
    ## Auxiliary parameters
    objective == :minTC ? N = @NLexpression(model, given[:N]) : N = @NLexpression(model, Cᶜ * m * f + Cᶜ′ * m′ * f′)                        # Market size
    η  = (1 - (1 + r)^(-Y))/r * d                                                                                                           # Amortization factor 
    :fₒ ∈ keys(given) ? fₒ = @NLexpression(model, given[:fₒ]) : fₒ = @NLexpression(model, f)                                                # Purchased fleet
    # ────────────────────────────────────────────────────────────────────────────────
    ## Operations
    # E-retailer operational parameters
    k  = 0.73                                                                                                                               # Last-mile coefficient
    nᵐʰ= @NLexpression(model, Ψᵐʰ * Nᵐʰ + (1 - Ψᵐʰ))                                                                                        # Adjusted number of micro-hubs
    nᶜᵖ= @NLexpression(model, Ψᶜᵖ * Nᶜᵖ + (1 - Ψᶜᵖ))                                                                                        # Adjusted number of collection-points
    Cᶠ = @NLexpression(model, Cᵐʰ * N * p/nᵐʰ + Cᶜᵖ * N * p/nᶜᵖ)                                                                            # Customers served through micro-hubs and/or collection-points in a delivery tour
    C  = @NLexpression(model, Cᶜ + Cᶠ)                                                                                                      # Total customers served in a delivery tour
    S  = @NLexpression(model, Sᶜ + Sᵐʰ + Sᶜᵖ)                                                                                               # Total stops served in a delivery tour
    δ  = @NLexpression(model, N/A * (1 - p)/θ + Nᵐʰ  * γ′/A + Nᶜᵖ  * γ′/A + ϵ)                                                              # Stop density
    ρ  = @NLexpression(model, (ρₓ + sqrt(A)/4) * ϕ + (ρₓ^2/sqrt(A) + sqrt(A)/2) * (1 - ϕ))                                                  # Long-haul length
    t  = @NLexpression(model, 1/φ * ((ρₓ/vₒ + sqrt(A)/4 * (3/vᵢ - 2/vₒ)) * ϕ + (ρ/vᵢ) * (1 - ϕ)))                                           # Long-haul travel time
    L  = @NLexpression(model, 2 * ρ + k * S/sqrt(δ))                                                                                        # Delivery tour length
    T  = @NLexpression(model, C * τᶠ + 2 * t + k * S/(vᵢ * φ * sqrt(δ)) + Cᶜ * τᶜ + Cᶠ * τᶠ)                                                # Delivery tour travel time
    # Outsource operational parameters
    k′ = 0.73                                                                                                                               # Last-mile coefficient
    δ′ = @NLexpression(model, N/A * p/θ + ϵ)                                                                                                # Stop density
    ρ′ = @NLexpression(model, Ψᶜˢ * ρ + Ψᵐʰ * (2/3 * sqrt(A/nᵐʰ)) + Ψᶜᵖ * (2/3 * sqrt(A/nᶜᵖ)))                                              # Long-haul length
    t′ = @NLexpression(model, 1/φ * (Ψᶜˢ * ((ρₓ/vₒ′ + sqrt(A)/4 * (3/vᵢ′ - 2/vₒ′)) * ϕ + (ρ/vᵢ′) * (1 - ϕ)) + (Ψᵐʰ + Ψᶜᵖ) * (ρ′/vᵢ′)))      # Long-haul travel time
    L′ = @NLexpression(model, 2 * ρ′ + (Ψᶜˢ + Ψᵐʰ) * k′ * Sᶜ′/sqrt(δ′))                                                                     # Delivery tour length
    T′ = @NLexpression(model, (Ψᶜˢ + Ψᵐʰ) * Cᶜ′ * τᶠ′ + 2 * t′ + (Ψᶜˢ + Ψᵐʰ) * k′ * Sᶜ′/(vᵢ′ * φ * sqrt(δ′)) + Cᶜ′ * τᶜ′)                   # Delivery tour travel time
    # ────────────────────────────────────────────────────────────────────────────────
    ## Costs
    # E-retailer costs
    FC = @NLexpression(model, (356.37 * ρₓ ^ (-0.231) * N * 5 + PC * fₒ)/(η * N))                                                           # Fixed cost
    OC = @NLexpression(model, m * f * (T * πᵈ + L * (πᵐ + rᶠ * πᶠ))/N)                                                                      # Operation cost per package
    EC = @NLexpression(model, m * f * L * (rᵉ * πᵉ)/N)                                                                                      # Emission cost per package
    # Outsource costs
    FC′= @NLexpression(model, ((458.648 * (sqrt(A) ^ (-0.231)) * N * p * 5 + (Nᵐʰ + Nᶜᵖ) * 3e5) * (Ψᵐʰ + Ψᶜᵖ) + PC′ * f′) * γ′/(η * N))     # Fixed cost per package
    OC′= @NLexpression(model, m′ * f′ * (T′ * πᵈ′ + L′ * (πᵐ′ + rᶠ′ * πᶠ′))/N)                                                              # Operation cost per package
    EC′= @NLexpression(model, m′ * f′ * L′ * (rᵉ′ * πᵉ′)/N)                                                                                 # Emission cost per package
    # Total costs
    TC = @NLexpression(model, (FC + EC + OC) + (FC′ + EC′ + OC′))                                                                           # Total cost per package
    # ────────────────────────────────────────────────────────────────────────────────
    ## Constriants
    # E-retailer constriants
    if objective == :maxN @NLconstraint(model, N >= ϵ) end                                             
    @NLconstraint(model, ρₓ - sqrt(A)/2 <= 1000 * ϕ)
    @NLconstraint(model, ρₓ - sqrt(A)/2 >= -1000 * (1-ϕ))
    @NLconstraint(model, C <= VC)                                                                       # Capacity constraint
    @NLconstraint(model, T * m <= W)                                                                    # Working hours constraint
    if objective == :maxN @NLconstraint(model, Cᶜ * m * f == N * (1 - p)) end                           # Level of service constraint
    @NLconstraint(model, Sᶜ >= Cᶜ/θ)                                                                    # Number of customer stops constraint
    # Outsource constriants
    @NLconstraint(model, T′ * m′ <= W)                                                                  # Working hours constraint
    if objective == :minTC && Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ == 1 @NLconstraint(model, Cᶜ′ * m′ * f′ == N * p) end     # Level of service constraint
    # ────────────────────────────────────────────────────────────────────────────────
    ## Optimize
    objective == :minTC ? @NLobjective(model, Min, TC) : @NLobjective(model, Max, N)
    JuMP.optimize!(model)
    TS, N⃰, ρₓ⃰, f⃰, f′⃰, TC⃰ = termination_status(model), value(N), value(ρₓ), value(f), value(f′), value(TC)
    # Lite mode
    if TS ≠ MOI.LOCALLY_SOLVED
        printstyled("Switching to Lite mode\n", color=:yellow)
        @NLconstraint(model, Cᶜ′ == VC′ * (Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ))
        JuMP.optimize!(model)
        TS, N⃰, ρₓ⃰, f⃰, f′⃰, TC⃰ = termination_status(model), value(N), value(ρₓ), value(f), value(f′), value(TC)
    end

    if TS ≠ MOI.LOCALLY_SOLVED
        printstyled("Termination Status: ", color=:yellow)
        println(TS)
    else 
        if !silent
            println("E-retailer parameters")
            println("   Facility location:              $(round(value(ρₓ), digits=2))")
            println("   Customers per tour:             $(round(value(Cᶜ)  * value(γ), digits=2))")
            println("   Micro-hubs per tour:            $(round(value(Cᵐʰ) * value(γ), digits=2))")
            println("   Collection-points per tour:     $(round(value(Cᶜᵖ) * value(γ), digits=2))")
            println("   Total customers per tour:       $(round(value(C)   * value(γ), digits=2))")
            println("   Total stops per tour:           $(round(value(S)   * value(γ), digits=2))")
            println("   Tours per vehicle:              $(round(value(m)   * value(γ), digits=2))")
            println("   Fleet size:                     $(round(value(f)   * value(γ), digits=2))")
            println("   Delivery tour length:           $(round(value(L)   * value(γ), digits=2))")
            println("   Delivery tour TT:               $(round(value(T)   * value(γ), digits=2))")
            println("   Fixed cost per package:         $(round(value(FC), digits=2))")
            println("   Operational cost per package:   $(round(value(OC), digits=2))")
            println("   Emission cost per package:      $(round(value(EC), digits=2))")
            println("   Outsourcing share:              $(round(value(p),  digits=2))")
            println("Outsourcing parameters")
            println("   Customers per tour:             $(round(value(Cᶜ′) * value(γ′), digits=2))")
            println("   Total stops per tour:           $(round(value(Sᶜ′) * value(γ′), digits=2))")
            println("   Tours per vehicle:              $(round(value(m′)  * value(γ′), digits=2))")
            println("   Fleet size:                     $(round(value(f′)  * value(γ′), digits=2))")
            println("   Delivery tour length:           $(round(value(L′)  * value(γ′), digits=2))")
            println("   Delivery tour TT:               $(round(value(T′)  * value(γ′), digits=2))")
            println("   Fixed cost per package:         $(round(value(FC′), digits=2))")
            println("   Operational cost per package:   $(round(value(OC′), digits=2))")
            println("   Emission cost per package:      $(round(value(EC′), digits=2))")
            println("Total cost per package:            $(round(value(TC),  digits=2))")
            println()
        end
    end
    return TS, N⃰, ρₓ⃰, f⃰, f′⃰, TC⃰
    # ────────────────────────────────────────────────────────────────────────────────
end