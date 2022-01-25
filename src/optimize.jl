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
    - `Ψᶜᵖ` (Collection-point binary)       : Optional. Defaulted at `Ψᶜᵖ = 0`. To outsource via collection-point `Ψᶜᵖ = 1` must be given.
    - `Ψᵐʰ` (Micro-hub binary)              : Optional. Defaulted at `Ψᵐʰ = 0`. To outsource via micro-hub `Ψᵐʰ = 1` must be given.
    - `Nᶜᵖ` (Number of collection-points)   : Optional. Defaulted at `Nᶜᵖ = Ψᶜᵖ`.
    - `Nᵐʰ` (Number of micro-hubs)          : Optional. Defaulted at `Nᵐʰ = Ψᵐʰ`.
    - `f̅` (Available truck drivers)         : Optional. Defaulted at `N`.
    - `f̅′` (Available outsource drivers)    : Optional. Defaulted at `f̅′ = Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ`.
    - `pₗ` (Min outsource demand share)     : Optional. Defaulted at `pₗ = 0`.
    - `pᵤ` (Max outsource demand share)     : Optional. Defaulted at `pᵤ = Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ`.
    ### To maximize distribution capacity:
    - `ρₓ` (Facility location)              : Optional. If given, facility location is treated as a parameter else as a decision variable.
    - `φ` (Congestion factor)               : Optional. Defaulted at `φ = 1.0`.
    - `Ψᶜˢ` (Crowd-shipping binary)         : Optional. Defaulted at `Ψᶜˢ = 0`. When outsourcing via crowd-shipping `Ψᶜˢ = 1` must be given.
    - `Ψᶜᵖ` (Collection-point binary)       : Optional. Defaulted at `Ψᶜᵖ = 0`. When outsourcing via collection-point `Ψᶜᵖ = 1` must be given.
    - `Ψᵐʰ` (Micro-hub binary)              : Optional. Defaulted at `Ψᵐʰ = 0`. When outsourcing via micro-hub `Ψᵐʰ = 1` must be given.
    - `Nᶜᵖ` (Number of collection-points)   : Optional. Defaulted at `Nᶜᵖ = Ψᶜᵖ`.
    - `Nᵐʰ` (Number of micro-hubs)          : Optional. Defaulted at `Nᵐʰ = Ψᵐʰ`.
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
    :Ψᶜᵖ ∈ keys(given) ? Ψᶜᵖ = given[:Ψᶜᵖ] : Ψᶜᵖ = 0.0                                  # Collection-point binary
    :Ψᵐʰ ∈ keys(given) ? Ψᵐʰ = given[:Ψᵐʰ] : Ψᵐʰ = 0.0                                  # Micro-hub binary
    :Nᶜᵖ ∈ keys(given) ? Nᶜᵖ = given[:Nᶜᵖ] : Nᶜᵖ = Ψᶜᵖ                                  # Number of collection-points
    :Nᵐʰ ∈ keys(given) ? Nᵐʰ = given[:Nᵐʰ] : Nᵐʰ = Ψᵐʰ                                  # Number of micro-hubs
    :f̅ ∈ keys(given) ? f̅ = given[:f̅] : f̅ = given[:N]                                    # Available fleet size
    :f̅′ ∈ keys(given) ? f̅′= given[:f̅′] : f̅′ = given[:N] * (Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ)             # Available outsource fleet
    :pₗ ∈ keys(given) ? pₗ = given[:pₗ] : pₗ = 0.0                                      # Min outsource demand share
    :pᵤ ∈ keys(given) ? pᵤ = given[:pᵤ] : pᵤ = Ψᶜˢ + Ψᵐʰ + Ψᶜᵖ                          # Max outsource demand share
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
    
    # Outsource fleet parameters (cs:pickup-truck, cp:personal car,  mh:cargo-bike)
    vᵢ′ = sum([    24,     24,      9] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ]) + 1                          # Vehicle free-flow speed inside the service region (mph)
    vₒ′ = sum([    59,     59,      9] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ]) + 1                          # Vehicle free-flow speed outside the service region (mph)
    VC′ = sum([    30,      1,     30] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                              # Vehicle capacity (number of customers)
    PC′ = sum([     0,      0,   9500] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                              # Vehcile purchase cost
    τᶜ′ = sum([0.5/60,   1/60, 0.5/60] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                              # Service time at customer (hours)
    τᶠ′ = sum([0.5/60,   0/60, 0.3/60] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                              # Service time at facility (hours)
    rᶠ′ = sum([  0.05,   0.03,   0.29] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                              # Rate of fuel consumption (gallon/mile)
    rᵉ′ = sum([     1,      1,      0] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                              # Rate of emission (grams/mile)
    
    # Outsource cost parameters (cs:pickup-truck, cp:personal car,  mh:cargo-bike)
    πᵈ′ = sum([    35,      0,   35] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                                # Hourly driver cost
    πᵐ′ = sum([   0.0,    0.0, 0.02] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                                # Maintenance cost per mile
    πᶠ′ = sum([   0.0,    0.0, 0.12] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                                # Fuel cost per gallon
    πᵉ′ = sum([0.0405, 0.0273,  0.0] .* [Ψᶜˢ, Ψᶜᵖ, Ψᵐʰ])                                # Emission cost per gram
    
    # Auxiliary parameters
    k = 0.73                                                                            # Last-mile coefficient
    k′= 0.73                                                                            # Last-mile coefficient
    V = 50                                                                              # Collection-point capacity
    η = (1 - (1 + r)^(-Y))/r * d                                                        # Amortization factor 
    
    
    # ────────────────────────────────────────────────────────────────────────────────
    ## Variables
    if objective == :minTC
        N = given[:N]

        # E-retailer
        m̄ = W * vᵢ/sqrt(A)
        ρₓ= :ρₓ ∈ keys(given) ? @NLexpression(model, given[:ρₓ]) : @variable(model, 0.001 ≤ ρₓ ≤ 50) 
        p = @variable(model, pₗ ≤ p ≤ pᵤ)
        f = :f̅ ∈ keys(given) ? @NLexpression(model, f̅) : @variable(model, 0 ≤ f ≤ f̅, Int, start = f̅)                                            
        m = @variable(model, 0 ≤ m ≤ m̄, Int, start = m̄)    
        Cᶜ= @NLexpression(model, N * (1 - p)/(m * f + ϵ))
        Sᶜ= @variable(model, 0 ≤ Sᶜ ≤ VC/θ, Int)
        
        γ = @NLexpression(model, ceil(1 - Ψᶜˢ * p))                                                                                             
        γ′= @NLexpression(model, ceil(p)) 
        ϕ = @variable(model, ϕ, Bin, start = 0)
        
        # Outsource
        if Ψᶜˢ == 1
            Cᶜ′= @variable(model, 0 ≤ Cᶜ′ ≤ VC′)
            m′ = @NLexpression(model, 1) 
            f′ = @NLexpression(model, N * p/(Cᶜ′ * m′ + ϵ))
            @NLconstraint(model, f′ ≤ f̅′)
        elseif Ψᶜᵖ == 1
            Cᶜ′= @NLexpression(model, 1) 
            m′ = @NLexpression(model, 1)
            f′ = @NLexpression(model, N * p)
        elseif Ψᵐʰ == 1
            mₗ′ = floor(W/((2 * 2/3 * sqrt(A/Nᵐʰ) + k′ * VC′/θ)/(vᵢ′ * 0.887) + VC′ * (τᶜ′ + τᶠ′)))
            mᵤ′= ceil(W /(2 * 2/3 * sqrt(A/Nᵐʰ)/(vᵢ′ * 0.982) + VC′ * (τᶜ′ + τᶠ′)))
            f′ = @variable(model, 0 ≤ f′ ≤ f̅′, Int)
            m′ = @variable(model, 0 ≤ m′ ≤ mᵤ′,Int) 
            Cᶜ′= @NLexpression(model, VC′)
            @NLconstraint(model, Cᶜ′ * m′ * f′ ==  N * p)
        else
            Cᶜ′= @NLexpression(model, 0)
            m′ = @NLexpression(model, 0)
            f′ = @NLexpression(model, 0)
        end
        Sᶜ′ = @NLexpression(model, ceil(Cᶜ′/θ)) 
        
        Cᵐʰ = @NLexpression(model, Nᵐʰ * γ′/(m * f + ϵ))                                                                                        
        Cᶜᵖ = @NLexpression(model, Cᶜ * (p / (1 - p)) / V * Ψᶜᵖ)                                                                                 
        Sᵐʰ = @NLexpression(model, ceil(Cᵐʰ))                                                                                                   
        Sᶜᵖ = @NLexpression(model, ceil(Cᶜᵖ))

        @NLconstraint(model, Sᶜ ≥ Cᶜ/θ)
    
    elseif objective == :maxN        
        # E-retailer
        m̄ = W * vᵢ/sqrt(A)
        p = @variable(model, pₗ ≤ p ≤ pᵤ)
        ρₓ= :ρₓ ∈ keys(given) ? @NLexpression(model, given[:ρₓ]) : @variable(model, 0.001 ≤ ρₓ ≤ 50) 
        f = @NLexpression(model, given[:f̅])                                        
        m = @variable(model, 0 ≤ m ≤ m̄, Int)
        Cᶜ= @variable(model, 0 ≤ Cᶜ ≤ VC)
        Sᶜ= @variable(model, 0 ≤ Sᶜ ≤ VC/θ, Int)
        
        γ = @NLexpression(model, ceil(1 - Ψᶜˢ * p))                                                                                             
        γ′= @NLexpression(model, ceil(p)) 
        ϕ = @variable(model, ϕ, Bin, start = 0)

        N = @NLexpression(model, Cᶜ * m * f/(1 - p))
        
        # Outsource
        if Ψᶜˢ == 1
            Cᶜ′= @variable(model, 0 ≤ Cᶜ′ ≤ VC′)
            m′ = @NLexpression(model, 1) 
            f′ = @NLexpression(model, f̅′)
            @NLconstraint(model, N * p == Cᶜ′ * f̅′)
        elseif Ψᶜᵖ == 1
            Cᶜ′= @NLexpression(model, 1) 
            m′ = @NLexpression(model, 1)
            f′ = @NLexpression(model, N * p)
        elseif Ψᵐʰ == 1
            mₗ′= floor(W/((2 * 2/3 * sqrt(A/Nᵐʰ) + k′ * VC′/θ)/(vᵢ′ * 0.887) + VC′ * (τᶜ′ + τᶠ′)))
            mᵤ′= ceil(W /(2 * 2/3 * sqrt(A/Nᵐʰ)/(vᵢ′ * 0.982) + VC′ * (τᶜ′ + τᶠ′)))
            f′ = @variable(model, 0 ≤ f′ ≤ f̅′, Int)
            m′ = @variable(model, mₗ′ ≤ m′ ≤ mᵤ′, start = mᵤ′, Int) 
            Cᶜ′= @NLexpression(model, VC′)
            @NLconstraint(model, Cᶜ′ * m′ * f′ * (1-p) == Cᶜ * m * f * p)
        else
            Cᶜ′= @NLexpression(model, 0)
            m′ = @NLexpression(model, 0)
            f′ = @NLexpression(model, 0)
        end
        Sᶜ′ = @NLexpression(model, ceil(Cᶜ′/θ)) 

        Cᵐʰ = @NLexpression(model, Nᵐʰ * γ′ /(m * f + ϵ))                                                                                        
        Cᶜᵖ = @NLexpression(model, Cᶜ * (p / (1 - p)) / V * Ψᶜᵖ)                                                                                 
        Sᵐʰ = @NLexpression(model, ceil(Cᵐʰ))                                                                                                   
        Sᶜᵖ = @NLexpression(model, ceil(Cᶜᵖ))

        @NLconstraint(model, Sᶜ ≥ Cᶜ/θ) 
    end
    # ────────────────────────────────────────────────────────────────────────────────
    ## Operations
    # E-retailer operational parameters
    nᵐʰ= @NLexpression(model, Ψᵐʰ * Nᵐʰ + (1 - Ψᵐʰ))                                                                                        # Adjusted number of micro-hubs
    nᶜᵖ= @NLexpression(model, Ψᶜᵖ * Nᶜᵖ + (1 - Ψᶜᵖ))                                                                                        # Adjusted number of collection-points               
    Cᶠ = @NLexpression(model, Cᵐʰ * N * p/nᵐʰ + Cᶜᵖ * V)
    C  = @NLexpression(model, Cᶜ + Cᶠ)                                                                                                      # Total customers served in a delivery tour
    S  = @NLexpression(model, Sᶜ + Sᵐʰ + Sᶜᵖ)                                                                                               # Total stops served in a delivery tour
    δ  = @NLexpression(model, N/A * (1 - p)/θ + Nᵐʰ * γ′/A + Nᶜᵖ * γ′/A + ϵ)                                                                # Stop density
    ρ  = @NLexpression(model, (ρₓ + sqrt(A)/4) * ϕ + (ρₓ^2/sqrt(A) + sqrt(A)/2) * (1 - ϕ))                                                  # Long-haul length
    t  = @NLexpression(model, 1/φ * ((ρₓ/vₒ + sqrt(A)/4 * (3/vᵢ - 2/vₒ)) * ϕ + (ρ/vᵢ) * (1 - ϕ)))                                           # Long-haul travel time
    L  = @NLexpression(model, 2 * ρ + k * S/sqrt(δ))                                                                                        # Delivery tour length
    T  = @NLexpression(model, C * τᶠ + 2 * t + k * S/(vᵢ * φ * sqrt(δ)) + Cᶜ * τᶜ + Cᶠ * τᶠ)                                                # Delivery tour travel time
    
    # Outsource operational parameters
    δ′ = @NLexpression(model, N/A * p/θ + ϵ)                                                                                                # Stop density
    ρ′ = @NLexpression(model, Ψᶜˢ * ρ + Ψᵐʰ * (2/3 * sqrt(A/nᵐʰ)) + Ψᶜᵖ * (2/3 * sqrt(A/nᶜᵖ)))                                              # Long-haul length
    t′ = @NLexpression(model, 1/φ * (Ψᶜˢ * ((ρₓ/vₒ′ + sqrt(A)/4 * (3/vᵢ′ - 2/vₒ′)) * ϕ + (ρ/vᵢ′) * (1 - ϕ)) + (Ψᵐʰ + Ψᶜᵖ) * (ρ′/vᵢ′)))      # Long-haul travel time
    L′ = @NLexpression(model, 2 * ρ′ + (Ψᶜˢ + Ψᵐʰ) * k′ * Sᶜ′/sqrt(δ′))                                                                     # Delivery tour length
    T′ = @NLexpression(model, (Ψᶜˢ + Ψᵐʰ) * Cᶜ′ * τᶠ′ + 2 * t′ + (Ψᶜˢ + Ψᵐʰ) * k′ * Sᶜ′/(vᵢ′ * φ * sqrt(δ′)) + Cᶜ′ * τᶜ′)                   # Delivery tour travel time
   
   
    # ────────────────────────────────────────────────────────────────────────────────
    ## Costs
    # E-retailer costs
    FC = @NLexpression(model, (356.37 * ρₓ ^ (-0.231) * 30000 * 5 + PC * f)/(η * N))                                                          # Fixed cost
    OC = @NLexpression(model, m * f * (T * πᵈ + L * (πᵐ + rᶠ * πᶠ))/N)                                                                      # Operation cost per package
    EC = @NLexpression(model, m * f * L * (rᵉ * πᵉ)/N)                                                                                      # Emission cost per package
    
    # Outsource costs
    FC′= @NLexpression(model, ((458.648 * (sqrt(A) ^ (-0.231)) * N * p * 5 + (3e5Nᵐʰ + 3.75e4Nᶜᵖ)) * (Ψᵐʰ + Ψᶜᵖ) + PC′ * f′) * γ′/(η * N))  # Fixed cost per package
    OC′= @NLexpression(model, m′ * f′ * (T′ * πᵈ′ + L′ * (πᵐ′ + rᶠ′ * πᶠ′))/N)                                                              # Operation cost per package
    EC′= @NLexpression(model, m′ * f′ * L′ * (rᵉ′ * πᵉ′)/N)                                                                                 # Emission cost per package
    
    # Total costs
    TC = @NLexpression(model, (FC + EC + OC) + (FC′ + EC′ + OC′))                                                                           # Total cost per package
    
    
    # ────────────────────────────────────────────────────────────────────────────────
    ## Constriants
    if objective == :maxN @NLconstraint(model, N ≥ ϵ) end                                             
    @NLconstraint(model, ρₓ - sqrt(A)/2 ≤ 1000 * ϕ)
    @NLconstraint(model, ρₓ - sqrt(A)/2 ≥ -1000 * (1-ϕ))
    @NLconstraint(model, C ≤ VC)                                                                        # Capacity constraint
    @NLconstraint(model, T * m ≤ W)                                                                     # Working hours constraint
    @NLconstraint(model, T′ * m′ - (Ψᵐʰ + Ψᶜᵖ) * 1/φ * ρ′/vᵢ′ ≤ W)                                      # Working hours constraint
    
    # ────────────────────────────────────────────────────────────────────────────────
    ## Optimize
    objective == :minTC ? @NLobjective(model, Min, TC) : @NLobjective(model, Max, N)
    JuMP.optimize!(model)
    TS, N⃰, ρₓ⃰, f⃰, p⃰, f′⃰, TC⃰ = termination_status(model), value(N), value(ρₓ), value(f), value(p), value(f′), value(TC)
   
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
        println("   Outsourcing share:              $(round(value(p),  digits=3))")
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
    return TS, N⃰, ρₓ⃰, f⃰, p⃰, f′⃰, TC⃰
    # ────────────────────────────────────────────────────────────────────────────────
end