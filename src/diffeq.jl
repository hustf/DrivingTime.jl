##############################################################
# Prepare and solve the differential equations for streamlines
##############################################################

function solve_guarded(rhs, Γᵢₙ, 𝐣::Journey, tspan, cbs; debug=false, odekws...)
    s = Unitful.s
    # Test the functions here. The internal __init will do this 
    # as well, but this keeps stacktraces shallow and checks
    # somewhat more strict.
    @inferred packout(Γᵢₙ)
    # Throw-away "derivative" w.r. time
    Γ´ = 0.01Γᵢₙ / s
    if debug
        with_logger(Logging.ConsoleLogger(stderr, Logging.Debug)) do
            x´, x´´ = packout(Γ´)
            @debug x´ x´´
        end
    end
    @inferred packin!(Γ´, 1.0u"m/s", 2.0u"m/s^2")
    @inferred rhs(Γ´, Γᵢₙ, 𝐣, nothing)
    # Define and solve
    prob = ODEProblem(rhs, Γᵢₙ, tspan, 𝐣, callback = cbs)
    sol = if debug
        with_logger(Logging.ConsoleLogger(stderr, Logging.Debug)) do
            solve(prob, Tsit5(); odekws...)
        end
    else
        solve(prob, Tsit5(), odekws...)
    end
    sol
end

"""
    solve_journey(rhs, Γᵢₙ, 𝐣::Journey; odekws...)
    solve_journey(Γᵢₙ::ArrayPartition, 𝐣::Journey;odekws...)
    solve_journey(𝐣::Journey;odekws...)
    ---> SciMLBase.ODESolution
"""
function solve_journey(rhs, Γᵢₙ, 𝐣::Journey; odekws...)
    tspan = make_tspan()
    cbs = callbacks_journey(𝐣; odekws...)
    # Drop the 'already spent' keywords comprising 'tspan'.
    # The remaining keywords will be passed on to the solver.
    remaining_kws = filter(odekws) do (kw, kwval)
        kw == :tstart && return false
        kw == :tstop && return false
        kw == :dtfloor && return false
        true
    end
    solve_guarded(rhs, Γᵢₙ, 𝐣, tspan, cbs; remaining_kws...)
end
solve_journey(Γᵢₙ::ArrayPartition, 𝐣::Journey;odekws...) = solve_journey(rhs!, Γᵢₙ, 𝐣; odekws...)
solve_journey(𝐣::Journey;odekws...) = solve_journey(ArrayPartition([0.0u"m"], [0.0u"m/s"]), 𝐣; odekws...)


"""
    make_tspan(; odekws...)
    --> Tuple{<:Time, <:Time}

Default: (0.0 s, 10800.0 s)

# Keyword arguments

- `tstart`
- `tstop`
- `tspan`
"""
function make_tspan(; odekws...)
    if :tspan ∈ keys(odekws)
        @assert :tstop ∉ keys(odekws) "Don't specify both tspan and tstop"
        tspan = odekws[:tspan]
    else
        t0 = get(odekws, :tstart, 0.0u"s")
        t1 = get(odekws, :tstop, u"s"(3.0u"hr"))
        tspan = (t0, t1)
    end
    @assert tspan isa Tuple{<:Time, <:Time}
    return tspan
end

function callbacks_journey(𝐣::Journey; odekws...)
    # Early end conditions
    vccb = ContinuousCallback[]
    #vccb = [ContinuousCallback(signed_distance_within_domain, affect!),
    #        ContinuousCallback(too_flat, affect!)]
    vdcb = DiscreteCallback[]
    if :dtfloor ∈ keys(odekws)
        push!(vdcb, let dtfloor = odekws[:dtfloor]
                DiscreteCallback(
                   (u, t, integrator ) -> get_proposed_dt(integrator) ≤ dtfloor, 
                       affect!)
            end)
    end
    CallbackSet(vccb..., vdcb...)
end



############################
# Callees for OrdinaryDiffEq
############################
packout(Γ::ArrayPartition) = Γ.x[1][1], Γ.x[2][1]
function packin!(Γ´::ArrayPartition, x´, x´´)
    Γ´.x[1][1] = x´
    Γ´.x[2][1] = x´´
    return Γ´
end

# ODE right-hand side
function rhs!(du, u, 𝐣::Journey, t)
    x, x´ = packout(u)
    @debug "u"    x    x´   maxlog = 2
    x´´ =     -g
    @debug "du "    x´´     maxlog = 2
    packin!(du, x´, x´´)
    du
end

# Our own function for potential debugging termination causes.
function affect!(integrator)
    terminate!(integrator)
end


# DEAD
#=
"""
    signed_distance_within_domain(u, t, integrator::ODEIntegrator)
    signed_distance_within_domain(saxy::SelectedVec2AtXY, x, y)
    signed_distance_within_domain(fxy::Journey, x, y)

Method extensions to Domain type argument.
"""
function signed_distance_within_domain(u, t, integrator::ODEIntegrator)
    @assert integrator.p isa Journey
    signed_distance_within_domain(integrator.p, u[1], u[2])
end
signed_distance_within_domain(saxy::SelectedVec2AtXY, x, y) =  signed_distance_within_domain(saxy.baxy.d, x, y)
signed_distance_within_domain(fxy::Journey, x, y) =  signed_distance_within_domain(fxy.d, x, y)

# Exit criterion. The hard coded criterion here is not
# intended for fine tuning. Instead, set a threshold 
# when defining the integrator functions (directional
# vector set to zero when under that threshold)
function too_flat(u, t, integrator::ODEIntegrator)
    out = max(-1.0, norm(u[1], u[2]) - 0.008)
    if out < 0f0
        @debug "Too flat ($out ) at (x, y) = $u  \n(xprev, yprev) = $(integrator.uprev) " maxlog = 120
    end
    out
end 



function condition_flip_bidirection(u, t, integrator::ODEIntegrator)
    @assert t > integrator.tprev
    is_close_to_opposite(integrator.p, integrator.uprev, u)
end
function is_close_to_opposite(saxy::SelectedVec2AtXY, u0, u1)
    d = dot_product_with_previous(saxy, u0, u1)
    is_close_to_opposite(d)
end

function affect_flip_bidirection!(integrator)
    integrator.p.flip[] = ! integrator.p.flip[] 
    integrator.p.v .=- integrator.p.v
    nothing
end

function condition_swap_major_minor(u, t, integrator::ODEIntegrator)
    @assert t > integrator.tprev
    dotprod = dot_product_with_previous(integrator.p, integrator.uprev, u)
    isit = is_close_to_perpendicular(dotprod)
    if isit
        @debug "Swap major minor at (x, y) = $u  \n(xprev, yprev) = $(integrator.uprev) " maxlog = 120
    else 
        if dotprod < 0.999
            @debug "In doubt: dotprod = $dotprod at (x, y) = $u  \n(xprev, yprev) = $(integrator.uprev) " maxlog = 120
            # Vector from last solution point
            Δu = (integrator.uprev - u)
            magΔ = norm(Δu)
            magΔ < MAG_EPS && return 1.0
            # ....normalized to unit length
            Δu ./= magΔ
            # Normalized move is zero? 
            du = integrator.p.v
            @debug "Δu = $(Δu)    du = $du"
        end
    end
    isit
end
function affect_swap_major_minor!(integrator)
    saxy = integrator.p
    saxy.baxy.major[] = ! saxy.baxy.major[]
    nothing
end

"""
    dot_product_with_previous(saxy::SelectedVec2AtXY, u0, u1)

Returns 1.0 if the dot product can't be evaluated.
"""
function dot_product_with_previous(saxy::SelectedVec2AtXY, u0, u1)
    # Vector from last solution point
    Δu = (u1 - u0)
    magΔ = norm(Δu)
    magΔ < MAG_EPS && return 1.0
    # ....normalized to unit length
    Δu ./= magΔ
    # The direction from here 
    du = saxy.v
    if norm(du) < 0.5
        @debug "norm(du) = $(norm(du)) at u0 = $u0 u1 = $u1"
        return 1.0
    end
    # Normalize the direction from here
    magd = norm(du)
    if magd < 0.92 && magd > 0
        throw("unexpected magd = $magd (temporarily at least) at u1 = $(u1)")
    end
    # Dot product    
    Δu[1] * du[1] + Δu[2] * du[2]
end

function add_discrete_callbacks!(vdcb, saxy::SelectedVec2AtXY)
    # Flip selected direction callback
    push!(vdcb, DiscreteCallback(condition_flip_bidirection, affect_flip_bidirection!, save_positions=(true,true)))
    # Flip major <--> minor
    push!(vdcb, DiscreteCallback(condition_swap_major_minor, affect_swap_major_minor!, save_positions=(true,true)))
end
add_discrete_callbacks!(vdcb, fxy) = vdcb 
function solve_ensemble(saxy, vu0, tspan, cbs; odekws...)
    throw("Nah, this about intitial con")
    u0 = MVector{2, Float64}(first(vu0))
    function prob_func(prob, i, repeat)
        reset!(prob.p)
        u0 = MVector{2, Float64}(vu0[i])
        remake(prob, u0 = u0)
    end
    prob = ODEProblem(rhs!, u0, tspan, saxy, callback = cbs)
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    # EnsembleThreads: Early results indicated tests varying between runs.
    # We changed to EnsembleSerial. This slows down the calculation.
    solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories = length(vu0); odekws...)
end


"""
    vu0_from_pts(baxy::BidirectionAtXY, pts)
    vu0_from_pts(saxy::SelectedVec2AtXY, pts)
    ---> 
"""
vu0_from_pts(fxy::T, pts) where T<:Journey = vu0_from_pts(fxy.negy, pts)
vu0_from_pts(saxy::SelectedVec2AtXY, pts) = vu0_from_pts(saxy.baxy.negy, pts)
function vu0_from_pts(negy::NegateY, pts)
    @assert eltype(pts) <: CartesianIndex{2}
    vx = map(pt -> float(pt.I[2]), pts)
    vy = map(pt -> negy(float(pt.I[1])), pts)
    [MVector{2, Float64}([x, y]) for (x, y) in zip(vx, vy)]
end



function get_solution(fxy::T, vu0; odekws...) where T<:Journey
    # We need to define a stopping point. 
    # Take it from keywords if supplied. 
    tspan = make_tspan(;odekws...)
    cbs = callbacks_journey(fxy; odekws...)
    # Drop the 'already spent' keywords comprising 'tspan'.
    # The remaining keywords will be passed on to the solver.
    remaining_kws = filter(odekws) do (kw, kwval)
        kw == :tstart && return false
        kw == :tstop && return false
        kw == :dtfloor && return false
        true
    end
    solve_ensemble(fxy, vu0, tspan, cbs; remaining_kws...)
end
=#
