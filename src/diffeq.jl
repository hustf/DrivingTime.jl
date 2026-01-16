##############################################################
# Prepare and solve the differential equations for streamlines
##############################################################

function solve_guarded(rhs, Γᵢₙ, 𝐣::Journey, tspan, cbs; debug=false, odekws...)
    s = Unitful.s
    # Test the functions here. The internal __init will do this 
    # as well, but this keeps stacktraces shallow and checks
    # somewhat more strict.
    p, p´ = @inferred packout(Γᵢₙ)
    # Throw-away "derivative" w.r. time
    Γ´ = 0.01Γᵢₙ / s
    if debug
        with_logger(Logging.ConsoleLogger(stderr, Logging.Debug)) do
            x´, x´´ = packout(Γ´)
            @debug x´ x´´
        end
    end
    @inferred packin!(Γ´, 1.0u"m/s", 2.0u"m/s^2")
    @inferred 𝐣.fslopeacc(p)
    @inferred 𝐣.fairacc(p´)
    @inferred 𝐣.fmotoracclim(p´)
    @inferred 𝐣.frollacc()
    @inferred rhs(Γ´, Γᵢₙ, 𝐣, nothing)
    # Define 
    prob = ODEProblem(rhs, Γᵢₙ, tspan, 𝐣, callback = cbs)
    integrator = init(prob, Tsit5(); odekws...)
    # More inferrence checking
    @inferred condition_toofast(Γᵢₙ, 1.0u"s", integrator)
    @inferred condition_reversing(Γᵢₙ, 1.0u"s", integrator)
    sol = if debug
        with_logger(Logging.ConsoleLogger(stderr, Logging.Debug)) do
            solve!(integrator)
        end
    else
        solve!(integrator)
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
    vccb = ContinuousCallback[]
    vdcb = DiscreteCallback[]
    #
    if :dtfloor ∈ keys(odekws)
        push!(vdcb, let dtfloor = odekws[:dtfloor]
                DiscreteCallback(
                   (u, t, integrator ) -> get_proposed_dt(integrator) ≤ dtfloor, 
                       affect!)
            end)
    end
    # Stopped vehicle, we're not interested in reversing
    push!(vdcb, DiscreteCallback(condition_reversing, affect_reversing!, save_positions=(true,true)))
    # Too fast vehicle, we're not interested in results for a bad model
    push!(vdcb, DiscreteCallback(condition_toofast, affect_toofast!, save_positions=(true,true)))
    # End of the journey geometry reached. Extrapolated ends makes 
    # overshooting the geometry not problematic. Capture pstop:
    pstop = 𝐣.pstop
    function condition_endprogress(u, t, integrator::ODEIntegrator)
        p, p´ = packout(u)
        p >= pstop
    end
    push!(vdcb, DiscreteCallback(condition_endprogress, affect_endprogress!, save_positions=(true,true)))
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
    p, p´ = packout(u)
    @debug "u"    p    p´   maxlog = 2
    p´´1 = 𝐣.fslopeacc(p)
    p´´2 = 𝐣.fairacc(p´)
    p´´3 = 𝐣.fmotoracclim(p´)
    p´´4 = 𝐣.frollacc()
    @debug "p´´"  p´´1   p´´2   p´´3  p´´4 maxlog = 2
    if p´ < 80u"km/hr" # Temp of course
        p´´ = p´´1  + p´´2 + p´´3  + p´´4
    else
        p´´ = p´´1  + p´´2 + p´´4
    end
    @debug "du "    p´´     maxlog = 2
    packin!(du, p´, p´´)
    du
end

# Our own function for potential debugging termination causes.
function affect!(integrator)
    @debug "Terminate"
    terminate!(integrator)
end

function condition_reversing(u, t, integrator::ODEIntegrator)
    @assert t > integrator.tprev
    p, p´ = packout(u)
    p > zero(p) && p´ < zero(p´)
end

function affect_reversing!(integrator)
    # For debugging
    @debug "Terminate due to reversing"
    terminate!(integrator)
end

function condition_toofast(u, t, integrator::ODEIntegrator)
    p, p´ = packout(u)
    p´ > 350.0u"km/hr"
end

function affect_toofast!(integrator)
    # For debugging
    @debug "Terminate due to driving too fast"
    terminate!(integrator)
end

function affect_endprogress!(integrator)
    # For debugging
    @debug "Terminate due to reaching end of journey"
    terminate!(integrator)
end
