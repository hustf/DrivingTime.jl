using Test
using DrivingTime
using DrivingTime: packin!, packout, rhs!
using Unitful
#using StaticArrays
using DrivingTime: make_tspan, callbacks_journey
import OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, ContinuousCallback, DiscreteCallback, CallbackSet
using OrdinaryDiffEq: terminate!, Tsit5, solve, get_proposed_dt, init, solve!
import OrdinaryDiffEqCore
using OrdinaryDiffEqCore: ODEIntegrator
import SciMLBase
using SciMLBase: ReturnCode.Success, ReturnCode.Terminated, ReturnCode.DtLessThanMin
using SciMLBase: successful_retcode
using ComponentArrays
using RecursiveArrayTools
using Logging
using Plots

na1 = "Eiksund"
ea1 = 27978
no1 = 6935574
na2 = "Sørestranda"
no2 = 6928773
ea2 = 33196
𝐣 = Journey(ea1, no1, ea2, no2);
@test sizeof(𝐣) == 72


Γᵢₙ = ArrayPartition([0.0u"m"], [0.0u"m/s"])
sol = solve_journey(Γᵢₙ, 𝐣);
sol.t[end]

function extract_from_solution(sol; length = 300)
    ts = range(sol.t[1], sol.t[end]; length)
    ps = map(t->sol(t).x[1][1], ts)
    vs = map(t->sol(t).x[2][1], ts)
    # Also retrieve the acceleration
    # The field .du is not available from the solver,
    # but we'll just map solution points to the rhs function
    rhs = sol.prob.f
    Γs = map(sol, ts)
    dummyΓ´ = Γs[1] / (ts[2]- ts[1])
    𝐣 = sol.prob.p
    acs = [rhs(dummyΓ´, Γ, 𝐣, nothing).x[2][1] for  Γ in Γs]
    ts, ps, vs, acs
end

function journey_plots(ts, ps, vs, acs; 
    plims = (minimum(ps), maximum(ps)),
    vlims = (minimum(vs), maximum(vs)),
    tlims = (minimum(ts), maximum(ts)),
    alims = (min(zero(acs[1]), minimum(acs)), max(zero(acs[1]), maximum(acs))),
    length = 300, 
    p_vp = plot(),
    p_ap = plot(),
    p_tp = plot(),
    plotkws...)
    #
    plot!(p_vp, ps, vs, xguide="p", yguide ="v", label = false, xlims = plims, ylims = vlims; plotkws...)
    areaplot!(p_ap, ps, acs, xguide="p", yguide ="a", label = false, xlims = plims, ylims = alims;
          plotkws...)
    plot!(p_tp, ps, ts, xguide="p", yguide ="t", label = false, xlims = plims, ylims = tlims; plotkws...)
    #
    p_vp, p_ap, p_tp
end
function plot_journey(sol; length = 300, kws...)
    #
    ts, ps, vs, acs = extract_from_solution(sol; length)
    t = Unitful.minute.(ts)
    p = Unitful.km.(ps)
    v = u"km/hr".(vs)
    p_vp, p_ap, p_tp = journey_plots(t, p, v, acs; kws...)
    plot(layout = (3, 1), p_vp, p_ap, p_tp, legend=false)
end

plot_journey(sol)
