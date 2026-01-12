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
    # Interpolate the slope
    ss = [slope_angle(sol.prob.p, p) for  p in ps]
    ts, ps, vs, acs, ss
end

function journey_plots(ts, ps, vs, acs, ss; 
    plims = (minimum(ps), maximum(ps)),
    vlims = (minimum(vs), maximum(vs)),
    tlims = (minimum(ts), maximum(ts)),
    alims = (min(zero(acs[1]), minimum(acs)), max(zero(acs[1]), maximum(acs))),
    slims = (min(zero(ss[1]), minimum(ss)), max(zero(ss[1]), maximum(ss))),
    length = 300, 
    p_vp = plot(),
    p_ap = plot(),
    p_sp = plot(),
    p_tp = plot(),
    plotkws...)
    #
    plot!(p_vp, ps, vs, xguide="p", yguide ="v", label = false, xlims = plims, ylims = vlims; plotkws...)
    areaplot!(p_ap, ps, acs, xguide="p", yguide ="a", label = false, xlims = plims, ylims = alims;
          plotkws...)
    areaplot!(p_sp, ps, ss, xguide="p", yguide ="s", label = false, xlims = plims, ylims = slims; plotkws...)
    plot!(p_tp, ps, ts, xguide="p", yguide ="t", label = false, xlims = plims, ylims = tlims; plotkws...)
    #
    p_vp, p_ap, p_sp, p_tp
end
function plot_journey(sol::SciMLBase.ODESolution; length = 300, tit = "", kws...)
    #
    ts, ps, vs, acs, ss = extract_from_solution(sol; length)
    t = Unitful.minute.(ts)
    p = Unitful.km.(ps)
    v = u"km/hr".(vs)
    p_vp, p_ap, p_sp, p_tp = journey_plots(t, p, v, acs, ss; kws...)
    pl = plot(layout = (4, 1), p_vp, p_ap, p_sp, p_tp)
    if tit !==""
        title!(pl[1], tit)
    end
    pl
end
