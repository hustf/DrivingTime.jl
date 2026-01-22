# Also see `exported.jl`. These are callees.

function extract_from_solution(sol, ts)
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

function journey_plots(ts, vxaxis, vs, acs, ss; 
    plims = (minimum(vxaxis), maximum(vxaxis)),
    vlims = (minimum(vs), maximum(vs)),
    tlims = (minimum(ts), maximum(ts)),
    alims = (min(zero(acs[1]), minimum(acs)), max(zero(acs[1]), maximum(acs))),
    slims = (min(zero(ss[1]), minimum(ss)), max(zero(ss[1]), maximum(ss))),
    p_vp = plot(),
    p_ap = plot(),
    p_sp = plot(),
    p_tp = plot(),
    plotkws...)
    #
    if eltype(vxaxis) <: Unitful.Time
        xguide = "t"
    else
        xguide = "p"
    end
    #
    plot!(p_vp, vxaxis, vs, xguide = xguide, yguide ="v", label = false, xlims = plims, ylims = vlims; plotkws...)
    areaplot!(p_ap, vxaxis, acs, xguide = xguide, yguide ="a", label = false, xlims = plims, ylims = alims;
          plotkws...)
    areaplot!(p_sp, vxaxis, ss, xguide = xguide, yguide ="s", label = false, xlims = plims, ylims = slims; plotkws...)
    plot!(p_tp, vxaxis, ts, xguide = xguide, yguide ="t", label = false, xlims = plims, ylims = tlims; plotkws...)
    #
    p_vp, p_ap, p_sp, p_tp
end

"""
    time_range_distributed_evenly_along_progression(sol, p_min, p_max, length)
    ---> Vector{Quantity{<:Time}}

Progress `p` must be monotoneously rising with time, and
we're ok with a little inexactness.
"""
function time_range_distributed_evenly_along_progression(sol, p_min, p_max, length; resolution = 0.1u"s")
    @assert p_min < p_max
    # A largish time vector, this is what we will pick times from
    td = collect(range(sol.t[1], sol.t[end], step = resolution))
    n = size(td,1) # because 'length' is a named parameter in this scope.
    # Output pre-allocs, values will be overwritten below
    ts = zeros(typeof(1.0u"s"), length)
    # The evenly distributed progress values
    ps = range(p_min, p_max; length)
    # This is a similarly large vector of rising progress values
    pd = [sol(t; idxs=1) for t in td]
    # 
    @inbounds for (k, p) in pairs(ps) 
        # k is index 1:length
        # p is the current progress value 
        # i is the index for a near match from pd
        i = searchsortedfirst(pd, p)
        # clamp to valid range
        if i <= 1
            ts[k] = td[1]
        elseif i > n
            ts[k] = td[end]
        else
            # We're sloppy here - doesn't matter for a debug plot.
            # We might interpolate, but it's easy to decrease `resolution`
            ts[k] = td[i]
        end
    end
    ts
end


function extract_acceleration_contributions(𝐣::Journey, p, p´)
    @assert length(p) == length(p´)
    motoracc = map(𝐣.fmotoracclim, p´)
    slopeacc = map(𝐣.fslopeacc, p)
    airacc = map(𝐣.fairacc, p´)
    rollacc = map(x -> 𝐣.frollacc(), p)
    motoracc, slopeacc, airacc, rollacc 
end

function plot_acceleration_components!(p_ap, vxax, motoracc, slopeacc, airacc, rollacc)
    c1 = :red
    c2 = :green
    c3 = :blue
    c4 = :black
    plot!(p_ap, vxax, motoracc, label = false, color = c1)
    areaplot!(p_ap, vxax, slopeacc, label = false, color = c2, fillalpha = 0.3)
    plot!(p_ap, vxax, airacc, label = false, color = c3)
    plot!(p_ap, vxax, rollacc, label = false, color = c4)
    xmi = vxax[1] + 0.25 * (vxax[end] - vxax[1])
    xma = vxax[end] * 0.75
    xra = range(xmi, xma, length = 4)
    i1 = searchsortedfirst(vxax, xra[1])
    i2 = searchsortedfirst(vxax, xra[2])
    i3 = searchsortedfirst(vxax, xra[3])
    i4 = searchsortedfirst(vxax, xra[4])
    y1 = motoracc[i1]
    y2 = slopeacc[i2]
    y3 = airacc[i3]
    y4 = rollacc[i4]
    t1 = text("Motor", 14, :hcenter, :vcenter, c1, rotation = 30)
    t2 = text("Slope", 14, :hcenter, :vcenter, c2, rotation = 15)
    t3 = text("Air", 14, :hcenter, :vcenter, c3, rotation = 15)
    t4 = text("Roll", 14, :hcenter, :vcenter, c4, rotation = 0)
    annotate!(p_ap, [(xra[1],  y1, t1),
                  (xra[2],  y2, t2),
                  (xra[3],  y3, t3),
                  (xra[4],  y4, t4)],
                  )
end

function extract_velocity_set_point_and_deviation(𝐣::Journey, p, p´)
    @assert length(p) == length(p´)
    vsetpoint = map(𝐣.itp_v, p)
    vdeviation = vsetpoint .- p´
    vsetpoint, vdeviation  
end


function plot_velocity_set_point_and_deviation!(p_vp, vxax, vsetpoint, vdeviation, p´)
    c1 = :red
    c2 = :green
    plot!(p_vp, vxax, vsetpoint, label = false, color = c1)
    areaplot!(p_vp, vxax, vdeviation, label = false, color = c2, fillalpha = 0.3)
    xmi = vxax[1] + 0.25 * (vxax[end] - vxax[1])
    xma = vxax[end] * 0.75
    xra = range(xmi, xma, length = 2)
    i1 = searchsortedfirst(vxax, xra[1])
    i2 = searchsortedfirst(vxax, xra[2])
    y1 = vsetpoint[i1]
    y2 = vdeviation[i2]
    t1 = text("Setpoint", 14, :hcenter, :vcenter, c1, rotation = 30)
    t2 = text("Deviation", 14, :hcenter, :vcenter, c2, rotation = 30)
    miy = Tvel(minimum(vdeviation)) / 1u"km/hr"
    may = Tvel(max(maximum(p´), maximum(vsetpoint))) / 1u"km/hr"
    ylims!(p_vp, (miy, may))
    annotate!(p_vp, [(xra[1],  y1, t1),
                  (xra[2],  y2, t2)])
end




