"""
    drivetime(ea1, no1, ea2, no2; default_fartsgrense = 50)
    drivetime(waypoints; default_fartsgrense = 50)
    drivetime(j::Journey)
    --> Dates.CompoundPeriod

Simulated time to drive a heavy vehicle the specified route

  - Starting velocity: 0 km/hr
  - End velocity: 1 km/hr

See `continuous_route_data` regarding waypoints.

The route, road geometry and speed limits are taken from Norsk Vegdatabase's API if not already cached locally. See `RouteSlopeDistance.jl`.

Prior to starting the integration, the model parameters are collected in a `Journey` object. For detail control
or for output resolution, see:

  - `Journey`
  - `solve_journey`
  - `plot_journey`
  - `slope`
  - `@u_str`

Driving time is found by time integration of a simple model involving:

 - distance (of course)
 - speed limits (slow down in advance of reduced speed limit)
 - slope (affects acceleration and sometimes max speed)
 - curvature (limits speed through a centripetal acceleration limit)
 - speed bumps (limits velocity to 15 km/hr below the local speed limit)
 - unofficial vehicle data for an electric intercity bus, Yutong IC12E
 - temperature, air pressure and resistance
 - roll resistance


# Example
```
julia> using DrivingTime

julia> drivetime(24062, 6939037, 28592, 6939504)
7 minutes, 48 seconds
```

"""
drivetime(ea1, no1, ea2, no2; default_fartsgrense = 50) = drivetime([(ea1, no1), (ea2, no2)]; default_fartsgrense)
drivetime(waypoints; default_fartsgrense = 50) = drivetime(Journey(waypoints; default_fartsgrense))
function drivetime(j::Journey)
    sol = solve_journey(j)
    t = sol.t[end]
    # Unitful.Time to Dates.Period
    nmi = Int(floor(t / u"minute"))
    nse = Int(floor((t / u"minute" - nmi) * 60))
    Dates.Minute(nmi) + Dates.Second(nse)
end


"""
    continuous_route_data(points::Vector{NTuple{2, Int64}}; default_fartsgrense = 50)
    continuous_route_data(fromtos::Vector{NTuple{4, Int64}}; default_fartsgrense = 50)
    continuous_route_data(fromto::NTuple{4, Int64}; default_fartsgrense = 50)
    ---> Dict{Symbol, Any}

Calls to 'RouteSlopeDistance.route_leg_data' between A and C will cache the result, while this will cache legs individually.

If several routes between `A` and `C` are of interest, caching one of those is not a good idea.

`continuous_route_data` will cache `dAB` and `dBC`, but not the returned `dABC`. Multiple waypoints OK.
"""
function continuous_route_data(points::Vector{NTuple{2, Int64}}; default_fartsgrense = 50)
    @assert length(points) > 1
    fromtos = Vector{NTuple{4, Int64}}()
    for i in eachindex(points)
        i == 1 && continue
        fromto = (points[i - 1]..., points[i]...)
        push!(fromtos, fromto)
    end
    continuous_route_data(fromtos; default_fartsgrense)
end
function continuous_route_data(fromtos::Vector{NTuple{4, Int64}}; default_fartsgrense = 50)
    @assert ! isempty(fromtos)
    # Check continuity
    for i in eachindex(fromtos)
        i == 1 && continue
        if fromtos[i - 1][3:4] !== fromtos[i][1:2]
            throw(ArgumentError("Expected start position of $i to match end position of $(i - 1): $(fromtos[i-1:i])"))
        end
    end
    # Join 
    d = route_leg_data(fromtos[1]...; default_fartsgrense)
    for i in 2:length(fromtos)
        dleg = route_leg_data(fromtos[i]...; default_fartsgrense)
        d = join_route_data(d, dleg) 
    end
    d
end
continuous_route_data(fromto::NTuple{4, Int64}; default_fartsgrense = 50) = continuous_route_data([fromto], default_fartsgrense)

"""
    plot_journey(sol::SciMLBase.ODESolution; xtime::Bool = false, 
        length = 300, tit = "", progress_max = nothing, kws...)
    ---> Plot

For checking. Full visualization belongs elsewhere.
 
Note here that the x-axis representing progress (position)
or time is selected with keyword argument `xtime`.
"""
function plot_journey(sol::SciMLBase.ODESolution; xtime::Bool = false, 
    length = 300, tit = "", progress_max = nothing, kws...)
    #
    if xtime
        # Time range distributed evenly along time
        @assert isnothing(progress_max) "progress_max can't be set when xtime is true"
        ts = range(sol.t[1], sol.t[end]; length)
    else
        # Time range distibuted evenly along progress
        if !isnothing(progress_max)
            @assert dimension(progress_max) == dimension(sol.u[1][1])
        else
            progress_max = sol.u[end][1]
        end
        progress_min = sol.u[1][1]
        ts = time_range_distributed_evenly_along_progression(sol, progress_min, progress_max, length)
    end
    # Time, progress, velocity, acceleration, slope angle
    ts, ps, vs, acs, ss = extract_from_solution(sol, ts)
    t = Unitful.minute.(ts)
    p = Unitful.km.(ps)
    v = u"km/hr".(vs)
    if xtime
        vxaxis = t
    else
        vxaxis = p
    end
    p_vp, p_ap, p_sp, p_tp = journey_plots(t, vxaxis, v, acs, ss; kws...)
    # Also plot contributions to acceleration.
    motoracc, slopeacc, airacc, rollacc = extract_acceleration_contributions(sol.prob.p, ps, vs)
    plot_acceleration_components!(p_ap, vxaxis, motoracc, slopeacc, airacc, rollacc )
    # Also plot velocity limit with reductions, and deviation
    vsetpoint, vdeviation = extract_velocity_set_point_and_deviation(sol.prob.p, ps, vs)
    plot_velocity_set_point_and_deviation!(p_vp, vxaxis, vsetpoint, vdeviation, vs)
    #
    # Assemble plots
    pl = plot(layout = (4, 1), p_vp, p_ap, p_sp, p_tp)
    if tit !==""
        title!(pl[1], tit)
    end
    pl
end
