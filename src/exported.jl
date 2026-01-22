"""
    drivetime(ea1, no1, ea2, no2)
    --> Dates.Minute

Starting point at Utm33 coordinates ("easting, northing") `ea1`, `no1`. Starting velocity: 0 km/hr.

End point at Utm33 coordinates ("easting, northing") `ea1`, `no1`. Starting velocity: 1 km/hr.

The route, road geometry and speed limits are taken from Norsk Vegdatabase's API if not already cached locally. 
Driving time is found by integration of a simple model involving:

 - distance (of course)
 - speed limits (slow down in advance of reduced speed limit)
 - slope (affects acceleration and sometimes max speed)
 - curvature (limits speed through a centripetal acceleration limit)
 - speed bumps (limits velocity to 15 km/hr below the local speed limit)
 - unofficial vehicle data for an electric intercity bus, Yutong IC12E
 - temperature, air pressure and resistance
 - roll resistance

Prior to starting the integration, the model parameters are collected in a `Journey` object. For detail control
or for output resolution, see:

  - `Journey`
  - `solve_journey`
  - `plot_journey`
  - `slope`
  - `@u_str`

# Example
```
julia> using DrivingTime

julia> drivetime(24062, 6939037, 28592, 6939504)
8 minutes
```

"""
function drivetime(ea1, no1, ea2, no2)
    sol = solve_journey(Journey(ea1, no1, ea2, no2))
    t = sol.t[end]
    # Unitful.Time to Dates.Period
    Minute(round(round(t / u"minute")))
end



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
