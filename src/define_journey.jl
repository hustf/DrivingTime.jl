# This constructs a Journey, which will be the integrator in 
# a differential equation. All the data pertinent to the calculation
# is stored as integrator fields, and also 


"""
    locally_smooth_interpolation(p::T, s) where T
    ---> interpolator(::T --> typeof(s)

For varying distance `p`. `p` must be rising. 
"""
function locally_smooth_interpolation(p::Vector{T}, s) where T
    @assert length(p) == length(s)
    # Equally-spaced progression points
    pgrid = range(p[1], p[end]; length = Int(round(p[end] / oneunit(T))))
    # An iterpolation functor for slope onto grid 
    itp_sgrid = extrapolate(interpolate((p, ), s, Gridded(Linear())), Throw())
    # Equally-spaced slope points
    sgrid = itp_sgrid.(pgrid)
    # Interpolation functor, extended tangentially past both ends
    extrapolate(scale(interpolate(sgrid, BSpline(Cubic(Line(OnGrid())))), pgrid), Line())
end

"""
    max_speed_adjusted_for_curvature(v, r)
    ---> Vector{Quantity{<:Velocity}}

Measurements from EcoSafe analyzed in `dev_acceptable_lateral_acceleration`. 
The maximum lateral acceleration should be 

    2.2 m/s² = a = v² / abs(r)
 => 
    v <= √(abs(r) * 2.2 m/s^2) 

where `r` is signed road radius of curvature, `a` is centripetal acceleration and 
`v` is the maxium speed for comfortable passenger transport.

This one is hard coded since it's well founded.
"""
function max_speed_adjusted_for_curvature(v, r)
    replace!(r, NaN * u"m"  => 0.0u"m")
    @assert length(v) == length(r)
    # Wanted velocity
    w = copy(v)
    for i in eachindex(v)
        if r[i] !== 0.0u"m"
            w =  √(abs(r[i]) * 2.2u"m/s^2")
            if v[i] > w
                v[i] = w
            end
        end
    end
    v
end

"""
    Journey((ea1, no1, ea2, no2; 
    default_fartsgrense = 50, 
    f_air_acc = AirAcceleration(),
    f_motor_acclim = MotorlimAcceleration(),
    f_roll_acc = RollRAcceleration()

An external constructor for the ODE right-hand side parameters.

# Arguments

- `ea1, no1, ea2, no2` Utm33 easting and northing coordinates. Start at 1, end at 2.
- `default_fartsgrense` is used in case the starting point has no defined speed limit, e.g. in bus terminals.
"""
function Journey(ea1, no1, ea2, no2; 
    default_fartsgrense = 50, 
    f_air_acc = AirAcceleration(),
    f_motor_acclim = MotorlimAcceleration(),
    f_roll_acc = RollRAcceleration())
    #
    # Slope and progression (and fartsgrense)
    d = route_leg_data(ea1, no1, ea2, no2; default_fartsgrense)
    # Gravity's component along the surface. Positive slope => negative component
    s = map(gravity_comp_along_surface, d[:slope])
    # Progression
    p = d[:progression] * u"m"
    # Interpolator progression -> slope
    itp_s = locally_smooth_interpolation(p, s)
    # End of journey
    pstop = Unitful.km(p[end])
    # Wanted speed
    v = max_speed_adjusted_for_curvature(d[:speed_limitation] * u"km/hr", d[:radius_of_curvature] * u"m")
    # v may have abrupt changes, but our simple control function doesn't look ahead.
    # Instead, we reduce the velocity gradient ahead of step changes downward.
    # TODO
    # Interpolator progression -> wanted speed
    itp_v = locally_smooth_interpolation(p, v)
    # Construct
    Journey(pstop, itp_s, f_air_acc, f_motor_acclim, f_roll_acc, itp_v)
end

# Utility
"""
    slope_angle(𝐣::Journey, p::T) where T<:Length
    --> Quantity{°}

We give the angle as a "unit", because the plotting recipe for 
quantities would be confused by mixing in a plain Float64.
""" 
function slope_angle(𝐣::Journey, p::T) where T<:Length
    # a = -g * sin(atan(slope))
    # <=>
    # a/-g = sin(atan(slope))
    # =>
    # asin(a/-g) = atan(slope) = slope_angle
    a = 𝐣.fslopeacc(p)
    Unitful.°(asin(a/-g))
end
