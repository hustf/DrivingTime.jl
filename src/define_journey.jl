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
    reduce_speed_looking_ahead!(v, p, i; consider_dist = 500.0u"m", acc_constant = -1.0u"m/s^2")
    reduce_speed_looking_ahead!(v, p; consider_dist = 500.0u"m", acc_constant = -1.0u"m/s^2")
    ---> Vector{Tvel}

This could have been done simpler by looping from end of `v` to start, all
the while keeping track of the last maximum, reducing it as we go. But remember, 
progress 'p' density varies widely, and time is unknown. So maybe this is not that bad anyway.
"""
function reduce_speed_looking_ahead!(v, p, i; consider_dist = 500.0u"m", acc_constant = -1.0u"m/s^2")
    #     
    n = length(v)
    # We're looking at points ahead for which we need to 
    # brake in advance. But there's no point in looking 
    # past our stopping distance at constant acceleration: s = v1² / (2a) 
    Δp_at_potential_stop = Unitful.m(v[i]^2 / (2 * -acc_constant))
    # Progress increase monotoneously, so it's easy to identify the 
    # last interesting index:
    coming_progress = view(p, (i + 1):n)
    Δi_after_potential_stop = searchsortedfirst(coming_progress, p[i] + Δp_at_potential_stop)
    # Function of indices i + 1, i + 2 etc.
    v_achieveable(i_ahead::Int64) = let a = acc_constant; v0² = v[i]^2; p0 = p[i]; p = p
        Δp = p[i_ahead] - p0
        @assert Δp < Δp_at_potential_stop "i_ahead = $i_ahead Δp = $Δp < Δp_at_potential_stop = $(Δp_at_potential_stop )"
        @assert v0² >= 2 * a * Δp 
        sqrt(v0² + 2 * a * Δp)
    end
    # Which point ahead is the one determining reduction in planned velocity here?
    i_determine = argmax((i + 1):(i + Δi_after_potential_stop - 1)) do i_ahead
        if v[i_ahead] < v[i]
            v_ach = v_achieveable(i_ahead)
            if v[i_ahead] < v_ach 
                v_ach - v[i_ahead]
            else
                0.0u"m/s"
            end 
        else 
            0.0u"m/s"
        end
    end
    Δv_positive = v_achieveable(i_determine) - v[i_determine]
    if Δv_positive > zero(Δv_positive)
        v[i] -= Δv_positive
    else
        # We can easily reach the speed needed at p[i_determine] without 
        # changing the speed (set point) at p[i]
    end
    v
end

function reduce_speed_looking_ahead!(v, p; consider_dist = 500.0u"m", acc_constant = -1.0u"m/s^2")
    @assert length(v) == length(p)
    @assert acc_constant < zero(acc_constant)
    # We're assuming constant acceleration
    # v1² = 2as + v0²
    for i in 1:(length(v) - 1)
        reduce_speed_looking_ahead!(v, p, i; consider_dist, acc_constant)
    end
    v
end


"""
    Journey(waypoints; 
        default_fartsgrense = 50, 
        f_air_acc = AirAcceleration(),
        f_motor_acclim = MotorlimAcceleration(),
        f_roll_acc = RollRAcceleration()

An external constructor for the ODE right-hand side parameters.

# Arguments

- `waypoints` Utm33 easting and northing coordinates, mininum two pairs. See `continuous_route_data`.
- `default_fartsgrense` is used in case the starting point has no defined speed limit, e.g. in bus terminals.
"""
function Journey(waypoints; 
    default_fartsgrense = 50, 
    f_air_acc = AirAcceleration(),
    f_motor_acclim = MotorlimAcceleration(),
    f_roll_acc = RollRAcceleration())
    #
    # Slope and progression (and fartsgrense)
    d = continuous_route_data(waypoints; default_fartsgrense)
    Journey(d; f_air_acc, f_motor_acclim, f_roll_acc)
end

function Journey(d::Dict{Symbol, Any}; 
    f_air_acc = AirAcceleration(),
    f_motor_acclim = MotorlimAcceleration(),
    f_roll_acc = RollRAcceleration())
    # Gravity's component along the surface. Positive slope => negative component
    s = map(gravity_comp_along_surface, d[:slope])
    # Progression
    p = d[:progression] * u"m"
    # Interpolator progression -> slope
    itp_s = locally_smooth_interpolation(p, s)
    # End of journey
    pstop = Unitful.km(p[end])
    # Speed limit with reductions based on curvature and speed bumps
    v = max_speed_adjusted_for_curvature(d[:speed_limitation] * u"km/hr", d[:radius_of_curvature] * u"m")
    # Additionaly, we're going to end the journey with speed ~0 km/hr.
    v[end] = 1.0u"km/hr"
    # v may have abrupt changes, but our simple motor control function doesn't look ahead.
    # Instead, we reduce the 'wanted velocity' ahead of step changes downward.
    reduce_speed_looking_ahead!(v, p)
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


