# This constructs a Journey, which will be the integrator in 
# a differential equation. All the data pertinent to the calculation
# is stored as integrator fields, and also 

"""
    Journey(ea1, no1, ea2, no2; ; default_fartsgrense = 50)

An external constructor for the right-hand side functor for 
this ordinary differential equation.

# Arguments

- `ea1, no1, ea2, no2` Utm33 easting and northing coordinates. Start at 1, end at 2.
- `default_fartsgrense` is used in case the starting point has no defined speed limit, e.g. in bus terminals.
"""
function Journey(ea1, no1, ea2, no2; 
    default_fartsgrense = 50, 
    mass = VEHICLE_DEFAULTS.mass,
    power = VEHICLE_DEFAULTS.power,
    motorlim = VEHICLE_DEFAULTS.motorlim,
    frontarea = VEHICLE_DEFAULTS.frontarea,
    shapecoeff = VEHICLE_DEFAULTS.shapecoeff,
    ρ = ρ_humid(ENVIRONMENT_DEFAULTS.T))
    # Acceleration from air resistance, function of velocity
    f_air_acc = let ρ = ρ; A = frontarea; Cs = shapecoeff; m = mass
        v -> -wind_force(Cs, A, ρ, v) / m
    end
    # Motor limit, function of velocity.
    # Rotary masses are neglected (no gearing for electric vehicle)
    f_motor_acclim = let P = power; Flim = motorlim; m = mass
        v -> max(0.0u"m/s^2", min(P / v, Flim) / m)
    end
    # Slope and progression
    d = route_leg_data(ea1, no1, ea2, no2; default_fartsgrense)
    p = d[:progression]
    s = d[:slope]
    # Equally-spaced progression points
    pgrid = range(0.0, p[end]; length = Int(round(p[end])) )
    # An iterpolation functor for slope onto grid 
    itp_sgrid = extrapolate(interpolate((p, ), s, Gridded(Linear())), Throw())
    # Equally-spaced slope points
    sgrid = itp_sgrid.(pgrid)
    # An interpolation functor for slope, extended tangentially past both ends
    itp_s = extrapolate(scale(interpolate(sgrid, BSpline(Cubic(Line(OnGrid())))), pgrid), Line())
    # Slope acceleration, function of progression
    f_slope_acc = let g = g; itp_s = itp_s
        p -> let
            p_m = NoUnits(p / Unitful.m)
            # Consider TODO:
            # Of course, we should have stored this as the extrapolation...
            comp_along_progression = sin(atan(itp_s(p_m)))
            -g * comp_along_progression
        end
    end
    Journey(itp_s, f_slope_acc, f_air_acc, f_motor_acclim)
end

# Utility
"""
    slope_angle(𝐣::Journey, p::T) where T<:Length
    --> Quantity{°}

We leave in this "unit" because the plotting recipe for quantities would be confused by
a plain Float64.
""" 
function slope_angle(𝐣::Journey, p::T) where T<:Length
    p_m = NoUnits(p / Unitful.m)
    Unitful.°(atan(𝐣.fslope(p_m)))
end
