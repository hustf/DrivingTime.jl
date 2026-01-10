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
function Journey(ea1, no1, ea2, no2; default_fartsgrense = 50)
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
    Journey(itp_s)
end

