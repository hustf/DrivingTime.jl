using Test
using DrivingTime
using RouteSlopeDistance
using Plots
import Interpolations
using Interpolations: extrapolate, interpolate, Gridded, Linear, Line, Cubic, OnGrid, BSpline, scale, gradient
using Interpolations: SteffenMonotonicInterpolation, FritschButlandMonotonicInterpolation, LinearMonotonicInterpolation
using Interpolations: Flat, Tangent, Throw, linear_interpolation, cubic_spline_interpolation
na1 = "Eiksund"
ea1 = 27978
no1 = 6935574
na2 = "Sørestranda"
no2 = 6928773
ea2 = 33196
d = route_leg_data(ea1, no1, ea2, no2)
x, y, z = unique_unnested_coordinates_of_multiline_string(d[:multi_linestring])
p = d[:progression]
s = d[:slope]
@test length(p) == length(x) == length(s) == length(unique(p))
@test p[1] == 0.0
function plotpath(x, y, z, p)
   pl = plot(x, y,  
       aspect_ratio = :equal,
       widen = 1.2, 
       ribbon = (-z, 0 .* z), 
       markersize = 2, markershape = :xcross,
       label = "east,north,(z)")
   zmin, izmin = findmin(z)
   zmax, izmax = findmax(z)
   t1 = text("0.0 km", 14, :left, :top, :green)
   t2 = text("$(round(p[end] / 1000; digits = 1)) km", 14, :hcenter, :vcenter, :green)
   t3 = text("\n\n\nz $(round(zmin; digits = 1)) m", 10, :hcenter, :vcenter, :blue, rotation=45)
   t4 = text("z $(round(zmax; digits = 1)) m\n\n\n", 10, :hcenter, :vcenter, :blue, rotation=45)
   annotate!(pl, [(x[1], y[1], t1),
                  (x[end], y[end], t2),
                  (x[izmin], y[izmin], t3),
                  (x[izmax], y[izmax], t4)],
                  )
   
end
plotpath(x, y, z, p)
rng = 1:60
plotpath(x[rng], y[rng], z[rng], p[rng])
let 
    pl = plot(p[1:10], label = false)
    xlabel!(pl, "Point no.")
    ylabel!(pl, "Progression [m]")
    pl
end

pl = let 
    rn = 1:72
    pl = plot(p[rn], s[rn], label = false, seriestype=:scatter, widen = 1.1)
    xlabel!(pl, "Progression [m]")
    ylabel!(pl, "Slope")
    pl
end

#######################################
# Define knots for smooth interpolation
#######################################
# A shorthand for 
# extrapolate(interpolate((p, ), s, Gridded(Linear())), Throw())
itp_sgrid = linear_interpolation(p, s)
itp_sgrid1 = extrapolate(interpolate((p, ), s, Gridded(Linear())), Throw())


itp_sgrid(0.0)
itp_sgrid(p[end])
itp_sgrid(1500)
plot!(pl, x->itp_sgrid(x), 0.0:1.0:1750, label = false, seriestype=:line)
# Equally-spaced progression points
pgrid = range(0.0, p[end]; length = Int(round(p[end])) )
# Equally-spaced slope points
sgrid = itp_sgrid.(pgrid)
plot(pgrid, sgrid, seriestype=:scatter, markersize=2)
pl = let 
    rn = 1:72
    pl = plot(pgrid[rn], sgrid[rn], label = false, seriestype=:scatter, widen = 1.1)
    xlabel!(pl, "Progression [m]")
    ylabel!(pl, "Slope")
    pl
end

#############################
# Define smooth interpolation
#############################
#
# We believe we need this step for effectively
# solving the differential equations.
# We want e.g. slope(progression) to be continuous
# and once-differentiable. Of course, we placed the 
# grid on straight lines betweeen given points, 
# so the 'smooth interpolation' is only smooth 
# in the differentiable respect.

itp_s = cubic_spline_interpolation(pgrid, sgrid; extrapolation_bc = Line())
itp_s1 = extrapolate(  
            scale(
                interpolate(sgrid, BSpline(Cubic(Line(OnGrid())))), 
                pgrid),
            Line())
@test typeof(itp_s) == typeof(itp_s1)

itp_s2 = extrapolate(
            interpolate(sgrid, BSpline(Cubic(Line(OnGrid())))),
            Line())
itp_s(-1)
itp_s(20000)
pl = plot(x->itp_s(x), xlims = (-1000, 440), label = "itp_s")
plot!(pl, x->itp_s1(x), label = "itp_s1")
plot!(pl, x->itp_s2(x), label = "itp_s2")
plot!(pl, p, s, seriestype = :scatter, label = "point")
plot!(pl, xlims = (150,180), ylims = (-0.02,0.002))

pl = plot(x->itp_s1(x), xlims = (166, 177), ylims = (0.001, 0.002), label = "Smooth slope", xlabel = "Progression [m]")
plot!(pl, p, s, seriestype = :scatter, label = "point")
pl = plot(x->itp_s1(x), xlims = (159, 179), ylims = (0.00, 0.002), label = "itp_s1", xlabel = "Progression [m]")
plot!(pl, x->itp_s2(x), label = "itp_s2")
plot!(pl, p, s, seriestype = :scatter, label = "point")
itp_s1(179)

# We could drop `scale` to improve speed, as in itp_s2. But not worth doing (yet).