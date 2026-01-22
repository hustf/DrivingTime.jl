using Test
using DrivingTime
using RouteSlopeDistance
using Plots
na1 = "Møre skule"
ea1 = 24062
no1 = 6939037
na2 = "Ringstaddalen"
ea2 = 28592
no2 = 6939504
tit = na1 * " -> " * na2
sol = solve_journey(Journey(ea1, no1, ea2, no2));
pl = plot_journey(sol; tit)
plot(pl[1], xlim =(5.5u"km", 6u"km"))
d = route_leg_data(ea1, no1, ea2, no2; default_fartsgrense = 50.0u"km")
using RouteSlopeDistance
plot_elevation_and_slope_vs_progression(d, na1, na2)