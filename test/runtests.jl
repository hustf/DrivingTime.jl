using Test
using DrivingTime
using Dates
na1 = "Møre skule"
ea1 = 24062
no1 = 6939037
na2 = "Ringstaddalen"
ea2 = 28592
no2 = 6939504
tit = na1 * " -> " * na2
@test drivetime(ea1, no1, ea2, no2) == Minute(8)
# A closer look
sol = solve_journey(Journey(ea1, no1, ea2, no2));
pl = plot_journey(sol; tit)
DrivingTime.plot(pl[1], xlim =(5.5u"km", 6u"km"))
# Some 'unitless' plots of the input.
using RouteSlopeDistance
d = route_leg_data(ea1, no1, ea2, no2; default_fartsgrense = 50.0u"km")
plot_elevation_and_slope_vs_progression(d, na1, na2)
plot_elevation_slope_speed_vs_progression(d, na1, na2)
@test true