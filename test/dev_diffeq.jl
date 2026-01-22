using Test
using DrivingTime
using DrivingTime: MotorlimAcceleration, AirAcceleration
using DrivingTime: ENVIRONMENT_DEFAULTS, VEHICLE_DEFAULTS
using DrivingTime: ρ_humid, psat, ρ
using Plots
import Unitful
using Unitful:uconvert
using Unitful: °C

# Velocity check
# Volda- Eiksund, 6°C, 4 passenger: 63 km/hr, stable
# Eiksund - Sørheim, 6°C, 1 passenger: 73 km/hr, stable
# Acceleration check: 10.4s to 50km/hr, near flat, "Fantehuset" to Hareid u.s.

na1 = "Eiksund"
ea1 = 27978
no1 = 6935574
na2 = "Sørestranda"
no2 = 6928773
ea2 = 33196
tit = na1 * " -> " * na2
𝐣 = Journey(ea1, no1, ea2, no2);
@test sizeof(𝐣) == 240



#  13.943960 seconds (70.56 M allocations: 3.432 GiB, 5.13% gc time, 188.06% compilation time: <1% of which was recompilation)
# 0.000584 seconds (1.23 k allocations: 50.570 KiB, 0.21% compilation time: 100% of which was recompilation)
# 0.009327 seconds (123.87 k allocations: 4.592 MiB, 238.85% compilation time: <1% of which was recompilation)
# 0.008081 seconds (95.79 k allocations: 3.181 MiB, 0.02% compilation time: 100% of which was recompilation)
# 0.000983 seconds (6.80 k allocations: 195.695 KiB, 0.14% compilation time: 100% of which was recompilation)
# 0.002936 seconds (30.47 k allocations: 1.013 MiB, 0.04% compilation time: 100% of which was recompilation)
# 0.003486 seconds (33.83 k allocations: 1.624 MiB, 0.04% compilation time: 100% of which was recompilation)
# 0.001245 seconds (1.79 k allocations: 76.148 KiB, 0.13% compilation time: 100% of which was recompilation)
# 0.000722 seconds (1.85 k allocations: 80.523 KiB, 0.22% compilation time: 100% of which was recompilation)
# 0.000638 seconds (4.51 k allocations: 207.805 KiB, 3472.15% compilation time: <1% of which was recompilation)
# 0.000797 seconds (1.98 k allocations: 89.211 KiB, 0.21% compilation time: 100% of which was recompilation)
# 0.001096 seconds (4.92 k allocations: 185.852 KiB, 0.23% compilation time: 100% of which was recompilation)
# 0.001203 seconds (4.92 k allocations: 186.383 KiB, 0.20% compilation time: 100% of which was recompilation)
# 0.002372 seconds (29.66 k allocations: 1.147 MiB, 1213.61% compilation time: <1% of which was recompilation)
# @time solve_journey(𝐣);
sol = solve_journey(𝐣);
sol.t[end] # 453s
pl = plot_journey(sol; tit)
plot(pl[2])
# We have not yet developed a braking feature, so it's more interesting to start elsewhere.
na1 = "Tunnelbotn"
ea1 = 28132
no1 = 6931884
na2 = "Sørheim"
ea2 = 32452
no2 = 6930519
tit = na1 * " -> " * na2

f_motor_acclim = MotorlimAcceleration(;rmp = 0.5u"km/hr", η = 0.9)
𝐣 = Journey(ea1, no1, ea2, no2; f_motor_acclim);
sol = solve_journey(𝐣);
pl = plot_journey(sol; tit, progress_max = 3u"km")
# Velocity check (tested: 73km/hr)
plot(pl[1], minorgrid = true, minorticks = 5, ylims = (60,80))

na1 = "Tunnelbotn før botn"
ea1 = 28687
no1 = 6931679
na1 = "Tunnelbotn"
ea1 = 28132
no1 = 6931884

na2 = "Eiksund"
ea2 = 27978
no2 = 6935574
f_motor_acclim = MotorlimAcceleration(;rmp = 0.5u"km/hr", η = 0.9)
𝐣 = Journey(ea1, no1, ea2, no2; f_motor_acclim);
sol = solve_journey(𝐣);
pl = plot_journey(sol; tit, progress_max = 4u"km")
# Velocity check (tested: 63 km/hr)
plot(pl[1], minorgrid = true, minorticks = 5, ylims = (60,80))


# Start on a flat stretch to calibrate 0-50km/hr in 10.4s.
na1 = "Fantehuset"
ea1 = 35513
no1 = 6947482
na2 = "Barnehage"
ea2 = 36004
no2 = 6947678
tit = na1 * " -> " * na2

vrng =  0.0u"km/hr": 0.1u"km/hr":120.0u"km/hr"
f_motor_acclim = MotorlimAcceleration(;rmp = 0.5u"km/hr", η = 0.9)
pl = plot(vrng, f_motor_acclim.(vrng), ylims = (0, 1.65), xlims=(0, 82))
𝐣 = Journey(ea1, no1, ea2, no2; f_motor_acclim);
sol = solve_journey(𝐣);
t_test = 10.4u"s"
uconvert(u"km/hr", sol(t_test; idxs = 2))



pl = plot_journey(sol; tit, xtime = false)
pl = plot_journey(sol; tit, xtime = true)
pl = plot_journey(sol; tit, xtime = true, xlim = (0, 0.17))

# Multiple velocity limit changes
using DrivingTime
using DrivingTime: route_leg_data, max_speed_adjusted_for_curvature, reduce_speed_looking_ahead!
using Plots
na1 = "Møre skule"
ea1 = 24062
no1 = 6939037
na2 = "Ringstaddalen"
ea2 = 28592
no2 = 6939504
tit = na1 * " -> " * na2
d = route_leg_data(ea1, no1, ea2, no2; default_fartsgrense = 50.0u"km")
p = d[:progression] * u"m"
v = max_speed_adjusted_for_curvature(d[:speed_limitation] * u"km/hr", d[:radius_of_curvature] * u"m")
vold = copy(v)
plot(p, [vold v])
plot(p[100:200], [vold[100:200] v[100:200]])
plot(100:125, [vold[100:125] v[100:125]])
plot(p[100:125], [vold[100:125] v[100:125]])


reduce_speed_looking_ahead!(v, p)
plot(p, [vold v])
plot(1:length(v), [vold v])

# Following 
plot(1:length(v), [vold v], xlims = (462, 470), ticks = 10)
reduce_speed_looking_ahead!(v, p, 467);
reduce_speed_looking_ahead!(v, p, 466);
reduce_speed_looking_ahead!(v, p);


@code_warntype reduce_speed_looking_ahead!(v, p)
# 0.003687 seconds (117.36 k allocations: 2.518 MiB)
@time reduce_speed_looking_ahead!(v, p)



𝐣 = Journey(ea1, no1, ea2, no2);
sol = solve_journey(𝐣);

pl = plot_journey(sol; tit, xtime = false)
pl = plot_journey(sol; tit, xtime = true)
plot(pl[1], size = (800, 600))