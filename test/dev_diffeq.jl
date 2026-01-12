using Test
using DrivingTime

na1 = "Eiksund"
ea1 = 27978
no1 = 6935574
na2 = "Sørestranda"
no2 = 6928773
ea2 = 33196
tit = na1 * " -> " * na2
𝐣 = Journey(ea1, no1, ea2, no2);
@test sizeof(𝐣) == 128

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
# @time solve_journey(𝐣);
sol = solve_journey(𝐣);
sol.t[end]
plot_journey(sol; tit)

DrivingTime.slope_angle(𝐣, 0.0u"km")
plot_journey(sol[1:11])
