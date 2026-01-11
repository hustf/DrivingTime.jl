using Test
using DrivingTime

na1 = "Eiksund"
ea1 = 27978
no1 = 6935574
na2 = "Sørestranda"
no2 = 6928773
ea2 = 33196
𝐣 = Journey(ea1, no1, ea2, no2);
@test sizeof(𝐣) == 72

Γᵢₙ = ArrayPartition([0.0u"m"], [0.0u"m/s"])
sol = solve_journey(𝐣);
sol.t[end]
plot_journey(sol)
DrivingTime.slope(𝐣, 10.0u"km")