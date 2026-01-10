using Test
using DrivingTime
using Unitful
using StaticArrays
na1 = "Eiksund"
ea1 = 27978
no1 = 6935574
na2 = "Sørestranda"
no2 = 6928773
ea2 = 33196
rhs = Journey(ea1, no1, ea2, no2);
@test sizeof(rhs) == 72



