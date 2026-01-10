using Test
using DrivingTime
using Unitful
using StaticArrays

vx = @MVector([1.0u"m"])
vv = @MVector([1.0u"m/s"])
va = @MVector([1.0u"m/s^2"])
# Inner constructor
@test State(vx, vv, va) isa State
# Convenience constructor
Γ =  State(1.0u"m", 2.0u"m/s", 3.0u"m/s^2")
Γ1 = State(; position = 100.0u"cm", velocity = 2.0u"m/s", acceleration = 3.0u"m/s^2")
@test Γ.x == Γ1.x
@test Γ.x´ == Γ1.x´
@test Γ.x´´ == Γ1.x´´
@test Γ == Γ1
@test repr(Γ) !== repr(Γ1)
# Keyword constructor
@test State(; position = 0.001u"km", velocity = 2.0u"m/s", acceleration = 3.0u"m/s^2") isa State
@test State(; position = 1.0u"km") isa State
@test State() isa State
Γ2 = State()
@test iszero(Γ2.x)
@test iszero(Γ2.x´)
@test iszero(Γ2.x´´)