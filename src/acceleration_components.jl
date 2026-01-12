
"""
    gravity_comp_along_surface(slope::Float64)
    --> Quantity{Float64, 𝐋 𝐓⁻²}
"""
gravity_comp_along_surface(slope::Float64) = -g * sin(atan(slope))

@kwdef struct AirAcceleration
    m::Tmass = VEHICLE_DEFAULTS.mass
    ρ::typeof(1.0u"kg/m^3") = ρ_humid(ENVIRONMENT_DEFAULTS.T)
    Cs::Float64 = VEHICLE_DEFAULTS.shapecoeff
    A::Tarea = VEHICLE_DEFAULTS.frontarea
    f::typeof(wind_force) = wind_force
end
# Callable with velocity
(a::AirAcceleration)(v) = -a.f(a.Cs, a.A, a.ρ, v) / a.m

"""
    motor_acceleration_limit(m, P, F, v)
    ---> Tacc

Acceleration from motor at full "throttle", limited
by vehicle defaults.

Rotary inertia is neglected (no gear shift for electric vehicle).
"""
function motor_acceleration_limit(m, P, F, v)
    max(zero(Tacc), Tacc(min(P / v, F) / m))
end
@kwdef struct MotorlimAcceleration
    m::Tmass = VEHICLE_DEFAULTS.mass
    P::Tpower = VEHICLE_DEFAULTS.power
    F::Tforce = VEHICLE_DEFAULTS.motorlim
    f::typeof(motor_acceleration_limit) = motor_acceleration_limit
end
# Callable with velocity
(a::MotorlimAcceleration)(v) = a.f(a.m, a.P, a.F, v)

"""
    roll_acceleration(Cr)

`Cr` is a "lumped road-load coefficient", higher than a rolling coefficient.

Output is intended mostly for estimating energy consumption in an electric heavy intercity bus in
a hilly region. 

Vehicle mass disappears from the equation since output is acceleration, not resistance.

The acceleration always works backwards, since we're not modelling reverse speed.

Hence, the Cr value takes into account:

- Tire deformation, energy is lost through imperfect elasticity
- Road deformation inelasticity, spring back is slower than depression
- Internal friction in roller bearings loaded by vehicle weight
- Friction between ground and tire from contact patch slippage (in all 
  horizontal directions) 
- Friction between ground and tire due to longitudinal slip (depends on applied torque)
- Friction between ground and tire due to cornering
- Sound emission in air and ground
- Sometimes wet road
"""
function roll_acceleration(Cr)
    Cr * g
end
# This contribution is a bit over-engineered
# considering the physical simplicity,
# but we try to make a consistent and easily extendable
# model.
@kwdef struct RollAcceleration
    Cr::Float64 = VEHICLE_DEFAULTS.rollcoeff
    f::typeof(roll_acceleration) = roll_acceleration
end
# Callable, no argument
(a::RollAcceleration)() =  roll_acceleration(a.Cr)
