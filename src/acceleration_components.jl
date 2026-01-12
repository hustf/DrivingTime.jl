
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
