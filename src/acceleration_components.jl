
"""
    gravity_comp_along_surface(slope::Float64)
    --> Quantity{Float64, 𝐋 𝐓⁻²}
"""
gravity_comp_along_surface(slope::Float64) = -g * sin(atan(slope))


"""
    @kwdef struct AirAcceleration
        m::Tmass = VEHICLE_DEFAULTS.mass
        ρ::typeof(1.0u"kg/m^3") = ρ_humid(ENVIRONMENT_DEFAULTS.T)
        Cs::Float64 = VEHICLE_DEFAULTS.shapecoeff
        A::Tarea = VEHICLE_DEFAULTS.frontarea
        f::typeof(wind_force) = wind_force
    end

Refer to `wind_force`.
"""
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

A rotary inertia equivalent should be included in the first argument.
It's a ~2% contribution.
"""
function motor_acceleration_limit(m, P, F, rmp, v)
    forcelim_power = P / v
    if v <= zero(v)
        rampfrac = 0.5
    elseif v < rmp
        rampfrac = 0.5 + 0.5 * v / rmp
    else
        rampfrac = 1.0
    end
    forcelim_ramp = rampfrac * F
    if v <= zero(v)
        forcelim_power = forcelim_ramp
    else
        forcelim_power = P / v
    end
    max(zero(Tacc), Tacc(min(forcelim_power, forcelim_ramp) / m))
end


"""
    @kwdef struct MotorlimAcceleration
        m::Tmass = VEHICLE_DEFAULTS.mass
        mr::Tmass = VEHICLE_DEFAULTS.massrot
        P::Tpower = VEHICLE_DEFAULTS.power
        F::Tforce = VEHICLE_DEFAULTS.motorlim
        rmp::Tvel = VEHICLE_DEFAULTS.rampvel
        f::typeof(motor_acceleration_limit) = motor_acceleration_limit
    end

Refer to `motor_acceleration_limit`.
"""
@kwdef struct MotorlimAcceleration
    m::Tmass = VEHICLE_DEFAULTS.mass
    mr::Tmass = VEHICLE_DEFAULTS.massrot
    P::Tpower = VEHICLE_DEFAULTS.power
    η::Float64 = VEHICLE_DEFAULTS.η
    F::Tforce = VEHICLE_DEFAULTS.motorlim
    rmp::Tvel = VEHICLE_DEFAULTS.rampvel
    f::typeof(motor_acceleration_limit) = motor_acceleration_limit
end
# Callable with velocity
(a::MotorlimAcceleration)(v) = a.f(a.m + a.mr, a.P * a.η, a.F, a.rmp, v)

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
"""
    @kwdef struct RollRAcceleration
        Cr::Float64 = VEHICLE_DEFAULTS.rollcoeff
        f::typeof(roll_acceleration) = roll_acceleration
    end

Refer to `roll_acceleration`.
"""
@kwdef struct RollRAcceleration
    Cr::Float64 = VEHICLE_DEFAULTS.rollcoeff
    f::typeof(roll_acceleration) = roll_acceleration
end
# Callable, no argument
(a::RollRAcceleration)() =  -roll_acceleration(a.Cr)
