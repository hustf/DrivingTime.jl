module DrivingTime
import Interpolations
using Interpolations: Extrapolation, ScaledInterpolation, BSplineInterpolation
using Interpolations: BSpline, Cubic, Line, OnGrid
using Interpolations: extrapolate, interpolate, Gridded, Linear, Throw, scale
import RouteSlopeDistance
using RouteSlopeDistance: route_leg_data
import Unitful
using Unitful: Length, Velocity, Acceleration, Time
using Unitful: @u_str, Quantity, NoUnits, dimension
using Unitful: uconvert
using RecursiveArrayTools
import Base: show
import Test
using Test: @inferred
import OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, ContinuousCallback, DiscreteCallback, CallbackSet
using OrdinaryDiffEq: terminate!, Tsit5, solve, get_proposed_dt, init, solve!
import OrdinaryDiffEqCore
using OrdinaryDiffEqCore: ODEIntegrator
import SciMLBase
using SciMLBase: ReturnCode.Success, ReturnCode.Terminated, ReturnCode.DtLessThanMin
using SciMLBase: successful_retcode
using RecursiveArrayTools
using Logging
using Plots
import Dates
using Dates: Minute
export drivetime

const g = 9.81u"m/s^2"

"""
# Calculation of equivalent translatory mass from one wheel rotational inertia:

```
julia> mw = 100kg
100 kg

julia> rw = 0.52m
0.52 m

julia> Iw = 0.8 * mw * rw^2
21.632 kg m²

julia> r_rim = (22.5 * 2.54 / 200)u"m"
0.28575 m

julia> r_outer = r_rim + 0.7 * 0.315u"m"
0.50625 m

julia> θ = 1u"m" / (r_outer) # When the vehicle moves a unit translation, the wheel rotates Θ:
1.9753086419753088

julia> m_eq = Iw / r_outer^2
84.4046944063405 kg

julia> # So for four wheels, the rotational moment of intertia, referred to a translational frame (i.e. the ground moving relative to the vehicle frame)

julia> 4 * m_eq
337.618777625362 kg
```
"""
const VEHICLE_DEFAULTS = (;mass = 15300.0u"kg",    # "Egenvekt med fører" 
                           power = 350.0u"kW",     # (merk: "280 kW per 30 min")
                           η = 0.9,                # Mechanical power efficiency, from Eiksund-Sørheim @73km/hr
                           motorlim = 25.0u"kN",   # On flat ground, power reaches 100% at 50 km/hr
                           rampvel = 0.5u"km/hr",  # Motor limited to 50% at 0, 100% at ramp velocity
                           frontarea = 2.55u"m" * 3.5u"m",
                           shapecoeff = 1.1,
                           rollcoeff = 0.009,    # Lump losses proportional to mass. Winter tires.
                           massrot = 337.6u"kg"  # Four wheels translatory mass equivalent.
                           )
const ENVIRONMENT_DEFAULTS = ( 
                              T = 278.15u"K", # 5°C
                              p    = 101.325u"kPa", # Sea level
                              Rₐᵢᵣ =  287.058u"J * kg^-1 * K^-1", # Specific gas constant, dry air 
                              Rᵥ   = 461.495u"J * kg^-1 * K^-1", # Specific gas constant for water vapor
                              pₐₜₘ = 101.325u"kPa" # Sea level standard pressure
                           )
const Tprog = typeof(1.0u"km")
const Tmass = typeof(VEHICLE_DEFAULTS.mass)
const Tpower = typeof(VEHICLE_DEFAULTS.power)
const Tarea = typeof(VEHICLE_DEFAULTS.frontarea)
const Tforce = typeof(VEHICLE_DEFAULTS.motorlim)
const Tvel = typeof(VEHICLE_DEFAULTS.rampvel)
const Tacc = typeof(g)
include("thermodynamics.jl")
include("acceleration_components.jl")


struct Journey{S <: Extrapolation, T<: Extrapolation}
   pstop::Tprog                       # End progression of the journey
   fslopeacc::S                       # Interpolation(progression)
   fairacc::AirAcceleration           # Callable with velocity
   fmotoracclim::MotorlimAcceleration # Callable with velocity
   frollacc::RollRAcceleration         # Callable with no argument
   itp_v::T                           # Interpolation(progression)
end

include("show.jl")
include("define_journey.jl")
include("diffeq.jl")
include("plot_journey.jl")
include("exported.jl")


end