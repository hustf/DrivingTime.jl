module DrivingTime
import Interpolations
using Interpolations: Extrapolation, ScaledInterpolation, BSplineInterpolation
using Interpolations: BSpline, Cubic, Line, OnGrid
using Interpolations: extrapolate, interpolate, Gridded, Linear, Throw, scale
import RouteSlopeDistance
using RouteSlopeDistance: route_leg_data
import Unitful
using Unitful: Length, Velocity, Acceleration, Time
using Unitful: @u_str, Quantity, NoUnits
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

export Journey, solve_journey, plot_journey, slope, @u_str

const g = 9.81u"m/s^2"
const VEHICLE_DEFAULTS = (;mass = 15300.0u"kg", # Egenvekt med fører 
                           power = 350.0u"kW",  # (merk: "280 kW per 30 min")
                           motorlim = 21.25u"kN", # This is unknown, but we can't have wheels spinning at low speed.
                           frontarea = 2.55u"m" * 3.5u"m",
                           shapecoeff = 1.1
                           )
const ENVIRONMENT_DEFAULTS = ( 
                              T = 278.15u"K", # 5°C
                              p    = 101.325u"kPa", # Sea level
                              Rₐᵢᵣ =  287.058u"J * kg^-1 * K^-1", # Specific gas constant, dry air 
                              Rᵥ   = 461.495u"J * kg^-1 * K^-1", # Specific gas constant for water vapor
                              pₐₜₘ = 101.325u"kPa" # Sea level standard pressure
                           )
include("thermodynamics.jl")
include("acceleration_components.jl")


struct Journey{S <:Extrapolation, F1, F2, F3}
   fslope::S
   fslopeacc::F1
   fairacc::F2
   fmotoracclim::F3
end

include("show.jl")
include("define_journey.jl")
include("diffeq.jl")
include("plot_journey.jl")

end