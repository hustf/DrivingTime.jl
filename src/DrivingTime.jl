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

pack(x::T) where {T <: Quantity } = MVector{1,T}([x])


struct Journey{S<:Extrapolation}
   fslope::S
end

include("show.jl")
include("define_journey.jl")
include("diffeq.jl")
include("plot_journey.jl")

end