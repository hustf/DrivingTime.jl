module DrivingTime
import Interpolations
using Interpolations: Extrapolation, ScaledInterpolation, BSplineInterpolation
using Interpolations: BSpline, Cubic, Line, OnGrid
using Interpolations: extrapolate, interpolate, Gridded, Linear, Throw, scale
import RouteSlopeDistance
using RouteSlopeDistance: route_leg_data

struct Journey{S}
   fslope::S
end

export Journey

include("define_rhs.jl")


end