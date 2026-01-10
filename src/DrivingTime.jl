module DrivingTime
import Interpolations
using Interpolations: Extrapolation, ScaledInterpolation, BSplineInterpolation
using Interpolations: BSpline, Cubic, Line, OnGrid
using Interpolations: extrapolate, interpolate, Gridded, Linear, Throw, scale
import RouteSlopeDistance
using RouteSlopeDistance: route_leg_data
import Unitful
using Unitful: Length, Velocity, Acceleration
import StaticArrays
using StaticArrays: MVector
import Base: show, ==
export Journey, State

const g = 9.81u"m/s^2"

pack(x::T) where T = MVector{1,T}([x])
# Type for the local tuple, Γ 
struct State{P,V,A}
   x::MVector{1,P}
   x´::MVector{1,V}
   x´´::MVector{1,A}
   function State(
        x  :: MVector{1,P},
        x´  :: MVector{1,V},
        x´´ :: MVector{1,A}
    ) where {
        P <: Length,
        V <: Velocity,
        A <: Acceleration
    }
        new{P,V,A}(x, x´, x´´)
    end
end
# Convenience constructors
function State(x::P, x´::V, x´´::A) where{P<:Quantity, V<:Quantity, A<:Quantity}
   State(pack(x), pack(x´), pack(x´´))
end
State(; position = 0.0u"m", velocity = 0.0u"km/hr", acceleration = 0.0u"m/s^2") = State(position, velocity, acceleration)
function ==(Γ1::T1, Γ2::T2) where {T1 <:State, T2 <: State}
   Γ1.x == Γ2.x && Γ1.x´ == Γ2.x´ && Γ1.x´´ == Γ2.x´´
end

struct Journey{S}
   fslope::S
end

include("show.jl")
include("define_rhs.jl")
include("diffeq.jl")

end