# Relevant to air resistance

"""
    psat(T)

Saturation vapor pressure, pressure when relative humidity is 100%.
[Ref](https://en.wikipedia.org/wiki/Density_of_air#cite_note-wahiduddin_02-18)

# Example

```
julia> DrivingTime.psat(100.0u"°C")
102.21236954338666 kPa

julia> DrivingTime.psat(00.0u"°C")
0.61078 kPa
```
"""
function psat(T)
    0.61078u"kPa" * exp(17.27 * (T - 273.15u"K")
                         / (T- 35.85u"K"))
end


"""
    ρ_humid(T, pₐᵢᵣ, pᵥ, Rₐᵢᵣ, Rᵥ)
    ρ_humid(T, pᵥ)
    ρ_humid(T)
    --> Quantity{𝐌 𝐋⁻³}

# Example
```
julia> DrivingTime.ρ_humid(0.0u"°C")
1.289303331662808 kg m⁻³

julia> DrivingTime.ρ_humid(18.0u"°C")
1.2030217649479937 kg m⁻³

julia> DrivingTime.ρ_humid(100.0u"°C")
0.5852599086392304 kg m⁻³
```
"""
function ρ_humid(T, pₐᵢᵣ, pᵥ, Rₐᵢᵣ, Rᵥ)
    ρ = pₐᵢᵣ / (Rₐᵢᵣ * T) + pᵥ / (Rᵥ * T)
    Unitful.upreferred(ρ)
end
function ρ_humid(T, pᵥ)
    pₐₜₘ = ENVIRONMENT_DEFAULTS.pₐₜₘ 
    Rₐᵢᵣ = ENVIRONMENT_DEFAULTS.Rₐᵢᵣ
    Rᵥ   = ENVIRONMENT_DEFAULTS.Rᵥ
    pₐᵢᵣ   = pₐₜₘ - pᵥ  # Dalton
    ρ_humid(T, pₐᵢᵣ, pᵥ, Rₐᵢᵣ, Rᵥ)
end
function ρ_humid(T)
    ρ_humid(u"K"(T), psat(T))
end

"""
wind_force(Cs, A, ρ, v)
--> Quantity{ 𝐋 𝐌 𝐓⁻² }
"""
wind_force(Cs, A, ρ, v) = Unitful.upreferred(Cs * A * 0.5 * ρ * v^2)
