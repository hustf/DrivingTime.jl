# Relevant to air resistance

"""
    psat(T)

Saturation vapor pressure, pressure when relative humidity is 100%.
[Ref](https://en.wikipedia.org/wiki/Density_of_air#cite_note-wahiduddin_02-18)

# Example

```
julia> psat(100.0u"°C")
102.21236954338666 kPa

julia> psat(5.0u"°C")
0.8722823974498973 kPa

julia> psat(0.0u"°C")
0.61078 kPa
```
"""
function psat(T)
    0.61078u"kPa" * exp(17.27 * (T - 273.15u"K")
                         / (T- 35.85u"K"))
end


"`ρ(p, R, T)`  # Density"
ρ(p, R, T) = p / (R * T)

"""
    ρ_humid(T, pₐᵢᵣ, pᵥ, Rₐᵢᵣ, Rᵥ)
    ρ_humid(T, pᵥ)
    ρ_humid(T)
    --> Quantity{𝐌 𝐋⁻³}

# Example
```
julia> ρ_humid(0.0u"°C")
1.289303331662808 kg m⁻³

julia> ρ_humid(18.0u"°C")
1.2030217649479937 kg m⁻³

julia> ρ_humid(100.0u"°C")
0.5852599086392304 kg m⁻³
```

Note that for extreme but realistic temperature variation, 
air resistance can vary by 18%. Humidity is neglectable,
but dry conditions are worse.

# Influence of humidity

Density, ρ [kg m⁻³]
```
| Temp.    | Dry air | 100% humid air |  Humid / dry |
|----------|---------|----------------|--------------|
|   0 °C   |  1.292  |  1.289         |   0.997      |
|   5 °C   |  1.269  |  1.264         |   0.996      |
|  18 °C   |  1.212  |  1.203         |   0.992      |
|  40 °C   |  1.127  |  1.096         |   0.972      |
| 100 °C   |  0.945  |  0.585         |   0.619      |
```

 1.2922476553778015 kg m⁻³
 1.2690183248838633 kg m⁻³
 1.2123559919850473 kg m⁻³
 1.1271832893707379 kg m⁻³
 0.9459398286652728 kg m⁻³

 1.289303331662808 kg m⁻³
 1.2648889910131416 kg m⁻³
 1.2030217649479937 kg m⁻³
 1.096171025292019 kg m⁻³
 0.5852599086392304 kg m⁻³
"""
function ρ_humid(T, pₐᵢᵣ, pᵥ, Rₐᵢᵣ, Rᵥ)
    ρh = ρ(pₐᵢᵣ, Rₐᵢᵣ, T) + ρ(pᵥ, Rᵥ, T)
    Unitful.upreferred(ρh)
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
