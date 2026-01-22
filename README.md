# DrivingTime.jl

Estimate driving time from start and end coordinates.

The route, road geometry and speed limits are supplied by Norsk Vegdatabase. The unofficial vehicle data is calibrated for an electric bus, Yutong IC12E 2025.

Given start and end coordinates, 
- retrieves route geometry 
- builds an ordinary differential equation
- solves the equation

The simple invocation:

```
julia> using DrivingTime

julia> drivetime(24062, 6939037, 28592, 6939504)
8 minutes
```

This comes with custom, unitful plots. See inline docs.