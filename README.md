# DrivingTime.jl

Estimate driving time from start and end coordinates.

The route, road geometry and speed limits are supplied by Norsk Vegdatabase. The unofficial vehicle data is calibrated for an electric bus, Yutong IC12E 2025. Speed in curves is calibrated against data from an installed acceleration monitoring and warning system - the driver aims for driving as fast as possible without triggering warnings.

Given start and end coordinates, 
- retrieves route geometry 
- builds an ordinary differential equation
- solves the equation

### Examples

The simple invocation:

```
julia> using DrivingTime

julia> drivetime(24062, 6939037, 28592, 6939504)
8 minutes
```

For verifying a route, make the route data explicit:

```
julia> using DrivingTime

julia> d = continuous_route_data([(24062, 6939037), (28592, 6939504)])
Dict{Symbol, Any} with 9 entries:
  :radius_of_curvature         => [NaN, NaN, -554.085, 416.561, -17327.8, NaN, NaN, N…
  :multi_linestring            => [[(24062.0, 6.93904e6, 19.668), (24072.0, 6.93904e6…
  :fartsgrense_tuples          => [(1.0, 50, 50), (1.0, 50, 50), (1.0, 50, 50), (1.0,…
  :prefixed_vegsystemreferanse => ["1515 KV3220 S1D1 m92-188", "1515 KV3220 S2D1 m366…
  :key                         => "(24062 6939037)-(28592 6939504)"
  :progression                 => [0.0, 9.96505, 26.1712, 42.0826, 58.6966, 73.9979, …
  :speed_limitation            => [50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50…
  :slope                       => [0.00810815, 0.00900124, 0.0142226, 0.0143325, 0.00…
  :progression_at_ends         => [0.0, 95.619, 380.259, 431.847, 745.964, 845.411, 8…

julia> j = Journey(d);

julia> drivetime(j)
7 minutes, 48 seconds
```


This package comes with custom, unitful plots focused on verifying the simulation:

```
julia> plot_elevation_and_slope_vs_progression(d, "from", "to")

julia> sol = solve_journey(j);

julia> plot_journey(sol; tit = "From -> To")
```

For birds-eye verification of the route travelled, take data from the dictionary out of 'continuous_route_data'. You can unpack `d[:multi_linestring]` by calling  `DrivingTime.RouteSlopeDistance.unique_unnested_coordinates_of_multiline_string`.

