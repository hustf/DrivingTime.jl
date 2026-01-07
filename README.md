# DrivingTime.jl
Estimate driving time from slope, curvature and vehicle data

Given start and end coordinates, retrieves route geometry and solves the ordinary differential equation.

`drivetime(ea1, no1, ea2, no2)`
--> Dates.Duration
