# DrivingTime.jl
Estimate driving time from slope, curvature and vehicle data

Given start and end coordinates, 
- retrieves route geometry 
- builds an ordinary differential equation
- solves the ordinary differential equation

The simple invocation:

`drivetime(ea1, no1, ea2, no2)`
--> Dates.Duration
