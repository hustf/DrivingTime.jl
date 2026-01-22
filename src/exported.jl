"""
    drivetime(ea1, no1, ea2, no2)
    --> Dates.Minute

Starting point at Utm33 coordinates ("easting, northing") `ea1`, `no1`. Starting velocity: 0 km/hr.

End point at Utm33 coordinates ("easting, northing") `ea1`, `no1`. Starting velocity: 1 km/hr.

The route, road geometry and speed limits are taken from Norsk Vegdatabase's API if not already cached locally. 
Driving time is found by integration of a simple model involving:

 - distance (of course)
 - speed limits (slow down in advance of reduced speed limit)
 - slope (affects acceleration and sometimes max speed)
 - curvature (limits speed through a centripetal acceleration limit)
 - speed bumps (limits velocity to 15 km/hr below the local speed limit)
 - unofficial vehicle data for an electric intercity bus, Yutong IC12E
 - temperature, air pressure and resistance
 - roll resistance

Prior to starting the integration, the model parameters are collected in a `Journey` object. For detail control
or for output resolution, see:

  - `Journey`
  - `solve_journey`
  - `plot_journey`
  - `slope`
  - `@u_str`

# Example
```
julia> using DrivingTime

julia> drivetime(24062, 6939037, 28592, 6939504)
8 minutes
```

"""
function drivetime(ea1, no1, ea2, no2)
    sol = solve_journey(Journey(ea1, no1, ea2, no2))
    t = sol.t[end]
    # Unitful.Time to Dates.Period
    Minute(round(round(t / u"minute")))
end