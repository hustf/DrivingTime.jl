using Test
using DrivingTime
using Dates

@testset "I" begin
    na1 = "Møre skule"
    ea1 = 24062
    no1 = 6939037
    na2 = "Ringstaddalen"
    ea2 = 28592
    no2 = 6939504
    tit = na1 * " -> " * na2
    @test drivetime(ea1, no1, ea2, no2) == Minute(7) + Second(49) 
    @test drivetime([(ea2, no2), (ea1, no1)]) == Minute(7) + Second(14) 
end

@testset "II" begin
    # A closer look
    j = Journey([(ea1, no1, ea2, no2)]);
    sol = solve_journey(j);
    pl = plot_journey(sol; tit)
    DrivingTime.plot(pl[1], xlim =(5.5u"km", 6u"km"))
    # Some 'unitless' plots of the input.
    using RouteSlopeDistance
    d = DrivingTime.continuous_route_data([(ea1, no1), (ea2, no2)]; default_fartsgrense = 50u"km/hr")
    plot_elevation_and_slope_vs_progression(d, na1, na2)
    plot_elevation_slope_speed_vs_progression(d, na1, na2)
    # Multiple waypoints (joined route data)
    @test drivetime([(ea1, no1), 
                    (27406, 6938058),
                    (ea2, no2)]) == Minute(9) + Second(12) 
    @test drivetime([(ea1, no1), 
                    (24740, 6939113),
                    (27406, 6938058),
                    (ea2, no2)]) == Minute(9) + Second(10)
end

@testset "III" begin
    # Some serious downhill - motor braking required.
    na1 = "Eiksundbrua"
    ea1 = 27978
    no1 = 6935574
    na2 = "Sørheim"
    ea2 = 32428
    no2 = 6930729
    tit = na1 * " -> " * na2
    j = Journey([(ea1, no1, ea2, no2)]);
    sol = solve_journey(j);
    size = (900, 600)
    pl = plot_journey(sol; tit, size)
    DrivingTime.plot(pl[1]; size)
    DrivingTime.plot(pl[2]; size)
    @test drivetime(j) == Minute(8) + Second(23)
    j = Journey([(ea2, no2, ea1, no1)]);
    @test drivetime(j) == Minute(8) + Second(48)
end
