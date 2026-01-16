using Test
using DrivingTime
import RouteSlopeDistance
using RouteSlopeDistance: route_leg_data, unique_unnested_coordinates_of_multiline_string
using Plots
import Unitful: upreferred
import Statistics
# Based on hand notes from Ecosafe
# The first step is to identify locations and corresponding radius of curvature
function plot_xy_and_r(d::Dict, tit)
    pl = plot(; layout = (2, 1), size = (1600, 800))
    vea, vno, _ = unique_unnested_coordinates_of_multiline_string(d[:multi_linestring])
    plxy = pl[1]
    plot!(plxy, vea, vno, label = false, markersize = 3, marker = :xcross)
    # Scale both axes equal for the y-x plot
    Δ =  max(abs(-(ylims(plxy)...)), abs(-(xlims(plxy)...))) * 1.25 
    xcen = +(xlims(plxy)...) / 2
    ycen = +(ylims(plxy)...) / 2
    xlims!(plxy, (xcen - Δ / 2, xcen + Δ / 2))
    ylims!(plxy, (ycen - Δ / 2, ycen + Δ / 2))
    title!(plxy, tit)
    # Mark start, end in top plot
    t = text("Start", 8, :left, :bottom, :green, rotation = -90)
    y = vno[1]
    x = vea[1]
    annotate!(plxy, [(x, y, t)])
    t = text("End", 8, :centre, :bottom, :green, rotation = -90)
    y = vno[end]
    x = vea[end]
    annotate!(plxy, [(x, y, t)])
    # Mark interesting radii
    r = d[:radius_of_curvature] * u"m"
    ra = abs.(replace(r, NaN * u"m" => 1000.0u"m"))
    rmi, imi = findmin(ra)
    t = text("$(round(u"m", rmi))", 12, :centre, :bottom, :black, rotation = 0)
    y = vno[imi]
    x = vea[imi]
    annotate!(plxy, [(x, y, t)])
    #
    plr = pl[2]
    p = d[:progression] * u"m"
    plot!(plr, p, r,  label = "Signed radius(progression)")
    y = r[imi]
    x = p[imi]
    annotate!(plr, [(x, y, t)])
    ylims!(plr, (-200.0, 200.0))
    hline!(plr, [0], label = "", linestyle = :dash, color = :red)
    display(plot(plxy, size = (800, 800)))
    pl
end



function endpoints(i) 
    @assert i <= 15
    kalibrering = ["https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=23593.713839066448,6942485.5900078565&slutt=23771.052968726202,6942714.9388697725&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/lsvegnett/rute?start=10929.721370896965,6932488.729146618&slutt=10906.172339663783,6932221.839157347&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=38744.875253753446,6946346.347827989&slutt=38856.09116223373,6946141.295902717&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=17233.49,6933166.44&slutt=17411.16742345906,6933361.249064187&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=19888.753066316945,6944574.509461373&slutt=20264.900856445194,6944368.952139658&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=26521.552530414192,6940216.455383167&slutt=26457.440704222478,6940123.788953421&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=26343.660416060884,6950023.898755312&slutt=26736.355686723255,6950277.5537265055&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=28758.084434766264,6945107.442628212&slutt=28872.074380539998,6945084.460240152&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=20957.20257390075,6939620.818590216&slutt=20637.109098991437,6939678.886214065&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=12406.263465148688,6933858.521046377&slutt=12658.589337749,6933842.480074977&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=23131.961055209278,6936838.823032396&slutt=23178.72879053885,6937309.664867436&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=36702.36624956445,6950291.6474007955&slutt=36778.50157994096,6950016.548198562&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=20853.204308743123,6939592.358717515&slutt=20638.11771304172,6939679.444927726&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=16112.167162967904,6946267.653663013&slutt=16176.453134912474,6946340.517977856&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true",
                    "https://nvdbapiles-v3.atlas.vegvesen.no/beta/vegnett/rute?start=38792.42584442743,6946386.567167263&slutt=38946.91037624149,6946291.330794491&maks_avstand=10&omkrets=100&konnekteringslenker=true&detaljerte_lenker=false&behold_trafikantgruppe=false&trafikantgruppe=K&pretty=true&kortform=true"]
    nametags =      [
        "Storesandvika Dimna r68m r120m", 
        "Larsnes Halsen r66m ", 
        "Hjørungavåg -> Hareid r94", 
        "Vågen-> Leikong ved Skoge r140m", 
        "Nord for Herøybrua mot sør r150m", 
        "Garneskrysset inn, høgre r8m ", 
        "Gåsnesbukta mot Flø, r114m", 
        "Inn på Smårisevadet r?m", 
        "Tjørvåg Tuftevatnet -> Nybøvegen r149m", 
        "Drageskaret Larsnes -> Vågen r103m", 
        "Leikong Budanesvågen mot nord r108m", 
        "Brandalsvegen mot sør, første nes r116m", 
        "Tjørvåg Tuftevatnet -> Nybøvegen r145m", 
        "Herøy vgs ned mot symjehall, venstre 20m høgre 55m", 
        "Nyevegen ved Ovra, mot Hjørungavåg r210m"]
    @assert length(kalibrering) == length(nametags) == 15
    s = kalibrering[i]
    args = split(split(s, '?')[2], '&')
    start = split(args[1], '=')[2]
    slutt = split(args[2], '=')[2]
    stea, stno = split(start, ',')
    slea, slno = split(slutt, ',')
    ea1 = Int(round(tryparse(Float64, stea)))
    no1 = Int(round(tryparse(Float64, stno)))
    ea2 = Int(round(tryparse(Float64, slea)))
    no2 = Int(round(tryparse(Float64, slno)))
    ea1, no1, ea2, no2, nametags[i]
end





ea1, no1, ea2, no2, nametag = endpoints(12)
d = route_leg_data(ea1, no1, ea2, no2);
"$(Int(round(ea1))) , $(Int(round(no1)))" |> clipboard
"$(Int(round(ea2))) , $(Int(round(no2)))" |> clipboard


pl = plot_xy_and_r(d, nametag)
# 
plot(pl[1], size = (800, 800))
yaxis!(pl[2], tick = 20)
plot(pl[2], size = (800, 800), ylims = (-300, 300))


# Observations, summary from notes. Dropping details like bus, right or left, units
lv1v    51  76   65  28  30 29 57  30 29  29   57 80   63  64  82   63   61   69   59  
lv1r    66  140  150  8   8  8 114  8 8    8  114 149 103 108  149  108  116  114  116 
lv2v    72   57  79   62 62  
lv2r   114   116 210  68 68   
lv3v   60  53  59  58  62    
lv3r   67  94  68  68  120

# Make every entry unique by sorting level 1, then adding 0.2 for duplicates.

lv1v = [28, 29, 29.2, 29.4,  30, 30.2, 51, 57, 57.2,  59,  61,  63,  63,  64,  65,  69,  76,  80, 82] * u"km/hr"
lv1r = [8 ,  8,    8,  8.1,   8,    8, 66, 114, 114, 116, 116, 103, 108, 108, 150, 114, 140, 149, 149] * u"m"  
lv2v = [  72 , 57 ,  79 , 62 , 62.1] * u"km/hr"
lv2r = [ 114 ,116 , 210 , 68 , 68] * u"m"
lv3v = [ 60  , 53 , 59  , 58 , 62] * u"km/hr"
lv3r = [ 67  , 94 , 68  , 68 , 120] * u"m"

pl = plot(lv1v, lv1r; seriestype = :scatter, label = "lv1", xlim = (0, 85))
plot!(pl, lv2v, lv2r; seriestype = :scatter, label = "lv2")
plot!(pl, lv3v, lv3r; seriestype = :scatter, label = "lv3")
# We would expect to see a parabolic trendline,
# where v^2/r. It's not that evident here.
# Still, we choose this
lv1_acc = map(eachindex(lv1v)) do i
   upreferred(lv1v[i]^2 / lv1r[i])
end

lv2_acc = map(eachindex(lv2v)) do i
   upreferred(lv2v[i]^2 / lv2r[i])
end

lv3_acc = map(eachindex(lv3v)) do i
   upreferred(lv3v[i]^2 / lv3r[i])
end
# To connect points visibly...
pl = plot(layout = (2,2), seriestype = :scatter)
plot!(pl[2], lv1v, lv1r; seriestype = :scatter, label = "lv1", xlim = (0, 85))
plot!(pl[1], lv1_acc, lv1r; seriestype = :scatter, label = "lv1")
plot!(pl[4], lv1v, lv1_acc; seriestype = :scatter, label = "lv1", xlim = (0, 85))

# Its' pretty clear that we can disregard the < 40 km/ hr data points entirely.
lv1v = lv1v[7:end]
lv1r = lv1r[7:end]
lv1_acc = lv1_acc[7:end]
# Now plot again
pl = plot(layout = (2,2), seriestype = :scatter)
plot!(pl[2], lv1v, lv1r; seriestype = :scatter, label = "lv1", xlim = (0, 85))
plot!(pl[1], lv1_acc, lv1r; seriestype = :scatter, label = "lv1")
plot!(pl[4], lv1v, lv1_acc; seriestype = :scatter, label = "lv1", xlim = (0, 85))
# These all look reasonable. 
lv1_acc_mean = Statistics.mean(lv1_acc)
# 2.7966112269858963 m s⁻²
hline!(pl[4], [lv1_acc_mean], label = "", linestyle = :dash, color = :red)
vline!(pl[1], [lv1_acc_mean], label = "", linestyle = :dash, color = :red)
# For reasonability checks, let's plot the other levels here too. 
plot!(pl[2], lv2v, lv2r; seriestype = :scatter, label = "lv2", xlim = (0, 85))
plot!(pl[1], lv2_acc, lv2r; seriestype = :scatter, label = "lv2")
plot!(pl[4], lv2v, lv2_acc; seriestype = :scatter, label = "lv2", xlim = (0, 85))
lv2_acc_mean = Statistics.mean(lv2_acc)
# 3.3401647538370502 m s⁻²
hline!(pl[4], [lv2_acc_mean], label = "", linestyle = :dash, color = :green)
vline!(pl[1], [lv2_acc_mean], label = "", linestyle = :dash, color = :green)
#
plot!(pl[2], lv3v, lv3r; seriestype = :scatter, label = "lv3", xlim = (0, 85)color = :black)
plot!(pl[1], lv3_acc, lv3r; seriestype = :scatter, label = "lv3", color = :black)
plot!(pl[4], lv3v, lv3_acc; seriestype = :scatter, label = "lv3", color = :black, xlim = (0, 85))
lv3_acc_mean = Statistics.mean(lv3_acc)
# 3.3381083340463276 m s⁻², lower than lv2 but we have few sample
hline!(pl[4], [lv3_acc_mean], label = "", linestyle = :dash, color = :black)
vline!(pl[1], [lv3_acc_mean], label = "", linestyle = :dash, color = :black)

# => A reasonable acc level for passenger comfort driving is 
round(u"m/s^2", minimum(lv1_acc), digits = 1)
# 2.2 m/s^2 = v² / r
# => A reasonable speed is  
#   v < √(r * 2.2 m/s^2) 