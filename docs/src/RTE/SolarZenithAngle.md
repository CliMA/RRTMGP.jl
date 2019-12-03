# Solar Zenith Angle

```@meta
CurrentModule = RRTMGP.SolarZenithAngle
```

# Examples

```@example ZenithAngle
include("plot_zenith_angle.jl")

#####
##### Example 2
#####

# 2018 constants
γ_0 = 23.44
π_0 = 282.95
e_0 = 0.017

# ndays, nlats
ndays = 365
nlats = 180

F0 = calc_day_lat_insolation(ndays, nlats, γ_0, π_0, e_0)
plot_day_lat_insolation(ndays, nlats, F0, "YlOrRd", "gamma_0 = $(γ_0), pi_0 = $(π_0), e = $(e_0)", "example2.png")

#####
##### Example 3
#####

# turn longitude of perihelion by 180°
γ_0 = 23.44
π_1 = 282.95 + 180.0
e_0 = 0.017

ndays = 365
nlats = 180

F1 = calc_day_lat_insolation(ndays,nlats,γ_0,π_1,e_0)
plot_day_lat_insolation(ndays,nlats,F1,"YlOrRd",  "γ_0 = $(γ_0)^∘, π_0 = $(π_0)^∘, e = $(e_0)", "example3a.png")
plot_day_lat_insolation(ndays,nlats,F1-F0,"PRGn", "ToA Insolation Difference: π_0_1 = 102.95^∘ - π_0_0 = 282.95^∘", "example3b.png")

#####
##### Example 4
#####

# perihelion back to normal. decrease γ to 22.0°
γ_1 = 22.0
π_0 = 282.95
e_0 = 0.017

# ndays, nlats
ndays = 365
nlats = 180

F2 = calc_day_lat_insolation(ndays,nlats,γ_1,π_0,e_0)
plot_day_lat_insolation(ndays,nlats,F2,"YlOrRd",  "γ_0 = $(γ_0)^∘, π_0 = $(π_0)^∘, e = $(e_0)", "example4a.png")
plot_day_lat_insolation(ndays,nlats,F2-F0,"PRGn", "ToA Insolation Difference: γ_1 = $(γ_1)^∘ - γ_0 = $(γ_0)^∘", "example4b.png")

# decrease γ to 18.0°
γ_2 = 18.0
π_0 = 282.95
e_0 = 0.017

# ndays, nlats
ndays = 365
nlats = 180

F3 = calc_day_lat_insolation(ndays,nlats,γ_2,π_0,e_0)
plot_day_lat_insolation(ndays,nlats,F3,"YlOrRd", "γ_0 = $(γ_0)^∘, π_0 = $(π_0)^∘, e = $(e_0)", "example4c.png")
plot_day_lat_insolation(ndays,nlats,F3-F0,"PRGn", "ToA Insolation Difference: γ_1 = $(γ_1)^∘ - γ_0 = $(γ_0)^∘", "example4d.png")

#####
##### Example 6
#####

# now change obliquity to 60.0°
γ_3 = 60.0
π_0 = 282.95
e_0 = 0.017

# ndays, nlats
ndays = 365
nlats = 180

F4 = calc_day_lat_insolation(ndays,nlats,γ_3,π_0,e_0)
plot_day_lat_insolation(ndays,nlats,F4,"YlOrRd", "γ_0 = $(γ_0)^∘, π_0 = $(π_0)^∘, e = $(e_0)", "example6a.png")
plot_day_lat_insolation(ndays,nlats,F4-F0,"PRGn", "ToA Insolation Difference: γ_1 = 60.0^∘ - γ_0 = 23.44^∘", "example6b.png")

# now change obliquity to 97.86°
γ = 97.86
π_0 = 282.95
e_0 = 0.017

# ndays, nlats
ndays = 365
nlats = 180

F5 = calc_day_lat_insolation(ndays,nlats,γ,π_0,e_0)
plot_day_lat_insolation(ndays,nlats,F5,"YlOrRd", "γ_0 = $(γ_0)^∘, π_0 = $(π_0)^∘, e = $(e_0)", "example6c.png")
plot_day_lat_insolation(ndays,nlats,F5-F0,"PRGn", "ToA Insolation Difference: γ_1 = 97.86^∘ - γ_0 = 23.44^∘", "example6d.png")

```
## Example 2
![](example2.png)
## Example 3
![](example3a.png)
![](example3b.png)
## Example 4
![](example4a.png)
![](example4b.png)
![](example4c.png)
![](example4d.png)
## Example 6
![](example6a.png)
![](example6b.png)
![](example6c.png)
![](example6d.png)

