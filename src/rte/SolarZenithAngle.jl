module SolarZenithAngle

export calc_point_insolation,
       calc_day_lat_insolation,
       linspace
"""
    linspace(start,stop; num = 100)

Wrapper for Matlab's linspace function.
"""
linspace(start,stop; num = 100) = collect(range(start,stop=stop,length=num))

"""
    calc_point_insolation(t, ϕ, γ::FT, π_, e) where FT

"""
function calc_point_insolation(t, ϕ, γ::FT, π_, e) where FT
  # convert inputs from degrees to radians
  ϕ  = ϕ  * 2*π/FT(360)
  γ  = γ  * 2*π/FT(360)
  π_ = π_ * 2*π/FT(360)

  # constants
  Ya = 365.26 # days
  t_VE = 76.0 # days since Jan 1
  S_0 = 1362.0 # W/m^2

  # step 1, calculate the mean anomaly at vernal equinox
  β = sqrt(1-e^2)
  M_VE = -π_ + (e + e^3/4)*(1 + β)*sin(π_)

  # step 2, calculate the mean anomaly
  M = (2*π*(t-t_VE))/(Ya) + M_VE

  # step 3, calculate the true anomaly
  A = M + (2*e - e^3/4)*sin(M)

  # step 4, calculate the distance to the sun
  d = (1-e^2)/(1+e*cos(A))

  # step 5, calculate the solar longitude
  L_s = A + π_

  # step 6, calculate the declination angle
  δ = asin(sin(γ) * sin(L_s))

  # step 7, calculate the sunrise/sunset angle
  T = tan(ϕ) * tan(δ)
  if T >= 1
    η_d = π
  elseif T <= -1
    η_d = FT(0)
  else
    η_d = acos(-1 * T)
  end

  # step 8, calculate the daily averaged cos(zenith angle)
  c1 = η_d*sin(ϕ)*sin(δ)
  c2 = cos(ϕ)*cos(δ)*sin(η_d)
  cosbar = (1/π)*(c1+c2)

  # step 9, calculate the flux
  F = S_0 * (1/d)^2 * cosbar
  return F
end

"""
    calc_day_lat_insolation(n_days::I,
                            n_lats::I,
                            γ::FT,
                            π_,
                            e) where {FT<:AbstractFloat,I<:Int}


"""
function calc_day_lat_insolation(n_days::I,
                                 n_lats::I,
                                 γ::FT,
                                 π_,
                                 e) where {FT<:AbstractFloat,I<:Int}
  d_arr = Array{I}(round.(linspace(0,365,num = n_days)))
  l_arr = Array{I}(round.(linspace(-90,90,num = n_lats)))
  F_arr = zeros(FT, n_days,n_lats)
  # loop over days
  for (i,day) in enumerate(d_arr)
    for (j,lat) in enumerate(l_arr)
      F_arr[i,j] = calc_point_insolation(day,lat,γ,π_,e)
    end
  end
  return F_arr
end

end
