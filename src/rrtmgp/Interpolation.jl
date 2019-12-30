#####
##### Interpolation
#####

"""
    extrap_lower(p_lay::Array{FT},
                 p_lev::FT,
                 ϕ_lay::Array{FT})

Extrapolate lower boundary (weighted by pressure)
"""
function extrap_lower(p_lay::Array{FT},
                      p_lev::FT,
                      ϕ_lay::Array{FT}) where {FT<:AbstractFloat}
  return ϕ_lay[1] + (p_lev-p_lay[1])* (ϕ_lay[2]-ϕ_lay[1])/ (p_lay[2]-p_lay[1])
end

"""
    extrap_upper(p_lay::Array{FT},
                 p_lev::FT,
                 ϕ_lay::Array{FT}) where {FT<:AbstractFloat,I<:Int}

Extrapolate upper boundary (weighted by pressure)
"""
function extrap_upper(p_lay::Array{FT},
                      p_lev::FT,
                      ϕ_lay::Array{FT}) where {FT<:AbstractFloat,I<:Int}
  return ϕ_lay[2] + (p_lev-p_lay[2])* (ϕ_lay[2]-ϕ_lay[1])/ (p_lay[2]-p_lay[1])
end

"""
    interpolate_var(p_lay::Array{FT},
                    p_lev::FT,
                    ϕ_lay::Array{FT}) where {FT<:AbstractFloat,I<:Int}

Interpolate variable ϕ (weighted by pressure)
"""
function interpolate_var(p_lay::Array{FT},
                         p_lev::FT,
                         ϕ_lay::Array{FT}) where {FT<:AbstractFloat,I<:Int}
   return (p_lay[1]*ϕ_lay[1]* (p_lev-p_lay[2]) +
           p_lay[2]*ϕ_lay[2]* (p_lay[1]-p_lev))/
          (p_lev*(p_lay[1] - p_lay[2]))
end
function interpolate_var(p_lay::Array{FT},
                         p_lev::Array{FT},
                         ϕ_lay::Array{FT},
                         ncol::I,
                         nlay::I) where {FT<:AbstractFloat,I<:Int}
  ϕ_lev = zeros(FT, ncol,nlay+1)
  # Interpolate temperature to levels if not provided
  #   Interpolation and extrapolation at boundaries is weighted by pressure
  for icol = 1:ncol
     ϕ_lev[icol,1] = extrap_lower(p_lay[icol,1:2],p_lev[icol,1],ϕ_lay[icol,1:2])
  end
  for ilay in 2:nlay
    for icol in 1:ncol
      ϕ_lev[icol,ilay] = interpolate_var(p_lay[icol,ilay-1:ilay],p_lev[icol,ilay],ϕ_lay[icol,ilay-1:ilay])
    end
  end
  for icol = 1:ncol
    ϕ_lev[icol,nlay+1] = extrap_upper(p_lay[icol,nlay-1:nlay],p_lev[icol,nlay+1],ϕ_lay[icol,nlay-1:nlay])
  end
  return ϕ_lev
end
