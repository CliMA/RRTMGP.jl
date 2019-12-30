#####
##### Numeric calculations for gas optics. Absorption and Rayleigh optical depths, source functions.
#####

PaTohPa(::Type{FT}) where FT = FT(0.01)

"""
    interpolation!(ics::InterpolationCoefficients{FT,I},
                   go::AbstractGasOptics{FT,I},
                   as::AtmosphericState{FT,I}) where {I<:Int,B<:Bool,FT<:AbstractFloat}

Compute interpolation coefficients

 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)

for calculations of major optical depths, minor optical depths, Rayleigh,
and Planck fractions given

 - `go` gas optics, see [`AbstractGasOptics`](@ref)
 - `as` atmospheric state, see [`AtmosphericState`](@ref)

# Local variables
 - `ftemp`, `fpress` interpolation fraction for temperature, pressure
 - `locpress` needed to find location in pressure grid
 - `ratio_η_half` ratio of vmrs of major species that defines η=0.5
                  for given flavor and reference temperature level
 - `η`, `f_η` binary_species_parameter, interpolation variable for η
 - `loc_η` needed to find location in η grid

!!! Reduce allocations
"""
function interpolation!(ics::InterpolationCoefficients{FT,I},
                        go::AbstractGasOptics{FT,I},
                        as::AtmosphericState{FT,I}) where {I<:Int,B<:Bool,FT<:AbstractFloat}

  nflav = get_nflav(go)
  neta  = get_neta(go)
  npres = get_npres(go)
  ntemp = get_ntemp(go)

  @unpack_fields as col_gas t_lay p_lay tropo ncol nlay
  @unpack_fields ics jtemp fmajor fminor jpress col_mix j_η
  @unpack_fields go ref flavor
  @unpack_fields ref press_log vmr press_log_delta temp temp_min temp_delta

  # input dimensions
  ftemp = Array{FT}(undef, ncol, nlay)
  fpress = Array{FT}(undef, ncol, nlay)

  @inbounds for ilay in 1:nlay
    @inbounds for icol in 1:ncol
      # index and factor for temperature interpolation
      jtemp[icol,ilay] = fint((t_lay[icol,ilay] - (temp_min - temp_delta)) / temp_delta)
      jtemp[icol,ilay] = min(ntemp - 1, max(1, jtemp[icol,ilay])) # limit the index range
      ftemp[icol,ilay] = (t_lay[icol,ilay] - temp[jtemp[icol,ilay]]) / temp_delta

      # index and factor for pressure interpolation
      locpress = FT(1) + (log(p_lay[icol,ilay]) - press_log[1]) / press_log_delta
      jpress[icol,ilay] = min(npres-1, max(1, fint(locpress)))
      fpress[icol,ilay] = locpress - FT(jpress[icol,ilay])
    end
  end

  @inbounds for ilay in 1:nlay
    @inbounds for icol in 1:ncol
      # itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
      itropo = fmerge(1,2,tropo[icol,ilay])
      # loop over implemented combinations of major species
      @inbounds for iflav in 1:nflav
        i_gases1 = flavor[1,iflav]+1
        i_gases2 = flavor[2,iflav]+1
        @inbounds for itemp in 1:2
          # compute interpolation fractions needed for lower, then upper reference temperature level
          # compute binary species parameter (η) for flavor and temperature and
          #  associated interpolation index and factors
          ij = jtemp[icol,ilay]+itemp-1
          vmr1 = vmr[itropo,i_gases1,ij]
          vmr2 = vmr[itropo,i_gases2,ij]
          ratio_η_half = vmr1 / vmr2
          col_mix[itemp,iflav,icol,ilay] = col_gas[icol,ilay,i_gases1] + ratio_η_half * col_gas[icol,ilay,i_gases2]
          η = fmerge(col_gas[icol,ilay,i_gases1] / col_mix[itemp,iflav,icol,ilay], FT(0.5),
                      col_mix[itemp,iflav,icol,ilay] > FT(2) * floatmin(FT))
          loc_η = η * FT(neta-1)
          j_η[itemp,iflav,icol,ilay] = min(fint(loc_η)+1, neta-1)
          f_η = mod(loc_η, 1)
          # compute interpolation fractions needed for minor species
          # ftemp_term = (FT(1)-ftemp(icol,ilay)) for itemp = 1, ftemp(icol,ilay) for itemp=2
          ftemp_term = (FT(2-itemp) + FT(2*itemp-3) * ftemp[icol,ilay])
          fminor[1,itemp,iflav,icol,ilay] = (1-f_η) * ftemp_term
          fminor[2,itemp,iflav,icol,ilay] =    f_η  * ftemp_term
          # compute interpolation fractions needed for major species
          fmajor[1,1,itemp,iflav,icol,ilay] = (FT(1)-fpress[icol,ilay]) * fminor[1,itemp,iflav,icol,ilay]
          fmajor[2,1,itemp,iflav,icol,ilay] = (FT(1)-fpress[icol,ilay]) * fminor[2,itemp,iflav,icol,ilay]
          fmajor[1,2,itemp,iflav,icol,ilay] =        fpress[icol,ilay]  * fminor[1,itemp,iflav,icol,ilay]
          fmajor[2,2,itemp,iflav,icol,ilay] =        fpress[icol,ilay]  * fminor[2,itemp,iflav,icol,ilay]
        end # reference temperatures
      end # iflav
    end # icol,ilay
  end
  return nothing
end
function interpolation!(ics::InterpolationCoefficientsPGP{FT,I},
                        go::AbstractGasOptics{FT,I},
                        as::AtmosphericStatePGP{FT,I}) where {I<:Int,B<:Bool,FT<:AbstractFloat}

  nflav = get_nflav(go)
  neta  = get_neta(go)
  npres = get_npres(go)
  ntemp = get_ntemp(go)

  @unpack_fields as col_gas t_lay p_lay tropo
  @unpack_fields ics jtemp fmajor fminor jpress col_mix j_η
  @unpack_fields go ref flavor
  @unpack_fields ref press_log vmr press_log_delta temp temp_min temp_delta

  # index and factor for temperature interpolation
  jtemp = fint((t_lay - (temp_min - temp_delta)) / temp_delta)
  jtemp = convert(I, min(ntemp - 1, max(1, jtemp))) # limit the index range
  ftemp = (t_lay - temp[jtemp]) / temp_delta

  # index and factor for pressure interpolation
  locpress = FT(1) + (log(p_lay) - press_log[1]) / press_log_delta
  jpress = convert(I, min(npres-1, max(1, fint(locpress))))
  fpress = locpress - FT(jpress)

  # determine if in lower or upper part of atmosphere

  # itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
  itropo = fmerge(1,2,tropo)
  # loop over implemented combinations of major species
  @inbounds for iflav in 1:nflav
    i_gases1 = flavor[1,iflav]+1
    i_gases2 = flavor[2,iflav]+1
    @inbounds for itemp in 1:2
      # compute interpolation fractions needed for lower, then upper reference temperature level
      # compute binary species parameter (η) for flavor and temperature and
      #  associated interpolation index and factors
      ij = jtemp+itemp-1
      vmr1 = vmr[itropo,i_gases1,ij]
      vmr2 = vmr[itropo,i_gases2,ij]
      ratio_η_half = vmr1 / vmr2
      col_mix[itemp,iflav] = col_gas[i_gases1] + ratio_η_half * col_gas[i_gases2]
      η = fmerge(col_gas[i_gases1] / col_mix[itemp,iflav], FT(0.5),
                  col_mix[itemp,iflav] > FT(2) * floatmin(FT))
      loc_η = η * FT(neta-1)
      j_η[itemp,iflav] = min(fint(loc_η)+1, neta-1)
      f_η = mod(loc_η, 1)
      # compute interpolation fractions needed for minor species
      # ftemp_term = (FT(1)-ftemp) for itemp = 1, ftemp for itemp=2
      ftemp_term = (FT(2-itemp) + FT(2*itemp-3) * ftemp)
      fminor[1,itemp,iflav] = (1-f_η) * ftemp_term
      fminor[2,itemp,iflav] =    f_η  * ftemp_term
      # compute interpolation fractions needed for major species
      fmajor[1,1,itemp,iflav] = (FT(1)-fpress) * fminor[1,itemp,iflav]
      fmajor[2,1,itemp,iflav] = (FT(1)-fpress) * fminor[2,itemp,iflav]
      fmajor[1,2,itemp,iflav] =        fpress  * fminor[1,itemp,iflav]
      fmajor[2,2,itemp,iflav] =        fpress  * fminor[2,itemp,iflav]
    end # reference temperatures
  end # iflav
  ics.jtemp = jtemp
  ics.jpress = jpress
  return nothing
end

"""
    compute_τ_absorption!(τ::Array{FT,3},
                          go::AbstractGasOptics{FT,I},
                          ics::InterpolationCoefficients{FT,I},
                          as::AtmosphericState{FT,I},
                          last_call=false) where {FT<:AbstractFloat,I<:Int}

Compute minor and major species optical depth

 - `τ` optical depth (ngpt,nlay,ncol)

given

 - `go` gas optics, see [`AbstractGasOptics`](@ref)
 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)
 - `as` atmospheric state, see [`AtmosphericState`](@ref)
"""
function compute_τ_absorption!(τ::Array{FT,3},
                               go::AbstractGasOptics{FT,I},
                               ics::InterpolationCoefficients{FT,I},
                               as::AtmosphericState{FT,I},
                               last_call=false) where {FT<:AbstractFloat,I<:Int}
  fill!(τ, FT(0))
  @timeit to_gor "go_major!" gas_optical_depths_major!(τ, go, as, ics) # Major Species
  @timeit to_gor "go_minor1!" gas_optical_depths_minor!(τ, go, as, ics, 1) # Minor Species - lower
  @timeit to_gor "go_minor2!" gas_optical_depths_minor!(τ, go, as, ics, 2) # Minor Species - upper
  return nothing
end
function compute_τ_absorption!(τ::Array{FT,1},
                               go::AbstractGasOptics{FT,I},
                               ics::InterpolationCoefficientsPGP{FT,I},
                               as::AtmosphericStatePGP{FT,I},
                               last_call=false) where {FT<:AbstractFloat,I<:Int}
  fill!(τ, FT(0))
  @timeit to_gor "go_major_pgp!" gas_optical_depths_major!(τ, go, as, ics) # Major Species
  @timeit to_gor "go_minor1_pgp!" gas_optical_depths_minor!(τ, go, as, ics, 1) # Minor Species - lower
  @timeit to_gor "go_minor2_pgp!" gas_optical_depths_minor!(τ, go, as, ics, 2) # Minor Species - upper
  return nothing
end


"""
    gas_optical_depths_major!(τ::Array{FT,3},
                              go::AbstractGasOptics{FT,I},
                              as::AtmosphericState{FT,I},
                              ics::InterpolationCoefficients{FT,I}) where {FT<:AbstractFloat,I<:Int}

compute minor species optical depths

 - `τ` optical depths (ngpt,nlay,ncol)

given

 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)
 - `go` gas optics, see [`AbstractGasOptics`](@ref)
"""
function gas_optical_depths_major!(τ::Array{FT,3},
                                   go::AbstractGasOptics{FT,I},
                                   as::AtmosphericState{FT,I},
                                   ics::InterpolationCoefficients{FT,I}) where {FT<:AbstractFloat,I<:Int}

  @unpack_fields as tropo
  @unpack_fields ics jtemp fmajor jpress col_mix j_η
  @unpack_fields go kmajor gpoint_flavor optical_props
  band_lims_gpt = get_band_lims_gpoint(optical_props)
  nbnd = get_nband(optical_props)
  ngpt,nlay,ncol = size(τ)

  @inbounds for icol in 1:ncol
    @inbounds for ilay in 1:nlay
      # itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
      itropo = fmerge(1,2,tropo[icol,ilay])
      # optical depth calculation for major species
      @inbounds for ibnd in 1:nbnd
        gptS = band_lims_gpt[1, ibnd]
        gptE = band_lims_gpt[2, ibnd]
        iflav = gpoint_flavor[itropo, gptS] #η interpolation depends on band's flavor
        # interpolation in temperature, pressure, and η

        fmajor_sm1 = SMatrix{2,2,FT}(fmajor[1,1,1,iflav,icol,ilay],
                                     fmajor[2,1,1,iflav,icol,ilay],
                                     fmajor[1,2,1,iflav,icol,ilay],
                                     fmajor[2,2,1,iflav,icol,ilay])
        fmajor_sm2 = SMatrix{2,2,FT}(fmajor[1,1,2,iflav,icol,ilay],
                                     fmajor[2,1,2,iflav,icol,ilay],
                                     fmajor[1,2,2,iflav,icol,ilay],
                                     fmajor[2,2,2,iflav,icol,ilay])

        j_η_sv = SVector{2,I}(j_η[1,iflav,icol,ilay],
                              j_η[2,iflav,icol,ilay])
        scaling = SVector{2,FT}(col_mix[1,iflav,icol,ilay],
                                col_mix[2,iflav,icol,ilay])

        @inbounds for igpt in gptS:gptE
          τ_major = interpolate3D(scaling,
                              fmajor_sm1,
                              fmajor_sm2,
                              kmajor,
                              j_η_sv,
                              igpt,
                              jtemp[icol,ilay],
                              jpress[icol,ilay]+itropo)
          τ[igpt,ilay,icol] += τ_major
        end
      end # igpt
    end
  end # ilay
  return nothing
end
function gas_optical_depths_major!(τ::Array{FT,1},
                                   go::AbstractGasOptics{FT,I},
                                   as::AtmosphericStatePGP{FT,I},
                                   ics::InterpolationCoefficientsPGP{FT,I}) where {FT<:AbstractFloat,I<:Int}

  @unpack_fields as tropo
  @unpack_fields ics jtemp fmajor jpress col_mix j_η
  @unpack_fields go kmajor gpoint_flavor optical_props
  band_lims_gpt = get_band_lims_gpoint(optical_props)
  nbnd = get_nband(optical_props)
  ngpt = length(τ)

  # itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
  itropo = fmerge(1,2,tropo)
  # optical depth calculation for major species
  @inbounds for ibnd in 1:nbnd
    gptS = band_lims_gpt[1, ibnd]
    gptE = band_lims_gpt[2, ibnd]
    iflav = gpoint_flavor[itropo, gptS] #η interpolation depends on band's flavor
    # interpolation in temperature, pressure, and η

    fmajor_sm1 = SMatrix{2,2,FT}(fmajor[1,1,1,iflav],
                                 fmajor[2,1,1,iflav],
                                 fmajor[1,2,1,iflav],
                                 fmajor[2,2,1,iflav])
    fmajor_sm2 = SMatrix{2,2,FT}(fmajor[1,1,2,iflav],
                                 fmajor[2,1,2,iflav],
                                 fmajor[1,2,2,iflav],
                                 fmajor[2,2,2,iflav])

    j_η_sv = SVector{2,I}(j_η[1,iflav],
                          j_η[2,iflav])
    scaling = SVector{2,FT}(col_mix[1,iflav],
                            col_mix[2,iflav])

    @inbounds for igpt in gptS:gptE
      τ_major = interpolate3D(scaling,
                          fmajor_sm1,
                          fmajor_sm2,
                          kmajor,
                          j_η_sv,
                          igpt,
                          jtemp,
                          jpress+itropo)
      τ[igpt] += τ_major
    end
  end # igpt
  return nothing
end


"""
    gas_optical_depths_minor!(τ::Array{FT,3},
                              go::AbstractGasOptics{FT,I},
                              as::AtmosphericState{FT,I},
                              ics::InterpolationCoefficients{FT,I},
                              itropo::I) where {FT<:AbstractFloat, I<:Int}

Compute minor species optical depths

 - `τ` optical depths (ngpt,nlay,ncol)

given

 - `go` gas optics, see [`AbstractGasOptics`](@ref)
 - `as` atmospheric state, see [`AtmosphericState`](@ref)
 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)
 - `itropo` integer indicating i-th level of troposphere

# Local variables
 - `vmr_fact`, `dry_fact` conversion from column abundance to dry vol. mixing ratio;
 - `scaling`, `kminor_loc` minor species absorption coefficient, optical depth
"""
function gas_optical_depths_minor!(τ::Array{FT,3},
                                   go::AbstractGasOptics{FT,I},
                                   as::AtmosphericState{FT,I},
                                   ics::InterpolationCoefficients{FT,I},
                                   itropo::I) where {FT<:AbstractFloat, I<:Int}
  @unpack_fields as col_gas t_lay p_lay tropo_lims gt_0_tropo_lims
  @unpack_fields ics jtemp fminor j_η
  @unpack_fields go gpoint_flavor lower lower_aux upper upper_aux gas_names
  idx_h2o = loc_in_array(h2o(), gas_names)
  atmos = itropo==1 ? lower : upper
  aux = itropo==1 ? lower_aux : upper_aux
  @unpack_fields atmos kminor kminor_start minor_scales_with_density minor_limits_gpt scale_by_complement
  @unpack_fields aux idx_minor idx_minor_scaling
  ngpt,nlay,ncol = size(τ)

  # number of minor contributors, total num absorption coeffs
  #
  # Guard against layer limits being 0 -- that means don't do anything i.e. there are no
  #   layers with pressures in the upper or lower atmosphere respectively
  # First check skips the routine entirely if all columns are out of bounds...
  #

  if gt_0_tropo_lims[itropo]
    @inbounds for imnr in 1:size(scale_by_complement,1) # loop over minor absorbers in each band
      kminor_start_imnr = kminor_start[imnr]
      @inbounds for icol in 1:ncol
        #
        # This check skips individual columns with no pressures in range
        #
        tropo_lims1 = tropo_lims[icol,1,itropo]
        tropo_lims2 = tropo_lims[icol,2,itropo]
        if tropo_lims1 > 0
          @inbounds for ilay in tropo_lims1:tropo_lims2
            #
            # Scaling of minor gas absorption coefficient begins with column amount of minor gas
            #
            scaling = col_gas[icol,ilay,idx_minor[imnr]+1]
            #
            # Density scaling (e.g. for h2o continuum, collision-induced absorption)
            #
            if minor_scales_with_density[imnr]
              #
              # NOTE: P needed in hPa to properly handle density scaling.
              #
              scaling = scaling * (PaTohPa(FT)*p_lay[icol,ilay]/t_lay[icol,ilay])
              if idx_minor_scaling[imnr] > 0  # there is a second gas that affects this gas's absorption
                vmr_fact = FT(1) / col_gas[icol,ilay,1]
                dry_fact = FT(1) / (FT(1) + col_gas[icol,ilay,idx_h2o+1] * vmr_fact)
                # scale by density of special gas
                if scale_by_complement[imnr] # scale by densities of all gases but the special one
                  scaling = scaling * (FT(1) - col_gas[icol,ilay,idx_minor_scaling[imnr]+1] * vmr_fact * dry_fact)
                else
                  scaling = scaling *          col_gas[icol,ilay,idx_minor_scaling[imnr]+1] * vmr_fact * dry_fact
                end
              end
            end
            #
            # Interpolation of absorption coefficient and calculation of optical depth
            #
            # Which g-point range does this minor gas affect?
            gptS = minor_limits_gpt[1,imnr]
            gptE = minor_limits_gpt[2,imnr]
            iflav = gpoint_flavor[itropo,gptS]
            fminor_sm = SMatrix{2,2,FT}(fminor[1,1,iflav,icol,ilay],
                                        fminor[2,1,iflav,icol,ilay],
                                        fminor[1,2,iflav,icol,ilay],
                                        fminor[2,2,iflav,icol,ilay])
            j_η_sv = SVector{2,I}(j_η[1,iflav,icol,ilay],
                                  j_η[2,iflav,icol,ilay])

            @inbounds for i_gpt in gptS:gptE
              igpt = kminor_start_imnr+i_gpt-gptS
              τ_minor = interpolate2D(fminor_sm, kminor, j_η_sv, jtemp[icol,ilay], igpt)
              τ[i_gpt,ilay,icol] += scaling*τ_minor
            end

          end
        end
      end
    end
  end
  return nothing
end
function gas_optical_depths_minor!(τ::Array{FT,1},
                                   go::AbstractGasOptics{FT,I},
                                   as::AtmosphericStatePGP{FT,I},
                                   ics::InterpolationCoefficientsPGP{FT,I},
                                   itropo::I) where {FT<:AbstractFloat, I<:Int}
  @unpack_fields as col_gas t_lay p_lay tropo_lims gt_0_tropo_lims
  @unpack_fields ics jtemp fminor j_η
  @unpack_fields go gpoint_flavor lower lower_aux upper upper_aux gas_names
  idx_h2o = loc_in_array(h2o(), gas_names)
  atmos = itropo==1 ? lower : upper
  aux = itropo==1 ? lower_aux : upper_aux
  @unpack_fields atmos kminor kminor_start minor_scales_with_density minor_limits_gpt scale_by_complement
  @unpack_fields aux idx_minor idx_minor_scaling
  ngpt = length(τ)

  # number of minor contributors, total num absorption coeffs
  #
  # Guard against layer limits being 0 -- that means don't do anything i.e. there are no
  #   layers with pressures in the upper or lower atmosphere respectively
  # First check skips the routine entirely if all columns are out of bounds...
  #

  if gt_0_tropo_lims[itropo]
    @inbounds for imnr in 1:size(scale_by_complement,1) # loop over minor absorbers in each band
      kminor_start_imnr = kminor_start[imnr]
      #
      # This check skips individual columns with no pressures in range
      #
      tropo_lims1 = tropo_lims[1,itropo]
      tropo_lims2 = tropo_lims[2,itropo]
      if tropo_lims1 > 0
        if as.ilay in tropo_lims1:tropo_lims2
          #
          # Scaling of minor gas absorption coefficient begins with column amount of minor gas
          #
          scaling = col_gas[idx_minor[imnr]+1]
          #
          # Density scaling (e.g. for h2o continuum, collision-induced absorption)
          #
          if minor_scales_with_density[imnr]
            #
            # NOTE: P needed in hPa to properly handle density scaling.
            #
            scaling = scaling * (PaTohPa(FT)*p_lay/t_lay)
            if idx_minor_scaling[imnr] > 0  # there is a second gas that affects this gas's absorption
              vmr_fact = FT(1) / col_gas[1]
              dry_fact = FT(1) / (FT(1) + col_gas[idx_h2o+1] * vmr_fact)
              # scale by density of special gas
              if scale_by_complement[imnr] # scale by densities of all gases but the special one
                scaling = scaling * (FT(1) - col_gas[idx_minor_scaling[imnr]+1] * vmr_fact * dry_fact)
              else
                scaling = scaling *          col_gas[idx_minor_scaling[imnr]+1] * vmr_fact * dry_fact
              end
            end
          end
          #
          # Interpolation of absorption coefficient and calculation of optical depth
          #
          # Which g-point range does this minor gas affect?
          gptS = minor_limits_gpt[1,imnr]
          gptE = minor_limits_gpt[2,imnr]
          iflav = gpoint_flavor[itropo,gptS]
          fminor_sm = SMatrix{2,2,FT}(fminor[1,1,iflav],
                                      fminor[2,1,iflav],
                                      fminor[1,2,iflav],
                                      fminor[2,2,iflav])
          j_η_sv = SVector{2,I}(j_η[1,iflav],
                                j_η[2,iflav])
          @inbounds for i_gpt in gptS:gptE
            igpt = kminor_start_imnr+i_gpt-gptS
            τ_minor = interpolate2D(fminor_sm, kminor, j_η_sv, jtemp, igpt)
            τ[i_gpt] += scaling*τ_minor
          end
        end
      end
    end
  end
  return nothing
end

"""
    compute_τ_Rayleigh!(τ_Rayleigh::Array{FT,3},
                        go::AbstractGasOptics{FT,I},
                        ics::InterpolationCoefficients{FT,I},
                        as::AtmosphericState{FT,I}) where {FT<:AbstractFloat, I<:Integer}

Compute Rayleigh scattering optical depths

 - `τ_Rayleigh` Rayleigh scattering optical depths

given

 - `go` gas optics, see [`AbstractGasOptics`](@ref)
 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)
 - `as` atmospheric state, see [`AtmosphericState`](@ref)

# Local variables
 - `k` Rayleigh scattering coefficient (ngpt)
"""
function compute_τ_Rayleigh!(τ_Rayleigh::Array{FT,3},
                             go::AbstractGasOptics{FT,I},
                             ics::InterpolationCoefficients{FT,I},
                             as::AtmosphericState{FT,I}) where {FT<:AbstractFloat, I<:Integer}

  @unpack_fields as col_gas col_dry tropo
  @unpack_fields ics jtemp fminor j_η
  @unpack_fields go krayl gpoint_flavor

  band_lims_gpt = get_band_lims_gpoint(go.optical_props)
  idx_h2o = loc_in_array(h2o(), go.gas_names)
  ngpt,nlay,ncol = size(τ_Rayleigh)

  nbnd = get_nband(go.optical_props)
  @inbounds for ilay in 1:nlay
    @inbounds for icol in 1:ncol
      itropo = fmerge(1,2,tropo[icol,ilay]) # itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
      @inbounds for ibnd in 1:nbnd
        gptS = band_lims_gpt[1, ibnd]
        gptE = band_lims_gpt[2, ibnd]
        iflav = gpoint_flavor[itropo, gptS] #η interpolation depends on band's flavor
        fminor_sm = SMatrix{2,2,FT}(fminor[1,1,iflav,icol,ilay],
                                    fminor[2,1,iflav,icol,ilay],
                                    fminor[1,2,iflav,icol,ilay],
                                    fminor[2,2,iflav,icol,ilay])
        j_η_sv = SVector{2,I}(j_η[1,iflav,icol,ilay],
                              j_η[2,iflav,icol,ilay])

        @inbounds for igpt in gptS:gptE
          k = interpolate2D(fminor_sm, krayl, j_η_sv, jtemp[icol,ilay], igpt, itropo)
          τ_Rayleigh[igpt,ilay,icol] = k * (col_gas[icol,ilay,idx_h2o+1]+col_dry[icol,ilay])
        end
      end
    end
  end
  return nothing
end
function compute_τ_Rayleigh!(τ_Rayleigh::Array{FT,1},
                             go::AbstractGasOptics{FT,I},
                             ics::InterpolationCoefficientsPGP{FT,I},
                             as::AtmosphericStatePGP{FT,I}) where {FT<:AbstractFloat, I<:Integer}

  @unpack_fields as col_gas col_dry tropo
  @unpack_fields ics jtemp fminor j_η
  @unpack_fields go krayl gpoint_flavor

  band_lims_gpt = get_band_lims_gpoint(go.optical_props)
  idx_h2o = loc_in_array(h2o(), go.gas_names)
  ngpt = length(τ_Rayleigh)

  nbnd = get_nband(go.optical_props)

  itropo = fmerge(1,2,tropo) # itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
  @inbounds for ibnd in 1:nbnd
    gptS = band_lims_gpt[1, ibnd]
    gptE = band_lims_gpt[2, ibnd]
    iflav = gpoint_flavor[itropo, gptS] #η interpolation depends on band's flavor
    fminor_sm = SMatrix{2,2,FT}(fminor[1,1,iflav],
                                fminor[2,1,iflav],
                                fminor[1,2,iflav],
                                fminor[2,2,iflav])
    j_η_sv = SVector{2,I}(j_η[1,iflav],
                          j_η[2,iflav])

    @inbounds for igpt in gptS:gptE
      k = interpolate2D(fminor_sm, krayl, j_η_sv, jtemp, igpt, itropo)
      τ_Rayleigh[igpt] = k * (col_gas[idx_h2o+1]+col_dry)
    end
  end

  return nothing
end

"""
    compute_Planck_source!(sources::SourceFuncLongWave{FT,I},
                           as::AtmosphericState{FT,I},
                           ics::InterpolationCoefficients{FT,I},
                           go::KDistributionLongwave{FT,I}) where {FT<:AbstractFloat,I<:Int}

Compute internal (Planck) source functions at layers and levels,
which depend on mapping from spectral space that creates k-distribution.

given

 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)
 - `go` gas optics, see [`KDistributionLongwave`](@ref)
 - `as` atmospheric state, see [`AtmosphericState`](@ref)

!!! Reduce allocations
"""
function compute_Planck_source!(sources::SourceFuncLongWave{FT,I},
                                as::AtmosphericState{FT,I},
                                ics::InterpolationCoefficients{FT,I},
                                go::KDistributionLongwave{FT,I}) where {FT<:AbstractFloat,I<:Int}

  @unpack_fields go ref totplnk_delta totplnk gpoint_flavor planck_frac optical_props
  @unpack_fields ref temp_min
  @unpack_fields sources sfc_source lay_source lev_source_inc lev_source_dec p_frac
  @unpack_fields as t_lay t_lev t_sfc mesh_orientation tropo
  @unpack_fields ics jtemp j_η jpress fmajor
  @unpack_fields mesh_orientation ilay_bot

  ncol  = get_ncol(sources)
  nlay  = get_nlay(sources)
  ngpt  = get_ngpt(sources)
  nbnd = get_nband(optical_props)
  band_lims_gpt = get_band_lims_gpoint(optical_props)
  gpoint_bands = get_gpoint_bands(optical_props)

  planck_function = Array{FT}(undef,nbnd,nlay+1,ncol)
  fill!(planck_function, FT(0))
  fill!(p_frac, FT(0))

  # Calculation of fraction of band's Planck irradiance associated with each g-point
  @inbounds for icol in 1:ncol
    @inbounds for ilay in 1:nlay
      # itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
      itropo = fmerge(1,2,tropo[icol,ilay])
      @inbounds for ibnd = 1:nbnd
        gptS = band_lims_gpt[1, ibnd]
        gptE = band_lims_gpt[2, ibnd]
        iflav = gpoint_flavor[itropo, gptS] #η interpolation depends on band's flavor

        fmajor_sm1 = SMatrix{2,2,FT}(fmajor[1,1,1,iflav,icol,ilay],
                                     fmajor[2,1,1,iflav,icol,ilay],
                                     fmajor[1,2,1,iflav,icol,ilay],
                                     fmajor[2,2,1,iflav,icol,ilay])
        fmajor_sm2 = SMatrix{2,2,FT}(fmajor[1,1,2,iflav,icol,ilay],
                                     fmajor[2,1,2,iflav,icol,ilay],
                                     fmajor[1,2,2,iflav,icol,ilay],
                                     fmajor[2,2,2,iflav,icol,ilay])

        j_η_sv = SVector{2,I}(j_η[1,iflav,icol,ilay],
                              j_η[2,iflav,icol,ilay])
        scaling = SVector{2,FT}(1,1)

        # interpolation in temperature, pressure, and η
        @inbounds for igpt in gpt_range(optical_props, ibnd)
          τ_local = interpolate3D(scaling,
                              fmajor_sm1,
                              fmajor_sm2,
                              planck_frac,
                              j_η_sv,
                              igpt,
                              jtemp[icol,ilay],
                              jpress[icol,ilay]+itropo)
          p_frac[icol,ilay,igpt] = τ_local
        end
      end # band
    end   # layer
  end     # column

  #
  # Planck function by band for the surface
  # Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
  #
  @inbounds for icol in 1:ncol
    planck_function[1:nbnd,1,icol] .= interpolate1D(t_sfc[icol], temp_min, totplnk_delta, totplnk)

    # Map to g-points
    @inbounds for ibnd in 1:nbnd
      @inbounds for igpt in gpt_range(optical_props, ibnd)
        sfc_source[icol, igpt] = p_frac[icol,ilay_bot,igpt] * planck_function[ibnd, 1, icol]
      end
    end
  end # icol

  @inbounds for icol in 1:ncol
    @inbounds for ilay in 1:nlay
      # Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
      planck_function[1:nbnd,ilay,icol] .= interpolate1D(t_lay[icol,ilay], temp_min, totplnk_delta, totplnk)

      # Map to g-points
      @inbounds for ibnd in 1:nbnd
        @inbounds for igpt in gpt_range(optical_props, ibnd)
          lay_source[icol,ilay,igpt] = p_frac[icol,ilay,igpt] * planck_function[ibnd,ilay,icol]
        end
      end
    end # ilay
  end # icol

  # compute level source irradiances for each g-point, one each for upward and downward paths
  @inbounds for icol in 1:ncol
    planck_function[1:nbnd,       1,icol] .= interpolate1D(t_lev[icol,     1], temp_min, totplnk_delta, totplnk)
    @inbounds for ilay in 1:nlay
      planck_function[1:nbnd,ilay+1,icol] .= interpolate1D(t_lev[icol,ilay+1], temp_min, totplnk_delta, totplnk)

      # Map to g-points
      @inbounds for ibnd in 1:nbnd
        @inbounds for igpt in gpt_range(optical_props, ibnd)
          lev_source_inc[icol,ilay,igpt] = p_frac[icol,ilay,igpt] * planck_function[ibnd,ilay+1,icol]
          lev_source_dec[icol,ilay,igpt] = p_frac[icol,ilay,igpt] * planck_function[ibnd,ilay,  icol]
        end
      end
    end # ilay
  end # icol
  return nothing
end
function compute_Planck_source!(sources::SourceFuncLongWavePGP{FT,I},
                                as::AtmosphericStatePGP{FT,I},
                                ics::InterpolationCoefficientsPGP{FT,I},
                                go::KDistributionLongwave{FT,I}) where {FT<:AbstractFloat,I<:Int}

  @unpack_fields go ref totplnk_delta totplnk gpoint_flavor planck_frac optical_props
  @unpack_fields ref temp_min
  @unpack_fields sources sfc_source lay_source lev_source_inc lev_source_dec p_frac
  @unpack_fields as t_lay t_lev t_sfc tropo
  @unpack_fields ics jtemp j_η jpress fmajor

  ngpt  = get_ngpt(sources)
  nbnd = get_nband(optical_props)
  band_lims_gpt = get_band_lims_gpoint(optical_props)
  gpoint_bands = get_gpoint_bands(optical_props)

  planck_function = Array{FT}(undef,nbnd,2)
  fill!(planck_function, FT(0))
  fill!(p_frac, FT(0))

  # Calculation of fraction of band's Planck irradiance associated with each g-point
  # itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
  itropo = fmerge(1,2,tropo)
  @inbounds for ibnd = 1:nbnd
    gptS = band_lims_gpt[1, ibnd]
    gptE = band_lims_gpt[2, ibnd]
    iflav = gpoint_flavor[itropo, gptS] #η interpolation depends on band's flavor

    fmajor_sm1 = SMatrix{2,2,FT}(fmajor[1,1,1,iflav],
                                 fmajor[2,1,1,iflav],
                                 fmajor[1,2,1,iflav],
                                 fmajor[2,2,1,iflav])
    fmajor_sm2 = SMatrix{2,2,FT}(fmajor[1,1,2,iflav],
                                 fmajor[2,1,2,iflav],
                                 fmajor[1,2,2,iflav],
                                 fmajor[2,2,2,iflav])

    j_η_sv = SVector{2,I}(j_η[1,iflav],
                          j_η[2,iflav])
    scaling = SVector{2,FT}(1,1)

    # interpolation in temperature, pressure, and η
    @inbounds for igpt in gptS:gptE
      τ_local = interpolate3D(scaling,
                          fmajor_sm1,
                          fmajor_sm2,
                          planck_frac,
                          j_η_sv,
                          igpt,
                          jtemp,
                          jpress+itropo)
      p_frac[igpt] = τ_local
    end
  end # band

  # Planck function by band for the surface
  # Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
  planck_function[1:nbnd,1] .= interpolate1D(t_sfc, temp_min, totplnk_delta, totplnk)

  # Map to g-points
  @inbounds for ibnd in 1:nbnd
    @inbounds for igpt in gpt_range(optical_props, ibnd)
      sfc_source[igpt] = p_frac[igpt] * planck_function[ibnd, 1]
    end
  end

  # Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
  planck_function[1:nbnd,1] .= interpolate1D(t_lay, temp_min, totplnk_delta, totplnk)

  # Map to g-points
  @inbounds for ibnd in 1:nbnd
    @inbounds for igpt in gpt_range(optical_props, ibnd)
      lay_source[igpt] = p_frac[igpt] * planck_function[ibnd,1]
    end
  end

  # compute level source irradiances for each g-point, one each for upward and downward paths
  planck_function[1:nbnd,1] .= interpolate1D(t_lev[1], temp_min, totplnk_delta, totplnk)
  planck_function[1:nbnd,2] .= interpolate1D(t_lev[2], temp_min, totplnk_delta, totplnk)

  # Map to g-points
  @inbounds for ibnd in 1:nbnd
    @inbounds for igpt in gpt_range(optical_props, ibnd)
      lev_source_inc[igpt] = p_frac[igpt] * planck_function[ibnd,2]
      lev_source_dec[igpt] = p_frac[igpt] * planck_function[ibnd,1]
    end
  end
  return nothing
end

"""
    interpolate1D(val::FT,
                  offset::FT,
                  delta::FT,
                  table::Array{FT,2}) where {FT<:AbstractFloat}

One dimensional interpolation -- return all values along second table dimension

given

 - `val` axis value at which to evaluate table
 - `offset` minimum of table axis
 - `delta` step size of table axis
 - `table` dimensions (axis, values)

# local
 - `val0` fraction index adjusted by offset and delta
 - `index` index term
 - `frac` fractional term
"""
function interpolate1D(val::FT,
                       offset::FT,
                       delta::FT,
                       table::Array{FT,2}) where {FT<:AbstractFloat}
  res = Vector{FT}(undef, size(table, 2))
  val0 = (val - offset) / delta
  f0 = fint(val0)
  frac = val0 - f0 # get fractional part
  index = Integer(min(size(table,1)-1, max(1, f0+1))) # limit the index range
  @inbounds for i in 1:size(table, 2)
    res[i] = table[index,i] + frac * (table[index+1,i] - table[index,i])
  end
  return res
end

"""
    interpolate2D(fminor::SMatrix{2,2,FT},
                  k::Array{FT,4},
                  j_η::SVector{2,I},
                  jtemp::I,
                  igpt::I[,
                  itropo::I]) where {I<:Int,FT<:AbstractFloat}

A single value from a subset (in gpoint) of the k table
"""
function interpolate2D(fminor::SMatrix{2,2,FT},
                       k::Array{FT,3},
                       j_η::SVector{2,I},
                       jtemp::I,
                       igpt::I) where {I<:Int,FT<:AbstractFloat}
  return fminor[1,1] * k[igpt, j_η[1]  , jtemp  ] +
         fminor[2,1] * k[igpt, j_η[1]+1, jtemp  ] +
         fminor[1,2] * k[igpt, j_η[2]  , jtemp+1] +
         fminor[2,2] * k[igpt, j_η[2]+1, jtemp+1]
end
function interpolate2D(fminor::SMatrix{2,2,FT},
                       k::Array{FT,4},
                       j_η::SVector{2,I},
                       jtemp::I,
                       igpt::I,
                       itropo::I) where {I<:Int,FT<:AbstractFloat}
  return fminor[1,1] * k[igpt, j_η[1]  , jtemp  ,itropo] +
         fminor[2,1] * k[igpt, j_η[1]+1, jtemp  ,itropo] +
         fminor[1,2] * k[igpt, j_η[2]  , jtemp+1,itropo] +
         fminor[2,2] * k[igpt, j_η[2]+1, jtemp+1,itropo]
end

"""
    interpolate3D(scaling::SVector{2,FT},
                  fmajor1::SMatrix{2,2,FT},
                  fmajor2::SMatrix{2,2,FT},
                  k::Array{FT,4},
                  j_η::SVector{2,I},
                  igpt::I,
                  jtemp::I,
                  jpress::I) where {I<:Int, FT<:AbstractFloat}

Interpolation in temperature, pressure, and η

 - `scaling` scaling
 - `fmajor1` interpolation fractions for major species for lower reference temperature level
             index(1) : reference η level (temperature dependent)
             index(2) : reference pressure level
 - `fmajor2` interpolation fractions for major species for upper reference temperature level
             index(1) : reference η level (temperature dependent)
             index(2) : reference pressure level
 - `k` absorption coefficient (gpt, η,temp,press)
 - `igpt` i-th index for g-point
 - `j_η` interpolation indexes for binary species parameter (η)
 - `jtemp`  interpolation index for temperature
 - `jpress` interpolation index for pressure
"""
function interpolate3D(scaling::SVector{2,FT},
                       fmajor1::SMatrix{2,2,FT},
                       fmajor2::SMatrix{2,2,FT},
                       k::Array{FT,4},
                       j_η::SVector{2,I},
                       igpt::I,
                       jtemp::I,
                       jpress::I) where {I<:Int, FT<:AbstractFloat}
  # each code block is for a different reference temperature
  return scaling[1] *
    ( fmajor1[1,1] * k[igpt, j_η[1]  , jpress-1, jtemp  ] +
      fmajor1[2,1] * k[igpt, j_η[1]+1, jpress-1, jtemp  ] +
      fmajor1[1,2] * k[igpt, j_η[1]  , jpress  , jtemp  ] +
      fmajor1[2,2] * k[igpt, j_η[1]+1, jpress  , jtemp  ] ) +
    scaling[2] *
    ( fmajor2[1,1] * k[igpt, j_η[2]  , jpress-1, jtemp+1] +
      fmajor2[2,1] * k[igpt, j_η[2]+1, jpress-1, jtemp+1] +
      fmajor2[1,2] * k[igpt, j_η[2]  , jpress  , jtemp+1] +
      fmajor2[2,2] * k[igpt, j_η[2]+1, jpress  , jtemp+1] )
end

"""
    combine_and_reorder!(op::AbstractOpticalPropsArry{FT},
                         τ_abs::Array{FT,3},
                         τ_Rayleigh::Array{FT,3},
                         has_Rayleigh::Bool) where {FT<:AbstractFloat}

Combine absorption and Rayleigh optical depths for total `τ`, `ssa`, `g`

!!! Reduce allocations
"""
function combine_and_reorder!(op::TwoStream{FT},
                              τ_abs::Array{FT,3},
                              τ_Rayleigh::Array{FT,3},
                              has_Rayleigh::Bool) where {FT<:AbstractFloat}
  fill!(op.g, FT(0))
  if has_Rayleigh
    @inbounds for icol in 1:size(op.τ,1)
      @inbounds for ilay in 1:size(op.τ,2)
        @inbounds for igpt in 1:get_ngpt(op)
          τ_Rayleigh′ = τ_Rayleigh[igpt,ilay,icol]
          τ_abs′ = τ_abs[igpt,ilay,icol]
          τ = τ_abs′+τ_Rayleigh′
          op.τ[icol,ilay,igpt] = τ
          op.ssa[icol,ilay,igpt] = τ > FT(2) * floatmin(FT) ? τ_Rayleigh′/τ : FT(0)
        end
      end
    end
  else
    permutedims!(op.τ, τ_abs, [3,2,1])
    fill!(op.ssa, FT(0))
  end
  return nothing
end
function combine_and_reorder!(op::OneScalar{FT},
                              τ::Array{FT,3},
                              τ_Rayleigh::Array{FT,3},
                              has_Rayleigh::Bool) where FT
  permutedims!(op.τ, τ, [3,2,1])
end
function combine_and_reorder!(op::TwoStreamPGP{FT},
                              τ_abs::Array{FT,1},
                              τ_Rayleigh::Array{FT,1},
                              has_Rayleigh::Bool) where {FT<:AbstractFloat}
  fill!(op.g, FT(0))
  if has_Rayleigh
    @inbounds for igpt in 1:get_ngpt(op)
      τ_Rayleigh′ = τ_Rayleigh[igpt]
      τ_abs′ = τ_abs[igpt]
      τ = τ_abs′+τ_Rayleigh′
      op.τ[igpt] = τ
      op.ssa[igpt] = τ > FT(2) * floatmin(FT) ? τ_Rayleigh′/τ : FT(0)
    end
  else
    op.τ .= τ_abs
    fill!(op.ssa, FT(0))
  end
  return nothing
end
function combine_and_reorder!(op::OneScalarPGP{FT},
                              τ::Array{FT,1},
                              τ_Rayleigh::Array{FT,1},
                              has_Rayleigh::Bool) where FT
  op.τ .= τ
end
