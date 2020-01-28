#=
Numeric calculations for radiative transfer solvers.
 Emission/absorption (no-scattering) calculations
   solver for multi-angle Gaussian quadrature
   solver for a single angle, calling
     source function computation (linear-in-τ)
     transport
 Extinction-only calculation (direct solar beam)
 Two-stream calculations
   solvers for LW and SW with different boundary conditions and source functions
     source function calculation for LW, SW
     two-stream calculations for LW, SW (using different assumptions about phase function)
     transport (adding)
 Application of boundary conditions
=#

LW_diff_sec(::Type{FT}) where FT = FT(1.66)  # 1./cos(diffusivity angle)

#####
##### Top-level long-wave kernels
#####

#
# LW fluxes, no scattering, μ (cosine of integration angle) specified by column
#  Does radiation calculation at user-supplied angles; converts radiances to flux
#  using user-supplied weights
#

"""
    lw_solver_noscat!(mo::MeshOrientation{I},
                      D::Array{FT,2},
                      w_μ::FT,
                      op::OneScalar{FT,I},
                      source::SourceFuncLongWave{FT,I},
                      sfc_emis::Array{FT,2},
                      radn_up::Array{FT,3},
                      radn_dn::Array{FT,3}) where {FT<:AbstractFloat,I<:Int}

 - `source` longwave source function, see [`SourceFuncLongWave`](@ref)
 - `mo` mesh orientation, see [`MeshOrientation`](@ref)
 - `op` optical properties, see [`OneScalar`](@ref)
 - `sfc_emis` - surface emissivity
 - `radn_up` - upward radiance [W/m2-str]
 - `radn_dn` - downward radiance [W/m2-str] Top level must contain incident flux boundary condition

 - `D` secant of propagation angle  []
 - `w_μ` quadrature weight

# Local variables, no g-point dependency
real(FT), dimension(ncol,nlay) :: τ_loc,   # path length (τ/μ)
                                  trans    # transmissivity  = exp(-τ)
real(FT), dimension(ncol,nlay) :: source_dn, source_up
real(FT), dimension(ncol     ) :: source_sfc, sfc_albedo

real(FT), dimension(:,:,:), pointer :: lev_source_up, lev_source_dn # Mapping increasing/decreasing indicies to up/down
"""
function lw_solver_noscat!(mo::MeshOrientation{I},
                           D::Array{FT,2},
                           w_μ::FT,
                           op::OneScalar{FT,I},
                           source::SourceFuncLongWave{FT,I},
                           sfc_emis::Array{FT,2},
                           radn_up::Array{FT,3},
                           radn_dn::Array{FT,3}) where {FT<:AbstractFloat,I<:Int}
  @unpack_fields op τ
  @unpack_fields source lay_source sfc_source
  ncol  = get_ncol(op)
  nlay  = get_nlay(op)
  ngpt  = get_ngpt(op)

  τ_loc    = Array{FT}(undef,ncol,nlay)
  trans      = Array{FT}(undef,ncol,nlay)
  source_up  = Array{FT}(undef,ncol,nlay)
  source_dn  = Array{FT}(undef,ncol,nlay)
  source_sfc = Array{FT}(undef,ncol)
  sfc_albedo = Array{FT}(undef,ncol)
  # Which way is up?
  # Level Planck sources for upward and downward radiation
  # When top_at_1, lev_source_up => lev_source_dec
  #                lev_source_dn => lev_source_inc, and vice-versa
  i_lev_top = ilev_top(mo)
  lev_source_up = mo.top_at_1 ? source.lev_source_dec : source.lev_source_inc
  lev_source_dn = mo.top_at_1 ? source.lev_source_inc : source.lev_source_dec

  @inbounds for igpt in 1:ngpt
    #
    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    #
    radn_dn[:,i_lev_top,igpt] ./= 2 * FT(π) * w_μ

    #
    # Optical path and transmission, used in source function and transport calculations
    #
    @inbounds for ilev in 1:nlay
      τ_loc[:,ilev] .= τ[:,ilev,igpt] .* D[:,igpt]
      trans[:,ilev] .= exp.(-τ_loc[:,ilev])
    end
    #
    # Source function for diffuse radiation
    #
    @inbounds for ilay in 1:nlay
      @inbounds for icol in 1:ncol
        source_dn[icol,ilay], source_up[icol,ilay] =
          lw_source_noscat(lay_source[icol,ilay,igpt],
                           lev_source_up[icol,ilay,igpt],
                           lev_source_dn[icol,ilay,igpt],
                           τ_loc[icol,ilay],
                           trans[icol,ilay])
      end
    end

    #
    # Surface albedo, surface source function
    #
    sfc_albedo .= 1 .- sfc_emis[:,igpt]
    source_sfc .= sfc_emis[:,igpt] .* sfc_source[:,igpt]
    #
    # Transport
    #
    lw_transport_noscat!(ncol, nlay, mo,
                         trans,
                         sfc_albedo,
                         source_dn,
                         source_up,
                         source_sfc,
                         @view(radn_up[:,:,igpt]),
                         @view(radn_dn[:,:,igpt]))
    #
    # Convert intensity to flux assuming azimuthal isotropy and quadrature weight
    #
    radn_dn[:,:,igpt] .*= 2 * FT(π) * w_μ
    radn_up[:,:,igpt] .*= 2 * FT(π) * w_μ
  end  # g point loop
  return nothing

end

"""
    solve!(rte::RTELongWaveNoScattering{FT,I}, op::OneScalar{FT,I}, mo::MeshOrientation{I}) where {FT,I}

Longwave transport, no scattering, multi-angle quadrature
Users provide a set of weights and quadrature angles
Routine sums over single-angle solutions for each sets of angles/weights

given

 - `rte` radiative transfer equation, see [`RTELongWaveNoScattering`](@ref)
 - `op` optical properties, see [`OneScalar`](@ref)
 - `mo` mesh orientation, see [`MeshOrientation`](@ref)
"""
function solve!(rte::RTELongWaveNoScattering{FT,I}, op::OneScalar{FT,I}, mo::MeshOrientation{I}) where {FT,I}

  # Upper boundary condition
  angle_disc = rte.angle_disc
  source = rte.sources
  sfc_emis = rte.sfc_emis_gpt
  flux_up = rte.base.gpt_flux_up
  flux_dn = rte.base.gpt_flux_dn
  i_lev_top = ilev_top(mo)
  inc_flux = rte.base.bcs.inc_flux
  apply_BC!(flux_dn, i_lev_top, inc_flux)

  ncol  = get_ncol(op)
  nlay  = get_nlay(op)
  ngpt  = get_ngpt(op)

  n_μ = angle_disc.n_gauss_angles
  Ds = angle_disc.gauss_Ds[1:n_μ,n_μ]
  w_μ = angle_disc.gauss_wts[1:n_μ,n_μ]

  # For the first angle output arrays store total flux
  Ds_ncol = Array{FT}(undef, ncol,ngpt)
  radn_up, radn_dn = ntuple(i->zeros(FT, ncol,nlay+1,ngpt), 2)

  @inbounds for i_μ in 1:n_μ
    i_μ==1 && apply_BC!(radn_dn, i_lev_top, inc_flux)
    Ds_ncol .= Ds[i_μ]
    lw_solver_noscat!(mo,
                      Ds_ncol,
                      w_μ[i_μ],
                      op,
                      source,
                      sfc_emis,
                      radn_up,
                      radn_dn)
    flux_up .+= radn_up
    flux_dn .+= radn_dn
  end
  return nothing
end

"""
    solve!(rte::RTELongWave{FT,I}, op::TwoStream{FT,I}, mo::MeshOrientation{I}) where {FT,I}

Longwave calculation:
 - combine RRTMGP-specific sources at levels
 - compute layer reflectance, transmittance
 - compute total source function at levels using linear-in-τ
 - transport

given

 - `rte` radiative transfer equation, see [`RTELongWave`](@ref)
 - `op` optical properties, see [`TwoStream`](@ref)
 - `mo` mesh orientation, see [`MeshOrientation`](@ref)

 - `Rdif`, `Tdif`, `γ_1`, `γ_2`
 - `sfc_albedo`
 - `lev_source`
 - `source_dn`, `source_up`
 - `source_sfc`
"""
function solve!(rte::RTELongWave{FT,I}, op::TwoStream{FT,I}, mo::MeshOrientation{I}) where {FT,I}
  @unpack_fields rte base sfc_emis_gpt sources
  @unpack_fields base gpt_flux_dn gpt_flux_up
  inc_flux = rte.base.bcs.inc_flux
  apply_BC!(gpt_flux_dn, ilev_top(mo), inc_flux)

  source = sources
  sfc_emis = sfc_emis_gpt
  flux_up = gpt_flux_up
  flux_dn = gpt_flux_dn

  @unpack_fields op τ ssa g
  @unpack_fields source lev_source_inc lev_source_dec sfc_source lay_source
  ncol  = get_ncol(op)
  nlay  = get_nlay(op)
  ngpt  = get_ngpt(op)

  lev_source = Array{FT}(undef, ncol,nlay+1)
  source_sfc = Array{FT}(undef, ncol)
  source_dn = Array{FT}(undef, ncol, nlay)
  source_up = Array{FT}(undef, ncol, nlay)
  sfc_albedo = Array{FT}(undef, ncol)
  Rdif, Tdif, γ_1, γ_2 = ntuple(i->Array{FT}(undef, ncol, nlay),4)
  @inbounds for igpt in 1:ngpt

    # RRTMGP provides source functions at each level using the spectral mapping
    #   of each adjacent layer. Combine these for two-stream calculations
    lw_combine_sources!(ncol, nlay,
                        lev_source_inc[:,:,igpt],
                        lev_source_dec[:,:,igpt],
                        lev_source)

    # Cell properties: reflection, transmission for diffuse radiation
    #   Coupling coefficients needed for source function
    @inbounds for ilay in 1:nlay
      @inbounds for icol in 1:ncol
        γ_1[icol,ilay],
        γ_2[icol,ilay],
        Rdif[icol,ilay],
        Tdif[icol,ilay] =
          lw_two_stream(τ[icol,ilay,igpt],
                        ssa[icol,ilay,igpt],
                        g[icol,ilay,igpt])
      end
    end

    # Source function for diffuse radiation
    lw_source_2str!(ncol, nlay, mo,
                    sfc_emis[:,igpt],
                    sfc_source[:,igpt],
                    lay_source[:,:,igpt],
                    lev_source,
                    γ_1,
                    γ_2,
                    Rdif,
                    Tdif,
                    τ[:,:,igpt],
                    source_dn,
                    source_up,
                    source_sfc)

    # Transport
    sfc_albedo .= FT(1) .- sfc_emis[:,igpt]
    adding!(ncol, nlay,
            mo,
            sfc_albedo,
            Rdif,
            Tdif,
            source_dn,
            source_up,
            source_sfc,
            @view(flux_up[:,:,igpt]),
            @view(flux_dn[:,:,igpt]))
  end
  return nothing

end

#####
#####   Top-level shortwave kernels
#####

"""
    solve!(rte::RTEShortWaveNoScattering{FT,I}, op::OneScalar{FT,I}, mo::MeshOrientation) where {FT,I}

Solve the RTE with extinction-only i.e. solar direct beam.

 - `rte` radiative transfer equation, see [`RTEShortWaveNoScattering`](@ref)
 - `op` optical properties, see [`OneScalar`](@ref)
 - `mo` mesh orientation, see [`MeshOrientation`](@ref)
"""
function solve!(rte::RTEShortWaveNoScattering{FT,I}, op::OneScalar{FT,I}, mo::MeshOrientation) where {FT,I}

  @unpack_fields rte base μ_0
  @unpack_fields base gpt_flux_dn bcs gpt_flux_dir
  flux_dir = gpt_flux_dir
  i_lev_top = ilev_top(mo)
  apply_BC!(gpt_flux_dir, i_lev_top, bcs.toa_flux, μ_0)
  apply_BC!(gpt_flux_dn, i_lev_top, bcs.inc_flux_diffuse)

  @unpack_fields op τ
  ngpt  = get_ngpt(op)

  μ_0_inv = 1 ./ μ_0
  # Indexing into arrays for upward and downward propagation depends on orientation
  n̂ = nhat(mo)
  b = binary(mo)

  # Downward propagation
  # For the flux at this level, what was the previous level, and which layer has the
  #   radiation just passed through?
  @inbounds for igpt in 1:ngpt
    @inbounds for ilev in lev_range(mo)
      flux_dir[:,ilev,igpt] .= flux_dir[:,ilev-n̂,igpt] .* exp.(-τ[:,ilev-b,igpt].*μ_0_inv)
    end
  end
  return nothing
end


"""
    solve!(rte::RTEShortWave{FT,I}, op::TwoStream{FT,I}, mo::MeshOrientation{I}) where {FT,I}

Shortwave two-stream calculation:
  compute layer reflectance, transmittance
  compute solar source function for diffuse radiation
  transport

 - `rte` radiative transfer equation, see [`RTEShortWave`](@ref)
 - `op` optical properties, see [`TwoStream`](@ref)
 - `mo` mesh orientation, see [`MeshOrientation`](@ref)

 - `Rdif`, `Tdif`, `Rdir`, `Tdir`, `Tnoscat`
 - `source_up`, `source_dn`
 - `source_srf`
"""
function solve!(rte::RTEShortWave{FT,I}, op::TwoStream{FT,I}, mo::MeshOrientation{I}) where {FT,I}
  @unpack_fields rte base
  @unpack_fields base bcs
  μ_0 = rte.μ_0
  sfc_alb_dir = rte.sfc_alb_dir_gpt
  sfc_alb_dif = rte.sfc_alb_dif_gpt
  flux_up = base.gpt_flux_up
  flux_dn = base.gpt_flux_dn
  flux_dir = base.gpt_flux_dir
  i_lev_top = ilev_top(mo)

  apply_BC!(flux_dir, i_lev_top, bcs.toa_flux, rte.μ_0)
  apply_BC!(flux_dn, i_lev_top, bcs.inc_flux_diffuse)


  @unpack_fields op τ ssa g
  ncol  = get_ncol(op)
  nlay  = get_nlay(op)
  ngpt  = get_ngpt(op)

  Rdif, Tdif, Rdir, Tdir, Tnoscat, source_up, source_dn = ntuple(i->zeros(FT, ncol,nlay), 7)
  source_srf = zeros(FT, ncol)

  @inbounds for igpt in 1:ngpt

    # Cell properties: transmittance and reflectance for direct and diffuse radiation
    @inbounds for ilay in 1:nlay
      @inbounds for icol in 1:ncol
        Rdif[icol,ilay],
        Tdif[icol,ilay],
        Rdir[icol,ilay],
        Tdir[icol,ilay],
        Tnoscat[icol,ilay] =
          sw_two_stream(μ_0[icol],
                        τ[icol,ilay,igpt],
                        ssa[icol,ilay,igpt],
                        g[icol,ilay,igpt])
      end
    end

    # Direct-beam and source for diffuse radiation
    sw_source_2str!(ncol, nlay, mo, Rdir, Tdir, Tnoscat, sfc_alb_dir[:,igpt],
                    source_up, source_dn, source_srf, @view(flux_dir[:,:,igpt]))

    # Transport
    adding!(ncol, nlay,
            mo,
            sfc_alb_dif[:,igpt],
            Rdif,
            Tdif,
            source_dn,
            source_up,
            source_srf,
            @view(flux_up[:,:,igpt]),
            @view(flux_dn[:,:,igpt]))

    # adding computes only diffuse flux; flux_dn is total
    flux_dn[:,:,igpt] .= flux_dn[:,:,igpt] .+ flux_dir[:,:,igpt]
  end
  return nothing
end

#####
#####   Lower-level longwave kernels
#####

"""
    lw_source_noscat(lay_source::FT,
                     lev_source_up::FT,
                     lev_source_dn::FT,
                     τ::FT,
                     trans::FT) where {I<:Int,FT<:AbstractFloat}

Compute LW source function for upward and downward emission at levels using linear-in-τ assumption
See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13

 - `lay_source` Planck source at layer center
 - `lev_source_up` Planck source at levels (layer edges), increasing/decreasing layer index
 - `lev_source_dn` Planck source at levels (layer edges), increasing/decreasing layer index
 - `τ` Optical path (τ/μ)
 - `trans` Transmissivity (exp(-τ))
 - `source_dn` source dn function at upper
 - `source_up` source up function at lower
"""
function lw_source_noscat(lay_source::FT,
                          lev_source_up::FT,
                          lev_source_dn::FT,
                          τ::FT,
                          trans::FT) where {I<:Int,FT<:AbstractFloat}

  τ_thresh = sqrt(eps(FT)) # or abs(eps(FT))?

  # Weighting factor. Use 2nd order series expansion when rounding error (~τ^2)
  #   is of order epsilon (smallest difference from 1. in working precision)
  #   Thanks to Peter Blossey
  B̄ = lay_source
  B_D = lev_source_dn
  B_U = lev_source_up
  trans = trans
  fact = τ > τ_thresh ? (1 - trans)/τ - trans : τ * (FT(1/2) - 1/3*τ)

  # Equation below is developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
  source_dn = (1 - trans) * B_D + 2 * fact * (B̄ - B_D)
  source_up = (1 - trans) * B_U + 2 * fact * (B̄ - B_U)

  return source_dn, source_up
end

#####
##### Longwave no-scattering transport
#####

"""
    lw_transport_noscat!(ncol::I, nlay::I,
                         mo::MeshOrientation{I},
                         τ::Array{FT,2},
                         trans::Array{FT,2},
                         sfc_albedo::Array{FT,1},
                         source_dn::Array{FT,2},
                         source_up::Array{FT,2},
                         source_sfc::Array{FT,1},
                         radn_up::AbstractArray{FT,2},
                         radn_dn::AbstractArray{FT,2}) where {I<:Int, FT<:AbstractFloat}

 - `mo` mesh orientation, see [`MeshOrientation`](@ref)

integer,                          intent(in   ) :: ncol, nlay # Number of columns, layers, g-points
real(FT), dimension(ncol,nlay  ), intent(in   ) :: τ,         # Absorption optical thickness, pre-divided by μ []
                                                   trans      # transmissivity = exp(-τ)
real(FT), dimension(ncol       ), intent(in   ) :: sfc_albedo # Surface albedo
real(FT), dimension(ncol,nlay  ), intent(in   ) :: source_dn, source_up  # Diffuse radiation emitted by the layer
real(FT), dimension(ncol       ), intent(in   ) :: source_sfc # Surface source function [W/m2]
real(FT), dimension(ncol,nlay+1), intent(  out) :: radn_up    # Radiances [W/m2-str]
real(FT), dimension(ncol,nlay+1), intent(inout) :: radn_dn    # Top level must contain incident flux boundary condition
"""
function lw_transport_noscat!(ncol::I, nlay::I,
                              mo::MeshOrientation{I},
                              trans::Array{FT,2},
                              sfc_albedo::Array{FT,1},
                              source_dn::Array{FT,2},
                              source_up::Array{FT,2},
                              source_sfc::Array{FT,1},
                              radn_up::AbstractArray{FT,2},
                              radn_dn::AbstractArray{FT,2}) where {I<:Int, FT<:AbstractFloat}
  n̂ = nhat(mo)
  b = binary(mo)
  n_sfc = ilev_bot(mo)

  # Downward propagation
  @inbounds for ilev in lev_range(mo)
    radn_dn[:,ilev] .= trans[:,ilev-b].*radn_dn[:,ilev-n̂] .+ source_dn[:,ilev-b]
  end

  # Surface reflection and emission
  radn_up[:,n_sfc] .= radn_dn[:,n_sfc].*sfc_albedo .+ source_sfc

  # Upward propagation
  @inbounds for ilev in lev_range_reversed(mo)
    radn_up[:,ilev] .= trans[:,ilev-1+b].*radn_up[:,ilev+n̂] .+ source_up[:,ilev-1+b]
  end
  return nothing
end

"""
    lw_two_stream(τ::FT, ssa::FT, g::FT) where {FT<:AbstractFloat}

Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
   with optical depth `τ`, single scattering albedo `ssa`, and asymmetery parameter `g`.

Equations are developed in Meador and Weaver, 1980,
   doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
"""
function lw_two_stream(τ::FT, ssa::FT, g::FT) where {FT<:AbstractFloat}
  # Coefficients differ from SW implementation because the phase function is more isotropic
  #   Here we follow Fu et al. 1997, doi:10.1175/1520-0469(1997)054<2799:MSPITI>2.0.CO;2
  γ_1 = LW_diff_sec(FT) * (FT(1) - FT(0.5) * ssa * (FT(1) + g)) # Fu et al. Eq 2.9
  γ_2 = LW_diff_sec(FT) *          FT(0.5) * ssa * (FT(1) - g)  # Fu et al. Eq 2.10

  # Written to encourage vectorization of exponential, square root
  # Eq 18;  k = SQRT(γ_1**2 - γ_2**2), limited below to avoid div by 0.
  #   k = 0 for isotropic, conservative scattering; this lower limit on k
  #   gives relative error with respect to conservative solution
  #   of < 0.1% in Rdif down to τ = 10^-9
  temp1 = γ_1 - γ_2 * γ_1 + γ_2

  temp2 = max(temp1, FT(1.e-12))
  k = sqrt(temp2)
  exp_minuskτ = exp(-τ*k)

  # Diffuse reflection and transmission
  exp_minus2kτ = exp_minuskτ * exp_minuskτ

  # Refactored to avoid rounding errors when k, γ_1 are of very different magnitudes
  RT_term = FT(1) / (k * (FT(1) + exp_minus2kτ) +
                   γ_1 * (FT(1) - exp_minus2kτ) )

  Rdif = RT_term * γ_2 * (FT(1) - exp_minus2kτ) # Equation 25
  Tdif = RT_term * FT(2) * k * exp_minuskτ # Equation 26

  return γ_1,γ_2,Rdif,Tdif
end


"""
    lw_combine_sources!(ncol::I, nlay::I,
                        lev_src_inc::Array{FT,2},
                        lev_src_dec::Array{FT,2},
                        lev_source::Array{FT,2}) where {I<:Int,FT<:AbstractFloat}

Source function combination
RRTMGP provides two source functions at each level
  using the spectral mapping from each of the adjacent layers.
  Need to combine these for use in two-stream calculation.

 - `lev_src_inc`, `lev_src_dec`
 - `lev_source`
"""
function lw_combine_sources!(ncol::I, nlay::I,
                             lev_src_inc::Array{FT,2},
                             lev_src_dec::Array{FT,2},
                             lev_source::Array{FT,2}) where {I<:Int,FT<:AbstractFloat}
  ilay = 1
  @inbounds for icol in 1:ncol
    lev_source[icol, ilay] =        lev_src_dec[icol, ilay]
  end
  @inbounds for ilay in 2:nlay
    @inbounds for icol in 1:ncol
      lev_source[icol, ilay] = sqrt(lev_src_dec[icol, ilay] * lev_src_inc[icol, ilay-1])
    end
  end
  ilay = nlay+1
  @inbounds for icol in 1:ncol
    lev_source[icol, ilay] = lev_src_inc[icol, ilay-1]
  end
  return nothing
end

"""
    lw_source_2str!(ncol::I, nlay::I,
                    mo::MeshOrientation{I},
                    sfc_emis::Array{FT,1},
                    sfc_src::Array{FT,1},
                    lay_source::Array{FT,2},
                    lev_source::Array{FT,2},
                    γ_1::Array{FT,2},
                    γ_2::Array{FT,2},
                    rdif::Array{FT,2},
                    tdif::Array{FT,2},
                    τ::Array{FT,2},
                    source_dn::Array{FT,2},
                    source_up::Array{FT,2},
                    source_sfc::Array{FT,1}) where {I<:Int,FT<:AbstractFloat}

Compute LW source function for upward and downward emission at levels using linear-in-τ assumption
  This version straight from ECRAD
  Source is provided as W/m2-str; factor of π converts to flux units

 - `mo` mesh orientation, see [`MeshOrientation`](@ref)

real(FT), dimension(ncol      ), intent(in) :: sfc_emis, sfc_src
real(FT), dimension(ncol, nlay), intent(in) :: lay_source,     # Planck source at layer center
                                                τ,             # Optical depth (τ)
                                                γ_1, γ_2,      # Coupling coefficients
                                                rdif, tdif     # Layer reflectance and transmittance
real(FT), dimension(ncol, nlay+1), target, intent(in)  :: lev_source # Planck source at layer edges
real(FT), dimension(ncol, nlay), intent(out) :: source_dn, source_up
real(FT), dimension(ncol      ), intent(out) :: source_sfc   # Source function for upward radiation at surface

real(FT)            :: Z, Zup_top, Zup_bottom, Zdn_top, Zdn_bottom
real(FT), dimension(:), pointer :: lev_source_bot, lev_source_top
"""
function lw_source_2str!(ncol::I, nlay::I,
                         mo::MeshOrientation{I},
                         sfc_emis::Array{FT,1},
                         sfc_src::Array{FT,1},
                         lay_source::Array{FT,2},
                         lev_source::Array{FT,2},
                         γ_1::Array{FT,2},
                         γ_2::Array{FT,2},
                         rdif::Array{FT,2},
                         tdif::Array{FT,2},
                         τ::Array{FT,2},
                         source_dn::Array{FT,2},
                         source_up::Array{FT,2},
                         source_sfc::Array{FT,1}) where {I<:Int,FT<:AbstractFloat}
  b = binary(mo)
  @inbounds for ilay in 1:nlay
    lev_source_top = lev_source[:,ilay+1-b]
    lev_source_bot = lev_source[:,ilay+b]
    @inbounds for icol in 1:ncol
      if τ[icol,ilay] > FT(1.0e-8)
        #
        # Toon et al. (JGR 1989) Eqs 26-27
        #
        Z = (lev_source_bot[icol]-lev_source_top[icol]) / (τ[icol,ilay]*(γ_1[icol,ilay]+γ_2[icol,ilay]))
        Zup_top        =  Z + lev_source_top[icol]
        Zup_bottom     =  Z + lev_source_bot[icol]
        Zdn_top        = -Z + lev_source_top[icol]
        Zdn_bottom     = -Z + lev_source_bot[icol]
        source_up[icol,ilay] = π * (Zup_top    - rdif[icol,ilay] * Zdn_top    - tdif[icol,ilay] * Zup_bottom)
        source_dn[icol,ilay] = π * (Zdn_bottom - rdif[icol,ilay] * Zup_bottom - tdif[icol,ilay] * Zdn_top)
      else
        source_up[icol,ilay] = FT(0)
        source_dn[icol,ilay] = FT(0)
      end
    end
  end
  @inbounds for icol in 1:ncol
    source_sfc[icol] = π * sfc_emis[icol] * sfc_src[icol]
  end
  return nothing
end

#####
#####   Lower-level shortwave kernels
#####

"""
    sw_two_stream(μ_0::FT, τ::FT, ssa::FT, g::FT) where {FT<:AbstractFloat}

Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
   with optical depth `τ`, single scattering albedo `ssa`, and asymmetry parameter `g`.

Equations are developed in Meador and Weaver, 1980,
   doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2

Compute

 - `Rdif`, `Tdif`, `Rdir`, `Tdir`, `Tnoscat`

given

 - `μ_0`
 - `τ`, `ssa`, `g`

# Local Variables used in Meador and Weaver
real(FT) :: γ_1(ncol), γ_2(ncol), γ_3(ncol), γ_4(ncol)
real(FT) :: α_1(ncol), α_2(ncol), k(ncol)
"""
function sw_two_stream(μ_0::FT, τ::FT, ssa::FT, g::FT) where {FT<:AbstractFloat}

  μ_0_inv = 1 / μ_0
  # Zdunkowski Practical Improved Flux Method "PIFM"
  #  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
  γ_1 = (8 - ssa * (5 + 3 * g)) * FT(.25)
  γ_2 =  3 *(ssa * (1 -     g)) * FT(.25)
  γ_3 = (2 - 3 * μ_0 *      g ) * FT(.25)
  γ_4 =  1 - γ_3
  α_1 = γ_1 * γ_4 + γ_2 * γ_3           # Eq. 16
  α_2 = γ_1 * γ_3 + γ_2 * γ_4           # Eq. 17

  # Written to encourage vectorization of exponential, square root
  # Eq 18;  k = SQRT(γ_1**2 - γ_2**2), limited below to avoid div by 0.
  #   k = 0 for isotropic, conservative scattering; this lower limit on k
  #   gives relative error with respect to conservative solution
  #   of < 0.1% in Rdif down to τ = 10^-9
  temp = (γ_1 - γ_2) * (γ_1 + γ_2)
  k = sqrt(max(temp, FT(1.e-12)))
  exp_minuskτ = exp(-τ * k)

  # Diffuse reflection and transmission
  exp_minus2kτ = exp_minuskτ * exp_minuskτ

  # Refactored to avoid rounding errors when k, γ_1 are of very different magnitudes
  RT_term = 1 / (   k * (1 + exp_minus2kτ)  +
                  γ_1 * (1 - exp_minus2kτ) )

  Rdif = RT_term * γ_2 * (1 - exp_minus2kτ) # Equation 25
  Tdif = RT_term * 2 * k * exp_minuskτ      # Equation 26

  # Transmittance of direct, non-scattered beam. Also used below
  Tnoscat = exp(-τ * μ_0_inv)

  # Direct reflect and transmission
  k_μ   = k * μ_0
  k_γ_3 = k * γ_3
  k_γ_4 = k * γ_4

  # Equation 14, multiplying top and bottom by exp(-k*τ)
  #   and rearranging to avoid div by 0.
  cond = abs(1 - k_μ*k_μ) >= eps(FT)
  denom = cond ? 1 - k_μ*k_μ : eps(FT)
  RT_term =  ssa * RT_term/denom

  Rdir = RT_term  *
      ((1 - k_μ) * (α_2 + k_γ_3) -
       (1 + k_μ) * (α_2 - k_γ_3) * exp_minus2kτ -
        2 * (k_γ_3 - α_2 * k_μ)  * exp_minuskτ * Tnoscat)

  # Equation 15, multiplying top and bottom by exp(-k*τ),
  #   multiplying through by exp(-τ/μ_0) to
  #   prefer underflow to overflow
  # Omitting direct transmittance
  Tdir = -RT_term *
              ((1 + k_μ) * (α_1 + k_γ_4)                * Tnoscat -
               (1 - k_μ) * (α_1 - k_γ_4) * exp_minus2kτ * Tnoscat -
                2 * (k_γ_4 + α_1 * k_μ)  * exp_minuskτ)
  return Rdif,Tdif,Rdir,Tdir,Tnoscat
end

#####
##### Direct beam source for diffuse radiation in layers and at surface;
#####   report direct beam as a byproduct
#####

"""
    sw_source_2str!(ncol::I, nlay::I,
                    mo::MeshOrientation{I},
                    Rdir::Array{FT,2},
                    Tdir::Array{FT,2},
                    Tnoscat::Array{FT,2},
                    sfc_albedo::Array{FT,1},
                    source_up::Array{FT,2},
                    source_dn::Array{FT,2},
                    source_sfc::Array{FT,1},
                    flux_dn_dir::AbstractArray{FT,2}) where {I<:Int,FT<:AbstractFloat}

 - `mo` mesh orientation, see [`MeshOrientation`](@ref)

real(FT), dimension(ncol, nlay  ), intent(in   ) :: Rdir, Tdir, Tnoscat # Layer reflectance, transmittance for diffuse radiation
real(FT), dimension(ncol        ), intent(in   ) :: sfc_albedo          # surface albedo for direct radiation
real(FT), dimension(ncol, nlay  ), intent(  out) :: source_dn, source_up
real(FT), dimension(ncol        ), intent(  out) :: source_sfc          # Source function for upward radation at surface
real(FT), dimension(ncol, nlay+1), intent(inout) :: flux_dn_dir # Direct beam flux
                                                                # intent(inout) because top layer includes incident flux
"""
function sw_source_2str!(ncol::I, nlay::I,
                         mo::MeshOrientation{I},
                         Rdir::Array{FT,2},
                         Tdir::Array{FT,2},
                         Tnoscat::Array{FT,2},
                         sfc_albedo::Array{FT,1},
                         source_up::Array{FT,2},
                         source_dn::Array{FT,2},
                         source_sfc::Array{FT,1},
                         flux_dn_dir::AbstractArray{FT,2}) where {I<:Int,FT<:AbstractFloat}

  b = binary(mo)
  i_lev_sfc = ilev_bot(mo)

  @inbounds for ilev in lay_range(mo)
    source_up[:,ilev]     .=    Rdir[:,ilev] .* flux_dn_dir[:,ilev+1-b]
    source_dn[:,ilev]     .=    Tdir[:,ilev] .* flux_dn_dir[:,ilev+1-b]
    flux_dn_dir[:,ilev+b] .= Tnoscat[:,ilev] .* flux_dn_dir[:,ilev+1-b]
  end
  source_sfc .= flux_dn_dir[:,i_lev_sfc] .* sfc_albedo
  return nothing
end

#####
##### Transport of diffuse radiation through a vertically layered atmosphere.
#####   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
#####   This routine is shared by longwave and shortwave
#####

"""
    adding!(ncol::I, nlay::I,
            mo::MeshOrientation{I},
            albedo_sfc::Array{FT,1},
            rdif::Array{FT,2},
            tdif::Array{FT,2},
            src_dn::Array{FT,2},
            src_up::Array{FT,2},
            src_sfc::Array{FT,1},
            flux_up::AbstractArray{FT,2},
            flux_dn::AbstractArray{FT,2}) where {I<:Int,FT<:AbstractFloat}

 - `mo` mesh orientation, see [`MeshOrientation`](@ref)

real(FT), dimension(ncol       ), intent(in   ) :: albedo_sfc
real(FT), dimension(ncol,nlay  ), intent(in   ) :: rdif, tdif
real(FT), dimension(ncol,nlay  ), intent(in   ) :: src_dn, src_up
real(FT), dimension(ncol       ), intent(in   ) :: src_sfc
real(FT), dimension(ncol,nlay+1), intent(  out) :: flux_up
# intent(inout) because top layer includes incident flux
real(FT), dimension(ncol,nlay+1), intent(inout) :: flux_dn
# ------------------
real(FT), dimension(ncol,nlay+1)  :: albedo,   # reflectivity to diffuse radiation below this level
                                                # α in SH08
                                     src        # source of diffuse upwelling radiation from emission or
                                                # scattering of direct beam
                                                # G in SH08
real(FT), dimension(ncol,nlay  )  :: denom      # β in SH08
# ------------------
"""
function adding!(ncol::I, nlay::I,
                 mo::MeshOrientation{I},
                 albedo_sfc::Array{FT,1},
                 rdif::Array{FT,2},
                 tdif::Array{FT,2},
                 src_dn::Array{FT,2},
                 src_up::Array{FT,2},
                 src_sfc::Array{FT,1},
                 flux_up::AbstractArray{FT,2},
                 flux_dn::AbstractArray{FT,2}) where {I<:Int,FT<:AbstractFloat}
  albedo = Array{FT}(undef, ncol,nlay+1)
  src    = Array{FT}(undef, ncol,nlay+1)
  denom  = Array{FT}(undef, ncol,nlay)

  b = binary(mo)
  n̂ = nhat(mo)
  i_lev_sfc = ilev_bot(mo)
  i_lev_top = ilev_top(mo)

  # Indexing into arrays for upward and downward propagation depends on the vertical
  #   orientation of the arrays (whether the domain top is at the first or last index)
  # We write the loops out explicitly so compilers will have no trouble optimizing them.

  ilev = i_lev_sfc
  # Albedo of lowest level is the surface albedo...
  albedo[:,ilev] .= albedo_sfc
  # ... and source of diffuse radiation is surface emission
  src[:,ilev] .= src_sfc

  # From bottom to top of atmosphere --
  #   compute albedo and source of upward radiation
  @inbounds for ilev in lay_range_reversed(mo)
    denom[:, ilev] .= FT(1) ./ (FT(1) .- rdif[:,ilev].*albedo[:,ilev+b]) # Eq 10
    albedo[:,ilev+1-b] .= rdif[:,ilev] .+
                          tdif[:,ilev] .*
                          tdif[:,ilev] .*
                          albedo[:,ilev+b] .*
                          denom[:,ilev] # Equation 9

    # Equation 11 -- source is emitted upward radiation at top of layer plus
    #   radiation emitted at bottom of layer,
    #   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
    src[:,ilev+1-b] .= src_up[:, ilev] .+
                     tdif[:,ilev] .*
                     denom[:,ilev] .*
                     (src[:,ilev+b] .+
                      albedo[:,ilev+b] .*
                      src_dn[:,ilev])
  end

  # Eq 12, at the top of the domain upwelling diffuse is due to ...
  ilev = i_lev_top
  flux_up[:,ilev] .= flux_dn[:,ilev] .* albedo[:,ilev] .+  # ... reflection of incident diffuse and
                    src[:,ilev]                            # emission from below/scattering by the direct beam below
  # From the top of the atmosphere downward -- compute fluxes
  @inbounds for ilev in lev_range(mo)
    flux_dn[:,ilev] .= (tdif[:,ilev-b].*flux_dn[:,ilev-n̂] +  # Equation 13
                       rdif[:,ilev-b].*src[:,ilev] +
                       src_dn[:,ilev-b]) .* denom[:,ilev-b]
    flux_up[:,ilev] .= flux_dn[:,ilev] .* albedo[:,ilev] .+  # Equation 12
                      src[:,ilev]
  end

  return nothing
end

#####
##### Upper boundary condition
#####

"""
    apply_BC!(flux_dn::Array{FT},
              ilay::I) where {I<:Integer,FT<:AbstractFloat}

 - `ilay` apply BC at the i-th layer
 - `flux_dn` Flux to be used as input to solvers below
"""
function apply_BC!(flux_dn::Array{FT,3},
                   ilay::I,
                   ::Nothing) where {I<:Integer,FT<:AbstractFloat}
  flux_dn[:,ilay, :] .= FT(0)
  return nothing
end

"""
    apply_BC!(flux_dn::Array{FT,3},
              ilay::I,
              inc_flux::Array{FT,2}) where {I<:Integer,B<:Bool,FT<:AbstractFloat}

 - `flux_dn` Flux to be used as input to solvers below
 - `ilay` apply BC at the i-th layer
 - `inc_flux` Flux at top of domain
"""
function apply_BC!(flux_dn::Array{FT,3},
                   ilay::I,
                   inc_flux::Array{FT,2}) where {I<:Integer,FT<:AbstractFloat}
  fill!(flux_dn, 0)
  flux_dn[:, ilay, :] .= inc_flux
  return nothing
end

"""
    apply_BC!(flux_dn::Array{FT},
              ilay::I,
              inc_flux::Array{FT,2},
              factor::Array{FT,1}) where {I<:Integer,B<:Bool,FT<:AbstractFloat}

 - `flux_dn` Flux to be used as input to solvers below
 - `ilay` apply BC at the i-th layer
 - `inc_flux` Flux at top of domain
 - `factor` Factor to multiply incoming flux
"""
function apply_BC!(flux_dn::Array{FT,3},
                   ilay::I,
                   inc_flux::Array{FT,2},
                   factor::Array{FT,1}) where {I<:Integer,FT<:AbstractFloat}
  fill!(flux_dn, 0)
  flux_dn[:, ilay, :]  .= inc_flux .* spread(factor, 2, size(inc_flux,2))
  return nothing
end
