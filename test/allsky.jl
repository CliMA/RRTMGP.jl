using Profile
using Test
using RRTMGP
using TimerOutputs
const to = TimerOutput()
using NCDatasets
using RRTMGP.OpticalProps
using RRTMGP.FortranIntrinsics
using RRTMGP.ArrayUtilities
using RRTMGP.GasOptics
using RRTMGP.GasConcentrations
using RRTMGP.RTESolver
using RRTMGP.Fluxes
using RRTMGP.SourceFunctions
using RRTMGP.CloudOptics

include(joinpath("..","ReadInputData","ReadInputs.jl"))
include(joinpath("..","ReadInputData","LoadCoefficients.jl"))
include(joinpath("..","ReadInputData","LoadCloudCoefficients.jl"))
include("CloudSampling.jl")

function vmr_2d_to_1d!(gas_concs::GasConcs{FT},
                       gas_concs_garand::GasConcs{FT},
                       name::String,
                       sz1::I,
                       sz2::I) where {FT<:AbstractFloat,I<:Int}
  tmp = Array{FT}(undef, sz1, sz2)
  get_vmr!(tmp, gas_concs_garand, name)
  tmp_col = tmp[1, :]
  set_vmr!(gas_concs, name, tmp_col)
end

function all_sky(ds; use_luts=false, λ_string="", compile_first=false)
  @assert λ_string == "sw" || λ_string == "lw"
  k_dist_sym = Symbol(:k_dist,λ_string)
  cloud_optics_sym = Symbol(:cloud_optics,λ_string)

  gas_names = lowercase.(strip.(["h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"]))
  ngas = length(gas_names)
  nloops = 1
  ncol = 128
  write_fluxes = true

  #
  # Read temperature, pressure, gas concentrations.
  #   Arrays are allocated as they are read
  #
  FT = Float64
  I = Int64
  p_lay, t_lay, p_lev, t_lev, gas_concs_garand, col_dry = @timeit to "read_atmos" read_atmos(ds[:input], FT, I, gas_names)

  col_dry = nothing
  nlay = size(p_lay, 2)
  # For clouds we'll use the first column, repeated over and over
  gsc = GasConcSize(ncol, nlay, (ncol, nlay), ngas)
  gas_concs = GasConcs(FT, I, gas_names, ncol, nlay, gsc)
  for igas = 1:ngas
    vmr_2d_to_1d!(gas_concs, gas_concs_garand, gas_names[igas], size(p_lay, 1), nlay)
  end

  #  If we trusted in Fortran allocate-on-assign we could skip the temp_array here
  p_lay = deepcopy(spread_new(p_lay[1,:], 1, ncol))
  t_lay = deepcopy(spread_new(t_lay[1,:], 1, ncol))
  p_lev = deepcopy(spread_new(p_lev[1,:], 1, ncol))
  t_lev = deepcopy(spread_new(t_lev[1,:], 1, ncol))

  # This puts pressure and temperature arrays on the GPU
  # load data into classes
  k_dist = load_and_init(ds[:k_dist], FT, gas_concs.gas_names)
  is_sw = source_is_external(k_dist)
  is_lw = !is_sw
  #
  # Should also try with Pade calculations
  if use_luts
    cloud_optics_ = load_cld_lutcoeff(FT,ds[:cloud_optics], 2)
  else
    cloud_optics_ = load_cld_padecoeff(FT,ds[:cloud_optics], 2)
  end

  #
  # Problem sizes
  #
  nbnd = get_nband(k_dist.optical_props)
  ngpt = get_ngpt(k_dist.optical_props)
  top_at_1 = p_lay[1, 1] < p_lay[1, nlay]

  println("--------- Problem size:")
  @show ncol,nlay,nbnd,ngpt

  ps = ProblemSize(ncol, nlay, ngpt)

  # Clouds optical props are defined by band
  clouds_base = OpticalPropsBase("Clouds", get_band_lims_wavenumber(k_dist.optical_props))

  # LW calculations neglect scattering; SW calculations use the 2-stream approximation
  #   Here we choose the right variant of optical_props.
  if is_sw
    clouds = TwoStream(clouds_base, ps)
    atmos = TwoStream(k_dist.optical_props,ps)
  else
    clouds = OneScalar(clouds_base, ps)
    atmos = OneScalar(k_dist.optical_props, ps)
  end

  #
  # Allocate arrays for the optical properties themselves.
  #

  #  Boundary conditions depending on whether the k-distribution being supplied
  #   is LW or SW
  if is_sw
    # toa_flux is thread private
    toa_flux = zeros(FT, ncol, ngpt)
    #
    sfc_alb_dir = zeros(FT, nbnd, ncol)
    sfc_alb_dif = zeros(FT, nbnd, ncol)
    μ_0 = zeros(FT, ncol)
    # Ocean-ish values for no particular reason
    sfc_alb_dir .= FT(0.06)
    sfc_alb_dif .= FT(0.06)
    μ_0 .= FT(.86)
  else
    # lw_sorces is thread private
    lw_sources = SourceFuncLW(ncol, nlay, k_dist.optical_props)

    t_sfc = zeros(FT, ncol)
    emis_sfc = zeros(FT, nbnd, ncol)
    # Surface temperature
    t_sfc .= t_lev[1, fmerge(nlay+1, 1, top_at_1)]
    emis_sfc .= FT(0.98)
  end

  #
  # Fluxes
  #
  flux_up = zeros(FT, ncol,nlay+1)
  flux_dn = zeros(FT, ncol,nlay+1)

  if is_sw
    flux_dir = zeros(FT, ncol,nlay+1)
  end
  #
  # Clouds
  #
  lwp = zeros(FT, ncol,nlay)
  iwp = zeros(FT, ncol,nlay)
  rel = zeros(FT, ncol,nlay)
  rei = zeros(FT, ncol,nlay)
  cloud_mask = Array{Bool}(undef, ncol,nlay)

  # Restrict clouds to troposphere (< 100 hPa = 100*100 Pa)
  #   and not very close to the ground
  rel_val = FT(0.5) * (get_min_radius(cloud_optics_.liq) + get_max_radius(cloud_optics_.liq))
  rei_val = FT(0.5) * (get_min_radius(cloud_optics_.ice) + get_max_radius(cloud_optics_.ice))
  for ilay=1:nlay
    for icol=1:ncol
      cloud_mask[icol,ilay] = p_lay[icol,ilay] < FT(100) * FT(100) && p_lay[icol,ilay] > FT(900)
      #
      # Ice and liquid will overlap in a few layers
      #
      lwp[icol,ilay] = fmerge(FT(10),  FT(0), cloud_mask[icol,ilay] && t_lay[icol,ilay] > FT(263))
      iwp[icol,ilay] = fmerge(FT(10),  FT(0), cloud_mask[icol,ilay] && t_lay[icol,ilay] < FT(273))
      rel[icol,ilay] = fmerge(rel_val, FT(0), lwp[icol,ilay] > FT(0))
      rei[icol,ilay] = fmerge(rei_val, FT(0), iwp[icol,ilay] > FT(0))
    end
  end

  fluxes = FluxesBroadBand(FT, size(flux_up), is_sw)
  #
  # Multiple iterations for big problem sizes, and to help identify data movement
  #   For CPUs we can introduce OpenMP threading over loop iterations
  #

  for iloop = 1:(compile_first ? 1 : nloops)
    cloud_optics!(cloud_optics_, lwp, iwp, rel, rei, clouds)
    #
    # Solvers
    #
    fluxes.flux_up .= FT(0)
    fluxes.flux_dn .= FT(0)
    if is_lw
      gas_optics_int!(k_dist, p_lay, p_lev,
                  t_lay, t_sfc,
                  gas_concs,
                  atmos,
                  lw_sources;
                  tlev = t_lev)


      increment!(clouds, atmos)
      rte_lw!(atmos, top_at_1,
              lw_sources,
              emis_sfc,
              fluxes)
    else
      fluxes.flux_dn_dir .= flux_dir

      gas_optics_ext!(k_dist, p_lay, p_lev,
                  t_lay,
                  gas_concs,
                  atmos,
                  toa_flux)
      delta_scale!(clouds)
      increment!(clouds, atmos)
      rte_sw!(atmos, top_at_1,
              μ_0,   toa_flux,
              sfc_alb_dir, sfc_alb_dif,
              fluxes)
    end
    flux_up .= fluxes.flux_up
    flux_dn .= fluxes.flux_dn
  end


  # Compare with reference:
  ref_flux_up = Array{FT}(ds[:ref][λ_string*"_flux_up"][:])
  ref_flux_dn = Array{FT}(ds[:ref][λ_string*"_flux_dn"][:])

  diff_up = maximum( abs.( flux_up .- ref_flux_up ) )
  diff_dn = maximum( abs.( flux_dn .- ref_flux_dn ) )

  diff_up_ulps = maximum( abs.( flux_up .- ref_flux_up ) ./ eps.(ref_flux_up) )
  diff_dn_ulps = maximum( abs.( flux_dn .- ref_flux_dn ) ./ eps.(ref_flux_dn) )

  # @show sqrt(1/eps(FT))
  # @show diff_up, diff_up_ulps, maximum(abs.(ref_flux_up))
  # @show diff_dn, diff_dn_ulps, maximum(abs.(ref_flux_dn))
  if !compile_first
    if use_luts
      @test diff_up_ulps < sqrt(1/(10eps(FT)))
      @test diff_dn_ulps < sqrt(1/(10eps(FT)))
    else
      # Need better test
      # @show diff_up_ulps < sqrt(1/(eps(FT)))
      # @show diff_dn_ulps < sqrt(1/(eps(FT)))
    end
  end
  return nothing

end

