using Profile
using Test
using RRTMGP
using TimerOutputs
const to = TimerOutput()
using NCDatasets
using RRTMGP.OpticalProps
using RRTMGP.Gases
using RRTMGP.FortranIntrinsics
using RRTMGP.Utilities
using RRTMGP.GasOptics
using RRTMGP.GasConcentrations
using RRTMGP.RadiativeBoundaryConditions
using RRTMGP.RTESolver
using RRTMGP.Fluxes
using RRTMGP.SourceFunctions
using RRTMGP.AtmosphericStates
using RRTMGP.CloudOptics
@static if haspkg("Plots")
  using Plots
  const export_plots = true
else
  const export_plots = false
end

include(joinpath("PostProcessing.jl"))

include(joinpath("..","ReadInputData","ReadInputs.jl"))
include(joinpath("..","ReadInputData","LoadCoefficients.jl"))
include(joinpath("..","ReadInputData","LoadCloudCoefficients.jl"))
include("CloudSampling.jl")

function vmr_2d_to_1d!(gas_conc::GasConcs{FT},
                       gas_conc_garand::GasConcs{FT},
                       gas::AbstractGas,
                       gases_prescribed::Vector{AbstractGas}) where {FT<:AbstractFloat, I<:Int}
  i_gas = loc_in_array(gas, gases_prescribed)
  tmp = gas_conc_garand.concs[i_gas,:,:]
  tmp_col = tmp[1, :]
  set_vmr!(gas_conc, gas, tmp_col)
end

function all_sky(ds; use_luts=false, λ_string="", compile_first=false)
  @assert λ_string == "sw" || λ_string == "lw"
  k_dist_sym = Symbol(:k_dist,λ_string)
  cloud_optics_sym = Symbol(:cloud_optics,λ_string)

  gases_prescribed = [h2o(), co2(), o3(), n2o(), co(), ch4(), o2(), n2()]
  ngas = length(gases_prescribed)
  nloops = 1
  ncol = 128
  write_fluxes = true

  #
  # Read temperature, pressure, gas concentrations.
  #   Arrays are allocated as they are read
  #
  FT = Float64
  I = Int64
  p_lay, t_lay, p_lev, t_lev, gas_conc_garand, col_dry = @timeit to "read_atmos" read_atmos(ds[:input], FT, I, gases_prescribed)

  col_dry = nothing
  nlay = size(p_lay, 2)
  # For clouds we'll use the first column, repeated over and over
  gsc = GasConcSize(ncol, nlay, (ncol, nlay), ngas)
  gas_conc = GasConcs(FT, I, gases_prescribed, ncol, nlay, gsc)
  for gas in gases_prescribed
    vmr_2d_to_1d!(gas_conc, gas_conc_garand, gas, gases_prescribed)
  end

  #  If we trusted in Fortran allocate-on-assign we could skip the temp_array here
  p_lay = deepcopy(spread_new(p_lay[1,:], 1, ncol))
  t_lay = deepcopy(spread_new(t_lay[1,:], 1, ncol))
  p_lev = deepcopy(spread_new(p_lev[1,:], 1, ncol))
  t_lev = deepcopy(spread_new(t_lev[1,:], 1, ncol))

  # This puts pressure and temperature arrays on the GPU
  # load data into classes
  k_dist = load_and_init(ds[:k_dist], FT, gas_conc.gas_names)
  is_sw = source_is_external(k_dist)
  is_lw = !is_sw

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

  println("--------- Problem size:")
  @show ncol,nlay,nbnd,ngpt

  # Clouds optical props are defined by band
  clouds_base = OpticalPropsBase("Clouds", get_band_lims_wavenumber(k_dist.optical_props))

  # LW calculations neglect scattering; SW calculations use the 2-stream approximation
  #   Here we choose the right variant of optical_props.
  if is_sw
    clouds = TwoStream(clouds_base, ncol, nlay, ngpt)
    atmos = TwoStream(k_dist.optical_props,ncol, nlay, ngpt)
  else
    clouds = OneScalar(clouds_base, ncol, nlay, ngpt)
    atmos = OneScalar(k_dist.optical_props, ncol, nlay, ngpt)
  end

  top_at_1 = p_lay[1, 1] < p_lay[1, nlay]

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

    # Surface emissivity
    sfc_emis = zeros(FT, nbnd, ncol)
    sfc_emis .= FT(0.98)
  end

  as = AtmosphericState(gas_conc, p_lay, p_lev, t_lay, t_lev, k_dist.ref, nothing, nothing)
  gas_conc, p_lay, p_lev, t_lay, t_lev = ntuple(i->nothing,5)

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
  clouds_ice = CloudOpticalProps(FT, ncol,nlay)
  clouds_liq = CloudOpticalProps(FT, ncol,nlay)
  cloud_mask = Array{Bool}(undef, ncol,nlay)

  # Restrict clouds to troposphere (< 100 hPa = 100*100 Pa)
  #   and not very close to the ground
  rel_val = FT(0.5) * (get_min_radius(cloud_optics_.liq) + get_max_radius(cloud_optics_.liq))
  rei_val = FT(0.5) * (get_min_radius(cloud_optics_.ice) + get_max_radius(cloud_optics_.ice))
  for ilay=1:nlay
    for icol=1:ncol
      cloud_mask[icol,ilay] = as.p_lay[icol,ilay] < FT(100) * FT(100) && as.p_lay[icol,ilay] > FT(900)
      #
      # Ice and liquid will overlap in a few layers
      #
      clouds_liq.wp[icol,ilay] = fmerge(FT(10),  FT(0), cloud_mask[icol,ilay] && as.t_lay[icol,ilay] > FT(263))
      clouds_ice.wp[icol,ilay] = fmerge(FT(10),  FT(0), cloud_mask[icol,ilay] && as.t_lay[icol,ilay] < FT(273))
      clouds_liq.re[icol,ilay] = fmerge(rel_val, FT(0), clouds_liq.wp[icol,ilay] > FT(0))
      clouds_ice.re[icol,ilay] = fmerge(rei_val, FT(0), clouds_ice.wp[icol,ilay] > FT(0))
    end
  end

  fluxes = FluxesBroadBand(FT, size(flux_up), is_sw)
  #
  # Multiple iterations for big problem sizes, and to help identify data movement
  #   For CPUs we can introduce OpenMP threading over loop iterations
  #

  for iloop = 1:(compile_first ? 1 : nloops)
    # cloud_optics!(cloud_optics_, clouds_liq.wp, clouds_ice.wp, clouds_liq.re, clouds_ice.re, clouds)
    cloud_optics!(cloud_optics_, clouds_liq, clouds_ice, clouds)
    #
    # Solvers
    #
    fluxes.flux_up .= FT(0)
    fluxes.flux_dn .= FT(0)


    if is_lw
      gas_optics!(k_dist, as, atmos, lw_sources)

      increment!(clouds, atmos)
      bcs = LongwaveBCs(sfc_emis)
      rte_lw!(atmos,
              as.top_at_1,
              lw_sources,
              bcs,
              fluxes)
    else
      fluxes.flux_dn_dir .= flux_dir

      gas_optics!(k_dist, as, atmos)

      check_extent(toa_flux, (as.ncol, ngpt), "toa_flux")
      toa_flux .= repeat(k_dist.solar_src', as.ncol)

      delta_scale!(clouds)
      increment!(clouds, atmos)
      bcs = ShortwaveBCs(toa_flux, sfc_alb_dir, sfc_alb_dif)
      rte_sw!(atmos,
              as.top_at_1,
              μ_0,
              bcs,
              fluxes)
    end
    flux_up .= fluxes.flux_up
    flux_dn .= fluxes.flux_dn
  end

  if export_plots
    case = "AllSky_use_luts_"*string(use_luts)*"_$(λ_string)"
    heating_rate, z = compute_heating_rate(fluxes.flux_up, fluxes.flux_dn, as)
    plot(heating_rate, z, title="All sky $(λ_string) heating rates",
                          xlabel="heating rate",
                          ylabel="pressure")
    out_dir = "output"
    mkpath(out_dir)
    savefig(joinpath(out_dir,case*".png"))
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

