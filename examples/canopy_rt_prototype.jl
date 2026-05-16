#=
Canopy RT Prototype for RRTMGP.jl

Demonstrates coupled atmosphere-canopy shortwave radiative transfer by treating
canopy layers as additional layers below the atmosphere. The RRTMGP 2-stream
solver is layer-agnostic: it takes arrays of (τ, ssa, g) and boundary conditions.

Layer ordering (RRTMGP convention): index 1 = bottom (soil), index nlay = top (TOA)
  - Layers 1:ncanlay         → canopy sublayers (soil at bottom, canopy top at top)
  - Layers ncanlay+1:nlay_total → atmospheric layers
  - Soil albedo replaces surface albedo as lower boundary condition

Setup:
  cd examples/
  julia --project=.
  ] dev ..                         # adds local RRTMGP
  ] dev ~/code/gitHub/CanopyOptics.jl  # adds local CanopyOptics
  ] instantiate
  include("canopy_rt_prototype.jl")
=#

using NCDatasets
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

using RRTMGP
using RRTMGP: RRTMGPGridParams
using RRTMGP.Vmrs
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.RTE
using RRTMGP.RTESolver
using RRTMGP.ArtifactPaths
import RRTMGP.Parameters.RRTMGPParameters
import ClimaParams as CP

using CanopyOptics
using CanopyOptics: G, spherical_leaves
using Distributions: pdf
using Unitful
using Statistics: mean

# Access inner solver function (not exported but needed for per-gpoint control)
import RRTMGP.RTESolver: rte_sw_2stream!, sw_2stream_coeffs
import RRTMGP.RTESolver: rte_lw_2stream!

# Include test helpers for atmospheric state setup
include(joinpath(dirname(@__DIR__), "test", "reference_files.jl"))
include(joinpath(dirname(@__DIR__), "test", "read_clear_sky.jl"))

# ============================================================================
# Configuration
# ============================================================================
const FT = Float64
const ncol = 1          # single column prototype
const ncanlay = 5       # canopy sublayers
const LAI_total = FT(3) # total LAI [m²/m²]
const soil_alb_val = FT(0.2)  # broadband soil albedo
const sza_deg = FT(30)  # solar zenith angle [degrees]

# LW configuration
const ε_leaf = FT(0.95)          # leaf thermal emissivity (generic, editable)
const soil_emis_val = FT(0.96)   # soil thermal emissivity
const ω_lw = FT(1) - ε_leaf     # LW single scattering albedo
const g_lw = FT(0)               # isotropic scattering (R_lw = T_lw)
const G_hemi = FT(0.5)           # hemispheric-mean G (for diffuse-only LW)

# PROSPECT-PRO leaf parameters
const leaf = LeafProspectProProperties{FT}(
    N     = FT(1.4),
    Ccab  = FT(40),
    Ccar  = FT(10),
    Canth = FT(0.5),
    Cbrown = FT(0),
    Cw    = FT(0.012),
    Cm    = FT(0),
    Cprot = FT(0.001),
    Ccbc  = FT(0.009),
)

# ============================================================================
# Section 1: Load RRTMGP SW lookup tables
# ============================================================================
println("="^60)
println("Canopy RT Prototype")
println("="^60)
println("\n--- Loading RRTMGP SW lookup tables ---")

context = ClimaComms.context()
device = ClimaComms.device(context)
DA = ClimaComms.array_type(device)

overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
param_set = RRTMGPParameters(FT, overrides)

sw_file = get_lookup_filename(:gas, :sw)
lookup_sw, idx_gases = Dataset(sw_file, "r") do ds
    LookUpSW(ds, FT, DA)
end

nbnd = LookUpTables.get_n_bnd(lookup_sw)
n_gpt = length(lookup_sw.solar_src_scaled)
bnd_lims_wn = Array(lookup_sw.band_data.bnd_lims_wn)  # (2, nbnd) cm⁻¹
major_gpt2bnd = Array(lookup_sw.band_data.major_gpt2bnd)

println("  SW bands: $nbnd, g-points: $n_gpt")

println("\n--- Loading RRTMGP LW lookup tables ---")
lw_file = get_lookup_filename(:gas, :lw)
lookup_lw, _ = Dataset(lw_file, "r") do ds
    LookUpLW(ds, FT, DA)
end
nbnd_lw = LookUpTables.get_n_bnd(lookup_lw)
n_gpt_lw = length(lookup_lw.band_data.major_gpt2bnd)
major_gpt2bnd_lw = Array(lookup_lw.band_data.major_gpt2bnd)
println("  LW bands: $nbnd_lw, g-points: $n_gpt_lw")

# Extract Planck lookup tables (CPU arrays for canopy T-profile scaling)
t_planck_cpu = Array(lookup_lw.planck.t_planck)
tot_planck_cpu = Array(lookup_lw.planck.tot_planck)  # (n_t_plnk, nbnd_lw)

# ============================================================================
# Section 2: Set up atmospheric state
# ============================================================================
println("\n--- Setting up atmospheric state ---")

input_file = get_input_filename(:gas, :lw)
ds_in = Dataset(input_file, "r")
(as, _, sfc_alb_orig, cos_zenith_orig, toa_flux_orig, bot_at_1) =
    setup_clear_sky_as(context, ds_in, idx_gases, 1, lookup_sw, ncol, FT, Vmrs.VmrGM, param_set)
close(ds_in)

nlay_atm, _ = AtmosphericStates.get_dims(as)
nlev_atm = nlay_atm + 1

# Override zenith angle and TOA flux for our prototype
cos_zenith = DA{FT, 1}([cosd(sza_deg)])
toa_flux = DA{FT, 1}(Array(toa_flux_orig)[1:ncol])

println("  Atmospheric layers: $nlay_atm")
println("  SZA: $(sza_deg)°, μ₀: $(cosd(sza_deg))")

# Canopy temperature profile: linear from T_soil (surface) to T_air (atm bottom)
T_sfc_val = FT(Array(as.t_sfc)[1])      # soil skin temperature
T_air_val = FT(Array(as.t_lev)[1, 1])   # atmospheric bottom-level temperature
println("  T_soil: $(round(T_sfc_val, digits=2)) K, T_air: $(round(T_air_val, digits=2)) K")

# ============================================================================
# Section 3: Compute PROSPECT-PRO leaf optics and band-average
# ============================================================================
println("\n--- Computing leaf optics (PROSPECT-PRO) ---")

# Create leaf optical properties at 1nm resolution (400–2500 nm)
opti = createLeafOpticalStruct((400.0:2500.0) * u"nm")
T_leaf, R_leaf = prospect(leaf, opti)  # (T, R) at each wavelength
λ_nm = ustrip.(u"nm", opti.λ)  # wavelength grid in nm (strip units)

println("  PROSPECT wavelength range: $(λ_nm[1])–$(λ_nm[end]) nm")
println("  Mean R: $(round(mean(R_leaf), digits=4)), Mean T: $(round(mean(T_leaf), digits=4))")

# Band-average R, T over each RRTMGP SW band
# For each band, average over PROSPECT wavelengths within the band
# Using uniform weighting (solar weighting can improve accuracy in Phase 2)
leaf_R_band = zeros(FT, nbnd)
leaf_T_band = zeros(FT, nbnd)
band_in_prospect = falses(nbnd)

for ibnd in 1:nbnd
    wn_lo, wn_hi = bnd_lims_wn[1, ibnd], bnd_lims_wn[2, ibnd]
    # Convert wavenumber (cm⁻¹) to wavelength (nm): λ = 1e7/ν
    λ_hi_nm = 1e7 / wn_lo  # low wavenumber → long wavelength
    λ_lo_nm = 1e7 / wn_hi  # high wavenumber → short wavelength

    # Find PROSPECT wavelengths within this band
    idx = findall(λ_lo_nm .≤ λ_nm .≤ λ_hi_nm)

    if !isempty(idx)
        leaf_R_band[ibnd] = mean(R_leaf[idx])
        leaf_T_band[ibnd] = mean(T_leaf[idx])
        band_in_prospect[ibnd] = true
    end
    # Bands outside PROSPECT range: R=T=0 → canopy fully absorbing
    # (but τ=0 for those bands, so no effect; see below)
end

println("\n  Band-averaged leaf properties:")
for ibnd in 1:nbnd
    wn_lo, wn_hi = bnd_lims_wn[1, ibnd], bnd_lims_wn[2, ibnd]
    λ_lo, λ_hi = round(1e7/wn_hi, digits=0), round(1e7/wn_lo, digits=0)
    flag = band_in_prospect[ibnd] ? "" : " [outside PROSPECT]"
    ω = leaf_R_band[ibnd] + leaf_T_band[ibnd]
    println("    Band $ibnd ($(λ_lo)–$(λ_hi) nm): R=$(round(leaf_R_band[ibnd],digits=4)), " *
            "T=$(round(leaf_T_band[ibnd],digits=4)), ω=$(round(ω,digits=4))$flag")
end

# ============================================================================
# Section 4: Canopy → TwoStream mapping
# ============================================================================
"""Interpolate band-integrated Planck B(T) from RRTMGP lookup tables."""
function planck_interp(T::FT, t_planck, totplnk_band) where {FT}
    T < t_planck[1] && return FT(totplnk_band[1])
    T > t_planck[end] && return FT(totplnk_band[end])
    dT = t_planck[2] - t_planck[1]
    loc = min(floor(Int, (T - t_planck[1]) / dT) + 1, length(t_planck) - 1)
    frac = (T - t_planck[loc]) / dT
    return totplnk_band[loc] * (FT(1) - frac) + totplnk_band[loc + 1] * frac
end

"""Linear canopy temperature profile from T_soil (level 1) to T_air (level ncanlay+1)."""
function canopy_temp_profile(ilev::Int, ncanlay::Int, T_soil::FT, T_air::FT) where {FT}
    frac = FT(ilev - 1) / FT(ncanlay)
    return T_soil + frac * (T_air - T_soil)
end

"""
    canopy_2stream_props(R_b, T_b, G_μ₀, LAI_layer, FT)

Map canopy leaf properties to TwoStream parameters for one layer and band.

Returns (τ, ssa, g_eff) where:
  - τ = G(μ₀) × LAI_layer (optical depth)
  - ssa = R + T (single scattering albedo)
  - g_eff maps PIFM backward-scatter fraction to canopy β = R/(R+T)
"""
function canopy_2stream_props(R_b::FT, T_b::FT, G_μ₀::FT, LAI_layer::FT) where {FT}
    τ = G_μ₀ * LAI_layer
    ω = R_b + T_b
    if ω > eps(FT)
        β = R_b / ω
        # Derived from PIFM: γ₂/(γ₁+γ₂) = β
        # γ₁ = (8 - ω(5+3g))/4, γ₂ = 3ω(1-g)/4
        # Solving: g = 1 - β(8-2ω)/(3ω)
        g_eff = clamp(FT(1) - β * (FT(8) - FT(2) * ω) / (FT(3) * ω), FT(-1), FT(1))
    else
        g_eff = FT(0)
    end
    return (τ, ω, g_eff)
end

"""
    compute_β_dir(μ₀, R_b, T_b, ld; nθ=90, nφ=360)

Direct-beam upscatter fraction β₀(μ₀) for a canopy layer.

Integrates over the leaf angle distribution to find the fraction of single-scattered
radiation from a collimated beam at zenith angle acos(μ₀) that goes upward.
For each leaf orientation (θ_L, φ_L), computes the fraction of Lambertian-reflected
and transmitted light entering the upper hemisphere.

Returns β₀ ∈ [0,1]: fraction of scattered radiation going upward.
  - β₀ → R/(R+T) when μ₀ → 1 (overhead sun, horizontal leaves dominate)
  - β₀ → 0.5 when μ₀ → 0 (grazing sun, up/down symmetry)
"""
function compute_β_dir(μ₀::FT, R_b::FT, T_b::FT, ld; nθ=90, nφ=360) where {FT}
    ω = R_b + T_b
    ω < eps(FT) && return FT(0.5)

    sin_θ₀ = sqrt(max(FT(0), FT(1) - μ₀^2))
    dθ = FT(π / 2) / nθ
    dφ = FT(2π) / nφ

    num = FT(0)  # ∫ upscatter × |cos γ| × f_L dΩ
    den = FT(0)  # ∫ |cos γ| × f_L dΩ

    for iθ in 1:nθ
        θ_L = (iθ - FT(0.5)) * dθ
        sin_θ_L = sin(θ_L)
        cos_θ_L = cos(θ_L)
        # Leaf angle distribution weight (CanopyOptics convention)
        f_L = FT(pdf(ld.LD, 2θ_L / π) * ld.scaling)

        for iφ in 1:nφ
            φ_L = (iφ - FT(0.5)) * dφ
            # Projection of beam onto leaf normal: cos γ = n̂·(-Ω̂₀)
            cos_γ = cos_θ_L * μ₀ - sin_θ_L * cos(φ_L) * sin_θ₀
            abs_cos_γ = abs(cos_γ)

            # Lambertian scattering from tilted leaf:
            # Fraction of emitted light going into upper hemisphere = (1 + cos α)/2
            # where α is the angle of the emitting surface normal from vertical.
            # For reflected light: emitted from the hit surface
            #   cos γ > 0 (upper surface hit): normal up, α = θ_L → (1+cos θ_L)/2
            #   cos γ < 0 (lower surface hit): normal down, α = π-θ_L → (1-cos θ_L)/2
            # For transmitted light: emitted from the opposite surface
            s = cos_γ > 0 ? FT(1) : FT(-1)
            up_R = (FT(1) + s * cos_θ_L) / FT(2)
            up_T = (FT(1) - s * cos_θ_L) / FT(2)
            upscatter = R_b * up_R + T_b * up_T

            w = abs_cos_γ * f_L * dθ * dφ
            num += upscatter * w
            den += w
        end
    end

    # β₀ = total upscatter / (ω × total intercepted)
    return den > eps(FT) ? num / (ω * den) : FT(0.5)
end

"""
    canopy_sw_2stream_coeffs(τ, ssa, g, μ₀, β_dir)

Two-stream coefficients for canopy layers with physically-computed γ₃.

Like `sw_2stream_coeffs` but uses `β_dir` (directional backscatter from leaf angle
geometry) as γ₃ instead of the PIFM approximation γ₃ = (2 - 3μ₀g)/4.
γ₁ and γ₂ still use g_eff (correct for diffuse hemispherically-averaged scattering).
"""
function canopy_sw_2stream_coeffs(τ::FT, ssa::FT, g::FT, μ₀::FT, β_dir::FT) where {FT}
    k_min = sqrt(eps(FT))
    # PIFM γ₁, γ₂ from g_eff (hemispheric, correct for diffuse)
    γ1 = (FT(8) - ssa * (FT(5) + FT(3) * g)) * FT(0.25)
    γ2 = FT(3) * (ssa * (FT(1) - g)) * FT(0.25)
    # Direct beam: use β_dir from leaf angle geometry instead of PIFM
    γ3 = β_dir
    γ4 = FT(1) - γ3
    α1 = γ1 * γ4 + γ2 * γ3
    α2 = γ1 * γ3 + γ2 * γ4
    k = sqrt(max((γ1 - γ2) * (γ1 + γ2), k_min))

    exp_minusktau = exp(-τ * k)
    exp_minus2ktau = exp_minusktau * exp_minusktau

    RT_term = FT(1) / (k * (FT(1) + exp_minus2ktau) + γ1 * (FT(1) - exp_minus2ktau))

    Rdif = RT_term * γ2 * (FT(1) - exp_minus2ktau)
    Tdif = RT_term * FT(2) * k * exp_minusktau

    T₀ = Tnoscat = exp(-τ / max(μ₀, eps(FT)))

    k_μ = k * μ₀
    k_γ3 = k * γ3
    k_γ4 = k * γ4

    if abs(FT(1) - k_μ * k_μ) ≥ eps(FT)
        RT_term = ssa * RT_term / (FT(1) - k_μ * k_μ)
    else
        RT_term = ssa * RT_term / eps(FT)
    end

    Rdir_unconstrained =
        RT_term * (
            (FT(1) - k_μ) * (α2 + k_γ3) - (FT(1) + k_μ) * (α2 - k_γ3) * exp_minus2ktau -
            FT(2) * (k_γ3 - α2 * k_μ) * exp_minusktau * T₀
        )
    Tdir_unconstrained =
        -RT_term * (
            (FT(1) + k_μ) * (α1 + k_γ4) * T₀ - (FT(1) - k_μ) * (α1 - k_γ4) * exp_minus2ktau * T₀ -
            FT(2) * (k_γ4 + α1 * k_μ) * exp_minusktau
        )
    Rdir = max(FT(0), min(Rdir_unconstrained, (FT(1) - T₀)))
    Tdir = max(FT(0), min(Tdir_unconstrained, (FT(1) - T₀ - Rdir)))
    return (Rdir, Tdir, Tnoscat, Rdif, Tdif)
end

"""
    rte_sw_2stream_canopy!(...)

Modified SW 2-stream solver that uses physically-computed γ₃ for canopy layers.
For layers 1:ncanlay, calls `canopy_sw_2stream_coeffs` with β_dir.
For layers ncanlay+1:nlay, calls the standard `sw_2stream_coeffs`.
"""
function rte_sw_2stream_canopy!(
    (; τ, ssa, g)::TwoStream,
    (; albedo, src)::SourceSW2Str,
    bcs_sw::SwBCs,
    (; flux_up, flux_dn, flux_dn_dir)::FluxSW,
    solar_frac::FT,
    igpt::Int,
    n_gpt::Int,
    ibnd::Int,
    nlev::Int,
    gcol::Int,
    ncanlay::Int,
    β_dir::FT,
) where {FT}
    nlay = nlev - 1
    @inbounds begin
        toa_flux = bcs_sw.toa_flux[gcol]
        sfc_alb_direct = bcs_sw.sfc_alb_direct[ibnd, gcol]
        μ₀ = bcs_sw.cos_zenith[gcol]
    end
    τ_sum = FT(0)
    for ilay in 1:nlay
        @inbounds τ_sum += τ[ilay, gcol]
    end
    flux_dn_dir_top = toa_flux * solar_frac * μ₀
    flux_dn_dir_bot = flux_dn_dir_top * exp(-τ_sum / max(μ₀, eps(FT)))
    @inbounds flux_dn_dir[1, gcol] = flux_dn_dir_bot
    sfc_source = flux_dn_dir_bot * sfc_alb_direct

    @inbounds flux_dn[nlev, gcol] = FT(0)
    @inbounds surface_albedo = albedo[1, gcol] = bcs_sw.sfc_alb_diffuse[ibnd, gcol]
    @inbounds src[1, gcol] = sfc_source
    τ_cum = τ_sum
    albedo_ilev, src_ilev = surface_albedo, sfc_source
    @inbounds for ilev in 1:nlay
        τ_ilev, ssa_ilev, g_ilev = τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        (Rdir, Tdir, _, Rdif, Tdif) = if ilev ≤ ncanlay
            canopy_sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀, β_dir)
        else
            sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀)
        end
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)
        albedo_ilevplus1 = Rdif + Tdif * Tdif * albedo_ilev * denom
        τ_cum -= τ_ilev
        τ_cum = max(τ_cum, FT(0))
        flux_dn_dir_ilevplus1 = flux_dn_dir_top * exp(-τ_cum / max(μ₀, eps(FT)))
        src_up_ilev = Rdir * flux_dn_dir_ilevplus1
        src_dn_ilev = Tdir * flux_dn_dir_ilevplus1
        src_ilevplus1 = src_up_ilev + Tdif * denom * (src_ilev + albedo_ilev * src_dn_ilev)
        albedo[ilev + 1, gcol], src[ilev + 1, gcol] = albedo_ilevplus1, src_ilevplus1
        albedo_ilev = albedo_ilevplus1
        src_ilev = src_ilevplus1
    end
    @inbounds flux_up[nlev, gcol] =
        flux_dn[nlev, gcol] * albedo[nlev, gcol] + src[nlev, gcol]

    @inbounds flux_dn_ilevplus1 = flux_dn[nlev, gcol]
    @inbounds flux_dn[nlev, gcol] += flux_dn_dir_top
    τ_cum = FT(0)
    ilev = nlay
    @inbounds while ilev ≥ 1
        τ_ilev, ssa_ilev, g_ilev = τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        albedo_ilev, src_ilev = albedo[ilev, gcol], src[ilev, gcol]
        (_, Tdir, _, Rdif, Tdif) = if ilev ≤ ncanlay
            canopy_sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀, β_dir)
        else
            sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀)
        end
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)
        src_dn_ilev = Tdir * flux_dn_dir_top * exp(-τ_cum / max(μ₀, eps(FT)))
        τ_cum += τ_ilev
        flux_dn_ilev = (Tdif * flux_dn_ilevplus1 + Rdif * src_ilev + src_dn_ilev) * denom
        flux_up[ilev, gcol] = flux_dn_ilev * albedo_ilev + src_ilev
        flux_dn[ilev, gcol] = flux_dn_ilev + flux_dn_dir_top * exp(-τ_cum / max(μ₀, eps(FT)))
        flux_dn_ilevplus1 = flux_dn_ilev
        ilev -= 1
    end
    return nothing
end

# Compute G(μ₀) using CanopyOptics leaf angle distribution
ld = spherical_leaves()  # G = 0.5 for spherical
μ₀ = cosd(sza_deg)
G_vals = G([FT(μ₀)], ld)
G_μ₀ = FT(G_vals[1])
LAI_layer = FT(LAI_total / ncanlay)

# Canopy temperature profile (levels 1:ncanlay, bottom-to-top within canopy)
T_canopy_lev = [canopy_temp_profile(ilev, ncanlay, T_sfc_val, T_air_val) for ilev in 1:ncanlay]

println("\n--- Canopy 2-stream mapping ---")
println("  Leaf angle distribution: spherical")
println("  G(μ₀=$(round(μ₀,digits=3))) = $(round(G_μ₀,digits=4))")
println("  LAI/layer = $(LAI_layer)")
println("\n  Canopy temperature profile (bottom → top):")
for ilev in 1:ncanlay
    println("    Level $ilev: $(round(T_canopy_lev[ilev], digits=2)) K")
end
println("    Level $(ncanlay+1): $(round(T_air_val, digits=2)) K (atmosphere boundary)")

# Compute directional backscatter fraction β₀(μ₀) per band
β_dir_band = zeros(FT, nbnd)
for ibnd in 1:nbnd
    if band_in_prospect[ibnd]
        β_dir_band[ibnd] = compute_β_dir(FT(μ₀), FT(leaf_R_band[ibnd]), FT(leaf_T_band[ibnd]), ld)
    else
        β_dir_band[ibnd] = FT(0.5)  # symmetric (bands are transparent anyway)
    end
end

println("\n  Directional vs hemispheric backscatter (β):")
println("    Band | β_hemi (R/ω) | β_dir(μ₀) |   Δβ")
println("    -----+--------------+-----------+-------")
for ibnd in 1:nbnd
    if band_in_prospect[ibnd]
        ω = leaf_R_band[ibnd] + leaf_T_band[ibnd]
        β_hemi = ω > 0 ? leaf_R_band[ibnd] / ω : 0.5
        Δβ = β_dir_band[ibnd] - β_hemi
        println("      $(lpad(ibnd, 2)) |     $(lpad(round(β_hemi, digits=4), 6)) |    $(lpad(round(β_dir_band[ibnd], digits=4), 6)) | $(lpad(round(Δβ, digits=4), 7))")
    end
end

# ============================================================================
# Section 5: Create expanded column data structures
# ============================================================================
println("\n--- Building expanded column ---")

nlay_total = nlay_atm + ncanlay
nlev_total = nlay_total + 1

println("  Total layers: $nlay_total (atm: $nlay_atm, canopy: $ncanlay)")

# Grid params for expanded and atmospheric columns
grid_exp = RRTMGPGridParams(FT; context, nlay = nlay_total, ncol)
grid_atm = RRTMGPGridParams(FT; context, nlay = nlay_atm, ncol)

# Expanded column optical properties, sources, fluxes
op_exp = TwoStream(grid_exp)
src_exp = SourceSW2Str(grid_exp)
flux_gpt = FluxSW(grid_exp)   # per-gpoint flux
flux_sw = FluxSW(grid_exp)    # accumulated flux

# Atmospheric column optical properties (for compute_optical_props!)
op_atm = TwoStream(grid_atm)

# Boundary conditions with soil albedo
soil_alb_direct = DA{FT, 2}(fill(soil_alb_val, nbnd, ncol))
soil_alb_diffuse = DA{FT, 2}(fill(soil_alb_val, nbnd, ncol))
bcs_sw = SwBCs(cos_zenith, toa_flux, soil_alb_direct, nothing, soil_alb_diffuse)

# Also create reference solver (atmosphere-only, original surface albedo)
# for validation comparison
sfc_alb_direct_ref = DA{FT, 2}(fill(soil_alb_val, nbnd, ncol))
sfc_alb_diffuse_ref = DA{FT, 2}(fill(soil_alb_val, nbnd, ncol))
bcs_ref = SwBCs(cos_zenith, toa_flux, sfc_alb_direct_ref, nothing, sfc_alb_diffuse_ref)
src_ref = SourceSW2Str(grid_atm)
op_ref = TwoStream(grid_atm)
flux_gpt_ref = FluxSW(grid_atm)
flux_ref = FluxSW(grid_atm)

# APAR accumulators
apar_direct = zeros(FT, ncanlay, ncol)
apar_diffuse = zeros(FT, ncanlay, ncol)

println("  Expanded TwoStream: $(size(op_exp.τ))")
println("  Expanded FluxSW: $(size(flux_sw.flux_up))")

# --- LW data structures ---
op_exp_lw = TwoStream(grid_exp)
sf_exp_lw = SourceLW2Str(grid_exp; params = param_set)
flux_gpt_lw = FluxLW(grid_exp)
flux_lw = FluxLW(grid_exp)

op_atm_lw = TwoStream(grid_atm)
sf_atm_lw = SourceLW2Str(grid_atm; params = param_set)

soil_emis = DA{FT, 2}(fill(soil_emis_val, nbnd_lw, ncol))
bcs_lw = LwBCs(soil_emis, nothing)

# Reference LW solve (no canopy)
op_ref_lw = TwoStream(grid_atm)
sf_ref_lw = SourceLW2Str(grid_atm; params = param_set)
flux_gpt_ref_lw = FluxLW(grid_atm)
flux_ref_lw = FluxLW(grid_atm)
bcs_ref_lw = LwBCs(DA{FT, 2}(fill(soil_emis_val, nbnd_lw, ncol)), nothing)

# Per-layer LW absorption accumulator
lw_abs_canopy = zeros(FT, ncanlay, ncol)

println("  LW TwoStream: $(size(op_exp_lw.τ))")
println("  LW FluxLW: $(size(flux_lw.flux_up))")

# ============================================================================
# Section 6: Main g-point solve loop
# ============================================================================
println("\n--- Running expanded column solver ---")

flux_up_sw = flux_sw.flux_up
flux_dn_sw = flux_sw.flux_dn
flux_dn_dir_sw = flux_sw.flux_dn_dir
flux_net_sw = flux_sw.flux_net

flux_up_ref_sw = flux_ref.flux_up
flux_dn_ref_sw = flux_ref.flux_dn
flux_dn_dir_ref_sw = flux_ref.flux_dn_dir
flux_net_ref_sw = flux_ref.flux_net

@inbounds for igpt in 1:n_gpt
    ibnd = major_gpt2bnd[igpt]
    solar_frac = lookup_sw.solar_src_scaled[igpt]

    for gcol in 1:ncol
        μ₀_col = cos_zenith[gcol]
        if μ₀_col ≤ 0
            Fluxes.set_flux_to_zero!(flux_sw, gcol)
            Fluxes.set_flux_to_zero!(flux_ref, gcol)
            continue
        end

        # --- Atmospheric optics (into op_atm, nlay_atm layers) ---
        compute_optical_props!(op_atm, as, gcol, igpt, lookup_sw, nothing, nothing)

        # === Reference solve (atmosphere only, no canopy) ===
        # Copy atmospheric optics into reference op
        for ilay in 1:nlay_atm
            op_ref.τ[ilay, gcol] = op_atm.τ[ilay, gcol]
            op_ref.ssa[ilay, gcol] = op_atm.ssa[ilay, gcol]
            op_ref.g[ilay, gcol] = op_atm.g[ilay, gcol]
        end
        rte_sw_2stream!(op_ref, src_ref, bcs_ref, flux_gpt_ref, solar_frac, igpt, n_gpt, ibnd, nlev_atm, gcol)

        # === Expanded solve (atmosphere + canopy) ===
        # Copy atmospheric optics to upper layers of expanded column
        for ilay in 1:nlay_atm
            op_exp.τ[ilay + ncanlay, gcol] = op_atm.τ[ilay, gcol]
            op_exp.ssa[ilay + ncanlay, gcol] = op_atm.ssa[ilay, gcol]
            op_exp.g[ilay + ncanlay, gcol] = op_atm.g[ilay, gcol]
        end

        # Fill canopy layers (bottom ncanlay layers)
        R_b = leaf_R_band[ibnd]
        T_b = leaf_T_band[ibnd]
        if band_in_prospect[ibnd]
            τ_c, ω_c, g_c = canopy_2stream_props(FT(R_b), FT(T_b), G_μ₀, LAI_layer)
        else
            # Outside PROSPECT range: transparent canopy
            τ_c, ω_c, g_c = FT(0), FT(0), FT(0)
        end
        for ilay in 1:ncanlay
            op_exp.τ[ilay, gcol] = τ_c
            op_exp.ssa[ilay, gcol] = ω_c
            op_exp.g[ilay, gcol] = g_c
        end

        # Solve expanded column (canopy layers use β_dir for direct beam geometry)
        rte_sw_2stream_canopy!(op_exp, src_exp, bcs_sw, flux_gpt, solar_frac, igpt, n_gpt, ibnd, nlev_total, gcol, ncanlay, β_dir_band[ibnd])

        # --- Accumulate fluxes ---
        if igpt == 1
            for ilev in 1:nlev_total
                flux_up_sw[ilev, gcol] = flux_gpt.flux_up[ilev, gcol]
                flux_dn_sw[ilev, gcol] = flux_gpt.flux_dn[ilev, gcol]
            end
            flux_dn_dir_sw[1, gcol] = flux_gpt.flux_dn_dir[1, gcol]
            # Reference
            for ilev in 1:nlev_atm
                flux_up_ref_sw[ilev, gcol] = flux_gpt_ref.flux_up[ilev, gcol]
                flux_dn_ref_sw[ilev, gcol] = flux_gpt_ref.flux_dn[ilev, gcol]
            end
            flux_dn_dir_ref_sw[1, gcol] = flux_gpt_ref.flux_dn_dir[1, gcol]
        else
            for ilev in 1:nlev_total
                flux_up_sw[ilev, gcol] += flux_gpt.flux_up[ilev, gcol]
                flux_dn_sw[ilev, gcol] += flux_gpt.flux_dn[ilev, gcol]
            end
            flux_dn_dir_sw[1, gcol] += flux_gpt.flux_dn_dir[1, gcol]
            # Reference
            for ilev in 1:nlev_atm
                flux_up_ref_sw[ilev, gcol] += flux_gpt_ref.flux_up[ilev, gcol]
                flux_dn_ref_sw[ilev, gcol] += flux_gpt_ref.flux_dn[ilev, gcol]
            end
            flux_dn_dir_ref_sw[1, gcol] += flux_gpt_ref.flux_dn_dir[1, gcol]
        end

        # --- APAR extraction for this g-point ---
        # Absorption per canopy layer = flux convergence
        # flux_dn[ilev] is at the TOP of layer ilev, flux_dn[ilev-1] at bottom (or ilev+1 going up)
        # In RRTMGP: level ilev is at TOP of layer ilev (going bottom-to-top)
        #   level 1 = bottom of layer 1 = soil surface
        #   level 2 = top of layer 1 = bottom of layer 2
        #   level ncanlay+1 = top of canopy
        for ilay in 1:ncanlay
            # Flux convergence: what comes in minus what goes out
            # Level ilay is at bottom of layer ilay, level ilay+1 is at top
            flux_in = flux_gpt.flux_dn[ilay + 1, gcol] + flux_gpt.flux_up[ilay, gcol]
            flux_out = flux_gpt.flux_dn[ilay, gcol] + flux_gpt.flux_up[ilay + 1, gcol]
            absorbed_total = flux_in - flux_out

            # Direct beam component: Beer-Lambert extinction through canopy
            # Cumulative τ from TOA to top of this layer
            τ_above = FT(0)
            for j in (ilay + 1):nlay_total
                τ_above += op_exp.τ[j, gcol]
            end
            # Also τ from top of this layer to bottom of this layer
            τ_this = op_exp.τ[ilay, gcol]

            flux_dn_dir_top_of_layer = toa_flux[gcol] * solar_frac * μ₀_col * exp(-τ_above / μ₀_col)
            flux_dn_dir_bot_of_layer = flux_dn_dir_top_of_layer * exp(-τ_this / μ₀_col)

            # Direct beam extinguished in this layer
            direct_extinguished = flux_dn_dir_top_of_layer - flux_dn_dir_bot_of_layer
            # Of that, fraction (1-ω) is absorbed, rest is scattered
            direct_absorbed = direct_extinguished * (FT(1) - ω_c)

            diffuse_absorbed = absorbed_total - direct_absorbed

            apar_direct[ilay, gcol] += direct_absorbed
            apar_diffuse[ilay, gcol] += diffuse_absorbed
        end
    end
end

# Compute net fluxes
for gcol in 1:ncol
    if cos_zenith[gcol] > 0
        for ilev in 1:nlev_total
            flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
        end
        for ilev in 1:nlev_atm
            flux_net_ref_sw[ilev, gcol] = flux_up_ref_sw[ilev, gcol] - flux_dn_ref_sw[ilev, gcol]
        end
    end
end

println("  Solver complete.")

# ============================================================================
# Section 7: Results
# ============================================================================
println("\n" * "="^60)
println("RESULTS")
println("="^60)

gcol = 1
flux_up_arr = Array(flux_up_sw)
flux_dn_arr = Array(flux_dn_sw)
flux_dn_dir_arr = Array(flux_dn_dir_sw)

# TOA = level nlev_total (top of expanded column)
toa_incoming = flux_dn_arr[nlev_total, gcol]
toa_reflected = flux_up_arr[nlev_total, gcol]

# Canopy top = level ncanlay+1
canopy_top_dn = flux_dn_arr[ncanlay + 1, gcol]
canopy_top_up = flux_up_arr[ncanlay + 1, gcol]

# Soil surface = level 1
soil_dn = flux_dn_arr[1, gcol]
soil_up = flux_up_arr[1, gcol]
soil_absorbed = soil_dn - soil_up

# Direct beam at soil
soil_dir = flux_dn_dir_arr[1, gcol]

println("\n--- Flux Profile (W/m²) ---")
println("  TOA incoming:       $(round(toa_incoming, digits=2))")
println("  TOA reflected:      $(round(toa_reflected, digits=2))")
println("  Canopy top ↓:       $(round(canopy_top_dn, digits=2))")
println("  Canopy top ↑:       $(round(canopy_top_up, digits=2))")
println("  Soil ↓:             $(round(soil_dn, digits=2))")
println("  Soil ↑:             $(round(soil_up, digits=2))")
println("  Soil ↓ (direct):    $(round(soil_dir, digits=2))")
println("  Soil absorbed:      $(round(soil_absorbed, digits=2))")

# APAR results
total_apar = sum(apar_direct[:, gcol]) + sum(apar_diffuse[:, gcol])
println("\n--- APAR by Canopy Layer (W/m²) ---")
println("  Layer | Direct  | Diffuse | Total")
println("  ------+---------+---------+--------")
for ilay in ncanlay:-1:1  # print top-to-bottom
    d = apar_direct[ilay, gcol]
    f = apar_diffuse[ilay, gcol]
    println("    $(lpad(ilay, 2))  | $(lpad(round(d, digits=2), 7)) | $(lpad(round(f, digits=2), 7)) | $(lpad(round(d+f, digits=2), 7))")
end
println("  ------+---------+---------+--------")
println("  Total | $(lpad(round(sum(apar_direct[:,gcol]), digits=2), 7)) " *
        "| $(lpad(round(sum(apar_diffuse[:,gcol]), digits=2), 7)) " *
        "| $(lpad(round(total_apar, digits=2), 7))")

# Atmospheric absorption
atm_absorbed = sum(
    (flux_dn_arr[ilay + 1, gcol] + flux_up_arr[ilay, gcol]) -
    (flux_dn_arr[ilay, gcol] + flux_up_arr[ilay + 1, gcol])
    for ilay in (ncanlay + 1):nlay_total
)

# ============================================================================
# Section 8: Validation
# ============================================================================
println("\n--- Energy Conservation Check ---")
residual = toa_incoming - toa_reflected - atm_absorbed - total_apar - soil_absorbed
println("  TOA incoming:         $(round(toa_incoming, digits=4))")
println("  TOA reflected:        $(round(toa_reflected, digits=4))")
println("  Atm absorbed:         $(round(atm_absorbed, digits=4))")
println("  Canopy APAR:          $(round(total_apar, digits=4))")
println("  Soil absorbed:        $(round(soil_absorbed, digits=4))")
println("  Residual:             $(round(residual, digits=6))")
println("  Energy conserved:     $(abs(residual) < 1e-4 ? "YES ✓" : "NO ✗")")

# Reference comparison (atmosphere only, same surface albedo)
println("\n--- Reference Comparison (no canopy) ---")
flux_up_ref_arr = Array(flux_up_ref_sw)
flux_dn_ref_arr = Array(flux_dn_ref_sw)

ref_toa_dn = flux_dn_ref_arr[nlev_atm, gcol]
ref_toa_up = flux_up_ref_arr[nlev_atm, gcol]
ref_sfc_dn = flux_dn_ref_arr[1, gcol]
ref_sfc_up = flux_up_ref_arr[1, gcol]

println("  Reference TOA ↓:     $(round(ref_toa_dn, digits=2))")
println("  Reference TOA ↑:     $(round(ref_toa_up, digits=2))")
println("  Reference Sfc ↓:     $(round(ref_sfc_dn, digits=2))")
println("  Reference Sfc ↑:     $(round(ref_sfc_up, digits=2))")
println("  With canopy TOA ↑:   $(round(toa_reflected, digits=2))")
println("  Canopy reduces TOA ↑ by: $(round(ref_toa_up - toa_reflected, digits=2)) W/m²")

# Fraction of incoming SW absorbed by canopy
frac_apar = total_apar / toa_incoming * 100
println("\n--- Summary ---")
println("  Canopy APAR fraction: $(round(frac_apar, digits=1))% of incoming SW")
println("  Direct/Diffuse APAR:  $(round(sum(apar_direct[:,gcol])/total_apar*100, digits=1))% / " *
        "$(round(sum(apar_diffuse[:,gcol])/total_apar*100, digits=1))%")
println("  Canopy albedo:        $(round(canopy_top_up / canopy_top_dn, digits=4))")

# ============================================================================
# Section 9: LW g-point solve loop
# ============================================================================
println("\n" * "="^60)
println("LONGWAVE RADIATIVE TRANSFER")
println("="^60)
println("\n--- Running LW expanded column solver ---")
println("  ε_leaf = $ε_leaf, soil_emis = $soil_emis_val")
println("  ω_lw = $ω_lw, g_lw = $g_lw, G_hemi = $G_hemi")

flux_up_lw = flux_lw.flux_up
flux_dn_lw = flux_lw.flux_dn
flux_net_lw = flux_lw.flux_net

flux_up_ref_lw = flux_ref_lw.flux_up
flux_dn_ref_lw = flux_ref_lw.flux_dn
flux_net_ref_lw = flux_ref_lw.flux_net

τ_lw_canopy = G_hemi * LAI_layer

@inbounds for igpt in 1:n_gpt_lw
    ibnd = major_gpt2bnd_lw[igpt]

    for gcol in 1:ncol
        # --- Atmospheric optics + Planck sources ---
        compute_optical_props!(op_atm_lw, as, sf_atm_lw, gcol, igpt, lookup_lw, nothing, nothing)

        # === Reference solve (atmosphere only, no canopy) ===
        for ilay in 1:nlay_atm
            op_ref_lw.τ[ilay, gcol] = op_atm_lw.τ[ilay, gcol]
            op_ref_lw.ssa[ilay, gcol] = op_atm_lw.ssa[ilay, gcol]
            op_ref_lw.g[ilay, gcol] = op_atm_lw.g[ilay, gcol]
        end
        for ilev in 1:nlev_atm
            sf_ref_lw.lev_source[ilev, gcol] = sf_atm_lw.lev_source[ilev, gcol]
        end
        sf_ref_lw.sfc_source[gcol] = sf_atm_lw.sfc_source[gcol]
        rte_lw_2stream!(op_ref_lw, flux_gpt_ref_lw, sf_ref_lw, bcs_ref_lw, gcol, igpt, ibnd, nlev_atm, ncol)

        # === Expanded solve (atmosphere + canopy) ===
        # Copy atmospheric optics to upper layers
        for ilay in 1:nlay_atm
            op_exp_lw.τ[ilay + ncanlay, gcol] = op_atm_lw.τ[ilay, gcol]
            op_exp_lw.ssa[ilay + ncanlay, gcol] = op_atm_lw.ssa[ilay, gcol]
            op_exp_lw.g[ilay + ncanlay, gcol] = op_atm_lw.g[ilay, gcol]
        end

        # Copy atmospheric Planck sources to expanded column
        for atm_ilev in 1:nlev_atm
            sf_exp_lw.lev_source[ncanlay + atm_ilev, gcol] = sf_atm_lw.lev_source[atm_ilev, gcol]
        end

        # Fill canopy LW optics (layers 1:ncanlay)
        for ilay in 1:ncanlay
            op_exp_lw.τ[ilay, gcol] = τ_lw_canopy
            op_exp_lw.ssa[ilay, gcol] = ω_lw
            op_exp_lw.g[ilay, gcol] = g_lw
        end

        # Fill canopy Planck sources with temperature gradient
        # Scale by B(T_canopy)/B(T_ref) to preserve planckfrac weighting
        totplnk_band = view(tot_planck_cpu, :, ibnd)
        T_ref = FT(as.t_lev[1, gcol])
        B_ref = planck_interp(T_ref, t_planck_cpu, totplnk_band)
        for ilev in 1:ncanlay
            B_can = planck_interp(T_canopy_lev[ilev], t_planck_cpu, totplnk_band)
            sf_exp_lw.lev_source[ilev, gcol] = sf_atm_lw.lev_source[1, gcol] * (B_can / B_ref)
        end
        # Level ncanlay+1 already set from atmospheric data (at T_ref), no scaling needed

        # Surface Planck source
        sf_exp_lw.sfc_source[gcol] = sf_atm_lw.sfc_source[gcol]

        # Solve expanded column
        rte_lw_2stream!(op_exp_lw, flux_gpt_lw, sf_exp_lw, bcs_lw, gcol, igpt, ibnd, nlev_total, ncol)

        # --- Accumulate fluxes ---
        if igpt == 1
            for ilev in 1:nlev_total
                flux_up_lw[ilev, gcol] = flux_gpt_lw.flux_up[ilev, gcol]
                flux_dn_lw[ilev, gcol] = flux_gpt_lw.flux_dn[ilev, gcol]
            end
            for ilev in 1:nlev_atm
                flux_up_ref_lw[ilev, gcol] = flux_gpt_ref_lw.flux_up[ilev, gcol]
                flux_dn_ref_lw[ilev, gcol] = flux_gpt_ref_lw.flux_dn[ilev, gcol]
            end
        else
            for ilev in 1:nlev_total
                flux_up_lw[ilev, gcol] += flux_gpt_lw.flux_up[ilev, gcol]
                flux_dn_lw[ilev, gcol] += flux_gpt_lw.flux_dn[ilev, gcol]
            end
            for ilev in 1:nlev_atm
                flux_up_ref_lw[ilev, gcol] += flux_gpt_ref_lw.flux_up[ilev, gcol]
                flux_dn_ref_lw[ilev, gcol] += flux_gpt_ref_lw.flux_dn[ilev, gcol]
            end
        end

        # --- Per-layer LW absorption (flux convergence) ---
        for ilay in 1:ncanlay
            flux_in = flux_gpt_lw.flux_dn[ilay + 1, gcol] + flux_gpt_lw.flux_up[ilay, gcol]
            flux_out = flux_gpt_lw.flux_dn[ilay, gcol] + flux_gpt_lw.flux_up[ilay + 1, gcol]
            lw_abs_canopy[ilay, gcol] += flux_in - flux_out
        end
    end
end

# Compute LW net fluxes
for gcol in 1:ncol
    for ilev in 1:nlev_total
        flux_net_lw[ilev, gcol] = flux_up_lw[ilev, gcol] - flux_dn_lw[ilev, gcol]
    end
    for ilev in 1:nlev_atm
        flux_net_ref_lw[ilev, gcol] = flux_up_ref_lw[ilev, gcol] - flux_dn_ref_lw[ilev, gcol]
    end
end

println("  LW solver complete.")

# ============================================================================
# Section 10: LW Results
# ============================================================================
println("\n" * "="^60)
println("LW RESULTS")
println("="^60)

gcol = 1
flux_up_lw_arr = Array(flux_up_lw)
flux_dn_lw_arr = Array(flux_dn_lw)

# TOA = level nlev_total
lw_toa_dn = flux_dn_lw_arr[nlev_total, gcol]
lw_toa_up = flux_up_lw_arr[nlev_total, gcol]  # OLR

# Canopy top = level ncanlay+1
lw_canopy_top_dn = flux_dn_lw_arr[ncanlay + 1, gcol]
lw_canopy_top_up = flux_up_lw_arr[ncanlay + 1, gcol]

# Soil surface = level 1
lw_soil_dn = flux_dn_lw_arr[1, gcol]
lw_soil_up = flux_up_lw_arr[1, gcol]
lw_soil_absorbed = lw_soil_dn - lw_soil_up

# Atmospheric absorption
lw_atm_absorbed = sum(
    (flux_dn_lw_arr[ilay + 1, gcol] + flux_up_lw_arr[ilay, gcol]) -
    (flux_dn_lw_arr[ilay, gcol] + flux_up_lw_arr[ilay + 1, gcol])
    for ilay in (ncanlay + 1):nlay_total
)

println("\n--- LW Flux Profile (W/m²) ---")
println("  TOA ↓ (incident):   $(round(lw_toa_dn, digits=2))")
println("  TOA ↑ (OLR):        $(round(lw_toa_up, digits=2))")
println("  Canopy top ↓:       $(round(lw_canopy_top_dn, digits=2))")
println("  Canopy top ↑:       $(round(lw_canopy_top_up, digits=2))")
println("  Soil ↓:             $(round(lw_soil_dn, digits=2))")
println("  Soil ↑:             $(round(lw_soil_up, digits=2))")
println("  Soil net absorbed:  $(round(lw_soil_absorbed, digits=2))")

# Per-layer LW absorption
total_lw_canopy = sum(lw_abs_canopy[:, gcol])
println("\n--- LW Absorbed per Canopy Layer (W/m²) ---")
println("  Layer | Net LW absorbed")
println("  ------+----------------")
for ilay in ncanlay:-1:1  # print top-to-bottom
    lw_abs = lw_abs_canopy[ilay, gcol]
    println("    $(lpad(ilay, 2))  | $(lpad(round(lw_abs, digits=2), 10))")
end
println("  ------+----------------")
println("  Total | $(lpad(round(total_lw_canopy, digits=2), 10))")

# ============================================================================
# Section 11: LW Validation
# ============================================================================
println("\n--- LW Energy Conservation ---")
lw_residual = lw_toa_dn - lw_toa_up - lw_atm_absorbed - total_lw_canopy - lw_soil_absorbed
println("  TOA ↓:                $(round(lw_toa_dn, digits=4))")
println("  TOA ↑ (OLR):          $(round(lw_toa_up, digits=4))")
println("  Atm absorbed:         $(round(lw_atm_absorbed, digits=4))")
println("  Canopy absorbed:      $(round(total_lw_canopy, digits=4))")
println("  Soil absorbed:        $(round(lw_soil_absorbed, digits=4))")
println("  Residual:             $(round(lw_residual, digits=6))")
println("  Energy conserved:     $(abs(lw_residual) < 1e-4 ? "YES ✓" : "NO ✗")")

# Reference comparison
println("\n--- LW Reference Comparison (no canopy) ---")
flux_up_ref_lw_arr = Array(flux_up_ref_lw)
flux_dn_ref_lw_arr = Array(flux_dn_ref_lw)

ref_lw_olr = flux_up_ref_lw_arr[nlev_atm, gcol]
ref_lw_sfc_dn = flux_dn_ref_lw_arr[1, gcol]
ref_lw_sfc_up = flux_up_ref_lw_arr[1, gcol]

println("  Reference OLR:        $(round(ref_lw_olr, digits=2))")
println("  Reference Sfc ↓:      $(round(ref_lw_sfc_dn, digits=2))")
println("  Reference Sfc ↑:      $(round(ref_lw_sfc_up, digits=2))")
println("  With canopy OLR:      $(round(lw_toa_up, digits=2))")
println("  OLR change:           $(round(lw_toa_up - ref_lw_olr, digits=2)) W/m²")

# ============================================================================
# Section 12: Combined SW + LW Summary
# ============================================================================
println("\n" * "="^60)
println("COMBINED SW + LW SUMMARY")
println("="^60)

total_sw_canopy = total_apar
println("\n--- Net Radiation per Canopy Layer (W/m²) ---")
println("  Layer | SW absorbed | LW absorbed | Net (SW+LW)")
println("  ------+-------------+-------------+------------")
for ilay in ncanlay:-1:1
    sw_abs = apar_direct[ilay, gcol] + apar_diffuse[ilay, gcol]
    lw_abs = lw_abs_canopy[ilay, gcol]
    net = sw_abs + lw_abs
    println("    $(lpad(ilay, 2))  | $(lpad(round(sw_abs, digits=2), 11)) | $(lpad(round(lw_abs, digits=2), 11)) | $(lpad(round(net, digits=2), 10))")
end
println("  ------+-------------+-------------+------------")
println("  Total | $(lpad(round(total_sw_canopy, digits=2), 11)) | $(lpad(round(total_lw_canopy, digits=2), 11)) | $(lpad(round(total_sw_canopy + total_lw_canopy, digits=2), 10))")

println("\n--- Column Summary ---")
println("  SW TOA incoming:      $(round(toa_incoming, digits=2))")
println("  SW TOA reflected:     $(round(toa_reflected, digits=2))")
println("  LW OLR:               $(round(lw_toa_up, digits=2))")
println("  Net radiation (TOA):  $(round(toa_incoming - toa_reflected - lw_toa_up, digits=2))")
