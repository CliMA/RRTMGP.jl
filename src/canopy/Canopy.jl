module Canopy

using Adapt
using DocStringExtensions
import ClimaComms

import ..RRTMGPGridParams
import ..Parameters as RP
using ..Sources
using ..Optics

export CanopyState,
    CanopyRadiationOutput,
    canopy_2stream_props,
    canopy_sw_2stream_coeffs,
    compute_beta_dir,
    planck_interp,
    fill_canopy_sw_optics!,
    fill_canopy_lw_optics!,
    fill_canopy_lw_sources!,
    compute_canopy_radiation!

"""
    CanopyState{FT, V1, V2, V2B}

State container for canopy radiative transfer within an expanded atmospheric column.

Canopy layers sit below the atmospheric layers (index 1 = soil surface).
All arrays are allocated on the appropriate device via `RRTMGPGridParams`.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CanopyState{FT, V1, V2, V2B}
    "number of canopy sublayers"
    ncanlay::Int
    "LAI per sublayer (ncanlay, ncol)"
    lai_layer::V2
    "band-averaged leaf reflectance (nbnd_sw, ncol)"
    leaf_r::V2B
    "band-averaged leaf transmittance (nbnd_sw, ncol)"
    leaf_t::V2B
    "LW leaf single-scattering albedo per band (nbnd_lw, ncol)"
    leaf_omega_lw::V2B
    "Ross G-function at cos(SZA) per column (ncol)"
    G_dir::V1
    "hemispheric-mean G per column (ncol)"
    G_hemi::V1
    "direct-beam backscatter fraction per SW band (nbnd_sw, ncol)"
    beta_dir::V2B
    "canopy level temperatures, bottom to top (ncanlay+1, ncol)"
    t_canopy_lev::V2
    "soil direct albedo per SW band (nbnd_sw, ncol)"
    soil_alb_direct::V2B
    "soil diffuse albedo per SW band (nbnd_sw, ncol)"
    soil_alb_diffuse::V2B
    "soil emissivity per LW band (nbnd_lw, ncol)"
    soil_emis::V2B
    "soil surface temperature (ncol)"
    t_soil::V1
    "chlorophyll absorption fraction per SW band (nbnd_sw, ncol)"
    f_cab_band::V2B
    "carotenoid absorption fraction per SW band (nbnd_sw, ncol)"
    f_car_band::V2B
end
Adapt.@adapt_structure CanopyState

"""
    CanopyState(grid_params::RRTMGPGridParams; ncanlay, nbnd_sw, nbnd_lw)

Construct a `CanopyState` with all arrays allocated on the correct device.
Arrays are initialized to zero and must be filled by the user before solving.
"""
function CanopyState(grid_params::RRTMGPGridParams; ncanlay::Int, nbnd_sw::Int, nbnd_lw::Int)
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    ncol = grid_params.ncol

    lai_layer = DA{FT, 2}(zeros(FT, ncanlay, ncol))
    leaf_r = DA{FT, 2}(zeros(FT, nbnd_sw, ncol))
    leaf_t = DA{FT, 2}(zeros(FT, nbnd_sw, ncol))
    leaf_omega_lw = DA{FT, 2}(zeros(FT, nbnd_lw, ncol))
    G_dir = DA{FT, 1}(zeros(FT, ncol))
    G_hemi = DA{FT, 1}(zeros(FT, ncol))
    beta_dir = DA{FT, 2}(zeros(FT, nbnd_sw, ncol))
    t_canopy_lev = DA{FT, 2}(zeros(FT, ncanlay + 1, ncol))
    soil_alb_direct = DA{FT, 2}(zeros(FT, nbnd_sw, ncol))
    soil_alb_diffuse = DA{FT, 2}(zeros(FT, nbnd_sw, ncol))
    soil_emis = DA{FT, 2}(zeros(FT, nbnd_lw, ncol))
    t_soil = DA{FT, 1}(zeros(FT, ncol))
    f_cab_band = DA{FT, 2}(zeros(FT, nbnd_sw, ncol))
    f_car_band = DA{FT, 2}(zeros(FT, nbnd_sw, ncol))

    return CanopyState{FT, typeof(G_dir), typeof(lai_layer), typeof(leaf_r)}(
        ncanlay,
        lai_layer,
        leaf_r,
        leaf_t,
        leaf_omega_lw,
        G_dir,
        G_hemi,
        beta_dir,
        t_canopy_lev,
        soil_alb_direct,
        soil_alb_diffuse,
        soil_emis,
        t_soil,
        f_cab_band,
        f_car_band,
    )
end

using ..Fluxes

include("canopy_utils.jl")
include("canopy_optics.jl")
include("canopy_diagnostics.jl")

# Stubs for extension functions (implemented in RRTMGPCanopyOpticsExt)
"""
    fill_leaf_optics!(canopy, leaf, opti, bnd_lims_wn)

Band-average PROSPECT leaf R/T onto RRTMGP bands. Requires `CanopyOptics.jl`.
"""
function fill_leaf_optics! end

"""
    fill_G!(canopy, ld, cos_zenith)

Fill G_dir and G_hemi from a leaf angle distribution. Requires `CanopyOptics.jl`.
"""
function fill_G! end

"""
    fill_beta_dir!(canopy, ld, cos_zenith; nθ=90, nφ=360)

Precompute directional backscatter per band. Requires `CanopyOptics.jl`.
"""
function fill_beta_dir! end

end # module
