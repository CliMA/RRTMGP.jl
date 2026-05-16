"""
    CanopyRadiationOutput{FT, V1, V2, V3}

Diagnostic output for canopy radiation absorption, fAPAR, and green-APAR.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CanopyRadiationOutput{FT, V1, V2, V3}
    "SW absorbed per canopy layer — total (ncanlay, ncol)"
    apar::V2
    "SW absorbed per canopy layer per band (nbnd_sw, ncanlay, ncol)"
    apar_by_band::V3
    "LW net absorbed per canopy layer (ncanlay, ncol)"
    lw_absorbed::V2
    "fAPAR: fraction of incoming PAR absorbed by canopy (ncol)"
    fAPAR::V1
    "green-APAR per canopy layer (ncanlay, ncol)"
    green_apar::V2
    "carotenoid energy transfer efficiency [0-1]"
    η_car::FT
end
Adapt.@adapt_structure CanopyRadiationOutput

"""
    CanopyRadiationOutput(grid_params; ncanlay, nbnd_sw, η_car=0.0)

Construct a `CanopyRadiationOutput` with arrays allocated on the correct device.
"""
function CanopyRadiationOutput(grid_params::RRTMGPGridParams; ncanlay::Int, nbnd_sw::Int, η_car = nothing)
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    ncol = grid_params.ncol

    η = isnothing(η_car) ? FT(0) : FT(η_car)
    apar = DA{FT, 2}(zeros(FT, ncanlay, ncol))
    apar_by_band = DA{FT, 3}(zeros(FT, nbnd_sw, ncanlay, ncol))
    lw_absorbed = DA{FT, 2}(zeros(FT, ncanlay, ncol))
    fAPAR = DA{FT, 1}(zeros(FT, ncol))
    green_apar = DA{FT, 2}(zeros(FT, ncanlay, ncol))

    return CanopyRadiationOutput{FT, typeof(fAPAR), typeof(apar), typeof(apar_by_band)}(
        apar, apar_by_band, lw_absorbed, fAPAR, green_apar, η,
    )
end

"""
    compute_canopy_radiation!(output, flux_sw, flux_lw, canopy, nlay_atm, bnd_lims_wn)

Compute fAPAR and green-APAR from accumulated flux profiles.

`flux_sw` and `flux_lw` are the accumulated fluxes from the expanded-column solve.
`bnd_lims_wn` is `(2, nbnd_sw)` band wavenumber limits used to identify PAR bands.
"""
function compute_canopy_radiation!(
    output::CanopyRadiationOutput,
    flux_sw::FluxSW,
    flux_lw::Union{FluxLW, Nothing},
    canopy::CanopyState,
    nlay_atm::Int,
    bnd_lims_wn::AbstractMatrix,
)
    FT = eltype(output.apar)
    ncanlay = canopy.ncanlay
    ncol = size(output.apar, 2)
    nbnd_sw = size(output.apar_by_band, 1)
    nlev_total = nlay_atm + ncanlay + 1

    flux_up_arr = flux_sw.flux_up
    flux_dn_arr = flux_sw.flux_dn

    # Compute total APAR per canopy layer (from broadband fluxes)
    for gcol in 1:ncol
        for ilay in 1:ncanlay
            @inbounds begin
                flux_in = flux_dn_arr[ilay + 1, gcol] + flux_up_arr[ilay, gcol]
                flux_out = flux_dn_arr[ilay, gcol] + flux_up_arr[ilay + 1, gcol]
                output.apar[ilay, gcol] = flux_in - flux_out
            end
        end
    end

    # Identify PAR bands (400-700 nm ≈ 14286-25000 cm⁻¹)
    par_wn_lo = FT(14286)  # 700 nm
    par_wn_hi = FT(25000)  # 400 nm
    is_par_band = falses(nbnd_sw)
    for ibnd in 1:nbnd_sw
        wn_lo = bnd_lims_wn[1, ibnd]
        wn_hi = bnd_lims_wn[2, ibnd]
        # Band overlaps PAR if its range intersects [par_wn_lo, par_wn_hi]
        if wn_hi ≥ par_wn_lo && wn_lo ≤ par_wn_hi
            is_par_band[ibnd] = true
        end
    end

    # Compute fAPAR from canopy-top incoming PAR
    # (Note: proper per-band fAPAR requires apar_by_band which is filled during the solve loop)
    for gcol in 1:ncol
        total_canopy_apar = FT(0)
        for ilay in 1:ncanlay
            @inbounds total_canopy_apar += output.apar[ilay, gcol]
        end
        canopy_top_dn = flux_dn_arr[ncanlay + 1, gcol]
        @inbounds output.fAPAR[gcol] = canopy_top_dn > 0 ? total_canopy_apar / canopy_top_dn : FT(0)
    end

    # Compute green-APAR from per-band APAR and absorption fractions
    η_car = output.η_car
    for gcol in 1:ncol
        for ilay in 1:ncanlay
            green = FT(0)
            for ibnd in 1:nbnd_sw
                @inbounds begin
                    apar_band = output.apar_by_band[ibnd, ilay, gcol]
                    f_cab = canopy.f_cab_band[ibnd, gcol]
                    f_car = canopy.f_car_band[ibnd, gcol]
                    green += apar_band * (f_cab + η_car * f_car)
                end
            end
            @inbounds output.green_apar[ilay, gcol] = green
        end
    end

    # LW absorption
    if !isnothing(flux_lw)
        flux_up_lw = flux_lw.flux_up
        flux_dn_lw = flux_lw.flux_dn
        for gcol in 1:ncol
            for ilay in 1:ncanlay
                @inbounds begin
                    flux_in = flux_dn_lw[ilay + 1, gcol] + flux_up_lw[ilay, gcol]
                    flux_out = flux_dn_lw[ilay, gcol] + flux_up_lw[ilay + 1, gcol]
                    output.lw_absorbed[ilay, gcol] = flux_in - flux_out
                end
            end
        end
    end

    return nothing
end
