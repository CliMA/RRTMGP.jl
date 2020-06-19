module GrayAtmos
using CLIMAParameters
using ..AtmosphericStates
using RRTMGP.MeshOrientations
using RRTMGP.RadiativeBoundaryConditions
using ..SourceFunctions
using RRTMGP.OpticalProps
using DocStringExtensions

export gen_grid_gray_atmos,
    source_func_longwave_gray_atmos,
    gas_optics_gray_atmos!,
    compute_gas_τs_gray_atmos!,
    gen_grid_gray_atmos_old,
    GrayAtmosphericState,
    tlay_to_tlev!,
    GrayRRTMGP

"""
    GrayAtmosphericState{FT}

Atmospheric conditions, used to compute optical properties with the gray atmosphere approximation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayAtmosphericState{FT,I} <: AbstractAtmosphericState{FT,I}
    "Layer pressures `[Pa, mb]`; (`ncol,nlay`)"
    p_lay::Array{FT,2}
    "Level pressures `[Pa, mb]`; (`ncol,nlay+1`)"
    p_lev::Array{FT,2}
    "Layer temperatures `[K]`; (`ncol,nlay`)"
    t_lay::Array{FT,2}
    "Level temperatures `[K]`; (`ncol,nlay+1`)"
    t_lev::Array{FT,2}
    "Surface temperatures `[K]`; (`ncol`)"
    t_sfc::Vector{FT}
    "lapse rate"
    α::FT
    "optical thickness parameter"
    d0::FT
    "Mesh orientation, see [`MeshOrientation`](@ref)"
    mesh_orientation::MeshOrientation{I}
    "Number of layers."
    nlay::I
    "Number of columns."
    ncol::I

    function GrayAtmosphericState(
        nlay::Int,
        lat::FT,
        p0::FT,
        pe::FT;
        step = "linear",
    ) where {FT<:AbstractFloat}
        nlev = Int(nlay + 1)
        ncol = Int(1)
        p_lay = Array{FT}(undef, 1, nlay)
        p_lev = Array{FT}(undef, 1, nlev)
        t_lay = Array{FT}(undef, 1, nlay)
        t_lev = Array{FT}(undef, 1, nlev)

        efac = log(p0 / pe) / nlay # multiplication factor for each pressure step
        Δp = (p0 - pe) / nlay     # linear pressure distribution
        te = 300               # global mean surface temperature (K)
        tt = 200               # skin temp at top of atmosphere (K)
        Δt = 60
        ts = te + Δt * (1.0 / 3.0 - sin(lat)^2) # surface temp at a given latitude (K)
        d0 = FT((ts / tt)^4.0 - 1.0) # optical depth
        α = FT(3.5)           # lapse rate of radiative equillibrium

        p_lev[1, nlev] = p0
        t_lev[1, nlev] = tt * (1 + d0 * (p_lev[1, nlev] / p0)^α)^0.25

        for i = nlay:-1:1
            #p_lev[1,i] = p_lev[1,i+1] * exp(-efac)
            if step == "linear"
                p_lev[1, i] = p_lev[1, i+1] - Δp
            else
                p_lev[1, i] = p_lev[1, i+1] * exp(-efac)
            end

            p_lay[1, i] = (p_lev[1, i] + p_lev[1, i+1]) * 0.5

            t_lev[1, i] = tt * (1 + d0 * (p_lev[1, i] / p0)^α)^0.25
            t_lay[1, i] = tt * (1 + d0 * (p_lay[1, i] / p0)^α)^0.25
        end

        t_sfc = [t_lev[1, end]]
        top_at_1 = p_lay[1, 1] < p_lay[1, nlay]

        return new{FT,Int}(
            p_lay,
            p_lev,
            t_lay,
            t_lev,
            t_sfc,
            α,
            d0,
            MeshOrientation(top_at_1, nlay),
            nlay,
            ncol,
        )
    end
end

function source_func_longwave_gray_atmos(
    ::Type{FT},
    ncol::I,
    nlay::I,
    ngpt::I,
) where {FT<:AbstractFloat,I<:Int}
    lay_source = zeros(FT, ncol, nlay, ngpt)
    lev_source_inc = zeros(FT, ncol, nlay, ngpt)
    lev_source_dec = zeros(FT, ncol, nlay, ngpt)
    sfc_source = zeros(FT, ncol, ngpt)
    p_frac = zeros(FT, ncol, nlay, ngpt)
    return SourceFuncLongWave{FT,I}(
        nothing,
        lay_source,
        lev_source_inc,
        lev_source_dec,
        sfc_source,
        p_frac,
    )
end

function gas_optics_gray_atmos!(
    as::GrayAtmosphericState{FT,I},
    optical_props::AbstractOpticalPropsArry{FT,I},
    source::SourceFuncLongWave{FT,I},
) where {FT<:AbstractFloat,I<:Int}

    compute_gas_τs_gray_atmos!(as, optical_props)
    compute_Planck_source!(source, as)
    return nothing
end

function compute_gas_τs_gray_atmos!(
    as::GrayAtmosphericState{FT,I},
    optical_props::AbstractOpticalPropsArry{FT,I},
) where {FT<:AbstractFloat,I<:Int}
    nlay = as.nlay
    α = as.α
    d0 = as.d0
    top_at_1 = as.mesh_orientation.top_at_1
    p_lay = as.p_lay
    p_lev = as.p_lev

    p0 = (top_at_1 ? as.p_lev[1, end] : as.p_lev[1, 1])
    for i = 1:nlay
        optical_props.τ[1, i, 1] =
            (α * d0 * (p_lay[1, i] ./ p0) .^ α ./ p_lay[1, i]) *
            (p_lev[1, i+1] - p_lev[1, i])
    end
    return nothing
end

function compute_Planck_source!(
    source::SourceFuncLongWave{FT,I},
    as::GrayAtmosphericState{FT,I},
) where {FT<:AbstractFloat,I<:Int}
    ncol = as.ncol
    nlay = as.nlay
    for col = 1:ncol
        source.sfc_source[col, 1] = Stefan() * as.t_sfc[col]^4.0 / π # computing sfc_source
        for lay = 1:nlay
            source.lay_source[col, lay, 1] =
                Stefan() * as.t_lay[col, lay]^4.0 / π # computing lay_source

            source.lev_source_inc[col, lay, 1] =
                Stefan() * as.t_lev[col, lay+1]^4.0 / π
            source.lev_source_dec[col, lay, 1] =
                Stefan() * as.t_lev[col, lay]^4.0 / π
        end
    end
    return nothing
end

function tlay_to_tlev!(p_lay, p_lev, t_lay, t_lev)
    ncol, nlay, nlev = size(t_lay, 1), size(t_lay, 2), size(t_lev, 2)
    #    for ilay in nlay:-1:2
    #        t_lev[1, ilay] = 2 * t_lay[1, ilay] - t_lev[1, ilay+1]
    #    end
    #    t_lev[1, 1] = 2 * t_lay[1, 1] - t_lev[1, 2]

    for j = 2:nlay-1
        t_lev[1, j] =
            (1.0 / 3.0) * t_lay[1, j-1] + (5.0 / 6.0) * t_lay[1, j] -
            (1.0 / 6.0) * t_lay[1, j+1]
    end
    t_lev[1, 1] = 2 * t_lay[1, 1] - t_lev[1, 2]
    t_lev[1, nlay] =
        (1.0 / 3.0) * t_lay[1, nlay] + (5.0 / 6.0) * t_lay[1, nlay-1] -
        (1.0 / 6.0) * t_lay[1, nlay-2]
    t_lev[1, nlay+1] = 2 * t_lay[1, nlay] - t_lev[1, nlay]

end

struct GrayRRTMGP{FT,I}
    as::GrayAtmosphericState{FT,I}
    op::AbstractOpticalPropsArry{FT,I}
    src::AbstractSourceFunc{FT,I}
    bcs::AbstractRadiativeBoundaryConditions{FT}

    function GrayRRTMGP(
        nlay::Int,
        lat::FT,
        p0::FT,
        pe::FT,
        optical_props_constructor,
        sfc_emis::Array{FT,2},
        stype,
    ) where {FT<:AbstractFloat}
        ncol, ngpt, nbnd = 1, 1, 1
        if stype == "lw"
            return new{FT,Int}(
                GrayAtmosphericState(nlay, lat, p0, pe),
                optical_props_constructor(FT, ncol, nlay, ngpt),
                source_func_longwave_gray_atmos(FT, ncol, nlay, ngpt),
                LongwaveBCs(sfc_emis),
            )
        else
            return nothing
        end
    end
end

end
