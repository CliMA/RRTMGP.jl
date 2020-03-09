"""
    CloudOptics

Provides cloud optical properties as a function of effective radius for the RRTMGP bands
  Based on Mie calculations for liquid
    and results from doi:10.1175/JAS-D-12-039.1 for ice with variable surface roughness
  Can use either look-up tables or Pade approximates according to which data has been loaded
  Mike Iacono (AER) is the original author
The class can be used as-is but is also intended as an example of how to extend the RTE framework
"""
module CloudOptics

using DocStringExtensions
using ..OpticalProps
using ..Utilities
using ..FortranIntrinsics

export CloudOpticalProps, CloudOpticalPropsPGP
export cloud_optics!
export get_min_radius, get_max_radius

export PadeMethod, LookUpTable
export CloudOpticsLUT
export CloudOpticsPade

abstract type AbstractInterpolationMethod{FT<:AbstractFloat} end

"""
    PadeMethod{FT} <: AbstractInterpolationMethod{FT}

Pade approximation coefficients

# Fields
$(DocStringExtensions.FIELDS)
"""
struct PadeMethod{FT} <: AbstractInterpolationMethod{FT}
    "Extinction coefficients"
    ext::Array{FT}
    "Single scattering albedo"
    ssa::Array{FT}
    "Asymmetry parameter"
    asy::Array{FT}
    "Extinction coefficient particle size regime boundaries"
    sizreg_ext::Array{FT}
    "Single scattering albedo particle size regime boundaries"
    sizreg_ssa::Array{FT}
    "Asymmetry parameter particle size regime boundaries"
    sizreg_asy::Array{FT}
    "particle size lower bound for interpolation"
    rad_lwr::FT
    "particle size upper bound for interpolation"
    rad_upr::FT
    function PadeMethod(
        ext::Array{FT},
        ssa::Array{FT},
        asy::Array{FT},
        sizreg_ext::Array{FT},
        sizreg_ssa::Array{FT},
        sizreg_asy::Array{FT},
    ) where {FT}

    # Error checking
        nbnd = size(ext, 1)
        nsizereg = size(ext, 2)
        ncoeff_ssa_g = size(ssa, 3)
        nbound = length(sizreg_ext)
        nrghice = size(ext, 4)

        @assert nsizereg == 3
        if length(size(ext)) == 4
            @assert all(size(ssa) .== (nbnd, nsizereg, ncoeff_ssa_g, nrghice))
            @assert all(size(asy) .== (nbnd, nsizereg, ncoeff_ssa_g, nrghice))
        else
            @assert all(size(ssa) .== (nbnd, nsizereg, ncoeff_ssa_g))
            @assert all(size(asy) .== (nbnd, nsizereg, ncoeff_ssa_g))
        end
        @assert length(sizreg_ssa) == nbound
        @assert length(sizreg_asy) == nbound
    # Consistency among size regimes
        rad_lwr = sizreg_ext[1]
        rad_upr = sizreg_ext[end]
        @assert !any([sizreg_ssa[1], sizreg_asy[1]] .< rad_lwr)
        @assert !any([sizreg_ssa[nbound], sizreg_asy[nbound]] .> rad_upr)

        return new{FT}(
            ext,
            ssa,
            asy,
            sizreg_ext,
            sizreg_ssa,
            sizreg_asy,
            rad_lwr,
            rad_upr,
        )
    end
end

"""
    LookUpTable{FT,I<:Int} <: AbstractInterpolationMethod{FT}

Lookup table interpolation constants and coefficients

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LookUpTable{FT,I<:Int} <: AbstractInterpolationMethod{FT}
    "particle size lower bound for interpolation"
    rad_lwr::FT
    "particle size upper bound for interpolation"
    rad_upr::FT
    "constant for calculating interpolation indexes"
    rad_fac::FT
    "Number of evenly-spaced elements for linear interpolation"
    nsteps::I
    "Size of step for linear interpolation"
    step_size::FT
    "Extinction coefficients"
    ext::Array{FT}
    "Single scattering albedo"
    ssa::Array{FT}
    "Asymmetry parameter"
    asy::Array{FT}
    function LookUpTable(
        rad_lwr::FT,
        rad_upr::FT,
        rad_fac::FT,
        ext::Array{FT},
        ssa::Array{FT},
        asy::Array{FT},
    ) where {FT<:AbstractFloat}

    # Error checking
    # LUT coefficient dimensions
        nsize = size(ext, 1)
        nbnd = size(ext, 2)

    # Error checking
    #   Can we check for consistency between table bounds and _fac?
        if length(size(ext)) == 3
            nrghice = size(ext, 3)
            @assert all(size(ssa) .== (nsize, nbnd, nrghice))
            @assert all(size(asy) .== (nsize, nbnd, nrghice))
        else
            @assert all(size(ssa) .== (nsize, nbnd))
            @assert all(size(asy) .== (nsize, nbnd))
        end

        nsteps = size(ext, 1)
        step_size = (rad_upr - rad_lwr) / FT(nsteps - 1)

        return new{FT,Int}(
            rad_lwr,
            rad_upr,
            rad_fac,
            nsteps,
            step_size,
            ext,
            ssa,
            asy,
        )
    end
end
get_min_radius(aim::AbstractInterpolationMethod) = aim.rad_lwr
get_max_radius(aim::AbstractInterpolationMethod) = aim.rad_upr

get_num_roughness_types(this::PadeMethod) = size(this.ext, 4)
get_num_roughness_types(this::LookUpTable) = size(this.ext, 3)

"""
    CloudOpticsLUT{FT, I} <: AbstractOpticalProps{FT, I}

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudOpticsLUT{FT,I} <: AbstractOpticalProps{FT,I}
    base::OpticalPropsBase{FT,I}
    "Ice surface roughness category - needed for Yang (2013) ice optics parameterization"
    icergh::I
    "Interpolation method for liquid properties"
    liq::LookUpTable{FT,I}
    "Interpolation method for ice properties"
    ice::LookUpTable{FT,I}
end

"""
    CloudOpticsPade{FT, I} <: AbstractOpticalProps{FT, I}

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudOpticsPade{FT,I} <: AbstractOpticalProps{FT,I}
    "Base optical properties"
    base::OpticalPropsBase{FT,I}
    "Ice surface roughness category - needed for Yang (2013) ice optics parameterization"
    icergh::I
    "Interpolation method for liquid properties"
    liq::PadeMethod{FT}
    "Interpolation method for ice properties"
    ice::PadeMethod{FT}
end

"""
    CloudOpticalProps{FT<:AbstractFloat}

Cloud properties (liquid or ice) used to combine with
averaged optical properties

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudOpticalProps{FT<:AbstractFloat}
    "Water path"
    wp::Array{FT}
    "Effective radius"
    re::Array{FT}
end
CloudOpticalProps(FT, ncol, nlay) =
    CloudOpticalProps{FT}(zeros(FT, ncol, nlay), zeros(FT, ncol, nlay))

"""
    CloudOpticalPropsPGP{FT<:AbstractFloat}

Cloud properties (liquid or ice) used to combine with
averaged optical properties per grid point

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct CloudOpticalPropsPGP{FT<:AbstractFloat}
    "Water path"
    wp::FT
    "Effective radius"
    re::FT
end
CloudOpticalPropsPGP(FT) = CloudOpticalPropsPGP{FT}(FT(0), FT(0))

function Base.convert(
    ::Type{CloudOpticalProps},
    data::Array{CloudOpticalPropsPGP{FT}},
) where {FT}
    s = size(data)
    return CloudOpticalProps{FT}(
        Array{FT}([data[i, j].wp for i in 1:s[1], j in 1:s[2]]),
        Array{FT}([data[i, j].re for i in 1:s[1], j in 1:s[2]]),
    )
end

Base.convert(
    ::Type{Array{CloudOpticalPropsPGP}},
    data::CloudOpticalProps{FT},
) where {FT} =
    [CloudOpticalPropsPGP{FT}(data.wp[i, j], data.re[i, j]) for i in 1:size(
        data.wp,
        1,
    ), j in 1:size(data.wp, 2)]


#####
##### Derive cloud optical properties from provided cloud physical properties
#####

"""
    combine_optical_props!(op::AbstractOpticalPropsArry, liq::TwoStream{FT}, ice::TwoStream{FT}) where {FT<:AbstractFloat}

Combine liquid and ice contributions into total cloud optical properties
   See also the `increment!` routines in `OpticalProps_kernels`
"""
function combine_optical_props!(
    op::OneScalar{FT},
    liq::TwoStream{FT},
    ice::TwoStream{FT},
) where {FT<:AbstractFloat}
  # Absorption optical depth  = (1-ssa) * τ = τ - τssa
    nbnd = size(liq.τ, 3)
    op.τ[:, :, 1:nbnd] .= (liq.τ .- liq.ssa) + (ice.τ .- ice.ssa)
end
function combine_optical_props!(
    op::TwoStream{FT},
    liq::TwoStream{FT},
    ice::TwoStream{FT},
) where {FT<:AbstractFloat}
    nbnd = size(liq.τ, 3)
    τ = liq.τ + ice.τ
    ssa = liq.ssa + ice.ssa
    op.g[:, :, 1:nbnd] .= (liq.g + ice.g) ./ max.(eps(FT), ssa)
    op.ssa[:, :, 1:nbnd] .= ssa ./ max.(eps(FT), τ)
    op.τ[:, :, 1:nbnd] .= τ
end
function combine_optical_props!(
    op::OneScalarPGP{FT},
    liq::TwoStreamPGP{FT},
    ice::TwoStreamPGP{FT},
) where {FT<:AbstractFloat}
  # Absorption optical depth  = (1-ssa) * τ = τ - τssa
    nbnd = size(liq.τ, 1)
    op.τ[1:nbnd] .= (liq.τ .- liq.ssa) + (ice.τ .- ice.ssa)
end
function combine_optical_props!(
    op::TwoStreamPGP{FT},
    liq::TwoStreamPGP{FT},
    ice::TwoStreamPGP{FT},
) where {FT<:AbstractFloat}
    nbnd = size(liq.τ, 1)
    τ = liq.τ + ice.τ
    ssa = liq.ssa + ice.ssa
    op.g[1:nbnd] .= (liq.g + ice.g) ./ max.(eps(FT), ssa)
    op.ssa[1:nbnd] .= ssa ./ max.(eps(FT), τ)
    op.τ[1:nbnd] .= τ
end

function validate_cloud_optics!(
    this::AbstractOpticalProps{FT},
    liq::CloudOpticalProps{FT},
    ice::CloudOpticalProps{FT},
    optical_props::AbstractOpticalPropsArry{FT},
) where {FT}
    liqmsk = BitArray(liq.wp .> FT(0))
    icemsk = BitArray(ice.wp .> FT(0))
  # Error checking
    @assert bands_are_equal(this, optical_props)
    @assert get_nband(optical_props) == get_ngpt(optical_props)
    @assert !(
        this.icergh < 1 || this.icergh > get_num_roughness_types(this.ice)
    )
    @assert !any_vals_outside(
        liq.re,
        liqmsk,
        this.liq.rad_lwr,
        this.liq.rad_upr,
    )
    @assert !any_vals_outside(
        ice.re,
        icemsk,
        this.ice.rad_lwr,
        this.ice.rad_upr,
    )
    any(liqmsk) && @assert !any_vals_less_than(liq.wp, liqmsk, FT(0))
    any(icemsk) && @assert !any_vals_less_than(ice.wp, icemsk, FT(0))
    return nothing
end
function validate_cloud_optics!(
    this::AbstractOpticalProps{FT},
    liq::CloudOpticalPropsPGP{FT},
    ice::CloudOpticalPropsPGP{FT},
    optical_props::AbstractOpticalPropsPGP{FT},
) where {FT}
    liqmsk = liq.wp > FT(0)
    icemsk = ice.wp > FT(0)
  # Error checking
    @assert bands_are_equal(this, optical_props)
    @assert get_nband(optical_props) == get_ngpt(optical_props)
    @assert !(
        this.icergh < 1 || this.icergh > get_num_roughness_types(this.ice)
    )
    @assert !any_vals_outside(
        liq.re,
        liqmsk,
        this.liq.rad_lwr,
        this.liq.rad_upr,
    )
    @assert !any_vals_outside(
        ice.re,
        icemsk,
        this.ice.rad_lwr,
        this.ice.rad_upr,
    )
    liqmsk && @assert !any_vals_less_than(liq.wp, liqmsk, FT(0))
    icemsk && @assert !any_vals_less_than(ice.wp, icemsk, FT(0))
    return nothing
end

"""
    cloud_optics!(this::AbstractOpticalProps{FT},
                  clouds_liq::CloudOpticalProps{FT},
                  clouds_ice::CloudOpticalProps{FT},
                  optical_props::AbstractOpticalPropsArry{FT}) where {FT<:AbstractFloat}

Compute single-scattering properties in

 - `optical_props` optical properties, see [`AbstractOpticalPropsArry`](@ref)

given

 - `this` cloud optics, see [`AbstractOpticalProps`](@ref)
 - `clouds_liq` cloud properties for liquid, see [`CloudOpticalProps`](@ref)
 - `clouds_ice` cloud properties for ice   , see [`CloudOpticalProps`](@ref)
"""
function cloud_optics!(
    this::CloudOpticsPade{FT},
    clouds_liq::CloudOpticalProps{FT},
    clouds_ice::CloudOpticalProps{FT},
    optical_props::AbstractOpticalPropsArry{FT},
) where {FT<:AbstractFloat}

    nbnd = get_nband(this)

    validate_cloud_optics!(this, clouds_liq, clouds_ice, optical_props)

    liq = TwoStream(FT, size(clouds_liq.wp)..., nbnd)
    ice = TwoStream(FT, size(clouds_ice.wp)..., nbnd)

  #### Compute cloud optical properties.

    compute_all_from_pade!(nbnd, clouds_liq, this.liq, liq)              # Liquid
    compute_all_from_pade!(nbnd, clouds_ice, this.ice, ice, this.icergh) # Ice

  # Copy total cloud properties onto outputs
    combine_optical_props!(optical_props, liq, ice)

end
function cloud_optics!(
    this::CloudOpticsLUT{FT},
    clouds_liq::CloudOpticalProps{FT},
    clouds_ice::CloudOpticalProps{FT},
    optical_props::AbstractOpticalPropsArry{FT},
) where {FT}

    nbnd = get_nband(this)

    validate_cloud_optics!(this, clouds_liq, clouds_ice, optical_props)

    liq = TwoStream(FT, size(clouds_liq.wp)..., nbnd)
    ice = TwoStream(FT, size(clouds_ice.wp)..., nbnd)

  #### Compute cloud optical properties.

    compute_all_from_table!(nbnd, clouds_liq, this.liq, liq)              # Liquid
    compute_all_from_table!(nbnd, clouds_ice, this.ice, ice, this.icergh) # Ice

  # Copy total cloud properties onto outputs
    combine_optical_props!(optical_props, liq, ice)

end
function cloud_optics!(
    this::CloudOpticsPade{FT},
    clouds_liq::CloudOpticalPropsPGP{FT},
    clouds_ice::CloudOpticalPropsPGP{FT},
    optical_props::AbstractOpticalPropsPGP{FT},
) where {FT<:AbstractFloat}

    nbnd = get_nband(this)

    validate_cloud_optics!(this, clouds_liq, clouds_ice, optical_props)

    liq = TwoStreamPGP(FT, nbnd)
    ice = TwoStreamPGP(FT, nbnd)

  #### Compute cloud optical properties.

    compute_all_from_pade!(nbnd, clouds_liq, this.liq, liq)              # Liquid
    compute_all_from_pade!(nbnd, clouds_ice, this.ice, ice, this.icergh) # Ice

  # Copy total cloud properties onto outputs
    combine_optical_props!(optical_props, liq, ice)

end
function cloud_optics!(
    this::CloudOpticsLUT{FT},
    clouds_liq::CloudOpticalPropsPGP{FT},
    clouds_ice::CloudOpticalPropsPGP{FT},
    optical_props::AbstractOpticalPropsPGP{FT},
) where {FT}

    nbnd = get_nband(this)

    validate_cloud_optics!(this, clouds_liq, clouds_ice, optical_props)

    liq = TwoStreamPGP(FT, nbnd)
    ice = TwoStreamPGP(FT, nbnd)

  #### Compute cloud optical properties.

    compute_all_from_table!(nbnd, clouds_liq, this.liq, liq)              # Liquid
    compute_all_from_table!(nbnd, clouds_ice, this.ice, ice, this.icergh) # Ice

  # Copy total cloud properties onto outputs
    combine_optical_props!(optical_props, liq, ice)

end

"""
    compute_all_from_table!(nbnd::I,
                            clouds::CloudOpticalProps{FT},
                            lut::LookUpTable{FT,I},
                            op::TwoStream{FT,I},
                            icergh::Union{Nothing,I}=nothing) where {FT<:AbstractFloat,I<:Int}

Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
elements starting at "offset." in

 - `op` optical properties, see [`TwoStream`](@ref)

given

 - `nbnd` number of bands
 - `clouds` cloud optical properties, see [`CloudOpticalProps`](@ref)
 - `lut` interpolation method (look-up table), see [`LookUpTable`](@ref)
 - `icergh` ice surface roughness category - needed for Yang (2013) ice optics parameterization

The table's second dimension is band. Returns 0 where `clouds.wp ≤ 0`.
"""
function compute_all_from_table!(
    nbnd::I,
    clouds::CloudOpticalProps{FT},
    lut::LookUpTable{FT,I},
    op::TwoStream{FT,I},
    icergh::Union{Nothing,I} = nothing,
) where {FT<:AbstractFloat,I<:Int}

    ncol, nlay = size(clouds.re)

    offset = lut.rad_lwr

    τ_table = icergh ≠ nothing ? lut.ext[:, :, icergh] : lut.ext
    ssa_table = icergh ≠ nothing ? lut.ssa[:, :, icergh] : lut.ssa
    asy_table = icergh ≠ nothing ? lut.asy[:, :, icergh] : lut.asy

    @inbounds for ibnd in 1:nbnd
        @inbounds for ilay in 1:nlay
            @inbounds for icol in 1:ncol
                if clouds.wp[icol, ilay] > FT(0)
                    index = convert(
                        Int,
                        min(
                            floor((clouds.re[icol, ilay] - offset) / lut.step_size) +
                            1,
                            lut.nsteps - 1,
                        ),
                    )
                    fint = (clouds.re[icol, ilay] - offset) / lut.step_size - (
                        index - 1
                    )
                    τ = clouds.wp[icol, ilay] * (
                        τ_table[index, ibnd] +
                        fint * (τ_table[index+1, ibnd] - τ_table[index, ibnd])
                    )
                    τs = τ * (
                        ssa_table[index, ibnd] +
                        fint * (
                            ssa_table[index+1, ibnd] - ssa_table[index, ibnd]
                        )
                    )
                    op.g[icol, ilay, ibnd] = τs * (
                        asy_table[index, ibnd] +
                        fint * (
                            asy_table[index+1, ibnd] - asy_table[index, ibnd]
                        )
                    )
                    op.ssa[icol, ilay, ibnd] = τs
                    op.τ[icol, ilay, ibnd] = τ
                else
                    op.τ[icol, ilay, ibnd] = FT(0)
                    op.ssa[icol, ilay, ibnd] = FT(0)
                    op.g[icol, ilay, ibnd] = FT(0)
                end
            end
        end
    end
end
function compute_all_from_table!(
    nbnd::I,
    clouds::CloudOpticalPropsPGP{FT},
    lut::LookUpTable{FT,I},
    op::TwoStreamPGP{FT,I},
    icergh::Union{Nothing,I} = nothing,
) where {FT<:AbstractFloat,I<:Int}

    offset = lut.rad_lwr

    τ_table = icergh ≠ nothing ? lut.ext[:, :, icergh] : lut.ext
    ssa_table = icergh ≠ nothing ? lut.ssa[:, :, icergh] : lut.ssa
    asy_table = icergh ≠ nothing ? lut.asy[:, :, icergh] : lut.asy

    @inbounds for ibnd in 1:nbnd
        if clouds.wp > FT(0)
            index = convert(
                Int,
                min(
                    floor((clouds.re - offset) / lut.step_size) + 1,
                    lut.nsteps - 1,
                ),
            )
            fint = (clouds.re - offset) / lut.step_size - (index - 1)
            τ = clouds.wp * (
                τ_table[index, ibnd] +
                fint * (τ_table[index+1, ibnd] - τ_table[index, ibnd])
            )
            τs = τ * (
                ssa_table[index, ibnd] +
                fint * (ssa_table[index+1, ibnd] - ssa_table[index, ibnd])
            )
            op.g[ibnd] = τs * (
                asy_table[index, ibnd] +
                fint * (asy_table[index+1, ibnd] - asy_table[index, ibnd])
            )
            op.ssa[ibnd] = τs
            op.τ[ibnd] = τ
        else
            op.τ[ibnd] = FT(0)
            op.ssa[ibnd] = FT(0)
            op.g[ibnd] = FT(0)
        end
    end
end

#####
##### Pade functions
#####

"""
    compute_all_from_pade!(nbnd::I,
                           clouds::CloudOpticalProps{FT},
                           pm::PadeMethod{FT},
                           op::TwoStream{FT,I},
                           icergh::Union{Nothing,I}=nothing) where {FT<:AbstractFloat,I<:Int}

 - `nbnd` number of bands
 - `clouds` cloud optical properties, see [`CloudOpticalProps`](@ref)
 - `pm` interpolation method (Pade), see [`PadeMethod`](@ref)
 - `op` optical properties, see [`TwoStream`](@ref)
 - `icergh` ice surface roughness category - needed for Yang (2013) ice optics parameterization
"""
function compute_all_from_pade!(
    nbnd::I,
    clouds::CloudOpticalProps{FT},
    pm::PadeMethod{FT},
    op::TwoStream{FT,I},
    icergh::Union{Nothing,I} = nothing,
) where {FT<:AbstractFloat,I<:Int}

    nsizes = size(pm.ext, 2)
  # Cloud optical properties from Pade coefficient method
  #   Hard coded assumptions: order of approximants, three size regimes
    m_ext, n_ext = 2, 3
    m_ssa, n_ssa = 2, 2
    m_asy, n_asy = 2, 2

    coeffs_ext = icergh ≠ nothing ? pm.ext[:, :, :, icergh] : pm.ext
    coeffs_ssa = icergh ≠ nothing ? pm.ssa[:, :, :, icergh] : pm.ssa
    coeffs_asy = icergh ≠ nothing ? pm.asy[:, :, :, icergh] : pm.asy

    ncol, nlay = size(clouds.re)
    @inbounds for ibnd in 1:nbnd
        @inbounds for ilay in 1:nlay
            @inbounds for icol in 1:ncol
                if clouds.wp[icol, ilay] > FT(0)
                    re = clouds.re[icol, ilay]
          #
          # Finds index into size regime table
          # This works only if there are precisely three size regimes (four bounds) and it's
          #   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
          #
                    irad = convert(
                        Int,
                        min(
                            floor((re - pm.sizreg_ext[2]) / pm.sizreg_ext[3]) +
                            2,
                            3,
                        ),
                    )
                    τ = clouds.wp[icol, ilay] * pade_eval(
                        ibnd,
                        nbnd,
                        nsizes,
                        m_ext,
                        n_ext,
                        irad,
                        re,
                        coeffs_ext,
                    )

                    irad = convert(
                        Int,
                        min(
                            floor((re - pm.sizreg_ssa[2]) / pm.sizreg_ssa[3]) +
                            2,
                            3,
                        ),
                    )
          # Pade approximants for co-albedo can sometimes be negative
                    τs = τ * (
                        FT(1) - max(
                            FT(0),
                            pade_eval(
                                ibnd,
                                nbnd,
                                nsizes,
                                m_ssa,
                                n_ssa,
                                irad,
                                re,
                                coeffs_ssa,
                            ),
                        )
                    )
                    irad = convert(
                        Int,
                        min(
                            floor((re - pm.sizreg_asy[2]) / pm.sizreg_asy[3]) +
                            2,
                            3,
                        ),
                    )
                    op.g[icol, ilay, ibnd] = τs * pade_eval(
                        ibnd,
                        nbnd,
                        nsizes,
                        m_asy,
                        n_asy,
                        irad,
                        re,
                        coeffs_asy,
                    )
                    op.ssa[icol, ilay, ibnd] = τs
                    op.τ[icol, ilay, ibnd] = τ
                else
                    op.τ[icol, ilay, ibnd] = FT(0)
                    op.ssa[icol, ilay, ibnd] = FT(0)
                    op.g[icol, ilay, ibnd] = FT(0)
                end
            end
        end
    end
end
function compute_all_from_pade!(
    nbnd::I,
    clouds::CloudOpticalPropsPGP{FT},
    pm::PadeMethod{FT},
    op::TwoStreamPGP{FT,I},
    icergh::Union{Nothing,I} = nothing,
) where {FT<:AbstractFloat,I<:Int}

    nsizes = size(pm.ext, 2)
  # Cloud optical properties from Pade coefficient method
  #   Hard coded assumptions: order of approximants, three size regimes
    m_ext, n_ext = 2, 3
    m_ssa, n_ssa = 2, 2
    m_asy, n_asy = 2, 2

    coeffs_ext = icergh ≠ nothing ? pm.ext[:, :, :, icergh] : pm.ext
    coeffs_ssa = icergh ≠ nothing ? pm.ssa[:, :, :, icergh] : pm.ssa
    coeffs_asy = icergh ≠ nothing ? pm.asy[:, :, :, icergh] : pm.asy

    @inbounds for ibnd in 1:nbnd
        if clouds.wp > FT(0)
            re = clouds.re
      #
      # Finds index into size regime table
      # This works only if there are precisely three size regimes (four bounds) and it's
      #   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
      #
            irad = convert(
                Int,
                min(floor((re - pm.sizreg_ext[2]) / pm.sizreg_ext[3]) + 2, 3),
            )
            τ = clouds.wp * pade_eval(
                ibnd,
                nbnd,
                nsizes,
                m_ext,
                n_ext,
                irad,
                re,
                coeffs_ext,
            )

            irad = convert(
                Int,
                min(floor((re - pm.sizreg_ssa[2]) / pm.sizreg_ssa[3]) + 2, 3),
            )
      # Pade approximants for co-albedo can sometimes be negative
            τs = τ * (
                FT(1) - max(
                    FT(0),
                    pade_eval(
                        ibnd,
                        nbnd,
                        nsizes,
                        m_ssa,
                        n_ssa,
                        irad,
                        re,
                        coeffs_ssa,
                    ),
                )
            )
            irad = convert(
                Int,
                min(floor((re - pm.sizreg_asy[2]) / pm.sizreg_asy[3]) + 2, 3),
            )
            op.g[ibnd] = τs * pade_eval(
                ibnd,
                nbnd,
                nsizes,
                m_asy,
                n_asy,
                irad,
                re,
                coeffs_asy,
            )
            op.ssa[ibnd] = τs
            op.τ[ibnd] = τ
        else
            op.τ[ibnd] = FT(0)
            op.ssa[ibnd] = FT(0)
            op.g[ibnd] = FT(0)
        end
    end
end

#####
##### Ancillary functions
#####

"""
    pade_eval(iband::I,
              nbnd::I,
              nrads::I,
              m::I,
              n::I,
              irad::I,
              re::FT,
              pade_coeffs::AbstractArray{FT}) where {FT<:AbstractFloat, I<:Int}

Evaluate Pade approximant of order [m/n], given

 - `pade_coeffs` pade coefficients [nbnd, nrads, 1:m+n+1]
"""
function pade_eval(
    iband::I,
    nbnd::I,
    nrads::I,
    m::I,
    n::I,
    irad::I,
    re::FT,
    pade_coeffs::AbstractArray{FT},
) where {FT<:AbstractFloat,I<:Int}

    denom = pade_coeffs[iband, irad, n + m + 1]
    @inbounds for i in n-1+m:-1:1+m
        denom = pade_coeffs[iband, irad, i+1] + re * denom
    end
    denom = 1 + re * denom

    numer = pade_coeffs[iband, irad, m+1]
    @inbounds for i in m-1:-1:1
        numer = pade_coeffs[iband, irad, i+1] + re * numer
    end
    numer = pade_coeffs[iband, irad, 1] + re * numer

    pade_eval_1 = numer / denom
    return pade_eval_1
end

end
