
include("OpticsUtils.jl")
include("GasOptics.jl")
include("CloudOptics.jl")
include("GrayOpticsKernels.jl")

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glay,
    gcol,
    sf::AbstractSourceLW{FT},
    igpt,
    lkp::LookUpLW{FT},
    lkp_cld::LookUpCld,
) where {FT <: AbstractFloat}
    ibnd = lkp.major_gpt2bnd[igpt]

    compute_optical_props_kernel!(op, as, glay, gcol, sf, igpt, lkp) # longwave gas optics
    add_cloud_optics_2stream(op, as, lkp, lkp_cld, glay, gcol, ibnd, igpt) # add cloud optics
    return nothing
end

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glay,
    gcol,
    sf::AbstractSourceLW{FT},
    igpt,
    lkp::LookUpLW{FT},
) where {FT <: AbstractFloat}
    t_sfc = as.t_sfc[gcol]
    t_lev = as.t_lev
    compute_optical_props_kernel!(op, as, glay, gcol, igpt, lkp, sf, t_lev, t_sfc)
    return nothing
end

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glay,
    gcol,
    igpt,
    lkp::LookUpSW{FT},
    lkp_cld::LookUpCld,
) where {FT <: AbstractFloat}
    ibnd = lkp.major_gpt2bnd[igpt]

    compute_optical_props_kernel!(op, as, glay, gcol, igpt, lkp) # shortwave gas optics
    add_cloud_optics_2stream(op, as, lkp, lkp_cld, glay, gcol, ibnd, igpt) # add cloud optics
    return nothing
end

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glay,
    gcol,
    igpt::Int,
    lkp::LookUpLW{FT},
    src_args...,
) where {FT <: AbstractFloat}

    vmr = as.vmr
    col_dry = as.col_dry[glay, gcol]
    p_lay = as.p_lay[glay, gcol]
    t_lay = as.t_lay[glay, gcol]
    ibnd = lkp.major_gpt2bnd[igpt]
    τ, ssa = compute_τ_ssa_lw_src!(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol, src_args...)

    @inbounds op.τ[glay, gcol] = τ

    if op isa TwoStream
        @inbounds op.ssa[glay, gcol] = ssa
        @inbounds op.g[glay, gcol] = FT(0) # initializing asymmetry parameter
    end
    return nothing
end

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glay,
    gcol,
    igpt::Int,
    lkp::LookUpSW{FT},
) where {FT <: AbstractFloat}

    vmr = as.vmr
    col_dry = as.col_dry[glay, gcol]
    p_lay = as.p_lay[glay, gcol]
    t_lay = as.t_lay[glay, gcol]
    ibnd = lkp.major_gpt2bnd[igpt]
    τ, ssa = compute_τ_ssa_lw_src!(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)

    @inbounds op.τ[glay, gcol] = τ

    if op isa TwoStream
        @inbounds op.ssa[glay, gcol] = ssa
        @inbounds op.g[glay, gcol] = FT(0) # initializing asymmetry parameter
    end
    return nothing
end
