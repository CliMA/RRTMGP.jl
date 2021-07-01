
include("OpticsUtils.jl")
include("GasOptics.jl")
include("CloudOptics.jl")
include("GrayOpticsKernels.jl")

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glaycol,
    sf::AbstractSourceLW{FT},
    igpt,
    lkp::LookUpLW{I,FT},
    lkp_cld::LookUpCld,
) where {I<:Int,FT<:AbstractFloat}
    ibnd = lkp.major_gpt2bnd[igpt]

    compute_optical_props_kernel!(op, as, glaycol, sf, igpt, lkp) # longwave gas optics
    add_cloud_optics_2stream(op, as, lkp, lkp_cld, glaycol, ibnd, igpt) # add cloud optics
    return nothing
end

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glaycol,
    sf::AbstractSourceLW{FT},
    igpt,
    lkp::LookUpLW{I,FT},
) where {I<:Int,FT<:AbstractFloat}
    t_sfc = as.t_sfc[glaycol[2]]
    t_lev = as.t_lev
    compute_optical_props_kernel!(op, as, glaycol, igpt, lkp, sf, t_lev, t_sfc)
    return nothing
end

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glaycol,
    igpt,
    lkp::LookUpSW{I,FT},
    lkp_cld::LookUpCld,
) where {I<:Int,FT<:AbstractFloat}
    ibnd = lkp.major_gpt2bnd[igpt]

    compute_optical_props_kernel!(op, as, glaycol, igpt, lkp) # shortwave gas optics
    add_cloud_optics_2stream(op, as, lkp, lkp_cld, glaycol, ibnd, igpt) # add cloud optics
    return nothing
end

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glaycol::Tuple{I,I},
    igpt::I,
    lkp::LookUpLW{I,FT},
    src_args...,
) where {I<:Int,FT<:AbstractFloat}

    vmr = as.vmr
    col_dry = as.col_dry[glaycol...]
    p_lay = as.p_lay[glaycol...]
    t_lay = as.t_lay[glaycol...]
    ibnd = lkp.major_gpt2bnd[igpt]
    τ, ssa = compute_τ_ssa_lw_src(
        lkp,
        vmr,
        col_dry,
        igpt,
        ibnd,
        p_lay,
        t_lay,
        glaycol,
        src_args...,
    )

    @inbounds op.τ[glaycol...] = τ

    if op isa TwoStream
        @inbounds op.ssa[glaycol...] = ssa
        @inbounds op.g[glaycol...] = FT(0) # initializing asymmetry parameter
    end
    return nothing
end

function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    glaycol::Tuple{I,I},
    igpt::I,
    lkp::LookUpSW{I,FT},
) where {I<:Int,FT<:AbstractFloat}

    vmr = as.vmr
    col_dry = as.col_dry[glaycol...]
    p_lay = as.p_lay[glaycol...]
    t_lay = as.t_lay[glaycol...]
    ibnd = lkp.major_gpt2bnd[igpt]
    τ, ssa = compute_τ_ssa_lw_src(
        lkp,
        vmr,
        col_dry,
        igpt,
        ibnd,
        p_lay,
        t_lay,
        glaycol,
    )

    @inbounds op.τ[glaycol...] = τ

    if op isa TwoStream
        @inbounds op.ssa[glaycol...] = ssa
        @inbounds op.g[glaycol...] = FT(0) # initializing asymmetry parameter
    end
    return nothing
end
