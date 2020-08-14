module GrayRTESolver

using ..Device: array_type
using ..GrayAngularDiscretizations
using ..GraySources
using ..GrayAtmos
using ..GrayOptics

export rte_lw_noscat_gray_solve!, rte_lw_2stream_gray_solve!, adding_gray!

using KernelAbstractions
using CUDA

array_device(::Union{Array}) = CPU()
array_device(::CuArray) = CUDADevice()


include("GrayRTESolverKernels.jl")
# -------------------------------------------------------------------------------------------------
#
# LW fluxes, no scattering, mu (cosine of integration angle) specified by column
#   Does radiation calculation at user-supplied angles; converts radiances to flux
#   using user-supplied weights
#
# ---------------------------------------------------------------    
function rte_lw_noscat_gray_solve!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D,FTA3D,B},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    I<:Int,
    B<:Bool,
}
    # Setting references 
    flux_up = gray_rrtmgp.flux.flux_up     # upward flux
    flux_dn = gray_rrtmgp.flux.flux_dn     # downward flux
    flux_net = gray_rrtmgp.flux.flux_net   # net flux
    top_at_1 = gray_rrtmgp.as.top_at_1     # is top-of-atmos located at point # 1?
    source = gray_rrtmgp.src               # Planck sources
    sfc_emis = gray_rrtmgp.bcs.sfc_emis    # Surface emissivity
    inc_flux = gray_rrtmgp.bcs.inc_flux    # Surface flux
    τ = gray_rrtmgp.op.τ                   # Optical thickness

    ncol = size(flux_up, 1)         # number of columns
    nlev = size(flux_up, 2)         # number of levels
    nlay = nlev - 1                 # number of layers
    ngpt = 1                        # number of g-points (only 1 for gray radiation model

    i_lev_top = top_at_1 ? 1 : nlev # index for top level of column
    n_sfc = top_at_1 ? nlev : 1     # index for bottom level of column
    n_μ = gray_rrtmgp.angle_disc.n_gauss_angles
    Ds = gray_rrtmgp.angle_disc.gauss_Ds
    w_μ = gray_rrtmgp.angle_disc.gauss_wts

    lev_source_up = top_at_1 ? source.lev_source_dec : source.lev_source_inc # Mapping increasing/decreasing indices to up/down
    lev_source_dn = top_at_1 ? source.lev_source_inc : source.lev_source_dec
    lay_source = source.lay_source # Planck source at average layer temperature [W/m^2]
    sfc_source = source.sfc_source # Surface source function [W/m^2]
    source_up = source.source_up
    source_dn = source.source_dn
    trans = source.trans

    #------Launching KA kernel-----------------------
    println("timing source kernel")
    @time begin
        max_threads = 256
        thr_x = min(32, ncol)
        thr_y = min(Int(floor(FT(max_threads / thr_x))), nlay)

        device = array_device(flux_up)
        workgroup = (thr_x, thr_y)
        ndrange = (ncol, nlay)

        kernel = rte_lw_noscat_gray_source_kernel!(device, workgroup)
        event = kernel(
            lay_source,
            lev_source_up,
            lev_source_dn,
            τ,
            trans,
            source_up,
            source_dn,
            Ds,
            Val(ncol),
            Val(nlay),
            ndrange = ndrange,
        )
        wait(event)
    end
    println("-----------------------------------")

    println("timing transport kernel")
    @time begin
        thr_x = 256

        device = array_device(flux_up)
        workgroup = (thr_x)
        ndrange = (ncol)

        kernel = rte_lw_noscat_gray_transport_kernel!(device, workgroup)

        event = kernel(
            flux_up,
            flux_dn,
            flux_net,
            trans,
            source_up,
            source_dn,
            sfc_source,
            sfc_emis,
            inc_flux,
            τ,
            n_μ,
            Ds,
            w_μ,
            top_at_1,
            Val(ncol),
            Val(nlev),
            Val(nlay),
            ndrange = ndrange,
        )
        wait(event)
    end
    println("-----------------------------------")
    #------------------------------------------------
    return nothing
end

function rte_lw_2stream_gray_solve!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D,FTA3D,B},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    I<:Int,
    B<:Bool,
}
    # Setting references
    flux_up = gray_rrtmgp.flux.flux_up
    flux_dn = gray_rrtmgp.flux.flux_dn
    flux_net = gray_rrtmgp.flux.flux_net
    top_at_1 = gray_rrtmgp.as.top_at_1
    source = gray_rrtmgp.src
    sfc_emis = gray_rrtmgp.bcs.sfc_emis
    inc_flux = gray_rrtmgp.bcs.inc_flux
    ncol, nlev, ngpt = size(flux_up, 1), size(flux_up, 2), 1
    nlay = nlev - 1
    lev_src_inc, lev_src_dec = source.lev_source_inc, source.lev_source_dec
    lev_source = source.lev_source
    lay_source, sfc_source = source.lay_source, source.sfc_source
    i_lev_top = top_at_1 ? 1 : nlev
    n_sfc = top_at_1 ? nlev : 1

    τ = gray_rrtmgp.op.τ
    g = gray_rrtmgp.op.g
    ssa = gray_rrtmgp.op.ssa

    Rdif, Tdif = source.Rdif, source.Tdif
    src_up, src_dn = source.src_up, source.src_dn # Source function for diffuse radiation
    albedo, src = source.albedo, source.src
    #----------------------------------------------------------------------------------
    fill!(flux_dn, 0)
    fill!(flux_up, 0)
    if inc_flux ≠ nothing
        flux_dn[1, i_lev_top] = inc_flux[1]
    end
    #-----------------------------------------------------------------------------------
    lw_diff_sec = FT(1.66)
    #-----------------------------------------------------------------------------------

    for icol = 1:ncol
        # RRTMGP provides source functions at each level using the spectral mapping
        # of each adjacent layer. Combine these for two-stream calculations
        # lw combine sources
        lev_source[icol, 1] = lev_src_dec[icol, 1]
        lev_source[icol, nlev] = lev_src_inc[icol, nlay]

        for ilay = 2:nlay
            lev_source[icol, ilay] =
                sqrt(lev_src_dec[icol, ilay] * lev_src_inc[icol, ilay-1])
        end
        #-----------------------------------------------
    end
    # Cell properties: reflection, transmission for diffuse radiation
    # Coupling coefficients needed for source function
    # -------------------------------------------------------------------------------------------------
    #
    # Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
    #    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
    #
    # Equations are developed in Meador and Weaver, 1980,
    #    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
    #
    # -------------------------------------------------------------------------------------------------    
    for icol = 1:ncol
        for ilay = 1:nlay
            γ1 =
                lw_diff_sec *
                (1.0 - 0.5 * ssa[icol, ilay] * (1.0 + g[icol, ilay]))
            γ2 = lw_diff_sec * 0.5 * ssa[icol, ilay] * (1.0 - g[icol, ilay])
            k = sqrt(max((γ1 + γ2) * (γ1 - γ2), 1e-12))

            # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
            RT_term =
                1.0 / (
                    k * (1.0 + exp(-2.0 * τ[icol, ilay] * k)) +
                    γ1 * (1.0 - exp(-2 * τ[icol, ilay] * k))
                )

            Rdif[ilay] = RT_term * γ2 * (1.0 - exp(-2 * τ[icol, ilay] * k)) # Equation 25
            Tdif[ilay] = RT_term * 2.0 * k * exp(-τ[icol, ilay] * k) # Equation 26

            # Source function for diffuse radiation
            # Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
            # This version straight from ECRAD
            # Source is provided as W/m2-str; factor of pi converts to flux units
            # lw_source_2str
            if top_at_1
                top = ilay
                bot = ilay + 1
            else
                top = ilay + 1
                bot = ilay
            end

            if τ[icol, ilay] > 1.0e-8
                # Toon et al. (JGR 1989) Eqs 26-27
                Z =
                    (lev_source[icol, bot] - lev_source[icol, top]) /
                    (τ[icol, ilay] * (γ1 + γ2))
                Zup_top = Z + lev_source[icol, top]
                Zup_bottom = Z + lev_source[icol, bot]
                Zdn_top = -Z + lev_source[icol, top]
                Zdn_bottom = -Z + lev_source[icol, bot]
                src_up[ilay] =
                    pi *
                    (Zup_top - Rdif[ilay] * Zdn_top - Tdif[ilay] * Zup_bottom)
                src_dn[ilay] =
                    pi * (
                        Zdn_bottom - Rdif[ilay] * Zup_bottom -
                        Tdif[ilay] * Zdn_top
                    )
            else
                src_up[ilay] = 0.0
                src_dn[ilay] = 0.0
            end
        end

        # Transport
        adding_gray!(
            top_at_1,
            sfc_emis,
            Rdif,
            Tdif,
            src_up,
            src_dn,
            sfc_source,
            flux_up,
            flux_dn,
            albedo,
            src,
            icol,
        )

        for ilev = 1:nlev
            flux_net[icol, ilev] = flux_up[icol, ilev] - flux_dn[icol, ilev]
        end
    end


    return nothing
end

# ---------------------------------------------------------------
#
# Transport of diffuse radiation through a vertically layered atmosphere.
#   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
#   This routine is shared by longwave and shortwave
#
# -------------------------------------------------------------------------------------------------
function adding_gray!(
    top_at_1::B,
    sfc_emis::FTA1D,
    Rdif::FTA1D,
    Tdif::FTA1D,
    src_up::FTA1D,
    src_dn::FTA1D,
    sfc_source::FTA2D,
    flux_up::FTA2D,
    flux_dn::FTA2D,
    albedo::FTA1D,
    src::FTA1D,
    icol::I,
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    I<:Int,
    B<:Bool,
}
    ncol = size(flux_up, 1)
    nlev = size(flux_up, 2)
    nlay = nlev - 1

    if top_at_1
        #--------------------------------------------------------------------------------------
        # Albedo of lowest level is the surface albedo...
        albedo[nlev] = FT(1) - sfc_emis[icol]
        # ... and source of diffuse radiation is surface emission
        src[nlev] = FT(π) * sfc_emis[icol] * sfc_source[icol]
        # From bottom to top of atmosphere --
        # compute albedo and source of upward radiation
        for ilev = nlay:-1:1
            denom = FT(1) / (FT(1) - Rdif[ilev] * albedo[ilev+1])
            albedo[ilev] =
                Rdif[ilev] + (Tdif[ilev] * Tdif[ilev] * albedo[ilev+1] * denom)
            # Equation 11 -- source is emitted upward radiation at top of layer plus
            # radiation emitted at bottom of layer,
            # transmitted through the layer and reflected from layers below (tdiff*src*albedo)
            src[ilev] =
                src_up[ilev] +
                Tdif[ilev] *
                denom *
                (src[ilev+1] + albedo[ilev+1] * src_dn[ilev])
        end

        # Eq 12, at the top of the domain upwelling diffuse is due to ...
        flux_up[icol, 1] = flux_dn[icol, 1] * albedo[1] + src[1] # reflection of incident diffuse and
        # emission from below

        # From the top of the atmosphere downward -- compute fluxes
        for ilev = 2:nlay+1
            denom = FT(1) / (FT(1) - Rdif[ilev-1] * albedo[ilev])
            flux_dn[icol, ilev] =
                (
                    Tdif[ilev-1] * flux_dn[icol, ilev-1] + # Equation 13
                    Rdif[ilev-1] * src[ilev] +
                    src_dn[ilev-1]
                ) * denom
            flux_up[icol, ilev] =
                flux_dn[icol, ilev] * albedo[ilev] + # Equation 12
                src[ilev]
        end
        #--------------------------------------------------------------------------------------
    else
        #--------------------------------------------------------------------------------------
        # Albedo of lowest level is the surface albedo...
        albedo[1] = FT(1) - sfc_emis[icol]
        # ... and source of diffuse radiation is surface emission
        src[1] = FT(π) * sfc_emis[icol] * sfc_source[icol]
        #--------------------------------------------------------
        # From bottom to top of atmosphere --
        #   compute albedo and source of upward radiation
        for ilev = 1:nlay
            denom = FT(1) / (FT(1) - Rdif[ilev] * albedo[ilev])  # Eq 10
            albedo[ilev+1] =
                Rdif[ilev] + Tdif[ilev] * Tdif[ilev] * albedo[ilev] * denom # Equation 9
            # 
            # Equation 11 -- source is emitted upward radiation at top of layer plus
            # radiation emitted at bottom of layer,
            # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
            src[ilev+1] =
                src_up[ilev] +
                Tdif[ilev] * denom * (src[ilev] + albedo[ilev] * src_dn[ilev])
        end

        # Eq 12, at the top of the domain upwelling diffuse is due to ...
        flux_up[icol, nlev] =
            flux_dn[icol, nlev] * albedo[nlev] + # ... reflection of incident diffuse and
            src[nlev]                          # scattering by the direct beam below

        # From the top of the atmosphere downward -- compute fluxes
        for ilev = nlay:-1:1
            denom = FT(1) / (FT(1) - Rdif[ilev] * albedo[ilev])  # Eq 10
            flux_dn[icol, ilev] =
                (
                    Tdif[ilev] * flux_dn[icol, ilev+1] + # Equation 13
                    Rdif[ilev] * src[ilev] +
                    src_dn[ilev]
                ) * denom
            flux_up[icol, ilev] =
                flux_dn[icol, ilev] * albedo[ilev] + # Equation 12
                src[ilev]

        end
        #--------------------------------------------------------
    end
    #--------------------------------------------------------------------------------------

    return nothing
end

end
