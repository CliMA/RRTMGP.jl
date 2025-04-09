module RTESolver

import ClimaComms
using Adapt

using ..AngularDiscretizations
using ..Vmrs
using ..AtmosphericStates
using ..Fluxes
using ..BCs
using ..Sources
using ..RTE
using ..Optics
using ..LookUpTables

export solve_lw!, solve_sw!

include("longwave1scalar.jl")
include("longwave2stream.jl")
include("shortwave1scalar.jl")
include("shortwave2stream.jl")


"""
    solve_lw!((; context, flux, src, bcs, op)::NoScatLWRTE, as::GrayAtmosphericState, metric_scaling::M = nothing)

Non-scattering RTE solver for the longwave problem, using gray optics. 
Additionally, takes an optional argument `metric_scaling` which scales the resultant fluxes by the 
corresponding factor to account for column expansion in `deep-atmosphere` configurations.
"""
function solve_lw!(
    (; context, flux, src, bcs, op, angle_disc)::NoScatLWRTE,
    as::GrayAtmosphericState,
    metric_scaling::M = nothing,
) where {M}
    rte_lw_noscat_solve!(context.device, flux, src, bcs, op, angle_disc, as)
    apply_metric_scaling!(flux, metric_scaling)
end

"""
    solve_lw!((; context, flux, src, bcs, op)::TwoStreamLWRTE, as::GrayAtmosphericState, metric_scaling::M = nothing)

`Two Stream` RTE solver for the longwave problem, using gray optics.
Additionally, takes an optional argument `metric_scaling` which scales the resultant fluxes by the 
corresponding factor to account for column expansion in `deep-atmosphere` configurations.
"""
function solve_lw!(
    (; context, flux, src, bcs, op)::TwoStreamLWRTE,
    as::GrayAtmosphericState,
    metric_scaling::M = nothing,
) where {M}
    rte_lw_2stream_solve!(context.device, flux, src, bcs, op, as)
    apply_metric_scaling!(flux, metric_scaling)
end

"""
    solve_lw!(
        (; context, fluxb, flux, src, bcs, op)::NoScatLWRTE,
        as::AtmosphericState,
        lookup_lw::LookUpLW,
        lookup_lw_cld::Union{LookUpCld, Nothing},
        lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
        metric_scaling = nothing
    )

Non-scattering RTE solver for the longwave problem, using RRTMGP optics.
Additionally, takes an optional argument `metric_scaling` which scales the resultant fluxes by the 
corresponding factor to account for column expansion in `deep-atmosphere` configurations.
"""
function solve_lw!(
    (; context, fluxb, flux, src, bcs, op, angle_disc)::NoScatLWRTE,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
    metric_scaling::M = nothing,
) where {M}
    rte_lw_noscat_solve!(
        context.device,
        fluxb,
        flux,
        src,
        bcs,
        op,
        angle_disc,
        as,
        lookup_lw,
        lookup_lw_cld,
        lookup_lw_aero,
    )
    apply_metric_scaling!(flux, metric_scaling)
end

"""
    solve_lw!(
        (; context, fluxb, flux, src, bcs, op)::TwoStreamLWRTE,
        as::AtmosphericState,
        lookup_lw::LookUpLW,
        lookup_lw_cld::Union{LookUpCld, Nothing},
        lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
        metric_scaling = nothing
    )

`Two Stream` RTE solver for the longwave problem, using RRTMGP optics.
Additionally, takes an optional argument `metric_scaling` which scales the resultant fluxes by the 
corresponding factor to account for column expansion in `deep-atmosphere` configurations.
"""
function solve_lw!(
    (; context, fluxb, flux, src, bcs, op)::TwoStreamLWRTE,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
    metric_scaling::M = nothing,
) where {M}
    rte_lw_2stream_solve!(context.device, fluxb, flux, src, bcs, op, as, lookup_lw, lookup_lw_cld, lookup_lw_aero)
    apply_metric_scaling!(flux, metric_scaling)
end

"""
    solve_sw!((; context, flux, bcs, op)::NoScatSWRTE, as::GrayAtmosphericState)

Non-scattering RTE solver for the shortwave problem, using gray optics.
Additionally, takes an optional argument `metric_scaling` which scales the resultant fluxes by the 
corresponding factor to account for column expansion in `deep-atmosphere` configurations.
"""
function solve_sw!(
    (; context, flux, bcs, op)::NoScatSWRTE,
    as::GrayAtmosphericState,
    metric_scaling::M = nothing,
) where {M}
    rte_sw_noscat_solve!(context.device, flux, op, bcs, as) # non-scattering solver, gray optics
    apply_metric_scaling!(flux, metric_scaling)
end

"""
    solve_sw!((; context, flux, src, bcs, op)::TwoStreamSWRTE, as::GrayAtmosphericState, metric_scaling = nothing)

`Two Stream` RTE solver for the shortwave problem, using gray optics.
Additionally, takes an optional argument `metric_scaling` which scales the resultant fluxes by the 
corresponding factor to account for column expansion in `deep-atmosphere` configurations.
"""
function solve_sw!(
    (; context, flux, src, bcs, op)::TwoStreamSWRTE,
    as::GrayAtmosphericState,
    metric_scaling::M = nothing,
) where {M}
    rte_sw_2stream_solve!(context.device, flux, op, bcs, src, as)
    apply_metric_scaling!(flux, metric_scaling)
end

"""
    solve_sw!(
        (; context, fluxb, flux, bcs, op)::NoScatSWRTE,
        as::AtmosphericState,
        lookup_sw::LookUpSW,
        metric_scaling = nothing,
    )

Non-scattering RTE solver for the shortwave problem, using RRTMGP optics.
Additionally, takes an optional argument `metric_scaling` which scales the resultant fluxes by the 
corresponding factor to account for column expansion in `deep-atmosphere` configurations.
"""
function solve_sw!(
    (; context, fluxb, flux, bcs, op)::NoScatSWRTE,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    metric_scaling::M = nothing,
) where {M}
    rte_sw_noscat_solve!(context.device, fluxb, flux, op, bcs, as, lookup_sw)
    apply_metric_scaling!(flux, metric_scaling)
end

"""
    solve_sw!(
        (; context, fluxb, flux, src, bcs, op)::TwoStreamSWRTE,
        as::AtmosphericState,
        lookup_sw::LookUpSW,
        lookup_sw_cld::Union{LookUpCld, Nothing},
        lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
        metric_scaling = nothing
    )

`Two Stream` RTE solver for the shortwave problem, using RRTMGP optics.
Additionally, takes an optional argument `metric_scaling` which scales the resultant fluxes by the 
corresponding factor to account for column expansion in `deep-atmosphere` configurations.
"""
function solve_sw!(
    (; context, fluxb, flux, src, bcs, op)::TwoStreamSWRTE,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
    metric_scaling::M = nothing,
) where {M}
    rte_sw_2stream_solve!(context.device, fluxb, flux, op, bcs, src, as, lookup_sw, lookup_sw_cld, lookup_sw_aero)
    apply_metric_scaling!(flux, metric_scaling)
end

end
