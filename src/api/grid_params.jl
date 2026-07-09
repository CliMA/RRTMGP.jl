using ClimaComms


"""
    RRTMGPGridParams(
        FT;
        context::ClimaComms.AbstractCommsContext,
        domain_nlay::Int,
        ncol::Int,
        isothermal_boundary_layer::Bool = false,
    )

Grid parameters for RRTMGP, parametrized on the float type `FT`.

Specify `domain_nlay`, the number of layers in your **physical** grid. When
`isothermal_boundary_layer = true`, RRTMGP adds one extra layer at the top of the
domain *internally*, so the stored field `nlay` — the total number of layers RRTMGP
works with — is `domain_nlay + 1` (and just `domain_nlay` when the flag is `false`).
Callers never add the extra layer themselves, so the boundary-layer bookkeeping
cannot be got wrong; the getters return domain-sized arrays and `heating_rate`
reports on the physical grid.

# Keyword Arguments
 - `context`: the `ClimaComms` context.
 - `domain_nlay`: the number of layers in the physical domain.
 - `ncol`: the number of columns.
 - `isothermal_boundary_layer = false`: whether RRTMGP adds an isothermal boundary
    layer/level at the top of the domain.

"""
struct RRTMGPGridParams{FT, C}
    context::C
    "total layer count RRTMGP works with = `domain_nlay + Int(isothermal_boundary_layer)`"
    nlay::Int
    ncol::Int
    isothermal_boundary_layer::Bool # for isothermal boundary layer
end
function RRTMGPGridParams(
    ::Type{FT};
    context::ClimaComms.AbstractCommsContext,
    domain_nlay::Int,
    ncol::Int,
    isothermal_boundary_layer::Bool = false,
) where {FT}
    # Add the isothermal boundary layer here so callers only ever pass the physical
    # `domain_nlay`; the stored `nlay` is the total and cannot be got wrong.
    nlay = domain_nlay + Int(isothermal_boundary_layer)
    return RRTMGPGridParams{FT, typeof(context)}(
        context,
        nlay,
        ncol,
        isothermal_boundary_layer,
    )
end

# Overload some methods for convenience:
Base.eltype(s::RRTMGPGridParams{FT}) where {FT} = FT
ClimaComms.device(s::RRTMGPGridParams) = ClimaComms.device(s.context)
ClimaComms.array_type(s::RRTMGPGridParams) =
    ClimaComms.array_type(ClimaComms.device(s))


