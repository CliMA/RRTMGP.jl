using Test
using RRTMGP

# The public API is the version-compatibility contract: `RRTMGP.PUBLIC_NAMES`
# lists it, and this file locks it down from both sides. A name added to or
# removed from the module surface fails here until it is classified, so the
# public surface cannot grow by accident.
#
# Internal names are excluded by convention: a leading `_`, a submodule, or a
# deprecated name.

const DEPRECATED_NAMES = (
    :clear_lw_flux,
    :clear_sw_flux,
    :toa_flux,
    :top_of_atmosphere_lw_flux_dn,
    :top_of_atmosphere_diffuse_sw_flux_dn,
)

# Every name the module exposes that is neither internal nor deprecated.
function exposed_names()
    out = Symbol[]
    for n in names(RRTMGP; all = true, imported = false)
        s = String(n)
        (startswith(s, "_") || startswith(s, "#")) && continue
        n in (:RRTMGP, :eval, :include, :PUBLIC_NAMES) && continue
        n in DEPRECATED_NAMES && continue
        isdefined(RRTMGP, n) || continue
        getfield(RRTMGP, n) isa Module && continue
        push!(out, n)
    end
    return Set(out)
end

@testset "public API list is complete and defined" begin
    @test allunique(RRTMGP.PUBLIC_NAMES)
    for n in RRTMGP.PUBLIC_NAMES
        @test isdefined(RRTMGP, n)
    end
    # Nothing public may hide behind an underscore.
    @test !any(n -> startswith(String(n), "_"), RRTMGP.PUBLIC_NAMES)
end

@testset "no undeclared names on the module surface" begin
    undeclared = sort(collect(setdiff(exposed_names(), RRTMGP.PUBLIC_NAMES)))
    if !isempty(undeclared)
        @info "Names exposed by RRTMGP but not in PUBLIC_NAMES. Either add them \
               to `src/api/public.jl` or make them internal (`_`-prefix)." undeclared
    end
    @test isempty(undeclared)
end

@testset "every public name has a docstring" begin
    # Check the module's documentation table by binding, not `Docs.doc`: the
    # binding-aware `doc` method needs REPL loaded, and resolving the value
    # would miss docstrings attached to a constant's binding (`check_values`).
    documented = Base.Docs.meta(RRTMGP)
    undocumented = filter(
        n -> !haskey(documented, Base.Docs.Binding(RRTMGP, n)),
        collect(RRTMGP.PUBLIC_NAMES),
    )
    isempty(undocumented) || @info "Public names without docstrings" undocumented
    @test isempty(undocumented)
end

@testset "deprecated names still resolve to their replacements" begin
    out = RRTMGP.solve_gray(Float64; nlay = 10, ncol = 2)
    s = out.solver
    @test (@test_deprecated RRTMGP.toa_flux(s)) === RRTMGP.toa_sw_flux_dn(s)
    @test (@test_deprecated RRTMGP.top_of_atmosphere_lw_flux_dn(s)) ===
          RRTMGP.toa_lw_flux_dn(s)
end

nothing
