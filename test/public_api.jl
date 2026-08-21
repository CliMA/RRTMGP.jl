using Test
using RRTMGP

# `RRTMGP.PUBLIC_NAMES` lists the public API; this file locks it from both
# sides, so a name added to or removed from the module surface fails here until
# it is classified. Internal names are excluded by convention: a leading `_` or
# a submodule.

# Every name the module exposes that is not internal.
function exposed_names()
    out = Symbol[]
    for n in names(RRTMGP; all = true, imported = false)
        s = String(n)
        (startswith(s, "_") || startswith(s, "#")) && continue
        n in (:RRTMGP, :eval, :include, :PUBLIC_NAMES) && continue
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
    # binding-aware `doc` method needs the REPL loaded, and resolving the value
    # would miss docstrings attached to a constant's binding (`check_values`).
    documented = Base.Docs.meta(RRTMGP)
    undocumented = filter(
        n -> !haskey(documented, Base.Docs.Binding(RRTMGP, n)),
        collect(RRTMGP.PUBLIC_NAMES),
    )
    isempty(undocumented) || @info "Public names without docstrings" undocumented
    @test isempty(undocumented)
end

nothing
