using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
import NCDatasets # loads the lookup-table extension

# LookupBundle round-trip through save_lookup_tables/load_lookup_tables: the
# cache must reproduce the NetCDF-built bundle exactly (structure and values),
# so the spectral methods can run from the cache without NCDatasets.
@testset "lookup-table serialization round-trip" begin
    context = ClimaComms.context()
    gp = RRTMGP.RRTMGPGridParams(
        Float64;
        context,
        domain_nlay = 10,
        ncol = 2,
    )
    method = RRTMGP.AllSkyRadiationWithClearSkyDiagnostics(true, false)
    bundle = RRTMGP.lookup_tables(gp, method)
    @test bundle isa RRTMGP.LookupBundle

    path = joinpath(mktempdir(), "lookups.jls")
    RRTMGP.save_lookup_tables(path, bundle)
    loaded = RRTMGP.load_lookup_tables(path, gp)

    @test loaded isa RRTMGP.LookupBundle
    @test loaded.nbnd_lw == bundle.nbnd_lw
    @test loaded.nbnd_sw == bundle.nbnd_sw
    @test loaded.ngas_lw == bundle.ngas_lw
    @test loaded.idx_gases_sw == bundle.idx_gases_sw
    @test loaded.idx_aerosol_sw == bundle.idx_aerosol_sw
    # spot-check table payloads (host copies must match the originals)
    @test Array(loaded.lookup_lw.kmajor) == Array(bundle.lookup_lw.kmajor)
    @test Array(loaded.lookup_sw.solar_src_scaled) ==
          Array(bundle.lookup_sw.solar_src_scaled)
    @test !isnothing(loaded.lookup_lw_cld) && !isnothing(loaded.lookup_sw_aero)
    # gray bundles round-trip too (all-nothing tables)
    gray = RRTMGP.lookup_tables(gp, RRTMGP.GrayRadiation())
    path2 = joinpath(mktempdir(), "gray.jls")
    RRTMGP.save_lookup_tables(path2, gray)
    gray2 = RRTMGP.load_lookup_tables(path2, gp)
    @test gray2.nbnd_lw == 1 && isnothing(gray2.lookup_lw)
end
