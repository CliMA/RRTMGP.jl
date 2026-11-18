using Test
import Adapt
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
import NCDatasets # loads the lookup-table extension
import Serialization

# LookupBundle round-trip through save_lookup_tables/load_lookup_tables: the
# cache must reproduce the NetCDF-built bundle exactly (structure and values),
# so the spectral methods can run from the cache without NCDatasets.
@testset "lookup-table serialization round-trip" begin
    context = ClimaComms.context()
    gp = RRTMGP.RRTMGPGridParams(Float64; context, domain_nlay = 10, ncol = 2)
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
    @test !isnothing(loaded.lookup_lw_precip) &&
          !isnothing(loaded.lookup_sw_precip)
    # cloud lookup tables remain pure liquid + ice (4 bounds)
    for cld_lkp in (loaded.lookup_lw_cld, loaded.lookup_sw_cld)
        @test cld_lkp isa RRTMGP.LookUpTables.LookUpCld
        @test length(Array(cld_lkp.bounds)) == 4
    end
    # verify standalone snow and rain lookup tables are populated and preserved
    for precip_lkp in (loaded.lookup_lw_precip, loaded.lookup_sw_precip)
        @test precip_lkp isa RRTMGP.LookUpTables.LookUpPrecip
        precip_h = Adapt.adapt(Array, precip_lkp)
        @test RRTMGP.LookUpTables.get_nsize_sno(precip_h) == 64
        @test RRTMGP.LookUpTables.get_nsize_rai(precip_h) == 64
        bnds = precip_h.bounds
        @test length(bnds) == 4
        @test bnds[1] ≈ 5.0 && bnds[2] ≈ 1100.0
        @test bnds[3] ≈ 5.0 && bnds[4] ≈ 5000.0
        extsno, ssasno, asysno =
            RRTMGP.LookUpTables.getview_snodata(precip_h, 1)
        extrai, ssarai, asyrai =
            RRTMGP.LookUpTables.getview_raidata(precip_h, 1)
        # Mass-specific extinction decreases with particle radius across the table
        @test extsno[end] < extsno[1]
        @test extrai[end] < extrai[1]
        @test all(0 .<= ssasno .<= 1)
        @test all(0 .<= asysno .<= 1)
        @test all(0 .<= ssarai .<= 1)
        @test all(0 .<= asyrai .<= 1)
    end
    # gray bundles round-trip too (all-nothing tables)
    gray = RRTMGP.lookup_tables(gp, RRTMGP.GrayRadiation())
    path2 = joinpath(mktempdir(), "gray.jls")
    RRTMGP.save_lookup_tables(path2, gray)
    gray2 = RRTMGP.load_lookup_tables(path2, gp)
    @test gray2.nbnd_lw == 1 &&
          isnothing(gray2.lookup_lw) &&
          isnothing(gray2.lookup_lw_precip)

    # an aerosol-only (clear-sky) bundle keeps its aerosol tables and index
    # maps through the round-trip, with no cloud or precipitation tables
    clear_aero = RRTMGP.lookup_tables(gp, RRTMGP.ClearSkyRadiation(true))
    path3 = joinpath(mktempdir(), "clear_aero.jls")
    RRTMGP.save_lookup_tables(path3, clear_aero)
    clear_aero2 = RRTMGP.load_lookup_tables(path3, gp)
    @test isnothing(clear_aero2.lookup_lw_cld) &&
          isnothing(clear_aero2.lookup_sw_cld)
    @test isnothing(clear_aero2.lookup_lw_precip) &&
          isnothing(clear_aero2.lookup_sw_precip)
    @test !isnothing(clear_aero2.lookup_lw_aero) &&
          !isnothing(clear_aero2.lookup_sw_aero)
    @test clear_aero2.idx_aerosol_lw == clear_aero.idx_aerosol_lw
    @test clear_aero2.idx_aerosize_sw == clear_aero.idx_aerosize_sw

    # a file that was not written by save_lookup_tables is rejected with an
    # actionable error
    path4 = joinpath(mktempdir(), "not_a_bundle.jls")
    Serialization.serialize(path4, (; a = 1))
    @test_throws "did not deserialize to a LookupBundle" RRTMGP.load_lookup_tables(
        path4,
        gp,
    )
end
