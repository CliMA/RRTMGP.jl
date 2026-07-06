using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
using NCDatasets

# The spectral half of the standalone front door: `solve(profile)` with the
# default `ClearSkyRadiation`. Runs late in the suite because it needs the
# NCDatasets extension (test/standalone.jl runs before NCDatasets loads, to
# keep exercising the no-weakdeps gray path).
@testset "clear-sky solve(profile) (spectral standalone)" begin
    FT = Float64
    nlay = 50
    prof = RRTMGP.standard_atmosphere(FT; kind = :midlatitude_summer, nlay)
    out = RRTMGP.solve(prof)
    @test out isa RRTMGP.RadiationOutput
    @test size(Array(out.lw_up)) == (nlay + 1, 1)
    @test size(Array(out.heating_rate)) == (nlay, 1)
    for f in (
        out.lw_up,
        out.lw_dn,
        out.sw_up,
        out.sw_dn,
        out.sw_direct_dn,
        out.net,
        out.heating_rate,
    )
        @test all(isfinite, Array(f))
    end
    @test Array(out.net) ≈ Array(out.lw_net) .+ Array(out.sw_net)
    # clear-sky OLR for an idealized midlatitude-summer column
    olr = Array(out.lw_up)[end, 1]
    @test 220 < olr < 320
    # TOA downwelling shortwave is the solar constant times cos(zenith)
    @test Array(out.sw_dn)[end, 1] ≈ 1361 * 0.5 rtol = 1e-3
    # the direct beam never exceeds the total downwelling flux
    @test all(Array(out.sw_direct_dn) .<= Array(out.sw_dn) .* (1 + eps(FT)))

    # reusing the lookup tables skips the NetCDF re-read; a warmer atmosphere
    # emits more
    lookups = out.solver.lookups
    olr_of(kind) =
        Array(
            RRTMGP.solve(RRTMGP.standard_atmosphere(FT; kind); lookups).lw_up,
        )[
            end,
            1,
        ]
    @test olr_of(:tropical) > olr_of(:subarctic_winter)

    # Float32 runs end-to-end (its own lookups: the tables are FT-typed)
    out32 = RRTMGP.solve(RRTMGP.standard_atmosphere(Float32; nlay = 40))
    @test eltype(Array(out32.net)) == Float32
    @test all(isfinite, Array(out32.net))

    # gas-name validation of `well_mixed_vmr`
    prof_bad = RRTMGP.standard_atmosphere(FT; nlay = 10)
    prof_bad.well_mixed_vmr["h2o"] = FT(0.01) # spatially varying, not well-mixed
    @test_throws ErrorException RRTMGP.solve(prof_bad; lookups)
    delete!(prof_bad.well_mixed_vmr, "h2o")
    prof_bad.well_mixed_vmr["kryptonite"] = FT(1e-6) # not an RRTMGP gas
    @test_throws ErrorException RRTMGP.solve(prof_bad; lookups)
end
