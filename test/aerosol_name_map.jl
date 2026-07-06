using Test
import RRTMGP

@testset "aerosol name <-> index map" begin
    idx = RRTMGP.aerosol_idx()
    names = RRTMGP.aerosol_names()

    # There are 15 MERRA aerosol species.
    @test length(idx) == 15
    @test length(names) == 15

    # `aerosol_names()` is ordered by index: names[i] is the species at index i,
    # every index 1..15 is covered exactly once, and the name set is consistent.
    for (i, name) in enumerate(names)
        @test idx[name] == i
    end
    @test sort(collect(values(idx))) == collect(1:15)
    @test Set(names) == Set(keys(idx))

    # Kernel contract: these indices are hard-coded in
    # src/optics/aerosol_optics.jl (compute_lookup_aerosol). If the canonical
    # map changes, that kernel must change with it.
    @test [idx["dust$i"] for i in 1:5] == [1, 8, 9, 10, 11]
    @test [idx["sea_salt$i"] for i in 1:5] == [2, 12, 13, 14, 15]
    @test idx["sulfate"] == 3
    @test idx["black_carbon_rh"] == 4
    @test idx["black_carbon"] == 5
    @test idx["organic_carbon_rh"] == 6
    @test idx["organic_carbon"] == 7

    @test RRTMGP.canonical_aerosol_name("dust3") == "dust3" # canonical passes through


    # The order ClimaAtmos fills `aero_mass`/`aero_size` positionally must map
    # 1:1 onto canonical indices 1..15. This documents and locks why the current
    # ClimaAtmos aerosol code is correct (see RRTMGPInterface.jl).
    climaatmos_order = [
        "dust1",
        "sea_salt1",
        "sulfate",
        "black_carbon_rh",
        "black_carbon",
        "organic_carbon_rh",
        "organic_carbon",
        "dust2",
        "dust3",
        "dust4",
        "dust5",
        "sea_salt2",
        "sea_salt3",
        "sea_salt4",
        "sea_salt5",
    ]
    @test [RRTMGP.aerosol_index(n) for n in climaatmos_order] == collect(1:15)

    # Unknown names error.
    @test_throws ErrorException RRTMGP.canonical_aerosol_name("not_an_aerosol")
    @test_throws ErrorException RRTMGP.aerosol_index("xyz")
end
