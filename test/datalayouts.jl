#=
julia --check-bounds=yes --project=test
using Revise; include("test/datalayouts.jl")
=#
using RRTMGP: RRTMGPGridParams
using RRTMGP
using ClimaComms
ClimaComms.@import_required_backends
using Test

@testset "Datalayouts" begin
    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    nlay = 10
    ncol = 20
    gp = RRTMGPGridParams(Float64; context, nlay, ncol)

    data_NVC = RRTMGP.NVCData(gp; N = 2, nlay)
    @test size(data_NVC) == (2, nlay, ncol)
    @test size(RRTMGP.domain_view(true, data_NVC), 2) + 1 == size(RRTMGP.domain_view(false, data_NVC), 2)

    data_VC = RRTMGP.VCData(gp; nlay)
    @test size(data_VC) == (nlay, ncol)
    @test size(RRTMGP.domain_view(true, data_VC), 1) + 1 == size(RRTMGP.domain_view(false, data_VC), 1)

    RRTMGP.set_domain!(data_VC, DA(rand(nlay, ncol)), gp)

    data_NC = RRTMGP.NCData(gp; N = 2)
    @test size(data_NC) == (2, ncol)
    @test RRTMGP.domain_view(true, data_NC) === RRTMGP.domain_view(false, data_NC)

    data_N = RRTMGP.NData(gp; N = 2)
    @test size(data_N) == (2,)
    @test RRTMGP.domain_view(true, data_N) === RRTMGP.domain_view(false, data_N)
end

@testset "set_cols!" begin
    context = ClimaComms.context()
    nlay = 10
    ncol = 20
    gp = RRTMGPGridParams(Float64; context, nlay, ncol)

    data_NVC = RRTMGP.NVCData(gp; N = 2, nlay)
    RRTMGP.set_cols!(data_NVC, rand(2, nlay, ncol))
    @test_throws ErrorException RRTMGP.set_cols!(data_NVC, rand(nlay, ncol))
    RRTMGP.set_cols!(data_NVC, rand())

    data_NVC = RRTMGP.NVCData(gp; N = 1, nlay)
    RRTMGP.set_cols!(data_NVC, rand(1, nlay, ncol))
    RRTMGP.set_cols!(data_NVC, rand(nlay, ncol))
    RRTMGP.set_cols!(data_NVC, rand())

    data_VC = RRTMGP.VCData(gp; nlay)
    RRTMGP.set_cols!(data_VC, rand(nlay, ncol))
    RRTMGP.set_cols!(data_VC, rand())

    data_NC = RRTMGP.NCData(gp; N = 2)
    RRTMGP.set_cols!(data_NC, rand(2, ncol))
    RRTMGP.set_cols!(data_NC, rand())

    data_NC = RRTMGP.NCData(gp; N = 1)
    RRTMGP.set_cols!(data_NC, rand(1, ncol))
    RRTMGP.set_cols!(data_NC, rand(ncol))
    RRTMGP.set_cols!(data_NC, rand())

    data_N = RRTMGP.NData(gp; N = 2)
    RRTMGP.set_cols!(data_N, rand(2))
    RRTMGP.set_cols!(data_N, rand())
end
