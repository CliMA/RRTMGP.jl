#=
julia --check-bounds=yes --project=test
using Revise; include("test/datalayouts.jl")
=#
using RRTMGP: RRTMGPGridParams
using RRTMGP
using ClimaComms
ClimaComms.@import_required_backends
using Test
using Random
Random.seed!(1234)

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

    arr = DA(rand(nlay, ncol))
    RRTMGP.set_domain!(data_VC, arr, gp)
    @test parent(data_VC) == arr

    data_NC = RRTMGP.NCData(gp; N = 2)
    @test size(data_NC) == (2, ncol)
    @test RRTMGP.domain_view(true, data_NC) === RRTMGP.domain_view(false, data_NC)

    data_N = RRTMGP.NData(gp; N = 2)
    @test size(data_N) == (2,)
    @test RRTMGP.domain_view(true, data_N) === RRTMGP.domain_view(false, data_N)
end

@testset "set_cols!" begin
    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    nlay = 10
    ncol = 20
    gp = RRTMGPGridParams(Float64; context, nlay, ncol)

    data = RRTMGP.NVCData(gp; N = 2, nlay)
    arr = DA(rand(2, nlay, ncol))
    RRTMGP.set_cols!(data, arr)
    @test parent(data) == arr
    @test_throws ErrorException RRTMGP.set_cols!(data, rand(nlay, ncol))
    let val = rand()
        RRTMGP.set_cols!(data, val)
        @test all(x -> x == val, parent(data))
    end

    data = RRTMGP.NVCData(gp; N = 1, nlay)
    let val = DA(rand(1, nlay, ncol))
        RRTMGP.set_cols!(data, val)
        @test parent(data) == val
    end
    let val = DA(rand(nlay, ncol))
        RRTMGP.set_cols!(data, val)
        @test parent(data[1, :, :]) == val
    end
    let val = rand()
        RRTMGP.set_cols!(data, val)
        @test all(x -> x == val, parent(data))
    end

    data = RRTMGP.VCData(gp; nlay)
    let val = DA(rand(nlay, ncol))
        RRTMGP.set_cols!(data, val)
        @test parent(data) == val
    end
    let val = rand()
        RRTMGP.set_cols!(data, val)
        @test all(x -> x == val, parent(data))
    end

    data = RRTMGP.NCData(gp; N = 2)
    let val = DA(rand(2, ncol))
        RRTMGP.set_cols!(data, val)
        @test parent(data) == val
    end
    let val = rand()
        RRTMGP.set_cols!(data, val)
        @test all(x -> x == val, parent(data))
    end

    data = RRTMGP.NCData(gp; N = 1)
    let val = DA(rand(1, ncol))
        RRTMGP.set_cols!(data, val)
        @test parent(data) == val
    end
    let val = DA(rand(ncol))
        RRTMGP.set_cols!(data, val)
        @test parent(data)[1, :] == val
    end
    let val = rand()
        RRTMGP.set_cols!(data, val)
        @test all(x -> x == val, parent(data))
    end

    data = RRTMGP.NData(gp; N = 2)
    let val = DA(rand(2))
        RRTMGP.set_cols!(data, val)
        @test parent(data) == val
    end
    let val = rand()
        RRTMGP.set_cols!(data, val)
        @test all(x -> x == val, parent(data))
    end
end

@testset "set_domain! with isothermal_boundary_layer=false" begin
    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    nlay = 10
    ncol = 20
    gp = RRTMGPGridParams(Float64; context, nlay, ncol, isothermal_boundary_layer = false)

    data = RRTMGP.NVCData(gp; N = 2, nlay)
    let arr = DA(rand(2, nlay, ncol))
        RRTMGP.set_domain!(data, arr, gp)
        @test parent(data) == arr
    end

    data = RRTMGP.VCData(gp; nlay)
    fill!(parent(data), NaN)
    let arr = DA(rand(nlay, ncol))
        RRTMGP.set_domain!(data, arr, gp)
        @test parent(data) == arr
    end

    data = RRTMGP.NCData(gp; N = 2)
    let arr = DA(rand(2, ncol))
        RRTMGP.set_domain!(data, arr, gp)
        @test parent(data) == arr
    end

    data = RRTMGP.NData(gp; N = 2)
    let arr = rand()
        RRTMGP.set_domain!(data, arr, gp)
        @test all(x -> x == arr, parent(data))
    end
end

@testset "set_domain! with isothermal_boundary_layer=true" begin
    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    nlay = 10
    ncol = 20
    gp = RRTMGPGridParams(Float64; context, nlay, ncol, isothermal_boundary_layer = true)

    data = RRTMGP.NVCData(gp; N = 2, nlay)
    fill!(parent(data), NaN)
    let arr = DA(rand(2, nlay - 1, ncol))
        RRTMGP.set_domain!(data, arr, gp)
        @test parent(data)[:, 1:(nlay - 1), :] == arr[:, 1:(nlay - 1), :]
        @test all(isnan, parent(data)[:, nlay, :])
    end

    data = RRTMGP.VCData(gp; nlay)
    fill!(parent(data), NaN)
    let arr = DA(rand(nlay - 1, ncol))
        RRTMGP.set_domain!(data, arr, gp)
        @test parent(data)[1:(nlay - 1), :] == arr[1:(nlay - 1), :]
        @test all(isnan, parent(data)[nlay, :])
    end

    data = RRTMGP.NCData(gp; N = 2)
    let arr = DA(rand(2, ncol))
        RRTMGP.set_domain!(data, arr, gp)
        @test parent(data) == arr
    end

    data = RRTMGP.NData(gp; N = 2)
    let arr = rand()
        RRTMGP.set_domain!(data, arr, gp)
        @test all(x -> x == arr, parent(data))
    end
end
