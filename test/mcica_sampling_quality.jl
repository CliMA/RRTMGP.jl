using Test
import RRTMGP
import RRTMGP.Optics: _keyed_uniform
import RRTMGP: _mcica_key, _rank_salt
import ClimaComms

_mean(x) = sum(x) / length(x)
_var(x) = sum(abs2, x .- _mean(x)) / (length(x) - 1)
function _cor(a, b)
    A, B = vec(a), vec(b)
    da, db = A .- _mean(A), B .- _mean(B)
    return sum(da .* db) / sqrt(sum(abs2, da) * sum(abs2, db))
end

"""
    mcica_sampling_quality_test(::Type{FT})

Check that the keyed McICA draw behaves like a uniform random sample.

The cloud mask is a hash of `(key, column, g-point, layer)` rather than a draw
from a generator's state, which buys reproducibility but puts the burden of
statistical quality on the mixing function. McICA averages over g-points, so
draws that were uniform but *correlated along g-point* would bias the broadband
flux while leaving every flux-level test happy: a fused-versus-separate
comparison uses the same sampler on both sides, and clear sky involves no
sampling at all. Hence a direct test.

Nothing here touches a GPU, but `runtests.jl` runs on CPU, threaded CPU and GPU,
and the values must be identical in all three -- the draw depends only on its
arguments.
"""
# A stand-in context, so the rank assertions need no MPI.
struct _FakeCtx
    pid::Int
end
ClimaComms.mypid(c::_FakeCtx) = c.pid

function mcica_sampling_quality_test(::Type{FT}) where {FT}
    ngpt, nlay, ncol = 256, 64, 32
    key = 0x5eed1234 % UInt32

    draws = [
        Float64(_keyed_uniform(FT, key, gcol, igpt, ilay)) for igpt in 1:ngpt,
        ilay in 1:nlay, gcol in 1:ncol
    ]
    n = length(draws)

    @testset "range and uniformity" begin
        @test all(0 .<= draws .< 1)
        # Standard error of the mean is 1/sqrt(12n); this is ~17 sigma of slack.
        @test abs(_mean(draws) - 0.5) < 5e-3
        @test abs(_var(draws) - 1 / 12) < 5e-3
        # Sixteen equal bins, each within 5% of its expected count.
        counts = zeros(Int, 16)
        for d in draws
            counts[min(floor(Int, d * 16) + 1, 16)] += 1
        end
        @test all(abs.(counts .- n / 16) .< 0.05 * n / 16)
    end

    # Correlation along each index separately. g-point is the one that matters --
    # it is McICA's sample index, the axis the broadband flux averages over --
    # but a hash that fails on any axis is suspect.
    lagcorr(a, b) = _cor(a, b)
    @testset "independence along each axis" begin
        @test abs(lagcorr(draws[1:(end - 1), :, :], draws[2:end, :, :])) < 0.01
        @test abs(lagcorr(draws[:, 1:(end - 1), :], draws[:, 2:end, :])) < 0.01
        @test abs(lagcorr(draws[:, :, 1:(end - 1)], draws[:, :, 2:end])) < 0.01
    end

    @testset "neighboring keys decorrelate" begin
        # Hosts pass consecutive timestep indices, so adjacent keys must not
        # produce similar samples.
        next = [
            Float64(_keyed_uniform(FT, key + UInt32(1), gcol, igpt, ilay)) for
            igpt in 1:ngpt, ilay in 1:nlay, gcol in 1:ncol
        ]
        @test abs(_cor(draws, next)) < 0.01
    end

    @testset "the two bands sample independently" begin
        # driver_utils.jl mixes a per-band constant into the key. Without it the
        # longwave and shortwave masks are identical for a given column and
        # g-point, which they must not be.
        sw_key = key ⊻ (0x9e3779b9 % UInt32)
        sw = [
            Float64(_keyed_uniform(FT, sw_key, gcol, igpt, ilay)) for
            igpt in 1:ngpt, ilay in 1:nlay, gcol in 1:ncol
        ]
        @test sw != draws
        @test abs(_cor(draws, sw)) < 0.01
    end

    @testset "the rank is mixed into the key" begin
        # Columns are indexed within a rank, so without this every subdomain
        # draws the same mask for its column i.
        salts = [_rank_salt(_FakeCtx(pid)) for pid in 1:8]
        @test allunique(salts)
        # One rank must keep the keys it had, or this changes serial results.
        @test _rank_salt(_FakeCtx(1)) == 0x00000000
        keyed = [
            Float64(_keyed_uniform(FT, key ⊻ s, 1, 1, 1)) for s in salts
        ]
        @test allunique(keyed)
    end

    @testset "seedval conversion" begin
        # `seedval` was ignored unless reset_rng_seed was set, so hosts have
        # passed floats into it; an integral one must not throw.
        @test _mcica_key(3) === _mcica_key(3.0) === 0x00000003
        @test _mcica_key(-1) === _mcica_key(-1.0)
        @test_throws ArgumentError _mcica_key(3.5)
    end

    @testset "determinism" begin
        @test _keyed_uniform(FT, key, 7, 11, 13) ==
              _keyed_uniform(FT, key, 7, 11, 13)
        @test _keyed_uniform(FT, key, 7, 11, 13) isa FT
    end
end
