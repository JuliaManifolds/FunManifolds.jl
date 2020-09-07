using FunManifolds
using Test
using StaticArrays

include("../utils.jl")

@testset "CurveWarpingSRSFSpace" begin
    for N in [21, 50]
        nodes = range(0.0, 1.0, length = N)
        M = CurveWarpingSpace(nodes)
        MS = CurveWarpingSRSFSpace(nodes)

        test_manifold(
            MS,
            [srsf(M, t -> t^2), srsf(M, t -> t^3), srsf(M, t -> t^(0.7))];
            exp_log_atol_multiplier = N % 2 == 1 ? 10 : 100,
            test_default_vector_transport = true,
            test_reverse_diff = false,
            test_forward_diff = false,
            test_injectivity_radius = false,
        )
    end

    @testset "reverse SRSF" begin
        N = 11
        nodes = range(0.0, 1.0, length = N)
        M = CurveWarpingSpace(nodes)
        MS = CurveWarpingSRSFSpace(nodes)

        γ = FunManifolds.make_warping(M, [t^2 for t in nodes])
        γSRSF = srsf(M, γ)
        @test isapprox(M, γ, reverse_srsf(M, γSRSF), atol = 1e-15)
    end
end
