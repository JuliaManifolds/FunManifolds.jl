using FunManifolds
using Test
using StaticArrays

include("../utils.jl")

@testset "DiscretizedCurves" begin
    N = 50
    s2 = Sphere(2)
    r2 = Euclidean(2)

    dcss2 = DiscretizedCurves(s2, range(0.0, 1.0, length=N))
    test_manifold(
        dcss2,
        [
            discretize(dcss2, t -> project(s2, @SVector [t + 0.1, t^2, sin(t)])),
            discretize(dcss2, t -> project(s2, @SVector [t - 1, sin(t^2), 1])),
            discretize(dcss2, t -> project(s2, @SVector [t^2, t, cos(t)])),
        ];
        exp_log_atol_multiplier=10.0,
        is_tangent_atol_multiplier=10,
        test_default_vector_transport=true,
        test_project_point=true,
        test_project_tangent=true,
    )

    dcsr2 = DiscretizedCurves(r2, range(0.0, 1.0, length=N))
    test_manifold(
        dcsr2,
        [
            discretize(dcsr2, t -> project(r2, @SVector [t^2, sin(t)])),
            discretize(dcsr2, t -> project(r2, @SVector [sin(t^2), 1])),
            discretize(dcsr2, t -> project(r2, @SVector [t, cos(t)])),
        ];
        is_tangent_atol_multiplier=10,
        test_default_vector_transport=true,
    )

    @testset "SRVF" begin
        r3 = Euclidean(3)
        dcsr1 = DiscretizedCurves(r3, range(0.0, 1.0, length=9))
        f = t -> (@SVector [sin(t), t^2, 2 * t + 1])
        c = discretize(dcsr1, f)
        q = tsrvf(
            dcsr1,
            c,
            (@SVector [0.0, 0.0, 0.0]),
            FunManifolds.ProjectedDifferenceBackend(nothing),
        )
        c_rev = reverse_srvf(dcsr1, q, f(0.0))
        @test c ≈ c_rev atol = 1e-10
    end
end
