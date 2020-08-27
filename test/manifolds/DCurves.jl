using FunManifolds
using Test
using StaticArrays

include("../utils.jl")

using Manifolds: get_iterator, _write, _read

function discretize(M::DCurves, f)
    rep_size = representation_size(M.manifold)
    p = similar(f(M.grid[1]), rep_size..., length(M.grid))
    for i in get_iterator(M)
        copyto!(_write(M, rep_size, p, i), f(M.grid[i]))
    end
    return p
end

@testset "DCurves" begin
    N = 50
    s2 = Sphere(2)
    r2 = Euclidean(2)

    dcss2 = DCurves(s2, range(0.0, 1.0, length = N))
    test_manifold(
        dcss2,
        [
            discretize(dcss2, t -> project(s2, @SVector [t + 0.1, t^2, sin(t)])),
            discretize(dcss2, t -> project(s2, @SVector [t - 1, sin(t^2), 1])),
            discretize(dcss2, t -> project(s2, @SVector [t^2, t, cos(t)])),
        ];
        is_tangent_atol_multiplier = 10,
        test_default_vector_transport = true,
    )

    dcsr2 = DCurves(r2, range(0.0, 1.0, length = N))
    test_manifold(
        dcsr2,
        [
            discretize(dcsr2, t -> project(r2, @SVector [t^2, sin(t)])),
            discretize(dcsr2, t -> project(r2, @SVector [sin(t^2), 1])),
            discretize(dcsr2, t -> project(r2, @SVector [t, cos(t)])),
        ];
        is_tangent_atol_multiplier = 10,
        test_default_vector_transport = true,
    )

    @testset "SRVF" begin
        r3 = Euclidean(3)
        dcsr1 = DCurves(r3, range(0.0, 1.0, length = 9))
        f = t -> (@SVector [sin(t), t^2, 2 * t + 1])
        c = discretize(dcsr1, f)
        q = tsrvf(dcsr1, c, (@SVector [0.0, 0.0, 0.0]), FunManifolds.ProjectedDifferenceBackend(nothing))
        c_rev = reverse_srvf(dcsr1, q, f(0.0))
        @test c ≈ c_rev atol = 1e-10
    end
end
