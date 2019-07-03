using FunManifolds
using Test
using StaticArrays

include("../utils.jl")

@testset "Tangent manifolds" begin
    s2 = Sphere(2)
    p_sphere = project_point_wrapped(s2, [0., 1., 0.])
    r2 = EuclideanSpace(2)
    p_euclidean = ambient2point(r2, [1.0, -1.0])
    p_euclidean_static = ambient2point(r2, @SVector [1.0, -1.0])

    generic_manifold_tests(TSpaceManifold(p_euclidean),
        [TSpaceManifoldPt(EuclideanTV(p_euclidean, [0.5, 0.])),
        TSpaceManifoldPt(EuclideanTV(p_euclidean, [0., 0.5])),
        TSpaceManifoldPt(EuclideanTV(p_euclidean, [0.5, 0.5]))],
        "Tangent space manifold (euclidean)",
        0.0)

    generic_manifold_tests(TSpaceManifold(p_sphere),
        [TSpaceManifoldPt(project_tv([0.5, 0., 0.], p_sphere)),
        TSpaceManifoldPt(project_tv([0., 0., 0.5], p_sphere)),
        TSpaceManifoldPt(project_tv([0.5, 0., 0.5], p_sphere))],
        "Tangent space manifold (sphere)",
        0.0)

    generic_manifold_tests(TSpaceManifold(p_euclidean_static),
        [TSpaceManifoldPt(EuclideanTV(p_euclidean_static, @MVector [0.5, 0.])),
        TSpaceManifoldPt(EuclideanTV(p_euclidean_static, @MVector [0., 0.5])),
        TSpaceManifoldPt(EuclideanTV(p_euclidean_static, @MVector [0.5, 0.5]))],
        "Tangent space manifold (euclidean) (static)",
        0.0)

    @testset "Other tests" begin
        p = project_point_wrapped(Sphere(2), [0., 1., 0.])
        tv = SphereTV(p, [0.5, 0., 0.])
        y = exp(tv)
        tvp = TSpaceManifoldPt(tv)
        @test manifold_dimension(gettype(tvp)) == 2
        ztv = zero_tangent_vector(tvp)
        ttvp = TSpaceManifoldPt(ztv)
        @test manifold_dimension(gettype(ttvp)) == 2
    end
end
