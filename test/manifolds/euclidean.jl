using FunManifolds
using StaticArrays
using Test

include("../utils.jl")

@testset "Euclidean manifolds tests" begin
    # generic tests
    generic_manifold_tests(
        EuclideanSpace(2),
        [EuclideanPt([1.0, 1.0]),
        EuclideanPt([-2.0, 3.0]),
        EuclideanPt([3.0, -2.0])],
        "Euclidean space",
        0.0)

    generic_manifold_tests(EuclideanSpace(2),
        [EuclideanPt(@SVector [1.0, 1.0]),
        EuclideanPt(@SVector [-2.0, 3.0]),
        EuclideanPt(@SVector [3.0, -2.0])],
        "Euclidean space (static)",
        0.0)

    p1 = EuclideanPt([1.0, 1.0])
    p2 = EuclideanPt([-2.0, 3.0])
    v1 = EuclideanTV(p1, [2.0, 3.0])
    v2 = EuclideanTV(p1, [4.0, -1.0])
    @test norm(v1) ≈ sqrt(13.0)
    @test innerproduct_amb(p1, p2) ≈ 1.0
    vs = v1 + v2
    v2mul = 2.0 * v2
    @test vs.v ≈ v1.v + v2.v
    @test v2mul.v ≈ 2.0 * v2.v
    p2 = expmap(v1)
    @test p2.x ≈ [3.0, 4.0]

    parv1 = parallel_transport_geodesic(v1, p2)
    @test parv1.v ≈ v1.v

    #curve differentiation
    testF = t::Real -> EuclideanPt([t, t*t])
    c = CurvePt(testF, EuclideanSpace(2))
    cV = velocity(c, Val(:continuous))
    @test at_point(cV(1.0).x) ≈ EuclideanPt([1.0, 1.0])
    @test cV(1.0).x.v ≈ [1.0, 2.0]
end
