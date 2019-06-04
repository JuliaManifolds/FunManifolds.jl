using FunManifolds
using Test
using StaticArrays

include("../utils.jl")

@testset "Product space" begin

    s2 = Sphere(2)

    r2 = (EuclideanSpace(2),
        [EuclideanPt([1.0, 1.0]),
        EuclideanPt([-2.0, 3.0]),
        EuclideanPt([3.0, -2.0])])

    sphere = (s2,
        [project_point_wrapped([0., 1., 0.], s2),
        project_point_wrapped([1., 0., 0.], s2),
        project_point_wrapped([0.2, 0., 1.], s2)])

    generic_manifold_tests(ProductSpace((r2[1], sphere[1])),
        [ProductPt((r2[2][i], sphere[2][i])) for i ∈ 1:3],
        "Product space",
        0.0)

    r2_s = (EuclideanSpace(2),
        [EuclideanPt(@SVector [1.0, 1.0]),
        EuclideanPt(@SVector [-2.0, 3.0]),
        EuclideanPt(@SVector [3.0, -2.0])])

    sphere_s = (s2,
        [project_point_wrapped((@SVector [0., 1., 0.]), s2),
        project_point_wrapped((@SVector [1., 0., 0.]), s2),
        project_point_wrapped((@SVector [0.2, 0., 1.]), s2)])

    generic_manifold_tests(ProductSpace((r2_s[1], sphere_s[1])),
        [ProductPt((r2_s[2][i], sphere_s[2][i])) for i ∈ 1:3],
        "Product space (static)",
        0.0)
end
