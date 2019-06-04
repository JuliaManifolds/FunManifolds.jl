using FunManifolds
using Test
using StaticArrays

include("../utils.jl")

@testset "Sphere tests" begin

    s2 = Sphere(2)
    generic_manifold_tests(s2,
        [project_point_wrapped([0., 1., 0.], s2),
        project_point_wrapped([1., 0., 0.], s2),
        project_point_wrapped([0.2, 0., 1.], s2)],
        "Sphere",
        0.0)

    generic_manifold_tests(s2,
        [project_point_wrapped((@SVector [0., 1., 0.]), s2),
        project_point_wrapped((@SVector [1., 0., 0.]), s2),
        project_point_wrapped((@SVector [0.2, 0., 1.]), s2)],
        "Sphere (static)",
        0.0)
end
