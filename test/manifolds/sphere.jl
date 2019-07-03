using FunManifolds
using Test
using StaticArrays

include("../utils.jl")

@testset "Sphere tests" begin

    s2 = Sphere(2)
    generic_manifold_tests(s2,
        [project_point_wrapped(s2, [0., 1., 0.]),
        project_point_wrapped(s2, [1., 0., 0.]),
        project_point_wrapped(s2, [0.2, 0., 1.])],
        "Sphere",
        0.0)

    generic_manifold_tests(s2,
        [project_point_wrapped(s2, @SVector [0., 1., 0.]),
        project_point_wrapped(s2, @SVector [1., 0., 0.]),
        project_point_wrapped(s2, @SVector [0.2, 0., 1.])],
        "Sphere (static)",
        0.0)
end
