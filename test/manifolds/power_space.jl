using FunManifolds
using Test

include("../utils.jl")

@testset "Power space" begin
    s2 = Sphere(2)

    sphere = (s2,
        [project_point_wrapped([0., 1., 0.], s2),
        project_point_wrapped([1., 0., 0.], s2),
        project_point_wrapped([0.2, 0., 1.], s2)])

    generic_manifold_tests(PowerSpace(sphere[1], 3),
        [PowerPt([sphere[2][mod1(i+j, 3)] for j ∈ 1:3]) for i ∈ 1:3],
        "Power space",
        0.0)
end
