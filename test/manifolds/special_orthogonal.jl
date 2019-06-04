using FunManifolds
using StaticArrays
using Test

include("../utils.jl")

@testset "Special orthogonal space" begin
    generic_manifold_tests(SpecialOrthogonalSpace(3),
        [rotation3d_from_yaw_pitch_roll(.0, π/2, .0),
        rotation3d_from_yaw_pitch_roll(π/2, π/2, .0),
        rotation3d_from_yaw_pitch_roll(.0, .0, 2.)],
        "Special orthogonal group",
        0.0)

    generic_manifold_tests(SpecialOrthogonalSpace(3),
        [rotation3d_from_yaw_pitch_roll_s(.0, π/2, .0),
        rotation3d_from_yaw_pitch_roll_s(π/2, π/2, .0),
        rotation3d_from_yaw_pitch_roll_s(.0, .0, 2.)],
        "Special orthogonal group (static)",
        0.0)
end
