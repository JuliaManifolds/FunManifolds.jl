include("../utils.jl")

@testset "FunctionCurve manifold" begin
    M = FunctionCurveSpace(Euclidean(2))
    @testset "Oblique manifold Basics" begin
        @test manifold_dimension(M) === Inf
    end

    x = t -> SA[t, t - 1]
    y = t -> SA[t^2, sin(t)]
    z = t -> SA[2 - 2, cos(t)]

    test_manifold(
        M,
        [x, y, z],
        test_forward_diff = false,
        test_reverse_diff = false,
        test_vector_spaces = true,
        test_project_tangent = false,
        test_musical_isomorphisms = false,
        test_vector_transport = false,
        test_tangent_vector_broadcasting = false,
        test_representation_size = false,
        is_mutating = false,
    )
end
