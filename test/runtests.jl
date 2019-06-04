using FunManifolds
using Test
using ForwardDiff
using StaticArrays

include("utils.jl")

@testset "FunManifolds tests" begin

    include("tuple_array.jl")

    include("manifolds/euclidean.jl")
    include("manifolds/sphere.jl")
    include("manifolds/product_space.jl")
    include("manifolds/power_space.jl")
    include("manifolds/special_orthogonal.jl")
    include("manifolds/tangent_manifold.jl")
    include("manifolds/tangent_bundle.jl")

    #include("optimization.jl")
    #include("mean_like_functions.jl")

end
