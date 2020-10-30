using FunManifolds
using Test
using ForwardDiff
using StaticArrays

include("utils.jl")

@testset "FunManifolds tests" begin
    include("manifolds/DiscretizedCurves.jl")
    include("manifolds/CurveWarpingSRSFSpace.jl")
    include("manifolds/FunctionCurve.jl")
    include("curve_warping.jl")
end
