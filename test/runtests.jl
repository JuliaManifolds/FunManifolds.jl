using FunManifolds
using Test
using ForwardDiff
using StaticArrays

include("utils.jl")

@testset "FunManifolds tests" begin
    include("Manifolds/DCurves.jl")
    include("manifolds/FunctionCurve.jl")
end
