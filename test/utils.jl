const TEST_FLOAT32 = false
const TEST_DOUBLE64 = false
const TEST_STATIC_SIZED = false

using Manifolds
using ManifoldsBase
using ManifoldsBase: number_of_coordinates
using FunManifolds

using LinearAlgebra
using Distributions
using DoubleFloats
using ForwardDiff
using Random
using ReverseDiff
using StaticArrays
using Test

using Manifolds.ManifoldTests
using Manifolds.ManifoldTests: find_eps

# local methods
function Manifolds.ManifoldTests.find_eps(f1::Function, fs::Function...)
    return find_eps(f1(0.0), map(g -> g(0.0), fs)...)
end
function Manifolds.ManifoldTests.find_eps(
    f1::FunManifolds.VectorizedFunction,
    fs::FunManifolds.VectorizedFunction...,
)
    return find_eps(f1(0.0), map(g -> g(0.0), fs)...)
end
