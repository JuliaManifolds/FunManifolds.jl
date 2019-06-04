using FunManifolds
using Test
using StaticArrays
using BenchmarkTools
using Statistics
using UnsafeArrays

include("utils.jl")

function mul2(array::TV) where TV<:FunManifolds.BNBArray
    FunManifolds.@condbc TV (array .= 2.0.*array)
end

function test_condbc(array::FunManifolds.BNBArray, allowed_time)
    array_ref = deepcopy(array)
    mul2(array)
    if size(array) == ()
        @test array[] ≈ 2.0 * array_ref
    else
        @test array ≈ 2.0 * array_ref
    end
    bench_results = @benchmark mul2($array) seconds=0.5
    #@test median(bench_results.times) < allowed_time
end

@testset "TupleArray tests" begin
    a = [1.0, 2.0, 3.0]
    b = [@SVector [1.0, 2.0, 3.0]]
    c = @MVector [1.0, 2.0, 3.0]
    d = [1.0 2.0 3.0; 4.0 5.0 6.0]
    dm = @MMatrix [1.0 2.0 3.0; 4.0 5.0 6.0]
    e = [FunManifolds.TupleArray((1.0, 2.0, 3.0))]
    em = @MVector [FunManifolds.TupleArray((1.0, 2.0, 3.0))]
    f = [FunManifolds.TupleArray(((@SVector [1.0, 2.0, 3.0]), (@SVector [3.0, 5.0])))]
    fm = @MVector [FunManifolds.TupleArray(((@SVector [1.0, 2.0, 3.0]), (@SVector [3.0, 5.0])))]
    # storing non-isbits types in TupleArray is *not* a good idea
    # for performance reasons and therefore it's not a fully supported mode
    # of operation
    g = [FunManifolds.TupleArray(((@MVector [1.0, 2.0, 3.0]), (@MVector [3.0, 5.0])))]
    test_condbc(a, 10)
    test_condbc(view(b, 1), 10)
    test_condbc(uview(b, 1), 10)
    test_condbc(c, 10)
    test_condbc(view(d, 1, :), 20)
    test_condbc(uview(d, 1, :), 20)
    test_condbc(view(dm, 1, :), 15)
    test_condbc(uview(dm, 1, :), 15)
    test_condbc(view(e, 1), 10)
    test_condbc(uview(e, 1), 10)
    test_condbc(view(em, 1), 10)
    test_condbc(uview(em, 1), 10)
    test_condbc(view(f, 1), 10)
    test_condbc(uview(f, 1), 10)
    test_condbc(view(fm, 1), 10)
    test_condbc(uview(fm, 1), 10)
    # these two tests are very slow and represent operations that should not
    # be performed in practice (but they do work)
    test_condbc(view(g, 1), 500)
    test_condbc(uview(g, 1), 500)
    test_condbc(g[1], 40)

    g_sim = similar(g[1])
    @test typeof(g_sim) == typeof(g[1])

    @test deepcopy(g[1]) ≈ g[1]

    tup_a = FunManifolds.TupleArray(((@MVector [1.0, 2.0, 3.0]), (@MVector [3.0, 5.0])))
    tup_b = FunManifolds.TupleArray(((@MVector [0.0, 2.0, 3.0]), (@MVector [3.0, 5.0])))
    FunManifolds.@condbc typeof(tup_a) (tup_a .+= tup_b) (tup_b,)
    @test tup_a[1][1] ≈ 1.0
    @test tup_a[1][2] ≈ 4.0
end
