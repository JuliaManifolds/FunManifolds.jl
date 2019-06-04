using FunManifolds
using Test
using StaticArrays

@testset "Function optimization" begin
    r2 = (EuclideanSpace(2),
        [EuclideanPt([1.0, 1.0]),
        EuclideanPt([-2.0, 3.0]),
        EuclideanPt([3.0, -2.0])])

    r3 = EuclideanSpace(3)
    f1(x) = norm(x-[1., 2., 3.])^2
    x, cond = optimize(f1, [-1., -1, -1.], r3)
    @test x ≈ EuclideanPt([1., 2., 3.])

    km = mean_karcher(r2[2])
    @test gettype(km) == gettype(r2[2][1])
    @test km ≈ mean_extrinsic(r2[2])

    r2s = (EuclideanSpace(2),
        [EuclideanPt(@SVector [1.0, 1.0]),
        EuclideanPt(@SVector [-2.0, 3.0]),
        EuclideanPt(@SVector [3.0, -2.0])])

    xs, cond = optimize(f1, (@SVector [-1., -1, -1.]), r3)
    @test x ≈ EuclideanPt(@SVector [1., 2., 3.])

    kms = mean_karcher(r2s[2])
    @test gettype(kms) == gettype(r2s[2][1])
    @test kms ≈ mean_extrinsic(r2s[2])

    sphere_pts = [
        SpherePt([1.0, 0.0, 0.0]),
        SpherePt([0.0, 1.0, 0.0]),
        SpherePt([0.0, 1.0, 0.0])
    ]

    @test mean_karcher(sphere_pts) ≈ SpherePt([0.5, sqrt(0.75), 0.0])
    @test mean_extrinsic(sphere_pts) ≈ SpherePt([1/sqrt(5.0), 2/sqrt(5.0), 0.0])
end
