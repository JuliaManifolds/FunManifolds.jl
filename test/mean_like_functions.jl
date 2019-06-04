using FunManifolds
using Test

@testset "Mean-like functions" begin

    pts_eucl = [EuclideanPt([0.0, 0.0]),
        EuclideanPt([0.3, 2.0]),
        EuclideanPt([-34.0, -2.0])]

    @test mean_karcher(pts_eucl) ≈ mean_extrinsic(pts_eucl)
end
