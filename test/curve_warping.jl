include("utils.jl")

using FunManifolds: CurveWarpingGroup, CurveWarpingSpace, CurveWarpingAction

@testset "Curve warping action" begin
    N = 20
    r2 = Euclidean(2)
    knots = range(0.0, 1.0, length = N)
    dcsr2 = DCurves(r2, knots)
    M = CurveWarpingSpace(knots)
    G = CurveWarpingGroup(M)
    A_left = CurveWarpingAction(dcsr2, G)

    #@test repr(A_left) == "CurveWarpingAction($(repr(dcsr2)), $(repr(G)))"

    @test g_manifold(A_left) == dcsr2
    @test base_group(A_left) == G
    @test isa(A_left, AbstractGroupAction{LeftAction})
    @test base_manifold(G) == M

    a_funs = [x -> x^2, x -> x^3, x -> sqrt(x)]
    a_pts = [FunManifolds.make_warping(M, map(f, knots)) for f in a_funs]
    m_pts = [
        discretize(dcsr2, t -> project(r2, @SVector [t^2, sin(t)])),
        discretize(dcsr2, t -> project(r2, @SVector [sin(t^2), 1])),
        discretize(dcsr2, t -> project(r2, @SVector [t, cos(t)])),
    ]

    test_action(
        A_left,
        a_pts,
        m_pts;
        atol = 2 / N,
        test_optimal_alignment = false,
        test_diff = false,
        test_mutating = false,
        test_switch_direction = false,
    )
end
