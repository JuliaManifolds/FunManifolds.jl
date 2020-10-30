include("utils.jl")

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


@testset "Curve SRSF warping action" begin
    N = 100
    r2 = Euclidean(2)
    knots = range(0.0, 1.0, length = N)
    dcsr2 = DCurves(r2, knots)
    M = CurveWarpingSRSFSpace(knots)
    M2 = CurveWarpingSpace(knots)
    G = CurveWarpingSRSFGroup(M)
    p = [0.0, 0.0]
    dcsr2t = DCurves(TangentSpaceAtPoint(r2, p), knots)
    A_left = CurveWarpingSRSFAction(dcsr2t, p, G)

    #@test repr(A_left) == "CurveWarpingSRSFAction($(repr(dcsr2)), $(repr(G)))"

    @test g_manifold(A_left) == dcsr2t
    @test base_group(A_left) == G
    @test isa(A_left, AbstractGroupAction{LeftAction})
    @test base_manifold(G) == M

    a_funs = [x -> x^2, x -> x^3, x -> sqrt(x)]
    a_pts = [srsf(M2, FunManifolds.make_warping(M2, map(f, knots))) for f in a_funs]
    m_pts = [
        tsrvf(
            dcsr2,
            discretize(dcsr2, t -> project(r2, @SVector [t^2, sin(t)])),
            p,
            FunManifolds.ProjectedDifferenceBackend(nothing),
        ),
        tsrvf(
            dcsr2,
            discretize(dcsr2, t -> project(r2, @SVector [sin(t^2), 1])),
            p,
            FunManifolds.ProjectedDifferenceBackend(nothing),
        ),
        tsrvf(
            dcsr2,
            discretize(dcsr2, t -> project(r2, @SVector [t, cos(t)])),
            p,
            FunManifolds.ProjectedDifferenceBackend(nothing),
        ),
    ]

    test_action(
        A_left,
        a_pts,
        m_pts;
        atol = 3 / sqrt(N),
        atol_ident_compose = 2 * sqrt(eps(Float64)),
        test_optimal_alignment = true,
        test_diff = false,
        test_mutating = false,
        test_switch_direction = false,
    )
end

@testset "Mean amplitude" begin
    grid_small = 0.0:0.125:1.0
    r1 = Euclidean(1)
    r2 = Euclidean(2)
    f1 = t -> begin
        [if t < 0.5
            0.0
        elseif t < 0.75
            t - 0.5
        else
            1.0 - t
        end]
    end
    f2 = t -> begin
        [if t < 0.25
            0.0
        elseif t < 0.5
            t - 0.25
        elseif t < 0.75
            0.75 - t
        else
            0.0
        end]
    end
    f3 = t -> begin
        [if t < 0.25
            0.0
        elseif t < 0.5
            -20 * (t - 0.25)
        elseif t < 0.75
            -20 * (0.75 - t)
        else
            0.0
        end]
    end
    dcsr1 = DCurves(r1, grid_small)
    M = CurveWarpingSRSFSpace(grid_small)
    M2 = CurveWarpingSpace(grid_small)
    G = CurveWarpingSRSFGroup(M)
    p = [0.0]
    dcsr1t = DCurves(TangentSpaceAtPoint(r1, p), grid_small)
    A_left = CurveWarpingSRSFAction(dcsr1t, p, G)

    f1qc = tsrvf(
        dcsr1,
        discretize(dcsr1, f1),
        p,
        FunManifolds.ProjectedDifferenceBackend(nothing),
    )
    f2qc = tsrvf(
        dcsr1,
        discretize(dcsr1, f2),
        p,
        FunManifolds.ProjectedDifferenceBackend(nothing),
    )
    f3qc = tsrvf(
        dcsr1,
        discretize(dcsr1, f3),
        p,
        FunManifolds.ProjectedDifferenceBackend(nothing),
    )
    μq = FunManifolds.karcher_mean_amplitude(A_left, [f1qc, f2qc])
    for pt in μq
        @test abs(pt) < 2.0
    end
    (p̃, γs, aligned_ps) = FunManifolds.phase_amplitude_separation(A_left, [f1qc, f2qc])
    @test is_manifold_point(dcsr1t, p̃)
    for γ in γs
        @test is_manifold_point(M, γ)
    end
    for p in aligned_ps
        @test is_manifold_point(dcsr1t, p)
    end
end
