using FunManifolds
using Test
using StaticArrays

include("../utils.jl")

@testset "DCurves" begin
    N = 50
    s2 = Sphere(2)
    r2 = Euclidean(2)

    dcss2 = DCurveSpace(s2, range(0.0, 1.0, length = N))
    generic_manifold_tests(dcss2,
        [discretize(t -> point2ambient(project_point_wrapped(s2, @SVector [t+0.1, t^2, sin(t)])), dcss2),
        discretize(t -> point2ambient(project_point_wrapped(s2, @SVector [t-1, sin(t^2), 1])), dcss2),
        discretize(t -> point2ambient(project_point_wrapped(s2, @SVector [t^2, t, cos(t)])), dcss2)],
        "Space of DCurves on a 2-sphere",
        0.05)

    dcsr2 = DCurveSpace(r2, range(0.0, 1.0, length = N))
    generic_manifold_tests(dcsr2,
        [discretize(t -> point2ambient(project_point_wrapped(r2, @SVector [t^2, sin(t)])), dcsr2),
        discretize(t -> point2ambient(project_point_wrapped(r2, @SVector [sin(t^2), 1])), dcsr2),
        discretize(t -> point2ambient(project_point_wrapped(r2, @SVector [t, cos(t)])), dcsr2)],
        "Space of DCurves on ℝ^2",
        0.05)

    rot = RotationAction(SpecialOrthogonalSpace(2))
    qm = QuotientManifold(r2, rot)
    dcsqm = DCurveSpace(qm, range(0.0, 1.0, length = N))
    generic_manifold_tests(dcsqm,
        [discretize(t -> point2ambient(project_point_wrapped(qm, @SVector [t^2, sin(t)])), dcsqm),
        discretize(t -> point2ambient(project_point_wrapped(qm, @SVector [sin(t^2), 1])), dcsqm),
        discretize(t -> point2ambient(project_point_wrapped(qm, @SVector [t, cos(t)])), dcsqm)],
        "Space of DCurves on ℝ^2 up to rotation",
        0.05,
        atol_vel = 1.0e-2,
        atol_velvelnorm = 5.0e-2)

    @testset "SRVF" begin
        r3 = Euclidean(3)
        dcsr1 = DCurveSpace(r3, range(0.0, 1.0, length = 9))
        c = discretize(t -> (@SVector [sin(t), t^2, 2*t+1]), dcsr1)
        q = tsrvf(c, Val(:forward_diff), ambient2point(r3, @SVector [0.0, 0.0, 0.0]))
        c_rev = reverse_srvf(q, c(0.0))
        c_sample = uniform_sample(c, 10)
        c_rev_sample = uniform_sample(c_rev, 10)
        for i in 1:length(c_sample)
            @test c_sample[i] ≈ c_rev_sample[i] atol = 1e-10
        end
    end
end
