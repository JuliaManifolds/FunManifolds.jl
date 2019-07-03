using FunManifolds
using Test
using StaticArrays

function generic_manifold_tests(space::Manifold, pts, name::String, atol::Real;
    atol_g1 = nothing, atol_vel = 1.e-10, atol_velvelnorm = 1.e-8,
    atol_innerprod = atol, atol_project_tv = 1.0e-15, check_space_equality = true,
    cds = false)
    @testset "Generic testset for $name" begin
        if check_space_equality
            for p in pts
                @test gettype(p) == space
            end
        end

        tv1 = log(pts[1], pts[2])
        tv2 = log(pts[1], pts[3])
        if atol > 0.0
            @test exp(tv1) ≈ pts[2] atol = atol
        else
            @test exp(tv1) ≈ pts[2]
        end
        @test dim(tv1) == dim(pts[1])
        @test dim(space) == dim(pts[2])
        PARAMS.quad_abs_tol = 1e-6
        @test norm(tv1) >= 0
        @test norm(tv1) ≈ geodesic_distance(pts[1], pts[2]) atol = 2*atol+1e-6
        @test norm(tv2) ≈ geodesic_distance(pts[1], pts[3]) atol = 2*atol+1e-6
        @test norm(tv1) ≈ sqrt(inner(tv1, tv1))
        @test norm(tv2) ≈ sqrt(inner(tv2, tv2))
        PARAMS.quad_abs_tol = nothing
        if dim(tv1) < Inf
            @test norm(tv1) >= 0
        end
        if dim(tv1) < Inf && name != "Space of discretized curves on a 2-sphere"
            tv3 = log(pts[1], pts[2])
            tv4 = log(pts[1], pts[3])
            tv5 = deepcopy(tv3)
            @test tv5 ≈ tv3
            mul_vec!(tv3, 2.0)
            @test tv3 ≈ 2.0 * tv1
            add_vec!(tv3, tv4)
            @test tv3 ≈ 2.0 * tv1 + tv2
            sub_vec!(tv3, tv4)
            @test tv3 ≈ 2.0 * tv1
            @test tv5 ≉ tv3
        end
        @test tv1 + tv1 ≈ 2.0 * tv1
        @test tv1 - tv1 ≈ zero_tv(pts[1])
        @test inner_amb(pts[1], pts[1]) ≥ 0.0

        PARAMS.quad_abs_tol = 1e-6
        PARAMS.quad_rel_tol = 1e-6
        @test riemannian_distortion(pts) ≥ 0.0
        PARAMS.quad_abs_tol = nothing
        PARAMS.quad_rel_tol = nothing

        tv1ptg = parallel_transport_geodesic(tv1, pts[2])
        tv2ptg = parallel_transport_geodesic(tv2, pts[2])
        #println(name, " | ", inner(tv1, tv1, quad_abs_tol = 1e-6),
        #    " | ", inner(tv1ptg, tv1ptg, quad_abs_tol = 1e-6))
        #println(name, " | ", inner(tv2, tv2, quad_abs_tol = 1e-6),
        #    " | ", inner(tv2ptg, tv2ptg, quad_abs_tol = 1e-6))
        #println(name, " | ", inner(tv1, tv2, quad_abs_tol = 1e-6),
        #    " | ", inner(tv1ptg, tv2ptg, quad_abs_tol = 1e-6))
        PARAMS.quad_abs_tol = 1e-6
        @test inner(tv1, tv2) ≈ inner(tv1ptg, tv2ptg) atol = atol_innerprod
        PARAMS.quad_abs_tol = nothing
        if dim(space) < Inf
            @test dim_ambient(space) == prod(ambient_shape(space))
            for i ∈ 1:3
                @test ambient2point(point2ambient(pts[i]), space) ≈ pts[i]
                @test project_point_wrapped(point2ambient(pts[i]), space) ≈ pts[i]
            end
            amb_some = point2ambient(pts[1]) + point2ambient(pts[2])
            amb_proj = project_point_wrapped(deepcopy(amb_some), space)
            #if amb_some isa SArray
            if isbits(amb_some)
                amb_wrapped = @MVector [amb_some]
                project_point!(view(amb_wrapped, 1), space)
                @test amb_proj ≈ ambient2point(amb_wrapped[1], space)
            else
                project_point!(amb_some, space)
                @test amb_proj ≈ ambient2point(amb_some, space)
            end
            @test ambient2tangent(tangent2ambient(tv1), at_point(tv1)) ≈ tv1 atol = 1.e-15
            @test ambient2tangent(tangent2ambient(tv2), at_point(tv2)) ≈ tv2 atol = 1.e-15
            @test project_tv(tangent2ambient(tv1), at_point(tv1)) ≈ tv1 atol = atol_project_tv
            @test project_tv(tangent2ambient(tv2), at_point(tv2)) ≈ tv2 atol = atol_project_tv
            if !cds
                # CurveDiscretizedSpace is deprecated and we don't really need these
                @test inner(tv1, tv2) ≈ inner(tangent2ambient(tv1), tangent2ambient(tv2), point2ambient(pts[1]), space)
                @test geodesic_distance(pts[1], pts[2]) ≈ geodesic_distance(point2ambient(pts[1]), point2ambient(pts[2]), space)
                @test inner(tv1, tv2) ≈ inner(tangent2ambient(tv1), tangent2ambient(tv2), point2ambient(at_point(tv1)), space)
            end
            #testing modifying actions
            if !cds
                tv1p = tangent2ambient(tv1) + tangent2ambient(tv2)
                @test typeof(tv1p) === typeof(tangent2ambient(tv1))
                #println(typeof(tv1p))
                project_tv!(tv1p, at_point(tv1))
                @test ambient2tangent(tv1p, at_point(tv1)) ≈ project_tv(tangent2ambient(tv1) + tangent2ambient(tv2), at_point(tv1)) atol=1.e-15

                log!(tv1p, point2ambient(pts[1]), point2ambient(pts[2]), space)
                @test tv1p ≈ tangent2ambient(tv1)

                parallel_transport_geodesic!(tv1p, tangent2ambient(tv1), point2ambient(at_point(tv1)), point2ambient(pts[2]), space)
                @test tv1p ≈ tangent2ambient(tv1ptg)

                ptest = FunManifolds._ensure_mutable(point2ambient(exp(tv2)))
                exp!(ptest, tangent2ambient(tv1), point2ambient(at_point(tv1)), space)
                @test ptest ≈ point2ambient(exp(tv1))

                zero_tv!(tv1p, point2ambient(pts[1]), space)
                @test tv1p ≈ tangent2ambient(zero_tv(pts[1]))

                tv1p = deepcopy(tangent2ambient(tv1))
                mul_vec!(tv1p, 2.0, point2ambient(pts[1]), space)
                @test tv1p ≈ tangent2ambient(2.0*tv1)

                tv1p = deepcopy(tangent2ambient(tv1))
                add_vec!(tv1p, tv1p, point2ambient(pts[1]), space)
                @test tv1p ≈ tangent2ambient(2.0*tv1)

                sub_vec!(tv1p, tangent2ambient(tv1), point2ambient(pts[1]), space)
                @test tv1p ≈ tangent2ambient(tv1)
            end

            #perfomance tests
            @test begin; (@inferred point2ambient(pts[1])); true; end
            @test begin; (@inferred ambient2point(point2ambient(pts[1]), space)); true; end
            @test begin; (@inferred project_point(point2ambient(pts[1]), space)); true; end
            @test begin; (@inferred project_point_wrapped(point2ambient(pts[1]), space)); true; end
            @test begin; (@inferred tangent2ambient(tv1)); true; end
            @test begin; (@inferred ambient2tangent(tangent2ambient(tv1), at_point(tv1))); true; end
            @test begin; (@inferred project_tv(tangent2ambient(tv1), at_point(tv1))); true; end
        end
        g = geodesic(pts[1], pts[2])
        @test g(0.0) ≈ pts[1]
        @test g(0.3) ≈ geodesic_at(0.3, pts[1], pts[2])
        if atol_g1 === nothing
            @test g(1.0) ≈ pts[2]
        else
            @test g(1.0) ≈ pts[2] atol = atol_g1
        end
        if dim(space) < Inf
            @test ambient2point(geodesic_at(0.3, point2ambient(pts[1]), point2ambient(pts[2]), space), space) ≈ g(0.3)

            gv = velocity(g, Val(:continuous))
            gvv = velocity(gv, Val(:continuous))
            vels = [norm(s.x) for s in uniform_sample(gv, 100)]
            velvelnorms = [norm(x.x.v_ts.v) for x in uniform_sample(gvv, 10)]
            # println(name, " | velvelnorms: ", velvelnorms)
            # println(name, " | vels: ", vels)
            @test all(isapprox(vels[2], vi, atol = atol_vel) for vi in vels[2:end])
            @test all(isapprox(x, 0.0, atol = atol_velvelnorm) for x in velvelnorms)
            #println(name, ": ", curve_length(g) - geodesic_distance(pts[1], pts[2]))
            PARAMS.quad_abs_tol = atol/2.0
            @test curve_length(g) ≈ geodesic_distance(pts[1], pts[2]) atol = atol
            PARAMS.quad_abs_tol = nothing

            # performance tests
            @inferred velocity(g, Val(:continuous))
        end
        PARAMS.quad_abs_tol = 1e-6
        geod_dist12 = geodesic_distance(pts[1], pts[2])
        geod_dist23 = geodesic_distance(pts[2], pts[3])
        geod_dist13 = geodesic_distance(pts[1], pts[3])
        @test geod_dist12 + geod_dist23 + 3e-6 ≥ geod_dist13

        @test ambient_distance(pts[1], pts[2]) ≥ 0.0
        @test ambient_distance(pts[2], pts[3]) ≥ 0.0
        @test ambient_distance(pts[1], pts[3]) ≥ 0.0
        PARAMS.quad_abs_tol = nothing

        # general performance tests
        if dim(space) < Inf
            # type stability of these functions is less critical this legacy manifold
            if !cds
                @test begin; (@inferred project_point(point2ambient(pts[1]), space)); true; end
                @test begin; (@inferred project_point_wrapped(point2ambient(pts[1]), space)); true; end
                @test begin; (@inferred project_tv(tangent2ambient(tv1), at_point(tv1))); true; end
                @test begin; (@inferred zero_tv(pts[1])); true; end
                @test begin; (@inferred log(pts[1], pts[2])); true; end
                @test begin; (@inferred exp(tv1)); true; end
                @test begin; (@inferred norm(tv1)); true; end
            end
        end
    end
end
