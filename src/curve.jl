
"""
    CurveSpace(m, grid)

Space of curves in manifold `m`. When needed for calculations, uses grid `grid`
but does not discretize the curve when not necessary.
They have no embedding in a finite-dimensional euclidean space
(use `DCurveSpace` if you need this feature).
"""
struct CurveSpace{M <: Manifold, R <: AbstractRange} <: Manifold
    manifold_type::M
    approx_grid::R
    CurveSpace{M, R}(manifold_type::M, r::R = 0.0:0.01:1.0) where {M <: Manifold, R <: AbstractRange} = new(manifold_type, r)
end

function values_in(cs::CurveSpace)
    return cs.manifold_type
end

function CurveSpace(manifold_type::M, r::R = 0.0:0.01:1.0) where {M <: Manifold, R <: AbstractRange}
    return CurveSpace{M, R}(manifold_type, r)
end

function dim(x::CurveSpace)
    return Inf
end

"""
    CurvePt(f, cs)

A continuous curve on a manifold from curve space `cs` defined by function `f`.
It has no embedding in a finite-dimensional euclidean space
(use `DCurvePt` if you need this feature).
"""
struct CurvePt{M <: Manifold, F} <: AbstractCurvePt
    f::F
    pt_type::CurveSpace{M}
end

CurvePt(f::F, m::M) where {F,M<:Manifold} = CurvePt{M, F}(f, CurveSpace(m))

function (fun::CurvePt)(x::Real)
    return fun.f(x)
end

function isapprox(v1::CurvePt, v2::CurvePt; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    for i in v1.pt_type.approx_grid
        if !(isapprox(v1(i), v2(i), atol = atol, rtol = rtol))
            println("*** $i, $(v1(i)), $(v2(i)).")
        end
    end
    return all(isapprox(v1(i), v2(i), atol = atol, rtol = rtol) for i ∈ v1.pt_type.approx_grid)
end

function gettype(x::CurvePt)
    return x.pt_type
end

function uniform_sample(c::CurvePt, N::Int)
    return [c(t) for t = 0.0:1/N:1.0]
end


@doc doc"
    CurveTV(x, v)

Tangent vector for a curve space over manifold $M$ in tangent space
at point `x`, represented by curve `v` in tangent bundle of $M$.
"
struct CurveTV{M <: Manifold, F1, F2} <: TangentVector
    pt::CurvePt{M, F1}
    v::CurvePt{TangentBundleSpace{M}, F2}
end

function at_point(v::CurveTV)
    return v.pt
end

function zero_tv(pt::CurvePt)
    curve_tv = CurvePt(TangentBundleSpace(pt.pt_type.manifold_type)) do t
        return TangentBundlePt(zero_tv(pt(t)))
    end
    return CurveTV(pt, curve_tv)
end

function isapprox(v1::CurveTV, v2::CurveTV; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    if !(isapprox(v1.pt, v2.pt, atol = atol, rtol = rtol))
        return false
    end
    #TODO add tolerance parameters to isapprox for this quadrature?
    I, err = quadgk(t -> norm(v1.v(t).x - v2.v(t).x), 0.0, 1.0)
    return I <= 2*err+atol && I >= -2*err-atol
end

function +(v1::CurveTV, v2::CurveTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return CurveTV(at_point(v1),
        CurvePt(TangentBundleSpace(v1.pt.pt_type.manifold_type)) do t
            return TangentBundlePt(v1.v(t).x + v2.v(t).x)
        end)
end

function -(v1::CurveTV, v2::CurveTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return CurveTV(at_point(v1),
        CurvePt(TangentBundleSpace(v1.pt.pt_type.manifold_type)) do t
            return TangentBundlePt(v1.v(t).x - v2.v(t).x)
        end)
end

function *(α::Real, v::CurveTV)
    return CurveTV(at_point(v),
        CurvePt(TangentBundleSpace(v.pt.pt_type.manifold_type)) do t
            return TangentBundlePt(α * v.v(t).x)
        end)
end

function inner(v1::CurveTV, v2::CurveTV)

    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end

    reltol, abstol = concretize_tols(v1, v2, reltol=PARAMS.quad_rel_tol, abstol=PARAMS.quad_abs_tol)
    v1inner = v1.v
    v2inner = v2.v
    I, err = QuadGK.quadgk(s -> inner(v1inner(s).x, v2inner(s).x), 0.0, 1.0, rtol = reltol, atol = abstol)
    return I
end

function geodesic_at(t::Number, x1::CurvePt, x2::CurvePt)
    return CurvePt(x1.pt_type.manifold_type) do s
        return geodesic_at(t, x1(s), x2(s))
    end
end

function geodesic_distance(x1::CurvePt, x2::CurvePt)
    # Calculates simple L2-like distance
    # There are also other reasonable Riemannian structures for this space
    reltol, abstol = concretize_tols(x1, x2, reltol=PARAMS.quad_rel_tol, abstol=PARAMS.quad_abs_tol)

    return sqrt(quadgk(t -> geodesic_distance(x1(t), x2(t))^2, 0.0, 1.0, rtol = reltol, atol = abstol)[1])
end

function ambient_distance(x1::CurvePt, x2::CurvePt)

    reltol, abstol = concretize_tols(x1, x2, reltol=PARAMS.quad_rel_tol, abstol=PARAMS.quad_abs_tol)

    return sqrt(quadgk(t -> ambient_distance(x1(t), x2(t))^2, 0.0, 1.0, rtol = reltol, atol = abstol)[1])
end

function inner_amb(x1::CurvePt, x2::CurvePt)

    reltol, abstol = concretize_tols(x1, x2, reltol=PARAMS.quad_rel_tol, abstol=PARAMS.quad_abs_tol)

    return sqrt(quadgk(t -> inner_amb(x1(t), x2(t))^2, 0.0, 1.0, rtol = reltol, atol = abstol)[1])
end

function exp(v::CurveTV)
    return CurvePt(v.pt.pt_type) do t
        return exp(v.v(t).x)
    end
end

function log(x::CurvePt, y::CurvePt)
    return CurveTV(x,
        CurvePt(TangentBundleSpace(x.pt_type.manifold_type)) do t
            return TangentBundlePt(log(x(t), y(t)))
        end)
end

function parallel_transport_geodesic(v::CurveTV, to_point::CurvePt)
    return CurveTV(to_point,
        CurvePt(TangentBundleSpace(to_point.pt_type.manifold_type)) do t
            return TangentBundlePt(parallel_transport_geodesic(v.v(t).x, to_point(t)))
        end)
end
