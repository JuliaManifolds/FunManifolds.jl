
"""
    FunctionCurveSpace(m, grid)

Space of curves in manifold `m`. When needed for calculations, uses grid `grid`
but does not discretize the curve when not necessary.
They have no embedding in a finite-dimensional euclidean space
(use `DCurveSpace` if you need this feature).
"""
struct FunctionCurveSpace{F, M<:Manifold{F}, R<:AbstractRange} <: Manifold{F}
    M::M
    approx_grid::R
end

function FunctionCurveSpace(M::TM) where {TM<:Manifold}
    approx_grid = 0:0.01:1.0
    return FunctionCurveSpace{TM,typeof(approx_grid)}(M, approx_grid)
end

function manifold_dimension(x::FunctionCurveSpace)
    return Inf
end

function isapprox(M::FunctionCurveSpace, p1, p2; atol = atoldefault(p1, p2), rtol = rtoldefault(p1, p2))
    for i in M.approx_grid
        if !(isapprox(M.M, p1(i), p2(i), atol = atol, rtol = rtol))
            println("*** $i, $(p1(i)), $(p2(i)).")
        end
    end
    return all(isapprox(M.M, p1(i), p2(i), atol = atol, rtol = rtol) for i ∈ M.approx_grid)
end


function zero_tangent_vector(M::FunctionCurveSpace, p)
    return t -> zero_tangent_vector(M.M, p(t))
end

function isapprox(M::FunctionCurveSpace, p, X1, X2; atol = atoldefault(X1, X2), rtol = rtoldefault(X1, X2))
    #TODO add tolerance parameters to isapprox for this quadrature?
    I, err = quadgk(t -> norm(M.M, p(t), X1(t) - X2(t)), 0.0, 1.0)
    return I <= 2*err+atol && I >= -2*err-atol
end

function inner(M::FunctionCurveSpace, p, v1, v2)
    reltol, abstol = concretize_tols(v1, v2, reltol=PARAMS.quad_rel_tol, abstol=PARAMS.quad_abs_tol)
    v1inner = v1.v
    v2inner = v2.v
    I, err = QuadGK.quadgk(s -> inner(M.M, p(s), v1(s), v2(s)), 0.0, 1.0, rtol = reltol, atol = abstol)
    return I
end

function geodesic(M::FunctionCurveSpace, p, X)
    return t -> geodesic(M.M, p, X, t)
end

function distance(M::FunctionCurveSpace, p1, p2)
    # Calculates simple L2-like distance
    # There are also other reasonable Riemannian structures for this space
    reltol, abstol = concretize_tols(p1, p2, reltol=PARAMS.quad_rel_tol, abstol=PARAMS.quad_abs_tol)

    return sqrt(quadgk(t -> distance(M.M, p1(t), p2(t))^2, 0.0, 1.0, rtol = reltol, atol = abstol)[1])
end

function exp(M::FunctionCurveSpace, p, X)
    return t -> exp(M.M, p(t), X(t))
end

function log(M::FunctionCurveSpace, p, q)
    return t -> log(M.M, p(t), q(t))
end

function vector_transport_direction(M::FunctionCurveSpace, p, X, d, method::AbstractVectorTransportMethod)
    return t -> vector_transport_direction(M.M, p(t), X(t), d(t), method)
end
