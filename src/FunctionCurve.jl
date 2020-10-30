
"""
    FunctionCurveSpace(m, grid)

Space of curves in manifold `m`. When needed for calculations, uses grid `grid`
but does not discretize the curve when not necessary.
They have no embedding in a finite-dimensional euclidean space
(use `DiscretizedCurves` if you need this feature).
"""
struct FunctionCurveSpace{F,TM<:Manifold{F},R<:AbstractRange} <: Manifold{F}
    manifold::TM
    approx_grid::R
end

function FunctionCurveSpace(M::TM) where {F,TM<:Manifold{F}}
    approx_grid = 0:0.01:1.0
    return FunctionCurveSpace{F,TM,typeof(approx_grid)}(M, approx_grid)
end

"""
    VectorizedFunction(f)

A Function with vector space operations defined.
"""
struct VectorizedFunction{F}
    f::F
end

*(a::Number, f::VectorizedFunction) = VectorizedFunction(t -> a * f.f(t))
function +(f1::VectorizedFunction, f2::VectorizedFunction)
    return VectorizedFunction(t -> f1.f(t) + f2.f(t))
end
function -(f1::VectorizedFunction, f2::VectorizedFunction)
    return VectorizedFunction(t -> f1.f(t) - f2.f(t))
end
-(f::VectorizedFunction) = VectorizedFunction(t -> -f.f(t))

(f::VectorizedFunction)(t) = f.f(t)

function distance(M::FunctionCurveSpace, p1, p2)
    # Calculates simple L2-like distance
    # There are also other reasonable Riemannian structures for this space
    reltol, abstol = concretize_tols(
        M,
        p1,
        p2,
        reltol = PARAMS.quad_rel_tol,
        abstol = PARAMS.quad_abs_tol,
    )

    return sqrt(quadgk(
        t -> distance(M.manifold, p1(t), p2(t))^2,
        0.0,
        1.0,
        rtol = reltol,
        atol = abstol,
    )[1])
end

function exp(M::FunctionCurveSpace, p, X)
    return t -> exp(M.manifold, p(t), X(t))
end

function geodesic(M::FunctionCurveSpace, p, X)
    return t -> (s -> geodesic(M.manifold, p(t), X(t), s))
end

function injectivity_radius(M::FunctionCurveSpace)
    # TODO: check this
    return injectivity_radius(M.manifold)
end

function inner(M::FunctionCurveSpace, p, X1, X2)
    reltol, abstol = concretize_tols(
        M,
        X1,
        X2,
        reltol = PARAMS.quad_rel_tol,
        abstol = PARAMS.quad_abs_tol,
    )
    I, err = QuadGK.quadgk(
        s -> inner(M.manifold, p(s), X1(s), X2(s)),
        0.0,
        1.0,
        rtol = reltol,
        atol = abstol,
    )
    return I
end

function inverse_retract(
    M::FunctionCurveSpace,
    p,
    q,
    method::AbstractInverseRetractionMethod,
)
    return VectorizedFunction(t -> inverse_retract(M.manifold, p(t), q(t), method))
end

function inverse_retract(M::FunctionCurveSpace, p, q)
    return VectorizedFunction(t -> inverse_retract(M.manifold, p(t), q(t)))
end

function isapprox(
    M::FunctionCurveSpace,
    p1,
    p2;
    atol = atoldefault(M, p1, p2),
    rtol = rtoldefault(M, p1, p2),
)
    #=for i in M.approx_grid
        if !(isapprox(M.manifold, p1(i), p2(i), atol = atol, rtol = rtol))
            println("*** $i, $(p1(i)), $(p2(i)).")
        end
    end=#
    return all(
        isapprox(M.manifold, p1(i), p2(i), atol = atol, rtol = rtol) for i in M.approx_grid
    )
end

function isapprox(
    M::FunctionCurveSpace,
    p,
    X1,
    X2;
    atol = atoldefault(M, X1, X2),
    rtol = rtoldefault(M, X1, X2),
)
    #TODO add tolerance parameters to isapprox for this quadrature?
    I, err = quadgk(t -> norm(M.manifold, p(t), X1(t) - X2(t)), 0.0, 1.0)
    return I <= 2 * err + atol && I >= -2 * err - atol
end

function log(M::FunctionCurveSpace, p, q)
    return VectorizedFunction(t -> log(M.manifold, p(t), q(t)))
end

function manifold_dimension(M::FunctionCurveSpace)
    return Inf
end

function mid_point(M::FunctionCurveSpace, p, q)
    return t -> mid_point(M.manifold, p(t), q(t))
end

function retract(M::FunctionCurveSpace, p, X)
    return VectorizedFunction(t -> retract(M.manifold, p(t), X(t)))
end

function retract(M::FunctionCurveSpace, p, X, method::AbstractRetractionMethod)
    return VectorizedFunction(t -> retract(M.manifold, p(t), X(t), method))
end

function transport_srvf(M::FunctionCurveSpace, c_p, c_X, p)
    return t -> vector_transport_to(M.manifold, c_p(t), c_X(t), p)
end

function vector_transport_direction(
    M::FunctionCurveSpace,
    p,
    X,
    d,
    method::AbstractVectorTransportMethod,
)
    return VectorizedFunction(
        t -> vector_transport_direction(M.manifold, p(t), X(t), d(t), method),
    )
end

function zero_tangent_vector(M::FunctionCurveSpace, p)
    return VectorizedFunction(t -> zero_tangent_vector(M.manifold, p(t)))
end

function zero_vector(M::Manifolds.TangentBundleFibers{<:FunctionCurveSpace}, p)
    return VectorizedFunction(t -> zero_tangent_vector(M.manifold.manifold, p(t)))
end
