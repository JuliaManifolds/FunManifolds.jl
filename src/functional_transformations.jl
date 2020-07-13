
@doc raw"""
    velocity_curve(p, c::Curve, backend)

Velocity curve for a given curve `p` of type `c` calculated using differentiation
defined by `backend`.
If $c$ is a function such that $c\colon [0, 1] \to M$, then
`velocity_curve(c)(t)` is a vector tangent to `c` at `t`.
"""
function velocity_curve(M::Manifold, p, backend::Manifolds.AbstractRiemannianDiffBackend)
    return t -> Manifolds.differential(M, p, t, backend)
end

struct ProjectedDifferenceBackend{TDT<:Number} <: Manifolds.AbstractRiemannianDiffBackend
    dt::TDT
end

ProjectedDifferenceBackend() = ProjectedDifferenceBackend{Float64}(1e-7)

function velocity_curve(M::Manifold, p, backend::ProjectedDifferenceBackend)
    return t -> project(M.M, p(t), (embed(p(t + backend.dt)) - embed(p(t))) / backend.dt)
end

@doc raw"""
    curve_length(p, c::FunctionCurveSpace, a = 0.0, b = 1.0, dtype = Val(:continuous))

Returns length of the curve `c` between parameters `a` and `b`. By default
`a = 0.0` and `b = 1.0`. Calculates velocity using method `dtype`.
"""
function curve_length(
    M::Manifold,
    p,
    a = 0.0,
    b = 1.0,
    backend = Manifolds.rdifferential_backend(),
)
    v = velocity_curve(M, p, backend)
    #TODO: add quadrature tolerances to curve_length?
    #TODO: add velocity calculation options
    i_val, error = quadgk(t -> norm(M, p(t), v(t)), a, b)
    return i_val
end

"""
    q_function(M::Manifold, p, X)

Rescales given tangent vector `X` at point `p` to square root of its original length.
"""
function q_function(M::Manifold, p, X)
    norm_X = norm(M, p, X)
    return if norm_X ≈ 0
        return zero_tangent_vector(M, p)
    else
        return (1 / sqrt(norm_X)) * X
    end
end

"""
    srvf(M, f, backend)

Calculate the Square Root Velocity Function of a given curve `f`.
Parameter `backend` describes type of differentiaion used. Supported values:
* `Val(:continuous)` -- symbolic differentiation resulting in a continuous
curve,
* `Val(:forward_diff)` -- forward finite differentiaion of a discretized curve
`c` based on its grid.
"""
function srvf(M::Manifold, f, backend::Val{:continuous})
    vel = velocity_curve(M, f, backend)
    return t -> q_function(M, f(t), vel(t))
end
