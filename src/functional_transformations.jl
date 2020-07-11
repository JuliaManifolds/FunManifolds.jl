
@doc raw"""
    velocity_curve(p, c::Curve, backend)

Velocity curve for a given curve `p` of type `c` calculated using differentiation
defined by `backend`.
If $c$ is a function such that $c\colon [0, 1] \to M$, then
`velocity_curve(c)(t)` is a vector tangent to `c` at `t`.
"""
function velocity_curve(
    p,
    c::Manifolds.Curve,
    backend::Manifolds.AbstractRiemannianDiffBackend,
)
    return t -> Manifolds.differential(p, c, t, backend)
end

struct ProjectedDifferenceBackend{TDT<:Number} <: Manifolds.AbstractRiemannianDiffBackend
    dt::TDT
end

ProjectedDifferenceBackend() = ProjectedDifferenceBackend{Float64}(1e-7)

function velocity_curve(p, M::Manifolds.Curve, backend::ProjectedDifferenceBackend)
    return t -> project(M.M, p(t), (embed(p(t + backend.dt)) - embed(p(t))) / backend.dt)
end

@doc raw"""
    curve_length(p, c::FunctionCurveSpace, a = 0.0, b = 1.0, dtype = Val(:continuous))

Returns length of the curve `c` between parameters `a` and `b`. By default
`a = 0.0` and `b = 1.0`. Calculates velocity using method `dtype`.
"""
function curve_length(
    p,
    c::FunctionCurveSpace,
    a = 0.0,
    b = 1.0,
    backend = Manifolds.rdifferential_backend(),
)
    v = velocity_curve(p, Manifolds.Curve(c.M), backend)
    #TODO: add quadrature tolerances to curve_length?
    #TODO: add velocity calculation options
    i_val, error = quadgk(t -> norm(c.M, p(t), v(t)), a, b)
    return i_val
end
