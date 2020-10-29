
@doc raw"""
    velocity_curve(p, c::Curve, backend)

Velocity curve for a given curve `p` of type `c` calculated using differentiation
defined by `backend`.
If $c$ is a function such that $c\colon [0, 1] \to M$, then
`velocity_curve(c)(t)` is a vector tangent to `c` at `t`.
"""
function velocity_curve(M::Manifold, p, backend::AbstractRiemannianDiffBackend)
    return t -> Manifolds.differential(M, p, t, backend)
end

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
function srvf(M::Manifold, f, backend::AbstractRiemannianDiffBackend)
    vel = velocity_curve(M, f, backend)
    return t -> q_function(M, f(t), vel(t))
end

"""
    transport_srvf(M::Manifold, c_p, c_X, p)

Transports all tangent vectors `c_X` of curve `c_p` to a given
point `p` along shortest geodesics.
`M` is the respective manifold of curves.
"""
function transport_srvf(M::Manifold, c_p, c_X, p)
    c_out = allocate_result(M, transport_srvf, c_X, c_p)
    transport_srvf!(M, c_out, c_p, c_X, p)
    return c_out
end

"""
    tsrvf(M::Manifold, c, p, backend::AbstractRiemannianDiffBackend)

Calculate SRVF of curve `c` from manifold `M` using method `backend` and transport to `p`.
"""
function tsrvf(M::Manifold, c, p, backend::AbstractRiemannianDiffBackend)
    c_X = srvf(M, c, backend)
    return transport_srvf(M, c, c_X, p)
end

function reverse_srvf(M::Manifold, c_X, initial_point)
    c_out = allocate_result(M, reverse_srvf, c_X)
    reverse_srvf!(M, c_out, c_X, initial_point)
    return c_out
end

"""
    srsf(M::CurveWarpingSpace, p)

Square Root Slope Function of curve warping `p`.
"""
function srsf(M::CurveWarpingSpace, p)
    grid = M.knots
    coeffs = map(p, grid)
    grads = [
        (coeffs[i + 1] - coeffs[i]) / (grid[i + 1] - grid[i]) for i in 1:(length(grid) - 1)
    ]
    push!(grads, grads[end])
    vals = map(sqrt, grads)
    return project(CurveWarpingSRSFSpace(grid), vals)
end

"""
    reverse_srsf(c)

Reverse Square Root Slope Function of a given SRSF of a curve warping.
"""
function reverse_srsf(M::CurveWarpingSpace, p)
    grid = M.knots
    steps = diff(grid)
    ys = accumulate(+, [zero(eltype(p)); steps .* p[1:(end - 1)] .^ 2])
    ys ./= ys[end]
    return make_warping(M, ys)
end
function reverse_srsf(M::CurveWarpingSpace, p::Identity)
    reverse_srsf(M, p.p)
end
