"""
    CurveWarpingSRSFSpace(knots, quadweights = get_quad_weights(knots))

Space of SRSFs of curve warpings on knots `knots`. Points are represented by vectors
in ambient space (ℝ^(N+1) if there are `N` knots). Given quadrature weights are used
to compute inner products.
"""
struct CurveWarpingSRSFSpace{TK<:AbstractVector,TW<:AbstractVector} <: Manifold{ℝ}
    knots::TK
    quadweights::TW
end

function CurveWarpingSRSFSpace(knots::AbstractVector)
    ws = get_quad_weights(knots)
    return CurveWarpingSRSFSpace(knots, ws)
end

function manifold_dimension(M::CurveWarpingSRSFSpace)
    return length(M.knots)
end

function representation_size(M::CurveWarpingSRSFSpace)
    return (length(M.knots)+1,)
end

function isapprox(M::CurveWarpingSRSFSpace, p, q; kwargs...)
    return isapprox(p, q; kwargs...)
end


function isapprox(M::CurveWarpingSRSFSpace, p, X, Y; kwargs...)
    return isapprox(X, Y; kwargs...)
end

"""
    inner_ambient(M::CurveWarpingSRSFSpace, p, q)

Inner product in ambient space of points `p` and `q` on `M`. Takes into account quadrature
weights.
"""
function inner_ambient(M::CurveWarpingSRSFSpace, p, q)
    return sum(M.quadweights .* p .* q)
end

function distance(M::CurveWarpingSRSFSpace, p, q)
    return acos(clamp(inner_ambient(M, p, q), -1, 1))
end

function embed!(::CurveWarpingSRSFSpace, q, p)
    copyto!(q, p)
    return q
end

function embed!(M::CurveWarpingSRSFSpace, Y, p, X)
    copyto!(Y, X)
    return Y
end

function project!(M::CurveWarpingSRSFSpace, q, p)
    scale = sqrt(inner_ambient(M, p, p))
    q .= p ./ scale
    return q
end

function project!(M::CurveWarpingSRSFSpace, Y, p, X)
    Y .= X .- inner_ambient(M, p, X) .* p
    return Y
end

function get_quad_weights(nodes::AbstractRange)
    n = length(nodes)
    dt = 1/(n-1)
    if n%2 == 1
        # Simpson's rule
        return dt/3 .* [1; [i%2==0 ? 2 : 4 for i in 1:n-2]; 1]
    else
        # trapezoidal rule
        return dt/2 .* [1; [2 for i in 1:n-2]; 1]
    end

    # Simpson's 3/8
    #n%3 == 1 || error("argument mod 3 must be 1")
    #return 3dt/8 .* [1; [i%3==0 ? 2 : 3 for i in 1:n-2]; 1]
end

function exp!(M::CurveWarpingSRSFSpace, q, p, X)
    θ = norm(M, p, X)
    q .= cos(θ) .* p .+ Manifolds.usinc(θ) .* X
    return q
end

function log!(M::CurveWarpingSRSFSpace, X, p, q)
    cosθ = clamp(inner_ambient(M, p, q), -1, 1)
    θ = acos(cosθ)
    X .= (q .- cosθ .* p) ./ Manifolds.usinc(θ)
    return project!(M, X, p, X)
end

function inner(M::CurveWarpingSRSFSpace, p, X, Y)
    return inner_ambient(M, X, Y)
end

function vector_transport_to!(M::CurveWarpingSRSFSpace, Y, p, X, q, ::ParallelTransport)
    copyto!(Y, X)
    X_pq = log(M, p, q)
    X1 = norm(M, p, X_pq)
    if X1 > 0
        sum_pq = p + q
        factor = 2 * inner_ambient(M, X, q) / inner_ambient(M, sum_pq, sum_pq)
        Y .-= factor .* sum_pq
    end
    return Y
end

function zero_tangent_vector!(M::CurveWarpingSRSFSpace, X, p)
    fill!(X, 0)
    return X
end
