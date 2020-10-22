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
    return (length(M.knots) + 1,)
end

function isapprox(M::CurveWarpingSRSFSpace, p, q; kwargs...)
    return isapprox(distance(M, p, q), 0; kwargs...)
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
    dt = 1 / (n - 1)
    return dt / 2 .* [1; [2 for i in 1:(n - 2)]; 1]
    # not that precise for some reason
    #=if n % 2 == 1
        # Simpson's rule
        return dt / 3 .* [1; [i % 2 == 0 ? 2 : 4 for i in 1:(n - 2)]; 1]
    else
        # trapezoidal rule
        return dt / 2 .* [1; [2 for i in 1:(n - 2)]; 1]
    end=#

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

const CurveWarpingSRSFGroup{TCWS<:CurveWarpingSRSFSpace} =
    GroupManifold{ℝ,TCWS,WarpingCompositionOperation}

function CurveWarpingSRSFGroup(cws::CurveWarpingSRSFSpace)
    return GroupManifold(cws, WarpingCompositionOperation())
end

function identity(cwg::CurveWarpingSRSFGroup, p)
    return ones(eltype(p), size(p))
end

function inv(cwg::CurveWarpingSRSFGroup, p)
    cws = CurveWarpingSpace(cwg.manifold.knots)
    pinv = inv(CurveWarpingGroup(cws), reverse_srsf(cws, p))
    return srsf(cws, pinv)
end
inv(cwg::CurveWarpingSRSFGroup, p::Identity) = p

function compose(cwg::CurveWarpingSRSFGroup, p1, p2)
    knots = cwg.manifold.knots
    cws = CurveWarpingSpace(knots)
    p1inv = reverse_srsf(cws, p1)
    p2w = make_warping(cws, p2)
    p2warped = map(t -> p2w(p1inv(t)), knots)
    return p2warped .* p1
end
function compose(cwg::TCWG, p1::Identity{TCWG}, p2) where {TCWG<:CurveWarpingSRSFGroup}
    return p2
end
function compose(cwg::TCWG, p1, p2::Identity{TCWG}) where {TCWG<:CurveWarpingSRSFGroup}
    return p1
end
function compose(
    cwg::TCWG,
    p1::Identity{TCWG},
    p2::Identity{TCWG},
) where {TCWG<:CurveWarpingSRSFGroup}
    return p1
end


"""
    CurveWarpingSRSFAction(M::Manifold, cwg::CurveWarpingSRSFGroup)

Space of left actions of the group of SRSFs of curve warpings `cwg` on the manifold `M`
of SRVFs of curves.
"""
struct CurveWarpingSRSFAction{TM<:Manifold,TCWG<:CurveWarpingSRSFGroup} <:
       AbstractGroupAction{LeftAction}
    manifold::TM
    cwg::TCWG
end


function Manifolds.g_manifold(A::CurveWarpingSRSFAction)
    return A.manifold
end

function Manifolds.base_group(A::CurveWarpingSRSFAction)
    return A.cwg
end

function apply!(A::CurveWarpingSRSFAction{<:DCurves}, q, a, p)
    itp = make_interpolant(A.manifold, p)
    a_rev = reverse_srsf(CurveWarpingSpace(A.cwg.manifold.knots), a)
    ts = map(a_rev, A.cwg.manifold.knots)
    rep_size = representation_size(A.manifold.manifold)
    for (i, t) in zip(get_iterator(A.manifold), ts)
        copyto!(_write(A.manifold, rep_size, q, i), itp(t) * a[i])
    end
    return q
end

function apply!(
    A::CurveWarpingSRSFAction{<:DCurves,TCWG},
    q,
    ::Identity{TCWG},
    p,
) where {TCWG<:CurveWarpingSRSFGroup}
    copyto!(q, p)
    return q
end

function Manifolds.inverse_apply(A::CurveWarpingSRSFAction, a, p)
    inva = inv(base_group(A), a)
    return apply(A, inva, p)
end
