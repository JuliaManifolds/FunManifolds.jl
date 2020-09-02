
"""
    CurveWarpingSpace(knots)

Space of warpings of the interval [0, 1] on knots `knots`.
"""
struct CurveWarpingSpace{TK<:AbstractVector,TIM<:Interpolations.MonotonicInterpolationType} <: Manifold{ℝ}
    knots::TK
    interpolation_method::TIM
end

function CurveWarpingSpace(knots::TK) where {TK<:AbstractVector}
    method = FritschButlandMonotonicInterpolation()
    return CurveWarpingSpace{TK,typeof(method)}(knots, method)
end

function manifold_dimension(x::CurveWarpingSpace)
    return Inf
end

function make_warping(M::CurveWarpingSpace, y::AbstractVector)
    return extrapolate(interpolate(convert(Vector, M.knots), convert(Vector, y), M.interpolation_method), Flat())
end

struct WarpingCompositionOperation <: AbstractGroupOperation end

const CurveWarpingGroup{TCWS<:CurveWarpingSpace} = GroupManifold{ℝ,TCWS,WarpingCompositionOperation}

function CurveWarpingGroup(cws::CurveWarpingSpace)
    return GroupManifold(cws, WarpingCompositionOperation())
end

function identity(cwg::CurveWarpingGroup, p)
    return make_warping(cwg.manifold, cwg.manifold.knots)
end

function inv(cwg::CurveWarpingGroup, p)
    ys = map(p, cwg.manifold.knots)
    ys[end] = 1 # makes things easier for autodiff
                # TODO: make it better somehow?
    m = [1/Interpolations.gradient1(p, t) for t in cwg.manifold.knots]
    return extrapolate(interpolate(ys, convert(Vector, cwg.manifold.knots), KnownDerivativesMonotonicInterpolation(m)), Flat())
end
inv(cwg::CurveWarpingGroup, p::Identity) = p

function compose(cwg::CurveWarpingGroup, p1, p2)
    return make_warping(cwg.manifold, map(t -> p1(p2(t)), cwg.manifold.knots))
end
function compose(cwg::TCWG, p1::Identity{TCWG}, p2) where {TCWG<:CurveWarpingGroup}
    return p2
end
function compose(cwg::TCWG, p1, p2::Identity{TCWG}) where {TCWG<:CurveWarpingGroup}
    return p1
end
function compose(cwg::TCWG, p1::Identity{TCWG}, p2::Identity{TCWG}) where {TCWG<:CurveWarpingGroup}
    return p1
end

"""
    CurveWarpingAction(M::Manifold, cwg::CurveWarpingGroup)

Space of left actions of the group of curve warpings `cwg` on the manifold `M` of curves.
"""
struct CurveWarpingAction{TM<:Manifold,TCWG<:CurveWarpingGroup} <: AbstractGroupAction{LeftAction}
    manifold::TM
    cwg::TCWG
end

function Manifolds.g_manifold(A::CurveWarpingAction)
    return A.manifold
end

function Manifolds.base_group(A::CurveWarpingAction)
    return A.cwg
end

function apply!(A::CurveWarpingAction{<:DCurves}, q, a, p)
    ts = map(a, A.cwg.manifold.knots)
    itp = make_interpolant(A.manifold, p)
    grid = A.manifold.grid
    rep_size = representation_size(A.manifold.manifold)
    for (i, t) in zip(get_iterator(A.manifold), ts)
        copyto!(_write(A.manifold, rep_size, q, i), itp(t))
    end
    return q
end

function apply!(A::CurveWarpingAction{<:DCurves,TCWG}, q, ::Identity{TCWG}, p) where {TCWG<:CurveWarpingGroup}
    copyto!(q, p)
    return q
end

function Manifolds.inverse_apply(A::CurveWarpingAction, a, p)
    inva = inv(base_group(A), a)
    return apply(A, inva, p)
end

function isapprox(M::CurveWarpingSpace, x, y; kwargs...)
    xc = map(x, M.knots)
    yc = map(y, M.knots)
    return isapprox(xc, yc; kwargs...)
end
