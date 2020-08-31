
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
    return interpolate(convert(Vector, M.knots), convert(Vector, y), M.interpolation_method)
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
    ys = Interpolations.coefficients(p)
    ys[end] = 1 # makes things easier for autodiff
                # TODO: make it better somehow?
    m = [1/derivative(cwg.manifold, t) for t in cwg.manifold.knots]
    return interpolate(ys, convert(Vector, cwg.manifold.knots), KnownDerivativesMonotonicInterpolation(m))
end
inv(cwg::CurveWarpingGroup, p::Identity) = p

function compose(cwg::CurveWarpingGroup, p1, p2)
    y2s = Interpolations.coefficients(p2)
    return make_warping(cwg.manifold, map(t -> p1(t), y2s))
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
    error("TODO")
end

function Manifolds.inverse_apply(A::CurveWarpingAction, a, p)
    inva = inv(base_group(A), a)
    return apply(A, inva, p)
end
