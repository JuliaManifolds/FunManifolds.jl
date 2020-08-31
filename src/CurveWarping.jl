
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

function make_warping(x::AbstractVector, y::Vector{<:Real}, interpolation_method)
    return interpolate(Vector(x), y, interpolation_method)
end

struct WarpingCompositionOperation <: AbstractGroupOperation end

const CurveWarpingGroup{TCWS<:CurveWarpingSpace} = GroupManifold{ℝ,TCWS,WarpingCompositionOperation}

function CurveWarpingGroup(cws::CurveWarpingSpace)
    return GroupManifold(cws, WarpingCompositionOperation())
end

function identity!(cws::CurveWarpingGroup, q, p)
    q.copyto!(q, cws.knots)
    return q
end

function inv!(cws::CurveWarpingGroup, q, p)
    ys = Interpolations.coefficients(cwp.interp)
    interp_grid = cws.knots
    ys[end] = 1 # makes things easier for autodiff
                # TODO: make it better somehow?
    m = [1/derivative(cwp, t) for t in interp_grid]
    newintp = interpolate(ys, interp_grid, KnownDerivativesMonotonicInterpolation(m))
    return newintp
end

function compose(cws::CurveWarpingGroup, p1, p2)
    return make_warping(interp_grid, [p1(p2(t)) for t in cws.knots], cws.interpolation_method)
end

"""
    CurveWarpingAction(M::Manifold, cwg::CurveWarpingGroup)

Space of left actions of the group of curve warpings `cwg` on the manifold `M` of curves.
"""
struct CurveWarpingAction{TM<:Manifold,TCWG<:CurveWarpingGroup} <: AbstractGroupAction{LeftAction}
    manifold::TM
    cwg::TCWG
end

function Manifolds.g_manifold(cwa::CurveWarpingAction)
    return cwa.manifold
end

function Manifolds.base_group(cwa::CurveWarpingAction)
    return cwa.cwg
end
