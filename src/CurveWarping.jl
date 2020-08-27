
"""
    CurveWarpingSpace(knots)

Space of warpings of the interval [0, 1] on knots `knots`.
"""
struct CurveWarpingSpace{TK<:AbstractVector,TIM<:Interpolations.MonotonicInterpolationType} <: Manifold{ℝ}
    knots::TK
    interpolation_method::TIM
end

function CurveWarpingSpace(knots::TK) where TK<:AbstractVector
    method = FritschButlandMonotonicInterpolation()
    return CurveWarpingSpace{TK,typeof(method)}(knots, method)
end

function manifold_dimension(x::CurveWarpingSpace)
    return Inf
end

function make_warping(x::AbstractVector, y::Vector{<:Real}, interpolation_method)
    return interpolate(Vector(x), y, interpolation_method)
end

function identity!(cws::CurveWarpingSpace, p)
    interpolate(p.knots, p.knots, cws.interpolation_method)
end

function inv(cws::CurveWarpingSpace, p::Interpolations.MonotonicInterpolationType)
    ys = Interpolations.coefficients(cwp.interp)
    interp_grid = cws.knots
    ys[end] = 1 # makes things easier for autodiff
                # TODO: make it better somehow?
    m = [1/derivative(cwp, t) for t in interp_grid]
    newintp = interpolate(ys, interp_grid, KnownDerivativesMonotonicInterpolation(m))
    return newintp
end

function compose(cws::CurveWarpingSpace, p1, p2)
    interp_grid = x2.interp.knots
    return make_warping(interp_grid, [x1(x2(t)) for t in interp_grid], cws.interpolation_method)
end

"""
    CurveWarpingAction(cws::CurveWarpingSpace)

Space of actions of the group of curve warpings on the manifold of curves.
"""
struct CurveWarpingAction{TCWS<:CurveWarpingSpace} <: AbstractGroupOperation
    cws::TCWS
end
