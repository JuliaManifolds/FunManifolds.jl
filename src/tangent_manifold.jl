
"""
    TSpaceManifold(pt)

Default manifold structure for tangent space at point `pt`.
"""
struct TSpaceManifold{Pt <: Point} <: Manifold
    pt::Pt
end

function ambient_shape(m::TSpaceManifold)
    return ambient_shape(gettype(m.pt))
end

"""
    TSpaceManifoldPt(x)

A point on TSpaceManifold represented by a tangent vector `x`.
"""
struct TSpaceManifoldPt{TV <: TangentVector} <: Point
    x::TV
end

function copyto!(x_to::TSpaceManifoldPt, x_from::TSpaceManifoldPt)
    copyto!(x_to.x, x_from.x)
    return x_to
end

function deepcopy(x::TSpaceManifoldPt)
    return TSpaceManifoldPt(deepcopy(x.x))
end

function +(v1::TSpaceManifoldPt, v2::TSpaceManifoldPt)
    DEBUG && if !(at_point(v1.x) ≈ at_point(v2.x))
        error("Can't add tangent vectors from different tangent spaces")
    end
    return TSpaceManifoldPt(v1.x + v2.x)
end

function add_vec!(v1::TSpaceManifoldPt, v2::TSpaceManifoldPt)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    add_vec!(v1.x, v2.x)
    return v1
end

@inline function add_vec!(v1::BNBArray, v2::AbstractArray, at_pt::AbstractArray, m::TSpaceManifold)
    add_vec!(v1, v2, point2ambient(m.pt), gettype(m.pt))
end

function -(v1::TSpaceManifoldPt, v2::TSpaceManifoldPt)
    DEBUG && if !(at_point(v1.x) ≈ at_point(v2.x))
        error("Can't subtract tangent vectors from different tangent spaces")
    end
    return TSpaceManifoldPt(v1.x - v2.x)
end

function sub_vec!(v1::TSpaceManifoldPt, v2::TSpaceManifoldPt)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    sub_vec!(v1.x, v2.x)
    return v1
end

@inline function sub_vec!(v1::BNBArray, v2::AbstractArray, at_pt::AbstractArray, m::TSpaceManifold)
    sub_vec!(v1, v2, point2ambient(m.pt), gettype(m.pt))
end

function *(α::Real, v::TSpaceManifoldPt)
    return TSpaceManifoldPt(α * v.x)
end

function mul_vec!(v::TSpaceManifoldPt, α::Real)
    mul_vec!(v.x, α)
    return v
end

@inline function mul_vec!(v::BNBArray, α::Real, at_pt::AbstractArray, m::TSpaceManifold)
    mul_vec!(v, α, point2ambient(m.pt), gettype(m.pt))
end

function gettype(x::TSpaceManifoldPt)
    return TSpaceManifold(at_point(x.x))
end

function manifold_dimension(m::TSpaceManifold)
    return manifold_dimension(m.pt)
end

function isapprox(v1::TSpaceManifoldPt, v2::TSpaceManifoldPt; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    return isapprox(v1.x, v2.x, atol = atol, rtol = rtol)
end

"""
    TSpaceManifoldTV(x, v)

Tangent vector to a TSpaceManifold from tangent space at point `x`
represented by a tangent vector `v`.
"""
struct TSpaceManifoldTV{TV <: TangentVector} <: TangentVector
    at_pt::TSpaceManifoldPt{TV}
    v::TV
end

function copyto!(v_to::TSpaceManifoldTV, v_from::TSpaceManifoldTV)
    copyto!(v_to.at_pt, v_from.at_pt)
    copyto!(v_to.v, v_from.v)
    return v_to
end

function deepcopy(v::TSpaceManifoldTV)
    return TSpaceManifoldTV(deepcopy(v.at_pt), deepcopy(v.v))
end

function +(v1::TSpaceManifoldTV, v2::TSpaceManifoldTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't add tangent vectors from different tangent spaces")
    end
    return TSpaceManifoldTV(v1.at_pt, v1.v + v2.v)
end

function add_vec!(v1::TSpaceManifoldTV{TV}, v2::TSpaceManifoldTV{TV}) where TV <: TangentVector
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    add_vec!(v1.v, v2.v)
    return v1
end

function -(v1::TSpaceManifoldTV, v2::TSpaceManifoldTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't subtract tangent vectors from different tangent spaces")
    end
    return TSpaceManifoldTV(v1.at_pt, v1.v - v2.v)
end

function sub_vec!(v1::TSpaceManifoldTV, v2::TSpaceManifoldTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    sub_vec!(v1.v, v2.v)
    return v1
end

function *(α::Real, v::TSpaceManifoldTV)
    return TSpaceManifoldTV(v.at_pt, α * v.v)
end

function mul_vec!(v::TSpaceManifoldTV, α::Real)
    mul_vec!(v.v, α)
    return v
end

function isapprox(v1::TSpaceManifoldTV, v2::TSpaceManifoldTV; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    if !(isapprox(v1.at_pt, v2.at_pt, atol = atol, rtol = rtol))
        return false
    end
    return isapprox(v1.v, v2.v, atol = atol, rtol = rtol)
end

function inner(v1::TSpaceManifoldTV, v2::TSpaceManifoldTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return inner(v1.v, v2.v)
end

function inner(v1::AbstractArray, v2::AbstractArray, p::AbstractArray, m::TSpaceManifold)
    return inner(v1, v2, point2ambient(m.pt), gettype(m.pt))
end

function point2ambient(p::TSpaceManifoldPt)
    return tangent2ambient(p.x)
end

function ambient2point(amb::AbstractArray, m::TSpaceManifold)
    return TSpaceManifoldPt(ambient2tangent(amb, m.pt))
end

function project_point(amb::AbstractArray, m::TSpaceManifold)
    return tangent2ambient(project_tv(amb, m.pt))
end

function project_point!(amb::AbstractArray, m::TSpaceManifold)
    return project_tv!(amb, m.pt)
end

function ambient2tangent(v::AbstractArray, p::TSpaceManifoldPt)
    return TSpaceManifoldTV(p, ambient2tangent(v, at_point(p.x)))
end

function project_tv(v::AbstractArray, p::TSpaceManifoldPt)
    return TSpaceManifoldTV(p, project_tv(v, at_point(p.x)))
end

function project_tv!(v::AbstractArray, p::AbstractArray, m::TSpaceManifold)
    project_tv!(v, point2ambient(m.pt), gettype(m.pt))
end

function tangent2ambient(v::TSpaceManifoldTV)
    return tangent2ambient(v.v)
end

function zero_tangent_vector(tspt::TSpaceManifoldPt)
    return TSpaceManifoldTV(tspt, zero_tangent_vector(at_point(tspt.x)))
end

function zero_tangent_vector!(v::BNBArray, p::AbstractArray, m::TSpaceManifold)
    zero_tangent_vector!(v, p, gettype(m.pt))
end

function geodesic_at(t::Number, x1::AbstractArray, x2::AbstractArray, m::TSpaceManifold)
    return (1-t).*x1 .+ t.*x2
end

function distance(x1::AbstractArray, x2::AbstractArray, m::TSpaceManifold)
    return norm(x1 - x2)
end

function exp(v::TSpaceManifoldTV)
    return TSpaceManifoldPt(at_point(v).x + v.v)
end

function exp!(p::TV, v::AbstractArray, at_pt::AbstractArray, m::TSpaceManifold) where TV<:BNBArray
    @condbc TV (p .= at_pt .+ v)
end

function log(x::TSpaceManifoldPt, y::TSpaceManifoldPt)
    return TSpaceManifoldTV(x, y.x - x.x)
end

function log!(v::TV, x::AbstractArray, y::AbstractArray, m::TSpaceManifold) where TV<:BNBArray
    @condbc TV (v .= y .- x)
end

function parallel_transport_geodesic(v::TSpaceManifoldTV, to_point::TSpaceManifoldPt)
    return TSpaceManifoldTV(to_point, v.v)
end

function parallel_transport_geodesic!(vout::TV, vin::AbstractArray, at_pt::AbstractArray, to_point::AbstractArray, m::TSpaceManifold) where TV<:BNBArray
    @condbc TV (vout .= vin)
end
