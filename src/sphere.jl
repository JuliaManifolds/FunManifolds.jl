
@doc doc"""
    Sphere(n)

An `n`-dimensional unit sphere embedded in $\mathbb{R}^{n+1}$.
"""
struct Sphere <: Manifold
    dim::Int64
end

function ambient_shape(m::Sphere)
    return m.dim+1
end

"""
    SpherePt(x)

Point on a sphere embedded represented in ambient space by vector `x`.
"""
struct SpherePt{TVector <: AbstractArray{<:Real, 1}} <: Point
    x::TVector
end

function copyto!(x_to::SpherePt, x_from::SpherePt)
    copyto!(x_to.x, x_from.x)
    return x_to
end

function deepcopy(x::SpherePt)
    return SpherePt(copy(x.x))
end

function gettype(x::SpherePt)
    return Sphere(length(x.x)-1)
end

function isapprox(p1::SpherePt, p2::SpherePt; atol = atoldefault(p1, p2), rtol = rtoldefault(p1, p2))
    return isapprox(p1.x, p2.x, atol = atol, rtol = rtol)
end

"""
    SphereTV(x, v)

Tangent vector to a sphere from tangent space at point `x` with ambient space
representation `v`.
"""
struct SphereTV{TVectorPt <: AbstractArray{<:Real, 1}, TVectorV <: AbstractArray{<:Real, 1}} <: TangentVector
    at_pt::SpherePt{TVectorPt}
    v::TVectorV
end

function copyto!(v_to::SphereTV, v_from::SphereTV)
    copyto!(v_to.at_pt, v_from.at_pt)
    copyto!(v_to.v, v_from.v)
    return v_to
end

function deepcopy(v::SphereTV)
    return SphereTV(deepcopy(v.at_pt), copy(v.v))
end

function isapprox(v1::SphereTV, v2::SphereTV; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    if !(isapprox(v1.at_pt, v2.at_pt, atol = atol, rtol = rtol))
        return false
    end
    return isapprox(v1.v, v2.v, atol = atol, rtol = rtol)
end

function +(v1::SphereTV, v2::SphereTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't add tangent vectors from different tangent spaces $(at_point(v1)) and $(at_point(v2)).")
    end
    return SphereTV(at_point(v1), v1.v + v2.v)
end

@inline function add_vec!(m::Sphere, v1::TV, v2::AbstractArray, at_pt::AbstractArray) where TV<:BNBArray
    @condbc TV (v1 .+= v2)
    return v1
end

function -(v1::SphereTV, v2::SphereTV)
    if !(at_point(v1) ≈ at_point(v2))
        error("Can't subtract tangent vectors from different tangent spaces $(at_point(v1)) and $(at_point(v2)).")
    end
    return SphereTV(at_point(v1), v1.v - v2.v)
end

@inline function sub_vec!(m::Sphere, v1::TV, v2::AbstractArray, at_pt::AbstractArray) where TV<:BNBArray
    @condbc TV (v1 .-= v2)
    return v1
end

function *(α::Real, v::SphereTV)
    return SphereTV(at_point(v), α * v.v)
end

function mul_vec!(m::Sphere, v::TV, α::Real, at_pt::AbstractArray) where TV<:BNBArray
    @condbc TV (v .*= α)
    return v
end

function zero_tangent_vector(pt::SpherePt)
    return SphereTV(pt, _ensure_mutable(zero(pt.x)))
end

function zero_tangent_vector!(m::Sphere, v::TV, at_pt::AbstractArray) where TV<:BNBArray
    @condbc TV (v .= zero(at_pt))
end

function inner(v1::AbstractArray, v2::AbstractArray, p::AbstractArray, m::Sphere)
    return dot(v1, v2)
end

function geodesic_at(t::Number, x1::AbstractArray, x2::AbstractArray, m::Sphere)
    η = distance(x1, x2, m)
    if η ≈ 0.0
        return x1
    else
        return 1/sin(η).*(sin(η*(1-t)).*x1 .+ sin(η*t).*x2)
    end
end

function distance(x1::AbstractVector, x2::AbstractVector, ::Sphere)
    # in some rare cases due rounding errors dot(...) may be slightly outside the [-1,1] interval
    # and the acos function doesn't like it -- and so it's clamped
    dotval = clamp(dot(x1, x2), -1.0, 1.0)
    return acos(dotval)
end

function point2ambient(p::SpherePt)
    return p.x
end

function ambient2point(m::Sphere, amb::AbstractVector{<:Real})
    return SpherePt(amb)
end

function project_point(m::Sphere, amb::AbstractVector{<:Real})
    return amb / norm(amb)
end

function project_point!(m::Sphere, amb::TP) where TP<:BNBArray
    @condbc TP (amb .= amb ./ norm(amb))
    return amb
end

function ambient2tangent(v::AbstractVector{<:Real}, p::SpherePt)
    return SphereTV(p, _ensure_mutable(v))
end

function project_tv!(::Sphere, v::TV, p::AbstractVector) where TV<:BNBArray
    @condbc TV (v .-= (dot(p, v)/dot(p, p)).*p)
    return v
end

function tangent2ambient(v::SphereTV)
    return v.v
end

function exp!(m::Sphere, p::TV, at_pt::AbstractArray, v::AbstractArray) where TV<:BNBArray
    nv = norm(v)
    if nv ≈ 0.0
        @condbc TV (p .= at_pt)
    else
        @condbc TV (p .= cos(nv).*at_pt .+ (sin(nv)/nv).*v)
    end
end

@inline function log!(tv::TV, x::AbstractArray, y::AbstractArray, m::Sphere) where TV<:BNBArray
    θ = acos(dot(x, y))
    if θ ≈ 0.0
        zero_tangent_vector!(m, tv, x)
    else
        @condbc TV (tv .= (θ/sin(θ)) .* (y .- cos(θ).*x))
    end
    tv
end

function parallel_transport_geodesic!(vout::TV, vin::AbstractArray, at_pt::AbstractArray, to_point::AbstractArray, m::Sphere) where TV<:BNBArray
    factor = 2*dot(vin, to_point)/norm(at_pt + to_point)^2
    @condbc TV (vout .= vin .- factor.*(at_pt .+ to_point))
    return vout
end
