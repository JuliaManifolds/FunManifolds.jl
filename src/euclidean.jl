"""
    EuclideanSpace(n)

Euclidean space ``\\mathbb{R}^{n}``.

# Arguments
- `n::Int`: dimension of the space.
"""
struct EuclideanSpace <: Manifold
    dim::Int
end

function ambient_shape(m::EuclideanSpace)
    return m.dim
end

"""
    EuclideanPt(x)

Point on an ``n``-dimensional euclidean space ``\\mathbb{R}^{n}``

# Arguments
- `x::Vector`: vector of coordinates of length ``n``.
"""
struct EuclideanPt{TVector <: AbstractVector{<:Real}} <: Point
    x::TVector
end

function copyto!(x_to::EuclideanPt, x_from::EuclideanPt)
    copyto!(x_to.x, x_from.x)
    return x_to
end

function deepcopy(x::EuclideanPt)
    return EuclideanPt(copy(x.x))
end

function gettype(x::EuclideanPt)
    return EuclideanSpace(length(x.x))
end

function isapprox(x1::EuclideanPt, x2::EuclideanPt; atol = atoldefault(x1, x2), rtol = rtoldefault(x1, x2))
    return isapprox(x1.x, x2.x, atol = atol, rtol = rtol)
end

"""
    EuclideanTV(x, v)

Tangent vector from an ``n``-dimensional euclidean space ``\\mathbb{R}^{n}``
from the tangent space at point `x` with ambient space representation `v`.
"""
struct EuclideanTV{TVectorPt <: AbstractVector{<:Real}, TVectorV <: Union{Array{<:Real,1}, MVector{N, <:Real} where N}} <: TangentVector
    at_pt::EuclideanPt{TVectorPt}
    v::TVectorV
end

function fastwarp!(vDest::EuclideanTV, a1::Real, v1::EuclideanTV,
    a2::Real, v2::EuclideanTV)
    vDest.v .= a1 .* v1.v .- a2 .* v2.v
    return vDest
end

function copyto!(v_to::EuclideanTV, v_from::EuclideanTV)
    copyto!(v_to.at_pt, v_from.at_pt)
    copyto!(v_to.v, v_from.v)
    return v_to
end

function deepcopy(v::EuclideanTV)
    return EuclideanTV(deepcopy(v.at_pt), copy(v.v))
end

function zero_tangent_vector(pt::EuclideanPt)
    return EuclideanTV(pt, _ensure_mutable(zero(pt.x)))
end

function zero_tangent_vector!(m::EuclideanSpace, v::TV, at_pt::AbstractArray) where TV<:BNBArray
    @condbc TV (v .= zero(at_pt))
end

function isapprox(v1::EuclideanTV, v2::EuclideanTV; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    DEBUG && if !(isapprox(v1.at_pt, v2.at_pt, atol = atol, rtol = rtol))
        return false
    end
    return isapprox(v1.v, v2.v, atol = atol, rtol = rtol)
end

function +(v1::EuclideanTV, v2::EuclideanTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return EuclideanTV(at_point(v1), v1.v + v2.v)
end

function add_vec!(v1::EuclideanTV, v2::EuclideanTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    v1.v .+= v2.v
    return v1
end

@inline function add_vec!(m::EuclideanSpace, v1::TV, v2::AbstractArray, at_pt::AbstractArray) where TV<:BNBArray
    @condbc TV (v1 .+= v2)
    return v1
end

function -(v1::EuclideanTV, v2::EuclideanTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return EuclideanTV(at_point(v1), _ensure_mutable(v1.v - v2.v))
end

function sub_vec!(v1::EuclideanTV, v2::EuclideanTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    v1.v .-= v2.v
    return v1
end

@inline function sub_vec!(m::EuclideanSpace, v1::TV, v2::AbstractArray, at_pt::AbstractArray) where TV<:BNBArray
    @condbc TV (v1 .-= v2)
    return v1
end

function *(α::Real, v::EuclideanTV)
    return EuclideanTV(at_point(v), α * v.v)
end

@inline function mul_vec!(m::EuclideanSpace, v::TV, α::Real, at_pt::AbstractArray) where TV<:BNBArray
    @condbc TV (v .*= α)
    return v
end

function point2ambient(p::EuclideanPt)
    return p.x
end

function ambient2point(m::EuclideanSpace, amb::AbstractVector{<:Real})
    return EuclideanPt(amb)
end

function project_point(m::EuclideanSpace, amb::AbstractVector{<:Real})
    return amb
end

function project_point!(m::EuclideanSpace, amb::BNBArray)
    return amb
end

function ambient2tangent(v::AbstractVector{<:Real}, p::EuclideanPt)
    return EuclideanTV(p, _ensure_mutable(v))
end

function project_tangent!(::EuclideanSpace, v::TV, p::AbstractVector) where TV<:BNBArray
    return v
end

function tangent2ambient(v::EuclideanTV)
    return v.v
end

function inner(m::EuclideanSpace, p::AbstractVector, v1::AbstractVector, v2::AbstractVector)
    return dot(v1, v2)
end

function geodesic_at(m::EuclideanSpace, t::Number, x1::AbstractVector, x2::AbstractVector)
    return (1-t)*x1 + t*x2
end

function distance(x1::AbstractArray, x2::AbstractArray, m::EuclideanSpace)
    return norm(x1 - x2)
end

function exp!(m::EuclideanSpace, p::TP, at_pt::AbstractVector, v::AbstractVector) where TP<:BNBArray
    @condbc TP (p .= at_pt .+ v)
    return p
end

function log!(m::EuclideanSpace, v::TV, x::AbstractVector, y::AbstractVector) where TV<:BNBArray
    @condbc TV (v .= y .- x)
    return v
end

function parallel_transport_geodesic!(m::EuclideanSpace, vout::TV, at_pt::AbstractArray, vin::AbstractArray, to_point::AbstractArray) where TV<:BNBArray
    @condbc TV (vout .= vin)
end
