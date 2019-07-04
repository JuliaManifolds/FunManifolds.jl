
"""
    PowerSpace(m, n)

`n`th power manifold of manifold `m`.
"""
struct PowerSpace{M <: Manifold} <: Manifold
    m::M
    n::Int64
end

function manifold_dimension(m::PowerSpace{M}) where M <: Manifold
    return manifold_dimension(m.m)^m.n
end

function ambient_shape(m::PowerSpace{M}) where M <: Manifold
    return prod(ambient_shape(m.m))*m.n
end

"""
    PowerPt(xs)

Point on a power space, composed of points from the vector `xs`.
"""
struct PowerPt{P <: Point} <: Point
    xs::Vector{P}
end

function copyto!(x_to::PowerPt, x_from::PowerPt)
    copyto!(x_to.xs, x_from.xs)
    return x_to
end

function deepcopy(x::PowerPt)
    return PowerPt([deepcopy(p) for p in x.xs])
end

function gettype(x::PowerPt)
    xLen = length(x.xs)
    if xLen == 0
        return PowerSpace(Union{}, 0)
    else
        return PowerSpace(gettype(x.xs[1]), xLen)
    end
end

function isapprox(x1::PowerPt, x2::PowerPt; atol = atoldefault(x1, x2), rtol = rtoldefault(x1, x2))
    if gettype(x1) != gettype(x2)
        return false
    end
    return all(isapprox(x1.xs[i], x2.xs[i], atol = atol, rtol = rtol) for i ∈ 1:gettype(x1).n)
end

"""
    PowerTV(x, vs)

Tangent vector to a power manifold from the tangent space at `x`,
where `vs` is a vector of tangent vectors from respective subspaces
tangent at `x`.
"""
struct PowerTV{TV <: TangentVector, PPt <: PowerPt} <: TangentVector
    at_pt::PPt
    vs::Vector{TV}
end

function copyto!(v_to::PowerTV, v_from::PowerTV)
    copyto!(v_to.at_pt, v_from.at_pt)
    copyto!(v_to.vs, v_from.vs)
    return v_to
end

function deepcopy(v::PowerTV)
    return PowerTV(deepcopy(v.at_pt), [deepcopy(tv) for tv ∈ v.vs])
end

function zero_tangent_vector(pt::PowerPt)
    return PowerTV(pt, [zero_tangent_vector(subpt) for subpt ∈ pt.xs])
end

function zero_tangent_vector!(m::PowerSpace, v::BNBArray, at_pt::AbstractArray)
    for i ∈ 1:m.n
        zero_tangent_vector!(m.m, v[i], at_pt[i])
    end
end

function isapprox(v1::PowerTV, v2::PowerTV; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    if !(isapprox(at_point(v1), at_point(v2), atol = atol, rtol = rtol))
        return false
    end
    return all(isapprox(v1.vs[i], v2.vs[i], atol = atol, rtol = rtol) for i ∈ 1:length(v1.vs))
end

function +(v1::PowerTV, v2::PowerTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return PowerTV(v1.at_pt, v1.vs .+ v2.vs)
end

function add_vec!(v1::PowerTV, v2::PowerTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    add_vec!.(v1.vs, v2.vs)
    return v1
end

function add_vec!(m::PowerSpace, v1::BNBArray, v2::BNBArray, at_pt::AbstractArray)
    for i in 1:m.n
        add_vec!(m.m, v1[i], v2[i], at_pt[i])
    end
    return v1
end

function -(v1::PowerTV, v2::PowerTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return PowerTV(v1.at_pt, v1.vs .- v2.vs)
end

function sub_vec!(v1::PowerTV, v2::PowerTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    sub_vec!.(v1.vs, v2.vs)
    return v1
end

function sub_vec!(m::PowerSpace, v1::BNBArray, v2::BNBArray, at_pt::AbstractArray)
    for i in 1:m.n
        sub_vec!(m.m, v1[i], v2[i], at_pt[i])
    end
    return v1
end

function *(α::Real, v::PowerTV)
    return PowerTV(v.at_pt, α .* v.vs)
end

function mul_vec!(v::PowerTV, α::Real)
    mul_vec!.(v.vs, α)
    return v
end

function mul_vec!(m::PowerSpace, v::BNBArray, α::Real, at_pt::AbstractArray)
    for i in 1:m.n
        mul_vec!(m.m, v[i], α, at_pt[i])
    end
    return v
end

function point2ambient(p::PowerPt)
    return [point2ambient(x) for x ∈ p.xs]
    #return vcat((point2ambient(x) for x ∈ p.xs)...)
end

function ambient2point(m::PowerSpace, amb::AbstractArray)
    return PowerPt([ambient2point(m.m, mi) for mi ∈ amb])
end

function project_point(m::PowerSpace, amb::AbstractArray)
    return [project_point(m.m, mi) for mi ∈ amb]
end

function project_point!(m::PowerSpace, amb::TP) where TP<:BNBArray
    if TP <: NoBroadcastArray
        amb[] = project_point(m, amb[])
    else
        for mi ∈ amb
            project_point!(m.m, mi)
        end
    end
    return amb
end

function ambient2tangent(amb::AbstractArray, p::PowerPt)
    return PowerTV(p, [ambient2tangent(amb[i], p.xs[i]) for i ∈ 1:length(amb)])
end

function project_tv(amb::AbstractArray, p::PowerPt)
    return PowerTV(p, [project_tv(amb[i], p.xs[i]) for i ∈ 1:length(amb)])
end

function project_tv!(m::PowerSpace, v::TV, p::AbstractArray) where TV<:BNBArray
    for i ∈ 1:length(v)
        project_tv!(m.m, v[i], p[i])
    end
end

function tangent2ambient(v::PowerTV)
    return [tangent2ambient(x) for x ∈ v.vs]
end

function inner_amb(x1::PowerPt, x2::PowerPt)
    return sum(inner_amb(x1.xs[i], x2.xs[i]) for i ∈ 1:length(x1.xs))
end

function inner(m::PowerSpace, p::AbstractArray, v1::AbstractArray, v2::AbstractArray)
    return sum(inner(m.m, p[i], v1[i], v2[i]) for i ∈ 1:m.n)
end

function geodesic(x1::PowerPt, x2::PowerPt)
    num = length(x1.xs)
    gs = [geodesic(x1.xs[i], x2.xs[i]) for i ∈ 1:num]
    return CurvePt(gettype(x1)) do t
        return PowerPt([gs[i](t) for i ∈ 1:num])
    end
end

function geodesic_at(m::PowerSpace, t::Number, x1::AbstractArray, x2::AbstractArray)
    return [geodesic_at(m.m, t, x1[i], x2[i]) for i ∈ 1:m.n]
end

function distance(x1::PowerPt, x2::PowerPt)
    return norm([distance(x1.xs[i], x2.xs[i]) for i in 1:length(x1.xs)])
end

function distance(x1::AbstractArray, x2::AbstractArray, m::PowerSpace)
    return norm([distance(x1[i], x2[i], m.m) for i in 1:m.n])
end

function exp(v::PowerTV)
    return PowerPt([exp(v) for v ∈ v.vs])
end

function exp!(m::PowerSpace, p::BNBArray, at_pt::AbstractArray, v::AbstractArray)
    for i ∈ 1:m.n
        exp!(m.m, p[i], at_pt[i], v[i])
    end
    return p
end

function log(x::PowerPt, y::PowerPt)
    return PowerTV(x, [log(x.xs[i], y.xs[i]) for i ∈ 1:length(x.xs)])
end

function log!(m::PowerSpace, tv::BNBArray, x::AbstractArray, y::AbstractArray)
    for i ∈ 1:m.n
        log!(m.m, tv[i], x[i], y[i])
    end
    return tv
end

function parallel_transport_geodesic(v::PowerTV, to_point::PowerPt)
    return PowerTV(to_point, [parallel_transport_geodesic(v.vs[i], to_point.xs[i]) for i ∈ 1:length(v.vs)])
end

function parallel_transport_geodesic!(m::PowerSpace, vout::BNBArray, at_pt::AbstractArray, vin::AbstractArray, to_point::AbstractArray)
    for i ∈ 1:m.n
        parallel_transport_geodesic!(m.m, vout[i], at_pt[i], vin[i], to_point[i])
    end
    return vout
end
