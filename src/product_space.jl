
"""
    ProductSpace((m1, m2, ..., mN))

Product manifold of manifolds `m1`, `m2`, ..., `mN`.
"""
struct ProductSpace{Ms <: Tuple} <: Manifold
    ms::Ms
end

(==)(s1::ProductSpace, s2::ProductSpace) = s1.ms == s2.ms

function manifold_dimension(x::ProductSpace)
    return sum(manifold_dimension(mi) for mi in x.ms)
end

function ambient_shape(m::ProductSpace)
    return sum(dim_ambient(mi) for mi in m.ms)
end

"""
    ProductPt((x1, x2, ..., xN), m)

Point on a product space composed of points `x1`,  `x2`, ..., `xN` in
respective multiplied manifolds represented by product space `m`.
"""
struct ProductPt{Pts<:NTuple{N,AbstractArray} where N, TPS<:ProductSpace} <: Point
    xs::Pts
    m::TPS
end

function ProductPt(x::NTuple{N,Point}) where N
    xs = map(x) do xi
        return point2ambient(xi)
    end
    m = ProductSpace(map(x) do xi
        return gettype(xi)
    end)
    return ProductPt(xs, m)
end

function deepcopy(x::ProductPt)
    return ProductPt(map(p -> deepcopy(p), x.xs), x.m)
end

function gettype(x::ProductPt)
    return x.m
end

function isapprox(x1::ProductPt, x2::ProductPt; atol = atoldefault(x1, x2), rtol = rtoldefault(x1, x2))
    return length(x1.xs) == length(x2.xs) &&
        all(isapprox(
            ambient2point(x1.m.ms[i], x1.xs[i]),
            ambient2point(x2.m.ms[i], x2.xs[i]),
            atol = atol, rtol = rtol) for i in 1:length(x1.xs))
end

"""
    ProductTV(pt, (v1, v2, ..., vN))

Tangent vector from tangent space at point `pt` to a product manifold.
Vectors in respective multiplied manifolds are given by `v1`, `v2`, ..., `vN`.
"""
struct ProductTV{PPt <: ProductPt, TVs <: NTuple{N, AbstractArray} where N} <: TangentVector
    at_pt::PPt
    vs::TVs
end

function deepcopy(v::ProductTV)
    return ProductTV(deepcopy(v.at_pt), map(tv -> deepcopy(tv), v.vs))
end

function zero_tangent_vector(pt::ProductPt)
    exs = enumeratetuple(pt.xs)
    tv = map(exs) do tup
        i = tup[1]
        xi = tup[2]
        mi = pt.m.ms[i]
        return zero_tangent_vector(xi, mi)
    end
    return ProductTV(pt, tv)
end

function zero_tangent_vector!(v::BNBArray, at_pt::AbstractArray, m::ProductSpace)
    exs = enumeratetuple(m.ms)
    map(exs) do tup
        i = tup[1]
        mi = tup[2]
        zero_tangent_vector!(v[i], at_pt[i], mi)
    end
end

function isapprox(v1::ProductTV, v2::ProductTV; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    if !(isapprox(v1.at_pt, v2.at_pt, atol = atol, rtol = rtol))
        return false
    end
    return all(isapprox(v1.vs[i], v2.vs[i], atol = atol, rtol = rtol) for i in 1:length(v1.vs))
end

function +(v1::ProductTV, v2::ProductTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return ProductTV(v1.at_pt, Tuple(v1.vs[i] + v2.vs[i] for i in 1:length(v1.vs)))
end

function add_vec!(v1::ProductTV, v2::ProductTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    add_vec!.(v1.vs, v2.vs, v1.at_pt.xs, v1.at_pt.m.ms)
    return v1
end

function add_vec!(v1::BNBArray, v2::BNBArray, at_pt::AbstractArray, m::ProductSpace)
    for i in 1:length(m.ms)
        add_vec!(v1[i], v2[i], at_pt[i], m.ms[i])
    end
    return v1
end

function -(v1::ProductTV, v2::ProductTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return ProductTV(v1.at_pt, Tuple(v1.vs[i] - v2.vs[i] for i in 1:length(v1.vs)))
end

function sub_vec!(v1::ProductTV, v2::ProductTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    sub_vec!.(v1.vs, v2.vs, v1.at_pt.xs, v1.at_pt.m.ms)
    return v1
end

function sub_vec!(v1::BNBArray, v2::BNBArray, at_pt::AbstractArray, m::ProductSpace)
    for i in 1:length(m.ms)
        sub_vec!(v1[i], v2[i], at_pt[i], m.ms[i])
    end
    return v1
end

function *(α::T, v::ProductTV) where T <: Real
    return ProductTV(v.at_pt, Tuple(α * vi for vi in v.vs))
end

function mul_vec!(v::ProductTV, α::Real)
    mul_vec!.(v.vs, α, v.at_pt.xs, v.at_pt.m.ms)
    return v
end

function mul_vec!(v::BNBArray, α::Real, at_pt::AbstractArray, m::ProductSpace)
    for i in 1:length(m.ms)
        mul_vec!(v[i], α, at_pt[i], m.ms[i])
    end
    return v
end

function point2ambient(p::ProductPt)
    return TupleArray(p.xs)
end

function ambient2point(m::ProductSpace, amb::AbstractArray)
    return ProductPt(amb.a, m)
end

function project_point(m::ProductSpace, amb::AbstractVector)
    return TupleArray(map(i -> project_point(i[1], i[2]), ziptuples(m.ms, amb.a)))
end

function project_point!(m::ProductSpace, amb::TP) where TP<:BNBArray
    if TP <: TupleArray
        map(i -> project_point!(i[1], i[2]), ziptuples(m.ms, amb.a))
    else
        amb[] = project_point(m, amb[])
    end
    return amb
end

function ambient2tangent(amb::AbstractVector, p::ProductPt)
    return ProductTV(p, amb.a)
end

function project_tv(amb::TupleArray, p::ProductPt)
    tv = map(i -> project_tv(i[1], i[2], i[3]), ziptuples(p.m.ms, amb.a, p.xs))
    return ProductTV(p, tv)
end

function project_tv!(m::ProductSpace, v::AbstractVector, p::AbstractArray)
    exs = enumeratetuple(m.ms)
    map(exs) do tup
        i = tup[1]
        mi = tup[2]
        project_tv!(mi, v[i], p[i])
    end
end

function tangent2ambient(v::ProductTV)
    return TupleArray(v.vs)
end

function inner(v1::ProductTV, v2::ProductTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    ms = v1.at_pt.m.ms
    return sum(inner(v1.vs[i], v2.vs[i], v1.at_pt.xs[i], ms[i]) for i in 1:length(v1.vs))
end

function inner(v1::AbstractArray, v2::AbstractArray, p::AbstractArray, m::ProductSpace)
    return sum(inner(v1[i], v2[i], p[i], m.ms[i]) for i in 1:length(m.ms))
end

function geodesic(x1::ProductPt, x2::ProductPt)
    gs = map(Tuple(1:length(x1.xs))) do i
        return geodesic(x1.xs[i], x2.xs[i], x1.m.ms[i])
    end
    return CurvePt(gettype(x1)) do t
        return ProductPt(map(g -> point2ambient(g(t)), gs), x1.m)
    end
end

function geodesic_at(t::Number, x1::AbstractArray, x2::AbstractArray, m::ProductSpace)
    gs = map(i -> geodesic_at(t, x1[i], x2[i], m.ms[i]), Tuple(1:length(m.ms)))
    return TupleArray(gs)
end

function distance(x1::ProductPt, x2::ProductPt)
    return sqrt(sum(distance(x1.xs[i], x2.xs[i], x1.m.ms[i])^2 for i in 1:length(x1.xs)))
end

function distance(x1::TupleArray, x2::TupleArray, m::ProductSpace)
    return sqrt(sum(distance(x1[i], x2[i], m.ms[i])^2 for i in 1:length(m.ms)))
end

function exp(v::ProductTV)
    x = map(ziptuples(v.vs, v.at_pt.xs, v.at_pt.m.ms)) do t
        exp(t[1], t[2], t[3])
    end
    return ProductPt(x, v.at_pt.m)
end

function exp!(p::TupleArray, v::TupleArray, at_pt::TupleArray, m::ProductSpace)
    for i ∈ 1:length(m.ms)
        exp!(p[i], v[i], at_pt[i], m.ms[i])
    end
    return p
end

function log(x::ProductPt, y::ProductPt)
    args = ziptuples(x.xs, y.xs, x.m.ms)
    tv = map(args) do a
        return log(a[1], a[2], a[3])
    end
    return ProductTV(x, tv)
end

function log!(v::TupleArray, x::TupleArray, y::TupleArray, m::ProductSpace)
    for i ∈ 1:length(m.ms)
        log!(v[i], x[i], y[i], m.ms[i])
    end
end

function parallel_transport_geodesic(v::ProductTV, to_point::ProductPt)
    ms = v.at_pt.m.ms
    tv = map(tuple(1:length(v.vs)...)) do i
        parallel_transport_geodesic(v.vs[i], v.at_pt.xs[i], to_point.xs[i], ms[i])
    end
    return ProductTV(to_point, tv)
end

function parallel_transport_geodesic!(vout::BNBArray, vin::AbstractArray, at_pt::AbstractArray, to_point::AbstractArray, m::ProductSpace)
    for i ∈ 1:length(m.ms)
        parallel_transport_geodesic!(vout[i], vin[i], at_pt[i], to_point[i], m.ms[i])
    end
end

function inner_amb(x1::AbstractArray, x2::AbstractArray, m::ProductSpace)
    return sum(enumeratetuple(m.ms)) do (i, m)
        inner_amb(x1[i], x2[i], m)
    end
end
