
"""
    Manifold

An abstract manifold.
"""
abstract type Manifold
end

"""
    AbstractCurveSpace

An abstract manifold of curves.
"""
abstract type AbstractCurveSpace <: Manifold
end

"""
    values_in(cs)

Returns space of values of curves from a given curve space `cs`.
"""
function values_in(cs::AbstractCurveSpace)
    error("Function values_in is not yet defined for type $(typeof(cs)).")
end

"""
    Point

Point on a manifold.
"""
abstract type Point
end

"""
    AbstractCurvePt

Point on a manifolds of curves (i.e. a curve).
"""
abstract type AbstractCurvePt <: Point
end

"""
    paramgrid(x::AbstractCurvePt)

Discretization grid of the given curve `x`.
"""
function paramgrid(x::AbstractCurvePt)
    error("Function paramgrid is not yet defined for type $(typeof(x)).")
end

"""
    gettype(x)

Extract type of a manifold `x` belongs to.
"""
function gettype(x::Point)
    error("Function gettype is not yet defined for type $(typeof(x)).")
end

function rtoldefault(x1::Point, x2::Point)
    return 1.e-8
end

function atoldefault(x1::Point, x2::Point)
    return 0.0
end

"""
    TangentVector

Vector tangent to a manifold at a given point.
"""
abstract type TangentVector
end

function rtoldefault(x1::TangentVector, x2::TangentVector)
    return 1.e-8
end

function atoldefault(x1::TangentVector, x2::TangentVector)
    return 0.0
end


function concretize_tols(x1::Union{TangentVector,Point}, x2::Union{TangentVector,Point};
    reltol::Union{Real,Nothing}=nothing, abstol::Union{Real,Nothing}=nothing)

    rtol = if reltol === nothing
        rtoldefault(x1, x2)
    else
        reltol
    end

    atol = if abstol === nothing
        atoldefault(x1, x2)
    else
        abstol
    end

    return (rtol, atol)
end

"""
    point2ambient(p)

Convert `p` to representation in the ambient space.
"""
function point2ambient(p::Point)
    error("Function point2ambient is not yet defined for type $(typeof(p)).")
end

"""
    ambient2point(amb, m)

Convert `amb` (representation of point in the ambient space) to a point in
manifold `m`. If `amb` does not represent a valid point on `m`, the function
either fails (in debug mode) or returns an incorrect point.
"""
function ambient2point(amb::AbstractArray, m::Manifold)
    error("Function ambient2point is not yet defined for types $(typeof(amb)), $(typeof(m)).")
end

"""
    project_point(amb::AbstractArray, m::Manifold)

The function find a point from manifold `m` closest to `amb` in the
ambient space.
"""
function project_point(amb::AbstractArray, m::Manifold)
    error("Function project_point is not yet defined for types $(typeof(amb)), $(typeof(m)).")
end

"""
    project_point_wrapped(amb::AbstractArray, m::Manifold)

Convert `amb` (representation of point in the ambient space) to a point in
manifold `m`. If `amb` does not represent a valid point on `m`, the function
tries to find a valid one closest to `amb` in the ambient space.
"""
function project_point_wrapped(amb::AbstractArray, m::Manifold)
    return ambient2point(project_point(amb, m), m)
end

"""
    project_point!(amb::BNBArray, m::Manifold)

Project `amb` (representation of point in the ambient space) to a point in
manifold `m` tha is the closest point to `amb` in the ambient space.
"""
function project_point!(amb::BNBArray, m::Manifold)
    error("Function project_point! is not yet defined for types $(typeof(amb)), $(typeof(m)).")
end

"""
    ambient2tangent(v, p)

Make a tangent vector at point `p` with ambient space representation `v`. Does
not project the given ambient space representation.
"""
function ambient2tangent(v::AbstractArray, p::Point)
    error("Function ambient2tangent is not yet defined for types $(typeof(v)), $(typeof(p)).")
end

"""
    project_tv(v::AbstractArray, p::Point)

Make a tangent vector at point `p` with ambient space representation `v`.
If `v` does not lie in the tangent space, it is appropriately projected.
"""
function project_tv(v::AbstractArray, p::Point)
    tvamb = project_tv(v, point2ambient(p), gettype(p))
    return ambient2tangent(tvamb, p)
end

"""
    project_tv(v::AbstractArray, p::AbstractArray, m::Manifold)

Make a tangent vector at point `p` from manifold `m` with ambient space
representation `v`. If `v` does not lie in the tangent space, it is
appropriately projected.
"""
function project_tv(v::TV, p::AbstractArray, m::Manifold) where TV<:AbstractArray
    if TV <: SArray
        v2 = @MVector [deepcopy(v)]
        project_tv!(view(v2, 1), p, m)
        return v2[1]
    else
        v2 = deepcopy(v)
        project_tv!(v2, p, m)
        return v2
    end
end

"""
    project_tv!(v, p::Point)

Project vector `v` in-place to a tangent vector at point `p`.
"""
function project_tv!(v::BNBArray, p::Point)
    project_tv!(v, point2ambient(p), gettype(p))
end

"""
    project_tv!(v, p, m::Manifold)

Project vector `v` in-place to a tangent vector at point with ambient space
representation `p` on manifold `m`.
"""
function project_tv!(v::BNBArray, p::AbstractArray, m::Manifold)
    error("Function project_tv! is not yet defined for types $(typeof(v)), $(typeof(p)) and $(typeof(m)).")
end

"""
    tangent2ambient(v)

Get the ambient space representation of tangent vector `v`.
"""
function tangent2ambient(v::TangentVector)
    error("Function tangent2ambient is not yet defined for type $(typeof(v)).")
end

"""
    zero_tv(pt)

Get the zero tangent vector at point `pt`.
"""
function zero_tv(at_pt::Point)
    error("Function zero_tv is not yet defined for type $(typeof(at_pt)).")
end

"""
    zero_tv(p, m)

Get the zero tangent vector at point `p` from manifold `m`.
"""
function zero_tv(p::AbstractArray, m::Manifold)
    return tangent2ambient(zero_tv(ambient2point(p, m)))
end

"""
    zero_tv!(v, pt)

Set `v`  to the zero tangent vector at point `pt`.
"""
function zero_tv!(v::BNBArray, at_pt::Point)
    zero_tv!(v, point2ambient(at_pt), gettype(at_pt))
end

"""
    zero_tv!(v, p, m)

Set `v`  to the zero tangent vector at point with ambient space
representation `p` from a manifold `m`.
"""
function zero_tv!(v::BNBArray, at_pt::AbstractArray, m::Manifold)
    error("Function zero_tv! is not yet defined for types $(typeof(v)), $(typeof(at_pt)) and $(typeof(m)).")
end

function +(v1::TangentVector, v2::TangentVector)
    if !(at_point(v1) ≈ at_point(v2))
        error("Can't add tangent vectors from different tangent spaces $(at_point(v1)) and $(at_point(v2)).")
    end
    return ambient2tangent(tangent2ambient(v1) + tangent2ambient(v2), at_point(v1))
end

"""
    add_vec!(v1, v2)

Add tangent vector `v2` to tangent vector `v1` in-place.
"""
@inline function add_vec!(v1::TangentVector, v2::TangentVector)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't add tangent vectors from different tangent spaces $(at_point(v1)) and $(at_point(v2)).")
    end
    add_vec!(tangent2ambient(v1), tangent2ambient(v2), at_point(v1))
end

"""
    add_vec!(v1, v2, at_pt)

Add tangent vector with ambient space representation `v2` to tangent vector
with ambient space representation `v1` in-place.
Both vectors at tangent at `at_pt`
"""
@inline function add_vec!(v1::BNBArray, v2::BNBArray, at_pt::Point)
    add_vec!(v1, v2, point2ambient(at_pt), gettype(at_pt))
end

"""
    add_vec!(v1, v2, at_pt, m)

Add tangent vector with ambient space representation `v2` to tangent vector
with ambient space representation `v1` in-place.
Both vectors at tangent at a point with ambient space representation `at_pt`
from a manifold `m`.
"""
function add_vec!(v1::BNBArray, v2::BNBArray, at_pt::AbstractArray, m::Manifold)
    error("Function add_vec! is not yet defined for for types $(typeof(v1)), $(typeof(v2)), $(typeof(at_pt)) and $(typeof(m)).")
end

"""
    add_vec(v1, v2, at_pt)

Add tangent vector `v2` to tangent vector `v1` making a new vector.
`v1` and `v2` are bare arrays representing vectors tangent at `at_pt`.
"""
function add_vec(v1::AbstractArray, v2::AbstractArray, at_pt::Point)
    v1c = deepcopy(v1)
    add_vec!(v1c, v2, at_pt)
    return v1c
end

@inline function -(v1::TangentVector, v2::TangentVector)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't subtract tangent vectors from different tangent spaces $(at_point(v1)) and $(at_point(v2)).")
    end
    return ambient2tangent(tangent2ambient(v1) - tangent2ambient(v2), at_point(v1))
end

"""
    sub_vec!(v1, v2)

Subtract tangent vector `v2` from tangent vector `v1` in-place.
"""
@inline function sub_vec!(v1::TangentVector, v2::TangentVector)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't subtract tangent vectors from different tangent spaces $(at_point(v1)) and $(at_point(v2)).")
    end
    sub_vec!(tangent2ambient(v1), tangent2ambient(v2), at_point(v1))
end

"""
    sub_vec!(v1, v2, at_pt)

Subtract tangent vector with ambient space representation `v2` from tangent
vector with ambient space representation `v1` in-place.
Both vectors at tangent at `at_pt`
"""
@inline function sub_vec!(v1::BNBArray, v2::BNBArray, at_pt::Point)
    sub_vec!(v1, v2, point2ambient(at_pt), gettype(at_pt))
end

"""
    sub_vec!(v1, v2, at_pt, m)

Subtract tangent vector with ambient space representation `v2` from tangent
vector with ambient space representation `v1` in-place.
Both vectors at tangent at a point with ambient space representation `at_pt`
from a manifold `m`.
"""
function sub_vec!(v1::BNBArray, v2::BNBArray, at_pt::AbstractArray, m::Manifold)
    error("Function sub_vec! is not yet defined for for types $(typeof(v1)), $(typeof(v2)), $(typeof(at_pt)) and $(typeof(m)).")
end

"""
    sub_vec(v1, v2, at_pt)

Subtract tangent vector `v2` to tangent vector `v1` making a new vector.
`v1` and `v2` are bare arrays representing vectors tangent at `at_pt`.
"""
function sub_vec(v1::AbstractArray, v2::AbstractArray, at_pt::Point)
    v1c = deepcopy(v1)
    sub_vec!(v1c, v2, at_pt)
    return v1c
end

@inline function *(α::Real, v::TangentVector)
    return ambient2tangent(α * tangent2ambient(v), at_point(v))
end

"""
    mul_vec!(v, α)

Multiply tangent vector `v` by number `α` in-place.
"""
@inline function mul_vec!(v::TangentVector, α::Real)
    mul_vec!(tangent2ambient(v), α, at_point(v))
end

"""
    mul_vec!(v, α, at_pt)

Multiply tangent vector `v` by `α` in-place. `v` is a bare array representing
a vector tangent at `at_pt`.
"""
@inline function mul_vec!(v::BNBArray, α::Real, at_pt::Point)
    mul_vec!(v, α, point2ambient(at_pt), gettype(at_pt))
end

"""
    mul_vec!(v, α, at_pt, m)

Multiply tangent vector `v` by `α` in-place. `v` is a bare array representing
a vector tangent at point with ambient space representation `at_pt`
from manifold `m`.
"""
function mul_vec!(v::BNBArray, α::Real, at_pt::AbstractArray, m::Manifold)
    error("Function mul_vec! is not yet defined for for types $(typeof(v)), $(typeof(α)), $(typeof(at_pt)) and $(typeof(m)).")
end

"""
    mul_vec(v, α, at_pt)

Multiply tangent vector `v` by `α` making a new vector.
`v` is a bare array representing a vector tangent at `at_pt`.
"""
function mul_vec(v::AbstractArray, α::Real, at_pt::Point)
    vc = deepcopy(v)
    mul_vec!(vc, α, at_pt)
    return vc
end

"""
    mul_vec(v, α, at_pt, m)

Multiply tangent vector `v` by `α` making a new vector.
`v` is a bare array representing a vector tangent at `at_pt` from manifold `m`.
"""
function mul_vec(v::TV, α::Real, at_pt::AbstractArray, m::Manifold) where TV<:AbstractArray

    if TV <: SArray
        v2 = @MVector [v]
        mul_vec!(view(v2, 1), α, at_pt, m)
        return v2[1]
    else
        v2 = deepcopy(v)
        mul_vec!(v2, α, at_pt, m)
        return v2
    end
end

"""
    fastwarp!(v_dest, a1, v1, a2, v2)

Sets vector part of `v_dest` to a1*v1 - a2*v2. Temporary optimization solution
until a proper expression rewriting is available.
"""
function fastwarp!(v_dest::TV, a1::Real, v1::TV,
    a2::Real, v2::TV) where TV <: TangentVector
    copyto!(v_dest, a1 * v1 - a2 * v2)
    return v_dest
end

"""
    at_point(v)

Return the point where the tangent vector `v` is attached.
"""
@inline function at_point(v::TV) where TV <: TangentVector
    return v.at_pt
end

"""
    dim(m)

Return dimension manifold `m`.
"""
function dim(x::Manifold)
    return x.dim
end

function dim(v::TangentVector)
    return dim(at_point(v))
end

function dim(p::Point)
    return dim(gettype(p))
end

"""
    ambient_shape(m)

Sizes of matrices representing point on manifold `m` in the ambient space.
"""
function ambient_shape(m::Manifold)
    error("Function ambient_shape is not yet defined for type $(typeof(m)).")
end

"""
    dim_ambient(m)

Dimension of the ambient space for manifold `m`.
"""
function dim_ambient(m::Manifold)
    return prod(ambient_shape(m))
end

function dim_ambient(p::Point)
    return dim_ambient(gettype(p))
end

function dim_ambient(tv::TangentVector)
    return dim_ambient(at_point(tv))
end


"""
    exp(v)

Compute exponential map of tangent vector `v`.
"""
function exp(v::TangentVector)
    p = at_point(v)
    m = gettype(p)
    return ambient2point(exp(tangent2ambient(v), p), m)
end

"""
    exp!(p, v, at_pt)

Compute exponential map of bare array `v`, a tangent vector at `at_pt`
and save to bare array `p` of the size and shape of `at_pt`.
"""
function exp!(p::BNBArray, v::AbstractArray, at_pt::Point)
    exp!(p, v, point2ambient(at_pt), gettype(at_pt))
end

"""
    exp!(p, v, at_pt, m)

Compute exponential map of bare array `v`, a tangent vector at a bare point
`at_pt` and save to bare array `p` of the size and shape of `at_pt`.
Underlying manifold is `m`.
"""
function exp!(p::BNBArray, v::AbstractArray, at_pt::AbstractArray, m::Manifold)
    error("Function exp! is not yet defined for types $(typeof(p)), $(typeof(v)), $(typeof(at_pt)) and $(typeof(m)).")
end

"""
    exp(v, at_pt)

Compute exponential map of bare array `v`, a tangent vector at `at_pt`,
returning a bare array.
"""
function exp(v::AbstractArray, at_pt::Point)
    p = similar_ambient(point2ambient(at_pt), v)
    exp!(p, v, at_pt)
    return p
end

"""
    exp(v, at_pt, m::Manifold)

Compute exponential map of bare array `v`, a tangent vector at `at_pt`,
on a manifold `m` returning a bare array.
"""
function exp(v::AbstractArray, at_pt::AbstractArray, m::Manifold)
    p = similar_ambient(at_pt, v)
    exp!(p, v, at_pt, m)
    return p
end

"""
    retract(v)

Computes an approximate exponential map of tangent vector `v`.
Good enough for small tangent vectors.
"""
function retract(v::TangentVector)
    return exp(v)
end

"""
    retract!(p, v, at_pt)

Compute exponential map of bare array `v`, a tangent vector at `at_pt`
and save to bare array `p` of the size and shape of `at_pt`.
"""
function retract!(p::BNBArray, v::AbstractArray, at_pt::Point)
    retract!(p, v, point2ambient(at_pt), gettype(at_pt))
end

"""
    retract!(p, v, at_pt, m)

Compute exponential map of bare array `v`, a tangent vector at a bare point
`at_pt` and save to bare array `p` of the size and shape of `at_pt`.
Underlying manifold is `m`.
"""
function retract!(p::BNBArray, v::AbstractArray, at_pt::AbstractArray, m::Manifold)
    exp!(p, v, at_pt, m)
end

"""
    retract(v, at_pt)

Compute exponential map of bare array `v`, a tangent vector at `at_pt`,
returning a bare array.
"""
function retract(v::AbstractArray, at_pt::Point)
    p = similar_ambient(point2ambient(at_pt), v)
    retract!(p, v, at_pt)
    return p
end

"""
    retract(v, at_pt, m::Manifold)

Compute exponential map of bare array `v`, a tangent vector at `at_pt`,
on a manifold `m` returning a bare array.
"""
function retract(v::AbstractArray, at_pt::AbstractArray, m::Manifold)
    p = similar_ambient(at_pt, v)
    retract!(p, v, at_pt, m)
    return p
end

"""
    inner(v1, v2)

Compute inner product of two given tangent vectors `v1` and `v2`.
For specifying tolerances use `PARAMS.quad_abs_tol` and `PARAMS.quad_rel_tol`.
"""
@inline function inner(v1::TangentVector, v2::TangentVector)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Given vectors are attached at different points $(at_point(v1)) and $(at_point(v2)).")
    end
    return inner(tangent2ambient(v1), tangent2ambient(v2), at_point(v1))
end

"""
    inner(v1, v2, p::Point)

Compute inner product of two given vectors tangent at `p`
with ambient representations `v1` and `v2`.
If integration is required then given tolerances
`PARAMS.quad_rel_tol` and `PARAMS.quad_abs_tol` are used.
"""
@inline function inner(v1::AbstractArray, v2::AbstractArray, p::Point)
    return inner(v1, v2, point2ambient(p), gettype(p))
end

"""
    inner(v1, v2, p, m::Manifold)

Compute inner product of two given vectors tangent at a point with tangent
space representation `p` from a manifold `m`,
where tangent vectors have ambient representations `v1` and `v2`.
If integration is required then given tolerances
`PARAMS.quad_rel_tol` and `PARAMS.quad_abs_tol` are used.
"""
function inner(v1::AbstractArray, v2::AbstractArray, p::AbstractArray, m::Manifold)
    error("Function inner is not yet defined for types $(typeof(v1)), $(typeof(v2)), $(typeof(p)) and $(typeof(m)).")
end

"""
    norm(v)

Norm of tangent vector `v`.
"""
function norm(v::TangentVector)
    return sqrt(inner(v, v))
end

"""
    norm(v, p::Point)

Norm of tangent vector at `p` with ambient space representation `v`.
"""
function norm(v::AbstractArray, p::Point)
    return sqrt(inner(v, v, p))
end

"""
    norm(v, p, m::Manifold)

Norm of tangent vector at `p` with ambient space representation `v`, on a
manifold `m`.
"""
function norm(v::AbstractArray, p::AbstractArray, m::Manifold)
    return sqrt(inner(v, v, p, m))
end

"""
    log(x, y)

Compute logarithmic map, that is a tangent vector at point `x`
pointing at the point `y`.
"""
function log(x::Point, y::Point)
    return ambient2tangent(log(point2ambient(x), point2ambient(y), gettype(x)), x)
end

"""
    log(m, x, y)

Compute logarithmic map, that is a tangent vector to manifold `m` at point
with ambient space representation `x` pointing at the point with ambient
space representation `y`.
Returns ambient space representation, as determined by input arguments.
"""
function log(x::AbstractArray, y::AbstractArray, m::Manifold)
    v = zero_tv(x, m)
    log!(v, x, y, m)
    return v
end

"""
    log!(tv, x, y, m)

Compute logarithmic map, that is a tangent vector to manifold `m` at point
with ambient space representation `x` pointing at the point with ambient
space representation `y`.
Store the result in `tv` (a reference or a zero-dimensional abstract array).
"""
function log!(tv::BNBArray, x::AbstractArray, y::AbstractArray, m::Manifold)
    error("Function log! is not yet defined for types $(typeof(tv)), $(typeof(x)), $(typeof(y)) and $(typeof(m)).")
end

"""
    parallel_transport_geodesic(v, to_point)

Compute parallel transport of tangent vector `v` along the shortest geodesic
to point `to_point`.
"""
function parallel_transport_geodesic(v::TangentVector, to_point::Point)
    vin = tangent2ambient(v)
    at_pt = point2ambient(at_point(v))
    vout = parallel_transport_geodesic(vin, at_pt, point2ambient(to_point), gettype(to_point))
    return ambient2tangent(vout, to_point)
end

"""
    parallel_transport_geodesic!(vout, vin, to_point)

Compute parallel transport of tangent vector `vin` along the shortest geodesic
to point `to_point` and save the result in `vout`.
"""
function parallel_transport_geodesic!(vout::TangentVector, vin::TangentVector, to_point::Point)
    parallel_transport_geodesic!(tangent2ambient(vout), tangent2ambient(vin), point2ambient(to_point), gettype(to_point))
end

"""
    parallel_transport_geodesic(vin, at_pt, to_point, m)

Compute parallel transport of tangent vector `vin` at `at_pt` along
the shortest geodesic to point `to_point` on manifold `m`.
"""
function parallel_transport_geodesic(vin::AbstractArray, at_pt::AbstractArray, to_point::AbstractArray, m::Manifold)
    vout = similar(vin)
    parallel_transport_geodesic!(vout, vin, at_pt, to_point, m)
    return vout
end

"""
    parallel_transport_geodesic!(vout, vin, at_pt, to_point, m)

Compute parallel transport of tangent vector `vin` attached at point `at_pt`
along the shortest geodesic to point `to_point` and save the result in `vout`.
All of this variables accept ambient space representations. Underlying manifold
is given in `m`.
"""
function parallel_transport_geodesic!(vout::BNBArray, vin::AbstractArray, at_pt::AbstractArray, to_point::AbstractArray, m::Manifold)
    error("Function parallel_transport_geodesic! is not yet defined for types $(typeof(vout)), $(typeof(vin)), $(typeof(at_pt)) $(typeof(to_point)) and $(typeof(m)).")
end

"""
    geodesic(x1, x2)

Compute geodesic on a manifold between `x1` and `x2`.
"""
function geodesic(x1::Point, x2::Point)
    return CurvePt(gettype(x1)) do t
        return geodesic_at(t, x1, x2)
    end
end

"""
    geodesic(x1::AbstractArray, x2::AbstractArray, m::Manifold)

Compute geodesic on a manifold `m` between `x1` and `x2` (bare vectors).
"""
function geodesic(x1::AbstractArray, x2::AbstractArray, m::Manifold)
    return CurvePt(m) do t
        return ambient2point(geodesic_at(t, x1, x2, m), m)
    end
end

"""
    geodesic_at(t, x1, x2)

Compute geodesic between `x1` and `x2` (points) at
coordinate `t` and return a point.
"""
function geodesic_at(t::Number, x1::Point, x2::Point)
    m = gettype(x1)
    return ambient2point(geodesic_at(t, point2ambient(x1), point2ambient(x2), m), m)
end

"""
    geodesic_at(t, x1, x2, m::Manifold)

Compute geodesic on a manifold `m` between `x1` and `x2` (bare vectors) at
coordinate `t` and return bare vector.
"""
function geodesic_at(t::Number, x1::AbstractArray, x2::AbstractArray, m::Manifold)
    error("Function geodesic_at is not yet defined for types $(typeof(t)), $(typeof(x1)), $(typeof(x2)) and $(typeof(m)).")
end

"""
    geodesic_distance(x1, x2)

Compute distance (along the shortest geodesic) between two given points
`x1` and `x2`. If integration is required then given tolerances
`PARAMS.quad_rel_tol` and `PARAMS.quad_abs_tol` are used.
"""
function geodesic_distance(x1::Point, x2::Point)
    return geodesic_distance(point2ambient(x1), point2ambient(x2), gettype(x1))
end

"""
    geodesic_distance(m::Manifold, x1, x2)

Compute distance (along the shortest geodesic) between two given points
on the manifold `m` with ambient space representations `x1` and `x2`.
If integration is required then given tolerances
`PARAMS.quad_rel_tol` and `PARAMS.quad_abs_tol` are used.
"""
function geodesic_distance(x1::AbstractArray, x2::AbstractArray, m::Manifold)
    error("Function geodesic_distance is not yet defined for types $(typeof(x1)), $(typeof(x2)) and $(typeof(m)).")
end

"""
    inner_amb(x1, x2)

Compute the inner product of two points in a given manifold
in an ambient space selected for the specific representation.
"""
function inner_amb(x1::Point, x2::Point)
    return inner_amb(point2ambient(x1), point2ambient(x2), gettype(x1))
end

"""
    inner_amb(x1, x2, m)

Compute the inner product of two points in a given manifold `m`
in an ambient space selected for the specific representation.
"""
function inner_amb(x1::AbstractArray, x2::AbstractArray, m::Manifold)
    return sum(x1 .* x2)
end

"""
    ambient_distance(x1, x2)

Compute the distance in the ambient space between `x1` and `x2`.
"""
function ambient_distance(x1::Point, x2::Point)
    return ambient_distance(point2ambient(x1), point2ambient(x2), gettype(x1))
end

"""
    ambient_distance(x1, x2, m)

Compute the distance in the ambient space between points on `m` with
ambient space representations `x1` and `x2`.
"""
function ambient_distance(x1::AbstractArray, x2::AbstractArray, ::Manifold)
    return norm(x1 - x2)
end

"""
    riemannian_distortion(pts)

Compute the distortion due to curvature for the given set of points `pts`.
Implements FSDA Eq. (7.8)
"""
function riemannian_distortion(pts::Vector{<:Point})

    geod_dists = [i₁ == i₂ ? 0.0 : geodesic_distance(p₁, p₂) for (i₁, p₁) ∈ enumerate(pts), (i₂, p₂) ∈ enumerate(pts)]
    num = sum((geod_dists[i₁, i₂] - ambient_distance(pts[i₁], pts[i₂]))^2 for i₁ ∈ 1:length(pts), i₂ ∈ 1:length(pts))
    den = sum(geod_dists[i₁, i₂]^2 for i₁ ∈ 1:length(pts), i₂ ∈ 1:length(pts))
    if abs(num) < 1.e-16 #this should be enough to pevent NaNs
        return 0.0;
    else
        return num/den
    end
end
