__precompile__()

"""
    FunManifolds

Main module for `FunManifolds.jl` -- a Julia package for functional differential geometry.
"""
module FunManifolds

using LinearAlgebra
using Interpolations
using Manifolds
using ManifoldsBase
using Markdown: @doc_str
using QuadGK
using StaticArrays

import Base: +, -, *, isapprox

import ManifoldsBase:
    allocate,
    allocate_result,
    allocate_result_type,
    allocation_promotion_function,
    array_value,
    base_manifold,
    check_point,
    check_vector,
    decorated_manifold,
    distance,
    embed,
    embed!,
    exp,
    exp!,
    geodesic,
    get_basis,
    get_coordinates,
    get_coordinates!,
    get_embedding,
    get_vector,
    get_vector!,
    get_vectors,
    injectivity_radius,
    inner,
    is_point,
    is_vector,
    inverse_retract,
    inverse_retract!,
    log,
    log!,
    manifold_dimension,
    mid_point,
    mid_point!,
    norm,
    number_eltype,
    number_of_coordinates,
    parallel_transport_direction,
    parallel_transport_direction!,
    parallel_transport_to,
    parallel_transport_to!,
    power_dimensions,
    project,
    project!,
    representation_size,
    retract,
    retract!,
    shortest_geodesic,
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_to,
    vector_transport_to!,
    zero_vector,
    zero_vector!

import Manifolds:
    apply,
    apply!,
    compose,
    compose!,
    get_iterator,
    identity_element,
    identity_element!,
    inv,
    inv!,
    optimal_alignment,
    optimal_alignment!,
    zero_vector

using Manifolds: _read, _write, AbstractRiemannianDiffBackend

mutable struct GeneralParams
    quad_rel_tol::Union{Real,Nothing}
    quad_abs_tol::Union{Real,Nothing}
end

const global PARAMS = GeneralParams(nothing, nothing)

function rtoldefault(M::AbstractManifold, x1, x2)
    return 1.e-8
end

function atoldefault(M::AbstractManifold, x1, x2)
    return 0.0
end

function concretize_tols(M::AbstractManifold, x1, x2; reltol=nothing, abstol=nothing)
    rtol = if reltol === nothing
        rtoldefault(M, x1, x2)
    else
        reltol
    end

    atol = if abstol === nothing
        atoldefault(M, x1, x2)
    else
        abstol
    end

    return (rtol, atol)
end

struct ProjectedDifferenceBackend{TDT<:Union{Number,Nothing}} <:
       AbstractRiemannianDiffBackend
    dt::TDT
end

ProjectedDifferenceBackend() = ProjectedDifferenceBackend{Float64}(1e-7)

include("interpolation.jl")

include("DiscretizedCurves.jl")
include("FunctionCurve.jl")
include("CurveWarping.jl")
include("CurveWarpingSRSF.jl")
include("curve_warping_alignment.jl")
include("functional_transformations.jl")

export CurveWarpingAction,
    CurveWarpingGroup,
    CurveWarpingSpace,
    CurveWarpingSRSFAction,
    CurveWarpingSRSFGroup,
    CurveWarpingSRSFSpace,
    DiscretizedCurves,
    FunctionCurveSpace,
    reverse_srsf,
    reverse_srvf,
    srsf,
    srvf,
    tsrvf,
    transport_srvf,
    transport_srvf!,
    UniformDiscretizedCurves

end #module
