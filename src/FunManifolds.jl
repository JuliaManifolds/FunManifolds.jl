__precompile__()

"""
Main module for `FunManifolds.jl` -- a Julia package for
differential geometry (also functional).
"""
module FunManifolds

using Interpolations
using Manifolds
using ManifoldsBase
using Markdown: @doc_str
using QuadGK

mutable struct GeneralParams
    quad_rel_tol::Union{Real,Nothing}
    quad_abs_tol::Union{Real,Nothing}
end

const global PARAMS = GeneralParams(nothing, nothing)

function rtoldefault(M::Manifold, x1, x2)
    return 1.e-8
end

function atoldefault(M::Manifold, x1, x2)
    return 0.0
end

function concretize_tols(M::Manifold, x1, x2; reltol=nothing, abstol=nothing)
    rtol = if reltol === nothing
        rtoldefault(M.M, x1, x2)
    else
        reltol
    end

    atol = if abstol === nothing
        atoldefault(M.M, x1, x2)
    else
        abstol
    end

    return (rtol, atol)
end


@doc doc"
    velocity(c, dtype)

Velocity curve for a given curve `c` calculated using differentiation
of type `dtype`.
If $c$ is a function such that $c\colon [0, 1] \to M$, then
`velocity(c)(t)` is a vector tangent to `c` at `t`.

Possible values of `dtype`:
* `Val(:continuous)` -- automatic differentiation resulting in a continuous
curve
* `Val(:discretized)` -- discretized derivative (TODO: use something for
configuration)
"
function velocity(M::Manifold, c, ::Val{:continuous})
    return nothing
end

function velocity(M::Manifold, c, ::Val{:discretized}; dt = 1.e-7)
    return nothing
end

@doc raw"
    curve_length(c, a = 0.0, b = 1.0, dtype = Val(:continuous))

Returns length of the curve `c` between parameters `a` and `b`. By default
`a = 0.0` and `b = 1.0`. Calculates velocity using method `dtype`.
"
function curve_length(M::Manifold, c, a = 0.0, b = 1.0, dtype = Val(:continuous))
    v = velocity(M, c, dtype)
    #TODO: add quadrature tolerances to curve_length?
    #TODO: add velocity calculation options
    i_val, error = quadgk(t -> norm(M, c(t), v(t)), a, b)
    return i_val
end

include("FunctionCurve.jl")

end #module
