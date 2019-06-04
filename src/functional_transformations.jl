
function vel_derivative(f, t::Real)
    value = f(t)
    if isa(value, TupleArray)
        nums = Tuple(1:length(value.a))
        return TupleArray(map(i -> vel_derivative(s -> f(s)[i], t), nums))
    else
        return ForwardDiff.derivative(f, t)
    end
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
function velocity(c::AbstractCurvePt, ::Val{:continuous})
    m = values_in(gettype(c))
    return CurvePt(TangentBundleSpace(m)) do t::Real
        return TangentBundlePt(project_tv(vel_derivative(s -> point2ambient(c(s)), t), c(t)))
    end
end

function velocity(c::CurvePt, ::Val{:discretized}; dt = 1.e-7)
    return CurvePt(TangentBundleSpace(values_in(gettype(c)))) do t::Real
        return TangentBundlePt(project_tv((point2ambient(c(t+dt)) - point2ambient(c(t)))/dt, c(t)))
    end
end

@doc doc"
    curve_length(c, a = 0.0, b = 1.0, dtype = Val(:continuous))

Returns length of the curve `c` between parameters `a` and `b`. By default
`a = 0.0` and `b = 1.0`. Calculates velocity using method `dtype`.
"
function curve_length(c::AbstractCurvePt, a = 0.0, b = 1.0, dtype = Val(:continuous))
    v = velocity(c, dtype)
    #TODO: add quadrature tolerances to curve_length?
    #TODO: add velocity calculation options
    i_val, error = quadgk(t -> norm(v(t).x), a, b)
    return i_val
end

function gradient(f, x::Point)
    return gradient(f, point2ambient(x), gettype(x))
end

function gradient(f, x::AbstractArray, m::Manifold)
    return project_tv(ForwardDiff.gradient(f, x), x, m)
end

function gradient!(f, v::BNBArray, x::AbstractArray, m::Manifold)
    copyto!(v, ForwardDiff.gradient(fwrap, x))
    project_tv!(v, x, m)
    return v
end
