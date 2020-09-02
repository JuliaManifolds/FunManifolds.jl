"""
   KnownDerivativesMonotonicInterpolation(m)

Cubic interpolations with given derivatives at nodes.
"""
struct KnownDerivativesMonotonicInterpolation{TEl} <: Interpolations.MonotonicInterpolationType
    m::Vector{TEl}
end

function Interpolations.calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::KnownDerivativesMonotonicInterpolation) where {TCoeffs, TEl}

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k in 1:n-1
        Δ[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
        m[k] = method.m[k]
    end
    m[n] = method.m[n]
    return (m, Δ)
end

function ManifoldsBase.number_eltype(itp::Interpolations.AbstractInterpolation)
    return eltype(itp)
end
