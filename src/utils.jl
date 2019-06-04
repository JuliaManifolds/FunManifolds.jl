
"""
    _ensure_mutable(a)

Converts given array `a` to the most similar mutable array.
Normal `Array`s and `MArray`s are not changed, `SArray`s are converted to
`MArray`s.
"""
_ensure_mutable(a::Array) = a
_ensure_mutable(a::MArray) = a
_ensure_mutable(a::SArray{S,T,N,L}) where {S<:Tuple,T,N,L} = MArray{S,T,N,L}(a)

"""
    ziptuples(a, b)

Zips tuples `a` and `b` in a fast, type-stable way. If they have different
lengths, the result is trimmed to the length of the shorter tuple.
"""
@generated function ziptuples(a::NTuple{N,Any}, b::NTuple{M,Any}) where {N,M}
    ex = Expr(:tuple)
    for i = 1:min(N, M)
        push!(ex.args, :((a[$i], b[$i])))
    end
    ex
end

"""
    ziptuples(a, b, c)

Zips tuples `a`, `b` and `c` in a fast, type-stable way. If they have different
lengths, the result is trimmed to the length of the shorter tuple.
"""
@generated function ziptuples(a::NTuple{N,Any}, b::NTuple{M,Any}, c::NTuple{L,Any}) where {N,M,L}
    ex = Expr(:tuple)
    for i = 1:min(N, M, L)
        push!(ex.args, :((a[$i], b[$i], c[$i])))
    end
    ex
end

"""
    enumeratetuple(a)

Type stable `enumerate` for tuples.
"""
@generated function enumeratetuple(a::Tuple)
    N = length(a.parameters)
    Expr(:tuple, [Expr(:tuple, k, Expr(:ref, :a, k)) for k in 1:N]...)
end

if get(ENV, "DEBUG", "0") != "0"
    macro check(a)
        return quote
            if !$(esc(a))
                throw(ErrorException("Debug check failed!"))
            end
        end
    end
else
    macro check(a) end
end

"""
    substitute_eltype(TCollection, TEl)

Returns collection type `TCollection` with element type changed to `TEl`.
"""
function substitute_eltype(TCollection::Type, TEl::Type)
    error("substitute_eltype is not defined for collection $TCollection and element type $TEl")
end


function substitute_eltype(::Type{Array{T, N}}, ::Type{Eltype}) where {T,N,Eltype}
    return Array{Eltype, N}
end

function substitute_eltype(::Type{SArray{S,T,A,B}}, ::Type{Eltype}) where {S,T,A,B,Eltype}
    return SArray{S, Eltype, A, B}
end

struct TupleArray{TT<:Tuple} <: AbstractVector{Any}
    a::TT
end

size(t::TupleArray) = (length(t.a),)

function getindex(t::TupleArray, i::Int)
    return t.a[i]
end

function copyto!(tout::TupleArray, tin::TupleArray)
    for i ∈ 1:size(tout)[1]
        copyto!(tout[i], tin[i])
    end
    return tout
end

function convert(::Type{T}, v::AbstractVector) where {N, T<:TupleArray{<:NTuple{N, Any}}}
    return T(Tuple(StaticArrays.unroll_tuple(v, StaticArrays.Length(N))))
end

import Base: BroadcastStyle, copy, similar
struct TupleArrayStyle{Length} <: Broadcast.AbstractArrayStyle{1} end

TupleArrayStyle{Length}(::Val{1}) where Length = TupleArrayStyle{Length}()

@generated function Base.BroadcastStyle(::Type{<:TupleArray{TT}}) where TT
    return TupleArrayStyle{length(TT.parameters)}()
end

@inline function similar(a::TupleArray)
    return TupleArray(map(i -> similar(i), a.a))
end

@inline function deepcopy(x::TupleArray)
    return TupleArray(map(i -> copy(i), x.a))
end

@inline function copy(bc::Base.Broadcast.Broadcasted{TupleArrayStyle{Length}}) where Length
    dim = axes(bc)
    length(dim) == 1 || throw(DimensionMismatch("TupleArray only supports one dimension"))
    length(dim[1]) == Length || throw(DimensionMismatch("TupleArray: how is this even possible?"))
    return TupleArray(ntuple(k -> @inbounds(Base.Broadcast._broadcast_getindex(bc, k)), Val(Length)))
end

function zero(::Type{TupleArray{T}}) where T<:Tuple
    return TupleArray(ntuple(k -> zero(T.types[k]), Val(length(T.types))))
end

function zero(::TupleArray{T}) where T<:Tuple
    return TupleArray(ntuple(k -> zero(T.types[k]), Val(length(T.types))))
end

(+)(t1::TupleArray, t2::TupleArray) = TupleArray(map(i -> i[1]+i[2], ziptuples(t1.a, t2.a)))
(-)(t1::TupleArray, t2::TupleArray) = TupleArray(map(i -> i[1]-i[2], ziptuples(t1.a, t2.a)))
#(*)(α::Number, t::TupleArray) = TupleArray(map(i -> α * i, t.a))

_ensure_mutable(a::TupleArray) = TupleArray(map(i -> _ensure_mutable(i), a.a))

const NoBroadcastArray = Union{AbstractArray{<:SArray},Ref{<:AbstractArray},AbstractArray{<:TupleArray},AbstractArray{<:Any,0}}
const BroadcastArray = AbstractArray
const BNBArray = Union{AbstractArray,Ref}

macro condbc(TV, e, additional_tuples = :(()))
    newhead = Symbol(string(e.head)[2:end])
    headsym = e.args[1]
    newbody = MacroTools.postwalk(t -> t === headsym ? Expr(:ref, headsym) : t, Expr(newhead, e.args...))
    return quote
        if $(esc(TV))<:TupleArray
            for i ∈ 1:length($(esc(headsym)))
                $(MacroTools.postwalk(t ->  if t === headsym || t ∈ additional_tuples.args
                                                Expr(:ref, esc(t), :i)
                                            elseif isa(t, Symbol)
                                                esc(t)
                                            else
                                                t
                                            end, e))
            end
        elseif $(esc(TV))<:NoBroadcastArray
            $(esc(newbody))
        else
            $(esc(e))
        end
    end
end

function similar_ambient(original::TupleArray, other::TupleArray)
    return TupleArray(map(i -> similar_ambient(original[i], other[i]), Tuple(1:length(original))))
end

function similar_ambient(original::AbstractArray{<:Number}, other::AbstractArray{<:Number})
    return similar(original, typeof(zero(eltype(original)) + zero(eltype(other))))
end

function similar_ambient(original::AbstractArray{<:AbstractArray}, other::AbstractArray{<:AbstractArray})
    return similar(original, typeof(similar_ambient(original[1], other[1])))
end
