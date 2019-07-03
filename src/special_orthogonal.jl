
## Destructive matrix exponential using algorithm from Higham, 2008,
## "Functions of Matrices: Theory and Computation", SIAM
function exp(Aorig::StridedMatrix{T}) where T<:ForwardDiff.Dual
    A = copy(Aorig)
    n = size(A, 1)
    #ilo, ihi, scale = LAPACK.gebal!('B', A)     # modifies A
    nA = opnorm(A, 1)
    Inn = Matrix{T}(I, n, n)
    ## For sufficiently small nA, use lower order Padé-Approximations
    if (nA <= 2.1)
        if nA > 0.95
            C = T[17643225600.,8821612800.,2075673600.,302702400.,
                     30270240.,   2162160.,    110880.,     3960.,
                           90.,         1.]
    elseif nA > 0.25
        C = T[17297280.,8648640.,1995840.,277200.,
                 25200.,   1512.,     56.,     1.]
    elseif nA > 0.015
        C = T[30240.,15120.,3360.,
                420.,   30.,   1.]
    else
        C = T[120.,60.,12.,1.]
    end
    A2 = A * A
    P  = copy(Inn)
    U  = C[2] * P
    V  = C[1] * P
    for k in 1:(div(size(C, 1), 2) - 1)
        k2 = 2 * k
        P *= A2
        U += C[k2 + 2] * P
        V += C[k2 + 1] * P
    end
    U = A * U
    X = V + U
    #LAPACK.gesv!(V-U, X)
    #replaces gesv!
    X = (V-U)\X
    else
        s  = log2(nA/5.4)               # power of 2 later reversed by squaring
        if s > 0
            si = ceil(Int,s)
            A /= convert(T,2^si)
        end
        CC = T[64764752532480000.,32382376266240000.,7771770303897600.,
                1187353796428800.,  129060195264000.,  10559470521600.,
                    670442572800.,      33522128640.,      1323241920.,
                        40840800.,           960960.,           16380.,
                             182.,                1.]
        A2 = A * A
        A4 = A2 * A2
        A6 = A2 * A4
        U  = A * (A6 * (CC[14].*A6 .+ CC[12].*A4 .+ CC[10].*A2) .+
                  CC[8].*A6 .+ CC[6].*A4 .+ CC[4].*A2 .+ CC[2].*Inn)
        V  = A6 * (CC[13].*A6 .+ CC[11].*A4 .+ CC[9].*A2) .+
                   CC[7].*A6 .+ CC[5].*A4 .+ CC[3].*A2 .+ CC[1].*Inn

        X = V + U
        #LAPACK.gesv!(V-U, X)
        #replaces gesv!
        X = (V-U)\X

        if s > 0            # squaring to reverse dividing by power of 2
            for t=1:si; X *= X end
        end
    end
    return X
end


@doc doc"""
    SpecialOrthogonalSpace(n)

Space of rotations of an $n$-d space in a rotation matrix representation.
"""
struct SpecialOrthogonalSpace <: Manifold
    dim::Int32
    ambientSize::Int32
    SpecialOrthogonalSpace(n::Int64) = new(div(n*(n-1), 2), n)
end

function ambient_shape(m::SpecialOrthogonalSpace)
    return (m.ambientSize, m.ambientSize)
end

@doc doc"""
    SpecialOrthogonalPt(x)

Rotation of an $n$-d euclidean space represented by rotation matrix `x`
of size $n \times n$.
"""
struct SpecialOrthogonalPt{TArray <: AbstractArray{<:Real, 2}} <: Point
    x::TArray
    SpecialOrthogonalPt{TArray}(m::TArray) where TArray <: AbstractArray{<:Real, 2} = begin
        DEBUG && if !(isapprox(m*m', TArray(I, size(m, 1), size(m, 1)), atol=1e-4))
            error("Orthogonality conditions not maintained")
        end
        return new(m)
    end
end

SpecialOrthogonalPt(m::TArray) where TArray <: AbstractArray{<:Real,2} = SpecialOrthogonalPt{TArray}(m)

function deepcopy(x::SpecialOrthogonalPt)
    return SpecialOrthogonalPt(copy(x.x))
end

function gettype(x::SpecialOrthogonalPt)
    return SpecialOrthogonalSpace(size(x.x, 1))
end

function isapprox(x1::SpecialOrthogonalPt, x2::SpecialOrthogonalPt; atol = atoldefault(x1, x2), rtol = rtoldefault(x1, x2))
    return isapprox(x1.x, x2.x, atol = atol, rtol = rtol)
end

function id_element(s::SpecialOrthogonalSpace)
    return SpecialOrthogonalPt(Array{Float64}(I, s.ambientSize, s.ambientSize))
end

function inverse(pt::SpecialOrthogonalPt)
    return SpecialOrthogonalPt(inv(pt.x))
end

"""
    rotation2d(φ)

Creates a 2D rotation represented by `SpecialOrthogonalPt` corresponding
to given angle.
"""
function rotation2d(φ::Real)
    return SpecialOrthogonalPt([cos(φ) -sin(φ); sin(φ) cos(φ)])
end

"""
    rotation2d(φ)

Creates a 2D rotation represented by `SpecialOrthogonalPt` corresponding
to given angle. Uses static array representation.
"""
function rotation2d_s(φ::Real)
    return SpecialOrthogonalPt(SMatrix{2,2}([cos(φ) -sin(φ); sin(φ) cos(φ)]))
end

"""
    rotation3d_from_yaw_pitch_roll(yaw, pitch, roll)

Creates a 3D rotation represented by `SpecialOrthogonalPt` corresponding
to given yaw, pitch and roll.
"""
function rotation3d_from_yaw_pitch_roll(yaw::Real, pitch::Real, roll::Real)
    return SpecialOrthogonalPt(
        [cos(yaw) -sin(yaw) 0.; sin(yaw) cos(yaw) 0.;0. 0. 1.] *
        [cos(pitch) 0. sin(pitch); 0. 1. 0.; -sin(pitch) 0. cos(pitch)] *
        [1. 0. 0.; 0. cos(roll) -sin(roll); 0. sin(roll) cos(roll)])
end

"""
    rotation3d_from_yaw_pitch_roll_s(yaw, pitch, roll)

Creates a 3D rotation represented by `SpecialOrthogonalPt` corresponding
to given yaw, pitch and roll. Uses static array representation.
"""
function rotation3d_from_yaw_pitch_roll_s(yaw::Real, pitch::Real, roll::Real)
    return SpecialOrthogonalPt(SMatrix{3,3}(rotation3d_from_yaw_pitch_roll(yaw, pitch, roll).x))
end

@doc doc"""
    SpecialOrthogonalTV(x, v)

Tangent vector to an $n$-d rotation from tangent space at point
(rotation matrix) `x` represented in ambient space by matrix `v`.
"""
struct SpecialOrthogonalTV{TVectorPt <: AbstractArray{<:Real,2}, TVectorV <: Union{AbstractArray{<:Real,2}, MMatrix{N,N,<:Real} where N}} <: TangentVector
    at_pt::SpecialOrthogonalPt{TVectorPt}
    v::TVectorV
end

function deepcopy(v::SpecialOrthogonalTV)
    return SpecialOrthogonalTV(deepcopy(v.at_pt), copy(v.v))
end

function zero_tangent_vector(pt::SpecialOrthogonalPt)
    return SpecialOrthogonalTV(pt, _ensure_mutable(zeros(eltype(pt.x), size(pt.x))))
end

function zero_tangent_vector!(m::SpecialOrthogonalSpace, v::TV, p::AbstractArray) where TV<:BNBArray
    @condbc TV (v .= zeros(eltype(p), size(p)))
    return v
end

function +(v1::SpecialOrthogonalTV, v2::SpecialOrthogonalTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't add tangent vectors from different tangent spaces")
    end
    return SpecialOrthogonalTV(v1.at_pt, v1.v + v2.v)
end

function add_vec!(v1::SpecialOrthogonalTV, v2::SpecialOrthogonalTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't add tangent vectors from different tangent spaces $(at_point(v1)) and $(at_point(v2)).")
    end
    v1.v .+= v2.v
    return v1
end

@inline function add_vec!(v1::TV, v2::AbstractArray, at_pt::AbstractArray, m::SpecialOrthogonalSpace) where TV<:BNBArray
    @condbc TV (v1 .+= v2)
    return v1
end

function -(v1::SpecialOrthogonalTV, v2::SpecialOrthogonalTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't subtract tangent vectors from different tangent spaces")
    end
    return SpecialOrthogonalTV(v1.at_pt, v1.v - v2.v)
end

function sub_vec!(v1::SpecialOrthogonalTV, v2::SpecialOrthogonalTV)
    DEBUG && if !(at_point(v1) ≈ at_point(v2))
        error("Can't subtract tangent vectors from different tangent spaces $(at_point(v1)) and $(at_point(v2)).")
    end
    v1.v .-= v2.v
    return v1
end

@inline function sub_vec!(v1::TV, v2::AbstractArray, at_pt::AbstractArray, m::SpecialOrthogonalSpace) where TV<:BNBArray
    @condbc TV (v1 .-= v2)
    return v1
end

function *(α::T, v::SpecialOrthogonalTV) where T <: Real
    return SpecialOrthogonalTV(v.at_pt, α * v.v)
end

function mul_vec!(v::SpecialOrthogonalTV, α::Real)
    v.v .*= α
    return v
end

@inline function mul_vec!(v::TV, α::Real, at_pt::AbstractArray, m::SpecialOrthogonalSpace) where TV<:BNBArray
    @condbc TV (v .*= α)
    return v
end

function isapprox(v1::SpecialOrthogonalTV, v2::SpecialOrthogonalTV; atol = atoldefault(v1, v2), rtol = rtoldefault(v1, v2))
    if !(isapprox(v1.at_pt, v2.at_pt, atol = atol, rtol = rtol))
        return false
    end
    return isapprox(v1.v, v2.v, atol = atol, rtol = rtol)
end

function point2ambient(p::SpecialOrthogonalPt)
    return p.x
end

function ambient2point(m::SpecialOrthogonalSpace, amb::AbstractMatrix{<:Real})
    return SpecialOrthogonalPt(amb)
end

function project_point(m::SpecialOrthogonalSpace, amb::TVectorPt) where TVectorPt <: AbstractMatrix{<:Real}
    (q, r) = qr(amb)
    dim = size(amb, 1)
    dfix = _ensure_mutable(TVectorPt(Diagonal([1.0 for i ∈ 1:dim])))
    for i ∈ 1:size(amb, 1)
        #println(amb[:,i], q[:,i])
        if norm(amb[:,i] + q[:,i]) < norm(amb[:,i] - q[:,i])
            dfix[i, i] = -1.0
        end
    end
    q = q * dfix
    if det(q) < 0
        error("Input data $amb is broken: qr factorization $q, $r gives q matrix with det $(det(q)).")
    end
    return q
end

function project_point!(m::SpecialOrthogonalSpace, amb::TP) where TP<:BNBArray
    if TP<:NoBroadcastArray
        q = project_point(m, amb[])
        amb[] = q
    else
        q = project_point(m, amb)
        amb .= q
    end
    return amb
end

function ambient2tangent(v::AbstractMatrix{<:Real}, p::SpecialOrthogonalPt)
    return SpecialOrthogonalTV(p, v)
end

function project_tv!(m::SpecialOrthogonalSpace, v::TV, p::AbstractArray) where TV<:BNBArray
    if TV<:NoBroadcastArray
        v_id = v[] * p'
    else
        v_id = v * p'
    end
    @condbc TV (v .= (v_id .- v_id')./2 * p)
    return v
end

function tangent2ambient(v::SpecialOrthogonalTV)
    return v.v
end

"""
    ∘(x1::SpecialOrthogonalPt, x2::SpecialOrthogonalPt)

Composes two given rotations `x1` and `x2`.
"""
function ∘(x1::SpecialOrthogonalPt, x2::SpecialOrthogonalPt)
    return SpecialOrthogonalPt(x1.x * x2.x)
end

function exp!(p::TV, v::AbstractArray, at_pt::AbstractArray, m::SpecialOrthogonalSpace) where TV<:BNBArray
    @condbc TV (p .= at_pt * exp(v * at_pt'))
end

function inner(v1::AbstractMatrix, v2::AbstractMatrix, p::AbstractMatrix, m::SpecialOrthogonalSpace)
    return dot(v1, v2)/2
end

function log!(v::TV, x::A1, y::AbstractMatrix, m::SpecialOrthogonalSpace) where {TV<:BNBArray, A1<:AbstractMatrix}
    logarray = A1(real.(log(Array(x' * y))))
    @condbc TV (v .= _ensure_mutable(logarray * x))
    return v
end

function parallel_transport_geodesic!(vout::TV, vin::AbstractArray, at_pt::AbstractArray, to_point::AbstractArray, m::SpecialOrthogonalSpace) where TV<:BNBArray
    @condbc TV (vout .= vin * at_pt' * to_point)
    return vout
end

function geodesic(x1::SpecialOrthogonalPt, x2::SpecialOrthogonalPt)
    tv = log(x1, x2)
    return CurvePt(t -> exp(t*tv), gettype(x1))
end

function geodesic_at(t::Number, x1::AbstractMatrix, x2::AbstractMatrix, m::SpecialOrthogonalSpace)
    tv = log(x1, x2, m)
    return exp(t*tv, x1, m)
end

function distance(x1::AbstractArray, x2::AbstractArray, m::SpecialOrthogonalSpace)
    return norm(log(x1, x2, m))/sqrt(2.0)
end
