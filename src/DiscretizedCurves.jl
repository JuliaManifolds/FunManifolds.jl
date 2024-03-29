
abstract type CurveInterpolationMethod end

"""
    ProjectionCurveInterpolation(embedding_interpolation_method::InterpolationType)

Interpolation in the embedding using embedding_interpolation_method and projection
onto the manifold.
"""
struct ProjectionCurveInterpolation{TIM<:Interpolations.InterpolationType} <:
       CurveInterpolationMethod
    embedding_interpolation_method::TIM
end

function ProjectionCurveInterpolation()
    im = Gridded(Linear())
    return ProjectionCurveInterpolation{typeof(im)}(im)
end

"""
    DiscretizedCurves(M::AbstractManifold, grid::AbstractVector)

Space of curves on manifold `M` discretized on the given `grid`.
"""
struct DiscretizedCurves{
    𝔽,
    TM<:AbstractManifold{𝔽},
    TG<:AbstractVector,
    TIM<:CurveInterpolationMethod,
} <: Manifolds.AbstractPowerManifold{𝔽,TM,Manifolds.ArrayPowerRepresentation}
    manifold::TM
    grid::TG
    interpolation_method::TIM
end

function DiscretizedCurves(M::AbstractManifold{𝔽}, grid::AbstractVector) where {𝔽}
    itpm = ProjectionCurveInterpolation()
    return DiscretizedCurves{𝔽,typeof(M),typeof(grid),typeof(itpm)}(M, grid, itpm)
end

struct ProjectionInterpolant{TM<:AbstractManifold,TEITP}
    M::TM
    embedding_itp::TEITP
end

function (pitp::ProjectionInterpolant)(t::Number)
    embedding_val = pitp.embedding_itp(t)
    return project(pitp.M, embedding_val)
end

function make_interpolant(
    M::DiscretizedCurves{<:Any,<:Any,<:Any,<:ProjectionCurveInterpolation},
    p,
)
    rep_size = representation_size(M.manifold)
    embedded_p = [embed(M.manifold, _read(M, rep_size, p, i)) for i in get_iterator(M)]
    embedding_itp = extrapolate(
        interpolate(
            (M.grid,),
            embedded_p,
            M.interpolation_method.embedding_interpolation_method,
        ),
        Flat(),
    )
    return ProjectionInterpolant(M.manifold, embedding_itp)
end

"""
    UniformDiscretizedCurves
"""
const UniformDiscretizedCurves{TM} =
    DiscretizedCurves{𝔽,TM,<:AbstractRange} where {𝔽,TM<:AbstractManifold{𝔽}}

embed!(M::DiscretizedCurves, q, p) = copyto!(q, p)
embed!(M::DiscretizedCurves, Y, p, X) = copyto!(Y, X)

get_iterator(M::DiscretizedCurves) = axes(M.grid, 1)

function manifold_dimension(M::DiscretizedCurves)
    return manifold_dimension(M.manifold) * length(M.grid)
end

function power_dimensions(M::DiscretizedCurves)
    return (length(M.grid),)
end

function representation_size(M::DiscretizedCurves)
    return (representation_size(M.manifold)..., length(M.grid))
end

function reverse_srvf!(M::UniformDiscretizedCurves, c_out, c_X, initial_point)
    dt = step(M.grid)
    tmp = Ref(allocate(initial_point))
    last_point = Ref(allocate(initial_point))
    copyto!(last_point[], initial_point)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        copyto!(_write(M, rep_size, c_out, i), last_point[])
        vector_transport_to!(
            M.manifold,
            tmp[],
            last_point[],
            _read(M, rep_size, c_X, i),
            last_point[],
        )
        scale = dt * norm(tmp[])
        exp!(M.manifold, last_point[], last_point[], scale * tmp[])
    end
    return c_out
end

function srvf(M::UniformDiscretizedCurves, c, backend::ProjectedDifferenceBackend)
    vel = velocity_curve(M, c, backend)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        copyto!(
            _write(M, rep_size, vel, i),
            q_function(M.manifold, _read(M, rep_size, c, i), _read(M, rep_size, vel, i)),
        )
    end
    return vel
end

function transport_srvf!(M::UniformDiscretizedCurves, c_out, c_p, c_X, p)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        vector_transport_to!(
            M.manifold,
            _write(M, rep_size, c_out, i),
            _read(M, rep_size, c_p, i),
            _read(M, rep_size, c_X, i),
            p,
        )
    end
    return c_out
end

function velocity_curve(M::UniformDiscretizedCurves, c, backend::ProjectedDifferenceBackend)
    N = length(M.grid)
    factor = N - 1
    c_out = allocate(c)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        wri = _write(M, rep_size, c_out, i)
        if i < N
            log!(M.manifold, wri, _read(M, rep_size, c, i), _read(M, rep_size, c, i + 1))
        else
            log!(M.manifold, wri, _read(M, rep_size, c, N - 1), _read(M, rep_size, c, N))
        end
        wri .*= factor
    end
    return c_out
end
