
"""
    DCurves(M::Manifold, grid::AbstractVector)

Space of curves on manifold `M` discretized on the given `grid`.
"""
struct DCurves{𝔽,TM<:Manifold{𝔽},TG<:AbstractVector} <:
       Manifolds.AbstractPowerManifold{𝔽,TM,Manifolds.ArrayPowerRepresentation}
    manifold::TM
    grid::TG
end

"""
    UniformDCurves
"""
const UniformDCurves{TM} = DCurves{𝔽,TM,<:AbstractRange} where {𝔽,TM<:Manifold{𝔽}}

get_iterator(M::DCurves) = axes(M.grid, 1)

function manifold_dimension(M::DCurves)
    return manifold_dimension(M.manifold) * length(M.grid)
end

function representation_size(M::DCurves)
    return (representation_size(M.manifold)..., length(M.grid))
end

function reverse_srvf!(M::UniformDCurves, c_out, c_X, initial_point)
    dt = step(M.grid)
    tmp = Ref(allocate(initial_point))
    last_point = Ref(allocate(initial_point))
    copyto!(last_point[], initial_point)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        copyto!(_write(M, rep_size, c_out, i), last_point[])
        vector_transport_to!(M.manifold, tmp[], last_point[], c_X[i], last_point[])
        scale = dt * norm(tmp[])
        exp!(M.manifold, last_point[], last_point[], scale * tmp[])
    end
    return c_out
end

function srvf(M::UniformDCurves, c, backend::ProjectedDifferenceBackend)
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

function transport_srvf!(M::UniformDCurves, c_out, c_p, c_X, p)
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

function velocity_curve(M::UniformDCurves, c, backend::ProjectedDifferenceBackend)
    N = length(M.grid)
    factor = N-1
    c_out = allocate(c)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        wri = _write(M, rep_size, c_out, i)
        if i < N
            log!(M.manifold, wri, c[i], c[i+1])
        else
            log!(M.manifold, wri, c[N-1], c[N])
        end
        wri .*= factor
    end
    return c_out
end
