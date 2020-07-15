
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
const UniformDCurves{TM} = DCurves{TM,<:AbstractRange} where {TM<:Manifold}

get_iterator(M::DCurves) = axes(M.grid, 1)

function manifold_dimension(M::DCurves)
    return manifold_dimension(M.manifold) * length(M.grid)
end

function representation_size(M::DCurves)
    return (representation_size(M.manifold)..., length(M.grid))
end

function reverse_srvf!(M::DCurves, c_out, c_p, c_X, initial_point)
    dt = step(M.grid)
    tmp = @MVector [zero(TResPoint)]
    tmpview = view(tmp, 1)
    last_point = @MVector [initial_point]
    last_point_view = view(last_point, 1)
    rep_size = representation_size(M.manifold)
    for i in get_iterator(M)
        copyto!(_write(M, rep_size, c_out, i), last_point[1])
        vector_transport_to!(M.manifold, tmpview, c_p[i], c_X[i], last_point[1])
        scale = dt * norm(tmp[1])
        exp!(M.manifold, last_point_view, last_point[1], scale * tmp[1])
    end
    return c_out
end

function srvf(M::UniformDCurves, c, backend::Val{:forward_diff})
    vel = velocity_curve(M.manifold, c, backend)
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
