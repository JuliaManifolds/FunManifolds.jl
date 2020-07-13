
"""
    DCurves(M::Manifold, grid::AbstractVector)

Space of curves on manifold `M` discretized on the given `grid`.
"""
struct DCurves{𝔽,TM<:Manifold{𝔽},TG<:AbstractVector} <:
       Manifolds.AbstractPowerManifold{𝔽,Manifolds.ArrayPowerRepresentation}
    manifold::TM
    grid::TG
end

"""
    UniformDCurves
"""
const UniformDCurves{TM} = DCurves{TM,<:AbstractRange} where {TM<:Manifold}

get_iterator(M::DCurves) = axes(M.grid, 1)

function representation_size(M::DCurves)
    return (representation_size(M.manifold)..., length(M.grid))
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
