"""
    mean_karcher(pts)

Karcher mean (also known as the Frechet mean) of a given set of points.
"""
function mean_karcher(pts::Vector{<:Point})
    m = gettype(pts[1])
    ambs = map(pts) do pt
        return point2ambient(pt)
    end
    ext = point2ambient(mean_extrinsic(pts))
    x, cond = optimize(ext, m) do p
        return sum(geodesic_distance(p, pi, m)^2 for pi ∈ ambs)
    end
    return x
end

"""
    mean_extrinsic(pts)

Extrinsic mean is a mean calculated in the ambient space
and then the result is cast into the closest point on the given manifold.
"""
function mean_extrinsic(pts::Vector{<:Point})
    space = gettype(pts[1])
    ambs = [point2ambient(p) for p in pts]
    amb_vec = map(enumerate(ambs[1])) do (i, _)
        return mean(x[i] for x in ambs)
    end
    return project_point_wrapped(convert(typeof(ambs[1]), amb_vec), space)
end
