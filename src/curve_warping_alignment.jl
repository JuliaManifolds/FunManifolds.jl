
abstract type WarpingIntegrator
end

struct GLWarpingIntegrator{TQNodes<:AbstractVector,TQWeights<:AbstractVector,GType<:AbstractRange,TC1,TC2} <: WarpingIntegrator
    quadnodes::TQNodes
    quadweights::TQWeights
    grid::GType
    c1::TC1
    c2::TC2
end

function GLWarpingIntegrator(c1, c2, node_num::Int, grid::AbstractRange)
    quadnodes, quadweights = gausslegendre( node_num )
    GLWarpingIntegrator(quadnodes, quadweights, grid, c1, c2)
end

function calc_val(wi::GLWarpingIntegrator, i, j, ki, kj)
    grid = wi.grid
    scale = (grid[j] - grid[j - kj])/(grid[i] - grid[i - ki])
    aa = (grid[i] - grid[i - ki])/2
    c1 = wi.c1
    c2 = wi.c2
    function integrand(t)
        st = grid[i-ki] + aa*(t + 1.0)
        return aa * norm(c1(st).x - sqrt(scale) * c2(scale * (st - grid[i - ki]) + grid[j - kj]).x)^2
    end
    return dot(wi.quadweights, integrand.(wi.quadnodes))
end

struct GLLinearIntegrator{T1<:AbstractVector,T2<:AbstractVector,GT<:AbstractRange,TTV} <: WarpingIntegrator
    c1_vals::T1
    c2_vals::T2
    grid::GT
    all_interval_splits::Matrix{Vector{Tuple{Float64,Char,Int32}}}
    tv::TTV
end

function GLLinearIntegrator(c1, c2, grid::AbstractRange, sigma)
    GLLinearIntegrator(map(c1, grid), map(c2, grid), grid, interval_splits(sigma),
        copy(c1(1.0)))
end

struct LinearDCurveIntegrator{T1,T2,TM1<:DCurves,TM2<:DCurves,TP,TTV} <: WarpingIntegrator
    c1_vals::T1
    c2_vals::T2
    M1::TM1
    M2::TM2
    p::TP
    all_interval_splits::Matrix{Vector{Tuple{Float64,Char,Int32}}}
    gamma_roots2::Matrix{Float64}
    tv::TTV
end

function generate_gamma_roots2(sigma)
    max_num_first = 1 + max_sigma_len(sigma, 1)
    max_num_second = 1 + max_sigma_len(sigma, 2)
    ret = Matrix{Float64}(undef, max_num_first, max_num_second)
    for ki ∈ 0:max_num_first-1
        for kj ∈ 0:max_num_second-1
            ret[ki+1, kj+1] = ki != 0 ? sqrt(kj / ki) : 1.0
        end
    end
    return ret
end

function make_linear_dcurve_integrator(M1::DCurves, M2::DCurves, dc1, dc2, p, sigma)

    return LinearDCurveIntegrator(
        dc1,
        dc2,
        M1,
        M2,
        p,
        interval_splits(sigma),
        generate_gamma_roots2(sigma),
        allocate(_read(M1, representation_size(M1.manifold), dc1, 1)),
    )
end

function max_sigma_len(sigma, elem_ind::Int)
    ret = 0
    for s ∈ sigma
        ret = max(ret, s[elem_ind])
    end
    return ret
end

function generate_interval_splits(ki, kj)
    ret = Vector{Tuple{Float64, Char, Int32}}()
    for i ∈ 1:ki
        push!(ret, ( i/ki, 'i', i ))
    end
    for j ∈ 1:(kj-1)
        push!(ret, ( j/kj, 'j', j ))
    end
    sort!(ret, by = x -> x[1])
    return ret
end

function interval_splits(sigma)
    max_num_first = 1 + max_sigma_len(sigma, 1)
    max_num_second = 1 + max_sigma_len(sigma, 2)

    all_interval_splits = Matrix{Vector{Tuple{Float64, Char, Int32}}}(undef, max_num_first, max_num_second)
    for ki ∈ 0:max_num_first-1
        for kj ∈ 0:max_num_second-1
            all_interval_splits[ki+1, kj+1] = generate_interval_splits(ki, kj)
        end
    end
    return all_interval_splits
end

function calc_val(wi::GLLinearIntegrator, i, j, ki, kj)
    false && if (ki == 0 && kj > 1) || (kj == 0 && ki > 1)
        throw("when one of the paths has no edges, the other must be of length 1")
    end

    rep_size = representation_size(wi.M1)
    integral = zero(eltype(wi.tv))
    if ki == 0
        integral = norm(wi.M2.manifold, wi.p, _read(wi.M2, rep_size, wi.c2_vals, j))^2
    elseif kj == 0
        integral = norm(wi.M1.manifold, wi.p, _read(wi.M1, rep_size, wi.c1_vals, i))^2
    else
        gamma_root1 = ki != 0 ? one(eltype(wi.tv)) : zero(eltype(wi.tv))
        gamma_root2 = ki != 0 ? convert(eltype(wi.tv), sqrt(kj / ki)) : one(eltype(wi.tv))
        dt = wi.grid[i] - wi.grid[i-ki]

        interval_splits = wi.all_interval_splits[ki+1, kj+1]
        last_t = .0
        cur_i = Int32(0)
        cur_j = Int32(0)
        for (new_t, ij, new_val) ∈ interval_splits
            # integrate
            fastwarp!(
                wi.tv,
                gamma_root1,
                _read(wi.M1, rep_size, wi.c1_vals, i-ki+cur_i),
                gamma_root2,
                _read(wi.M2, rep_size, wi.c2_vals, j-kj+cur_j),
            )
            norm_sq = norm(wi.M1, wi.p, wi.tv)^2
            integral += (new_t - last_t) * norm_sq
            # update vars
            last_t = new_t
            if ij == 'i'
                cur_i = new_val
            elseif ij == 'j'
                cur_j = new_val
            else
                throw("incorrect value of ij")
            end
        end
        integral *= dt
    end
    return integral
end

function calc_val(wi::LinearDCurveIntegrator, i, j, ki, kj)
    false && if (ki == 0 && kj > 1) || (kj == 0 && ki > 1)
        throw("when one of the paths has no edges, the other must be of length 1")
    end

    integral = .0
    if ki == 0
        integral = dot(wi.c2_vals[j], wi.c2_vals[j])
    elseif kj == 0
        integral = dot(wi.c1_vals[i], wi.c1_vals[i])
    else
        gamma_root1 = ki != 0 ? one(eltype(wi.tv)) : zero(eltype(wi.tv))
        gamma_root2 = wi.gamma_roots2[ki+1, kj+1]
        dt = wi.M1.grid[i] - wi.M1.grid[i-ki]

        interval_splits = wi.all_interval_splits[ki+1, kj+1]
        last_t = zero(eltype(wi.tv))
        cur_i = Int32(0)
        cur_j = Int32(0)
        for (new_t, ij, new_val) ∈ interval_splits
            # integrate
            wi.tv .= gamma_root1.*wi.c1_vals[i-ki+cur_i] .- gamma_root2 .* wi.c2_vals[j-kj+cur_j]
            norm_sq = norm(wi.M1.manifold, wi.p, wi.tv)^2
            integral += (new_t - last_t) * norm_sq
            # update vars
            last_t = new_t
            if ij == 'i'
                cur_i = new_val
            elseif ij == 'j'
                cur_j = new_val
            else
                throw("incorrect value of ij")
            end
        end
        integral *= dt
    end
    return integral
end

mutable struct PairwiseOptimalWarpingConfiguration
    sigma::Vector{Tuple{Int64,Int64}}
end

const POWC = PairwiseOptimalWarpingConfiguration(
    [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (4, 1), (3, 4), (4, 3)])


function _pairwise_optimal_warping(integrator::WarpingIntegrator, grid1, grid2, sigma)
    cost_matrix = Float64[Inf for i ∈ grid1, j ∈ grid2]
    prev_matrix = Tuple{Int64, Int64}[(0, 0) for i ∈ grid1, j ∈ grid2]
    cost_matrix[1, 1] = .0
    for i ∈ 1:length(grid1)
        for j ∈ 1:length(grid2)
            for (ki, kj) ∈ sigma
                if i-ki >= 1 && j-kj >= 1
                    I = calc_val(integrator, i, j, ki, kj)
                    cur_cost = cost_matrix[i-ki, j-kj] + I
                    if cur_cost < cost_matrix[i, j]
                        cost_matrix[i, j] = cur_cost
                        prev_matrix[i, j] = (i-ki, j-kj)
                    end
                end
            end
        end
    end
    xs = Float64[]
    ys = Float64[]
    curi, curj = length(grid1), length(grid2)
    while prev_matrix[curi, curj] != (0, 0)
        push!(xs, grid1[curi])
        push!(ys, grid2[curj])
        curi, curj = prev_matrix[curi, curj]
    end
    push!(xs, 0.0)
    push!(ys, 0.0)

    reverse!(xs)
    reverse!(ys)
    # println("pairwise | ", xs, ys)
    warp_interp = interpolate(xs, ys, FritschButlandMonotonicInterpolation())
    warp_curve = extrapolate(interpolate(Vector(grid1), [warp_interp(t) for t ∈ grid1], FritschButlandMonotonicInterpolation()), Flat())
    final_cost = cost_matrix[end, end]
    return (warp_curve, final_cost, xs, ys)
end

"""
    pairwise_optimal_warping(M1::DCurves, M2::DCurves, c1, c2, p, sigma = nothing)

Return a warping `γ` and elastic distance such that when `c2` is warped by `γ`
it is the closest to `c1`. The warping is determined on the discretization
grids. Default `sigma` is equal to ```[(1, 1), (1, 2), (2, 1), (1, 3), (3, 1),
(2, 3), (3, 2), (1, 4), (4, 1), (3, 4), (4, 3)]```

`p` is the point to which vectors from `c1` and `c2` are tangent.

Warning: for functions with vastly different amplitudes the returned
reparametrization may not be useful. If you know how to make this statement
more precise, please let us know. One example of this effect is shown
in function `testSRVFAcc` in `FunManifoldsExamples`.
"""
function pairwise_optimal_warping(M1::DCurves, M2::DCurves, c1, c2, p,
    sigma_arg = nothing::Union{Nothing, Vector{Tuple{Int64,Int64}}})

    if isa(sigma_arg, Nothing)
        sigma = POWC.sigma
    else
        sigma = sigma_arg
    end

    integrator = make_linear_dcurve_integrator(M1, M2, c1, c2, p, sigma)
    grid1 = M1.grid
    grid2 = M2.grid

    return _pairwise_optimal_warping(integrator, grid1, grid2, sigma)
end
