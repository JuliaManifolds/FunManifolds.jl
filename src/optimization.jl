
"""
    optimize(f, x0, m; max_iter, step_eps, grad_eps)

Minimize objective function `f` on manifold `m`, starting at point `x0`.
Implements a variant of the steepest descent algorithm.
"""
function optimize(f, x0::AbstractArray{TF}, m::Manifold;
    max_iter = 1000, step_eps = 1.e-8, grad_eps = 1.e-8) where TF

    cur_x = _ensure_mutable(deepcopy(x0))
    desc_dir = -1. * gradient(f, cur_x, m)

    tmp_v = deepcopy(desc_dir)
    tmp_v2 = deepcopy(desc_dir)
    tmp_pt = deepcopy(cur_x)

    function s_retr!(x::AbstractArray, α)
        copyto!(tmp_v, desc_dir)
        mul_vec!(tmp_v, α, cur_x, m)
        retract!(x, tmp_v, cur_x, m)
        return x
    end

    function ϕ(α)
        s_retr!(tmp_pt, α)
        return f(tmp_pt)
    end
    function dϕ(α)
        s_retr!(tmp_pt, α)
        gradient!(f, tmp_v2, tmp_pt, m)
        return dot(tmp_v2, s)
    end
    function ϕdϕ(α)
        s_retr!(tmp_pt, α)
        gradient!(f, tmp_v2, tmp_pt, m)
        phi = f(tmp_pt)
        dphi = dot(tmp_v2, s)
        return (phi, dphi)
    end

    iter = 0
    cur_grad_norm = norm(desc_dir, cur_x, m)
    cur_step_size = 1.0
    lineseach = LineSearches.BackTracking{TF}()
    while iter < max_iter && cur_grad_norm > grad_eps && cur_step_size > step_eps
        fx = f(cur_x)
        desc_dir = -1. * gradient(f, cur_x, m)
        cur_grad_norm = norm(desc_dir, cur_x, m)
        iter += 1
        dϕ_0 = -cur_grad_norm^2

        #cur_x, step_size = line_search_adaptive(objective, desc_dir, cost, -grad_norm^2, lss)
        cur_step_size, fx = lineseach(ϕ, dϕ, ϕdϕ, 1.0, fx, dϕ_0)
        s_retr!(cur_x, cur_step_size)
    end

    return ambient2point(cur_x, m), :hopefully_ok
end
