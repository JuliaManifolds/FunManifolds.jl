# [Manifolds](@id manifolds)

This library supports a number of manifolds.

```@contents
Pages = ["manifolds.md"]
Depth = 4
```

## Manifold of continuous curves on a manifold

Curves (functions $f \colon [0,1] \to M$ where $M$ is a manifold) constitute a manifold themselves. This manifold is represented by the following types:

- [`CurveSpace`](@ref)
- [`CurvePt`](@ref)
- [`CurveTV`](@ref)

Their usage is restricted by the fact that they don't have an ambient space (functions [`point2ambient`](@ref), [`ambient2point`](@ref), [`tangent2ambient`](@ref) and [`ambient2tangent`](@ref) are not defined for this space). As a result, certain operations can't be performed:
- optimisation on the space of curves, that is finding optimal curves with respect to a certain criterion,
- differentiation of functions defined on the space of curves or taking values in it.
If you need one of these functionalities, use manifolds of discretized curves.

### Common sources of data

- Geodesics on other manifolds (see [`geodesic`](@ref)).
