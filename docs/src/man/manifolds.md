# [Manifolds](@id manifolds)

This library supports a number of manifolds.

```@contents
Pages = ["manifolds.md"]
Depth = 4
```

## Euclidean Space

Euclidean space is described by the following types (click to go to the documentation):
- [`EuclideanSpace`](@ref)
- [`EuclideanPt`](@ref)
- [`EuclideanTV`](@ref)

### Common sources of data

These types can be used to represent either points or free vectors (not attached to any point), for example translation vectors.

## Sphere

Unit spheres $\mathrm{S}^n$ embedded in $\mathbb{R}^{n+1}$ are described by the following types:
- [`Sphere`](@ref)
- [`SpherePt`](@ref)
- [`SphereTV`](@ref)

### Common sources of data

Spherical data often arises when:
- data is rescaled to have a given norm,
- geographical positions are analysed,
- representing orientation.

## Special Orthogonal Space

Special orthogonal spaces $\mathrm{SO}(n, \mathbb{R})$ represent rotation of $n$-dimension euclidean space $\mathbb{R}^n$. The selected embedding is described by rotation matrices of size $n\times n$. They are described by the following types:
- [`SpecialOrthogonalSpace`](@ref)
- [`SpecialOrthogonalPt`](@ref)
- [`SpecialOrthogonalTV`](@ref)

Additionally, composing two rotations from the same space represented by `SpecialOrthogonalPt` can be performed using a method of the `Base.∘` function:
- [`∘(::SpecialOrthogonalPt, ::SpecialOrthogonalPt)`](@ref)

### Common sources of data

- Rotations of euclidean spaces.

## Tangent Space

Representation of tangent spaces as manifolds. Mostly used internally. Since tangent spaces have a simple linear structure, they are not strictly necessary and only useful when one wants to explicitly mark something as a tangent space manifold.

- [`TSpaceManifold`](@ref)
- [`TSpaceManifoldPt`](@ref)
- [`TSpaceManifoldTV`](@ref)

## Tangent Bundle

Tangent bundle is the manifold of tangent vectors. In purely mathematical descriptions tangent vectors usually don't explicitly hold the point they are tangent at but it is useful to combine them as it was done in this library.

- [`TangentBundleSpace`](@ref)
- [`TangentBundlePt`](@ref)
- [`TangentBundleTV`](@ref)

### Common sources of data

Derivatives of curves and certain transformations of them.

## Product Manifold

Product manifolds combine data from two simpler manifolds.

- [`ProductSpace`](@ref)
- [`ProductPt`](@ref)
- [`ProductTV`](@ref)

### Common sources of data

Heterogeneous data (for example tuples of values of different origin).

## Power Manifold

Power manifolds are a special case of product manifolds where all data comes from the same manifold, which is a common use case. They are faster than deeply nested product manifolds.

- [`PowerSpace`](@ref)
- [`PowerPt`](@ref)
- [`PowerTV`](@ref)

### Common sources of data

- Sampling of a curve.
- Multiple sensors producing the same kind of data.

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
