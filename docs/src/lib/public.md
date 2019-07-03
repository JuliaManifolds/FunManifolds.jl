# Public Documentation

Documentation for `FunManifolds.jl`'s public interface.

## Contents

```@contents
Pages = ["public.md"]
Depth = 4
```

## Public Interface

Modules:
```@docs
FunManifolds.FunManifolds
```

### Basic Types

```@docs
Manifold
Point
TangentVector
```

### Manifolds

#### Euclidean Space

```@docs
EuclideanSpace
EuclideanPt
EuclideanTV
```
#### Sphere

```@docs
Sphere
SpherePt
SphereTV
```

#### Special Orthogonal Space

```@docs
SpecialOrthogonalSpace
SpecialOrthogonalPt
SpecialOrthogonalTV
∘(::SpecialOrthogonalPt, ::SpecialOrthogonalPt)
rotation2d
rotation2d_s
rotation3d_from_yaw_pitch_roll
rotation3d_from_yaw_pitch_roll_s
```

#### Tangent Space

```@docs
TSpaceManifold
TSpaceManifoldPt
TSpaceManifoldTV
```

#### Tangent Bundle

```@docs
TangentBundleSpace
TangentBundlePt
TangentBundleTV
```

#### Product Manifold

```@docs
ProductSpace
ProductPt
ProductTV
```

#### Power Manifold

```@docs
PowerSpace
PowerPt
PowerTV
```

#### Abstract manifolds of curves

```@docs
AbstractCurveSpace
AbstractCurvePt
values_in
paramgrid
```

#### Manifold of continuous curves on a manifold

```@docs
CurveSpace
CurvePt
CurveTV
```

### General operations on manifolds

```@docs
manifold_dimension
dim_ambient
ambient_shape
gettype
zero_tangent_vector
zero_tangent_vector!
add_vec
add_vec!
sub_vec
sub_vec!
mul_vec
mul_vec!
at_point
exp
exp!
retract
log
log!
inner
norm
geodesic
geodesic_at
distance
parallel_transport_geodesic
parallel_transport_geodesic!
ambient_distance
inner_amb
ambient2point
point2ambient
ambient2tangent
project_point
project_tv
project_tv!
tangent2ambient
riemannian_distortion
```

### Real-valued functions on manifold and optimisation

```@docs
optimize
```

### Mean-like functions

```@docs
mean_karcher
mean_extrinsic
```

## Index

```@index
Pages = ["public.md"]
```
