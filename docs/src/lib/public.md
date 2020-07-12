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
project_tangent
project_tangent!
tangent2ambient
riemannian_distortion
```

## Index

```@index
Pages = ["public.md"]
```
