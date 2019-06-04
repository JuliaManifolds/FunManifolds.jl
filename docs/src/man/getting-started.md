# Getting started

This library provides high-level primitives for constructing algorithms dealing with functional and manifold-valued data. With growing interest in applying functional and geometrical data analysis methods, in particular in Machine Learning and Computer Vision, the problem of growing complexity of proposed algorithms arose. Often many parts of algorithms are similar (for example function optimization) but the state-of-the-art approach is to reinvent (or at least reimplement) the wheel for each particular case. This library solves this problem by providing high-level primitive components that facilitate faster prototyping of algorithms dealing with functional and geometric data. In particular the following tasks can be performed:

- data linearization,
- optimization of functions defined on manifolds,
- interpolation and approximation on manifolds.

Another goal for this package is reducing the barrier to entry for developers. Functional analysis and Riemannian geometry are advanced mathematical topics that are rarely presented to non-mathematicians. It is our belief that even a non-rigorous understanding of their basic concepts can be helpful. By having a better grasp of underlying structure of complex data one can use more appropriate processing methods, leading to improved results.

Content overview:

- Section [Introduction to geometry](@ref geometry-intro) introduces a few basic geometric concepts necessary for good understanding of mathematics underlying this library.
- Section [Examples](@ref examples) shows a variety of applications of this library.
- Section [Manifolds](@ref manifolds) describes different types of manifolds available in this library and how to use them.
