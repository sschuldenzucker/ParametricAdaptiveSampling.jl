[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://todo.com)

# ParametricAdaptiveSampling.jl

A little library for adaptive sampling of one-dimensional parametric curves, mainly for plotting.

A small package that exposes a function to adaptively sample points from a parametric function with one-dimensional real domain (e.g., `t::Float64 -> 1/t .* (cos(t), sin(t))`). The main use case is plotting but we do support an arbitrary number of dimensions and other use cases are certainly possible.

The main algorithm works by successive bipartition in `t` space until linear approximation is good enough at the midpoint (or we run out of points).

See [`sample_adaptive_parametric()`](https://todo.com) for details.

As a bonus, we also get non-parametric adaptive sampling (of `sin`, say). See [`sample_adaptive()`](https://todo.com).

## Installation

```
] add ParametricAdaptiveSampling
```

## Usage

```julia
using ParametricAdaptiveSampling
using Plots

f(t) = 1/t .* (cos(t), sin(t))  # for example

let _ts, points = sample_adaptive_parametric(f, 1, 10*pi)
  # `ts` is the set of t values that were used as samples. If you don't need them, just ignore them.
  plot(points)
end
```

`sample_adaptive_parametric()` has a number of keyword arguments to control its behavior. See [the docs](https://todo.com).

### Non-parametric version

```julia
using ParametricAdaptiveSampling
using Plots

# This does not return any t values because they are simply the x values.
plot(sample_adaptive(sin, 0, 2*pi))
```

Note that in a non-parametric use case, adaptive sampling is often not needed. See [the docs](https://todo.com)..

## Features

We do _not_ require (and actually do currently not make use of) differentiability of f. Therefore, this even works when `f` is some numerical method. This might actually be a quite relevant use case because adaptive sampling can reduce the number of data points required, which matters when `f` is slow (configuration matters here, though).

## Limitations

Some limitations and potential improvements are mentioned in the docs.

### Higher-dimensional parameter domain

One big limitation is that we only support one-dimensional parameter domains (t). Going multi-dimensional in t would be really cool but also make things _much_ more complicated because we now have to deal with meshes: for single-dimensional t, meshes are trivial because you just go from one data point to the next higher one. But in higher dimensions, it's suddenly not obvious anymore which point is connected to which. I _guess_ a recursive simplex mesh would be the proper higher-dimensional generalization of the method used here and it seems ok to implement. Requires some thought though because there are many ways to split a simplex (e.g., 4 obvious ones for 2-dimensional t) and there is some accounting to do for how they are connected. Create an issue if you wanna collaborate.

People gave some pointers in [this](https://discourse.julialang.org/t/library-hunt-adaptive-sampling-for-parametric-plots/128448) forum thread.

## Related Packages

- [AdaptiveSampling.jl](https://github.com/iuliancioarca/AdaptiveSampling.jl) does something related but seems to be more about compressing oversampled signals than sampling functions directly. The approach could potentially be adapted to achieve the same thing, though.
- [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl) is adaptive but computes a _global_ approximation of a function, which is not desired here.
- [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) is a library for differential equation simulations and [can](https://trixi-framework.github.io/TrixiDocumentation/stable/tutorials/adaptive_mesh_refinement/#adaptive_mesh_refinement) perform adaptive mesh refinement. This seems related in spirit to this package, but of course we do not care about differential equations here and (ab-)using Trixi for this would be overkill.
- [Mamba](https://github.com/brian-j-smith/Mamba.jl) was an Markov-Chain Monte Carlo framework that supposedly does some adaptive sampling under the hood.
- If we ever go down the route of higher-dimensional t, mesh libraries like [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) or [DistMesh](https://github.com/precise-simulation/distmesh-julia?tab=readme-ov-file) might be related, as well as Trixi's adaptive meshes.

## License

MIT

