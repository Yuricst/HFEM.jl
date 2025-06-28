# `HFEM.jl`: High-Fidelity Ephemeris Model for Astrodynamics

![test workflow](https://github.com/Yuricst/HFEM.jl/actions/workflows/runtest.yml/badge.svg)

`HFEM.jl` is a minimal implementation of high-fidelity ephemeris model dynamics compatible with the [`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl) ecosystem (i.e. its solvers, parallelism, etc.).

What `HFEM.jl` contains:
- full-ephemeris equations of motion relevant for astrodynamics
- callback conditions for common astrodynamics events (e.g. detection of osculating true anomaly)
- ephemeris interpolation, to define equations of motion compatible with `EnsembleThreads` & automatic differentiation, e.g. `ForwardDiff`

What `HFEM.jl` is *not*:
- not an integrator, i.e. there are no integration schemes (e.g. Runge-Kutta algorithms, step-correction, event detection features, etc.) impemented (at least for now)

We strive for minimal dependencies (listed in `Project.toml`), consisting of: `Dierckx`, `ForwardDiff`, `LinearAlgebra`, `OrdinaryDiffEq`, `SPICE`, `Symbolics`.


## Quick start

1. `git clone` this repositiory
2. In your project directory, add:

```julia-repl
pkg> dev ./path/to/HFEM.jl
```

3. To run tests, `cd` to the root of this repository, then

```julia-repl
(@v1.10) pkg> activate .
(HFEM) pkg> test
```
