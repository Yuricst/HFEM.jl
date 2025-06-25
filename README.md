# `HFEM.jl`: High-Fidelity Ephemeris Model for Astrodynamics

`HFEM.jl` is a minimal implementation of high-fidelity ephemeris model dynamics compatible with the [`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl) ecosystem (i.e. its solvers, parallelism, etc.).

What `HFEM.jl` contains:
-  full-ephemeris (and other - TBD) equations of motion relevant for astrodynamics
- commonly used event functions (TODO)

What `HFEM.jl` is *not*:
- it is not an integrator, i.e. there are no integration schemes (e.g. Runge-Kutta algorithms, step-correction, event detection features, etc.) impemented (at least for now)

We strive for minimal dependencies (listed in `Project.toml`), consisting of: `Dierckx`, `LinearAlgebra`, `OrdinaryDiffEq`, `SPICE`, `Symbolics`.


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

## Examples

### N-body Dynamics

```julia
using HFEM
using OrdinaryDiffEq

# define parameters
GMs = [
    4.9028000661637961E+03,
    3.9860043543609598E+05,
    1.3271244004193938E+11,                 # GM's in km^3/s^2
]
naif_ids = ["301", "399", "10"]
naif_frame = "J2000"
abcorr = "NONE"
DU = 3000.0                                 # distance canonical scale, in km
et0 = str2et("2020-01-01T00:00:00")         # reference epoch (i.e. epoch when t = 0 within the eom)
parameters = HFEM.HFEMParameters(et0, DU, GMs, naif_ids, naif_frame, abcorr)

# construct & solve ODE problem
x0 = [1.0, 0.0, 0.3, 0.5, 1.0, 0.0]         # initial state in DU & DU/TU
tspan = (0.0, 7*86400/parameters.TU)
prob = ODEProblem(HFEM.eom_Nbody_SPICE!, x0, tspan, parameters)
sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)

# propagate both state & STM
x0_stm = [x0; reshape(I(6),36)]
prob = ODEProblem(HFEM.eom_stm_Nbody_SPICE!, x0_stm, tspan, parameters)
sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)
```

To parallelize with multi-thread, use instead interpolated ephemerides:

```julia
et0 = str2et("2020-01-01T00:00:00")
etf = et0 + 30 * 86400.0
interpolate_ephem_span = [et0, etf]
parameters = HFEM.HFEMParameters(
    et0, DU, GMs, naif_ids, naif_frame, abcorr;
    interpolate_ephem_span = interpolate_ephem_span
)

# construct & solve ODE problem
x0 = [1.0, 0.0, 0.3, 0.5, 1.0, 0.0]         # initial state in DU & DU/TU
tspan = (0.0, 7*86400/parameters.TU)
prob = ODEProblem(HFEM.eom_Nbody_Interp!, x0, tspan, parameters)
sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)

# propagate both state & STM
x0_stm = [x0; reshape(I(6),36)]
prob = ODEProblem(HFEM.eom_stm_Nbody_Interp!, x0_stm, tspan, parameters)
sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)
```