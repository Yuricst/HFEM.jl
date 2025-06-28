# Overview

There are a number of equations of motion implemented in `HFEM.jl`.

- functions starting with `eom_` integrates the translational state (`[x,y,z,vx,vy,vz]`)
- functions starting with `eom_stm_` integrates both the translational state and the flattened 6-by-6 STM.

!!! note

    The STM is flattened row-wise, so to extract the state & STM from the `ODESolution`, make sure to reshape then tranapose; for example,

    ```julia
    x_stm_tf = sol.u[end]                       # concatenated state & STM (flattened)
    x_tf     = x_stm_tf[1:6]                    # final state [x,y,z,vx,vy,vz]
    STM_tf   = reshape(sol.u[end][7:42],6,6)'   # final 6-by-6 STM
    ```

## Dynamics model

In `HFEM.jl`, the dynamcis consists of the central gravitational term, together with the following perturbations:

- third-body perturbations
- spherical harmonics
- solar radiation pressure (todo)
- drag (todo)

```math
\dot{\boldsymbol{x}}(t) = 
\begin{bmatrix}
    \dot{\boldsymbol{r}}(t) \\ \dot{\boldsymbol{v}}(t)
\end{bmatrix} = 
\begin{bmatrix}
    \boldsymbol{v}(t)
    \\ -\dfrac{\mu}{\| \boldsymbol{r}(t) \|_2^3}\boldsymbol{r}(t)
\end{bmatrix}
+
\sum_{i} 
\begin{bmatrix}
    \boldsymbol{0}_{3 \times 1} 
    \\ \boldsymbol{a}_{\mathrm{3bd},i}(t)
\end{bmatrix}
+ 
\begin{bmatrix}
    \boldsymbol{0}_{3 \times 1}
    \\ \boldsymbol{a}_{\mathrm{SH},n_{\max}}(t)
\end{bmatrix}
+ 
\begin{bmatrix}
    \boldsymbol{0}_{3 \times 1}
    \\ \boldsymbol{a}_{\mathrm{SRP}}(t)
\end{bmatrix}
```

### Third-body perturbation

The third-body perturbation due to body $i$, $\boldsymbol{a}_{\mathrm{3bd},i}$, is given by

```math
\boldsymbol{a}_{\mathrm{3bd},i}(t)
= -\mu_i \left(
    \dfrac{\boldsymbol{r}(t) - \boldsymbol{r}_i(t)}{\| \boldsymbol{r}(t) - \boldsymbol{r}_i(t) \|_2^3}
    +
    \dfrac{\boldsymbol{r}_i(t)}{\| \boldsymbol{r}_i(t) \|_2^3}
\right)
```

where $\boldsymbol{r}_i$ is the position vector of the perturbing body.
In `HEFM.jl`, this term is implemented using Battin's $F(q)$ function:

```math
\boldsymbol{a}_{\mathrm{3bd},i}(t) =
-\dfrac{\mu_i}{\| \boldsymbol{r}(t) - \boldsymbol{r}_i(t) \|_2^3} (\boldsymbol{r}(t) + F(q_i)\boldsymbol{r}_i(t))
```

where $F(q_i)$ is given by

```math
F(q_i) = q_i \left( \dfrac{3 + 3q_i + q_i^2}{1 + (\sqrt{1 + q_i})^3} \right)
,\quad
q_i = \dfrac{\boldsymbol{r}(t)^T (\boldsymbol{r}(t) - 2\boldsymbol{r}_i(t))}{\boldsymbol{r}_i(t)^T \boldsymbol{r}_i(t)}
```

### Spherical Harmonics

The spherical harmonics perturbation $\boldsymbol{a}_{\mathrm{SH},n_{\max}}$ is given by

```math
\boldsymbol{a}_{\mathrm{SH},n_{\max}} = 
\sum_{n=2}^{n_{\max}} \sum_{m=0}^n \boldsymbol{a}_{\mathrm{SH},nm}
```

where $\boldsymbol{a}_{\mathrm{SH},nm}$ is given by

```math
\boldsymbol{a}_{\mathrm{SH},nm} = 
\begin{bmatrix}
    \ddot{x}_{n m} \\ \ddot{y}_{n m} \\ \ddot{z}_{n m}
\end{bmatrix}
```

where

```math
\begin{aligned}
& \ddot{x}_{n m} =
\begin{cases}
    \frac{G M}{R_{\oplus}^2} \cdot\left\{-C_{n 0} V_{n+1,1}\right\} & m = 0 \\[1.0em]
    \frac{G M}{R_{\oplus}^2} \cdot \frac{1}{2} \cdot\left\{\left(-C_{n m} V_{n+1, m+1}-S_{n m} W_{n+1, m+1}\right) + \frac{(n-m+2)!}{(n-m)!} \cdot\left(+C_{n m} V_{n+1, m-1}+S_{n m} W_{n+1, m-1}\right)\right\} & m > 0
\end{cases}
\\[2.5em]
& \ddot{y}_{n m} = 
\begin{cases}
    \frac{G M}{R_{\oplus}^2} \cdot\left\{-C_{n 0} W_{n+1,1}\right\} & m = 0 \\[1.0em]
    \frac{G M}{R_{\oplus}^2} \cdot \frac{1}{2} \cdot\left\{\left(-C_{n m} \cdot W_{n+1, m+1}+S_{n m} \cdot V_{n+1, m+1}\right) + \frac{(n-m+2)!}{(n-m)!} \cdot\left(-C_{n m} W_{n+1, m-1}+S_{n m} V_{n+1, m-1}\right)\right\} & m > 0
\end{cases}
\\[2.5em]
& \ddot{z}_{n m} = \frac{G M}{R_{\oplus}^2} \cdot\left\{(n-m+1) \cdot\left(-C_{n m} V_{n+1, m}-S_{n m} W_{n+1, m}\right)\right\}
\end{aligned}
```

and 

```math
V_{n m}=\left(\frac{R_{\oplus}}{r}\right)^{n+1} \cdot P_{n m}(\sin \phi) \cdot \cos m \lambda ,
\quad
W_{n m}=\left(\frac{R_{\oplus}}{r}\right)^{n+1} \cdot P_{n m}(\sin \phi) \cdot \sin m \lambda
```

(c.f. Montenbruck & Gill Chapter 3.2)


### Solar Radiation Pressure 

TODO


## List of equations of motion in `HFEM.jl`

The table below summarizes the equations of motion. Note: 

- `Nbody`: central gravity term + third-body perturbations ($\boldsymbol{a}_{\mathrm{3bd},i}$)
- `NbodySH`: central gravity term + third-body perturbations + spherical harmonics perturbations up to `nmax` degree ($\boldsymbol{a}_{\mathrm{SH},n_{\max}}$)
- The STM is integrated with the Jacobian, which is computed either analytically (using symbolic derivative) or via `ForwardDiff` (functions containing `_fd`)

| eom                   | eom + STM (analytical)  | eom + STM (ForwardDiff)      | `EnsembleThreads` compatibility |
|-----------------------|-------------------------|------------------------------|---------------------------------|
| `eom_Nbody_SPICE!`    | `eom_stm_Nbody_SPICE!`  |                              | no                              |
| `eom_Nbody_Interp!`   | `eom_stm_Nbody_Interp!` | `eom_stm_Nbody_Interp_fd!`   | yes                             |
| `eom_NbodySH_SPICE!`  |                         | `eom_stm_NbodySH_SPICE_fd!`  | no                              |
| `eom_NbodySH_Interp!` |                         | `eom_stm_NbodySH_Interp_fd!` | yes                             |


!!! note

    In order to use Julia's dual numbers, make sure to use a function that does not contain SPICE calls (i.e. use ones with `_interp` in the name); this is enabled by interpolating ahead of time ephemerides/transformation matrices.

!!! warning

    The accuracy of interpolated equations of motion (with `_interp` in the name) depends on the `interpolation_time_step`; if high-accuracy integration is required, it is advised to directly use the equations of motion that internally call SPICE (i.e. with `_SPICE` in the name)


## Initializing the parameter

We first need to define the parameter struct to be parsed as argument to the equations of motion.

Below is the most general example compatible with `eom_NbodySH_Interp!`/`eom_stm_NbodySH_Interp_fd!`:

```julia
using OrdinaryDiffEq
using HFEM

naif_ids = ["301", "399", "10"]        # NAIF IDs of bodies to be included; first ID is of the central body
GMs = [4.9028000661637961E+03, 3.9860043543609598E+05, 1.3271244004193938E+11]   # in km^3/s^2
naif_frame = "J2000"
abcorr = "NONE"
DU = 1e5                               # canonical distance unit, in km

nmax = 4                               # using up to 4-by-4 spherical harmonics
filepath_spherical_harmonics = "HFEM.jl/data/luna/gggrx_1200l_sha_20x20.tab"

et0 = str2et("2026-01-05T00:00:00")    # reference epoch
etf = et0 + 30 * 86400.0
interpolate_ephem_span = [et0, etf]    # range of epoch to interpolate ephemeris
interpolation_time_step = 1000.0       # time-step to sample ephemeris for interpolation

parameters = HFEM.HFEMParameters(
    et0, DU, GMs, naif_ids, naif_frame, abcorr;
    interpolate_ephem_span=interpolate_ephem_span,
    interpolation_time_step = interpolation_time_step,
    filepath_spherical_harmonics = filepath_spherical_harmonics,
    nmax = nmax,
    frame_PCPF = "MOON_PA",
)
```

Note:

- NAIF body IDs are defined according to: [https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html)
- if using `_SPICE` equations of motion, you do not need to parse `interpolate_ephem_span` and `interpolation_time_step`
- if using `Nbody` dynamics instead of `NbodySH`, you do not need to parse `filepath_spherical_harmonics`, `nmax`, and `frame_PCPF`


## Solving an Initial Value Problem

The integration is done with the `OrdinaryDiffEq.jl` library (or equivalently with `DifferentialEquations.jl`).

```julia
# initial state (in canonical scale)
x0 = [1.05, 0.0, 0.3, 0.5, 1.0, 0.0]

# time span (in canonical scale)
tspan = (0.0, 6 * 3600/parameters.TU)

# solve with SPICE
prob_spice = ODEProblem(HFEM.eom_NbodySH_SPICE!, x0, tspan, parameters)
sol_spice = solve(prob_spice, Vern8(), reltol=1e-14, abstol=1e-14)

# solve with interpolation
prob_interp = ODEProblem(HFEM.eom_NbodySH_Interp!, x0, tspan, parameters)
sol_interp = solve(prob_interp, Vern8(), reltol=1e-14, abstol=1e-14)
```
