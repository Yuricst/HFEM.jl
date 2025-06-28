# Choosing the equations of motion

There are a number of equations of motion implemented in `HFEM.jl`.

- functions starting with `eom_` integrates the translational state (`[x,y,z,vx,vy,vz]`)
- functions starting with `eom_stm_` integrates both the translational state and the flattened 6-by-6 STM.

The dynamics models are:

- `Nbody`: central gravity term + third-body perturbations
- `NbodySH`: central gravity term + third-body perturbations + spherical harmonics perturbations up to `nmax` degree

The STM is integrated with the Jacobian, which is computed either analytically (using symbolic derivative) or via `ForwardDiff` (functions containing `_fd`)

| eom                   | eom + STM (analytical)  | eom + STM (ForwardDiff)      | `EnsembleThreads` compatibility |
|-----------------------|-------------------------|------------------------------|---------------------------------|
| `eom_Nbody_SPICE!`    | `eom_stm_Nbody_SPICE!`  |                              | no                              |
| `eom_Nbody_Interp!`   | `eom_stm_Nbody_Interp!` | `eom_stm_Nbody_Interp_fd!`   | yes                             |
| `eom_NbodySH_SPICE!`  |                         | `eom_stm_NbodySH_SPICE_fd!`  | no                              |
| `eom_NbodySH_Interp!` |                         | `eom_stm_NbodySH_Interp_fd!` | yes                             |


> [!TIP]
> In order to use Julia's dual numbers, make sure to use a function that does not contain SPICE calls (i.e. use ones with `_interp` in the name); this is enabled by interpolating ahead of time ephemerides/transformation matrices.

> [!WARNING]  
> The accuracy of interpolated equations of motion (with `_interp` in the name) depends on the `interpolation_time_step`; if high-accuracy integration is required, it is advised to directly use the equations of motion that internally call SPICE (i.e. with `_SPICE` in the name)