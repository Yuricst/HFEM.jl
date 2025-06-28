# Core routines

## Parameters & Interpolations

```@autodocs
Modules = [HFEM]
Order   = [:function, :type, :struct]
Pages   = [
  "parameters.jl",
  "ephemeris_interpolation.jl",
  "transformation_interpolation.jl",
]
```

## Equations of motion

```@autodocs
Modules = [HFEM]
Order   = [:function, :type]
Pages   = [
  "eoms/eom_Nbody_Interp.jl",
  "eoms/eom_Nbody_SPICE.jl",
  "eoms/eom_NbodySH_Interp.jl",
  "eoms/eom_NbodySH_SPICE.jl"
]
```

## Perturbations

```@autodocs
Modules = [HFEM]
Order   = [:function, :type]
Pages   = [
  "perturbations/spherical_harmonics.jl",
  "perturbations/third_body.jl"
]
```

## Events

```@autodocs
Modules = [HFEM]
Order   = [:function, :type]
Pages   = [
  "events.jl",
]
```