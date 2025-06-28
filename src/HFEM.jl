module HFEM

using Dierckx
using ForwardDiff
using LinearAlgebra
using Printf
using SPICE
import Symbolics

include("perturbations/third_body.jl")
include("perturbations/spherical_harmonics.jl")

include("jacobians_symbolic.jl")
include("ephemeris_interpolation.jl")
include("transformation_interpolation.jl")
include("parameters.jl")

include("eoms/eom_Nbody_SPICE.jl")
include("eoms/eom_Nbody_Interp.jl")
include("eoms/eom_NbodySH_SPICE.jl")

include("events.jl")

export InterpolatedEphemeris
export InterpolatedTransformation
export HFEMParameters

export eom_Nbody_SPICE!, eom_stm_Nbody_SPICE!, dfdx_Nbody_SPICE
export eom_Nbody_Interp!, eom_Nbody_Interp, eom_stm_Nbody_Interp!, dfdx_Nbody_Interp
export eom_stm_Nbody_Interp_fd!, dfdx_Nbody_Interp_fd

export get_trueanomaly_event

end # module HFEM
