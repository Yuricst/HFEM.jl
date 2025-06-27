module HFEM

using Dierckx
using LinearAlgebra
using Printf
using SPICE
import Symbolics

include("perturbations/third_body.jl")
include("jacobians_symbolic.jl")
include("ephemeris_interpolation.jl")
include("parameters.jl")
include("eoms/eom_Nbody_SPICE.jl")
include("eoms/eom_Nbody_Interp.jl")
include("events.jl")

export InterpolatedEphemeris
export HFEMParameters
export eom_Nbody_SPICE!, eom_stm_Nbody_SPICE!, dfdx_Nbody_SPICE 
export get_trueanomaly_event

end # module HFEM
