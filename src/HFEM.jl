module HFEM

using LinearAlgebra
using SPICE
import Symbolics

include("jacobians_symbolic.jl")
include("ephemeris_interpolation.jl")
include("parameters.jl")
include("eoms/eom_Nbody_SPICE.jl")
include("eoms/eom_Nbody_Interp.jl")

export HFEMParameters
export eom_Nbody_SPICE!, eom_stm_Nbody_SPICE!, dfdx_Nbody_SPICE 

end # module HFEM
