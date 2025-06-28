"""Run tests"""

using SPICE
using Test

include("utils.jl")
include(joinpath(@__DIR__, "../src/HFEM.jl"))

# furnish spice kernels
furnsh_kernels()

@testset "Ephemeris interpolation" begin
    include("test_interpolate_ephem.jl")
    include("test_interpolate_transformation.jl")
end

@testset "N-body ODE             " begin
    include("test_thirdbody.jl")
    include("test_Nbody_SPICE.jl")
    include("test_Nbody_Interp.jl")
    include("test_Nbody_ensemble.jl")
end

@testset "Spherical harmonics    " begin
    include("test_spherical_harmonics.jl")
    include("test_NbodySH_SPICE.jl")
    include("test_NbodySH_Interp.jl")
    include("test_NbodySH_ensemble.jl")
end

@testset "Callbacks              " begin
    include("test_callback.jl")
end