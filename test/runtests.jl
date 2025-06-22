"""Run tests"""

using Test

include(joinpath(@__DIR__, "../src/HFEM.jl"))

@testset "Jacobian" begin
    include("test_Nbody_SPICE.jl")
end