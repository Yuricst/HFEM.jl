"""Run tests"""

using SPICE
using Test

include("utils.jl")
include(joinpath(@__DIR__, "../src/HFEM.jl"))

# furnish spice kernels
if !haskey(ENV, "SPICE")
    spice_dir = joinpath(@__DIR__, "../spice/test")
    furnsh(joinpath(spice_dir, "naif0012.tls"))
    furnsh(joinpath(spice_dir, "de440.bsp"))
    furnsh(joinpath(spice_dir, "gm_de440.tpc"))
else
    spice_dir = ENV["SPICE"]
    furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
    furnsh(joinpath(spice_dir, "spk", "de440.bsp"))
    furnsh(joinpath(spice_dir, "pck", "gm_de440.tpc"))
end


@testset "Ephemeris interpolation" begin
    include("test_interpolate_ephem.jl")
end

@testset "N-body ODE             " begin
    include("test_Nbody_SPICE.jl")
    include("test_Nbody_Interp.jl")
end