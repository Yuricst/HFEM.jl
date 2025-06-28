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
    furnsh(joinpath(spice_dir, "moon_pa_de440_200625.bpc"))
    furnsh(joinpath(spice_dir, "moon_de440_250416.tf"))
    furnsh(joinpath(spice_dir, "receding_horiz_3189_1burnApo_DiffCorr_15yr.bsp"))
else
    spice_dir = ENV["SPICE"]
    furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
    furnsh(joinpath(spice_dir, "spk", "de440.bsp"))
    furnsh(joinpath(spice_dir, "pck", "gm_de440.tpc"))
    furnsh(joinpath(spice_dir, "pck", "moon_pa_de440_200625.bpc"))
    furnsh(joinpath(spice_dir, "fk", "moon_de440_250416.tf"))
    furnsh(joinpath(spice_dir, "misc", "dsg_naif", "receding_horiz_3189_1burnApo_DiffCorr_15yr.bsp"))
end


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