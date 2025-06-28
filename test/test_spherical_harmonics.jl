"""
Test spherical harmonics
"""

using LinearAlgebra
using OrdinaryDiffEq
using SPICE
using Test

if !@isdefined(HFEM)
    include(joinpath(@__DIR__, "../src/HFEM.jl"))
end

function test_spherical_harmonics()
    # load gggrd_20x20.tab file
    nmax = 8
    denormalize = true

    filepath = joinpath(@__DIR__, "../data/luna/gggrx_1200l_sha_20x20.tab")
    spherical_harmonics_data = HFEM.load_spherical_harmonics(filepath, nmax, denormalize)

    rvec = [0.0, 0.0, 1838.0]
    lmb, phi, r = HFEM.cart2sph(rvec)

    # evaluate a single acceleration term in the infinite series
    a_nm = HFEM.spherical_harmonics_nm_accel_PCPF(
        phi, lmb, r,
        spherical_harmonics_data["Cnm"],
        spherical_harmonics_data["Snm"],
        spherical_harmonics_data["GM"],
        spherical_harmonics_data["REFERENCE RADIUS"],
        2,0
    )

    # evaluate the cumulative acceleration due to spherical harmonics
    a_full = HFEM.spherical_harmonics_accel_PCPF(
        rvec,
        spherical_harmonics_data["Cnm"],
        spherical_harmonics_data["Snm"],
        spherical_harmonics_data["GM"],
        spherical_harmonics_data["REFERENCE RADIUS"],
        nmax
    )
    a_full_check = [
        3.0575156328702336e-7,
        -1.812482284272607e-8,
        4.304294506669716e-7
    ]
    @test all(isapprox.(a_full, a_full_check, atol=1e-10))
end


test_eom_NbodySH_SPICE = function()
    # load gggrd_20x20.tab file
    nmax = 8
    denormalize = true

    filepath_spherical_harmonics = joinpath(@__DIR__, "../data/luna/gggrx_1200l_sha_20x20.tab")
    #spherical_harmonics_data = HFEM.load_spherical_harmonics(filepath, nmax, denormalize)

    # define parameters
    GMs = [
        4.9028000661637961E+03,
        0.0, #3.9860043543609598E+05,
        0.0, #1.3271244004193938E+11,
    ]
    naif_ids = ["301", "399", "10"]
    naif_frame = "J2000"
    abcorr = "NONE"
    DU = 1737.4

    et0 = str2et("2020-01-01T00:00:00")
    parameters = HFEM.HFEMParameters(
        et0, DU, GMs, naif_ids, naif_frame, abcorr;
        filepath_spherical_harmonics = filepath_spherical_harmonics,
        nmax = nmax,
        frame_PCPF = "MOON_PA"
    )
    # @show parameters.DU, parameters.TU, parameters.VU
    # @show parameters.mus

    # initial state (in canonical scale)
    u0 = [1.0, 0.0, 0.3, 0.5, 1.0, 0.0]

    # time span (in canonical scale)
    tspan = (0.0, 3*86400/parameters.TU)

    # solve
    prob = ODEProblem(HFEM.eom_NbodySH_SPICE!, u0, tspan, parameters)
    sol = solve(prob, Vern7(), reltol=1e-14, abstol=1e-14)
    u_check = [-1.3008005902886173, 1.0821476922891953, -0.568881995188118, -0.13115294012566467, -0.6582795434237702, 0.0612542261511961]
    @test norm(sol.u[end] - u_check) < 1e-12

    # also solve the two-body problem for plotting
    # prob_twobody = ODEProblem(HFEM.eom_Nbody_SPICE!, u0, tspan, parameters)
    # sol_twobody = solve(prob_twobody, Vern7(), reltol=1e-14, abstol=1e-14)

    # # plot
    # fig = Figure(size=(600,600))
    # ax3d = Axis3(fig[1,1])
    # lines!(ax3d, Array(sol_twobody)[1,:], Array(sol_twobody)[2,:], Array(sol_twobody)[3,:], color=:blue)
    # lines!(ax3d, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:], color=:red, linewidth=0.5)
    # fig
end


test_spherical_harmonics()
test_eom_NbodySH_SPICE()