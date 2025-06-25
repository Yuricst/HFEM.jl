"""
Test integrating N-body dynamics with SPICE call within eom.
Uses low-level API
"""

using LinearAlgebra
using OrdinaryDiffEq
using SPICE
using Printf
using Test

# if !@isdefined(HFEM)
    include(joinpath(@__DIR__, "../src/HFEM.jl"))
# end


# furnish spice kernels
spice_dir = ENV["SPICE"]
furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
furnsh(joinpath(spice_dir, "spk", "de440.bsp"))
furnsh(joinpath(spice_dir, "pck", "gm_de440.tpc"))


function print_matrix(A)
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            @printf("% 1.6e  ", A[i,j])
        end
        println()
    end
end


test_eom_Nbody_Interp = function()
    # define parameters 
    GMs = [
        4.9028000661637961E+03,
        3.9860043543609598E+05,
        1.3271244004193938E+11,
    ]
    naif_ids = ["301", "399", "10"]
    naif_frame = "J2000"
    abcorr = "NONE"
    DU = 3000.0

    et0 = str2et("2020-01-01T00:00:00")
    etf = et0 + 30 * 86400.0
    interpolate_ephem_span = [et0, etf]
    parameters = HFEM.HFEMParameters(
        et0, DU, GMs, naif_ids, naif_frame, abcorr;
        interpolate_ephem_span=interpolate_ephem_span)
    @show parameters.DU, parameters.TU, parameters.VU
    @show parameters.mus

    # initial state (in canonical scale)
    u0 = [1.0, 0.0, 0.3, 0.5, 1.0, 0.0]

    # time span (in canonical scale)
    tspan = (0.0, 7*86400/parameters.TU)

    # solve
    prob = ODEProblem(HFEM.eom_Nbody_Interp!, u0, tspan, parameters)
    sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)
    u_check = [0.5223150961449969, 2.096145052759142, -0.1636600693576513,
               -0.40936130128708903, 0.2538626713075729, -0.16005651999775147]
    # @show sol.u[end]
    # @show norm(sol.u[end] - u_check)
    @test norm(sol.u[end] - u_check) < 1e-11
end


test_eom_stm_Nbody_Interp = function(;verbose::Bool = false)
    # define parameters
    GMs = [
        4.9028000661637961E+03,
        3.9860043543609598E+05,
        1.3271244004193938E+11,
    ]
    naif_ids = ["301", "399", "10"]
    naif_frame = "J2000"
    abcorr = "NONE"
    DU = 3000.0

    et0 = str2et("2020-01-01T00:00:00")
    etf = et0 + 30 * 86400.0
    interpolate_ephem_span = [et0, etf]
    parameters = HFEM.HFEMParameters(et0, DU, GMs, naif_ids, naif_frame, abcorr;
        interpolate_ephem_span=interpolate_ephem_span)

    # initial state (in canonical scale)
    x0 = [1.0, 0.0, 0.3, 0.5, 1.0, 0.0]
    x0_stm = [x0; reshape(I(6),36)]

    # evaluate Jacobian
    jac_analytical = HFEM.dfdx_Nbody_SPICE(x0, x0, parameters, 0.0)
    @show jac_analytical

    f_eval = zeros(6)
    HFEM.eom_Nbody_Interp!(f_eval, x0, parameters, 0.0)
    jac_numerical = zeros(6,6)
    h = 1e-8
    for i = 1:6
        x0_copy = copy(x0)
        x0_copy[i] += h
        _f_eval = zeros(6)
        HFEM.eom_Nbody_Interp!(_f_eval, x0_copy, parameters, 0.0)
        jac_numerical[:,i] = (_f_eval - f_eval) / h
    end
    println("Analytical Jacobian:")
    print_matrix(jac_analytical)
    println()
    println("Numerical Jacobian:")
    print_matrix(jac_numerical)
    println()
    println("Diff:")
    print_matrix(jac_analytical - jac_numerical)
    @test maximum(abs.(jac_analytical - jac_numerical)) < 1e-6


    # time span (in canonical scale)
    tspan = (0.0, 1.0)

    # solve
    prob = ODEProblem(HFEM.eom_stm_Nbody_Interp!, x0_stm, tspan, parameters)
    sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)
    
    # construct STM
    STM_analytical = reshape(sol.u[end][7:42],6,6)'
    STM_numerical = zeros(6,6)
    h = 1e-8
    for i = 1:6
        x0_copy = copy(x0)
        x0_copy[i] += h
        sol_ptrb = solve(ODEProblem(HFEM.eom_Nbody_Interp!, x0_copy, tspan, parameters), Vern7(), reltol=1e-12, abstol=1e-12)
        STM_numerical[:,i] = (sol_ptrb.u[end][1:6] - sol.u[end][1:6]) / h
    end
    # println("Analytical STM:")
    # print_matrix(STM_analytical)
    # println()
    # println("Numerical STM:")
    # print_matrix(STM_numerical)
    # println()
    # println("Diff:")
    # print_matrix(STM_analytical - STM_numerical)
    @test maximum(abs.(STM_analytical - STM_numerical)) < 1e-6
end


#test_eom_Nbody_Interp()
test_eom_stm_Nbody_Interp(verbose = true)