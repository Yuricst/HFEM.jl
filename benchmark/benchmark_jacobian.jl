"""Benchmark the Jacobian of the N-body problem"""

using BenchmarkTools
using ForwardDiff
using LinearAlgebra
using OrdinaryDiffEq
using SPICE
using Test

include(joinpath(@__DIR__, "../src/HFEM.jl"))


benchmark_jacobian = function(;verbose::Bool = false)
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
    jac_analytical = HFEM.dfdx_Nbody_SPICE(x0, 0.0, parameters, 0.0)

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
    
    eom_x_only = function (x)
        return HFEM.eom_Nbody_Interp(x, parameters, 0.0)
    end
    jac_numerical_fd = ForwardDiff.jacobian(x -> HFEM.eom_Nbody_Interp(x, parameters, 0.0), x0)

    # jacobian via ForwardDiff
    println("Analytical Jacobian:")
    print_matrix(jac_analytical)
    println()
    println("Numerical Jacobian:")
    print_matrix(jac_numerical)
    println()
    println("ForwardDiff Jacobian:")
    print_matrix(jac_numerical_fd)
    println()

    # benchmark
    # println("Benchmarking analytical Jacobian:")
    # @benchmark HFEM.dfdx_Nbody_SPICE($x0, 0.0, $parameters, 0.0)
    
    println("Benchmarking ForwardDiff Jacobian:")
    @benchmark ForwardDiff.jacobian(x -> HFEM.eom_Nbody_Interp(x, $parameters, 0.0), $x0)
end


benchmark_jacobian()