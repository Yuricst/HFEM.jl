"""Equations of motion for N-body problem with spherical harmonics using interpolated structs"""


"""Right-hand side of N-body equations of motion compatible with `DifferentialEquations.jl`"""
function eom_NbodySH_Interp!(dx, x, params, t)
    dx[1:3] = x[4:6]
    dx[4:6] = -params.mus[1] / norm(x[1:3])^3 * x[1:3]

    for i = 2:length(params.mus)
        pos_3body = HFEM.get_pos(params.interpolated_ephems[i-1], params.et0 + t*params.TU) /params.DU
        dx[4:6] += third_body_accel(x[1:3], pos_3body, params.mus[i])
    end

    if !isnothing(params.spherical_harmonics_data)
        T_inr2pcpf = pxform(params.interpolated_transformation, params.et0 + t*params.TU)
        a_SH = spherical_harmonics_accel(
            T_inr2pcpf,
            x[1:3] * params.DU,
            params.spherical_harmonics_data["Cnm"],
            params.spherical_harmonics_data["Snm"],
            params.spherical_harmonics_data["GM"],
            params.spherical_harmonics_data["REFERENCE RADIUS"],
            params.spherical_harmonics_data["nmax"]
        )
        dx[4:6] += a_SH / (params.VU/params.TU)
    end

    return nothing
end


"""Right-hand side of N-body equations of motion compatible with `DifferentialEquations.jl`"""
function eom_NbodySH_Interp(x, params, t)
    dx = [x[4:6]; -params.mus[1] / norm(x[1:3])^3 * x[1:3]]

    for i = 2:length(params.mus)
        pos_3body = HFEM.get_pos(params.interpolated_ephems[i-1], params.et0 + t*params.TU) /params.DU
        dx[4:6] += third_body_accel(x[1:3], pos_3body, params.mus[i])
    end

    if !isnothing(params.spherical_harmonics_data)
        T_inr2pcpf = pxform(params.interpolated_transformation, params.et0 + t*params.TU)
        a_SH = spherical_harmonics_accel(
            T_inr2pcpf,
            x[1:3] * params.DU,
            params.spherical_harmonics_data["Cnm"],
            params.spherical_harmonics_data["Snm"],
            params.spherical_harmonics_data["GM"],
            params.spherical_harmonics_data["REFERENCE RADIUS"],
            params.spherical_harmonics_data["nmax"]
        )
        dx[4:6] += a_SH / (params.VU/params.TU)
    end
    return dx
end


"""Evaluate Jacobian of N-body problem"""
function dfdx_NbodySH_Interp_fd(x, u, params, t)
    return ForwardDiff.jacobian(x -> HFEM.eom_NbodySH_Interp(x, params, t), x)
end



"""Right-hand side of N-body equations of motion with STM compatible with `DifferentialEquations.jl`"""
function eom_stm_NbodySH_Interp_fd!(dx_stm, x_stm, params, t)
    dx_stm[1:6] = eom_NbodySH_Interp(x_stm[1:6], params, t)
    A = dfdx_NbodySH_Interp_fd(x_stm[1:6], 0.0, params, t)
    dx_stm[7:42] = reshape((A * reshape(x_stm[7:42],6,6)')', 36)
    return nothing
end