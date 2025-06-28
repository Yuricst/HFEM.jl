"""Interpolated ephemeris-based N-body equations of motion"""


"""
    eom_Nbody_Interp!(dx, x, params, t)

Right-hand side of N-body equations of motion compatible with `DifferentialEquations.jl`"""
function eom_Nbody_Interp!(dx, x, params, t)
    dx[1:3] = x[4:6]
    dx[4:6] = -params.mus[1] / norm(x[1:3])^3 * x[1:3]

    for i = 2:length(params.mus)
        pos_3body = get_pos(params.interpolated_ephems[i-1], params.et0 + t*params.TU) /params.DU
        dx[4:6] += third_body_accel(x[1:3], pos_3body, params.mus[i])
    end
    return nothing
end


"""
    eom_Nbody_Interp(x, params, t)
    
Right-hand side of N-body equations of motion compatible with `DifferentialEquations.jl`"""
function eom_Nbody_Interp(x, params, t)
    dx = [x[4:6]; -params.mus[1] / norm(x[1:3])^3 * x[1:3]]
    for i = 2:length(params.mus)
        pos_3body = get_pos(params.interpolated_ephems[i-1], params.et0 + t*params.TU) /params.DU
        dx[4:6] += third_body_accel(x[1:3], pos_3body, params.mus[i])
    end
    return dx
end


"""
    eom_stm_Nbody_Interp!(dx_stm, x_stm, params, t)
    
Right-hand side of N-body equations of motion with STM compatible with `DifferentialEquations.jl`"""
function eom_stm_Nbody_Interp!(dx_stm, x_stm, params, t)
    dx_stm[1:3] = x_stm[4:6]
    dx_stm[4:6] = -params.mus[1] / norm(x_stm[1:3])^3 * x_stm[1:3]

    for i = 2:length(params.mus)
        pos_3body = get_pos(params.interpolated_ephems[i-1], params.et0 + t*params.TU) /params.DU
        params.Rs[1+3(i-2):3(i-1)] = pos_3body
        dx_stm[4:6] += third_body_accel(x_stm[1:3], pos_3body, params.mus[i])
    end
    A = params.f_jacobian(x_stm[1:6], params.mus, params.Rs)
    dx_stm[7:42] = reshape((A * reshape(x_stm[7:42],6,6)')', 36)
    return nothing
end


"""
    dfdx_Nbody_Interp(x, u, params, t)
    
Evaluate Jacobian of N-body problem"""
function dfdx_Nbody_Interp(x, u, params, t)
    for i = 2:length(params.mus)
        pos_3body = get_pos(params.interpolated_ephems[i-1], params.et0 + t*params.TU) /params.DU
        params.Rs[1+3(i-2):3(i-1)] = pos_3body
    end
    return params.f_jacobian(x[1:6], params.mus, params.Rs)
end


"""
    dfdx_Nbody_Interp_fd(x, u, params, t)
    
Evaluate Jacobian of N-body problem"""
function dfdx_Nbody_Interp_fd(x, u, params, t)
    return ForwardDiff.jacobian(x -> eom_Nbody_Interp(x, params, t), x)
end


"""
    eom_stm_Nbody_Interp_fd!(dx_stm, x_stm, params, t)
    
Right-hand side of N-body equations of motion with STM compatible with `DifferentialEquations.jl`"""
function eom_stm_Nbody_Interp_fd!(dx_stm, x_stm, params, t)
    dx_stm[1:6] = eom_Nbody_Interp(x_stm[1:6], params, t)
    A = dfdx_Nbody_Interp_fd(x_stm[1:6], 0.0, params, t)
    A[1:3,4:6] .= I(3)   # force identity for linear map
    dx_stm[7:42] = reshape((A * reshape(x_stm[7:42],6,6)')', 36)
    return nothing
end