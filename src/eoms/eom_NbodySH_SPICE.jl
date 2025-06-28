"""Equations of motion for N-body problem with spherical harmonics"""


"""Right-hand side of N-body equations of motion compatible with `DifferentialEquations.jl`"""
function eom_NbodySH_SPICE!(dx, x, params, t)
    dx[1:3] = x[4:6]
    dx[4:6] = -params.mus[1] / norm(x[1:3])^3 * x[1:3]

    for i = 2:length(params.mus)
        pos_3body, _ = spkpos(
            params.naif_ids[i],
            params.et0 + t*params.TU,
            params.naif_frame,
            params.abcorr,
            params.naif_ids[1]
        )
        pos_3body /= params.DU
        dx[4:6] += third_body_accel(x[1:3], pos_3body, params.mus[i])
    end

    if !isnothing(params.spherical_harmonics_data)
        T_inr2pcpf = pxform(params.naif_frame, params.frame_PCPF, params.et0 + t*params.TU)
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


# """Right-hand side of N-body equations of motion with STMcompatible with `DifferentialEquations.jl`"""
# function eom_stm_Nbody_SPICE!(dx_stm, x_stm, params, t)
#     dx_stm[1:3] = x_stm[4:6]
#     dx_stm[4:6] = -params.mus[1] / norm(x_stm[1:3])^3 * x_stm[1:3]

#     for i = 2:length(params.mus)
#         pos_3body, _ = spkpos(
#             params.naif_ids[i],
#             params.et0 + t*params.TU,
#             params.naif_frame,
#             params.abcorr,
#             params.naif_ids[1]
#         )
#         pos_3body /= params.DU
#         params.Rs[1+3(i-2):3(i-1)] = pos_3body
#         dx_stm[4:6] += third_body_accel(x_stm[1:3], pos_3body, params.mus[i])
#     end
#     A = params.f_jacobian(x_stm[1:6], params.mus, params.Rs)
#     dx_stm[7:42] = reshape((A * reshape(x_stm[7:42],6,6)')', 36)
#     return nothing
# end


# """Evaluate Jacobian of N-body problem"""
# function dfdx_Nbody_SPICE(x, u, params, t)
#     for i = 2:length(params.mus)
#         pos_3body, _ = spkpos(
#             params.naif_ids[i],
#             params.et0 + t*params.TU,
#             params.naif_frame,
#             params.abcorr,
#             params.naif_ids[1]
#         )
#         pos_3body /= params.DU
#         params.Rs[1+3(i-2):3(i-1)] = pos_3body
#     end
#     return params.f_jacobian(x[1:6], params.mus, params.Rs)
# end