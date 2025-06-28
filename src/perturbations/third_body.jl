"""Perturbations due to third-body"""


"""
Compute third-body acceleration via Battin's formula
"""
function third_body_accel(r_spacecraft, r_3body, mu_3body)
    s = r_spacecraft - r_3body
    q = dot(r_spacecraft, r_spacecraft - 2s)/dot(s, s)
    F = q * (3 + 3q + q^2)/(1 + sqrt(1+q)^3)
    return -mu_3body/norm(r_3body)^3 * (r_spacecraft + F*s)
end


function third_body_accel2(r_spacecraft, r_3body, mu_3body)
    d = r_spacecraft - r_3body
    q = dot(r_spacecraft, r_spacecraft - 2r_3body)/dot(r_3body, r_3body)
    F = q * (3 + 3q + q^2)/(1 + sqrt(1+q)^3)
    return -mu_3body/norm(d)^3 * (r_spacecraft + F*r_3body)
end


function third_body_accel_classical(r_spacecraft, r_3body, mu_3body)
    dr = r_spacecraft - r_3body
    return -mu_3body * (dr/norm(dr)^3 + r_3body/norm(r_3body)^3)
end