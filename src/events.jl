"""Event functions"""


"""
    mod_custom(a, n)

Custom modulo function
Ref: https://stackoverflow.com/questions/1878907/how-can-i-find-the-difference-between-two-angles
"""
function mod_custom(a, n)
    return a - floor(a / n) * n
end


"""
    angle_difference(ϕ_fwd::Real, ϕ_bck::Real)

Compute angle difference for periodic angles between 0 and 2π
"""
function angle_difference(ϕ_fwd::Real, ϕ_bck::Real)
    # modulo based
    dϕ = mod_custom((ϕ_bck - ϕ_fwd + π), 2π) - π
    return dϕ
end



"""
Compute osculating true anomaly from state and mu
    
# Arguments
- `state::Vector`: state vector in Cartesian coordinates, in order [x, y, z, vx, vy, vz]
- `mu::Float64`: gravitational parameter
- `to2pi::Bool`: if true, return the true anomaly in the range [0, 2π]
"""
function cart2trueanomaly(state::Vector, mu::Float64; to2pi::Bool=false)
    r = state[1:3]
    v = state[4:6]
    h = cross(r,v)
    hnorm = norm(h)
    vr = dot(v,r)/norm(r)
    ta = atan(hnorm*vr, hnorm^2/norm(r) - mu)
    if to2pi == true
        return mod(ta, 2π)
    else
        return ta
    end
end


function get_trueanomaly_event(
    θ_target::Real,
    t_bounds::Tuple{Real,Real},
    radius_bounds::Tuple{Real,Real},
    mu::Float64;
    θ_check_range::Float64 = deg2rad(60),
)
    if θ_target >= 0.9 * π
        to2pi = true
    else
        to2pi = false
    end

    function _condition(x, t, integrator)
        rnorm = norm(x[1:3])
        if (radius_bounds[1] <= rnorm <= radius_bounds[2]) && (t_bounds[1] <= t <= t_bounds[2])
            θ = cart2trueanomaly(x[1:6], mu; to2pi = to2pi)
            if abs(θ - θ_target) < θ_check_range
                #@printf("t = %2.4f, θ = %2.4f, θ_target = %2.4f\n", t, rad2deg(θ), rad2deg(θ_target))
                return angle_difference(θ, θ_target)
            else
                return NaN
            end
        else
            return NaN
        end
    end
    return _condition
end