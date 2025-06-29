"""Utility functions"""


function vector_hessian_forwarddiff(f::Function, x)
    out = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(f, x), x)
    return reshape(out, (length(x), length(x), length(x)))
end


function eom_jacobian_fd(eom::Function, x, u, params, t)
    return ForwardDiff.jacobian(x -> eom(x, params, t), x)
end


function eom_hessian_fd(eom::Function, x, u, params, t)
    return vector_hessian_forwarddiff(x -> eom(x, params, t), x)
end