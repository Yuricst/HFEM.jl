"""Parameters struct"""


mutable struct HFEMParameters
    et0::Float64
    DU::Real
    TU::Real
    VU::Real

    GMs::Vector{Float64}
    mus::Vector{Float64}
    naif_ids::Vector{String}
    naif_frame::String
    abcorr::String

    f_jacobian::Union{Nothing,Function}
    Rs::Vector{Float64}
end


function Base.show(io::IO, params::HFEMParameters)
    println("HFEMParameters struct")
    @printf("    et0        : %s (et = %1.8f)\n", et2utc(params.et0, "ISOC", 3), params.et0)
    @printf("    DU         : %1.8f\n", params.DU)
    @printf("    TU         : %1.8f\n", params.TU)
    @printf("    VU         : %1.8f\n", params.VU)
end


function HFEMParameters(
    et0::Float64,
    DU::Real,
    GMs::Vector{Float64},
    naif_ids::Vector{String},
    naif_frame::String,
    abcorr::String;
    get_jacobian_func::Bool = true,
)
    VU = sqrt(GMs[1]/DU)
    TU = DU/VU
    mus = GMs / GMs[1]         # scaled GM's

    # Jacobian function
    if get_jacobian_func
        f_jacobian = symbolic_Nbody_jacobian(length(GMs))
    else
        f_jacobian = nothing
    end
    Rs = zeros(3 * (length(mus)-1))  # storage for third-body positions

    return HFEMParameters(
        et0, DU, TU, VU,
        GMs, mus, naif_ids, naif_frame, abcorr,
        f_jacobian, Rs
    )
end