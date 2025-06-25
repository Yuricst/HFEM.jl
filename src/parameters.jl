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
    interpolated_ephems::Union{Nothing,Vector{InterpolatedEphemeris}}

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
    interpolate_ephem_span::Union{Nothing,Vector{Float64}} = nothing,
    interpolation_time_step::Real = 3600.0,
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

    if isnothing(interpolate_ephem_span)
        interpolated_ephems = nothing
    else
        N_interp = Int(ceil((interpolate_ephem_span[2] - interpolate_ephem_span[1]) / interpolation_time_step))
        ets_interp = range(interpolate_ephem_span[1], interpolate_ephem_span[2], N_interp)
        interpolated_ephems = []
        for ID in naif_ids[2:end]
            rvs_interp = hcat([spkezr(ID, et, naif_frame, abcorr, naif_ids[1])[1] for et in ets_interp]...)
            push!(interpolated_ephems, HFEM.InterpolatedEphemeris(ID, ets_interp, rvs_interp, false, TU))
        end
    end

    return HFEMParameters(
        et0, DU, TU, VU,
        GMs, mus, naif_ids, naif_frame, abcorr,
        interpolated_ephems,
        f_jacobian, Rs
    )
end