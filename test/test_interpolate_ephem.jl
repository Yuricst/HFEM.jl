"""Test ephemeris interpolation"""

using LinearAlgebra
using SPICE

if !@isdefined(HFEM)
    include(joinpath(@__DIR__, "../src/HFEM.jl"))
end


test_interpolate_ephem = function ()
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
    parameters = HFEM.HFEMParameters(et0, DU, GMs, naif_ids, naif_frame, abcorr)

    # query states to be interpolated
    ets = range(et0, et0 + 30 * 86400.0, 1000)
    rvs = hcat([spkezr(naif_ids[2], et, naif_frame, abcorr, naif_ids[1])[1] for et in ets]...)

    ephem_interp = HFEM.InterpolatedEphemeris(
        naif_ids[2], 
        ets,
        rvs,
        false,
        parameters.TU
    )

    # evaluate position
    ets_test = range(et0 + 1e-4, et0 + 30 * 86400.0 - 1e-4, 200)

    state_interp = zeros(6, length(ets_test))
    state_spice = zeros(6, length(ets_test))
    diff = zeros(6, length(ets_test))
    for (idx,et_test) in enumerate(ets_test)
        state_interp[:,idx] = HFEM.get_state(ephem_interp, et_test)

        state_spice[:,idx], _ = spkezr(
            naif_ids[2],
            et_test,
            naif_frame,
            abcorr,
            naif_ids[1]
        )
        diff[:,idx] = state_spice[:,idx] - state_interp[:,idx]
    end
    # @show state_spice[:,end]
    # @show state_interp[:,end]
    # @show maximum(abs.(diff[1:3,:]))
    # @show maximum(abs.(diff[4:6,:]))
    @test maximum(abs.(diff[1:3,:])) < 1e-5     # this is in km
    @test maximum(abs.(diff[4:6,:])) < 1e-10    # this is in km/s
end

test_interpolate_ephem()