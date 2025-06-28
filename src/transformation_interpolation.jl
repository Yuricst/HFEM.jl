"""Interpolate transformation matrix"""


function rotmat1(phi::Float64)
    return [1.0  0.0      0.0;
            0.0  cos(phi) sin(phi);
            0.0 -sin(phi) cos(phi)]
end


function rotmat2(phi::Float64)
    return [cos(phi) 0.0 -sin(phi);
            0.0      1.0  0.0;
            sin(phi) 0.0  cos(phi)]
end


function rotmat3(phi::Float64)
    return [ cos(phi) sin(phi) 0.0;
            -sin(phi) cos(phi) 0.0;
             0.0      0.0      1.0]
end


struct InterpolatedTransformation
    et_range::Tuple{Float64, Float64}
    frame_from::String
    frame_to::String
    axis_sequence::Tuple{Int, Int, Int}
    splines::Array{Spline1D, 1}
    rescale_epoch::Bool
    TU::Float64

    function InterpolatedTransformation(
        ets,
        frame_from::String,
        frame_to::String,
        rescale_epoch::Bool,
        TU::Float64;
        axis_sequence::Tuple{Int, Int, Int} = (3, 1, 3),
        spline_order::Int = 3,
    )
        @assert 1 <= spline_order <= 5
        if rescale_epoch
            @warn "rescale_epoch == true is buggy"
            times_input = (ets .- ets[1]) / TU
        else
            times_input = ets
        end
        euler_angles = zeros(3, length(ets))
        for (idx,et) in enumerate(ets)
            T = SPICE.pxform(frame_from, frame_to, et)
            euler_angles[:,idx] .= m2eul(T, axis_sequence...)
        end
        splines = [
            Spline1D(times_input, euler_angles[1,:]; k=spline_order, bc="error"),
            Spline1D(times_input, euler_angles[2,:]; k=spline_order, bc="error"),
            Spline1D(times_input, euler_angles[3,:]; k=spline_order, bc="error"),
        ]
        new((ets[1], ets[end]), frame_from, frame_to, axis_sequence, splines, rescale_epoch, TU)
    end
end


"""
Overload method for showing InterpolatedTransformation
"""
function Base.show(io::IO, transformation::InterpolatedTransformation)
    println("Interpolated transformation struct")
    @printf("    et0           : %s (%1.8f)\n", et2utc(transformation.et_range[1], "ISOC", 3), transformation.et_range[1])
    @printf("    etf           : %s (%1.8f)\n", et2utc(transformation.et_range[2], "ISOC", 3), transformation.et_range[2])
    @printf("    frame from    : %s\n", transformation.frame_from)
    @printf("    frame to      : %s\n", transformation.frame_to)
    @printf("    axis sequence : %s\n", transformation.axis_sequence)
end


function get_euler_angles(transformation::InterpolatedTransformation, et::Float64)
    if transformation.rescale_epoch
        et_eval = et * transformation.TU + transformation.et_range[1]
        @assert transformation.et_range[1] <= et <= transformation.et_range[2]
    else
        et_eval = et
        @assert transformation.et_range[1] <= et <= transformation.et_range[2]
    end
    euler_angles = [Dierckx.evaluate(transformation.splines[1], et_eval),
                    Dierckx.evaluate(transformation.splines[2], et_eval),
                    Dierckx.evaluate(transformation.splines[3], et_eval)]
    return euler_angles
end


"""Interpolate transformation matrix at a given epoch"""
function pxform(transformation::InterpolatedTransformation, et::Float64)
    euler_angles = get_euler_angles(transformation, et)
    T = rotmat3(euler_angles[1]) * rotmat1(euler_angles[2]) * rotmat3(euler_angles[3])
    return T
end