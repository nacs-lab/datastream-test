#!/usr/bin/julia

using YAML

include("cl-utils.jl")
include("utils.jl")

function output_result(prefix, res, native)
    dev = get_device_id(res)
    suffix = dev
    if native
        suffix = suffix * "-native"
    end
    if endswith(prefix, "/")
        name = "$(prefix)$(suffix)"
    else
        name = "$(prefix)-$(suffix)"
    end
    open(name, "w") do io
        write(io, Float64(res["start"]))
        write(io, Float64(res["step"]))
        write(io, UInt32(res["nsteps"]))
        write_swapsign_leb128(io, res[native ? "diff_native" : "diff"])
    end
end

const prefix = ARGS[2]

for res in YAML.load_file(ARGS[1])
    output_result(prefix, res, true)
    output_result(prefix, res, false)
end
