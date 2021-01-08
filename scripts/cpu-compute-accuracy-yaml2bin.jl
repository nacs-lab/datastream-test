#!/usr/bin/julia

using YAML

function encode_leb128(io, vu)
    while true
        b = (vu & 0x7f) % UInt8
        vu = vu >> 7
        if vu == 0
            write(io, b)
            return
        end
        write(io, (b | 0x80)::UInt8)
    end
end

function encode_array(io, array)
    for v in array
        if v >= 0
            vu = UInt32(v) << 1
        else
            vu = UInt32(-v - 1) << 1 | 0x1
        end
        encode_leb128(io, vu)
    end
end

function output_result(prefix, key, res)
    suffix = key
    if endswith(prefix, "/")
        name = "$(prefix)$(suffix)"
    else
        name = "$(prefix)-$(suffix)"
    end
    open(name, "w") do io
        write(io, Float64(res["start"]))
        write(io, Float64(res["step"]))
        write(io, UInt32(res["nsteps"]))
        encode_array(io, res["diff"])
    end
end

const prefix = ARGS[2]

for (key, res) in YAML.load_file(ARGS[1])
    output_result(prefix, key, res)
end
