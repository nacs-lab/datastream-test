#!/usr/bin/julia

using YAML

include("cl-utils.jl")
include("utils.jl")

const all_res = Dict{String,Vector{Any}}()

function collect_results(dir, all_res)
    for file in readdir(dir, join=true)
        if isdir(file)
            collect_results(file, all_res)
            continue
        end
        endswith(file, ".yaml") || continue
        for res in YAML.load_file(file)
            dev_res = get!(()->[], all_res, get_device_id(res))
            push!(dev_res, (nrep=res["nrep"], nele=res["nele"], nconcurrent=res["nconcurrent"],
                            queue_times=res["queue_times"], submit_times=res["submit_times"],
                            start_times=res["start_times"],
                            complete_times=res["complete_times"]))
        end
    end
end

collect_results(ARGS[1], all_res)

function write_time_array(io, times)
    prev = 0
    for t in times
        write_leb128(io, UInt(t - prev))
        prev = t
    end
end

function write_times(io, times)
    for ts in times
        write_time_array(io, ts)
    end
end

function output_result(dir, dev, res)
    mkpath(dir, mode=0o755)
    name = "$(dev)-$(res.nele)-$(res.nconcurrent)"
    name = joinpath(dir, name)
    open(name, "w") do io
        write_leb128(io, UInt64(res.nele))
        write_leb128(io, UInt64(res.nconcurrent))
        write_leb128(io, UInt64(res.nrep))
        write_times(io, res.queue_times)
        write_times(io, res.submit_times)
        write_times(io, res.start_times)
        write_times(io, res.complete_times)
    end
end

const dir = ARGS[2]

for (dev, dev_res) in all_res
    for res in dev_res
        output_result(dir, dev, res)
    end
end
