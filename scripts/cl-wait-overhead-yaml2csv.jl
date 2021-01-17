#!/usr/bin/julia

using YAML

include("cl-utils.jl")

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
            push!(dev_res, (t=res["t"], nrep=res["nrep"], nele=res["nele"],
                            do_wait=res["do_wait"] ? 1 : 0))
        end
    end
end

collect_results(ARGS[1], all_res)

const dir = ARGS[2]
mkpath(dir, mode=0o755)

for (dev, dev_res) in all_res
    name = joinpath(dir, "$(dev).csv")
    open(name, "w") do io
        println(io, "Size,Repeat,Do Wait,Time (ns)")
        sort!(dev_res, by=x->(x.nele, x.nrep, x.do_wait, x.t))
        for res in dev_res
            println(io, ("$(res.nele),$(res.nrep),$(res.do_wait),$(res.t)"))
        end
    end
end
