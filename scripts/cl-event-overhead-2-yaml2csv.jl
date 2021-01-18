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
                            nbuffer=res["nbuffer"], nworker=res["nworker"],
                            do_wait=get(res, "do_wait", false) ? 1 : 0,
                            complete_cb=res["complete_cb"] ? 1 : 0,
                            do_unmap=get(res, "do_unmap", false) ? 1 : 0))
        end
    end
end

collect_results(ARGS[1], all_res)

const dir = ARGS[2]
mkpath(dir, mode=0o755)

for (dev, dev_res) in all_res
    name = joinpath(dir, "$(dev).csv")
    open(name, "w") do io
        println(io, "Size,Repeat,Workers,Buffer Ring,Do Wait,Complete CB,Unmap,Time (ns)")
        sort!(dev_res, by=x->(x.nele, x.nrep, x.nworker, x.nbuffer,
                              x.do_wait, x.complete_cb, x.do_unmap, x.t))
        for res in dev_res
            println(io, ("$(res.nele),$(res.nrep),$(res.nworker),$(res.nbuffer)," *
                         "$(res.do_wait),$(res.complete_cb),$(res.do_unmap),$(res.t)"))
        end
    end
end
