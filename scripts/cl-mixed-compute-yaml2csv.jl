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
            push!(dev_res, (tdummy=res["tdummy"], tcompute_in_order=res["tcompute_in_order"],
                            tcompute_out_of_order=res["tcompute_out_of_order"],
                            nrep=res["nrep"], nele=res["nele"],
                            ncalc_single=res["ncalc_single"],
                            ncalc_double=res["ncalc_double"]))
        end
    end
end

collect_results(ARGS[1], all_res)

const dir = ARGS[2]
mkpath(dir, mode=0o755)

for (dev, dev_res) in all_res
    name = joinpath(dir, "$(dev).csv")
    open(name, "w") do io
        println(io, "Size,Single Eval per Kernel,Double Eval per Kernel,Repeat," *
                "Dummy Time (ns),Compute Time In Order (ns),Compute Time Out of Order (ns)")
        sort!(dev_res, by=x->(x.nele, x.ncalc_single, x.ncalc_double, x.nrep))
        for res in dev_res
            println(io, "$(res.nele),$(res.ncalc_single),$(res.ncalc_double),$(res.nrep)," *
                    "$(res.tdummy),$(res.tcompute_in_order),$(res.tcompute_out_of_order)")
        end
    end
end
