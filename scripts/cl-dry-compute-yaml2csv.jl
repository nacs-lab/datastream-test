#!/usr/bin/julia

using YAML

function get_intel_device_id(res)
    name = res["name"]
    m = match(r"Core\(TM\) (i[-0-9A-Za-z]*) ", name)
    if m !== nothing
        return "intel-core-" * lowercase(m[1])
    end
    m = match(r"Graphics \[0x([0-9A-Fa-f]*)\]", name)
    if m !== nothing
        return "intel-gpu-" * lowercase(m[1])
    end
    return name
end

function get_device_id(res)
    if res["vendor"] == "Intel(R) Corporation"
        return get_intel_device_id(res)
    elseif res["vendor"] == "Advanced Micro Devices, Inc."
        return "amd-$(res["type"])-$(res["name"])"
    end
    return res["name"]
end

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
            push!(dev_res, (tdummy=res["tdummy"], tcompute=res["tcompute"],
                            tcompute_native=get(res, "tcompute_native", NaN),
                            ooo=res["ooo"] ? 1 : 0, nrep=res["nrep"], nele=res["nele"],
                            ncalc=get(res, "ncalc", 1)))
        end
    end
end

collect_results(ARGS[1], all_res)

const prefix = ARGS[2]

for (dev, dev_res) in all_res
    if endswith(prefix, "/")
        name = "$(prefix)$(dev).csv"
    else
        name = "$(prefix)-$(dev).csv"
    end
    open(name, "w") do io
        println(io, "Size,Eval per Kernel,Repeat,Out of order," *
                "Dummy Time (ns),Compute Time (ns),Native Compute Time (ns)")
        sort!(dev_res, by=x->(x.nele, x.ncalc, x.nrep, x.ooo))
        for res in dev_res
            println(io, "$(res.nele),$(res.ncalc),$(res.nrep),$(res.ooo)," *
                    "$(res.tdummy),$(res.tcompute),$(res.tcompute_native)")
        end
    end
end
