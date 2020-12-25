#!/usr/bin/julia

using JSON

lines = readlines(ARGS[1])

all_res = []

for i in 1:2:length(lines)
    cores, buff = parse.(Int, split(lines[i], '-'))
    cores += 1
    res = JSON.parse(lines[i + 1])
    push!(all_res, (cores=cores, size=buff, t=res["t"], inst=res["inst"], cycle=res["cycle"],
                    cache_ref=res["cache_ref"], cache_miss=res["cache_miss"]))
end

sort!(all_res, by=x->(x.cores, x.size))

println("Cores,Size (B),Time (ns),Inst,Cycle,Cache Ref,Cache Miss"); for res in all_res
    println("$(res.cores),$(res.size),$(res.t),$(res.inst),$(res.cycle),$(res.cache_ref),$(res.cache_miss)")
end
