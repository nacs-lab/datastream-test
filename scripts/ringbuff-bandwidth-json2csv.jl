#!/usr/bin/julia

using JSON

lines = readlines(ARGS[1])

all_res = []

for i in 1:3:length(lines)
    buff = parse.(Int, lines[i])
    datapipe = JSON.parse(lines[i + 1])
    blockring = JSON.parse(lines[i + 2])
    push!(all_res, (size=buff, datapipe=datapipe, blockring=blockring))

    "cycle" => 1.05632, "cache_ref" => 0.135487, "pipe_stall" => 2.22455e-7, "cache_miss" => 0.104959, "pipe_rw" => 2.38419e-7, "t" => 0.390561, "inst" => 0.500752, "pipe_sync" => 3.21039e-6
end

sort!(all_res, by=x->x.size)

function print_res(io, all_res, getter)
    println(io, "Size (B),Time (ns),Inst,Cycle,Cache Ref,Cache Miss,Pipe RW,Pipe Stall,Pipe Sync")
    for res in all_res
        data = getter(res)
        println(io, "$(res.size),$(data["t"]),$(data["inst"]),$(data["cycle"])," *
                "$(data["cache_ref"]),$(data["cache_miss"])," *
                "$(data["pipe_rw"]),$(data["pipe_stall"]),$(data["pipe_sync"])")
    end
end

const prefix = ARGS[2]

open(io->print_res(io, all_res, res->res.datapipe["read"]),
     "$(prefix)_datapipe_read.csv", "w")
open(io->print_res(io, all_res, res->res.datapipe["write"]),
     "$(prefix)_datapipe_write.csv", "w")
open(io->print_res(io, all_res, res->res.blockring["read"]),
     "$(prefix)_blockring_read.csv", "w")
open(io->print_res(io, all_res, res->res.blockring["write"]),
     "$(prefix)_blockring_write.csv", "w")
