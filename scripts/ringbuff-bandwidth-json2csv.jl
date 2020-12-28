#!/usr/bin/julia

using JSON

lines = readlines(ARGS[1])

all_res = []

for i in 1:3:length(lines)
    buff, blksz = parse.(Int, split(lines[i], '-'))
    datapipe = JSON.parse(lines[i + 1])
    blockring = JSON.parse(lines[i + 2])
    push!(all_res, (size=buff, blksz=blksz, datapipe=datapipe, blockring=blockring))
end

sort!(all_res, by=x->(x.size, x.blksz))

function print_res(io, all_res, getter)
    println(io, "Size (B),Block Size (B),Time (ns),Inst,Cycle,Cache Ref,Cache Miss,Pipe RW,Pipe Stall,Pipe Sync")
    for res in all_res
        data = getter(res)
        println(io, "$(res.size),$(res.blksz),$(data["t"]),$(data["inst"]),$(data["cycle"])," *
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
