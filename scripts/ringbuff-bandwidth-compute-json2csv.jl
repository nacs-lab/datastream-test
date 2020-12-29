#!/usr/bin/julia

using JSON

lines = readlines(ARGS[1])

all_res = []

for i in 1:2:length(lines)
    _size, blksz, nchn_wr, nchn_rd, localbuff = parse.(Int, split(lines[i], '-'))
    blockring = JSON.parse(lines[i + 1])
    push!(all_res, (size=_size, blksz=blksz, nchn_wr=nchn_wr, nchn_rd=nchn_rd,
                    localbuff=localbuff, read=blockring["read"], write=blockring["write"]))
end

sort!(all_res, by=x->(x.size, x.blksz, x.nchn_wr, x.nchn_rd, x.localbuff))

function print_res(io, all_res, getter)
    println(io, "Size (B),Block Size (B),Write Chn,Read Chn,Local Buffer Size (B)," *
            "Time (ns),Inst,Cycle,Cache Ref,Cache Miss,Pipe RW,Pipe Stall,Pipe Sync")
    for res in all_res
        data = getter(res)
        println(io, "$(res.size),$(res.blksz),$(res.nchn_wr),$(res.nchn_rd),$(res.localbuff)," *
                "$(data["t"]),$(data["inst"]),$(data["cycle"])," *
                "$(data["cache_ref"]),$(data["cache_miss"])," *
                "$(data["pipe_rw"]),$(data["pipe_stall"]),$(data["pipe_sync"])")
    end
end

const prefix = ARGS[2]

open(io->print_res(io, all_res, res->res.read), "$(prefix)_read.csv", "w")
open(io->print_res(io, all_res, res->res.write), "$(prefix)_write.csv", "w")
