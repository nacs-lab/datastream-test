#!/usr/bin/julia

using NaCsPlot
using PyPlot
using DelimitedFiles

include("utils.jl")

function preprocess_data(data)
    sizes = Int[]
    nblk = Int[]
    nchn_wr = Int[]
    nchn_rd = Int[]
    localbuff = Int[]
    byte_ns = Float64[]
    byte_cyl = Float64[]
    byte_ref = Float64[]
    miss_perc = Float64[]
    pipe_rw = Float64[]
    pipe_stall_perc = Float64[]
    pipe_sync = Float64[]
    for i in 1:size(data, 1)
        push!(sizes, data[i, 1])
        push!(nblk, data[i, 1] / data[i, 2])
        push!(nchn_wr, data[i, 3])
        push!(nchn_rd, data[i, 4])
        push!(localbuff, data[i, 5])
        push!(byte_ns, 4 / data[i, 6])
        push!(byte_cyl, 4 / data[i, 8])
        push!(byte_ref, 4 / data[i, 9])
        push!(miss_perc, 100 * data[i, 10] / data[i, 9])
        push!(pipe_rw, data[i, 11] / 4 * data[i, 2])
        push!(pipe_stall_perc, 100 * data[i, 12] / data[i, 11])
        push!(pipe_sync, data[i, 13] / data[i, 12])
    end
    return (size=sizes, nblk=nblk, nchn_wr=nchn_wr, nchn_rd=nchn_rd, localbuff=localbuff,
            byte_ns=byte_ns, byte_cyl=byte_cyl, byte_ref=byte_ref, miss_perc=miss_perc,
            pipe_rw=pipe_rw, pipe_stall_perc=pipe_stall_perc, pipe_sync=pipe_sync)
end

function load_data(rdname, wrname)
    return (read=preprocess_data(readdlm(rdname, ',', skipstart=1)),
            write=preprocess_data(readdlm(wrname, ',', skipstart=1)))
end

function filter_data(cb, data)
    sizes = Int[]
    nblk = Int[]
    nchn_wr = Int[]
    nchn_rd = Int[]
    localbuff = Int[]
    byte_ns = Float64[]
    byte_cyl = Float64[]
    byte_ref = Float64[]
    miss_perc = Float64[]
    pipe_rw = Float64[]
    pipe_stall_perc = Float64[]
    pipe_sync = Float64[]
    for i in 1:length(data.size)
        point = (size=data.size[i], nblk=data.nblk[i],
                 nchn_wr=data.nchn_wr[i], nchn_rd=data.nchn_rd[i], localbuff=data.localbuff[i],
                 byte_ns=data.byte_ns[i], byte_cyl=data.byte_cyl[i],
                 byte_ref=data.byte_ref[i], miss_perc=data.miss_perc[i],
                 pipe_rw=data.pipe_rw[i], pipe_stall_perc=data.pipe_stall_perc[i],
                 pipe_sync=data.pipe_sync[i])
        cb(point) || continue
        push!(sizes, point.size)
        push!(nblk, point.nblk)
        push!(nchn_wr, point.nchn_wr)
        push!(nchn_rd, point.nchn_rd)
        push!(localbuff, point.localbuff)
        push!(byte_ns, point.byte_ns)
        push!(byte_cyl, point.byte_cyl)
        push!(byte_ref, point.byte_ref)
        push!(miss_perc, point.miss_perc)
        push!(pipe_rw, point.pipe_rw)
        push!(pipe_stall_perc, point.pipe_stall_perc)
        push!(pipe_sync, point.pipe_sync)
    end
    return (size=sizes, nblk=nblk, nchn_wr=nchn_wr, nchn_rd=nchn_rd, localbuff=localbuff,
            byte_ns=byte_ns, byte_cyl=byte_cyl, byte_ref=byte_ref, miss_perc=miss_perc,
            pipe_rw=pipe_rw, pipe_stall_perc=pipe_stall_perc, pipe_sync=pipe_sync)
end

function plot_perf(data, nblks, localbuffs, nchn_rd, nchn_wr)
    # Note that this is very specific to how we took the data and isn't very generic...
    @assert length(localbuffs) == 2
    nchn = max(nchn_rd, nchn_wr)
    for i in 1:length(nblks)
        nblk = nblks[i]
        for localbuff in localbuffs
            line_rd = filter_data(x->(x.nblk == nblk && x.localbuff == localbuff &&
                                      x.nchn_rd == nchn_rd && x.nchn_wr == nchn_wr), data.read)
            line_wr = filter_data(x->(x.nblk == nblk && x.localbuff == localbuff &&
                                      x.nchn_rd == nchn_rd && x.nchn_wr == nchn_wr), data.write)
            @assert line_rd.size == line_wr.size
            isempty(line_rd.size) && continue
            # Make the color consistent with the previous measurements
            if localbuff == 0
                plot(line_rd.size, (line_rd.byte_ns .+ line_wr.byte_ns) ./ 2 ./ 4 .* nchn,
                     label="$nblk blocks", color="C$(i)", ls="-")
            else
                plot(line_rd.size, (line_rd.byte_ns .+ line_wr.byte_ns) ./ 2 ./ 4 .* nchn,
                     color="C$(i)", ls="--")
            end
        end
    end
    grid()
    legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
    title("Write $(nchn_wr)CH + Read $(nchn_rd)CH")
    xscale("log")
    xlabel("Buffer Size")
    ylabel("Throughput (S/ns)")
    ax = gca()
    xticks(2.0.^(18:2:22), size_to_str.(2.0.^(18:2:22)))
    ax.set_xticks(2.0.^(17:22), minor=true)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
end

function plot_stall(data, nblks, localbuffs, nchn_rd, nchn_wr)
    # Note that this is very specific to how we took the data and isn't very generic...
    @assert length(localbuffs) == 2
    nchn = max(nchn_rd, nchn_wr)
    for i in 1:length(nblks)
        nblk = nblks[i]
        for localbuff in localbuffs
            line_rd = filter_data(x->(x.nblk == nblk && x.localbuff == localbuff &&
                                      x.nchn_rd == nchn_rd && x.nchn_wr == nchn_wr), data.read)
            line_wr = filter_data(x->(x.nblk == nblk && x.localbuff == localbuff &&
                                      x.nchn_rd == nchn_rd && x.nchn_wr == nchn_wr), data.write)
            @assert line_rd.size == line_wr.size
            isempty(line_rd.size) && continue
            # Make the color consistent with the previous measurements
            if localbuff == 0
                plot(line_rd.size, line_rd.pipe_stall_perc,
                     label="$nblk blocks", color="C$(i)", ls="-")
                plot(line_rd.size, line_wr.pipe_stall_perc,
                     color="C$(i)", ls=":")
            else
                plot(line_rd.size, line_rd.pipe_stall_perc,
                     color="C$(i)", ls="--")
                plot(line_rd.size, line_wr.pipe_stall_perc,
                     color="C$(i)", ls="-.")
            end
        end
    end
    grid()
    legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
    title("Write $(nchn_wr)CH + Read $(nchn_rd)CH")
    xscale("log")
    xlabel("Buffer Size")
    ylabel("Pipe Stall (%)")
    ylim([0, 100])
    ax = gca()
    xticks(2.0.^(18:2:22), size_to_str.(2.0.^(18:2:22)))
    ax.set_xticks(2.0.^(17:22), minor=true)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
end
