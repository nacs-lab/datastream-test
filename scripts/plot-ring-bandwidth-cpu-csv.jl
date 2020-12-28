#!/usr/bin/julia

using NaCsPlot
using PyPlot
using DelimitedFiles
using PyCall
using Printf

@pyimport copy as pycopy

include("utils.jl")

function preprocess_data(data)
    cpu_wr = Int[]
    cpu_rd = Int[]
    byte_ns = Float64[]
    byte_cyl = Float64[]
    byte_ref = Float64[]
    miss_perc = Float64[]
    pipe_rw = Float64[]
    pipe_stall_perc = Float64[]
    pipe_sync = Float64[]
    for i in 1:size(data, 1)
        push!(cpu_wr, data[i, 1])
        push!(cpu_rd, data[i, 2])
        push!(byte_ns, 4 / data[i, 3])
        push!(byte_cyl, 4 / data[i, 5])
        push!(byte_ref, 4 / data[i, 6])
        push!(miss_perc, 100 * data[i, 7] / data[i, 6])
        push!(pipe_rw, data[i, 8] / 4 * data[i, 2])
        push!(pipe_stall_perc, 100 * data[i, 9] / data[i, 8])
        push!(pipe_sync, data[i, 10] / data[i, 9])
    end
    return (cpu_wr=cpu_wr, cpu_rd=cpu_rd, byte_ns=byte_ns, byte_cyl=byte_cyl,
            byte_ref=byte_ref, miss_perc=miss_perc,
            pipe_rw=pipe_rw, pipe_stall_perc=pipe_stall_perc, pipe_sync=pipe_sync)
end

function load_data(rdname, wrname)
    return (read=preprocess_data(readdlm(rdname, ',', skipstart=1)),
            write=preprocess_data(readdlm(wrname, ',', skipstart=1)))
end

function filter_data(cb, data)
    cpu_wr = Int[]
    cpu_rd = Int[]
    byte_ns = Float64[]
    byte_cyl = Float64[]
    byte_ref = Float64[]
    miss_perc = Float64[]
    pipe_rw = Float64[]
    pipe_stall_perc = Float64[]
    pipe_sync = Float64[]
    for i in 1:length(data.cpu_wr)
        point = (cpu_wr=data.cpu_wr[i], cpu_rd=data.cpu_rd[i],
                 byte_ns=data.byte_ns[i], byte_cyl=data.byte_cyl[i],
                 byte_ref=data.byte_ref[i], miss_perc=data.miss_perc[i],
                 pipe_rw=data.pipe_rw[i], pipe_stall_perc=data.pipe_stall_perc[i],
                 pipe_sync=data.pipe_sync[i])
        cb(point) || continue
        push!(cpu_wr, point.cpu_wr)
        push!(cpu_rd, point.cpu_rd)
        push!(byte_ns, point.byte_ns)
        push!(byte_cyl, point.byte_cyl)
        push!(byte_ref, point.byte_ref)
        push!(miss_perc, point.miss_perc)
        push!(pipe_rw, point.pipe_rw)
        push!(pipe_stall_perc, point.pipe_stall_perc)
        push!(pipe_sync, point.pipe_sync)
    end
    return (cpu_wr=cpu_wr, cpu_rd=cpu_rd, byte_ns=byte_ns, byte_cyl=byte_cyl,
            byte_ref=byte_ref, miss_perc=miss_perc,
            pipe_rw=pipe_rw, pipe_stall_perc=pipe_stall_perc, pipe_sync=pipe_sync)
end

function get_cpu2cpu(data, ncores)
    @assert data.read.cpu_wr == data.write.cpu_wr
    @assert data.read.cpu_rd == data.write.cpu_rd
    res = fill(NaN, ncores, ncores)
    for i in 1:length(data.read.cpu_wr)
        cpu_wr = data.read.cpu_wr[i] + 1
        cpu_rd = data.read.cpu_rd[i] + 1
        @assert isnan(res[cpu_wr, cpu_rd])
        res[cpu_wr, cpu_rd] = (1024 / data.read.byte_ns[i] + 1024 / data.write.byte_ns[i]) / 2
    end
    return res
end

function plot_cpu2cpu(data)
    cmap = pycopy.copy(matplotlib.cm.get_cmap("viridis"))
    cmap.set_bad(color="lightgray")
    imshow(data, cmap=cmap)
    ncores = size(data, 1)
    xticks(0:(ncores - 1))
    yticks(0:(ncores - 1))
    colorbar()
    xlabel("Reader Core")
    ylabel("Writer Core")
    fs = 110 / ncores
    for x in 1:ncores
        for y in 1:ncores
            v = data[x, y]
            if isnan(v)
                text(x - 1, y - 1, "n/a",
                     ha="center", va="center", color="w", fontsize=fs)
            else
                text(x - 1, y - 1, @sprintf("%.1f", data[y, x]),
                     ha="center", va="center", color="magenta", fontsize=fs)
            end
        end
    end
end
