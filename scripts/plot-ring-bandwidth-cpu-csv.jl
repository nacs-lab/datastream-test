#!/usr/bin/julia

using NaCsPlot
using PyPlot
using DelimitedFiles

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
