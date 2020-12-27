#!/usr/bin/julia

using NaCsPlot
using PyPlot
using DelimitedFiles

include("utils.jl")

function filter_data_core(data, core)
    sizes = Int[]
    byte_ns = Float64[]
    byte_cyl = Float64[]
    byte_ref = Float64[]
    miss_perc = Float64[]
    for i in 1:size(data, 1)
        if data[i, 1] != core
            continue
        end
        push!(sizes, data[i, 2])
        push!(byte_ns, 4 / data[i, 3])
        push!(byte_cyl, 4 / data[i, 5])
        push!(byte_ref, 4 / data[i, 6])
        push!(miss_perc, 100 * data[i, 7] / data[i, 6])
    end
    return (core=core, size=sizes, byte_ns=byte_ns, byte_cyl=byte_cyl,
            byte_ref=byte_ref, miss_perc=miss_perc)
end

function filter_data(cb, data)
    sizes = Int[]
    byte_ns = Float64[]
    byte_cyl = Float64[]
    byte_ref = Float64[]
    miss_perc = Float64[]
    for i in 1:length(data.size)
        point = (core=data.core, size=data.size[i], byte_ns=data.byte_ns[i],
                 byte_cyl=data.byte_cyl[i], byte_ref=data.byte_ref[i],
                 miss_perc=data.miss_perc[i])
        cb(point) || continue
        push!(sizes, point.size)
        push!(byte_ns, point.byte_ns)
        push!(byte_cyl, point.byte_cyl)
        push!(byte_ref, point.byte_ref)
        push!(miss_perc, point.miss_perc)
    end
    return (core=data.core, size=sizes, byte_ns=byte_ns, byte_cyl=byte_cyl,
            byte_ref=byte_ref, miss_perc=miss_perc)
end

function load_data(name)
    data = readdlm(name, ',', skipstart=1)
    cores = sort!(Int.(unique(data[:, 1])))
    return [filter_data_core(data, core) for core in cores]
end
