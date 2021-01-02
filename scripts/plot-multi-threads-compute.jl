#!/usr/bin/julia

using NaCsPlot
using PyPlot
using LibArchive

include("utils.jl")

struct VLInt
    v::Vector{UInt8}
end

function Base.iterate(vl::VLInt, state=1)
    if state > length(vl.v)
        return
    end
    v16 = vl.v[state] | (UInt16(vl.v[state + 1]) << 8)
    state += 2
    if (v16 & 0x8000) == 0
        return v16 % UInt64, state
    end
    v::UInt64 = v16 & 0x7fff
    shift = 15
    while true
        v8 = vl.v[state]
        state += 1
        if v8 & 0x80 == 0
            return v | UInt64(v8) << shift, state
        end
        v = v | UInt64(v8 & 0x7f) << shift
        shift += 7
    end
end

Base.IteratorSize(::VLInt) = Base.SizeUnknown()

read_compressed_res(fname) = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    vs = collect(UInt64, VLInt(read(reader)))
    v0 = vs[1]
    vmin = vs[2]
    for i in 2:length(vs) - 1
        dv = vmin + vs[i + 1]
        v0 += dv
        vs[i] = v0
    end
    resize!(vs, length(vs) - 1)
    return vs
end

get_times(data, ghz) = @view(data[end ÷ 16:end]) ./ ghz

function load_all_times(dir, ghz)
    res = Dict{String,Vector{Vector{Float64}}}()
    for d in readdir(dir, join=true)
        isdir(d) || continue
        for f in readdir(d, join=true)
            isfile(f) || continue
            endswith(f, ".zst") || continue
            r = get!(()->Vector{Float64}[], res, basename(f)[1:end - 4])
            push!(r, get_times(read_compressed_res(f), ghz))
        end
    end
    return res
end

function get_speed_limit(times, samples)
    tmax = zero(eltype(times))
    for i in 1:length(times) - samples
        dt = times[i + samples] - times[i]
        if dt > tmax
            tmax = dt
        end
    end
    return samples / tmax
end

function plot_line(times, xscale, yscale; kws...)
    if length(times) <= 1025
        xs = [1:(length(times) - 1);]
    else
        xs = Int[1:1024;]
        x = length(times) - 1
        while x > 1024
            insert!(xs, 1025, x)
            x = floor(Int, x * 0.9)
        end
    end
    plot(xs .* xscale, get_speed_limit.(Ref(times), xs) .* yscale; kws...)
end

function plot_lines(all_times, xscale, yscale, color, name; kws...)
    first = true
    for times in all_times
        if first
            first = false
            plot_line(times, xscale, yscale, color=color, label=name; kws...)
        else
            plot_line(times, xscale, yscale, color=color; kws...)
        end
    end
end
