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

@inline _max(v1, v2) = v1 > v2 ? v1 : v2

# Somehow these are actually slower
# Also LLVM can't vectorize this right now...
# @inline _max(v1::Float32, v2::Float32) =
#     ccall("llvm.maxnum.f32", llvmcall, Float32, (Float32, Float32), v1, v2)
# @inline _max(v1::Float64, v2::Float64) =
#     ccall("llvm.maxnum.f64", llvmcall, Float64, (Float64, Float64), v1, v2)
# @inline function get_speed_limit(times::Vector{Float64}, samples)
#     tmaxv = (VecElement(0.0), VecElement(0.0), VecElement(0.0), VecElement(0.0))
#     imax = length(times) - samples
#     @inbounds for i in 1:4:(imax - 3)
#         dt1 = times[i + samples] - times[i]
#         dt2 = times[i + samples + 1] - times[i + 1]
#         dt3 = times[i + samples + 2] - times[i + 2]
#         dt4 = times[i + samples + 3] - times[i + 3]
#         tmaxv = ccall("llvm.maxnum.v4f64", llvmcall, NTuple{4,VecElement{Float64}},
#                       (NTuple{4,VecElement{Float64}}, NTuple{4,VecElement{Float64}}),
#                       tmaxv, (VecElement(dt1), VecElement(dt2),
#                               VecElement(dt3), VecElement(dt4)))
#     end
#     tmax = ccall("llvm.experimental.vector.reduce.fmax.v4f64", llvmcall, Float64,
#                  (NTuple{4,VecElement{Float64}},), tmaxv)
#     # tmax = _max(_max(tmaxv[1].value, tmaxv[2].value), _max(tmaxv[3].value, tmaxv[4].value))
#     @inbounds @simd ivdep for i in (imax - 2):imax
#         dt = times[i + samples] - times[i]
#         tmax = _max(dt, tmax)
#     end
#     return samples / tmax
# end

@inline function get_speed_limit(times, samples)
    tmax1 = zero(eltype(times))
    tmax2 = zero(eltype(times))
    tmax3 = zero(eltype(times))
    tmax4 = zero(eltype(times))
    imax = length(times) - samples
    @inbounds for i in 1:4:(imax - 3)
        dt1 = times[i + samples] - times[i]
        dt2 = times[i + samples + 1] - times[i + 1]
        dt3 = times[i + samples + 2] - times[i + 2]
        dt4 = times[i + samples + 3] - times[i + 3]
        tmax1 = _max(dt1, tmax1)
        tmax2 = _max(dt2, tmax2)
        tmax3 = _max(dt3, tmax3)
        tmax4 = _max(dt4, tmax4)
    end
    tmax = _max(_max(tmax1, tmax2), _max(tmax3, tmax4))
    @inbounds @simd ivdep for i in _max(imax - 2, 1):imax
        dt = times[i + samples] - times[i]
        tmax = _max(dt, tmax)
    end
    return samples / tmax
end

# get_speed_limits(times, samples) = get_speed_limit.(Ref(times), samples)

function plot_line(times, xscale, yscale; kws...)
    xs = Int[]
    x = length(times) - 1
    while x > 0
        insert!(xs, 1, x)
        xnew = floor(Int, x * 0.95)
        if xnew == x
            xnew = x - 1
        end
        x = xnew
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
