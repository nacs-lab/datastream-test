#!/usr/bin/julia

using NaCsCalc
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

struct TimeInfo
    times::Vector{Float64}
    max_time::Dict{Int,Float64}
    TimeInfo(times) = new(times, Dict{Int,Float64}())
end

get_times(data, ghz) = TimeInfo(@view(data[end ÷ 16:end]) ./ ghz)

function load_all_times(dir, ghz)
    res = Dict{String,Vector{TimeInfo}}()
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

@inline get_max_time(timeinfo, samples) = get!(timeinfo.max_time, samples) do
    times = timeinfo.times
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
    return tmax
end

@inline get_speed_limit(timeinfo, samples) = samples / get_max_time(timeinfo, samples)

# get_speed_limits(times, samples) = get_speed_limit.(Ref(times), samples)

function plot_line(times, xscale, yscale; kws...)
    xs = Int[]
    x = length(times.times) - 1
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

struct HullPoint
    idx::Int
    val::Float64
    slope::Float64
end

function find_convex_hull(times)
    upper = [HullPoint(1, times[1], 0.0), HullPoint(2, times[2], times[2] - times[1])]
    lower = [HullPoint(1, times[1], 0.0), HullPoint(2, times[2], times[2] - times[1])]
    @inbounds for i in 3:length(times)
        t = times[i]
        new_slope_upper = (t - upper[end].val) / (i - upper[end].idx)
        while length(upper) > 1 && new_slope_upper > upper[end].slope
            pop!(upper)
            new_slope_upper = (t - upper[end].val) / (i - upper[end].idx)
        end
        push!(upper, HullPoint(i, t, new_slope_upper))
        new_slope_lower = (t - lower[end].val) / (i - lower[end].idx)
        while length(lower) > 1 && new_slope_lower < lower[end].slope
            pop!(lower)
            new_slope_lower = (t - lower[end].val) / (i - lower[end].idx)
        end
        push!(lower, HullPoint(i, t, new_slope_lower))
    end
    return upper, lower
end

# Find the max speed at which the convex hulls are still overlapping.
# Above this speed we are only checking the global properties
# and not the local properties anymore
function find_max_speed(times)
    upper, lower = find_convex_hull(times)
    # apart from the first and the last point there should be no points with the same index
    # in the upper and lower hull since each hull point must be a datapoint in the original set
    # and each point have unique index and can't appear in both hulls.
    idx_up = 2
    idx_lo = 2
    pt_up = upper[2]
    pt_lo = lower[2]
    @assert pt_up.slope > pt_lo.slope
    while true
        if pt_up.val < pt_lo.val
            idx_up += 1
            pt_up = upper[idx_up]
            if pt_up.slope > pt_lo.slope
                continue
            end
            return 1 / pt_lo.slope
        else
            idx_lo += 1
            pt_lo = lower[idx_lo]
            if pt_up.slope > pt_lo.slope
                continue
            end
            return 1 / pt_up.slope
        end
    end
end

function find_min_runahead(times, speed)
    max_output = 1 - speed * times[1]
    required_runahead = 0.0
    for i in 2:length(times)
        output = i - speed * times[i]
        if output > max_output
            max_output = output
        else
            runahead = max_output - output
            if runahead > required_runahead
                required_runahead = runahead
            end
        end
    end
    return required_runahead
end

function get_min_runaheads(timeinfo)
    speed_lb = get_speed_limit(timeinfo, 1)
    speed_ub = find_max_speed(timeinfo.times)
    # speed_ub = get_speed_limit(timeinfo, length(timeinfo.times) - 1)
    speeds = linspace(speed_lb, speed_ub, 201)
    return speeds, find_min_runahead.(Ref(timeinfo.times), speeds)
end

function plot_runahead(timeinfo, xscale, yscale; kws...)
    xs, ys = get_min_runaheads(timeinfo)
    plot(xs .* xscale, ys .* yscale; kws...)
end

function plot_runaheads(all_times, xscale, yscale, color, name; kws...)
    first = true
    for times in all_times
        if first
            first = false
            plot_runahead(times, xscale, yscale, color=color, label=name; kws...)
        else
            plot_runahead(times, xscale, yscale, color=color; kws...)
        end
    end
end
