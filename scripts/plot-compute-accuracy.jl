#!/usr/bin/julia

using NaCsCalc
using NaCsPlot
using PyPlot
using LibArchive

include("utils.jl")

const prefix = joinpath(@__DIR__, "../imgs/compute-accuracy")

struct LEB128
    v::Vector{UInt8}
end

function Base.iterate(vl::LEB128, state=1)
    if state > length(vl.v)
        return
    end
    v::UInt64 = 0
    shift = 0
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

Base.IteratorSize(::LEB128) = Base.SizeUnknown()

function decode_sign(v::UInt64)
    if v & 1 == 0
        return Int(v >> 1)
    else
        return -Int(v >> 1) - 1
    end
end

struct DiffData
    start::Float64
    step::Float64
    nsteps::Int
    diff::Vector{Int}
end

read_compressed_res(fname) = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    start = read(reader, Float64)
    step = read(reader, Float64)
    nsteps = Int(read(reader, UInt32))
    return DiffData(start, step, nsteps, Int[decode_sign(v) for v in LEB128(read(reader))])
end

const datadir = joinpath(@__DIR__, "../data/")
const cpudir = joinpath(datadir, "cpu-compute-accuracy")
const cldir = joinpath(datadir, "cl-compute-accuracy")

const cpunames = ["aarch64-scalar", "aarch64-asimd",
                  "x86_64-scalar", "x86_64-sse2", "x86_64-avx",
                  "x86_64-avx2", "x86_64-avx512"]

const clnames = ["amd-gpu-gfx1012-native", "amd-gpu-gfx1012",
                 "intel-core-i7-6700k-native", "intel-core-i7-6700k",
                 "intel-core-i9-10885h-native", "intel-core-i9-10885h",
                 "intel-gpu-1912-native", "intel-gpu-1912",
                 "intel-gpu-9bc4-native", "intel-gpu-9bc4"]

const cpudata = Dict(name=>read_compressed_res(joinpath(cpudir, "$(name).zst"))
                     for name in cpunames)
const cldata = Dict(name=>read_compressed_res(joinpath(cldir, "$(name).zst"))
                    for name in clnames)

function plot_diff(data::DiffData; xscale=1, kws...)
    xs = range(data.start, step=data.step .* xscale, length=data.nsteps)
    plot(xs, data.diff; kws...)
end

figure(figsize=[12.6, 16.8])

ax = subplot(3, 2, 1)
plot_diff(cpudata["aarch64-scalar"]; xscale=π, label="A64 Scalar")
plot_diff(cpudata["aarch64-asimd"]; xscale=π, label="A64 ASIMD")
plot_diff(cpudata["x86_64-scalar"]; xscale=π, label="x64 Scalar")
plot_diff(cpudata["x86_64-sse2"]; xscale=π, label="x64 SSE2")
plot_diff(cpudata["x86_64-avx"]; xscale=π, label="x64 AVX")
plot_diff(cpudata["x86_64-avx2"]; xscale=π, label="x64 AVX2")
plot_diff(cpudata["x86_64-avx512"]; xscale=π, label="x64 AVX512")
legend(fontsize="xx-small", ncol=4, columnspacing=0.8, handlelength=0.8,
       loc="upper center")
grid()
title("CPU")

ax = subplot(3, 2, 2)
plot_diff(cldata["amd-gpu-gfx1012-native"]; label="native_sin")
plot_diff(cldata["amd-gpu-gfx1012"]; label="sin")
legend(fontsize="xx-small", ncol=2, columnspacing=1, handlelength=1,
       loc="upper left")
grid()
title("OpenCL AMD RX 5500 XT")

ax = subplot(3, 2, 3)
plot_diff(cldata["intel-gpu-1912-native"]; label="native_sin")
plot_diff(cldata["intel-gpu-1912"]; label="sin")
legend(fontsize="xx-small", ncol=2, columnspacing=1, handlelength=1,
       loc="upper left")
grid()
title("OpenCL Intel UHD 530")

ax = subplot(3, 2, 4)
plot_diff(cldata["intel-gpu-9bc4-native"]; label="native_sin")
plot_diff(cldata["intel-gpu-9bc4"]; label="sin")
legend(fontsize="xx-small", ncol=2, columnspacing=1, handlelength=1,
       loc="upper left")
grid()
title("OpenCL Intel UHD 630")

ax = subplot(3, 2, 5)
plot_diff(cldata["intel-core-i7-6700k-native"]; label="native_sin")
plot_diff(cldata["intel-core-i7-6700k"]; label="sin")
legend(fontsize="xx-small", ncol=2, columnspacing=1, handlelength=1,
       loc="upper left")
grid()
title("OpenCL Intel CPU i7-6700K")

ax = subplot(3, 2, 6)
plot_diff(cldata["intel-core-i9-10885h-native"]; label="native_sin")
plot_diff(cldata["intel-core-i9-10885h"]; label="sin")
legend(fontsize="xx-small", ncol=2, columnspacing=1, handlelength=1,
       loc="upper left")
grid()
title("OpenCL Intel CPU i9-10885H")

tight_layout(pad=0.3)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
