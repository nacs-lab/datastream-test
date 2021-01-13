#!/usr/bin/julia

using NaCsCalc
using NaCsData
using NaCsData.Fitting: fit_data
using NaCsPlot
using PyPlot
using DelimitedFiles
using Statistics

include("utils.jl")

function preprocess_data(data)
    nele = Int.(@view data[:, 1])
    nrep = Int.(@view data[:, 2])
    nstream = Int.(@view data[:, 3])
    t = data[:, 4]
    return (nele=nele, nrep=nrep, nstream=nstream, t=t)
end

load_data(name) = preprocess_data(readdlm(name, ',', skipstart=1))

function filter_data(cb, data)
    nele = Int[]
    nrep = Int[]
    nstream = Int[]
    t = Float64[]
    for i in 1:length(data.nele)
        point = (nele=data.nele[i], nrep=data.nrep[i], nstream=data.nstream[i], t=data.t[i])
        cb(point) || continue
        push!(nele, point.nele)
        push!(nrep, point.nrep)
        push!(nstream, point.nstream)
        push!(t, point.t)
    end
    return (nele=nele, nrep=nrep, nstream=nstream, t=t)
end

sort_u(vals) = sort(unique(vals))

const prefix = joinpath(@__DIR__, "../imgs/cl-parallel-rw")
const datadir = joinpath(@__DIR__, "../data/cl-parallel-rw")

const devices = ["amd-gpu-gfx1012", "intel-gpu-1912", "intel-gpu-9bc4"]
const data = Dict(dev=>load_data(joinpath(datadir, "$(dev).csv"))
                  for dev in devices)
const max_bandwidths = Dict("amd-gpu-gfx1012"=>224e9,
                            "intel-gpu-1912"=>34.1e9,
                            "intel-gpu-9bc4"=>45.8e9)

function plot_nrep(data, max_bw, nstream)
    neles = 2 .^ (14:23)
    cnt_mid = 4
    cnt_big = 0
    xmin = Inf
    xmax = 0.0
    ymin = max_bw
    ymax = max_bw
    for nele in neles
        line = filter_data(x->x.nstream == nstream && x.nele == nele, data)
        if nele <= 2^16
            color = "C$(cnt_mid)"
            linestyle = "--"
            cnt_mid += 1
        else
            color = "C$(cnt_big)"
            linestyle = "-"
            cnt_big += 1
        end
        perm = sortperm(line.nrep)
        x = line.nrep[perm]
        y = 4e9 ./ line.t[perm] .* nstream

        x_u = sort_u(x)
        y_u = [mean((y[xi] for xi in 1:length(x) if x[xi] == vx)) for vx in x_u]

        xmin = min(xmin, minimum(x))
        xmax = max(xmax, maximum(x))
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
        plot(x, y, ".", color=color)
        plot(x_u, y_u, linestyle, color=color)
        plot([], [], ".", ls=linestyle, color=color, label="$(size_to_str(nele * 4))")
    end
    axhline(max_bw, color="r", ls="dotted", linewidth=2)
    axhline(max_bw / 2, color="g", ls="dotted", linewidth=2)
    xscale("log")
    yscale("log")
    grid()
    xlabel("Repetition")
    ylabel("Memory Throughput (B/s)")
    legend(fontsize="xx-small", ncol=5, columnspacing=0.8, handlelength=1,
           loc="lower center")

    ax = gca()

    x_lo = floor(Int, log2(xmin) + 0.25)
    x_hi = ceil(Int, log2(xmax) - 0.25)
    xticks(2.0.^(x_lo:x_hi), size_to_str.(2.0.^(x_lo:x_hi)))
    ax.set_xticks(2.0.^(x_lo:x_hi), minor=true)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    y_lo = floor(Int, log2(ymin) + 0.25)
    y_hi = ceil(Int, log2(ymax) - 0.25)
    yticks(2.0.^(y_lo:2:y_hi), size_to_str.(2.0.^(y_lo:2:y_hi)))
    ax.set_yticks(2.0.^(y_lo:y_hi), minor=true)
    ax.tick_params(axis="y", right=true, which="minor", direction="in")
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
end

function plot_all_nrep(d, max_bw)
    nstreams = sort_u(d.nstream)
    ncol = 2
    nrow = (length(nstreams) - 1) ÷ ncol + 1
    figure(figsize=[6.3 * ncol, 5.6 * nrow])
    for i in 1:length(nstreams)
        ax = subplot(nrow, ncol, i)
        plot_nrep(d, max_bw, nstreams[i])
        title("$(nstreams[i]) Streams")
    end
    tight_layout(pad=0.6)
end

function plot_nstream(data, max_bw, nrep)
    neles = 2 .^ (14:23)
    cnt_mid = 4
    cnt_big = 0
    xmin = Inf
    xmax = 0.0
    ymin = max_bw
    ymax = max_bw
    for nele in neles
        line = filter_data(x->x.nrep == nrep && x.nele == nele, data)
        if nele <= 2^16
            color = "C$(cnt_mid)"
            linestyle = "--"
            cnt_mid += 1
        else
            color = "C$(cnt_big)"
            linestyle = "-"
            cnt_big += 1
        end
        perm = sortperm(line.nstream)
        x = line.nstream[perm]
        y = 4e9 ./ line.t[perm] .* x

        x_u = sort_u(x)
        y_u = [mean((y[xi] for xi in 1:length(x) if x[xi] == vx)) for vx in x_u]

        xmin = min(xmin, minimum(x))
        xmax = max(xmax, maximum(x))
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
        plot(x, y, ".", color=color)
        plot(x_u, y_u, linestyle, color=color)
        plot([], [], ".", ls=linestyle, color=color, label="$(size_to_str(nele * 4))")
    end
    axhline(max_bw, color="r", ls="dotted", linewidth=2)
    axhline(max_bw / 2, color="g", ls="dotted", linewidth=2)
    xscale("log")
    yscale("log")
    grid()
    xlabel("Streams")
    ylabel("Memory Throughput (B/s)")
    legend(fontsize="xx-small", ncol=5, columnspacing=0.8, handlelength=1,
           loc="lower center")

    ax = gca()

    x_lo = floor(Int, log2(xmin) + 0.25)
    x_hi = ceil(Int, log2(xmax) - 0.25)
    xticks(2.0.^(x_lo:x_hi), size_to_str.(2.0.^(x_lo:x_hi)))
    ax.set_xticks(2.0.^(x_lo:x_hi), minor=true)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    y_lo = floor(Int, log2(ymin) + 0.25)
    y_hi = ceil(Int, log2(ymax) - 0.25)
    yticks(2.0.^(y_lo:2:y_hi), size_to_str.(2.0.^(y_lo:2:y_hi)))
    ax.set_yticks(2.0.^(y_lo:y_hi), minor=true)
    ax.tick_params(axis="y", right=true, which="minor", direction="in")
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
end

function plot_all_nstream(d, max_bw)
    nreps = sort_u(d.nrep)
    ncol = 2
    nrow = (length(nreps) - 1) ÷ ncol + 1
    figure(figsize=[6.3 * ncol, 5.6 * nrow])
    for i in 1:length(nreps)
        ax = subplot(nrow, ncol, i)
        plot_nstream(d, max_bw, nreps[i])
        title("$(nreps[i]) Repetitions")
    end
    tight_layout(pad=0.6)
end

for dev in devices
    d = data[dev]
    max_bw = max_bandwidths[dev]

    plot_all_nrep(d, max_bw)
    NaCsPlot.maybe_save("$(prefix)_$(dev)_nrep")

    plot_all_nstream(d, max_bw)
    NaCsPlot.maybe_save("$(prefix)_$(dev)_nstream")
end

NaCsPlot.maybe_show()
