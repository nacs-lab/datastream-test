#!/usr/bin/julia

using NaCsPlot
using PyPlot
using DelimitedFiles

include("utils.jl")

function preprocess_data(data)
    nele = Int.(@view data[:, 1])
    nrep = Int.(@view data[:, 2])
    nwrite = Int.(@view data[:, 3])
    t = data[:, 4]
    return (nele=nele, nrep=nrep, nwrite=nwrite, t=t)
end

load_data(name) = preprocess_data(readdlm(name, ',', skipstart=1))

function filter_data(cb, data)
    nele = Int[]
    nwrite = Int[]
    nrep = Int[]
    t = Float64[]
    for i in 1:length(data.nele)
        point = (nele=data.nele[i], nwrite=data.nwrite[i],
                 nrep=data.nrep[i], t=data.t[i])
        cb(point) || continue
        push!(nele, point.nele)
        push!(nwrite, point.nwrite)
        push!(nrep, point.nrep)
        push!(t, point.t)
    end
    return (nele=nele, nwrite=nwrite, nrep=nrep, t=t)
end

sort_u(vals) = sort(unique(vals))

const prefix = joinpath(@__DIR__, "../imgs/cl-write-bandwidth")
const datadir = joinpath(@__DIR__, "../data/cl-write-bandwidth")

const devices = ["amd-gpu-gfx1012", "intel-gpu-1912", "intel-gpu-9bc4"]
const data = Dict(dev=>load_data(joinpath(datadir, "$(dev).csv"))
                  for dev in devices)

function plot_dummy(data)
    neles = 2 .^ (10:23)
    cnt_mid = 0
    cnt_big = 0
    xmin = Inf
    xmax = 0.0
    ymin = Inf
    ymax = 0.0
    for nele in neles
        line = filter_data(x->nele - 2 <= x.nele <= nele + 2, data)
        if nele <= 2^16
            style = "C$(cnt_mid)--"
            cnt_mid += 1
        else
            style = "C$(cnt_big)-"
            cnt_big += 1
        end
        perm = sortperm(line.nrep)
        x = line.nrep[perm]
        y = 4e9 ./ line.t[perm]
        xmin = min(xmin, minimum(x))
        xmax = max(xmax, maximum(x))
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
        plot(x, y, style, label="$(size_to_str(nele * 4))")
    end
    xscale("log")
    yscale("log")
    grid()
    xlabel("Repetition")
    ylabel("Memory Throughput (B/s)")
    legend(fontsize="xx-small", ncol=6, columnspacing=0.8, handlelength=1,
           loc="lower center")
    # title(ooo ? "Out of Order dummy" : "In Order dummy")

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

function plot_compute(data)
    neles = 2 .^ (10:23)
    cnt_mid = 0
    cnt_big = 0
    # xmin = Inf
    # xmax = 0.0
    ymin = Inf
    ymax = 0.0
    for nele in neles
        line = filter_data(x->nele - 2 <= x.nele <= nele + 2, data)
        if nele <= 2^16
            style = "C$(cnt_mid)--"
            cnt_mid += 1
        else
            style = "C$(cnt_big)-"
            cnt_big += 1
        end
        perm = sortperm(line.nwrite)
        x = line.nwrite[perm] .* 4
        y = 4e9 ./ line.t[perm]
        # xmin = min(xmin, minimum(x))
        # xmax = max(xmax, maximum(x))
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
        plot(x, y, style, label="$(size_to_str(nele * 4))")
    end
    yscale("log")
    grid()
    xlabel("Memory Write per Kernel (B)")
    ylabel("Memory Throughput (B/s)")
    legend(fontsize="xx-small", ncol=6, columnspacing=0.8, handlelength=1,
           loc="lower center")
    # title((ooo ? "Out of Order " : "In Order ") * (native ? "native_sin" : "sin"))

    ax = gca()

    # n_lo = floor(Int, log2(xmin) + 0.25)
    # n_hi = ceil(Int, log2(xmax) - 0.25)
    # xticks(2.0.^(n_lo:2:n_hi), size_to_str.(2.0.^(n_lo:2:n_hi)))
    # ax.set_xticks(2.0.^(n_lo:n_hi), minor=true)
    # ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    y_lo = floor(Int, log2(ymin) + 0.25)
    y_hi = ceil(Int, log2(ymax) - 0.25)
    yticks(2.0.^(y_lo:2:y_hi), size_to_str.(2.0.^(y_lo:2:y_hi)))
    ax.set_yticks(2.0.^(y_lo:y_hi), minor=true)
    ax.tick_params(axis="y", right=true, which="minor", direction="in")
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
end

for dev in devices
    figure(figsize=[12.6, 5.6])

    d = data[dev]

    ax = subplot(1, 2, 1)
    plot_dummy(d)

    ax = subplot(1, 2, 2)
    plot_compute(d)

    tight_layout(pad=0.6)
    NaCsPlot.maybe_save("$(prefix)_$(dev)")
end

NaCsPlot.maybe_show()
