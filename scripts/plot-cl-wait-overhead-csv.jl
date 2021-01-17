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
    return (nele=Int.(@view data[:, 1]),
            nrep=Int.(@view data[:, 2]),
            do_wait=(@view data[:, 3]) .!= 0,
            t=data[:, 4])
end

load_data(name) = preprocess_data(readdlm(name, ',', skipstart=1))

function filter_data(cb, data)
    res = (nele=Int[], nrep=Int[], t=Float64[], do_wait=Bool[])
    for i in 1:length(data.nele)
        point = (nele=data.nele[i], nrep=data.nrep[i], t=data.t[i], do_wait=data.do_wait[i])
        cb(point) || continue
        push!(res.nele, point.nele)
        push!(res.nrep, point.nrep)
        push!(res.do_wait, point.do_wait)
        push!(res.t, point.t)
    end
    return res
end

sort_u(vals) = sort(unique(vals))

const prefix = joinpath(@__DIR__, "../imgs/cl-event-overhead")
const datadir = joinpath(@__DIR__, "../data/cl-event-overhead")

const devices = ["amd-gpu-gfx1012", "intel-gpu-1912", "intel-gpu-9bc4"]
const data = Dict(dev=>load_data(joinpath(datadir, "$(dev).csv"))
                  for dev in devices)
const max_bandwidths = Dict("amd-gpu-gfx1012"=>224e9,
                            "intel-gpu-1912"=>34.1e9,
                            "intel-gpu-9bc4"=>45.8e9)

function plot_nele(data)
    xmin = Inf
    xmax = 0.0
    ymin = Inf
    ymax = 0.0
    for i in 0:7
        do_wait = (i & 1) != 0
        line = filter_data(x->x.do_wait == do_wait, data)
        x = line.nele .* 4
        y = line.t .* 1e-9 .* line.nele

        x_u = sort_u(x)
        y_u = [mean((y[xi] for xi in 1:length(x) if x[xi] == vx)) for vx in x_u]

        xmin = min(xmin, minimum(x))
        xmax = max(xmax, maximum(x))
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))

        color = "C$i"

        plot(x, y, ".", color=color)
        plot(x_u, y_u, color=color)
        plot([], [], ".", ls="-", color=color, label=do_wait ? "wait" : "no wait")
    end
    # axhline(max_bw, color="r", ls="dotted", linewidth=2)
    # axhline(max_bw / 2, color="g", ls="dotted", linewidth=2)
    grid()
    legend(fontsize="xx-small", ncol=2, columnspacing=0.8, handlelength=1,
           loc="upper left")

    ax = gca()

    xscale("log")
    xlabel("Buffer Size (B)")
    x_lo = floor(Int, log2(xmin) + 0.25)
    x_hi = ceil(Int, log2(xmax) - 0.25)
    xticks(2.0.^(x_lo:2:x_hi), size_to_str.(2.0.^(x_lo:2:x_hi)))
    ax.set_xticks(2.0.^(x_lo:x_hi), minor=true)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    if ymin > 0
        yscale("log")
        ylabel("Time per Run")
        y_lo = floor(Int, log10(ymin) + 0.25)
        y_hi = ceil(Int, log10(ymax) - 0.25)
        yticks(10.0.^(y_lo:y_hi), num_to_si.(10.0.^(y_lo:y_hi), Ref("s")))
        # ax.set_yticks(10.0.^(y_lo:y_hi), minor=true)
        ax.tick_params(axis="y", right=true, which="minor", direction="in")
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else
        ylabel("Time per Run (s)")
    end
end

for dev in devices
    d = data[dev]
    max_bw = max_bandwidths[dev]

    figure()
    plot_nele(d)
    tight_layout(pad=0.6)
    NaCsPlot.maybe_save("$(prefix)_$(dev)")
end

NaCsPlot.maybe_show()
