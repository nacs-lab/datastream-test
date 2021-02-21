#!/usr/bin/julia

using NaCsPlot
using PyPlot
using DelimitedFiles
using Statistics

include("utils.jl")

function preprocess_data(data)
    return (nele=Int.(@view data[:, 1]),
            ncalc_single=Int.(@view data[:, 2]),
            ncalc_double=Int.(@view data[:, 3]),
            nrep=Int.(@view data[:, 4]),
            tdummy=data[:, 5],
            tcompute_in_order=data[:, 6],
            tcompute_out_of_order=data[:, 7])
end

load_data(name) = preprocess_data(readdlm(name, ',', skipstart=1))

function filter_data(cb, data)
    res = (nele=Int[],
           ncalc_single=Int[],
           ncalc_double=Int[],
           nrep=Int[],
           tdummy=Float64[],
           tcompute_in_order=Float64[],
           tcompute_out_of_order=Float64[])
    for i in 1:length(data.nele)
        point = (nele=data.nele[i], ncalc_single=data.ncalc_single[i],
                 ncalc_double=data.ncalc_double[i], nrep=data.nrep[i],
                 tdummy=data.tdummy[i], tcompute_in_order=data.tcompute_in_order[i],
                 tcompute_out_of_order=data.tcompute_out_of_order[i])
        cb(point) || continue
        push!(res.nele, point.nele)
        push!(res.ncalc_single, point.ncalc_single)
        push!(res.ncalc_double, point.ncalc_double)
        push!(res.nrep, point.nrep)
        push!(res.tdummy, point.tdummy)
        push!(res.tcompute_in_order, point.tcompute_in_order)
        push!(res.tcompute_out_of_order, point.tcompute_out_of_order)
    end
    return res
end

sort_u(vals) = sort(unique(vals))

const prefix = joinpath(@__DIR__, "../imgs/cl-mixed-compute")
const datadir = joinpath(@__DIR__, "../data/cl-mixed-compute")

const devices = ["amd-gpu-gfx1012", "intel-gpu-1912", "intel-gpu-9bc4"]
const data = Dict(dev=>load_data(joinpath(datadir, "$(dev).csv"))
                  for dev in devices)

function plot_compute(data, ooo)
    ncalc_doubles = sort_u(data.ncalc_double)
    cnt = 0
    xmin = Inf
    xmax = 0.0
    ymin = Inf
    ymax = 0.0
    for ncalc_double in ncalc_doubles
        line = filter_data(x->x.ncalc_double == ncalc_double, data)
        color = "C$(cnt)"
        cnt += 1

        x = line.ncalc_single
        y = ooo ? line.tcompute_out_of_order : line.tcompute_in_order

        x_u = sort_u(x)
        y_u = [mean((y[xi] for xi in 1:length(x) if x[xi] == vx)) for vx in x_u]

        xmin = min(xmin, minimum(x))
        xmax = max(xmax, maximum(x))
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))

        plot(x, y, ".", color=color)
        plot(x_u, y_u, color=color)
        plot([], [], ".", ls="-", color=color, label="$(ncalc_double)")
    end
    grid()
    xlabel("Single Evaluations per Kernel")
    ylabel("Time per Evaluation (ns)")
    legend(fontsize="xx-small", ncol=3, columnspacing=0.8, handlelength=0.8)
    title(ooo ? "Interleaved " : "Separated")
end

for dev in devices
    figure(figsize=[12.6, 5.6])

    d = data[dev]

    ax = subplot(1, 2, 1)
    plot_compute(d, false)

    ax = subplot(1, 2, 2)
    plot_compute(d, true)

    tight_layout(pad=0.6)
    NaCsPlot.maybe_save("$(prefix)_$(dev)")
end

NaCsPlot.maybe_show()
