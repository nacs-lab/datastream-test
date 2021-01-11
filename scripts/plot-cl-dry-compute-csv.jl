#!/usr/bin/julia

using NaCsPlot
using PyPlot
using DelimitedFiles

include("utils.jl")

function preprocess_data(data)
    nele = Int.(@view data[:, 1])
    ncalc = Int.(@view data[:, 2])
    nrep = Int.(@view data[:, 3])
    ooo = @view(data[:, 4]) .!= 0
    tdummy = data[:, 5]
    tcompute = data[:, 6]
    tcompute_native = data[:, 7]
    return (nele=nele, ncalc=ncalc, nrep=nrep, ooo=ooo,
            tdummy=tdummy, tcompute=tcompute, tcompute_native=tcompute_native)
end

load_data(name) = preprocess_data(readdlm(name, ',', skipstart=1))

function filter_data(cb, data)
    nele = Int[]
    ncalc = Int[]
    nrep = Int[]
    ooo = Bool[]
    tdummy = Float64[]
    tcompute = Float64[]
    tcompute_native = Float64[]
    for i in 1:length(data.nele)
        point = (nele=data.nele[i], ncalc=data.ncalc[i],
                 nrep=data.nrep[i], ooo=data.ooo[i],
                 tdummy=data.tdummy[i], tcompute=data.tcompute[i],
                 tcompute_native=data.tcompute_native[i])
        cb(point) || continue
        push!(nele, point.nele)
        push!(ncalc, point.ncalc)
        push!(nrep, point.nrep)
        push!(ooo, point.ooo)
        push!(tdummy, point.tdummy)
        push!(tcompute, point.tcompute)
        push!(tcompute_native, point.tcompute_native)
    end
    return (nele=nele, ncalc=ncalc, nrep=nrep, ooo=ooo,
            tdummy=tdummy, tcompute=tcompute, tcompute_native=tcompute_native)
end

sort_u(vals) = sort(unique(vals))

const prefix = joinpath(@__DIR__, "../imgs/cl-dry-compute")
const datadir = joinpath(@__DIR__, "../data/cl-dry-compute")

const devices = ["amd-gpu-gfx1012", "intel-core-i7-6700k", "intel-core-i9-10885h",
                 "intel-gpu-1912", "intel-gpu-9bc4"]
const data = Dict(dev=>load_data(joinpath(datadir, "$(dev).csv"))
                  for dev in devices)

function plot_dummy(data, ooo)
    data = filter_data(x->x.ooo == ooo, data)
    neles = sort_u(data.nele)
    cnt_small = 0
    cnt_mid = 0
    cnt_big = 0
    for nele in neles
        line = filter_data(x->x.nele == nele, data)
        if nele <= 2^7
            style = "C$(cnt_small)-"
            cnt_small += 1
        elseif nele <= 2^15
            style = "C$(cnt_mid)^"
            cnt_mid += 1
        else
            style = "C$(cnt_big)o"
            cnt_big += 1
        end
        plot(line.nrep, line.tdummy ./ 1e9, style, label="$(size_to_str(nele))")
    end
    xscale("log")
    yscale("log")
    grid()
    xlabel("Repetition")
    ylabel("Time per Element")
    legend(fontsize="xx-small", ncol=3, columnspacing=0.8, handlelength=0.8,
           loc="lower right")
    title(ooo ? "Out of Order dummy" : "In Order dummy")

    ax = gca()

    nele_lo = floor(Int, log2(minimum(data.nrep)) + 0.25)
    nele_hi = ceil(Int, log2(maximum(data.nrep)) - 0.25)
    xticks(2.0.^(nele_lo:3:nele_hi), size_to_str.(2.0.^(nele_lo:3:nele_hi)))
    ax.set_xticks(2.0.^(nele_lo:nele_hi), minor=true)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(axis="y", right=true, which="minor", direction="in")

    t_lo = floor(Int, log10(minimum(data.tdummy)) + 0.25) - 9
    t_hi = ceil(Int, log10(maximum(data.tdummy)) - 0.25) - 9
    yticks(10.0.^(t_lo:t_hi), num_to_si.(10.0.^(t_lo:t_hi), Ref("s")))
end

function plot_compute(data, ooo, native)
    data = filter_data(x->x.ooo == ooo, data)
    neles = sort_u(data.nele)
    cnt_mid = 0
    cnt_big = 0
    xmin = Inf
    xmax = 0.0
    ymin = Inf
    ymax = 0.0
    for nele in neles
        if nele <= 2^7
            continue
        end
        line = filter_data(x->x.nele == nele, data)
        if nele <= 2^15
            style = "C$(cnt_mid)^"
            cnt_mid += 1
        else
            style = "C$(cnt_big)o"
            cnt_big += 1
        end
        t = native ? line.tcompute_native : line.tcompute
        t = t ./ line.ncalc
        xmin = min(xmin, minimum(line.ncalc))
        xmax = max(xmax, maximum(line.ncalc))
        ymin = min(ymin, minimum(t))
        ymax = max(ymax, maximum(t))
        plot(line.ncalc, t ./ 1e9, style, label="$(size_to_str(nele))")
    end
    xscale("log")
    yscale("log")
    grid()
    xlabel("Evaluations per Kernel")
    ylabel("Time per Evaluation")
    legend(fontsize="xx-small", ncol=2, columnspacing=0.8, handlelength=0.8,
           loc="upper right")
    title((ooo ? "Out of Order " : "In Order ") * (native ? "native_sin" : "sin"))

    ax = gca()

    n_lo = floor(Int, log2(xmin) + 0.25)
    n_hi = ceil(Int, log2(xmax) - 0.25)
    xticks(2.0.^(n_lo:2:n_hi), size_to_str.(2.0.^(n_lo:2:n_hi)))
    ax.set_xticks(2.0.^(n_lo:n_hi), minor=true)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(axis="y", right=true, which="minor", direction="in")

    t_lo = floor(Int, log10(ymin) + 0.25) - 9
    t_hi = ceil(Int, log10(ymax) - 0.25) - 9
    yticks(10.0.^(t_lo:t_hi), num_to_si.(10.0.^(t_lo:t_hi), Ref("s")))
end

for dev in devices
    figure(figsize=[12.6, 16.8])

    d = data[dev]

    ax = subplot(3, 2, 1)
    plot_dummy(d, false)

    ax = subplot(3, 2, 2)
    plot_dummy(d, true)

    ax = subplot(3, 2, 3)
    plot_compute(d, false, false)

    ax = subplot(3, 2, 4)
    plot_compute(d, true, false)

    ax = subplot(3, 2, 5)
    plot_compute(d, false, true)

    ax = subplot(3, 2, 6)
    plot_compute(d, true, true)

    tight_layout(pad=0.6)
    NaCsPlot.maybe_save("$(prefix)_$(dev)")
end

NaCsPlot.maybe_show()
