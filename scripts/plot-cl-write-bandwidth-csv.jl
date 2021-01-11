#!/usr/bin/julia

using NaCsCalc
using NaCsData
using NaCsData.Fitting: fit_data
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
const max_bandwidths = Dict("amd-gpu-gfx1012"=>224e9,
                            "intel-gpu-1912"=>34.1e9,
                            "intel-gpu-9bc4"=>45.8e9)

function plot_nrep(data, max_bw)
    neles = 2 .^ (10:23)
    cnt_mid = 0
    cnt_big = 0
    xmin = Inf
    xmax = 0.0
    ymin = max_bw
    ymax = max_bw
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
    axhline(max_bw, color="r", ls="dotted", linewidth=2)
    xscale("log")
    yscale("log")
    grid()
    xlabel("Repetition")
    ylabel("Memory Throughput (B/s)")
    legend(fontsize="xx-small", ncol=6, columnspacing=0.8, handlelength=1,
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

function plot_nwrite(data, max_bw)
    neles = 2 .^ (10:23)
    cnt_mid = 0
    cnt_big = 0
    # xmin = Inf
    # xmax = 0.0
    ymin = max_bw
    ymax = max_bw
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
    axhline(max_bw, color="r", ls="dotted", linewidth=2)
    yscale("log")
    grid()
    xlabel("Memory Write per Kernel (B)")
    ylabel("Memory Throughput (B/s)")
    legend(fontsize="xx-small", ncol=6, columnspacing=0.8, handlelength=1,
           loc="lower center")

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
    plot_nrep(d, max_bandwidths[dev])

    ax = subplot(1, 2, 2)
    plot_nwrite(d, max_bandwidths[dev])

    tight_layout(pad=0.6)
    NaCsPlot.maybe_save("$(prefix)_$(dev)")
end

function fit_nrep(data)
    neles = 2 .^ (19:23)
    cnt_big = 0
    xmin = Inf
    xmax = 0.0
    ymin = Inf
    ymax = 0.0
    data = filter_data(x->(any((nele - 2 <= x.nele <= nele + 2 for nele in neles))), data)

    function real_model(nrep, buf_sz, p)
        max_bw, init_time, run_time = p
        perrun_time = init_time / nrep + run_time + buf_sz / max_bw
        return buf_sz / perrun_time
    end

    scalar_model(i, p) = real_model(data.nrep[i], data.nele[i] * 4, p)
    model(i, p) = scalar_model.(i, Ref(p))

    fit = fit_data(model, 1:length(data.t), 4e9 ./ data.t, [200e9, 3e-3, 0]; plotx=false)
    @show fit.uncs

    for nele in neles
        line = filter_data(x->nele - 2 <= x.nele <= nele + 2, data)
        color = "C$(cnt_big)"
        cnt_big += 1
        perm = sortperm(line.nrep)
        x = line.nrep[perm]
        y = 4e9 ./ line.t[perm]
        xmin = min(xmin, x[1])
        xmax = max(xmax, x[end])
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
        plot(x, y, "o", color=color, label="$(size_to_str(nele * 4))")
        plotx = linspace(x[1] * 0.9, x[end] * 1.1, 10000)
        plot(plotx, real_model.(plotx, nele * 4, Ref(fit.param)), color=color)
    end
    xscale("log")
    yscale("log")
    grid()
    xlabel("Repetition")
    ylabel("Memory Throughput (B/s)")
    legend(fontsize="xx-small", ncol=6, columnspacing=0.8, handlelength=1,
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

figure()
fit_nrep(data["amd-gpu-gfx1012"])
tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_amd-gpu-gfx1012-large_buff")

NaCsPlot.maybe_show()
