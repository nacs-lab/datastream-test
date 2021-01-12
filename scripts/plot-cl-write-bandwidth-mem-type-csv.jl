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
    nalloc = Int.(@view data[:, 3])
    readable = @view(data[:, 4]) .!= 0
    host_access = @view(data[:, 5]) .!= 0
    host_write = @view(data[:, 6]) .!= 0
    host_ptr = @view(data[:, 7]) .!= 0
    t = data[:, 8]
    return (nele=nele, nrep=nrep, nalloc=nalloc, readable=readable,
            host_access=host_access, host_write=host_write, host_ptr=host_ptr, t=t)
end

load_data(name) = preprocess_data(readdlm(name, ',', skipstart=1))

function filter_data(cb, data)
    nele = Int[]
    nrep = Int[]
    nalloc = Int[]
    readable = Bool[]
    host_access = Bool[]
    host_write = Bool[]
    host_ptr = Bool[]
    t = Float64[]
    for i in 1:length(data.nele)
        point = (nele=data.nele[i], nrep=data.nrep[i], nalloc=data.nalloc[i],
                 readable=data.readable[i], host_access=data.host_access[i],
                 host_write=data.host_write[i], host_ptr=data.host_ptr[i], t=data.t[i])
        cb(point) || continue
        push!(nele, point.nele)
        push!(nrep, point.nrep)
        push!(nalloc, point.nalloc)
        push!(readable, point.readable)
        push!(host_access, point.host_access)
        push!(host_write, point.host_write)
        push!(host_ptr, point.host_ptr)
        push!(t, point.t)
    end
    return (nele=nele, nrep=nrep, nalloc=nalloc, readable=readable,
            host_access=host_access, host_write=host_write, host_ptr=host_ptr, t=t)
end

sort_u(vals) = sort(unique(vals))

const prefix = joinpath(@__DIR__, "../imgs/cl-write-bandwidth-mem-type")
const datadir = joinpath(@__DIR__, "../data/cl-write-bandwidth-mem-type")

const devices = ["amd-gpu-gfx1012", "intel-gpu-1912", "intel-gpu-9bc4"]
const data = Dict(dev=>load_data(joinpath(datadir, "$(dev).csv"))
                  for dev in devices)
const max_bandwidths = Dict("amd-gpu-gfx1012"=>224e9,
                            "intel-gpu-1912"=>34.1e9,
                            "intel-gpu-9bc4"=>45.8e9)

function data_filter(min_nalloc, readable, host_access, host_write, host_ptr)
    return function (x, nele)
        if x.nele != nele
            return false
        end
        nalloc = max(min_nalloc, nele)
        if x.nalloc != nalloc
            return false
        end
        return (x.readable == readable && x.host_access == host_access &&
                x.host_write == host_write && x.host_ptr == host_ptr)
    end
end

function plot_nrep(data, max_bw, filter)
    neles = 2 .^ (12:23)
    cnt_mid = 2
    cnt_big = 0
    xmin = Inf
    xmax = 0.0
    ymin = max_bw
    ymax = max_bw
    for nele in neles
        line = filter_data(x->filter(x, nele), data)
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
        y = 4e9 ./ line.t[perm]

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

function plot_minsizes(d, max_bw, sp_offset, readable, host_access, host_write, host_ptr)
    if host_ptr
        suffix = "Host Ptr"
    elseif host_write
        suffix = "Host RW"
    elseif host_access
        suffix = "Host RO"
    else
        suffix = "Host NX"
    end
    if readable
        suffix = "Device RW/$suffix"
    else
        suffix = "Device WO/$suffix"
    end
    ax = subplot(8, 3, 1 + sp_offset)
    plot_nrep(d, max_bw, data_filter(0, readable, host_access, host_write, host_ptr))
    title("Min Buff: 0 B/$suffix", fontsize=20)
    ax = subplot(8, 3, 2 + sp_offset)
    plot_nrep(d, max_bw, data_filter(2^20, readable, host_access, host_write, host_ptr))
    title("Min Buff: $(size_to_str(2^20, "B"))/$suffix", fontsize=20)
    ax = subplot(8, 3, 3 + sp_offset)
    plot_nrep(d, max_bw, data_filter(2^23, readable, host_access, host_write, host_ptr))
    title("Min Buff: $(size_to_str(2^23, "B"))/$suffix", fontsize=20)
end

function plot_readable(d, max_bw, sp_offset, host_access, host_write, host_ptr)
    plot_minsizes(d, max_bw, sp_offset, false, host_access, host_write, host_ptr)
    plot_minsizes(d, max_bw, sp_offset + 3, true, host_access, host_write, host_ptr)
end

for dev in devices
    figure(figsize=[6.3 * 3, 5.6 * 8])

    d = data[dev]
    max_bw = max_bandwidths[dev]

    plot_readable(d, max_bw, 0, false, false, false)
    plot_readable(d, max_bw, 6, true, false, false)
    plot_readable(d, max_bw, 12, true, true, false)
    plot_readable(d, max_bw, 18, true, true, true)

    tight_layout(pad=0.6)
    NaCsPlot.maybe_save("$(prefix)_$(dev)")
end

NaCsPlot.maybe_show()
