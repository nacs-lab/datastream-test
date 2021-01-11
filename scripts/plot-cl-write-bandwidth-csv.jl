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

# function plot_perf(data, nblks, localbuffs, nchn_rd, nchn_wr)
#     # Note that this is very specific to how we took the data and isn't very generic...
#     @assert length(localbuffs) == 2
#     nchn = max(nchn_rd, nchn_wr)
#     for i in 1:length(nblks)
#         nblk = nblks[i]
#         for localbuff in localbuffs
#             line_rd = filter_data(x->(x.nblk == nblk && x.localbuff == localbuff &&
#                                       x.nchn_rd == nchn_rd && x.nchn_wr == nchn_wr), data.read)
#             line_wr = filter_data(x->(x.nblk == nblk && x.localbuff == localbuff &&
#                                       x.nchn_rd == nchn_rd && x.nchn_wr == nchn_wr), data.write)
#             @assert line_rd.size == line_wr.size
#             isempty(line_rd.size) && continue
#             # Make the color consistent with the previous measurements
#             if localbuff == 0
#                 plot(line_rd.size, (line_rd.byte_ns .+ line_wr.byte_ns) ./ 2 ./ 4 .* nchn,
#                      label="$nblk blocks", color="C$(i)", ls="-")
#             else
#                 plot(line_rd.size, (line_rd.byte_ns .+ line_wr.byte_ns) ./ 2 ./ 4 .* nchn,
#                      color="C$(i)", ls="--")
#             end
#         end
#     end
#     grid()
#     legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
#     title("Write $(nchn_wr)CH + Read $(nchn_rd)CH")
#     xscale("log")
#     xlabel("Buffer Size")
#     ylabel("Throughput (S/ns)")
#     ax = gca()
#     xticks(2.0.^(18:2:22), size_to_str.(2.0.^(18:2:22)))
#     ax.set_xticks(2.0.^(17:22), minor=true)
#     ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# end

# function plot_stall(data, nblks, localbuffs, nchn_rd, nchn_wr)
#     # Note that this is very specific to how we took the data and isn't very generic...
#     @assert length(localbuffs) == 2
#     nchn = max(nchn_rd, nchn_wr)
#     for i in 1:length(nblks)
#         nblk = nblks[i]
#         for localbuff in localbuffs
#             line_rd = filter_data(x->(x.nblk == nblk && x.localbuff == localbuff &&
#                                       x.nchn_rd == nchn_rd && x.nchn_wr == nchn_wr), data.read)
#             line_wr = filter_data(x->(x.nblk == nblk && x.localbuff == localbuff &&
#                                       x.nchn_rd == nchn_rd && x.nchn_wr == nchn_wr), data.write)
#             @assert line_rd.size == line_wr.size
#             isempty(line_rd.size) && continue
#             # Make the color consistent with the previous measurements
#             if localbuff == 0
#                 plot(line_rd.size, line_rd.pipe_stall_perc,
#                      label="$nblk blocks", color="C$(i)", ls="-")
#                 plot(line_wr.size, line_wr.pipe_stall_perc,
#                      color="C$(i)", ls=":")
#             else
#                 plot(line_rd.size, line_rd.pipe_stall_perc,
#                      color="C$(i)", ls="--")
#                 plot(line_wr.size, line_wr.pipe_stall_perc,
#                      color="C$(i)", ls="-.")
#             end
#         end
#     end
#     grid()
#     legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
#     title("Write $(nchn_wr)CH + Read $(nchn_rd)CH")
#     xscale("log")
#     xlabel("Buffer Size")
#     ylabel("Pipe Stall (%)")
#     ylim([0, 100])
#     ax = gca()
#     xticks(2.0.^(18:2:22), size_to_str.(2.0.^(18:2:22)))
#     ax.set_xticks(2.0.^(17:22), minor=true)
#     ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# end
