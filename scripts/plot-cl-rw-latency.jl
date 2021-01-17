#!/usr/bin/julia

using NaCsCalc
using NaCsPlot
using PyPlot
using LibArchive

include("utils.jl")

const prefix = joinpath(@__DIR__, "../imgs/cl-rw-latency")

function decode_sign(v::UInt64)
    if v & 1 == 0
        return Int(v >> 1)
    else
        return -Int(v >> 1) - 1
    end
end

struct TimeData
    nele::Int
    nconcurrent::Int
    queue_times::Matrix{Float64}
    submit_times::Matrix{Float64}
    start_times::Matrix{Float64}
    complete_times::Matrix{Float64}
end

read_compressed_res(fname, freq) = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    data = Int[v for v in LEB128(read(reader))]
    nele = data[1]
    nconcurrent = data[2]
    nrep = data[3]
    res = TimeData(nele, nconcurrent, Matrix{Float64}(undef, nrep, nconcurrent),
                   Matrix{Float64}(undef, nrep, nconcurrent),
                   Matrix{Float64}(undef, nrep, nconcurrent),
                   Matrix{Float64}(undef, nrep, nconcurrent))
    idx = 4
    for times in (res.queue_times, res.submit_times, res.start_times, res.complete_times)
        for ic in 1:nconcurrent
            t = 0
            for i in 1:nrep
                t = t + data[idx]
                times[i, ic] = t / freq
                idx += 1
            end
        end
    end
    return res
end

const datadir = joinpath(@__DIR__, "../data/cl-rw-latency")

const devices = ["amd-gpu-gfx1012","intel-gpu-1912", "intel-gpu-9bc4"]
const freqs = Dict("amd-gpu-gfx1012"=>4.0,"intel-gpu-1912"=>4.0, "intel-gpu-9bc4"=>2.4)

const max_bandwidths = Dict("amd-gpu-gfx1012"=>224e9,
                            "intel-gpu-1912"=>34.1e9,
                            "intel-gpu-9bc4"=>45.8e9)

const data = Dict(dev=>TimeData[] for dev in devices)

for file in readdir(datadir)
    for dev in devices
        if startswith(file, dev)
            res = read_compressed_res(joinpath(datadir, file), freqs[dev])
            push!(data[dev], res)
        end
    end
end
for res in values(data)
    sort!(res, by=x->(x.nele, x.nconcurrent))
end

function average_diff(d1, d2, offset=0)
    total = 0.0
    @assert size(d1) == size(d2)
    trange = 128:(size(d1, 1) - 64)
    for j in 1:size(d1, 2)
        for i in trange
            total += d1[i + offset, j] - d2[i, j]
        end
    end
    return total / length(trange) / size(d1, 2)
end

function compute_summary(data::Vector{TimeData})
    res = (nele=Int[], nconcurrent=Int[], submit_delay=Float64[], start_delay=Float64[],
           complete_delay=Float64[], queue_rate=Float64[], complete_rate=Float64[],
           next_queue=Float64[], next_submit=Float64[], next_start=Float64[])
    for td in data
        push!(res.nele, td.nele)
        push!(res.nconcurrent, td.nconcurrent)

        push!(res.submit_delay, average_diff(td.submit_times, td.queue_times))
        push!(res.start_delay, average_diff(td.start_times, td.submit_times))
        push!(res.complete_delay, average_diff(td.complete_times, td.start_times))

        push!(res.queue_rate, average_diff(td.queue_times, td.queue_times, 1))
        push!(res.complete_rate, average_diff(td.complete_times, td.complete_times, 1))

        push!(res.next_queue, average_diff(td.queue_times, td.complete_times, 1))
        push!(res.next_submit, average_diff(td.submit_times, td.complete_times, 1))
        push!(res.next_start, average_diff(td.start_times, td.complete_times, 1))
    end
    return res
end

function filter_data(cb, data)
    res = (nele=Int[], nconcurrent=Int[], submit_delay=Float64[], start_delay=Float64[],
           complete_delay=Float64[], queue_rate=Float64[], complete_rate=Float64[],
           next_queue=Float64[], next_submit=Float64[], next_start=Float64[])
    for i in 1:length(data.nele)
        point = (nele=data.nele[i], nconcurrent=data.nconcurrent[i],
                 submit_delay=data.submit_delay[i], start_delay=data.start_delay[i],
                 complete_delay=data.complete_delay[i],
                 queue_rate=data.queue_rate[i], complete_rate=data.complete_rate[i],
                 next_queue=data.next_queue[i], next_submit=data.next_submit[i],
                 next_start=data.next_start[i])
        cb(point) || continue
        push!(res.nele, point.nele)
        push!(res.nconcurrent, point.nconcurrent)

        push!(res.submit_delay, point.submit_delay)
        push!(res.start_delay, point.start_delay)
        push!(res.complete_delay, point.complete_delay)

        push!(res.queue_rate, point.queue_rate)
        push!(res.complete_rate, point.complete_rate)

        push!(res.next_queue, point.next_queue)
        push!(res.next_submit, point.next_submit)
        push!(res.next_start, point.next_start)
    end
    return res
end

sort_u(vals) = sort(unique(vals))

function plot_nele(data, field)
    nconcurrents = sort_u(data.nconcurrent)
    xmin = Inf
    xmax = 0.0
    ymin = Inf
    ymax = 0.0
    for nconcurrent in nconcurrents
        line = filter_data(x->x.nconcurrent == nconcurrent, data)
        x = line.nele .* 4
        y = getfield(line, field)::Vector{Float64} .* 1e-9

        xmin = min(xmin, minimum(x))
        xmax = max(xmax, maximum(x))
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
        plot(x, y, ".-", label="$(nconcurrent)")
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
        ylabel("Time")
        y_lo = floor(Int, log10(ymin) + 0.25)
        y_hi = ceil(Int, log10(ymax) - 0.25)
        yticks(10.0.^(y_lo:y_hi), num_to_si.(10.0.^(y_lo:y_hi), Ref("s")))
        # ax.set_yticks(10.0.^(y_lo:y_hi), minor=true)
        ax.tick_params(axis="y", right=true, which="minor", direction="in")
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else
        ylabel("Time (s)")
    end
end

function plot_maxbw(d, max_bw)
    xlim(xlim()) # Lock in X limit
    ylim(ylim()) # Lock in Y limit
    cidx = 0
    for nconcurrent in sort_u(d.nconcurrent)
        plotx = linspace(xlim()..., 1000)
        ploty = plotx ./ (max_bw / nconcurrent)
        plot(plotx, ploty, color="C$cidx", linestyle="dotted", alpha=0.8)
        cidx += 1
    end
end

for dev in devices
    d = compute_summary(data[dev])
    max_bw = max_bandwidths[dev] / 2

    figure(figsize=[6.3 * 2, 5.6 * 4])
    subplot(4, 2, 1)
    plot_nele(d, :queue_rate)
    title("Queue Rate")

    subplot(4, 2, 2)
    plot_nele(d, :complete_rate)
    title("Complete Rate")
    plot_maxbw(d, max_bw)

    subplot(4, 2, 3)
    plot_nele(d, :submit_delay)
    title("Queue -> Submit")

    subplot(4, 2, 4)
    plot_nele(d, :start_delay)
    title("Submit -> Start")

    subplot(4, 2, 5)
    plot_nele(d, :complete_delay)
    title("Start -> Complete")

    subplot(4, 2, 6)
    plot_nele(d, :next_queue)
    title("Complete -> Next Queue")

    subplot(4, 2, 7)
    plot_nele(d, :next_submit)
    title("Complete -> Next Submit")

    subplot(4, 2, 8)
    plot_nele(d, :next_start)
    title("Complete -> Next Start")

    tight_layout(pad=0.6)
    NaCsPlot.maybe_save("$(prefix)_$(dev)")
end

NaCsPlot.maybe_show()
