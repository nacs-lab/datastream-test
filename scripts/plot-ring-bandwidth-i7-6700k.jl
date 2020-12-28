#!/usr/bin/julia

include("plot-ring-bandwidth-csv.jl")
include("props_i7-6700k.jl")

const datadir = joinpath(@__DIR__, "../data/")

data_block = load_data(joinpath(datadir, "ring-bandwidth-i7-6700k_blockring_read.csv"),
                       joinpath(datadir, "ring-bandwidth-i7-6700k_blockring_write.csv"))
data_block_empty = load_data(joinpath(datadir,
                                      "ring-bandwidth-empty-i7-6700k_blockring_read.csv"),
                             joinpath(datadir,
                                      "ring-bandwidth-empty-i7-6700k_blockring_write.csv"))

data_pipe = load_data(joinpath(datadir, "ring-bandwidth-i7-6700k_datapipe_read.csv"),
                      joinpath(datadir, "ring-bandwidth-i7-6700k_datapipe_write.csv"))
data_pipe_empty = load_data(joinpath(datadir,
                                     "ring-bandwidth-empty-i7-6700k_datapipe_read.csv"),
                            joinpath(datadir,
                                     "ring-bandwidth-empty-i7-6700k_datapipe_write.csv"))

const prefix = joinpath(@__DIR__, "../imgs/ring-bandwidth-i7-6700k")

const nblks = [4, 8, 16, 32, 64]

figure(figsize=[12.6, 5.6])

ax = subplot(1, 2, 1)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_block.read)
    line_wr = filter_data(x->x.nblk == nblk, data_block.write)
    @assert line_rd.size == line_wr.size
    plot(line_rd.size, (1024 ./ line_rd.byte_ns .+ 1024 ./ line_wr.byte_ns) ./ 2,
         label="$nblk blocks", color="C$(i - 1)", ls="-")
end
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_block_empty.read)
    line_wr = filter_data(x->x.nblk == nblk, data_block_empty.write)
    @assert line_rd.size == line_wr.size
    plot(line_rd.size, (1024 ./ line_rd.byte_ns .+ 1024 ./ line_wr.byte_ns) ./ 2,
         color="C$(i - 1)", ls="--")
end
grid()
axhline(1024 / 36.8, color="r", linewidth=2, ls="dotted")
text(1.6 * 2.0^26, 1024 / 36.8, "36.8 B/ns", va="bottom", ha="center", color="r", fontsize=19)
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 75, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
text(L3_total / ncores_phys, 75, "L3/core", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
ylim([0, 145])
xlabel("Buffer Size")
ylabel("Time (ns/KiB)")
title("BlockRing")
xticks(2.0.^(15:4:27), size_to_str.(2.0.^(15:4:27)))
ax.set_xticks(2.0.^(14:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(1, 2, 2)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_pipe.read)
    line_wr = filter_data(x->x.nblk == nblk, data_pipe.write)
    @assert line_rd.size == line_wr.size
    plot(line_rd.size, (1024 ./ line_rd.byte_ns .+ 1024 ./ line_wr.byte_ns) ./ 2,
         label="$nblk blocks", color="C$(i - 1)", ls="-")
end
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_pipe_empty.read)
    line_wr = filter_data(x->x.nblk == nblk, data_pipe_empty.write)
    @assert line_rd.size == line_wr.size
    plot(line_rd.size, (1024 ./ line_rd.byte_ns .+ 1024 ./ line_wr.byte_ns) ./ 2,
         color="C$(i - 1)", ls="--")
end
grid()
axhline(1024 / 36.8, color="r", linewidth=2, ls="dotted")
text(1.6 * 2.0^26, 1024 / 36.8, "36.8 B/ns", va="bottom", ha="center", color="r", fontsize=19)
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 75, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
text(L3_total / ncores_phys, 75, "L3/core", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
ylim([0, 135])
xlabel("Buffer Size")
ylabel("Time (ns/KiB)")
title("DataPipe")
xticks(2.0.^(15:4:27), size_to_str.(2.0.^(15:4:27)))
ax.set_xticks(2.0.^(14:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_perf")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_block.read)
    line_wr = filter_data(x->x.nblk == nblk, data_block.write)
    plot(line_rd.size, line_rd.miss_perc, label="$nblk blocks", color="C$(i - 1)", ls="-")
    plot(line_wr.size, line_wr.miss_perc, color="C$(i - 1)", ls="--")
end
grid()
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 55, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
text(L3_total / ncores_phys, 50, "L3/core", rotation=90, va="center", ha="left", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
xlim([1.4 * 2^20, 1.4 * 2^28])
ylim([0, 100])
xlabel("Buffer Size")
ylabel("Cache Miss (%)")
xticks(2.0.^(21:2:27), size_to_str.(2.0.^(21:2:27)))
ax.set_xticks(2.0.^(21:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 2)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_block.read)
    line_wr = filter_data(x->x.nblk == nblk, data_block.write)
    plot(line_rd.size, line_rd.pipe_rw, label="$nblk blocks", color="C$(i - 1)", ls="-")
    plot(line_wr.size, line_wr.pipe_rw, color="C$(i - 1)", ls="--")
end
grid()
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
ylim([0.9, 1.1])
xlabel("Buffer Size")
ylabel("R/W per Block")
xticks(2.0.^(15:4:27), size_to_str.(2.0.^(15:4:27)))
ax.set_xticks(2.0.^(14:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 3)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_block.read)
    line_wr = filter_data(x->x.nblk == nblk, data_block.write)
    plot(line_rd.size, line_rd.pipe_stall_perc, label="$nblk blocks", color="C$(i - 1)", ls="-")
    plot(line_wr.size, line_wr.pipe_stall_perc, color="C$(i - 1)", ls="--")
end
grid()
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 50, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
text(L3_total / ncores_phys, 50, "L3/core", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
ylim([0, 100])
xlabel("Buffer Size")
ylabel("Pipe Stall (%)")
xticks(2.0.^(15:4:27), size_to_str.(2.0.^(15:4:27)))
ax.set_xticks(2.0.^(14:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 4)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_block.read)
    line_wr = filter_data(x->x.nblk == nblk, data_block.write)
    plot(line_rd.size, line_rd.pipe_sync, label="$nblk blocks", color="C$(i - 1)", ls="-")
    plot(line_wr.size, line_wr.pipe_sync, color="C$(i - 1)", ls="--")
end
grid()
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 2.5, "L3", rotation=90, va="center", ha="left", color="C2", fontsize=19)
axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
text(L3_total / ncores_phys, 34, "L3/core", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
ylim([0, 46])
xlabel("Buffer Size")
ylabel("Sync per Stall")
xticks(2.0.^(15:4:27), size_to_str.(2.0.^(15:4:27)))
ax.set_xticks(2.0.^(14:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_blockring")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_pipe.read)
    line_wr = filter_data(x->x.nblk == nblk, data_pipe.write)
    plot(line_rd.size, line_rd.miss_perc, label="$nblk blocks", color="C$(i - 1)", ls="-")
    plot(line_wr.size, line_wr.miss_perc, color="C$(i - 1)", ls="--")
end
grid()
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 55, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
text(L3_total / ncores_phys, 50, "L3/core", rotation=90, va="center", ha="left", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
xlim([1.4 * 2^20, 1.4 * 2^28])
ylim([0, 100])
xlabel("Buffer Size")
ylabel("Cache Miss (%)")
xticks(2.0.^(21:2:27), size_to_str.(2.0.^(21:2:27)))
ax.set_xticks(2.0.^(21:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 2)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_pipe.read)
    line_wr = filter_data(x->x.nblk == nblk, data_pipe.write)
    plot(line_rd.size, line_rd.pipe_rw, label="$nblk blocks", color="C$(i - 1)", ls="-")
    plot(line_wr.size, line_wr.pipe_rw, color="C$(i - 1)", ls="--")
end
grid()
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
ylim([1, 1.54])
xlabel("Buffer Size")
ylabel("R/W per Block")
xticks(2.0.^(15:4:27), size_to_str.(2.0.^(15:4:27)))
ax.set_xticks(2.0.^(14:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 3)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_pipe.read)
    line_wr = filter_data(x->x.nblk == nblk, data_pipe.write)
    plot(line_rd.size, line_rd.pipe_stall_perc, label="$nblk blocks", color="C$(i - 1)", ls="-")
    plot(line_wr.size, line_wr.pipe_stall_perc, color="C$(i - 1)", ls="--")
end
grid()
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 50, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
text(L3_total / ncores_phys, 50, "L3/core", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
ylim([0, 100])
xlabel("Buffer Size")
ylabel("Pipe Stall (%)")
xticks(2.0.^(15:4:27), size_to_str.(2.0.^(15:4:27)))
ax.set_xticks(2.0.^(14:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 4)
for i in 1:length(nblks)
    nblk = nblks[i]
    line_rd = filter_data(x->x.nblk == nblk, data_pipe.read)
    line_wr = filter_data(x->x.nblk == nblk, data_pipe.write)
    plot(line_rd.size, line_rd.pipe_sync, label="$nblk blocks", color="C$(i - 1)", ls="-")
    plot(line_wr.size, line_wr.pipe_sync, color="C$(i - 1)", ls="--")
end
grid()
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 25, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
text(L3_total / ncores_phys, 23, "L3/core", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
ylim([0, 46])
xlabel("Buffer Size")
ylabel("Sync per Stall")
xticks(2.0.^(15:4:27), size_to_str.(2.0.^(15:4:27)))
ax.set_xticks(2.0.^(14:28), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_datapipe")

NaCsPlot.maybe_show()