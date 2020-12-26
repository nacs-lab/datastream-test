#!/usr/bin/julia

include("plot-bandwidth-csv.jl")

const datadir = joinpath(@__DIR__, "../data/")

const L1d = 32 * 1024
const L2 = 256 * 1024
const L3_total = 8 * 1024^2

data_normal = load_data(joinpath(datadir, "write-bandwidth-i7-6700k.csv"));
data_nt = load_data(joinpath(datadir, "write-bandwidth-nt-i7-6700k.csv"));
data_no_prefetch = load_data(joinpath(datadir, "write-bandwidth-no-prefetch-i7-6700k.csv"));

const prefix = joinpath(@__DIR__, "../imgs/bandwidth-i7-6700k")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
for line in data_normal
    line = filter_data(x->x.size <= L2, line)
    plot(line.size, line.byte_ns, label="$(line.core) cores")
end
axvline(L1d, color="C3", linewidth=3)
text(L1d, 100, "L1d", rotation=90, va="center", ha="right", color="C3", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
ylim([58, 145])
grid()
title("Normal")
xlabel("Buffer Size")
ylabel("Throughput (B/ns)")
xticks(2.0.^(14:2:16), size_to_str.(2.0.^(14:2:16)))
ax.set_xticks(2.0.^(13:17), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 2)
for line in data_normal
    line = filter_data(x->x.size <= L2, line)
    plot(line.size, line.byte_cyl, label="$(line.core) cores")
end
axvline(L1d, color="C3", linewidth=3)
text(L1d, 21, "L1d", rotation=90, va="center", ha="right", color="C3", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
ylim([12.5, 32])
grid()
title("Normal")
xlabel("Buffer Size")
ylabel("Throughput (B/cycle)")
xticks(2.0.^(14:2:16), size_to_str.(2.0.^(14:2:16)))
ax.set_xticks(2.0.^(13:17), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 3)
for line in data_no_prefetch
    line = filter_data(x->x.size <= L2, line)
    plot(line.size, line.byte_ns, label="$(line.core) cores")
end
axvline(L1d, color="C3", linewidth=3)
text(L1d, 100, "L1d", rotation=90, va="center", ha="right", color="C3", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
ylim([58, 145])
grid()
title("No Prefetcher")
xlabel("Buffer Size")
ylabel("Throughput (B/ns)")
xticks(2.0.^(14:2:16), size_to_str.(2.0.^(14:2:16)))
ax.set_xticks(2.0.^(13:17), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 4)
for line in data_no_prefetch
    line = filter_data(x->x.size <= L2, line)
    plot(line.size, line.byte_cyl, label="$(line.core) cores")
end
axvline(L1d, color="C3", linewidth=3)
text(L1d, 21, "L1d", rotation=90, va="center", ha="right", color="C3", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
ylim([12.5, 32])
grid()
title("No Prefetcher")
xlabel("Buffer Size")
ylabel("Throughput (B/cycle)")
xticks(2.0.^(14:2:16), size_to_str.(2.0.^(14:2:16)))
ax.set_xticks(2.0.^(13:17), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_small")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
for line in data_normal
    line = filter_data(x->x.size >= L2, line)
    plot(line.size .* line.core, line.byte_ns .* line.core, label="$(line.core) cores")
end
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 18, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
yscale("log")
ylim([11, 300])
grid()
title("Normal")
xlabel("Total Buffer Size")
ylabel("Total Throughput (B/ns)")
xticks(2.0.^(20:4:28), size_to_str.(2.0.^(20:4:28)))
ax.set_xticks(2.0.^(18:29), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
yticks([30, 100, 300], ["30", "100", "300"])
ax.set_yticks([20:10:100; 200], minor=true)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 2)
for line in data_normal
    line = filter_data(x->x.size >= L2, line)
    plot(line.size .* line.core, line.byte_cyl .* line.core, label="$(line.core) cores")
end
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 4.5, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
yscale("log")
ylim([2.8, 65])
grid()
title("Normal")
xlabel("Total Buffer Size")
ylabel("Total Throughput (B/cycle)")
xticks(2.0.^(20:4:28), size_to_str.(2.0.^(20:4:28)))
ax.set_xticks(2.0.^(18:29), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
yticks([3, 10, 30], ["3", "10", "30"])
ax.set_yticks([3:10; 20:10:60], minor=true)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 3)
for line in data_no_prefetch
    line = filter_data(x->x.size >= L2, line)
    plot(line.size .* line.core, line.byte_ns .* line.core, label="$(line.core) cores")
end
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 18, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
yscale("log")
ylim([11, 300])
grid()
title("No Prefetcher")
xlabel("Total Buffer Size")
ylabel("Total Throughput (B/ns)")
xticks(2.0.^(20:4:28), size_to_str.(2.0.^(20:4:28)))
ax.set_xticks(2.0.^(18:29), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
yticks([30, 100, 300], ["30", "100", "300"])
ax.set_yticks([20:10:100; 200], minor=true)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(2, 2, 4)
for line in data_no_prefetch
    line = filter_data(x->x.size >= L2, line)
    plot(line.size .* line.core, line.byte_cyl .* line.core, label="$(line.core) cores")
end
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 4.5, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=1, columnspacing=1, handlelength=1)
xscale("log")
yscale("log")
ylim([2.8, 65])
grid()
title("No Prefetcher")
xlabel("Total Buffer Size")
ylabel("Total Throughput (B/cycle)")
xticks(2.0.^(20:4:28), size_to_str.(2.0.^(20:4:28)))
ax.set_xticks(2.0.^(18:29), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
yticks([3, 10, 30], ["3", "10", "30"])
ax.set_yticks([3:10; 20:10:60], minor=true)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_big")

figure(figsize=[12.6, 5.6])

ax = subplot(1, 2, 1)
for line in data_normal
    line = filter_data(x->(x.byte_ref < 150), line)
    plot(line.size .* line.core, line.miss_perc, label="$(line.core) cores")
end
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 70, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
xlim([3 * 1024^2, 0.7 * 1024^3])
ylim([0, 100])
grid()
title("Normal")
xlabel("Total Buffer Size")
ylabel("Cache Miss (%)")
xticks(2.0.^(24:4:28), size_to_str.(2.0.^(24:4:28)))
ax.set_xticks(2.0.^(23:29), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(1, 2, 2)
for line in data_no_prefetch
    line = filter_data(x->(x.byte_ref < 150), line)
    plot(line.size .* line.core, line.miss_perc, label="$(line.core) cores")
end
axvline(L3_total, color="C2", linewidth=3)
text(L3_total, 70, "L3", rotation=90, va="center", ha="right", color="C2", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
xlim([3 * 1024^2, 0.7 * 1024^3])
ylim([0, 100])
grid()
title("No Prefetcher")
xlabel("Total Buffer Size")
ylabel("Cache Miss (%)")
xticks(2.0.^(24:4:28), size_to_str.(2.0.^(24:4:28)))
ax.set_xticks(2.0.^(23:29), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_cache_miss")

figure(figsize=[12.6, 5.6])

ax = subplot(1, 2, 1)
for line in data_nt
    plot(line.size, line.byte_ns .* line.core, label="$(line.core) cores")
end
axvline(L1d, color="C3", linewidth=3)
text(L1d, 33.5, "L1d", rotation=90, va="center", ha="right", color="C3", fontsize=19)
axvline(L2, color="C1", linewidth=3)
text(L2, 33.5, "L2", rotation=90, va="center", ha="right", color="C1", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
grid()
title("Non-Temporal Store")
xlabel("Buffer Size")
ylabel("Total Throughput (B/ns)")
xticks(2.0.^(12:4:24), size_to_str.(2.0.^(12:4:24)))
ax.set_xticks(2.0.^(12:27), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax = subplot(1, 2, 2)
for line in data_nt
    plot(line.size, line.byte_cyl .* line.core, label="$(line.core) cores")
end
axvline(L1d, color="C3", linewidth=3)
text(L1d, 7.4, "L1d", rotation=90, va="center", ha="right", color="C3", fontsize=19)
axvline(L2, color="C1", linewidth=3)
text(L2, 7.4, "L2", rotation=90, va="center", ha="right", color="C1", fontsize=19)
legend(fontsize="x-small", ncol=2, columnspacing=1, handlelength=1)
xscale("log")
grid()
title("Non-Temporal Store")
xlabel("Buffer Size")
ylabel("Total Throughput (B/cycle)")
xticks(2.0.^(12:4:24), size_to_str.(2.0.^(12:4:24)))
ax.set_xticks(2.0.^(12:27), minor=true)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_nt")

NaCsPlot.maybe_show()
