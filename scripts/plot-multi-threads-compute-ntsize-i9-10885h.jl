#!/usr/bin/julia

include("plot-multi-threads-compute.jl")
include("props_i9-10885h.jl")

const datadir = joinpath(@__DIR__, "../data/cpu-multi-threads-compute-ntsize/i9-10885h")

const data = load_all_times(datadir, nominal_ghz)

const prefix = joinpath(@__DIR__, "../imgs/multi-threads-compute-ntsize-i9-10885h")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
for (i, size) in enumerate([262144, 1048576, 4194304, 16777216, 67108864, 268435456])
    plot_lines(data["8chn-$size"], 16384 / 625e3, 16384, "C$i",
               size_to_str(size, "B"), alpha=0.6)
end
axhline(0.625, color="r")
text(3300, 0.625, "625 MS/s", va="bottom", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="g", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Output Length (ms)")
ylabel("Minimum Speed (S/ns)")
title("8 Channels")

ax = subplot(2, 2, 2)
for (i, size) in enumerate([262144, 1048576, 4194304, 16777216, 67108864, 268435456])
    plot_lines(data["16chn-$size"], 16384 / 625e3, 16384, "C$i",
               size_to_str(size, "B"), alpha=0.6)
end
axhline(0.625, color="r")
text(3300, 0.625, "625 MS/s", va="bottom", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="g", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Output Length (ms)")
ylabel("Minimum Speed (S/ns)")
title("16 Channels")

ax = subplot(2, 2, 3)
for (i, size) in enumerate([262144, 1048576, 4194304, 16777216, 67108864, 268435456])
    plot_lines(data["32chn-$size"], 16384 / 625e3, 16384, "C$i",
               size_to_str(size, "B"), alpha=0.6)
end
axhline(0.625, color="r")
text(0.16, 0.625, "625 MS/s", va="bottom", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="g", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Output Length (ms)")
ylabel("Minimum Speed (S/ns)")
title("32 Channels")

ax = subplot(2, 2, 4)
for (i, size) in enumerate([262144, 1048576, 4194304, 16777216, 67108864, 268435456])
    plot_lines(data["64chn-$size"], 16384 / 625e3, 16384, "C$i",
               size_to_str(size, "B"), alpha=0.6)
end
axhline(0.625, color="r")
text(330, 0.615, "625 MS/s", va="top", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="g", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Output Length (ms)")
ylabel("Minimum Speed (S/ns)")
title("64 Channels")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_minspeed")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
for (i, size) in enumerate([262144, 1048576, 4194304, 16777216, 67108864, 268435456])
    plot_runaheads(data["8chn-$size"], 16384, 16384 / 625e3, "C$i",
                   size_to_str(size, "B"), alpha=0.6)
end
axvline(0.625, color="r")
text(0.625, 0.15, "625 MS/s", va="center", ha="right", rotation=90, color="r", fontsize=19)
grid()
yscale("log")
ax.set_xlim(left=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Sustained Throughput (S/ns)")
ylabel("Required Delay (ms)")
title("8 Channels")

ax = subplot(2, 2, 2)
for (i, size) in enumerate([262144, 1048576, 4194304, 16777216, 67108864, 268435456])
    plot_runaheads(data["16chn-$size"], 16384, 16384 / 625e3, "C$i",
                   size_to_str(size, "B"), alpha=0.6)
end
axvline(0.625, color="r")
text(0.625, 0.15, "625 MS/s", va="center", ha="right", rotation=90, color="r", fontsize=19)
grid()
yscale("log")
ax.set_xlim(left=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Sustained Throughput (S/ns)")
ylabel("Required Delay (ms)")
title("16 Channels")

ax = subplot(2, 2, 3)
for (i, size) in enumerate([262144, 1048576, 4194304, 16777216, 67108864, 268435456])
    plot_runaheads(data["32chn-$size"], 16384, 16384 / 625e3, "C$i",
                   size_to_str(size, "B"), alpha=0.6)
end
axvline(0.625, color="r")
text(0.625, 0.1, "625 MS/s", va="center", ha="right", rotation=90, color="r", fontsize=19)
grid()
yscale("log")
ax.set_xlim(left=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Sustained Throughput (S/ns)")
ylabel("Required Delay (ms)")
title("32 Channels")

ax = subplot(2, 2, 4)
for (i, size) in enumerate([262144, 1048576, 4194304, 16777216, 67108864, 268435456])
    plot_runaheads(data["64chn-$size"], 16384, 16384 / 625e3, "C$i",
                   size_to_str(size, "B"), alpha=0.6)
end
axvline(0.625, color="r")
text(0.625, 0.33, "625 MS/s", va="center", ha="right", rotation=90, color="r", fontsize=19)
grid()
yscale("log")
ax.set_xlim(left=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Sustained Throughput (S/ns)")
ylabel("Required Delay (ms)")
title("64 Channels")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_delay")

NaCsPlot.maybe_show()
