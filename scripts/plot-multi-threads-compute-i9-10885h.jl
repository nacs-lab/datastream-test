#!/usr/bin/julia

include("plot-multi-threads-compute.jl")
include("props_i9-10885h.jl")

const datadir = joinpath(@__DIR__, "../data/cpu-multi-threads-compute/i9-10885h")

const data = load_all_times(datadir, nominal_ghz)

const prefix = joinpath(@__DIR__, "../imgs/multi-threads-compute-i9-10885h")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
plot_lines(data["16chn-nolocal"], 16384 / 625e3, 16384, "C0", "No Buffer", ls="--")
plot_lines(data["16chn-local"], 16384 / 625e3, 16384, "C1", "16 KiB Buffer")
axhline(0.625, color="r")
text(3300, 0.625, "625 MS/s", va="bottom", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="C2", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Output Length (ms)")
ylabel("Minimum Speed (S/ns)")
title("16 Channels")

ax = subplot(2, 2, 2)
plot_lines(data["40chn-nolocal"], 16384 / 625e3, 16384, "C0", "No Buffer", ls="--")
plot_lines(data["40chn-local"], 16384 / 625e3, 16384, "C1", "16 KiB Buffer")
axhline(0.625, color="r")
text(15, 0.625, "625 MS/s", va="bottom", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="C2", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Output Length (ms)")
ylabel("Minimum Speed (S/ns)")
title("40 Channels")

ax = subplot(2, 2, 3)
plot_lines(data["80chn-nolocal"], 16384 / 625e3, 16384, "C0", "No Buffer", ls="--")
plot_lines(data["80chn-local"], 16384 / 625e3, 16384, "C1", "16 KiB Buffer")
axhline(0.625, color="r")
text(330, 0.615, "625 MS/s", va="top", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="C2", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Output Length (ms)")
ylabel("Minimum Speed (S/ns)")
title("80 Channels")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_minspeed")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
plot_runaheads(data["16chn-nolocal"], 16384, 16384 / 625e3, "C0", "No Buffer", ls="--")
plot_runaheads(data["16chn-local"], 16384, 16384 / 625e3, "C1", "16 KiB Buffer")
axvline(0.625, color="r")
text(0.625, 0.15, "625 MS/s", va="center", ha="right", rotation=90, color="r", fontsize=19)
grid()
yscale("log")
ax.set_xlim(left=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Sustained Throughput (S/ns)")
ylabel("Required Delay (ms)")
title("16 Channels")

ax = subplot(2, 2, 2)
plot_runaheads(data["40chn-nolocal"], 16384, 16384 / 625e3, "C0", "No Buffer", ls="--")
plot_runaheads(data["40chn-local"], 16384, 16384 / 625e3, "C1", "16 KiB Buffer")
axvline(0.625, color="r")
text(0.625, 0.1, "625 MS/s", va="center", ha="right", rotation=90, color="r", fontsize=19)
grid()
yscale("log")
ax.set_xlim(left=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Sustained Throughput (S/ns)")
ylabel("Required Delay (ms)")
title("40 Channels")

ax = subplot(2, 2, 3)
plot_runaheads(data["80chn-nolocal"], 16384, 16384 / 625e3, "C0", "No Buffer", ls="--")
plot_runaheads(data["80chn-local"], 16384, 16384 / 625e3, "C1", "16 KiB Buffer")
axvline(0.625, color="r")
text(0.625, 0.33, "625 MS/s", va="center", ha="right", rotation=90, color="r", fontsize=19)
grid()
yscale("log")
ax.set_xlim(left=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Sustained Throughput (S/ns)")
ylabel("Required Delay (ms)")
title("80 Channels")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_delay")

NaCsPlot.maybe_show()
