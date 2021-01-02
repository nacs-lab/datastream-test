#!/usr/bin/julia

include("plot-multi-threads-compute.jl")
include("props_i9-7900x.jl")

const datadir = joinpath(@__DIR__, "../data/multi-threads-compute/i9-7900x-2")

const data = load_all_times(datadir, nominal_ghz)

const prefix = joinpath(@__DIR__, "../imgs/multi-threads-compute-i9-7900x-2")

figure(figsize=[12.6, 11.2])

ax = subplot(2, 2, 1)
plot_lines(data["20chn-nolocal"], 16384 / 625e3, 16384,
           "C0", "No Buffer", ls="--", alpha=0.6)
plot_lines(data["20chn-local"], 16384 / 625e3, 16384,
           "C1", "16 KiB Buffer", alpha=0.3)
plot_lines(data["20chn-local2"], 16384 / 625e3, 16384,
           "C3", "32 KiB Buffer", alpha=0.3)
axhline(0.625, color="r")
text(1000, 0.625, "625 MS/s", va="bottom", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="C2", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Delay (ms)")
ylabel("Minimum Throughput (S/ns)")
title("20 Channels")

ax = subplot(2, 2, 2)
plot_lines(data["50chn-nolocal"], 16384 / 625e3, 16384,
           "C0", "No Buffer", ls="--", alpha=0.6)
plot_lines(data["50chn-local"], 16384 / 625e3, 16384,
           "C1", "16 KiB Buffer", alpha=0.3)
plot_lines(data["50chn-local2"], 16384 / 625e3, 16384,
           "C3", "32 KiB Buffer", alpha=0.3)
axhline(0.625, color="r")
text(3300, 0.625, "625 MS/s", va="bottom", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="C2", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Delay (ms)")
ylabel("Minimum Throughput (S/ns)")
title("50 Channels")

ax = subplot(2, 2, 3)
plot_lines(data["100chn-nolocal"], 16384 / 625e3, 16384,
           "C0", "No Buffer", ls="--", alpha=0.6)
plot_lines(data["100chn-local"], 16384 / 625e3, 16384,
           "C1", "16 KiB Buffer", alpha=0.3)
plot_lines(data["100chn-local2"], 16384 / 625e3, 16384,
           "C3", "32 KiB Buffer", alpha=0.3)
axhline(0.625, color="r")
text(0.13, 0.625, "625 MS/s", va="bottom", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="C2", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Delay (ms)")
ylabel("Minimum Throughput (S/ns)")
title("100 Channels")

ax = subplot(2, 2, 4)
plot_lines(data["200chn-nolocal"], 16384 / 625e3, 16384,
           "C0", "No Buffer", ls="--", alpha=0.6)
plot_lines(data["200chn-local"], 16384 / 625e3, 16384,
           "C1", "16 KiB Buffer", alpha=0.3)
plot_lines(data["200chn-local2"], 16384 / 625e3, 16384,
           "C3", "32 KiB Buffer", alpha=0.3)
axhline(0.625, color="r")
text(1000, 0.615, "625 MS/s", va="top", ha="center", color="r", fontsize=19)
axvline(1024^2 / 625e3, color="C2", linewidth=2)
grid()
xscale("log")
ax.set_ylim(bottom=0)
legend(fontsize="x-small", columnspacing=1, handlelength=1)
xlabel("Delay (ms)")
ylabel("Minimum Throughput (S/ns)")
title("200 Channels")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
