#!/usr/bin/julia

include("plot-ring-bandwidth-compute-csv.jl")
include("props_i9-7900x.jl")

const datadir = joinpath(@__DIR__, "../data/")

data_block = load_data(joinpath(datadir, "ring-bandwidth-compute-i9-7900x_read.csv"),
                       joinpath(datadir, "ring-bandwidth-compute-i9-7900x_write.csv"))

const prefix = joinpath(@__DIR__, "../imgs/ring-bandwidth-compute-i9-7900x")

const nblks = [8, 16, 32]
const localbuffs = [0, 32 * 1024]

function draw_l3_percore(y)
    axvline(L3_total / ncores_phys, color="C2", linewidth=3, ls="dotted")
    text(L3_total / ncores_phys, y, "L3/core", rotation=90,
         va="center", ha="right", color="C2", fontsize=19)
end

function set_ylim()
    ylim([7.7, 9.65])
end

function draw_max(y, x)
    axhline(y, color="r", linewidth=2, ls="dotted")
    text(x, y, "$y", va="top", ha="center", color="r", fontsize=19)
end

figure(figsize=[12.6, 16.8])

ax = subplot(3, 2, 1)
plot_perf(data_block, nblks, localbuffs, 2, 2)
draw_l3_percore(9.0)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 2)
plot_stall(data_block, nblks, localbuffs, 2, 2)
draw_l3_percore(20)

ax = subplot(3, 2, 3)
plot_perf(data_block, nblks, localbuffs, 4, 4)
draw_l3_percore(8.3)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 4)
plot_stall(data_block, nblks, localbuffs, 4, 4)
draw_l3_percore(78)

ax = subplot(3, 2, 5)
plot_perf(data_block, nblks, localbuffs, 8, 8)
draw_l3_percore(8.5)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 6)
plot_stall(data_block, nblks, localbuffs, 8, 8)
draw_l3_percore(65)

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_r2w2")

figure(figsize=[12.6, 16.8])

ax = subplot(3, 2, 1)
plot_perf(data_block, nblks, localbuffs, 3, 1)
draw_l3_percore(8.25)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 2)
plot_stall(data_block, nblks, localbuffs, 3, 1)
draw_l3_percore(50)

ax = subplot(3, 2, 3)
plot_perf(data_block, nblks, localbuffs, 6, 2)
draw_l3_percore(8.4)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 4)
plot_stall(data_block, nblks, localbuffs, 6, 2)
draw_l3_percore(50)

ax = subplot(3, 2, 5)
plot_perf(data_block, nblks, localbuffs, 12, 4)
draw_l3_percore(8.5)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 6)
plot_stall(data_block, nblks, localbuffs, 12, 4)
draw_l3_percore(50)

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_r3w1")

figure(figsize=[12.6, 16.8])

ax = subplot(3, 2, 1)
plot_perf(data_block, nblks, localbuffs, 1, 3)
draw_l3_percore(8.25)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 2)
plot_stall(data_block, nblks, localbuffs, 1, 3)
draw_l3_percore(20)

ax = subplot(3, 2, 3)
plot_perf(data_block, nblks, localbuffs, 2, 6)
draw_l3_percore(8.4)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 4)
plot_stall(data_block, nblks, localbuffs, 2, 6)
draw_l3_percore(50)

ax = subplot(3, 2, 5)
plot_perf(data_block, nblks, localbuffs, 4, 12)
draw_l3_percore(8.5)
set_ylim()
draw_max(9.62, 2.0^19)

ax = subplot(3, 2, 6)
plot_stall(data_block, nblks, localbuffs, 4, 12)
draw_l3_percore(50)

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_r1w3")

NaCsPlot.maybe_show()
