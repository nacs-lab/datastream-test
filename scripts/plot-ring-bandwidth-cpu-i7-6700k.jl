#!/usr/bin/julia

include("plot-ring-bandwidth-cpu-csv.jl")
include("props_i7-6700k.jl")

const datadir = joinpath(@__DIR__, "../data/cpu-ring-bandwidth-cpu/i7-6700k")

data_block = load_data(joinpath(datadir, "blockring_read.csv"),
                       joinpath(datadir, "blockring_write.csv"))

data_pipe = load_data(joinpath(datadir, "datapipe_read.csv"),
                      joinpath(datadir, "datapipe_write.csv"))

const prefix = joinpath(@__DIR__, "../imgs/ring-bandwidth-cpu-i7-6700k")

fig = figure(figsize=[12.6, 5.6])

ax = subplot(1, 2, 1)
plot_cpu2cpu(ax, get_cpu2cpu(data_block, ncores_phys))
ax.set_title("BlockRing")

ax = subplot(1, 2, 2)
plot_cpu2cpu(ax, get_cpu2cpu(data_pipe, ncores_phys))
ax.set_title("DataPipe")

fig.subplots_adjust(wspace=0.45)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
