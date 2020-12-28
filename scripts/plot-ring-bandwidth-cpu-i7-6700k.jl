#!/usr/bin/julia

include("plot-ring-bandwidth-cpu-csv.jl")
include("props_i7-6700k.jl")

const datadir = joinpath(@__DIR__, "../data/")

data_block = load_data(joinpath(datadir, "ring-bandwidth-cpu-i7-6700k_blockring_read.csv"),
                       joinpath(datadir, "ring-bandwidth-cpu-i7-6700k_blockring_write.csv"))

data_pipe = load_data(joinpath(datadir, "ring-bandwidth-cpu-i7-6700k_datapipe_read.csv"),
                      joinpath(datadir, "ring-bandwidth-cpu-i7-6700k_datapipe_write.csv"))

const prefix = joinpath(@__DIR__, "../imgs/ring-bandwidth-cpu-i7-6700k")

figure(figsize=[12.6, 5.6])

ax = subplot(1, 2, 1)
plot_cpu2cpu(get_cpu2cpu(data_block, ncores_phys))
title("BlockRing")

ax = subplot(1, 2, 2)
plot_cpu2cpu(get_cpu2cpu(data_pipe, ncores_phys))
title("DataPipe")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
