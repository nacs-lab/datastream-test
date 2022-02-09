#!/usr/bin/julia

using YAML

function create_worker_tree(cpus)
    ncpus = length(cpus)
    if ncpus == 1
        return Dict(cpus[1]=>Int[])
    elseif ncpus == 2
        return Dict(cpus[1]=>Int[], cpus[2]=>Int[cpus[1]])
    elseif ncpus == 3
        return Dict(cpus[1]=>Int[], cpus[2]=>Int[], cpus[3]=>Int[cpus[1], cpus[2]])
    elseif ncpus == 4
        return Dict(cpus[1]=>Int[], cpus[2]=>Int[], cpus[3]=>Int[],
                    cpus[4]=>Int[cpus[1], cpus[2], cpus[3]])
    end
    half_ncpus = ncpus ÷ 2
    return merge!(Dict(cpus[end]=>Int[cpus[half_ncpus], cpus[ncpus - 1]]),
                  create_worker_tree(cpus[1:half_ncpus]),
                  create_worker_tree(cpus[half_ncpus + 1:ncpus - 1]))
end

const dir = ARGS[1]
const config = YAML.load_file(ARGS[2])
const ncores = config["ncores"]
const worker_tree = sort(collect(create_worker_tree(1:ncores)))
const total_block = config["total_block"]
const block_size = config["block_size"]
const num_block = config["num_block"]
const max_num_block = config["max_num_block"]

function write_file(nchn, num_block_final)
    name = "$(nchn * ncores)chn-$(num_block_final * block_size)"
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, "total_block: $total_block")
        println(io, "block_size: $block_size")
        println(io, "num_block: $num_block")
        println(io, "num_block_final: $num_block_final")
        println(io, "num_channel: $nchn")
        println(io, "localbuff: 0")
        println(io, "workers:")
        for (cpu, inputs) in worker_tree
            println(io, "  - cpu: $cpu")
            if !isempty(inputs)
                println(io, "    input_cpu: $inputs")
            end
        end
        println(io, "    final: true")
    end
end

for nchn in [1, 2, 4, 8]
    num_block_final = max_num_block
    while num_block_final >= num_block
        write_file(nchn, num_block_final)
        num_block_final ÷= 4
    end
end
