#!/usr/bin/julia

const worker_config_short = """
workers:
  - # 0
  - # 1
  - # 2
  - # 3
  - # 4
    input: [0, 1]
  - # 5
    input: [2, 3]
  - # 6
    input: [4, 5]
    final: true"""

const worker_config_long = """
workers:
  - # 0
  - # 1
  - # 2
  - # 3
  - # 4
  - # 5
  - # 6
  - # 7
  - # 8
  - # 9
  - # 10
  - # 11
  - # 12
  - # 13
  - # 14
  - # 15
  - # 16
  - # 17
  - # 18
    input: [0, 1, 2]
  - # 19
    input: [3, 4, 5]
  - # 20
    input: [6, 7, 8]
  - # 21
    input: [9, 10, 11]
  - # 22
    input: [12, 13, 14]
  - # 23
    input: [15, 16, 17]
  - # 24
    input: [18, 19]
  - # 25
    input: [20, 21]
  - # 26
    input: [22, 23]
  - # 27
    input: [24, 25, 26]
    final: true"""

function write_file(dir, nele, long, do_wait, complete_cb)
    name = "$nele"
    if long
        name = "$name-long"
    else
        name = "$name-short"
    end
    if do_wait
        name = "$name-wait"
    end
    if complete_cb
        name = "$name-cb"
    end
    mkpath(dir, mode=0o755)
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, """
                    devices:
                      - cpu: false""")
        println(io, long ? worker_config_long : worker_config_short)
        println(io, "nrep: $(long ? 256 : 512)")
        println(io, "nele: $nele")
        println(io, "nbuffer: 4")
        println(io, "do_wait: $do_wait")
        println(io, "complete_cb: $complete_cb")
    end
end

const dir = ARGS[1]

for i in 10:23
    nele = 2^i
    for long in (false, true)
        for do_wait in (false, true)
            for complete_cb in (false, true)
                write_file(dir, nele, long, do_wait, complete_cb)
            end
        end
    end
end
