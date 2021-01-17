#!/usr/bin/julia

function write_file(dir, nele, do_wait)
    name = "$nele"
    if do_wait
        name = "$name-wait"
    end
    mkpath(dir, mode=0o755)
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, """
                    devices:
                      - cpu: false
                    nrep: 1024""")
        println(io, "nele: $nele")
        println(io, "do_wait: $do_wait")
    end
end

const dir = ARGS[1]

for i in 10:23
    nele = 2^i
    for do_wait in (false, true)
        write_file(dir, nele, do_wait)
    end
end
