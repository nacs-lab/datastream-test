#!/usr/bin/julia

function write_file(dir, nele, nconcurrent)
    name = "$nele-$nconcurrent"
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, """
                    devices:
                      - cpu: false
                    nrep: 1024""")
        println(io, "nele: $nele")
        println(io, "nconcurrent: $nconcurrent")
    end
end

const dir = ARGS[1]

for i in 14:23
    nele = 2^i
    for nconcurrent in (1, 2, 3, 4, 8, 16)
        write_file(dir, nele, nconcurrent)
    end
end
