#!/usr/bin/julia

function write_file(dir, nrep, nele, nstream)
    name = "$nele-$nrep-$nstream"
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, """
                    devices:
                      - cpu: false""")
        println(io, "nele: $nele")
        println(io, "nrep: $nrep")
        println(io, "nstream: $nstream")
    end
end

const dir = ARGS[1]

for i in 14:23
    nele = 2^i
    for nstream in (1, 2, 3, 4, 8, 16)
        for nrep in (32, 64, 128, 256)
            write_file(dir, nrep, nele, nstream)
        end
    end
end
