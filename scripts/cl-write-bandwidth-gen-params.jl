#!/usr/bin/julia

function write_file(dir, nrep, nele, nwrite)
    name = "$nele-$nrep-$nwrite"
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, """
                    devices:
                      - cpu: false""")
        println(io, "nele: $nele")
        println(io, "nrep: $nrep")
        println(io, "nwrite: $nwrite")
    end
end

const dir = ARGS[1]

for i in 10:23
    nele = 2^i
    for nwrite in (1, 2, 3, 4)
        for nrep in (32, 64, 128, 256)
            write_file(dir, nrep, nele, nwrite)
        end
    end
end
