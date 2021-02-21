#!/usr/bin/julia

function write_file(dir, nrep, nele, ncalc_single, ncalc_double)
    name = "$nele-$ncalc_single-$ncalc_double-$nrep"
    mkpath(dir, mode=0o755)
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, """
                    devices:
                      - cpu: false""")
        println(io, "nele: $nele")
        println(io, "ncalc_single: $ncalc_single")
        println(io, "ncalc_double: $ncalc_double")
        println(io, "nrep: $nrep")
    end
end

const dir = ARGS[1]

const nele = 2^23
for ncalc_single in (0, 64, 128, 256, 512)
    for ncalc_double in (0, 16, 32, 64, 128)
        if ncalc_double == 0 && ncalc_single == 0
            continue
        end
        write_file(dir, 64, nele, ncalc_single, ncalc_double)
    end
end
