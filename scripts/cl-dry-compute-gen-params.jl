#!/usr/bin/julia

function write_file(dir, ooo, nrep, nele, ncalc)
    name = "$nele-$ncalc-$nrep"
    if ooo
        name = name * "-ooo"
    end
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, "ooo: $(ooo ? "true" : "false")")
        println(io, "nele: $nele")
        println(io, "ncalc: $ncalc")
        println(io, "nrep: $nrep")
    end
end

const dir = ARGS[1]

# For smallest sizes (1 - 128), both run a few times and run larger number of times
# to measure the overhead
for i in 1:7
    nele = 2^i
    ncalc = 1
    for nrep in (1, 2, 4, 8, 16, 2^15, 2^16, 2^17)
        write_file(dir, true, nrep, nele, ncalc)
        write_file(dir, false, nrep, nele, ncalc)
    end
end

# Medium sizes (256 - 32k)
for i in 8:15
    nele = 2^i
    for ncalc in (1, 2, 4, 8)
        for nrep in (2^17, 2^18, 2^19) .÷ nele
            write_file(dir, true, nrep, nele, ncalc)
            write_file(dir, false, nrep, nele, ncalc)
        end
    end
end

# Large sizes (64k - 8M)
for i in 16:23
    nele = 2^i
    for ncalc in (1, 2, 3, 4, 6, 8, 16, 32, 64)
        for nrep in (32, 64, 128)
            write_file(dir, true, nrep, nele, ncalc)
            write_file(dir, false, nrep, nele, ncalc)
        end
    end
end
