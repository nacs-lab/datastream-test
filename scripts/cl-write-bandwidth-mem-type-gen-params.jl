#!/usr/bin/julia

function write_file(dir, nrep, nele, nalloc, readable::Bool,
                    host_access::Bool, host_write::Bool, host_ptr::Bool)
    name = "$nele-$nrep-$nalloc"
    if readable
        name = name * "-dev_read"
    end
    if host_ptr
        name = name * "-host_ptr"
    elseif host_write
        name = name * "-host_write"
    elseif host_access
        name = name * "-host_access"
    end
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, """
                    devices:
                      - cpu: false""")
        println(io, "nele: $nele")
        println(io, "nrep: $nrep")
        println(io, "nalloc: $nalloc")
        println(io, "readable: $readable")
        println(io, "host_access: $host_access")
        println(io, "host_write: $host_write")
        println(io, "host_ptr: $host_ptr")
    end
end

const dir = ARGS[1]

for i in 12:23
    nele = 2^i
    if nele < 2^20
        nallocs = (nele, 2^20, 2^23)
    elseif nele < 2^23
        nallocs = (nele, 2^23)
    else
        nallocs = (nele,)
    end
    for nalloc in nallocs
        for nrep in (32, 64, 128, 256)
            for readable in (false, true)
                write_file(dir, nrep, nele, nalloc, readable, false, false, false)
                write_file(dir, nrep, nele, nalloc, readable, true, false, false)
                write_file(dir, nrep, nele, nalloc, readable, true, true, false)
                write_file(dir, nrep, nele, nalloc, readable, true, true, true)
            end
        end
    end
end
