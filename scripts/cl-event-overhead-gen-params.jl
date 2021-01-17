#!/usr/bin/julia

function write_file(dir, nele, submit_cb, start_cb, complete_cb)
    name = "$nele"
    if submit_cb
        name = "$name-submit"
    end
    if start_cb
        name = "$name-start"
    end
    if complete_cb
        name = "$name-complete"
    end
    mkpath(dir, mode=0o755)
    open(joinpath(dir, name * ".yaml"), "w") do io
        println(io, """
                    devices:
                      - cpu: false
                    nrep: 1024""")
        println(io, "nele: $nele")
        println(io, "submit_cb: $submit_cb")
        println(io, "start_cb: $start_cb")
        println(io, "complete_cb: $complete_cb")
    end
end

const dir = ARGS[1]

for i in 10:23
    nele = 2^i
    for submit_cb in (false, true)
        for start_cb in (false, true)
            for complete_cb in (false, true)
                write_file(dir, nele, submit_cb, start_cb, complete_cb)
            end
        end
    end
end
