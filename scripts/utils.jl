#!/usr/bin/julia

function size_to_str(size)
    size = Int(size)
    str(num) = isinteger(num) ? string(Int(num)) : string(num)
    if size < 1024
        return str(size)
    else
        size /= 1024
    end
    if size < 1024
        return "$(str(size)) k"
    else
        size /= 1024
    end
    if size < 1024
        return "$(str(size)) M"
    else
        size /= 1024
    end
    if size < 1024
        return "$(str(size)) G"
    else
        size /= 1024
    end
    if size < 1024
        return "$(str(size)) T"
    else
        size /= 1024
    end
    return "$(str(size)) P"
end
