#!/usr/bin/julia

function size_to_str(size, unit="")
    size = Int(size)
    str(num) = isinteger(num) ? string(Int(num)) : string(num)
    if isempty(unit)
        _unit = ""
    else
        _unit = "i$unit"
    end
    if size < 1024
        if isempty(unit)
            return str(size)
        end
        return "$(str(size)) $unit"
    else
        size /= 1024
    end
    if size < 1024
        return "$(str(size)) k$_unit"
    else
        size /= 1024
    end
    if size < 1024
        return "$(str(size)) M$_unit"
    else
        size /= 1024
    end
    if size < 1024
        return "$(str(size)) G$_unit"
    else
        size /= 1024
    end
    if size < 1024
        return "$(str(size)) T$_unit"
    else
        size /= 1024
    end
    return "$(str(size)) P$_unit"
end

function num_to_si(num, unit)
    str(num) = isinteger(num) ? string(Int(num)) : string(num)
    if num >= 1e24
        return "$(str(num / 1e24)) Y$unit"
    elseif num >= 1e21
        return "$(str(num / 1e21)) Z$unit"
    elseif num >= 1e18
        return "$(str(num / 1e18)) E$unit"
    elseif num >= 1e15
        return "$(str(num / 1e15)) P$unit"
    elseif num >= 1e12
        return "$(str(num / 1e12)) T$unit"
    elseif num >= 1e9
        return "$(str(num / 1e9)) G$unit"
    elseif num >= 1e6
        return "$(str(num / 1e6)) M$unit"
    elseif num >= 1e3
        return "$(str(num / 1e3)) k$unit"
    elseif num >= 1
        if isempty(unit)
            return str(num)
        else
            return "$(str(num)) $unit"
        end
    elseif num >= 1e-3
        return "$(str(num * 1e3)) m$unit"
    elseif num >= 1e-6
        return "$(str(num * 1e6)) μ$unit"
    elseif num >= 1e-9
        return "$(str(num * 1e9)) n$unit"
    elseif num >= 1e-12
        return "$(str(num * 1e12)) p$unit"
    elseif num >= 1e-15
        return "$(str(num * 1e15)) f$unit"
    elseif num >= 1e-18
        return "$(str(num * 1e18)) a$unit"
    elseif num >= 1e-21
        return "$(str(num * 1e21)) z$unit"
    else # if num >= 1e-24
        return "$(str(num * 1e24)) y$unit"
    end
end
