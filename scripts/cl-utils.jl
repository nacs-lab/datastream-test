#!/usr/bin/julia

function get_intel_device_id(res)
    name = res["name"]
    m = match(r"Core\(TM\) (i[-0-9A-Za-z]*) ", name)
    if m !== nothing
        return "intel-core-" * lowercase(m[1])
    end
    m = match(r"Graphics \[0x([0-9A-Fa-f]*)\]", name)
    if m !== nothing
        return "intel-gpu-" * lowercase(m[1])
    end
    return name
end

function get_device_id(res)
    if res["vendor"] == "Intel(R) Corporation"
        return get_intel_device_id(res)
    elseif res["vendor"] == "Advanced Micro Devices, Inc."
        return "amd-$(res["type"])-$(res["name"])"
    end
    return res["name"]
end
