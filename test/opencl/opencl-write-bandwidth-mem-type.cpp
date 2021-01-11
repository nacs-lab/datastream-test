/*************************************************************************
 *   Copyright (c) 2020 - 2021 Yichao Yu <yyc1992@gmail.com>             *
 *                                                                       *
 *   This library is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU Lesser General Public          *
 *   License as published by the Free Software Foundation; either        *
 *   version 3.0 of the License, or (at your option) any later version.  *
 *                                                                       *
 *   This library is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *
 *   Lesser General Public License for more details.                     *
 *                                                                       *
 *   You should have received a copy of the GNU Lesser General Public    *
 *   License along with this library. If not,                            *
 *   see <http://www.gnu.org/licenses/>.                                 *
 *************************************************************************/

#include "helpers/test.h"
#include "helpers/opencl.h"

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <assert.h>

static YAML::Node test_device(cl::Device &dev, size_t nrep, size_t nele, size_t nalloc,
                              bool readable, bool host_access, bool host_write, bool host_ptr)
{
    assert(nalloc >= nele);
    assert((!host_ptr) || (host_access && host_write));
    assert((!host_write) || host_access);
    cl::Context ctx({dev});
    cl::CommandQueue queue(ctx, dev);
    // Use C++11 raw string literals for kernel source code
    std::string source{R"CLC(
        kernel void write(global int *res, int val)
        {
            res[get_global_id(0)] = val;
        }
    )CLC"};
    cl::Program prog(ctx, {source}, true);
    NaCs::Timer timer;
    auto flag = readable ? CL_MEM_READ_WRITE : CL_MEM_WRITE_ONLY;
    if (!host_access) {
        flag |= CL_MEM_HOST_NO_ACCESS;
    }
    else if (!host_write) {
        flag |= CL_MEM_HOST_READ_ONLY;
    }
    if (host_ptr)
        flag |= CL_MEM_ALLOC_HOST_PTR;
    cl::Buffer buff(ctx, flag, sizeof(float) * nalloc);

    {
        // Warm up
        cl::Kernel kernel(prog, "write");
        kernel.setArg(0, buff);
        kernel.setArg(1, int(42));
        cl::Event evt;
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evt);
        evt.wait();
    }

    timer.restart();
    cl::Event evt;
    for (size_t i = 0; i < nrep; i++) {
        cl::Kernel kernel(prog, "write");
        kernel.setArg(0, buff);
        kernel.setArg(1, int(42));
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evt);
    }
    evt.wait();
    auto t1 = (double)timer.elapsed() / double(nele) / double(nrep);

    auto res = OCL::get_device_ids(dev);
    res["t"] = t1;
    res["nrep"] = nrep;
    res["nele"] = nele;
    res["nalloc"] = nalloc;
    res["readable"] = readable;
    res["host_access"] = host_access;
    res["host_write"] = host_write;
    res["host_ptr"] = host_ptr;
    return res;
}

struct Config {
    YAML::Node dev_filter;
    size_t nrep;
    size_t nele;
    size_t nalloc;
    bool readable;
    bool host_access;
    bool host_write;
    bool host_ptr;
    static Config loadYAML(const char *fname)
    {
        Config conf;
        auto file = YAML::LoadFile(fname);
        auto required_key = [&] (auto name) {
            if (auto node = file[name])
                return node;
            throw std::runtime_error(std::string("Required key '") + name + "' missing.");
        };
        conf.nrep = required_key("nrep").as<size_t>();
        conf.nele = required_key("nele").as<size_t>();
        conf.nalloc = required_key("nalloc").as<size_t>();
        conf.readable = required_key("readable").as<bool>();
        conf.host_access = required_key("host_access").as<bool>();
        conf.host_write = required_key("host_write").as<bool>();
        conf.host_ptr = required_key("host_ptr").as<bool>();
        conf.dev_filter = file["devices"];
        return conf;
    }
};

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Needs at least one argument\n");
        exit(1);
    }
    auto config = Config::loadYAML(argv[1]);
    std::vector<cl::Device> devices = OCL::all_ocl2_devices(&config.dev_filter);
    if (devices.empty())
        throw std::runtime_error("Unable to find OpenCL devices");

    std::vector<YAML::Node> res;
    for (auto &dev: devices) {
        OCL::catch_error([&] {
            res.push_back(test_device(dev, config.nrep, config.nele, config.nalloc,
                                      config.readable, config.host_access,
                                      config.host_write, config.host_ptr));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
