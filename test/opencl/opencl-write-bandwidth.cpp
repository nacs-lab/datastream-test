/*************************************************************************
 *   Copyright (c) 2020 - 2020 Yichao Yu <yyc1992@gmail.com>             *
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

static YAML::Node test_device(cl::Device &dev, size_t nrep, size_t nele, int nwrite)
{
    cl::Context ctx({dev});
    cl::CommandQueue queue(ctx, dev);
    // Use C++11 raw string literals for kernel source code
    std::string source{R"CLC(
        kernel void write(global int *res, int nwrite, int val)
        {
            size_t base = get_global_id(0) * nwrite;
            for (int i = 0; i < nwrite; i++) {
                res[base + i] = val;
            }
        }
    )CLC"};
    auto gsize = nele / nwrite;
    nele = gsize * nwrite;
    cl::Program prog(ctx, {source}, true);
    NaCs::Timer timer;
    cl::Buffer buff(ctx, CL_MEM_WRITE_ONLY, sizeof(float) * nele);

    {
        // Warm up
        cl::Kernel kernel(prog, "write");
        kernel.setArg(0, buff);
        kernel.setArg(1, nwrite);
        kernel.setArg(2, int(42));
        cl::Event evt;
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(gsize),
                                   cl::NDRange(), nullptr, &evt);
        evt.wait();
    }

    timer.restart();
    cl::Event evt;
    for (size_t i = 0; i < nrep; i++) {
        cl::Kernel kernel(prog, "write");
        kernel.setArg(0, buff);
        kernel.setArg(1, nwrite);
        kernel.setArg(2, int(42));
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(gsize),
                                   cl::NDRange(), nullptr, &evt);
    }
    evt.wait();
    auto t1 = (double)timer.elapsed() / double(nele) / double(nrep);

    auto res = OCL::get_device_ids(dev);
    res["t"] = t1;
    res["nrep"] = nrep;
    res["nele"] = nele;
    res["nwrite"] = nwrite;
    return res;
}

struct Config {
    YAML::Node dev_filter;
    size_t nrep;
    size_t nele;
    int nwrite = 1;
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
        if (auto node = file["nwrite"])
            conf.nwrite = node.as<int>();
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
            res.push_back(test_device(dev, config.nrep, config.nele, config.nwrite));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
