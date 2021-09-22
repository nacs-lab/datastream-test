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

static YAML::Node test_device(cl::Device &dev, bool ooo, size_t nrep, size_t nele, int ncalc)
{
    cl::Context ctx({dev});
    cl::CommandQueue queue(ctx, dev, ooo ? CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE : 0);
    // Use C++11 raw string literals for kernel source code
    std::string source{R"CLC(
        kernel void dummy()
        {
        }

        kernel void compute(global float *res, float amp, float freq, int ncalc)
        {
            // The performance of sine function may depend on the input value.
            // We make sure the range is small in the real code and we should do that here
            // as well to make the performance more realistic.
            float val = get_global_id(0) % 512;
            for (int i = 0; i < ncalc; i++)
                val = amp * sin(freq * val);
            if (val > 1024) {
                // We pass in an `amp < 1024` so this will never happen
                // and we'll not be affected by any memory access bandwidth.
                // However, the compiler doesn't know about this
                // and will still run the function.
                res[get_global_id(0)] = val;
            }
        }

        kernel void compute_native(global float *res, float amp, float freq, int ncalc)
        {
            float val = get_global_id(0) % 512;
            for (int i = 0; i < ncalc; i++)
                val = amp * native_sin(freq * val);
            if (val > 1024) {
                res[get_global_id(0)] = val;
            }
        }
    )CLC"};
    cl::Program prog(ctx, {source}, true);
    std::vector<cl::Event> evts(nrep);
    auto marker_wait = &evts;
    NaCs::Timer timer;

    auto res = OCL::get_device_ids(dev);
    // Somehow the NVIDIA OpenCL driver (470.74)
    // doesn't actually wait for any of the provides events.
    // Pass in null and let it wait for everything previously queued works though...
    if (res["platform_name"].as<std::string>() == "NVIDIA CUDA")
        marker_wait = nullptr;

    {
        // Warm up
        cl::Kernel kernel(prog, "dummy");
        cl::Event evt;
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evt);
        evt.wait();
    }

    timer.restart();
    for (size_t i = 0; i < nrep; i++) {
        cl::Kernel kernel(prog, "dummy");
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evts[i]);
    }
    {
        cl::Event evt;
        queue.enqueueMarkerWithWaitList(marker_wait, &evt);
        evt.wait();
    }
    auto t0 = (double)timer.elapsed() / double(nele) / double(nrep);

    {
        // Warm up
        cl::Kernel kernel(prog, "compute");
        kernel.setArg(0, nullptr);
        kernel.setArg(1, 512.0f);
        kernel.setArg(2, 0.8f / 256);
        kernel.setArg(3, ncalc);
        cl::Event evt;
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evt);
        evt.wait();
    }

    timer.restart();
    for (size_t i = 0; i < nrep; i++) {
        cl::Kernel kernel(prog, "compute");
        kernel.setArg(0, nullptr);
        kernel.setArg(1, 512.0f);
        kernel.setArg(2, 0.8f / 256);
        kernel.setArg(3, ncalc);
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evts[i]);
    }
    {
        cl::Event evt;
        queue.enqueueMarkerWithWaitList(marker_wait, &evt);
        evt.wait();
    }
    auto t1 = (double)timer.elapsed() / double(nele) / double(nrep);

    {
        // Warm up
        cl::Kernel kernel(prog, "compute_native");
        kernel.setArg(0, nullptr);
        kernel.setArg(1, 512.0f);
        kernel.setArg(2, 0.8f / 256);
        kernel.setArg(3, ncalc);
        cl::Event evt;
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evt);
        evt.wait();
    }

    timer.restart();
    for (size_t i = 0; i < nrep; i++) {
        cl::Kernel kernel(prog, "compute_native");
        kernel.setArg(0, nullptr);
        kernel.setArg(1, 512.0f);
        kernel.setArg(2, 0.8f / 256);
        kernel.setArg(3, ncalc);
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evts[i]);
    }
    {
        cl::Event evt;
        queue.enqueueMarkerWithWaitList(marker_wait, &evt);
        evt.wait();
    }
    auto t2 = (double)timer.elapsed() / double(nele) / double(nrep);

    res["tdummy"] = t0;
    res["tcompute"] = t1;
    res["tcompute_native"] = t2;
    res["ooo"] = ooo;
    res["nrep"] = nrep;
    res["nele"] = nele;
    res["ncalc"] = ncalc;
    return res;
}

struct Config {
    YAML::Node dev_filter;
    bool ooo;
    size_t nrep;
    size_t nele;
    int ncalc = 1;
    static Config loadYAML(const char *fname)
    {
        Config conf;
        auto file = YAML::LoadFile(fname);
        auto required_key = [&] (auto name) {
            if (auto node = file[name])
                return node;
            throw std::runtime_error(std::string("Required key '") + name + "' missing.");
        };
        conf.ooo = required_key("ooo").as<bool>();
        conf.nrep = required_key("nrep").as<size_t>();
        conf.nele = required_key("nele").as<size_t>();
        if (auto node = file["ncalc"])
            conf.ncalc = node.as<int>();
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
            res.push_back(test_device(dev, config.ooo, config.nrep, config.nele, config.ncalc));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
