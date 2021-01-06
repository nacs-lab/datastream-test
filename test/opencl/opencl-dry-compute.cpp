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

static YAML::Node test_device(cl::Device &dev, bool ooo, size_t nrep, size_t nele)
{
    cl::Context ctx({dev});
    cl::CommandQueue queue(ctx, dev, ooo ? CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE : 0);
    // Use C++11 raw string literals for kernel source code
    std::string source{R"CLC(
        kernel void dummy()
        {
        }

        kernel void compute(global float *res, float amp, float freq)
        {
            float val = amp * sin(freq * get_global_id(0));
            if (val > 2) {
                // We pass in an `amp < 2` so this will never happen
                // and we'll not be affected by any memory access bandwidth.
                // However, the compiler doesn't know about this
                // and will still run the function.
                res[get_global_id(0)] = val;
            }
        }
    )CLC"};
    cl::Program prog(ctx, {source}, true);
    std::vector<cl::Event> evts(nrep);
    NaCs::Timer timer;

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
        queue.enqueueMarkerWithWaitList(&evts, &evt);
        evt.wait();
    }
    auto t0 = (double)timer.elapsed() / double(nele) / double(nrep);

    {
        // Warm up
        cl::Kernel kernel(prog, "compute");
        kernel.setArg(0, nullptr);
        kernel.setArg(1, 0.2f);
        kernel.setArg(2, 0.002f);
        cl::Event evt;
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evt);
        evt.wait();
    }

    timer.restart();
    for (size_t i = 0; i < nrep; i++) {
        cl::Kernel kernel(prog, "compute");
        kernel.setArg(0, nullptr);
        kernel.setArg(1, 0.2f);
        kernel.setArg(2, 0.002f);
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evts[i]);
    }
    {
        cl::Event evt;
        queue.enqueueMarkerWithWaitList(&evts, &evt);
        evt.wait();
    }
    auto t1 = (double)timer.elapsed() / double(nele) / double(nrep);
    auto res = OCL::get_device_ids(dev);
    res["tdummy"] = t0;
    res["tcompute"] = t1;
    res["ooo"] = ooo;
    res["nrep"] = nrep;
    res["nele"] = nele;
    return res;
}

int main(void)
{
    std::vector<cl::Device> devices = OCL::all_ocl2_devices();
    if (devices.empty())
        throw std::runtime_error("Unable to find OpenCL 2.0 devices");

    std::vector<YAML::Node> res;
    for (auto &dev: devices) {
        OCL::catch_error([&] {
            res.push_back(test_device(dev, false, 1024, 1024 * 1024));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
