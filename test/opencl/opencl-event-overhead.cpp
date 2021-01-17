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

#include <atomic>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <assert.h>

static YAML::Node test_device(cl::Device &dev, size_t nrep, size_t nele,
                              bool submit_cb, bool start_cb, bool complete_cb)
{
    cl::Context ctx({dev});
    cl::CommandQueue queue(ctx, dev, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
    // Use C++11 raw string literals for kernel source code
    std::string source{R"CLC(
        kernel void write(global int *res)
        {
            res[get_global_id(0)] = (int)get_global_id(0);
        }
        kernel void copy(global int *res, global const int *in)
        {
            res[get_global_id(0)] = in[get_global_id(0)];
        }
    )CLC"};
    cl::Program prog(ctx, {source}, true);
    NaCs::Timer timer;
    auto flag = CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS;

    std::vector<cl::Buffer> buffs;
    for (int i = 0; i < 2; i++)
        buffs.emplace_back(ctx, flag, sizeof(float) * nele);

    cl::Event evt;

    {
        // Init
        cl::Kernel init(prog, "write");
        init.setArg(0, buffs[0]);
        queue.enqueueNDRangeKernel(init, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &evt);
        // Warm up
        cl::Kernel copy(prog, "copy");
        copy.setArg(0, buffs[1]);
        copy.setArg(1, buffs[0]);
        std::vector<cl::Event> wait{evt};
        queue.enqueueNDRangeKernel(copy, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), &wait, &evt);
        queue.finish();
    }

    cl::Kernel copy(prog, "copy");
    auto dummy_cb = [] (cl_event, cl_int) {};
    timer.restart();
    for (size_t j = 0; j < nrep; j++) {
        if (j % 2 == 0) {
            // Forward
            copy.setArg(0, buffs[1]);
            copy.setArg(1, buffs[0]);
        }
        else {
            // Backward
            copy.setArg(0, buffs[0]);
            copy.setArg(1, buffs[1]);
        }
        std::vector<cl::Event> wait{evt};
        queue.enqueueNDRangeKernel(copy, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), &wait, &evt);
        if (submit_cb)
            OCL::set_event_callback(evt, CL_SUBMITTED, dummy_cb);
        if (start_cb)
            OCL::set_event_callback(evt, CL_RUNNING, dummy_cb);
        if (complete_cb)
            OCL::set_event_callback(evt, CL_COMPLETE, dummy_cb);
    }
    queue.finish();
    auto t1 = (double)timer.elapsed() / double(nele) / double(nrep);

    auto res = OCL::get_device_ids(dev);
    res["t"] = t1;
    res["nrep"] = nrep;
    res["nele"] = nele;
    res["submit_cb"] = submit_cb;
    res["start_cb"] = start_cb;
    res["complete_cb"] = complete_cb;
    return res;
}

struct Config {
    YAML::Node dev_filter;
    size_t nrep;
    size_t nele;
    bool submit_cb;
    bool start_cb;
    bool complete_cb;
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
        conf.submit_cb = required_key("submit_cb").as<bool>();
        conf.start_cb = required_key("start_cb").as<bool>();
        conf.complete_cb = required_key("complete_cb").as<bool>();
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
            res.push_back(test_device(dev, config.nrep, config.nele, config.submit_cb,
                                      config.start_cb, config.complete_cb));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
