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

static YAML::Node test_device(cl::Device &dev, size_t nrep, size_t nele, size_t nconcurrent)
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
    auto flag = CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS;

    std::vector<cl::Buffer> buffs;
    for (int i = 0; i < nconcurrent * 2; i++)
        buffs.emplace_back(ctx, flag, sizeof(float) * nele);

    std::vector<cl::Event> evts(nconcurrent);

    {
        // Init
        for (size_t i = 0; i < nconcurrent; i++) {
            cl::Kernel init(prog, "write");
            init.setArg(0, buffs[i * 2]);
            queue.enqueueNDRangeKernel(init, cl::NDRange(), cl::NDRange(nele),
                                       cl::NDRange(), nullptr, &evts[i]);
        }
        // Warm up
        for (size_t i = 0; i < nconcurrent; i++) {
            cl::Kernel copy(prog, "copy");
            copy.setArg(0, buffs[i * 2 + 1]);
            copy.setArg(1, buffs[i * 2]);
            std::vector<cl::Event> wait{evts[i]};
            queue.enqueueNDRangeKernel(copy, cl::NDRange(), cl::NDRange(nele),
                                       cl::NDRange(), &wait, &evts[i]);
        }
        queue.finish();
    }

    std::vector<std::vector<uint64_t>> queue_times;
    std::vector<std::vector<uint64_t>> submit_times;
    std::vector<std::vector<uint64_t>> start_times;
    std::vector<std::vector<uint64_t>> complete_times;
    for (size_t i = 0; i < nconcurrent; i++) {
        queue_times.emplace_back(nrep, 0);
        submit_times.emplace_back(nrep, 0);
        start_times.emplace_back(nrep, 0);
        complete_times.emplace_back(nrep, 0);
    }

    cl::Kernel copy(prog, "copy");
    auto start_time = Test::cycleclock();
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < nconcurrent; i++) {
            if (j % 2 == 0) {
                // Forward
                copy.setArg(0, buffs[i * 2 + 1]);
                copy.setArg(1, buffs[i * 2]);
            }
            else {
                // Backward
                copy.setArg(0, buffs[i * 2]);
                copy.setArg(1, buffs[i * 2 + 1]);
            }
            std::vector<cl::Event> wait{evts[i]};
            queue_times[i][j] = Test::cycleclock() - start_time;
            queue.enqueueNDRangeKernel(copy, cl::NDRange(), cl::NDRange(nele),
                                       cl::NDRange(), &wait, &evts[i]);
            OCL::set_event_callback(evts[i], CL_SUBMITTED, [=, &submit_times] (cl_event,
                                                                               cl_int) {
                submit_times[i][j] = Test::cycleclock() - start_time;
            });
            OCL::set_event_callback(evts[i], CL_RUNNING, [=, &start_times] (cl_event, cl_int) {
                start_times[i][j] = Test::cycleclock() - start_time;
            });
            OCL::set_event_callback(evts[i], CL_COMPLETE, [=, &complete_times] (cl_event,
                                                                                cl_int) {
                complete_times[i][j] = Test::cycleclock() - start_time;
            });
        }
    }
    queue.finish();

    auto res = OCL::get_device_ids(dev);
    res["nrep"] = nrep;
    res["nele"] = nele;
    res["nconcurrent"] = nconcurrent;
    res["queue_times"] = queue_times;
    res["submit_times"] = submit_times;
    res["start_times"] = start_times;
    res["complete_times"] = complete_times;
    return res;
}

struct Config {
    YAML::Node dev_filter;
    size_t nrep;
    size_t nele;
    size_t nconcurrent;
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
        conf.nconcurrent = required_key("nconcurrent").as<size_t>();
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
            res.push_back(test_device(dev, config.nrep, config.nele, config.nconcurrent));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
