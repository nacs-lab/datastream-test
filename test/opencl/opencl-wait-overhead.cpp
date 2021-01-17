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
#include <condition_variable>
#include <deque>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <assert.h>

static YAML::Node test_device(cl::Device &dev, size_t nrep, size_t nele, bool do_wait)
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


    std::vector<cl::Event> wait(1);
    {
        // Init
        cl::Kernel init(prog, "write");
        init.setArg(0, buffs[0]);
        queue.enqueueNDRangeKernel(init, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), nullptr, &wait[0]);
        // Warm up
        cl::Kernel copy(prog, "copy");
        copy.setArg(0, buffs[1]);
        copy.setArg(1, buffs[0]);
        queue.enqueueNDRangeKernel(copy, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), &wait, &wait[0]);
        queue.finish();
    }

    std::mutex lock;
    std::condition_variable cond_var;
    std::deque<cl::Event> evts;
    bool done = false;

    std::thread worker{[&] {
        if (!do_wait)
            return;
        std::unique_lock<std::mutex> locker(lock);
        while (!done) {
            cond_var.wait(locker, [&] {
                return done || !evts.empty();
            });
            while (!evts.empty()) {
                auto evt = evts.back();
                evts.pop_back();
                locker.unlock();
                evt.wait();
                locker.lock();
            }
        }
    }};

    cl::Kernel copy(prog, "copy");
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
        queue.enqueueNDRangeKernel(copy, cl::NDRange(), cl::NDRange(nele),
                                   cl::NDRange(), &wait, &wait[0]);
        {
            std::lock_guard<std::mutex> locker(lock);
            evts.push_front(wait[0]);
        }
        cond_var.notify_one();
    }
    queue.finish();
    auto t1 = (double)timer.elapsed() / double(nele) / double(nrep);
    {
        std::lock_guard<std::mutex> locker(lock);
        done = true;
    }
    cond_var.notify_one();
    worker.join();

    auto res = OCL::get_device_ids(dev);
    res["t"] = t1;
    res["nrep"] = nrep;
    res["nele"] = nele;
    res["do_wait"] = do_wait;
    return res;
}

struct Config {
    YAML::Node dev_filter;
    size_t nrep;
    size_t nele;
    bool do_wait;
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
        conf.do_wait = required_key("do_wait").as<bool>();
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
            res.push_back(test_device(dev, config.nrep, config.nele, config.do_wait));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
