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

struct WorkerConfig {
    unsigned nins = 0;
    unsigned ins[3];
    bool is_final = false;
    bool check(unsigned self) const
    {
        if (nins == 0)
            return true;
        // Require ordering to make sure there's no duplicate or loop
        if (nins == 2)
            return (ins[0] < ins[1] && ins[1] < self);
        if (nins == 3)
            return (ins[0] < ins[1] && ins[1] < ins[2] && ins[2] < self);
        return false;
    }
};

struct TestConfig {
    size_t nrep;
    size_t nele;
    size_t nbuffer;
    bool complete_cb;
    bool do_wait;
    std::vector<WorkerConfig> workers;
    bool check() const
    {
        if (nrep == 0 || nele == 0 || nbuffer == 0 || workers.empty())
            return false;
        auto nworkers = workers.size();
        std::vector<bool> consumed(nworkers, false);
        for (size_t i = 0; i < nworkers; i++) {
            const auto &wc = workers[i];
            if (!wc.check(unsigned(i)))
                return false;
            if (i == nworkers - 1) {
                if (!wc.is_final)
                    return false;
                consumed[i] = true;
            }
            else if (wc.is_final) {
                return false;
            }
            for (unsigned j = 0; j < wc.nins; j++) {
                if (consumed[wc.ins[j]])
                    return false;
                consumed[wc.ins[j]] = 1;
            }
        }
        for (auto v: consumed) {
            if (!v) {
                return false;
            }
        }
        return true;
    }
};

static YAML::Node test_device(cl::Device &dev, const TestConfig &config)
{
    cl::Context ctx({dev});
    cl::CommandQueue queue(ctx, dev, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
    // Use C++11 raw string literals for kernel source code
    std::string source{R"CLC(
        kernel void generate(global float *res)
        {
            res[get_global_id(0)] = (float)(get_global_id(0) % 4096);
        }
        kernel void merge(global float *res, global const float *in1, global const float *in2)
        {
            res[get_global_id(0)] = in1[get_global_id(0)] + in2[get_global_id(0)];
        }
        kernel void merge3(global float *res, global const float *in1,
                           global const float *in2, global const float *in3)
        {
            res[get_global_id(0)] = (in1[get_global_id(0)] + in2[get_global_id(0)] +
                                     in3[get_global_id(0)]);
        }
    )CLC"};
    cl::Program prog(ctx, {source}, true);
    NaCs::Timer timer;
    auto normal_flag = CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS;
    auto final_flag = CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR;

    auto nworkers = config.workers.size();

    std::vector<std::vector<cl::Buffer>> all_buffs(config.nbuffer);
    for (auto &buffs: all_buffs) {
        for (size_t i = 0; i < nworkers; i++) {
            auto flag = i == nworkers - 1 ? final_flag : normal_flag;
            buffs.emplace_back(ctx, flag, sizeof(float) * config.nele);
        }
    }
    std::vector<std::vector<cl::Event>> all_evts(config.nbuffer);
    for (auto &evts: all_evts)
        evts.resize(nworkers);

    cl::Kernel generate(prog, "generate");
    cl::Kernel merge(prog, "merge");
    cl::Kernel merge3(prog, "merge3");

    auto run_once = [&] (unsigned i, unsigned prev_i) {
        auto &buffs = all_buffs[i];
        auto &evts = all_evts[i];
        auto &prev_evts = all_evts[prev_i];

        for (size_t wi = 0; wi < nworkers; wi++) {
            auto &wc = config.workers[wi];
            std::vector<cl::Event> wait_vec;
            auto add_wait = [&] (const cl::Event &evt) {
                if (evt.get()) {
                    wait_vec.push_back(evt);
                }
            };
            add_wait(prev_evts[wi]);
            for (unsigned j = 0; j < wc.nins; j++)
                add_wait(evts[wc.ins[j]]);
            auto wait = wait_vec.empty() ? nullptr : &wait_vec;
            if (wc.nins == 0) {
                generate.setArg(0, buffs[wi]);
                queue.enqueueNDRangeKernel(generate, cl::NDRange(), cl::NDRange(config.nele),
                                           cl::NDRange(), wait, &evts[wi]);
            }
            else if (wc.nins == 2) {
                merge.setArg(0, buffs[wi]);
                merge.setArg(1, buffs[wc.ins[0]]);
                merge.setArg(2, buffs[wc.ins[1]]);
                queue.enqueueNDRangeKernel(merge, cl::NDRange(), cl::NDRange(config.nele),
                                           cl::NDRange(), wait, &evts[wi]);
            }
            else {
                merge3.setArg(0, buffs[wi]);
                merge3.setArg(1, buffs[wc.ins[0]]);
                merge3.setArg(2, buffs[wc.ins[1]]);
                merge3.setArg(3, buffs[wc.ins[2]]);
                queue.enqueueNDRangeKernel(merge3, cl::NDRange(), cl::NDRange(config.nele),
                                           cl::NDRange(), wait, &evts[wi]);
            }
        }
        std::vector<cl::Event> wait_vec{evts.back()};
        cl::Event evt;
        queue.enqueueMapBuffer(buffs.back(), false, CL_MAP_READ,
                               0, sizeof(float) * config.nele, &wait_vec, &evt);
        return evt;
    };

    // Warm up
    run_once(0, unsigned(config.nbuffer - 1));
    queue.finish();

    auto dummy_cb = [] (cl_event, cl_int) {};

    std::mutex lock;
    std::condition_variable cond_var;
    std::deque<cl::Event> wait_queue;
    bool done = false;

    std::thread worker{[&] {
        if (!config.do_wait)
            return;
        std::unique_lock<std::mutex> locker(lock);
        while (!done) {
            cond_var.wait(locker, [&] {
                return done || !wait_queue.empty();
            });
            while (!wait_queue.empty()) {
                auto evt = wait_queue.back();
                wait_queue.pop_back();
                locker.unlock();
                evt.wait();
                locker.lock();
            }
        }
    }};

    timer.restart();
    unsigned prevbuff_idx = 0;
    for (size_t j = 0; j < config.nrep; j++) {
        auto buff_idx = prevbuff_idx + 1;
        if (buff_idx >= config.nbuffer)
            buff_idx = 0;
        auto evt = run_once(buff_idx, prevbuff_idx);
        prevbuff_idx = buff_idx;
        if (config.complete_cb)
            OCL::set_event_callback(evt, CL_COMPLETE, dummy_cb);
        if (config.do_wait) {
            {
                std::lock_guard<std::mutex> locker(lock);
                wait_queue.push_front(evt);
            }
            cond_var.notify_one();
        }
    }
    queue.finish();
    auto t1 = (double)timer.elapsed() / double(config.nele) / double(config.nrep);
    {
        std::lock_guard<std::mutex> locker(lock);
        done = true;
    }
    cond_var.notify_one();
    worker.join();

    auto res = OCL::get_device_ids(dev);
    res["t"] = t1;
    res["nrep"] = config.nrep;
    res["nele"] = config.nele;
    res["nbuffer"] = config.nbuffer;
    res["nworker"] = config.workers.size();
    res["do_wait"] = config.do_wait;
    res["complete_cb"] = config.complete_cb;
    return res;
}

struct Config {
    YAML::Node dev_filter;
    TestConfig test;
    static Config loadYAML(const char *fname)
    {
        Config conf;
        auto file = YAML::LoadFile(fname);
        auto required_key = [&] (auto name) {
            if (auto node = file[name])
                return node;
            throw std::runtime_error(std::string("Required key '") + name + "' missing.");
        };
        conf.test.nrep = required_key("nrep").as<size_t>();
        conf.test.nele = required_key("nele").as<size_t>();
        conf.test.nbuffer = required_key("nbuffer").as<size_t>();
        conf.test.do_wait = required_key("do_wait").as<bool>();
        conf.test.complete_cb = required_key("complete_cb").as<bool>();
        auto workers = required_key("workers");
        if (!workers.IsSequence())
            throw std::runtime_error("Invalid workers config.");
        for (const auto &worker: workers) {
            if (worker.IsNull()) {
                conf.test.workers.emplace_back();
                continue;
            }
            if (!worker.IsMap())
                throw std::runtime_error("Invalid workers config.");
            auto &wc = conf.test.workers.emplace_back();
            if (auto node = worker["final"])
                wc.is_final = node.as<bool>();

            if (auto ins_node = worker["input"]) {
                if (!ins_node.IsSequence())
                    throw std::runtime_error("Invalid worker input config.");
                auto ins = ins_node.as<std::vector<int>>();
                if (ins.size() == 0) {
                    wc.nins = 0;
                }
                else if (ins.size() == 2) {
                    wc.nins = 2;
                    wc.ins[0] = ins[0];
                    wc.ins[1] = ins[1];
                }
                else if (ins.size() == 3) {
                    wc.nins = 3;
                    wc.ins[0] = ins[0];
                    wc.ins[1] = ins[1];
                    wc.ins[2] = ins[2];
                }
                else {
                    throw std::runtime_error("Invalid worker input config.");
                }
            }
        }

        conf.dev_filter = file["devices"];
        if (!conf.test.check())
            throw std::runtime_error("Invalid test config");
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
            res.push_back(test_device(dev, config.test));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
