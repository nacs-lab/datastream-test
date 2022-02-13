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

#include "threads.h"

#include <pthread.h>
#include <sched.h>
#include <stdlib.h>

#include <memory>
#include <stdexcept>
#include <system_error>

namespace Thread {

using namespace NaCs;

NACS_EXPORT() std::vector<int> parse_cpulist(const char *cpulist)
{
    if (!cpulist || !*cpulist)
        return {};

    std::vector<int> cpus;

    for (const char *p = cpulist; *p;) {
        if (*p < '0' || *p > '9')
            throw std::runtime_error(std::string("Invalid cpulist: ") + cpulist);
        char *p2;
        int cpu = (int)strtol(p, &p2, 10);
        if (p2 == p)
            throw std::runtime_error(std::string("Invalid cpulist: ") + cpulist);
        if (*p2 == 0) {
            cpus.push_back(cpu);
            break;
        }
        else if (*p2 == ',') {
            cpus.push_back(cpu);
            p = p2 + 1;
            continue;
        }
        else if (*p2 != '-') {
            throw std::runtime_error(std::string("Invalid cpulist: ") + cpulist);
        }

        p = p2 + 1;
        int cpu2 = (int)strtol(p, &p2, 10);
        if (p2 == p)
            throw std::runtime_error(std::string("Invalid cpulist: ") + cpulist);
        if (*p2 == 0) {
            for (int c = cpu; c <= cpu2; c++)
                cpus.push_back(c);
            break;
        }
        else if (*p2 == ',') {
            for (int c = cpu; c <= cpu2; c++)
                cpus.push_back(c);
            p = p2 + 1;
            continue;
        }
        throw std::runtime_error(std::string("Invalid cpulist: ") + cpulist);
    }

    return cpus;
}

namespace {

struct thread_data {
    std::function<void(int)> cb;
    int id;
    int cpu;
    static void *thread_fun(void *_data)
    {
        std::unique_ptr<thread_data> data{(thread_data*)_data};
        if (data->cpu >= 0)
            NaCs::Thread::pin(data->cpu);
        data->cb(data->id);
        return nullptr;
    }
    static void start(std::function<void(int)> cb, int id, int cpu=-1)
    {
        auto data = new thread_data{cb, id, cpu};
        pthread_t thread;
        int ret = pthread_create(&thread, nullptr, thread_fun, data);
        if (ret != 0) {
            delete data;
            throw std::system_error(-ret, std::system_category(),
                                    "Failed to start thread.");
        }
    }
};

}

NACS_EXPORT() void start(int nthreads, std::function<void(int)> cb)
{
    for (int i = 0; i < nthreads; i++) {
        thread_data::start(cb, i);
    }
}

NACS_EXPORT() void start(const std::vector<int> &cpus, std::function<void(int)> cb)
{
    for (int i = 0; i < (int)cpus.size(); i++) {
        thread_data::start(cb, i, cpus[i]);
    }
}

}
