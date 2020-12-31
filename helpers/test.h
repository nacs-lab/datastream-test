/*************************************************************************
 *   Copyright (c) 2019 - 2020 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef HELPERS_TEST_H
#define HELPERS_TEST_H

#include <nacs-utils/utils.h>
#include <nacs-utils/timer.h>

#include <map>
#include <string>

namespace Test {

using namespace NaCs;

// For background check and subtraction
NACS_EXPORT(ds_helper) extern bool empty;
NACS_EXPORT(ds_helper) extern bool output_json;

struct Timer {
    NaCs::Timer timer;
    PerfCounter insts{PerfCounter::CPUInsts};
    PerfCounter cycles{PerfCounter::CPUCycles};
    PerfCounter cacherefs{PerfCounter::CacheRefs};
    PerfCounter cachemisses{PerfCounter::CacheMisses};
    PerfCounter stall_fe{PerfCounter::CPUStallFrontend};
    PerfCounter stall_be{PerfCounter::CPUStallBackend};

    NACS_EXPORT(ds_helper) Timer();
    NACS_EXPORT(ds_helper) void enable_cache(bool on=true);
    NACS_EXPORT(ds_helper) void enable_stall(bool on=true);

    NACS_EXPORT(ds_helper) void restart();
    NACS_EXPORT(ds_helper) std::map<std::string,double> get_res(size_t nrep, size_t nele);
    NACS_EXPORT(ds_helper) void print(size_t nrep=1, size_t nele=1);

private:
    bool m_cache_on{false};
    bool m_stall_on{false};
};

// number of cycles since power-on
static NACS_INLINE uint64_t cycleclock()
{
#if NACS_CPU_X86_64
    uint64_t low, high;
    __asm__ volatile("rdtsc" : "=a"(low), "=d"(high));
    return (high << 32) | low;
#elif NACS_CPU_X86
    unt64_t ret;
    __asm__ volatile("rdtsc" : "=A"(ret));
    return ret;
#elif NACS_CPU_AARCH64
    // System timer of ARMv8 runs at a different frequency than the CPU's.
    // The frequency is fixed, typically in the range 1-50MHz.  It can be
    // read at CNTFRQ special register.  We assume the OS has set up
    // the virtual timer properly.
    uint64_t virtual_timer_value;
    __asm__ volatile("mrs %0, cntvct_el0" : "=r"(virtual_timer_value));
    return virtual_timer_value;
#elif NACS_CPU_PPC64
    // This returns a time-base, which is not always precisely a cycle-count.
    // https://reviews.llvm.org/D78084
    uint64_t tb;
    asm volatile("mfspr %0, 268" : "=r" (tb));
    return tb;
#else
    #warning No cycleclock() definition for your platform
    // copy from https://github.com/google/benchmark/blob/v1.5.0/src/cycleclock.h
    return 0;
#endif
}

}

#endif
