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

}

#endif
