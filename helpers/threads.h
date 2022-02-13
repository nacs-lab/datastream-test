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

#ifndef HELPERS_THREAD_H
#define HELPERS_THREAD_H

#include <nacs-utils/utils.h>
#include <nacs-utils/thread.h>

#include <functional>
#include <vector>

namespace Thread {

NACS_EXPORT(ds_helper) std::vector<int> parse_cpulist(const char *cpulist);

NACS_EXPORT(ds_helper) void start(int nthreads, std::function<void(int)> cb);
NACS_EXPORT(ds_helper) void start(const std::vector<int> &cpus, std::function<void(int)> cb);

}

#endif
