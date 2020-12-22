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

#include "test.h"

#include <iostream>

namespace Test {

// For background check and subtraction
NACS_EXPORT() bool empty = getenv("TEST_EMPTY");

NACS_EXPORT() Timer::Timer()
{
}

NACS_EXPORT() void Timer::enable_cache(bool on)
{
    m_cache_on = on;
}

NACS_EXPORT() void Timer::enable_stall(bool on)
{
    m_stall_on = on;
}

NACS_EXPORT() void Timer::restart()
{
    insts.reset();
    cycles.reset();
    if (m_cache_on) {
        cachemisses.reset();
        cacherefs.reset();
    }
    if (m_stall_on) {
        stall_fe.reset();
        stall_be.reset();
    }

    timer.restart();
    if (m_stall_on) {
        stall_be.start(false);
        stall_fe.start(false);
    }
    if (m_cache_on) {
        cacherefs.start(false);
        cachemisses.start(false);
    }
    cycles.start(false);
    insts.start(false);
}

NACS_EXPORT() void Timer::print(size_t nrep, size_t ncalc)
{
    insts.stop();
    cycles.stop();
    if (m_cache_on) {
        cachemisses.stop();
        cacherefs.stop();
    }
    if (m_stall_on) {
        stall_fe.stop();
        stall_be.stop();
    }
    auto tdry = (double)timer.elapsed() / (double)ncalc / (double)nrep;
    auto ninsts = (double)insts.finish(false) / (double)ncalc / (double)nrep;
    auto ncycles = (double)cycles.finish(false) / (double)ncalc / (double)nrep;
    double ncachemisses = 0;
    double ncacherefs = 0;
    double nstall_fe = 0;
    double nstall_be = 0;
    if (m_cache_on) {
        ncachemisses = (double)cachemisses.finish(false) / (double)ncalc / (double)nrep;
        ncacherefs = (double)cacherefs.finish(false) / (double)ncalc / (double)nrep;
    }
    if (m_stall_on) {
        nstall_fe = (double)stall_fe.finish(false) / (double)ncalc / (double)nrep;
        nstall_be = (double)stall_be.finish(false) / (double)ncalc / (double)nrep;
    }

    std::cout << tdry << " ns, " << ninsts << " insts, "
              << ncycles << " cycle," << std::endl;
    if (m_cache_on) {
        std::cout << "  " << ncacherefs << " cache refs, " << ncachemisses << " cache misses,"
                  << std::endl;
    }
    if (m_stall_on) {
        std::cout << "  " << nstall_fe << " frontend stall, " << nstall_be << " backend stall,"
                  << std::endl;
    }
}

}
