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
#include "print.h"

#include <iostream>

namespace Test {

static bool check_env(const char *name, bool def=false)
{
    auto env = getenv(name);
    if (!env || !*env)
        return def;
    if (strcasecmp(env, "0") == 0 || strcasecmp(env, "off") == 0 ||
        strcasecmp(env, "false") == 0 || strcasecmp(env, "f") == 0)
        return false;
    if (strcasecmp(env, "1") == 0 || strcasecmp(env, "on") == 0 ||
        strcasecmp(env, "true") == 0 || strcasecmp(env, "t") == 0)
        return true;
    return def;
}

// For background check and subtraction
NACS_EXPORT() bool empty = check_env("TEST_EMPTY");
NACS_EXPORT() bool output_json = check_env("TEST_OUTPUT_JSON");

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

NACS_EXPORT() std::map<std::string,double> Timer::get_res(size_t nrep, size_t nele)
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
    std::map<std::string,double> res;
    res["t"] = (double)timer.elapsed() / (double)nele / (double)nrep;
    res["inst"] = (double)insts.finish(false) / (double)nele / (double)nrep;
    res["cycle"] = (double)cycles.finish(false) / (double)nele / (double)nrep;
    if (m_cache_on) {
        res["cache_miss"] = (double)cachemisses.finish(false) / (double)nele / (double)nrep;
        res["cache_ref"] = (double)cacherefs.finish(false) / (double)nele / (double)nrep;
    }
    if (m_stall_on) {
        res["stall_fe"] = (double)stall_fe.finish(false) / (double)nele / (double)nrep;
        res["stall_be"] = (double)stall_be.finish(false) / (double)nele / (double)nrep;
    }
    return res;
}

NACS_EXPORT() void Timer::print(size_t nrep, size_t nele)
{
    auto res = get_res(nrep, nele);

    if (output_json) {
        Print::json(std::cout, res);
    }
    else {
        std::cout << res["t"] << " ns, " << res["inst"] << " insts, "
                  << res["cycle"] << " cycle," << std::endl;
        if (m_cache_on) {
            std::cout << "  " << res["cache_ref"] << " cache refs, "
                      << res["cache_miss"] << " cache misses,"
                      << std::endl;
        }
        if (m_stall_on) {
            std::cout << "  " << res["stall_fe"] << " frontend stall, "
                      << res["stall_be"] << " backend stall,"
                      << std::endl;
        }
    }
}

}
