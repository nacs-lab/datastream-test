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

#ifndef HELPERS_GEN_DATA_H
#define HELPERS_GEN_DATA_H

#include <random>

namespace Gen {

namespace {
static std::random_device rd; // Will be used to obtain a seed for the random number engine
static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

static inline float rand_single(float lb=0, float ub=1)
{
    std::uniform_real_distribution<float> dis(lb, ub);
    return dis(gen);
}

template<typename T>
static inline void rand_fill(T &out, float lb=0, float ub=1)
{
    std::uniform_real_distribution<float> dis(lb, ub);
    for (auto &d: out) {
        d = dis(gen);
    }
}

}

}

#endif
