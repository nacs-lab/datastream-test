/*************************************************************************
 *   Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "helpers/cpu_kernel.h"
#include "helpers/gen_data.h"

#include <nacs-utils/processor.h>
#include <nacs-utils/timer.h>
#include <nacs-utils/mem.h>

#include <iostream>
#include <vector>

using namespace NaCs;
using namespace CPUKernel;

// For background subtraction
static bool empty_run = getenv("TEST_EMPTY");

template<typename Kernel>
static void time_run(size_t nrep, size_t ncalc)
{
    float t = Gen::rand_single(0, 1);
    float freq = Gen::rand_single(0, 1);
    float amp = Gen::rand_single(0, 1000);

    // Warm-up
    Timer timer;
    Kernel::calc_dry(1, 1, t, freq, amp);

    timer.restart();
    if (!empty_run)
        Kernel::calc_dry(nrep, ncalc, t, freq, amp);
    auto tdry = (double)timer.elapsed() / (double)ncalc / (double)nrep;

    std::cout << "Dry: " << tdry << " ns / ele" << std::endl;
}

static void runtests()
{
    auto &host NACS_UNUSED = CPUInfo::get_host();

#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (host.test_feature(X86::Feature::avx512f) &&
        host.test_feature(X86::Feature::avx512dq)) {
        std::cout << "AVX512:" << std::endl;
        time_run<avx512::Kernel>(32 * 1024, 1024 * 1024);
        return;
    }
    if (host.test_feature(X86::Feature::avx2) && host.test_feature(X86::Feature::fma)) {
        std::cout << "AVX2:" << std::endl;
        time_run<avx2::Kernel>(16 * 1024, 1024 * 1024);
        return;
    }
    if (host.test_feature(X86::Feature::avx)) {
        std::cout << "AVX:" << std::endl;
        time_run<avx::Kernel>(12 * 1024, 1024 * 1024);
        return;
    }
    std::cout << "SSE2:" << std::endl;
    time_run<sse2::Kernel>(6 * 1024, 1024 * 1024);
    return;
#endif

    std::cout << "Scalar:" << std::endl;
    time_run<scalar::Kernel>(2 * 1024, 1024 * 1024);
}

int main()
{
    runtests();
    return 0;
}
