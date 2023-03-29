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

#include "helpers/cpu_kernel.h"
#include "helpers/gen_data.h"
#include "helpers/test.h"

#include <nacs-utils/mem.h>

#include <iostream>
#include <vector>

using namespace NaCs;
using namespace CPUKernel;

constexpr int max_params = 20;

template<typename Kernel>
static void time_run(unsigned nrep, unsigned nsteps)
{
    CPUKernel::ChnParamFixed params[max_params];
    for (int i = 0; i < max_params; i++) {
        params[i].freq = Gen::rand_single(0, 1);
        params[i].amp = Gen::rand_single(0, 1000);
    }

    std::vector<float> buff(nsteps);

    // Warm-up
    Test::Timer timer;

    std::cout << "single:" << std::endl;
    Kernel::sin_single(&buff[0], nsteps, 1, params[0]);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_single(&buff[0], nsteps, nrep, params[0]);
    timer.print(nrep, nsteps);

    std::cout << "multi-chn-loop 1:" << std::endl;
    Kernel::sin_multi_chn_loop(&buff[0], nsteps, 1, params, 1);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_chn_loop(&buff[0], nsteps, nrep, params, 1);
    timer.print(nrep, nsteps);

    std::cout << "multi-block-loop 1:" << std::endl;
    Kernel::sin_multi_block_loop(&buff[0], nsteps, 1, params, 1);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_block_loop(&buff[0], nsteps, nrep, params, 1);
    timer.print(nrep, nsteps);

    std::cout << "multi-chn-loop 2:" << std::endl;
    Kernel::sin_multi_chn_loop(&buff[0], nsteps, 1, params, 2);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_chn_loop(&buff[0], nsteps, nrep / 2, params, 2);
    timer.print(nrep / 2, nsteps);

    std::cout << "multi-block-loop 2:" << std::endl;
    Kernel::sin_multi_block_loop(&buff[0], nsteps, 1, params, 2);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_block_loop(&buff[0], nsteps, nrep / 2, params, 2);
    timer.print(nrep / 2, nsteps);

    std::cout << "multi-chn-loop 5:" << std::endl;
    Kernel::sin_multi_chn_loop(&buff[0], nsteps, 1, params, 5);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_chn_loop(&buff[0], nsteps, nrep / 4, params, 5);
    timer.print(nrep / 4, nsteps);

    std::cout << "multi-block-loop 5:" << std::endl;
    Kernel::sin_multi_block_loop(&buff[0], nsteps, 1, params, 5);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_block_loop(&buff[0], nsteps, nrep / 4, params, 5);
    timer.print(nrep / 4, nsteps);

    std::cout << "multi-chn-loop 10:" << std::endl;
    Kernel::sin_multi_chn_loop(&buff[0], nsteps, 1, params, 10);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_chn_loop(&buff[0], nsteps, nrep / 8, params, 10);
    timer.print(nrep / 8, nsteps);

    std::cout << "multi-block-loop 10:" << std::endl;
    Kernel::sin_multi_block_loop(&buff[0], nsteps, 1, params, 10);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_block_loop(&buff[0], nsteps, nrep / 8, params, 10);
    timer.print(nrep / 8, nsteps);

    std::cout << "multi-chn-loop 20:" << std::endl;
    Kernel::sin_multi_chn_loop(&buff[0], nsteps, 1, params, 20);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_chn_loop(&buff[0], nsteps, nrep / 16, params, 20);
    timer.print(nrep / 16, nsteps);

    std::cout << "multi-block-loop 20:" << std::endl;
    Kernel::sin_multi_block_loop(&buff[0], nsteps, 1, params, 20);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_multi_block_loop(&buff[0], nsteps, nrep / 16, params, 20);
    timer.print(nrep / 16, nsteps);
}

static void runtests()
{
#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cout << "AVX512:" << std::endl;
        time_run<avx512::Kernel>(8 * 1024 * 1024, 2 * 1024);
    }
    if (CPUKernel::hasavx2()) {
        std::cout << "AVX2:" << std::endl;
        time_run<avx2::Kernel>(4 * 1024 * 1024, 2 * 1024);
    }
    if (CPUKernel::hasavx()) {
        std::cout << "AVX:" << std::endl;
        time_run<avx::Kernel>(4 * 1024 * 1024, 2 * 1024);
    }
    std::cout << "SSE2:" << std::endl;
    time_run<sse2::Kernel>(2 * 1024 * 1024, 2 * 1024);
#endif

#if NACS_CPU_AARCH64
    if (CPUKernel::hassve()) {
        std::cout << "SVE:" << std::endl;
        time_run<sve::Kernel>(2 * 1024 * 1024, 2 * 1024);
    }
    std::cout << "ASIMD:" << std::endl;
    time_run<asimd::Kernel>(2 * 1024 * 1024, 2 * 1024);
#endif

    std::cout << "Scalar:" << std::endl;
    time_run<scalar::Kernel>(512 * 1024, 2 * 1024);
}

int main()
{
    runtests();
    return 0;
}
