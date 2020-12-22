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

#include "helpers/cpu_kernel.h"
#include "helpers/gen_data.h"
#include "helpers/test.h"

#include <nacs-utils/mem.h>

#include <iostream>
#include <vector>

using namespace NaCs;
using namespace CPUKernel;

template<typename Kernel>
static void time_run(size_t nrep, size_t ncalc, int *buff)
{
    int v = (int)Gen::rand_single(10, 1000000);

    // Warm-up
    Test::Timer timer;
    timer.enable_cache();
    Kernel::fill(1, ncalc, buff, v);

    timer.restart();
    if (!Test::empty)
        Kernel::fill(nrep, ncalc, buff, v);
    timer.print(nrep, ncalc);
}

static void runtests()
{
    auto &host NACS_UNUSED = CPUInfo::get_host();

    std::unique_ptr<int,void(*)(int*)> buff(
        (int*)mapAnonPage(16 * 1024 * 1024 * 4, Prot::RW),
        [] (int *ptr) { unmapPage((void*)ptr, 16 * 1024 * 1024 * 4); });

#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cout << "AVX512:" << std::endl;
        time_run<avx512::Kernel>(32 * 4 * 1024 * 1024, 1024, buff.get());
        time_run<avx512::Kernel>(32 * 32 * 1024, 128 * 1024, buff.get());
        time_run<avx512::Kernel>(32 * 4 * 1024, 1024 * 1024, buff.get());
        time_run<avx512::Kernel>(32 * 256, 16 * 1024 * 1024, buff.get());
        return;
    }
    if (CPUKernel::hasavx2()) {
        std::cout << "AVX2:" << std::endl;
        time_run<avx2::Kernel>(16 * 4 * 1024 * 1024, 1024, buff.get());
        time_run<avx2::Kernel>(16 * 32 * 1024, 128 * 1024, buff.get());
        time_run<avx2::Kernel>(16 * 4 * 1024, 1024 * 1024, buff.get());
        time_run<avx2::Kernel>(16 * 256, 16 * 1024 * 1024, buff.get());
        return;
    }
    if (CPUKernel::hasavx()) {
        std::cout << "AVX:" << std::endl;
        time_run<avx::Kernel>(12 * 4 * 1024 * 1024, 1024, buff.get());
        time_run<avx::Kernel>(12 * 32 * 1024, 128 * 1024, buff.get());
        time_run<avx::Kernel>(12 * 4 * 1024, 1024 * 1024, buff.get());
        time_run<avx::Kernel>(12 * 256, 16 * 1024 * 1024, buff.get());
        return;
    }
    std::cout << "SSE2:" << std::endl;
    time_run<sse2::Kernel>(6 * 4 * 1024 * 1024, 1024, buff.get());
    time_run<sse2::Kernel>(6 * 32 * 1024, 128 * 1024, buff.get());
    time_run<sse2::Kernel>(6 * 4 * 1024, 1024 * 1024, buff.get());
    time_run<sse2::Kernel>(6 * 256, 16 * 1024 * 1024, buff.get());
    return;
#endif

#if NACS_CPU_AARCH64
    std::cout << "ASIMD:" << std::endl;
    time_run<asimd::Kernel>(6 * 4 * 1024 * 1024, 1024, buff.get());
    time_run<asimd::Kernel>(6 * 32 * 1024, 128 * 1024, buff.get());
    time_run<asimd::Kernel>(6 * 4 * 1024, 1024 * 1024, buff.get());
    time_run<asimd::Kernel>(6 * 256, 16 * 1024 * 1024, buff.get());
    return;
#endif

    std::cout << "Scalar:" << std::endl;
    time_run<scalar::Kernel>(4 * 1024 * 1024, 1024, buff.get());
    time_run<scalar::Kernel>(32 * 1024, 128 * 1024, buff.get());
    time_run<scalar::Kernel>(4 * 1024, 1024 * 1024, buff.get());
    time_run<scalar::Kernel>(256, 16 * 1024 * 1024, buff.get());
}

int main()
{
    runtests();
    return 0;
}
