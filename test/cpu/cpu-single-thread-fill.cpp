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

template<typename Kernel>
static void time_run(size_t nrep, size_t ncalc, float *buff)
{
    float t = Gen::rand_single(0, 1);
    float freq = Gen::rand_single(0, 1);
    float amp = Gen::rand_single(0, 1000);

    // Warm-up
    Test::Timer timer;
    timer.enable_cache();
    Kernel::calc_fill(1, ncalc, buff, t, freq, amp);

    timer.restart();
    if (!Test::empty)
        Kernel::calc_fill(nrep, ncalc, buff, t, freq, amp);
    timer.print(nrep, ncalc);
}

template<typename Kernel>
static void time_run_nt(size_t nrep, size_t ncalc, float *buff)
{
    float t = Gen::rand_single(0, 1);
    float freq = Gen::rand_single(0, 1);
    float amp = Gen::rand_single(0, 1000);

    // Warm-up
    Test::Timer timer;
    timer.enable_cache();
    Kernel::calc_fill_nt(1, ncalc, buff, t, freq, amp);

    timer.restart();
    if (!Test::empty)
        Kernel::calc_fill_nt(nrep, ncalc, buff, t, freq, amp);
    timer.print(nrep, ncalc);
}

static void runtests()
{
    auto &host NACS_UNUSED = CPUInfo::get_host();

    std::unique_ptr<float,void(*)(float*)> buff(
        (float*)mapAnonPage(16 * 1024 * 1024 * 4, Prot::RW),
        [] (float *ptr) { unmapPage((void*)ptr, 16 * 1024 * 1024 * 4); });

#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cout << "AVX512:" << std::endl;
        std::cout << "Normal" << std::endl;
        time_run<avx512::Kernel>(32 * 2 * 1024 * 1024, 1024, buff.get());
        time_run<avx512::Kernel>(32 * 16 * 1024, 128 * 1024, buff.get());
        time_run<avx512::Kernel>(32 * 2 * 1024, 1024 * 1024, buff.get());
        time_run<avx512::Kernel>(32 * 128, 16 * 1024 * 1024, buff.get());
        std::cout << "Non-temporal" << std::endl;
        time_run_nt<avx512::Kernel>(32 * 2 * 1024 * 1024, 1024, buff.get());
        time_run_nt<avx512::Kernel>(32 * 16 * 1024, 128 * 1024, buff.get());
        time_run_nt<avx512::Kernel>(32 * 2 * 1024, 1024 * 1024, buff.get());
        time_run_nt<avx512::Kernel>(32 * 128, 16 * 1024 * 1024, buff.get());
        return;
    }
    if (CPUKernel::hasavx2()) {
        std::cout << "AVX2:" << std::endl;
        std::cout << "Normal" << std::endl;
        time_run<avx2::Kernel>(16 * 2 * 1024 * 1024, 1024, buff.get());
        time_run<avx2::Kernel>(16 * 16 * 1024, 128 * 1024, buff.get());
        time_run<avx2::Kernel>(16 * 2 * 1024, 1024 * 1024, buff.get());
        time_run<avx2::Kernel>(16 * 128, 16 * 1024 * 1024, buff.get());
        std::cout << "Non-temporal" << std::endl;
        time_run_nt<avx2::Kernel>(16 * 2 * 1024 * 1024, 1024, buff.get());
        time_run_nt<avx2::Kernel>(16 * 16 * 1024, 128 * 1024, buff.get());
        time_run_nt<avx2::Kernel>(16 * 2 * 1024, 1024 * 1024, buff.get());
        time_run_nt<avx2::Kernel>(16 * 128, 16 * 1024 * 1024, buff.get());
        return;
    }
    if (CPUKernel::hasavx()) {
        std::cout << "AVX:" << std::endl;
        std::cout << "Normal" << std::endl;
        time_run<avx::Kernel>(12 * 2 * 1024 * 1024, 1024, buff.get());
        time_run<avx::Kernel>(12 * 16 * 1024, 128 * 1024, buff.get());
        time_run<avx::Kernel>(12 * 2 * 1024, 1024 * 1024, buff.get());
        time_run<avx::Kernel>(12 * 128, 16 * 1024 * 1024, buff.get());
        std::cout << "Non-temporal" << std::endl;
        time_run_nt<avx::Kernel>(12 * 2 * 1024 * 1024, 1024, buff.get());
        time_run_nt<avx::Kernel>(12 * 16 * 1024, 128 * 1024, buff.get());
        time_run_nt<avx::Kernel>(12 * 2 * 1024, 1024 * 1024, buff.get());
        time_run_nt<avx::Kernel>(12 * 128, 16 * 1024 * 1024, buff.get());
        return;
    }
    std::cout << "SSE2:" << std::endl;
    std::cout << "Normal" << std::endl;
    time_run<sse2::Kernel>(6 * 2 * 1024 * 1024, 1024, buff.get());
    time_run<sse2::Kernel>(6 * 16 * 1024, 128 * 1024, buff.get());
    time_run<sse2::Kernel>(6 * 2 * 1024, 1024 * 1024, buff.get());
    time_run<sse2::Kernel>(6 * 128, 16 * 1024 * 1024, buff.get());
    std::cout << "Non-temporal" << std::endl;
    time_run_nt<sse2::Kernel>(6 * 2 * 1024 * 1024, 1024, buff.get());
    time_run_nt<sse2::Kernel>(6 * 16 * 1024, 128 * 1024, buff.get());
    time_run_nt<sse2::Kernel>(6 * 2 * 1024, 1024 * 1024, buff.get());
    time_run_nt<sse2::Kernel>(6 * 128, 16 * 1024 * 1024, buff.get());
    return;
#endif

#if NACS_CPU_AARCH64
    std::cout << "ASIMD:" << std::endl;
    std::cout << "Normal" << std::endl;
    time_run<asimd::Kernel>(6 * 2 * 1024 * 1024, 1024, buff.get());
    time_run<asimd::Kernel>(6 * 16 * 1024, 128 * 1024, buff.get());
    time_run<asimd::Kernel>(6 * 2 * 1024, 1024 * 1024, buff.get());
    time_run<asimd::Kernel>(6 * 128, 16 * 1024 * 1024, buff.get());
    std::cout << "Non-temporal" << std::endl;
    time_run_nt<asimd::Kernel>(6 * 2 * 1024 * 1024, 1024, buff.get());
    time_run_nt<asimd::Kernel>(6 * 16 * 1024, 128 * 1024, buff.get());
    time_run_nt<asimd::Kernel>(6 * 2 * 1024, 1024 * 1024, buff.get());
    time_run_nt<asimd::Kernel>(6 * 128, 16 * 1024 * 1024, buff.get());
    return;
#endif

    std::cout << "Scalar:" << std::endl;
    std::cout << "Normal" << std::endl;
    time_run<scalar::Kernel>(2 * 1024 * 1024, 1024, buff.get());
    time_run<scalar::Kernel>(16 * 1024, 128 * 1024, buff.get());
    time_run<scalar::Kernel>(2 * 1024, 1024 * 1024, buff.get());
    time_run<scalar::Kernel>(128, 16 * 1024 * 1024, buff.get());
    std::cout << "Non-temporal" << std::endl;
    time_run_nt<scalar::Kernel>(2 * 1024 * 1024, 1024, buff.get());
    time_run_nt<scalar::Kernel>(16 * 1024, 128 * 1024, buff.get());
    time_run_nt<scalar::Kernel>(2 * 1024, 1024 * 1024, buff.get());
    time_run_nt<scalar::Kernel>(128, 16 * 1024 * 1024, buff.get());
}

int main()
{
    runtests();
    return 0;
}
