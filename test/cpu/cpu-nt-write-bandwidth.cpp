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
#include "helpers/threads.h"

#include <nacs-utils/mem.h>
#include <nacs-utils/number.h>

#include <atomic>
#include <iostream>
#include <vector>

using namespace NaCs;
using namespace CPUKernel;

template<typename Kernel>
static void time_run_real(size_t nrep, size_t ncalc, int *buff, std::atomic<bool> &done)
{
    int v = (int)Gen::rand_single(10, 1000000);

    // Warm-up
    Test::Timer timer;
    timer.enable_cache();
    Kernel::fill_nt(1, ncalc, buff, v);

    timer.restart();
    if (!Test::empty)
        Kernel::fill_nt(nrep, ncalc, buff, v);
    timer.print(nrep, ncalc);
    done.store(true, std::memory_order_release);
}

template<typename Kernel>
static void time_run_dummy(size_t nrep, size_t ncalc,
                           std::atomic<int> &ready, const std::atomic<bool> &done)
{
    int v = (int)Gen::rand_single(10, 1000000);
    nrep /= 16;
    if (nrep < 1)
        nrep = 1;
    auto buff = (int*)mapAnonPage(ncalc * 4, Prot::RW);
    Kernel::fill_nt(nrep, ncalc, buff, v);
    // Signal that we are ready.
    ready.fetch_add(1, std::memory_order_release);
    CPU::wake();
    do {
        Kernel::fill_nt(nrep, ncalc, buff, v);
    } while (!done.load(std::memory_order_acquire));
    unmapPage((void*)buff, ncalc * 4);
    ready.fetch_add(-1, std::memory_order_release);
    CPU::wake();
}

template<typename Kernel>
static void time_run(size_t nrep, size_t ncalc, int *buff, const std::vector<int> &cpus)
{
    std::atomic<bool> done{false};
    std::atomic<int> ready{0};
    Thread::start(cpus, [=, &done, &ready] (int) {
        time_run_dummy<Kernel>(nrep, ncalc, ready, done);
    });
    while (ready.load(std::memory_order_acquire) != cpus.size()) {
        CPU::pause();
    }
    time_run_real<Kernel>(nrep, ncalc, buff, done);
    while (ready.load(std::memory_order_acquire) != 0) {
        CPU::pause();
    }
}

static void runtests(const std::vector<int> &cpus, long _size)
{
    size_t size = alignTo(_size, page_size);
    auto freepage = [=] (int *ptr) { unmapPage((void*)ptr, size); };
    std::unique_ptr<int,decltype(freepage)> buff((int*)mapAnonPage(size, Prot::RW), freepage);

#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cout << "AVX512:" << std::endl;
        time_run<avx512::Kernel>(size_t(32 * 16ull * 1024 * 1024 * 1024 / size),
                                 size / 4, buff.get(), cpus);
        return;
    }
    if (CPUKernel::hasavx2()) {
        std::cout << "AVX2:" << std::endl;
        time_run<avx2::Kernel>(size_t(16 * 16ull * 1024 * 1024 * 1024 / size),
                               size / 4, buff.get(), cpus);
        return;
    }
    if (CPUKernel::hasavx()) {
        std::cout << "AVX:" << std::endl;
        time_run<avx::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                              size / 4, buff.get(), cpus);
        return;
    }
    std::cout << "SSE2:" << std::endl;
    time_run<sse2::Kernel>(size_t(6 * 16ull * 1024 * 1024 * 1024 / size),
                           size / 4, buff.get(), cpus);
    return;
#endif

#if NACS_CPU_AARCH64
    std::cout << "ASIMD:" << std::endl;
    time_run<asimd::Kernel>(size_t(6 * 16ull * 1024 * 1024 * 1024 / size),
                            size / 4, buff.get(), cpus);
    return;
#endif

    std::cout << "Scalar:" << std::endl;
    time_run<scalar::Kernel>(size_t(16ull * 1024 * 1024 * 1024 / size),
                             size / 4, buff.get(), cpus);
}

static inline long parse_int(const char *s)
{
    if (!*s) {
        fprintf(stderr, "Error parsing integer '%s'\n", s);
        exit(1);
    }
    char *ep;
    long res = strtol(s, &ep, 10);
    if (s == ep || *ep != 0) {
        fprintf(stderr, "Error parsing integer '%s'\n", s);
        exit(1);
    }
    return res;
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "Needs at least three arguments\n");
        exit(1);
    }
    Thread::pin((int)parse_int(argv[1]));
    runtests(Thread::parse_cpulist(argv[2]), parse_int(argv[3]));
    return 0;
}
