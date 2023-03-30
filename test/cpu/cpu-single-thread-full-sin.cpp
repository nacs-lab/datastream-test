/*************************************************************************
 *   Copyright (c) 2023 - 2023 Yichao Yu <yyc1992@gmail.com>             *
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

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <vector>

using namespace NaCs;
using namespace CPUKernel;

constexpr int max_params = 50;

static inline void check_value(float expect, float got, unsigned idx)
{
    auto diff = abs(expect - got);
    auto avg = abs(expect + got) / 2;
    if (diff <= 2e-3)
        return;
    if (diff <= 1e-3 * avg)
        return;
    printf("%d: expect: %f, got: %f\n", idx, expect, got);
}

template<typename Kernel>
#if NACS_CPU_X86 || NACS_CPU_X86_64
__attribute__((target_clones("default,avx,fma,avx512f")))
#endif
static void test_values(unsigned nrep, unsigned nsteps, int nparams,
                        const CPUKernel::ChnParamFixed *params)
{
    std::vector<float> expects[nparams];
    for (int c = 0; c < nparams; c++) {
        expects[c] = std::vector<float>(nsteps);
        auto &expect = expects[c];
        auto param = params[c];
        double phase = param.freq * nsteps * (nrep - 1);
        while (phase > 32) {
            phase -= 64;
        }
#pragma omp parallel for simd
        for (unsigned i = 0; i < nsteps; i++) {
            expect[i] = float(param.amp * sin(M_PI * (phase + i * param.freq)) / M_PI);
        }
    }

    std::vector<float> buff(nsteps);
    std::vector<float> expect2(nsteps, 0);
    for (int c = 0; c < nparams; c++) {
        auto &expect = expects[c];
        Kernel::sin_single(&buff[0], nsteps, nrep, params[c]);
        for (unsigned i = 0; i < nsteps; i++) {
            expect2[i] += buff[i];
            check_value(expect[i], buff[i], i);
        }
    }
    Kernel::sin_multi_chn_loop(&buff[0], nsteps, nrep, params, nparams);
    for (unsigned i = 0; i < nsteps; i++) {
        check_value(expect2[i], buff[i], i);
    }
    Kernel::sin_multi_block_loop(&buff[0], nsteps, nrep, params, nparams);
    for (unsigned i = 0; i < nsteps; i++) {
        check_value(expect2[i], buff[i], i);
    }
}

template<typename Kernel>
static YAML::Node time_run(unsigned nrep, unsigned nsteps)
{
    YAML::Node res(YAML::NodeType::Map);
    CPUKernel::ChnParamFixed params[max_params];
    for (int i = 0; i < max_params; i++) {
        params[i].freq = Gen::rand_single(0, 1);
        params[i].amp = Gen::rand_single(0, 1000);
    }
    // Use a smaller repetition to make things go faster.
    // This is a very rough value test anyway since we aren't
    // paying much attention to the phase accumulation
    test_values<Kernel>(1, nsteps, max_params, params);
    test_values<Kernel>(3, nsteps, max_params, params);

    std::vector<float> buff(nsteps);

    Test::Timer timer;

    Kernel::sin_single(&buff[0], nsteps, 1, params[0]);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_single(&buff[0], nsteps, nrep, params[0]);
    res["single"] = timer.get_res(nrep, nsteps);

    YAML::Node chn_loop(YAML::NodeType::Map);
    YAML::Node blk_loop(YAML::NodeType::Map);

    for (int nparam: {1, 2, 5, 10, 20, 50}) {
        Kernel::sin_multi_chn_loop(&buff[0], nsteps, 1, params, nparam);
        timer.restart();
        if (!Test::empty)
            Kernel::sin_multi_chn_loop(&buff[0], nsteps, nrep / nparam,
                                       params, nparam);
        chn_loop[std::to_string(nparam)] = timer.get_res(nrep / nparam, nsteps);

        Kernel::sin_multi_block_loop(&buff[0], nsteps, 1, params, nparam);
        timer.restart();
        if (!Test::empty)
            Kernel::sin_multi_block_loop(&buff[0], nsteps, nrep / nparam,
                                         params, nparam);
        blk_loop[std::to_string(nparam)] = timer.get_res(nrep / nparam, nsteps);
    }

    res["chn_loop"] = chn_loop;
    res["blk_loop"] = blk_loop;

    return res;
}

static YAML::Node runtests()
{
    YAML::Node res(YAML::NodeType::Map);

#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cerr << "AVX512:" << std::endl;
        auto r = time_run<avx512::Kernel>(8 * 1024 * 1024, 2 * 1024);
        std::cerr << r << std::endl;
        res["avx512"] = r;
    }
    if (CPUKernel::hasavx2()) {
        std::cerr << "AVX2:" << std::endl;
        auto r = time_run<avx2::Kernel>(4 * 1024 * 1024, 2 * 1024);
        std::cerr << r << std::endl;
        res["avx2"] = r;
    }
    if (CPUKernel::hasavx()) {
        std::cerr << "AVX:" << std::endl;
        auto r = time_run<avx::Kernel>(4 * 1024 * 1024, 2 * 1024);
        std::cerr << r << std::endl;
        res["avx"] = r;
    }
    {
        std::cerr << "SSE2:" << std::endl;
        auto r = time_run<sse2::Kernel>(2 * 1024 * 1024, 2 * 1024);
        std::cerr << r << std::endl;
        res["sse2"] = r;
    }
#endif

#if NACS_CPU_AARCH64
    if (CPUKernel::hassve()) {
        std::cerr << "SVE:" << std::endl;
        auto r = time_run<sve::Kernel>(2 * 1024 * 1024, 2 * 1024);
        std::cerr << r << std::endl;
        res["sve"] = r;
    }
    {
        std::cerr << "ASIMD:" << std::endl;
        auto r = time_run<asimd::Kernel>(2 * 1024 * 1024, 2 * 1024);
        std::cerr << r << std::endl;
        res["asimd"] = r;
    }
#endif

    {
        std::cerr << "Scalar:" << std::endl;
        auto r = time_run<scalar::Kernel>(512 * 1024, 2 * 1024);
        std::cerr << r << std::endl;
        res["scalar"] = r;
    }

    return res;
}

int main()
{
    std::cout << runtests() << std::endl;
    return 0;
}
