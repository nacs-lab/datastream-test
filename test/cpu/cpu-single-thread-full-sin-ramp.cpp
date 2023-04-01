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
    if (diff <= 2e-5)
        return;
    if (diff <= 1e-5 * avg)
        return;
    printf("%d: expect: %f, got: %f\n", idx, expect, got);
}

template<typename Kernel>
#if NACS_CPU_X86 || NACS_CPU_X86_64
__attribute__((target_clones("default,avx,fma,avx512f")))
#endif
static void test_values(unsigned nrep, unsigned nsteps, int nparams,
                        const CPUKernel::ChnParamMod *amp_params,
                        const CPUKernel::ChnParamMod *freq_params)
{
    std::vector<float> expects[nparams];
    std::vector<float> expect2(nsteps, 0);
    for (int c = 0; c < nparams; c++) {
        expects[c] = std::vector<float>(nsteps);
        auto &expect = expects[c];
        scalar::Kernel::sin_ramp_single(&expect[0], nsteps, nrep,
                                        amp_params[c], freq_params[c]);
        for (unsigned i = 0; i < nsteps; i++) {
            expect2[i] += expect[i];
        }
    }

    Test::aligned_vector<float, 64> buff(nsteps);
    for (int c = 0; c < nparams; c++) {
        auto &expect = expects[c];
        Kernel::sin_ramp_single(&buff[0], nsteps, nrep,
                                amp_params[c], freq_params[c]);
        for (unsigned i = 0; i < nsteps; i++) {
            check_value(expect[i], buff[i], i);
        }
        Kernel::sin_ramp_single_pbuf(&buff[0], nsteps, nrep,
                                     amp_params[c], freq_params[c]);
        for (unsigned i = 0; i < nsteps; i++) {
            check_value(expect[i], buff[i], i);
        }
    }
    Kernel::sin_ramp_multi_chn_loop(&buff[0], nsteps, nrep,
                                    amp_params, freq_params, nparams);
    for (unsigned i = 0; i < nsteps; i++) {
        check_value(expect2[i], buff[i], i);
    }
    Kernel::sin_ramp_multi_chnblk_loop(&buff[0], nsteps, nrep,
                                       amp_params, freq_params, nparams);
    for (unsigned i = 0; i < nsteps; i++) {
        check_value(expect2[i], buff[i], i);
    }
    Kernel::sin_ramp_multi_block_loop(&buff[0], nsteps, nrep,
                                      amp_params, freq_params, nparams);
    for (unsigned i = 0; i < nsteps; i++) {
        check_value(expect2[i], buff[i], i);
    }
    Kernel::sin_ramp_multi_block_loop_pbuf(&buff[0], nsteps, nrep,
                                           amp_params, freq_params, nparams);
    for (unsigned i = 0; i < nsteps; i++) {
        check_value(expect2[i], buff[i], i);
    }
}

template<typename Kernel>
static YAML::Node time_run(unsigned nrep, unsigned nsteps)
{
    YAML::Node res(YAML::NodeType::Map);
    CPUKernel::ChnParamMod amp_params[max_params];
    CPUKernel::ChnParamMod freq_params[max_params];
    for (int i = 0; i < max_params; i++) {
        amp_params[i].slope = Gen::rand_single(-1, 1);
        amp_params[i].v0 = Gen::rand_single(0, 2);
        amp_params[i].v1 = Gen::rand_single(0, 2);
        amp_params[i].v2 = Gen::rand_single(0, 2);

        freq_params[i].slope = Gen::rand_single(-1, 1);
        freq_params[i].v0 = Gen::rand_single(0, 0.1f);
        freq_params[i].v1 = Gen::rand_single(0, 0.05f);
        freq_params[i].v2 = Gen::rand_single(0, 0.01f);
    }
    // Use a smaller repetition to make things go faster.
    // This is a very rough value test anyway since we aren't
    // paying much attention to the phase accumulation
    test_values<Kernel>(1, nsteps, max_params, amp_params, freq_params);
    test_values<Kernel>(3, nsteps, max_params, amp_params, freq_params);

    Test::aligned_vector<float, 64> buff(nsteps);

    Test::Timer timer;

    Kernel::sin_ramp_single(&buff[0], nsteps, 1, amp_params[0], freq_params[0]);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_ramp_single(&buff[0], nsteps, nrep, amp_params[0], freq_params[0]);
    res["single"] = timer.get_res(nrep, nsteps);

    Kernel::sin_ramp_single_pbuf(&buff[0], nsteps, 1, amp_params[0], freq_params[0]);
    timer.restart();
    if (!Test::empty)
        Kernel::sin_ramp_single_pbuf(&buff[0], nsteps, nrep,
                                     amp_params[0], freq_params[0]);
    res["single_pbuf"] = timer.get_res(nrep, nsteps);

    YAML::Node chn_loop(YAML::NodeType::Map);
    YAML::Node chnblk_loop(YAML::NodeType::Map);
    YAML::Node blk_loop(YAML::NodeType::Map);
    YAML::Node blk_loop_pbuf(YAML::NodeType::Map);

    for (int nparam: {1, 2, 5, 10, 20, 50}) {
        Kernel::sin_ramp_multi_chn_loop(&buff[0], nsteps, 1,
                                        amp_params, freq_params, nparam);
        timer.restart();
        if (!Test::empty)
            Kernel::sin_ramp_multi_chn_loop(&buff[0], nsteps, nrep / nparam,
                                            amp_params, freq_params, nparam);
        chn_loop[std::to_string(nparam)] = timer.get_res(nrep / nparam, nsteps);

        Kernel::sin_ramp_multi_chnblk_loop(&buff[0], nsteps, 1,
                                           amp_params, freq_params, nparam);
        timer.restart();
        if (!Test::empty)
            Kernel::sin_ramp_multi_chnblk_loop(&buff[0], nsteps, nrep / nparam,
                                               amp_params, freq_params, nparam);
        chnblk_loop[std::to_string(nparam)] = timer.get_res(nrep / nparam, nsteps);

        Kernel::sin_ramp_multi_block_loop(&buff[0], nsteps, 1,
                                          amp_params, freq_params, nparam);
        timer.restart();
        if (!Test::empty)
            Kernel::sin_ramp_multi_block_loop(&buff[0], nsteps, nrep / nparam,
                                              amp_params, freq_params, nparam);
        blk_loop[std::to_string(nparam)] = timer.get_res(nrep / nparam, nsteps);

        Kernel::sin_ramp_multi_block_loop_pbuf(&buff[0], nsteps, 1,
                                               amp_params, freq_params, nparam);
        timer.restart();
        if (!Test::empty)
            Kernel::sin_ramp_multi_block_loop_pbuf(&buff[0], nsteps, nrep / nparam,
                                                   amp_params, freq_params, nparam);
        blk_loop_pbuf[std::to_string(nparam)] = timer.get_res(nrep / nparam, nsteps);
    }

    res["chn_loop"] = chn_loop;
    res["chnblk_loop"] = chnblk_loop;
    res["blk_loop"] = blk_loop;
    res["blk_loop_pbuf"] = blk_loop_pbuf;

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
