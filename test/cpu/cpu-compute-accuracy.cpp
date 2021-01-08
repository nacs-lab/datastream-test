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
#include "helpers/test.h"

#include <yaml-cpp/yaml.h>

#include <stdlib.h>

#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace CPUKernel;

#if NACS_CPU_X86 || NACS_CPU_X86_64
__attribute__((target_clones("default,avx,fma,avx512f")))
#endif
static void compute_baseline(double *out, double start, double step, unsigned nsteps)
{
#pragma omp parallel for simd
    for (unsigned i = 0; i < nsteps; i++) {
        out[i] = sin(M_PI * (start + step * i)) / M_PI;
    }
}

#if NACS_CPU_X86 || NACS_CPU_X86_64
__attribute__((target_clones("default,avx,fma,avx512f")))
#endif
static void compute_difference(int32_t *diff, const float *res,
                               const double *baseline, unsigned nsteps)
{
#pragma omp parallel for simd
    for (unsigned i = 0; i < nsteps; i++) {
        diff[i] = (int32_t)round((res[i] - baseline[i]) * M_PI /
                                 std::numeric_limits<float>::epsilon());
    }
}

template<typename Kernel>
static YAML::Node test_kernel(double start, double step, unsigned nsteps, const double *baseline)
{
    std::unique_ptr<float> buff{(float*)aligned_alloc(64, size_t(nsteps) * sizeof(float))};
    std::vector<int32_t> diff(nsteps);
    Kernel::sin_range(buff.get(), start, step, nsteps);
    compute_difference(&diff[0], buff.get(), baseline, nsteps);

    YAML::Node res(YAML::NodeType::Map);
    res["start"] = start;
    res["step"] = step;
    res["nsteps"] = nsteps;
    res["diff"] = diff;
    return res;
}

static YAML::Node runtests(double start, double step, unsigned nsteps, const double *baseline)
{
    YAML::Node res(YAML::NodeType::Map);
#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        res["avx512"] = test_kernel<avx512::Kernel>(start, step, nsteps, baseline);
    }
    if (CPUKernel::hasavx2()) {
        res["avx2"] = test_kernel<avx2::Kernel>(start, step, nsteps, baseline);
    }
    if (CPUKernel::hasavx()) {
        res["avx"] = test_kernel<avx::Kernel>(start, step, nsteps, baseline);
    }
    res["sse2"] = test_kernel<sse2::Kernel>(start, step, nsteps, baseline);
#endif
#if NACS_CPU_AARCH64
    res["asimd"] = test_kernel<asimd::Kernel>(start, step, nsteps, baseline);
#endif
    res["scalar"] = test_kernel<scalar::Kernel>(start, step, nsteps, baseline);
    return res;
}

struct Config {
    double start;
    double step;
    unsigned nsteps;
    static Config loadYAML(const char *fname)
    {
        Config conf;
        auto file = YAML::LoadFile(fname);
        auto required_key = [&] (auto name) {
            if (auto node = file[name])
                return node;
            throw std::runtime_error(std::string("Required key '") + name + "' missing.");
        };
        conf.start = required_key("start").as<double>();
        conf.step = required_key("step").as<double>();
        conf.nsteps = required_key("nsteps").as<unsigned>();
        return conf;
    }
};

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Needs at least one argument\n");
        exit(1);
    }
    auto config = Config::loadYAML(argv[1]);

    std::vector<double> baseline(config.nsteps);
    compute_baseline(&baseline[0], config.start, config.step, config.nsteps);

    auto res = runtests(config.start, config.step, config.nsteps, &baseline[0]);
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
