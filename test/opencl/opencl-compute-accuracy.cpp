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

#include "helpers/test.h"
#include "helpers/opencl.h"

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#if NACS_CPU_X86 || NACS_CPU_X86_64
__attribute__((target_clones("default,avx,fma,avx512f")))
#endif
static void compute_baseline(double *out, double start, double step, unsigned nsteps)
{
#pragma omp parallel for simd
    for (unsigned i = 0; i < nsteps; i++) {
        out[i] = sin(start + step * i);
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
        diff[i] = (int32_t)round((res[i] - baseline[i]) /
                                 std::numeric_limits<float>::epsilon());
    }
}

static YAML::Node test_device(cl::Device &dev, double start, double step, unsigned nsteps,
                              const double *baseline)
{
    cl::Context ctx({dev});
    cl::CommandQueue queue(ctx, dev);
    std::string source{R"CLC(
        kernel void compute_sin(global float *res, double start, double step)
        {
            res[get_global_id(0)] = sin(start + get_global_id(0) * step);
        }

        kernel void compute_sin_native(global float *res, double start, double step)
        {
            res[get_global_id(0)] = native_sin(start + get_global_id(0) * step);
        }
    )CLC"};
    cl::Program prog(ctx, {source}, true);
    std::vector<float> hostbuff(nsteps);
    cl::Buffer buff(ctx, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
                    sizeof(float) * nsteps, &hostbuff[0]);

    std::vector<int32_t> diff(nsteps);
    {
        cl::Kernel kernel(prog, "compute_sin");
        kernel.setArg(0, buff);
        kernel.setArg(1, start);
        kernel.setArg(2, step);
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nsteps));
        queue.enqueueMapBuffer(buff, true, CL_MAP_READ, 0, sizeof(float) * nsteps);
        compute_difference(&diff[0], &hostbuff[0], baseline, nsteps);
    }

    std::vector<int32_t> diff_native(nsteps);
    {
        cl::Kernel kernel(prog, "compute_sin_native");
        kernel.setArg(0, buff);
        kernel.setArg(1, start);
        kernel.setArg(2, step);
        queue.enqueueNDRangeKernel(kernel, cl::NDRange(), cl::NDRange(nsteps));
        queue.enqueueMapBuffer(buff, true, CL_MAP_READ, 0, sizeof(float) * nsteps);
        compute_difference(&diff_native[0], &hostbuff[0], baseline, nsteps);
    }

    auto res = OCL::get_device_ids(dev);
    res["start"] = start;
    res["step"] = step;
    res["nsteps"] = nsteps;
    res["diff"] = diff;
    res["diff_native"] = diff_native;
    return res;
}

struct Config {
    YAML::Node dev_filter;
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
        conf.dev_filter = file["devices"];
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
    std::vector<cl::Device> devices = OCL::all_ocl2_devices(&config.dev_filter);
    if (devices.empty())
        throw std::runtime_error("Unable to find OpenCL devices");

    std::vector<double> baseline(config.nsteps);
    compute_baseline(&baseline[0], config.start, config.step, config.nsteps);

    std::vector<YAML::Node> res;
    for (auto &dev: devices) {
        OCL::catch_error([&] {
            res.push_back(test_device(dev, config.start, config.step,
                                      config.nsteps, &baseline[0]));
        });
    }
    YAML::Emitter out;
    out << res;
    std::cout << out.c_str() << std::endl;
    return 0;
}
