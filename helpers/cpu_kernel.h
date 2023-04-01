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

#ifndef HELPERS_CPU_KERNEL_H
#define HELPERS_CPU_KERNEL_H

#include <nacs-utils/utils.h>
#include <nacs-utils/number.h>
#include <nacs-utils/processor.h>

#if NACS_CPU_X86 || NACS_CPU_X86_64
#  include <immintrin.h>
#elif NACS_CPU_AARCH64
#  include <arm_neon.h>
#  include <arm_sve.h>
#endif

namespace CPUKernel {

NACS_EXPORT(ds_helper) extern const NaCs::CPUInfo &host;
NACS_EXPORT(ds_helper) bool hasavx();
NACS_EXPORT(ds_helper) bool hasavx2();
NACS_EXPORT(ds_helper) bool hasavx512();
NACS_EXPORT(ds_helper) bool hassve();

struct ChnParamFixed {
    double freq;
    double amp;
};

struct ChnParamMod {
    double slope;
    double v0;
    double v1;
    double v2;
};

namespace scalar {

struct Kernel {
    // Calculate `sin(pi * d) / pi`
    // By computing `sin(pi * d)` instead of `sin(d)`, we don't need to do
    // additional scaling of `d` before converting it to integer.
    // Then by computing `sin(pi * d) / pi` the first order
    // expansion term is kept as `d` which saves another multiplication at the end.
    // We still need to generate the correct output from the correct input in the end
    // so we still need the correct input/output scaling.
    // However, the input needs to be scaled anyway so we can fold the input scaling in there,
    // for the output, we can scale it once after all the channels are computed instead
    // of doing it once per channel.
    static NACS_INLINE float sinpif_pi(float d)
    {
        int q = NaCs::round<int>(d);
        // Now `d` is the fractional part in the range `[-0.5, 0.5]`
        d = d - (float)q;
        auto s = d * d;

        // For the original `d` in the range (0.5, 1.5), (2.5, 3.5) etc
        // their value is the same as (-0.5, 0.5) with a sign flip in the input.
        if (q & 1)
            d = -d;

        // These coefficients are numerically optimized to
        // give the smallest maximum error over [0, 4] / [-4, 4]
        // The maximum error is ~4.08e-7.
        // When only limited to [0, 0.5], the error could be reduced to ~2.9e-7
        // Either should be good enough since the output only has 16bit resolution.
        auto u = -0.17818783f * s + 0.8098674f;
        u = u * s - 1.6448531f;
        return (s * d) * u + d;
    }
    static NACS_INLINE double ramp_func(double x, ChnParamMod param)
    {
        x = x * param.slope;
        x = x - NaCs::round<int>(x);
        return param.v0 + x * (param.v1 + x * param.v2);
    }
    NACS_EXPORT(ds_helper) static
    void sin_range(float *out, double start, double step, unsigned nsteps);
    NACS_EXPORT(ds_helper) static
    void calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static void fill(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void fill_nt(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void copy(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void copy_nt(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void fill1(size_t nele, int *buff, int v);
    NACS_EXPORT(ds_helper) static void read1(size_t nele, const int *buff);
    NACS_EXPORT(ds_helper) static void sum(size_t nele, const float *buff1, const float *buff2);

    NACS_EXPORT(ds_helper) static
    void calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_fill_nt(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void read_calc_multi(size_t nele, int nchn, const float *buff,
                         float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                    float *out, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                       float *out, float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void sum_multi(size_t nele, int nins, const float *ins[], float *out);
    NACS_EXPORT(ds_helper) static
    void sum_multi_nt(size_t nele, int nins, const float *ins[], float *out);

    NACS_EXPORT(ds_helper) static
    void sin_single(float *out, unsigned nsteps, unsigned nrep, ChnParamFixed param);
    NACS_EXPORT(ds_helper) static
    void sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                            const ChnParamFixed *params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                              const ChnParamFixed *params, unsigned nparams);

    NACS_EXPORT(ds_helper) static
    void sin_ramp_single(float *out, unsigned nsteps, unsigned nrep,
                         ChnParamMod amp_param, ChnParamMod freq_param);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                 const ChnParamMod *amp_params,
                                 const ChnParamMod *freq_params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_multi_chnblk_loop(float *out, unsigned nsteps, unsigned nrep,
                                    const ChnParamMod *amp_params,
                                    const ChnParamMod *freq_params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                   const ChnParamMod *amp_params,
                                   const ChnParamMod *freq_params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_single_pbuf(float *out, unsigned nsteps, unsigned nrep,
                              ChnParamMod amp_param, ChnParamMod freq_param);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_multi_block_loop_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                        const ChnParamMod *amp_params,
                                        const ChnParamMod *freq_params,
                                        unsigned nparams);
};

} // namespace scalar

#if NACS_CPU_AARCH64

namespace asimd {
struct Kernel {
    static NACS_INLINE float32x4_t sinpif_pi(float32x4_t d)
    {
        auto q = vcvtnq_s32_f32(d);
        // Now `d` is the fractional part in the range `[-0.5, 0.5]`
        d = d - vcvtq_f32_s32(q);
        auto s = d * d;

        // Shift the last bit of `q` to the sign bit
        // and therefore flip the sign of `d` if `q` is odd
        d = float32x4_t(vshlq_n_s32(q, 31) ^ int32x4_t(d));

        auto u = -0.17818783f * s + 0.8098674f;
        u = u * s - 1.6448531f;
        return (s * d) * u + d;
    }
    static NACS_INLINE float64x2_t ramp_func(float64x2_t x, ChnParamMod param)
    {
        x = x * param.slope;
        x = x - vcvtq_f64_s64(vcvtnq_s64_f64(x));
        return param.v0 + x * (param.v1 + x * param.v2);
    }
    NACS_EXPORT(ds_helper) static
    void sin_range(float *out, double start, double step, unsigned nsteps);
    NACS_EXPORT(ds_helper) static
    void calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static void fill(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void fill_nt(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void copy(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void copy_nt(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void fill1(size_t nele, int *buff, int v);
    NACS_EXPORT(ds_helper) static void read1(size_t nele, const int *buff);
    NACS_EXPORT(ds_helper) static void sum(size_t nele, const float *buff1, const float *buff2);

    NACS_EXPORT(ds_helper) static
    void calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_fill_nt(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void read_calc_multi(size_t nele, int nchn, const float *buff,
                         float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                    float *out, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                       float *out, float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void sum_multi(size_t nele, int nins, const float *ins[], float *out);
    NACS_EXPORT(ds_helper) static
    void sum_multi_nt(size_t nele, int nins, const float *ins[], float *out);

    NACS_EXPORT(ds_helper) static
    void sin_single(float *out, unsigned nsteps, unsigned nrep, ChnParamFixed param);
    NACS_EXPORT(ds_helper) static
    void sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                            const ChnParamFixed *params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                              const ChnParamFixed *params, unsigned nparams);

    NACS_EXPORT(ds_helper) static
    void sin_ramp_single(float *out, unsigned nsteps, unsigned nrep,
                         ChnParamMod amp_param, ChnParamMod freq_param);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                 const ChnParamMod *amp_params,
                                 const ChnParamMod *freq_params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_multi_chnblk_loop(float *out, unsigned nsteps, unsigned nrep,
                                    const ChnParamMod *amp_params,
                                    const ChnParamMod *freq_params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                   const ChnParamMod *amp_params,
                                   const ChnParamMod *freq_params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_single_pbuf(float *out, unsigned nsteps, unsigned nrep,
                              ChnParamMod amp_param, ChnParamMod freq_param);
    NACS_EXPORT(ds_helper) static
    void sin_ramp_multi_block_loop_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                        const ChnParamMod *amp_params,
                                        const ChnParamMod *freq_params,
                                        unsigned nparams);
};

} // namespace asimd

namespace sve {
struct Kernel {
    static NACS_INLINE __attribute__((target("+sve")))
    svfloat32_t sinpif_pi(svfloat32_t d)
    {
        auto ptrue = svptrue_b32();
        auto qf = svrinta_x(ptrue, d);
        auto q = svcvt_s32_x(ptrue, qf);
        // Now `d` is the fractional part in the range `[-0.5, 0.5]`
        d = svsub_x(ptrue, d, qf);
        auto s = svmul_x(ptrue, d, d);

        // Shift the last bit of `q` to the sign bit
        // and therefore flip the sign of `d` if `q` is odd
        d = svreinterpret_f32(sveor_x(ptrue, svlsl_x(ptrue, q, 31),
                                      svreinterpret_s32(d)));

        auto u = svmad_x(ptrue, svdup_f32(-0.17818783f), s, svdup_f32(0.8098674f));
        u = svmad_x(ptrue, u, s, svdup_f32(-1.6448531f));
        return svmad_x(ptrue, svmul_x(ptrue, s, d), u, d);
    }
    static NACS_INLINE __attribute__((target("+sve")))
    svfloat64_t ramp_func(svfloat64_t x, ChnParamMod param)
    {
        auto ptrue = svptrue_b32();
        x = svmul_x(ptrue, x, svdup_f64(param.slope));
        x = svsub_x(ptrue, x, svrinta_x(ptrue, x));
        return svmad_x(ptrue, svmad_x(ptrue, x, svdup_f64(param.v2),
                                      svdup_f64(param.v1)),
                       x, svdup_f64(param.v0));
    }
    NACS_EXPORT(ds_helper) static
    void sin_range(float *out, double start, double step, unsigned nsteps);
    NACS_EXPORT(ds_helper) static
    void calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static void fill(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void fill_nt(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void copy(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void copy_nt(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void fill1(size_t nele, int *buff, int v);
    NACS_EXPORT(ds_helper) static void read1(size_t nele, const int *buff);
    NACS_EXPORT(ds_helper) static void sum(size_t nele, const float *buff1, const float *buff2);

    NACS_EXPORT(ds_helper) static
    void calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_fill_nt(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void read_calc_multi(size_t nele, int nchn, const float *buff,
                         float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                    float *out, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                       float *out, float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void sum_multi(size_t nele, int nins, const float *ins[], float *out);
    NACS_EXPORT(ds_helper) static
    void sum_multi_nt(size_t nele, int nins, const float *ins[], float *out);

    NACS_EXPORT(ds_helper) static
    void sin_single(float *out, unsigned nsteps, unsigned nrep, ChnParamFixed param);
    NACS_EXPORT(ds_helper) static
    void sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                            const ChnParamFixed *params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                              const ChnParamFixed *params, unsigned nparams);

    // I'm too lazy to implement most of these unless we actually need to use them...
    static inline
    void sin_ramp_single(float*, unsigned, unsigned, ChnParamMod, ChnParamMod)
    {}
    static inline
    void sin_ramp_multi_chn_loop(float*, unsigned, unsigned,
                                 const ChnParamMod*, const ChnParamMod*, unsigned)
    {}
    static inline
    void sin_ramp_multi_chnblk_loop(float*, unsigned, unsigned,
                                 const ChnParamMod*, const ChnParamMod*, unsigned)
    {}
    static inline
    void sin_ramp_multi_block_loop(float*, unsigned, unsigned,
                                   const ChnParamMod*, const ChnParamMod*, unsigned)
    {}
    static inline
    void sin_ramp_single_pbuf(float*, unsigned, unsigned,
                              ChnParamMod, ChnParamMod)
    {}
    static inline
    void sin_ramp_multi_block_loop_pbuf(float*, unsigned, unsigned, const ChnParamMod*,
                                        const ChnParamMod*, unsigned)
    {}
};

} // namespace sve
#endif

#if NACS_CPU_X86 || NACS_CPU_X86_64

namespace sse2 {

struct Kernel {
    static NACS_INLINE __attribute__((target("sse2")))
    __m128 sinpif_pi(__m128 d)
    {
        __m128i q = _mm_cvtps_epi32(d);
        d = d - _mm_cvtepi32_ps(q);

        __m128 s = d * d;

        // Shift the last bit of `q` to the sign bit
        // and therefore flip the sign of `d` if `q` is odd
        d = __m128(_mm_slli_epi32(q, 31) ^ __m128i(d));

        auto u = -0.17818783f * s + 0.8098674f;
        u = u * s - 1.6448531f;
        return (s * d) * u + d;
    }
    NACS_EXPORT(ds_helper) static
    void sin_range(float *out, double start, double step, unsigned nsteps);
    NACS_EXPORT(ds_helper) static
    void calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static void fill(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void fill_nt(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void copy(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void copy_nt(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void fill1(size_t nele, int *buff, int v);
    NACS_EXPORT(ds_helper) static void read1(size_t nele, const int *buff);
    NACS_EXPORT(ds_helper) static void sum(size_t nele, const float *buff1, const float *buff2);

    NACS_EXPORT(ds_helper) static
    void calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_fill_nt(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void read_calc_multi(size_t nele, int nchn, const float *buff,
                         float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                    float *out, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                       float *out, float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void sum_multi(size_t nele, int nins, const float *ins[], float *out);
    NACS_EXPORT(ds_helper) static
    void sum_multi_nt(size_t nele, int nins, const float *ins[], float *out);

    NACS_EXPORT(ds_helper) static
    void sin_single(float *out, unsigned nsteps, unsigned nrep, ChnParamFixed param);
    NACS_EXPORT(ds_helper) static
    void sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                            const ChnParamFixed *params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                              const ChnParamFixed *params, unsigned nparams);
};

} // namespace sse2

namespace avx {

struct Kernel {
    static NACS_INLINE __attribute__((target("avx")))
    __m256 sinpif_pi(__m256 d)
    {
        __m256i q = _mm256_cvtps_epi32(d);
        d = d - _mm256_cvtepi32_ps(q);

        __m256 s = d * d;

        // Shift the last bit of `q` to the sign bit
        // and therefore flip the sign of `d` if `q` is odd
        __m128i tmp[2] = {_mm256_castsi256_si128(q), _mm256_extractf128_si256(q, 1)};
        tmp[0] = _mm_slli_epi32(tmp[0], 31);
        tmp[1] = _mm_slli_epi32(tmp[1], 31);
        auto mask = _mm256_castsi128_si256(tmp[0]);
        mask = _mm256_insertf128_si256(mask, tmp[1], 1);
        d = __m256(mask ^ __m256i(d));

        auto u = -0.17818783f * s + 0.8098674f;
        u = u * s - 1.6448531f;
        return (s * d) * u + d;
    }
    NACS_EXPORT(ds_helper) static
    void sin_range(float *out, double start, double step, unsigned nsteps);
    NACS_EXPORT(ds_helper) static
    void calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static void fill(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void fill_nt(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void copy(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void copy_nt(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void fill1(size_t nele, int *buff, int v);
    NACS_EXPORT(ds_helper) static void read1(size_t nele, const int *buff);
    NACS_EXPORT(ds_helper) static void sum(size_t nele, const float *buff1, const float *buff2);

    NACS_EXPORT(ds_helper) static
    void calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_fill_nt(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void read_calc_multi(size_t nele, int nchn, const float *buff,
                         float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                    float *out, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                       float *out, float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void sum_multi(size_t nele, int nins, const float *ins[], float *out);
    NACS_EXPORT(ds_helper) static
    void sum_multi_nt(size_t nele, int nins, const float *ins[], float *out);

    NACS_EXPORT(ds_helper) static
    void sin_single(float *out, unsigned nsteps, unsigned nrep, ChnParamFixed param);
    NACS_EXPORT(ds_helper) static
    void sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                            const ChnParamFixed *params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                              const ChnParamFixed *params, unsigned nparams);
};

} // namespace avx

namespace avx2 {

struct Kernel {
    static NACS_INLINE __attribute__((target("avx2,fma")))
    __m256 sinpif_pi(__m256 d)
    {
        __m256i q = _mm256_cvtps_epi32(d);
        d = d - _mm256_cvtepi32_ps(q);

        __m256 s = d * d;

        // Shift the last bit of `q` to the sign bit
        // and therefore flip the sign of `d` if `q` is odd
        d = __m256(_mm256_slli_epi32(q, 31) ^ __m256i(d));

        auto u = -0.17818783f * s + 0.8098674f;
        u = u * s - 1.6448531f;
        return (s * d) * u + d;
    }
    NACS_EXPORT(ds_helper) static
    void sin_range(float *out, double start, double step, unsigned nsteps);
    NACS_EXPORT(ds_helper) static
    void calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static void fill(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void fill_nt(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void copy(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void copy_nt(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void fill1(size_t nele, int *buff, int v);
    NACS_EXPORT(ds_helper) static void read1(size_t nele, const int *buff);
    NACS_EXPORT(ds_helper) static void sum(size_t nele, const float *buff1, const float *buff2);

    NACS_EXPORT(ds_helper) static
    void calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_fill_nt(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void read_calc_multi(size_t nele, int nchn, const float *buff,
                         float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                    float *out, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                       float *out, float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void sum_multi(size_t nele, int nins, const float *ins[], float *out);
    NACS_EXPORT(ds_helper) static
    void sum_multi_nt(size_t nele, int nins, const float *ins[], float *out);

    NACS_EXPORT(ds_helper) static
    void sin_single(float *out, unsigned nsteps, unsigned nrep, ChnParamFixed param);
    NACS_EXPORT(ds_helper) static
    void sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                            const ChnParamFixed *params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                              const ChnParamFixed *params, unsigned nparams);
};

} // namespace avx2

namespace avx512 {

struct Kernel {
    static NACS_INLINE __attribute__((target("avx512f,avx512dq")))
    __m512 sinpif_pi(__m512 d)
    {
        __m512i q = _mm512_cvtps_epi32(d);
        d = d - _mm512_cvtepi32_ps(q);

        __m512 s = d * d;

        // Shift the last bit of `q` to the sign bit
        // and therefore flip the sign of `d` if `q` is odd
        d = __m512(_mm512_slli_epi32(q, 31) ^ __m512i(d));

        auto u = -0.17818783f * s + 0.8098674f;
        u = u * s - 1.6448531f;
        return (s * d) * u + d;
    }
    NACS_EXPORT(ds_helper) static
    void sin_range(float *out, double start, double step, unsigned nsteps);
    NACS_EXPORT(ds_helper) static
    void calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static void fill(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void fill_nt(size_t nrep, size_t ncalc, int *buff, int v);
    NACS_EXPORT(ds_helper) static void copy(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void copy_nt(size_t nele, const int *in, int *out);
    NACS_EXPORT(ds_helper) static void fill1(size_t nele, int *buff, int v);
    NACS_EXPORT(ds_helper) static void read1(size_t nele, const int *buff);
    NACS_EXPORT(ds_helper) static void sum(size_t nele, const float *buff1, const float *buff2);

    NACS_EXPORT(ds_helper) static
    void calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_fill_nt(size_t nele, int nchn, float *buff, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void read_calc_multi(size_t nele, int nchn, const float *buff,
                         float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                    float *out, float t, float freq, float amp);
    NACS_EXPORT(ds_helper) static
    void calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                       float *out, float t, float freq, float amp);

    NACS_EXPORT(ds_helper) static
    void sum_multi(size_t nele, int nins, const float *ins[], float *out);
    NACS_EXPORT(ds_helper) static
    void sum_multi_nt(size_t nele, int nins, const float *ins[], float *out);

    NACS_EXPORT(ds_helper) static
    void sin_single(float *out, unsigned nsteps, unsigned nrep, ChnParamFixed param);
    NACS_EXPORT(ds_helper) static
    void sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                            const ChnParamFixed *params, unsigned nparams);
    NACS_EXPORT(ds_helper) static
    void sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                              const ChnParamFixed *params, unsigned nparams);
};

} // namespace avx512
#endif

}

#endif
