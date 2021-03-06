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

#include "cpu_kernel.h"

#if NACS_CPU_X86 || NACS_CPU_X86_64
#  include <immintrin.h>
#elif NACS_CPU_AARCH64
#  include <arm_neon.h>
#endif

namespace CPUKernel {

NACS_EXPORT() const NaCs::CPUInfo &host = NaCs::CPUInfo::get_host();

#if NACS_CPU_X86 || NACS_CPU_X86_64

NACS_EXPORT() bool hasavx()
{
    return host.test_feature(NaCs::X86::Feature::avx);
}

NACS_EXPORT() bool hasavx2()
{
    return host.test_feature(NaCs::X86::Feature::avx2) &&
        host.test_feature(NaCs::X86::Feature::fma);
}

NACS_EXPORT() bool hasavx512()
{
    return host.test_feature(NaCs::X86::Feature::avx512f) &&
        host.test_feature(NaCs::X86::Feature::avx512dq);
}

#endif

namespace scalar {

// Calculate `sin(pi * d) / pi`
// By computing `sin(pi * d)` instead of `sin(d)`, we don't need to do additional scaling of `d`
// before converting it to integer. Then by computing `sin(pi * d) / pi` the first order
// expansion term is kept as `d` which saves another multiplication at the end.
// We still need to generate the correct output from the correct input in the end
// so we still need the correct input/output scaling.
// However, the input needs to be scaled anyway so we can fold the input scaling in there,
// for the output, we can scale it once after all the channels are computed instead
// of doing it once per channel.
static NACS_INLINE float sinpif_pi(float d)
{
#if NACS_CPU_X86 || NACS_CPU_X86_64
    int q = _mm_cvtss_si32(_mm_set_ss(d));
#elif NACS_CPU_AARCH64
    int q = vcvtns_s32_f32(d);
#else
    int q = d < 0 ? (int)(d - 0.5f) : (int)(d + 0.5f);
#endif
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

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_range(float *out, double start, double step, unsigned nsteps)
{
    for (unsigned i = 0; i < nsteps; i++) {
        out[i] = sinpif_pi(float(start + step * i));
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp)
{
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+r"(t), "+r"(freq), "+r"(amp) :: "memory");
            auto res = amp * sinpif_pi(t * freq);
            asm volatile ("" :: "r"(res) : "memory");
        }
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+r"(t), "+r"(freq), "+r"(amp) :: "memory");
            buff[i] = amp * sinpif_pi(t * freq);
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    calc_fill(nrep, ncalc, buff, t, freq, amp);
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::fill(size_t nrep, size_t ncalc, int *buff, int v)
{
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+r"(v) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            buff[i] = v;
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::fill_nt(size_t nrep, size_t ncalc, int *buff, int v)
{
    fill(nrep, ncalc, buff, v);
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::copy(size_t nele, const int *in, int *out)
{
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        out[i] = in[i];
    asm volatile ("" :: "r"(out) : "memory");
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::copy_nt(size_t nele, const int *in, int *out)
{
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        out[i] = in[i];
    asm volatile ("" :: "r"(out) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::fill1(size_t nele, int *buff, int v)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(v) :: "memory");
    for (size_t i = 0; i < nele; i++)
        buff[i] = v;
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::read1(size_t nele, const int *buff)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto v = buff[i];
        asm volatile ("" :: "r"(v) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sum(size_t nele, const float *buff1, const float *buff2)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" :: "r"(buff1) : "memory");
    asm volatile ("" :: "r"(buff2) : "memory");
    for (size_t i = 0; i < nele; i++) {
        float res = buff1[i] + buff2[i];
        asm volatile ("" :: "r"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        float res = 0;
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+r"(t), "+r"(freq), "+r"(amp) :: "memory");
            res += amp * sinpif_pi(t * freq);
        }
        buff[i] = res;
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_multi_fill_nt(size_t nele, int nchn, float *buff,
                                float t, float freq, float amp)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        float res = 0;
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+r"(t), "+r"(freq), "+r"(amp) :: "memory");
            res += amp * sinpif_pi(t * freq);
        }
        buff[i] = res;
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::read_calc_multi(size_t nele, int nchn, const float *buff,
                             float t, float freq, float amp)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    asm volatile ("" :: "r"(buff) : "memory");
    for (size_t i = 0; i < nele; i++) {
        float res = buff[i];
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+r"(t), "+r"(freq), "+r"(amp) :: "memory");
            res += amp * sinpif_pi(t * freq);
        }
        asm volatile ("" :: "r"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                        float *out, float t, float freq, float amp)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        float res = ins[0][i];
        for (int in = 1; in < nins; in++)
            res += ins[in][i];
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+r"(t), "+r"(freq), "+r"(amp) :: "memory");
            res += amp * sinpif_pi(t * freq);
        }
        out[i] = res;
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                           float *out, float t, float freq, float amp)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        float res = ins[0][i];
        for (int in = 1; in < nins; in++)
            res += ins[in][i];
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+r"(t), "+r"(freq), "+r"(amp) :: "memory");
            res += amp * sinpif_pi(t * freq);
        }
        out[i] = res;
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sum_multi(size_t nele, int nins, const float *ins[], float *out)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        float res = ins[0][i];
        for (int in = 1; in < nins; in++)
            res += ins[in][i];
        out[i] = res;
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sum_multi_nt(size_t nele, int nins, const float *ins[], float *out)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        float res = ins[0][i];
        for (int in = 1; in < nins; in++)
            res += ins[in][i];
        out[i] = res;
    }
}

} // namespace scalar

#if NACS_CPU_AARCH64

namespace asimd {

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

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_range(float *out, double start, double step, unsigned nsteps)
{
    uint64x2_t inc{0, 1};
    for (unsigned i = 0; i < nsteps; i += 4) {
        auto phase_lo = vcvt_f32_f64(start + vcvtq_f64_u64(inc + i) * step);
        auto phase = vcvt_high_f32_f64(phase_lo, start + vcvtq_f64_u64(inc + i + 2) * step);
        vst1q_f32(&out[i], sinpif_pi(phase));
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp)
{
    ncalc = ncalc / 4;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = vdupq_n_f32(t);
    auto fp = vdupq_n_f32(freq);
    auto ap = vdupq_n_f32(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res = ap * sinpif_pi(tp * fp);
            asm volatile ("" :: "w"(res) : "memory");
        }
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 4;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = vdupq_n_f32(t);
    auto fp = vdupq_n_f32(freq);
    auto ap = vdupq_n_f32(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            vst1q_f32(&buff[i * 4], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    calc_fill(nrep, ncalc, buff, t, freq, amp);
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::fill(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 4;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = vdupq_n_s32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+w"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            vst1q_s32(&buff[i * 4], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::fill_nt(size_t nrep, size_t ncalc, int *buff, int v)
{
    // TODO I think this can use `stnp`
    fill(nrep, ncalc, buff, v);
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::copy(size_t nele, const int *in, int *out)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        vst1q_s32(&out[i * 4], vld1q_s32(&in[i * 4]));
    asm volatile ("" :: "r"(out) : "memory");
}
NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::copy_nt(size_t nele, const int *in, int *out)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        vst1q_s32(&out[i * 4], vld1q_s32(&in[i * 4]));
    asm volatile ("" :: "r"(out) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::fill1(size_t nele, int *buff, int v)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    auto vp = vdupq_n_s32(v);
    for (size_t i = 0; i < nele; i++)
        vst1q_s32(&buff[i * 4], vp);
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::read1(size_t nele, const int *buff)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto v = vld1q_s32(&buff[i * 4]);
        asm volatile ("" : "+w"(v) :: "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sum(size_t nele, const float *buff1, const float *buff2)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" :: "r"(buff1) : "memory");
    asm volatile ("" :: "r"(buff2) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = vld1q_f32(&buff1[i * 4]) + vld1q_f32(&buff2[i * 4]);
        asm volatile ("" :: "w"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = vdupq_n_f32(t);
    auto fp = vdupq_n_f32(freq);
    auto ap = vdupq_n_f32(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = vdupq_n_f32(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        vst1q_f32(&buff[i * 4], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_multi_fill_nt(size_t nele, int nchn, float *buff,
                                float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = vdupq_n_f32(t);
    auto fp = vdupq_n_f32(freq);
    auto ap = vdupq_n_f32(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = vdupq_n_f32(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        vst1q_f32(&buff[i * 4], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::read_calc_multi(size_t nele, int nchn, const float *buff,
                             float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = vdupq_n_f32(t);
    auto fp = vdupq_n_f32(freq);
    auto ap = vdupq_n_f32(amp);
    asm volatile ("" :: "r"(buff) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = vld1q_f32(&buff[i * 4]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        asm volatile ("" :: "w"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                        float *out, float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = vdupq_n_f32(t);
    auto fp = vdupq_n_f32(freq);
    auto ap = vdupq_n_f32(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = vld1q_f32(&ins[0][i * 4]);
        for (int in = 1; in < nins; in++)
            res += vld1q_f32(&ins[in][i * 4]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        vst1q_f32(&out[i * 4], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                           float *out, float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = vdupq_n_f32(t);
    auto fp = vdupq_n_f32(freq);
    auto ap = vdupq_n_f32(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = vld1q_f32(&ins[0][i * 4]);
        for (int in = 1; in < nins; in++)
            res += vld1q_f32(&ins[in][i * 4]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        vst1q_f32(&out[i * 4], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sum_multi(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = vld1q_f32(&ins[0][i * 4]);
        for (int in = 1; in < nins; in++)
            res += vld1q_f32(&ins[in][i * 4]);
        vst1q_f32(&out[i * 4], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sum_multi_nt(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = vld1q_f32(&ins[0][i * 4]);
        for (int in = 1; in < nins; in++)
            res += vld1q_f32(&ins[in][i * 4]);
        vst1q_f32(&out[i * 4], res);
    }
}

} // namespace asimd
#endif

#if NACS_CPU_X86 || NACS_CPU_X86_64

namespace sse2 {

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

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_range(float *out, double start, double step, unsigned nsteps)
{
    int32_t inc_buff[] = {0, 1, 0, 1};
    __m128i inc = _mm_loadu_si128((const __m128i*)&inc_buff);
    for (unsigned i = 0; i < nsteps; i += 4) {
        auto phase_lo = _mm_cvtpd_ps(start + _mm_cvtepi32_pd(
                                         _mm_add_epi32(inc, _mm_set1_epi32(i))) * step);
        auto phase_hi = _mm_cvtpd_ps(start + _mm_cvtepi32_pd(
                                         _mm_add_epi32(inc, _mm_set1_epi32(i + 2))) * step);
        auto phase = _mm_shuffle_ps(phase_lo, phase_hi, 0x44);
        _mm_store_ps(&out[i], sinpif_pi(phase));
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp)
{
    ncalc = ncalc / 4;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm_set1_ps(t);
    auto fp = _mm_set1_ps(freq);
    auto ap = _mm_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            auto res = ap * sinpif_pi(tp * fp);
            asm volatile ("" :: "x"(res) : "memory");
        }
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 4;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm_set1_ps(t);
    auto fp = _mm_set1_ps(freq);
    auto ap = _mm_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            _mm_store_ps(&buff[i * 4], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 4;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm_set1_ps(t);
    auto fp = _mm_set1_ps(freq);
    auto ap = _mm_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            _mm_stream_ps(&buff[i * 4], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::fill(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 4;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = _mm_set1_epi32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+x"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            _mm_store_si128((__m128i*)&buff[i * 4], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::fill_nt(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 4;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = _mm_set1_epi32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+x"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            _mm_stream_si128((__m128i*)&buff[i * 4], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::copy(size_t nele, const int *in, int *out)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        _mm_store_si128((__m128i*)&out[i * 4], _mm_load_si128((const __m128i*)&in[i * 4]));
    asm volatile ("" :: "r"(out) : "memory");
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::copy_nt(size_t nele, const int *in, int *out)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        _mm_stream_si128((__m128i*)&out[i * 4], _mm_load_si128((const __m128i*)&in[i * 4]));
    asm volatile ("" :: "r"(out) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::fill1(size_t nele, int *buff, int v)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    auto vp = _mm_set1_epi32(v);
    for (size_t i = 0; i < nele; i++)
        _mm_store_si128((__m128i*)&buff[i * 4], vp);
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::read1(size_t nele, const int *buff)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto v = _mm_load_si128((const __m128i*)&buff[i * 4]);
        asm volatile ("" : "+x"(v) :: "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sum(size_t nele, const float *buff1, const float *buff2)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" :: "r"(buff1) : "memory");
    asm volatile ("" :: "r"(buff2) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm_load_ps(&buff1[i * 4]) + _mm_load_ps(&buff2[i * 4]);
        asm volatile ("" :: "x"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm_set1_ps(t);
    auto fp = _mm_set1_ps(freq);
    auto ap = _mm_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm_set1_ps(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm_store_ps(&buff[i * 4], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::calc_multi_fill_nt(size_t nele, int nchn, float *buff,
                                float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm_set1_ps(t);
    auto fp = _mm_set1_ps(freq);
    auto ap = _mm_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm_set1_ps(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm_stream_ps(&buff[i * 4], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::read_calc_multi(size_t nele, int nchn, const float *buff,
                             float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm_set1_ps(t);
    auto fp = _mm_set1_ps(freq);
    auto ap = _mm_set1_ps(amp);
    asm volatile ("" :: "r"(buff) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm_load_ps(&buff[i * 4]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        asm volatile ("" :: "x"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                        float *out, float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm_set1_ps(t);
    auto fp = _mm_set1_ps(freq);
    auto ap = _mm_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm_load_ps(&ins[0][i * 4]);
        for (int in = 1; in < nins; in++)
            res += _mm_load_ps(&ins[in][i * 4]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm_store_ps(&out[i * 4], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                           float *out, float t, float freq, float amp)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm_set1_ps(t);
    auto fp = _mm_set1_ps(freq);
    auto ap = _mm_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm_load_ps(&ins[0][i * 4]);
        for (int in = 1; in < nins; in++)
            res += _mm_load_ps(&ins[in][i * 4]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm_stream_ps(&out[i * 4], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sum_multi(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm_load_ps(&ins[0][i * 4]);
        for (int in = 1; in < nins; in++)
            res += _mm_load_ps(&ins[in][i * 4]);
        _mm_store_ps(&out[i * 4], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sum_multi_nt(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 4;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm_load_ps(&ins[0][i * 4]);
        for (int in = 1; in < nins; in++)
            res += _mm_load_ps(&ins[in][i * 4]);
        _mm_stream_ps(&out[i * 4], res);
    }
}

} // namespace sse2

namespace avx {

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

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_range(float *out, double start, double step, unsigned nsteps)
{
    int32_t inc_buff[] = {0, 1, 2, 3};
    __m128i inc = _mm_loadu_si128((const __m128i*)&inc_buff);
    for (unsigned i = 0; i < nsteps; i += 8) {
        auto phase_lo = _mm256_cvtpd_ps(
            start + _mm256_cvtepi32_pd(_mm_add_epi32(inc, _mm_set1_epi32(i))) * step);
        auto phase_hi = _mm256_cvtpd_ps(
            start + _mm256_cvtepi32_pd(_mm_add_epi32(inc, _mm_set1_epi32(i + 4))) * step);
        auto phase = _mm256_insertf128_ps(_mm256_castps128_ps256(phase_lo), phase_hi, 1);
        _mm256_store_ps(&out[i], sinpif_pi(phase));
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            auto res = ap * sinpif_pi(tp * fp);
            asm volatile ("" :: "x"(res) : "memory");
        }
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            _mm256_store_ps(&buff[i * 8], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            _mm256_stream_ps(&buff[i * 8], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::fill(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = _mm256_set1_epi32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+x"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            _mm256_store_si256((__m256i*)&buff[i * 8], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::fill_nt(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = _mm256_set1_epi32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+x"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            _mm256_stream_si256((__m256i*)&buff[i * 8], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::copy(size_t nele, const int *in, int *out)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        _mm256_store_si256((__m256i*)&out[i * 8],
                           _mm256_load_si256((const __m256i*)&in[i * 8]));
    asm volatile ("" :: "r"(out) : "memory");
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::copy_nt(size_t nele, const int *in, int *out)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        _mm256_stream_si256((__m256i*)&out[i * 8],
                            _mm256_load_si256((const __m256i*)&in[i * 8]));
    asm volatile ("" :: "r"(out) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::fill1(size_t nele, int *buff, int v)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    auto vp = _mm256_set1_epi32(v);
    for (size_t i = 0; i < nele; i++)
        _mm256_store_si256((__m256i*)&buff[i * 8], vp);
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::read1(size_t nele, const int *buff)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto v = _mm256_load_si256((const __m256i*)&buff[i * 8]);
        asm volatile ("" : "+x"(v) :: "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sum(size_t nele, const float *buff1, const float *buff2)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" :: "r"(buff1) : "memory");
    asm volatile ("" :: "r"(buff2) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&buff1[i * 8]) + _mm256_load_ps(&buff2[i * 8]);
        asm volatile ("" :: "x"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_set1_ps(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm256_store_ps(&buff[i * 8], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::calc_multi_fill_nt(size_t nele, int nchn, float *buff,
                                float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_set1_ps(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm256_stream_ps(&buff[i * 8], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::read_calc_multi(size_t nele, int nchn,
                             const float *buff, float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    asm volatile ("" :: "r"(buff) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&buff[i * 8]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        asm volatile ("" :: "x"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                        float *out, float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&ins[0][i * 8]);
        for (int in = 1; in < nins; in++)
            res += _mm256_load_ps(&ins[in][i * 8]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm256_store_ps(&out[i * 8], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                           float *out, float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&ins[0][i * 8]);
        for (int in = 1; in < nins; in++)
            res += _mm256_load_ps(&ins[in][i * 8]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm256_stream_ps(&out[i * 8], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sum_multi(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&ins[0][i * 8]);
        for (int in = 1; in < nins; in++)
            res += _mm256_load_ps(&ins[in][i * 8]);
        _mm256_store_ps(&out[i * 8], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sum_multi_nt(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&ins[0][i * 8]);
        for (int in = 1; in < nins; in++)
            res += _mm256_load_ps(&ins[in][i * 8]);
        _mm256_stream_ps(&out[i * 8], res);
    }
}

} // namespace avx

namespace avx2 {

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

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_range(float *out, double start, double step, unsigned nsteps)
{
    int32_t inc_buff[] = {0, 1, 2, 3};
    __m128i inc = _mm_loadu_si128((const __m128i*)&inc_buff);
    for (unsigned i = 0; i < nsteps; i += 8) {
        auto phase_lo = _mm256_cvtpd_ps(
            start + _mm256_cvtepi32_pd(_mm_add_epi32(inc, _mm_set1_epi32(i))) * step);
        auto phase_hi = _mm256_cvtpd_ps(
            start + _mm256_cvtepi32_pd(_mm_add_epi32(inc, _mm_set1_epi32(i + 4))) * step);
        auto phase = _mm256_insertf128_ps(_mm256_castps128_ps256(phase_lo), phase_hi, 1);
        _mm256_store_ps(&out[i], sinpif_pi(phase));
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            auto res = ap * sinpif_pi(tp * fp);
            asm volatile ("" :: "x"(res) : "memory");
        }
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            _mm256_store_ps(&buff[i * 8], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            _mm256_stream_ps(&buff[i * 8], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::fill(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = _mm256_set1_epi32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+x"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            _mm256_store_si256((__m256i*)&buff[i * 8], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::fill_nt(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 8;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = _mm256_set1_epi32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+x"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            _mm256_stream_si256((__m256i*)&buff[i * 8], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::copy(size_t nele, const int *in, int *out)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        _mm256_store_si256((__m256i*)&out[i * 8],
                           _mm256_load_si256((const __m256i*)&in[i * 8]));
    asm volatile ("" :: "r"(out) : "memory");
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::copy_nt(size_t nele, const int *in, int *out)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        _mm256_stream_si256((__m256i*)&out[i * 8],
                            _mm256_load_si256((const __m256i*)&in[i * 8]));
    asm volatile ("" :: "r"(out) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::fill1(size_t nele, int *buff, int v)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    auto vp = _mm256_set1_epi32(v);
    for (size_t i = 0; i < nele; i++)
        _mm256_store_si256((__m256i*)&buff[i * 8], vp);
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::read1(size_t nele, const int *buff)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto v = _mm256_load_si256((const __m256i*)&buff[i * 8]);
        asm volatile ("" : "+x"(v) :: "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sum(size_t nele, const float *buff1, const float *buff2)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" :: "r"(buff1) : "memory");
    asm volatile ("" :: "r"(buff2) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&buff1[i * 8]) + _mm256_load_ps(&buff2[i * 8]);
        asm volatile ("" :: "x"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_set1_ps(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm256_store_ps(&buff[i * 8], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::calc_multi_fill_nt(size_t nele, int nchn, float *buff,
                                float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_set1_ps(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm256_stream_ps(&buff[i * 8], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::read_calc_multi(size_t nele, int nchn, const float *buff,
                             float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    asm volatile ("" :: "r"(buff) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&buff[i * 8]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        asm volatile ("" :: "x"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                        float *out, float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&ins[0][i * 8]);
        for (int in = 1; in < nins; in++)
            res += _mm256_load_ps(&ins[in][i * 8]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm256_store_ps(&out[i * 8], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                           float *out, float t, float freq, float amp)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm256_set1_ps(t);
    auto fp = _mm256_set1_ps(freq);
    auto ap = _mm256_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&ins[0][i * 8]);
        for (int in = 1; in < nins; in++)
            res += _mm256_load_ps(&ins[in][i * 8]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm256_stream_ps(&out[i * 8], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sum_multi(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&ins[0][i * 8]);
        for (int in = 1; in < nins; in++)
            res += _mm256_load_ps(&ins[in][i * 8]);
        _mm256_store_ps(&out[i * 8], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sum_multi_nt(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 8;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm256_load_ps(&ins[0][i * 8]);
        for (int in = 1; in < nins; in++)
            res += _mm256_load_ps(&ins[in][i * 8]);
        _mm256_stream_ps(&out[i * 8], res);
    }
}

} // namespace avx2

namespace avx512 {

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

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_range(float *out, double start, double step, unsigned nsteps)
{
    int32_t inc_buff[] = {0, 1, 2, 3, 4, 5, 6, 7};
    __m256i inc = _mm256_loadu_si256((const __m256i*)&inc_buff);
    for (unsigned i = 0; i < nsteps; i += 16) {
        auto phase_lo = _mm512_cvtpd_ps(
            start + _mm512_cvtepi32_pd(_mm256_add_epi32(inc, _mm256_set1_epi32(i))) * step);
        auto phase_hi = _mm512_cvtpd_ps(
            start + _mm512_cvtepi32_pd(_mm256_add_epi32(inc, _mm256_set1_epi32(i + 8))) * step);
        auto phase = _mm512_insertf32x8(_mm512_castps256_ps512(phase_lo), phase_hi, 1);
        _mm512_store_ps(&out[i], sinpif_pi(phase));
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp)
{
    ncalc = ncalc / 16;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm512_set1_ps(t);
    auto fp = _mm512_set1_ps(freq);
    auto ap = _mm512_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            auto res = ap * sinpif_pi(tp * fp);
            asm volatile ("" :: "x"(res) : "memory");
        }
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 16;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm512_set1_ps(t);
    auto fp = _mm512_set1_ps(freq);
    auto ap = _mm512_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            _mm512_store_ps(&buff[i * 16], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    ncalc = ncalc / 16;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = _mm512_set1_ps(t);
    auto fp = _mm512_set1_ps(freq);
    auto ap = _mm512_set1_ps(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            _mm512_stream_ps(&buff[i * 16], ap * sinpif_pi(tp * fp));
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::fill(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 16;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = _mm512_set1_epi32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+x"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            _mm512_store_si512((__m512i*)&buff[i * 16], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::fill_nt(size_t nrep, size_t ncalc, int *buff, int v)
{
    ncalc = ncalc / 16;
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto vp = _mm512_set1_epi32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+x"(vp) :: "memory");
        for (size_t i = 0; i < ncalc; i++)
            _mm512_stream_si512((__m512i*)&buff[i * 16], vp);
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::copy(size_t nele, const int *in, int *out)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        _mm512_store_si512((__m512i*)&out[i * 16],
                           _mm512_load_si512((const __m512i*)&in[i * 16]));
    asm volatile ("" :: "r"(out) : "memory");
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::copy_nt(size_t nele, const int *in, int *out)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    for (size_t i = 0; i < nele; i++)
        _mm512_stream_si512((__m512i*)&out[i * 16],
                           _mm512_load_si512((const __m512i*)&in[i * 16]));
    asm volatile ("" :: "r"(out) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::fill1(size_t nele, int *buff, int v)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    auto vp = _mm512_set1_epi32(v);
    for (size_t i = 0; i < nele; i++)
        _mm512_store_si512((__m512i*)&buff[i * 16], vp);
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::read1(size_t nele, const int *buff)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto v = _mm512_load_si512((const __m512i*)&buff[i * 16]);
        asm volatile ("" : "+x"(v) :: "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sum(size_t nele, const float *buff1, const float *buff2)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" :: "r"(buff1) : "memory");
    asm volatile ("" :: "r"(buff2) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm512_load_ps(&buff1[i * 16]) + _mm512_load_ps(&buff2[i * 16]);
        asm volatile ("" :: "x"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm512_set1_ps(t);
    auto fp = _mm512_set1_ps(freq);
    auto ap = _mm512_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm512_set1_ps(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm512_store_ps(&buff[i * 16], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::calc_multi_fill_nt(size_t nele, int nchn, float *buff,
                                float t, float freq, float amp)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm512_set1_ps(t);
    auto fp = _mm512_set1_ps(freq);
    auto ap = _mm512_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm512_set1_ps(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm512_stream_ps(&buff[i * 16], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::read_calc_multi(size_t nele, int nchn, const float *buff,
                             float t, float freq, float amp)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm512_set1_ps(t);
    auto fp = _mm512_set1_ps(freq);
    auto ap = _mm512_set1_ps(amp);
    asm volatile ("" :: "r"(buff) : "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm512_load_ps(&buff[i * 16]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        asm volatile ("" :: "x"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                        float *out, float t, float freq, float amp)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm512_set1_ps(t);
    auto fp = _mm512_set1_ps(freq);
    auto ap = _mm512_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm512_load_ps(&ins[0][i * 16]);
        for (int in = 1; in < nins; in++)
            res += _mm512_load_ps(&ins[in][i * 16]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm512_store_ps(&out[i * 16], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                           float *out, float t, float freq, float amp)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = _mm512_set1_ps(t);
    auto fp = _mm512_set1_ps(freq);
    auto ap = _mm512_set1_ps(amp);
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm512_load_ps(&ins[0][i * 16]);
        for (int in = 1; in < nins; in++)
            res += _mm512_load_ps(&ins[in][i * 16]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+x"(tp), "+x"(fp), "+x"(ap) :: "memory");
            res += ap * sinpif_pi(tp * fp);
        }
        _mm512_stream_ps(&out[i * 16], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sum_multi(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm512_load_ps(&ins[0][i * 16]);
        for (int in = 1; in < nins; in++)
            res += _mm512_load_ps(&ins[in][i * 16]);
        _mm512_store_ps(&out[i * 16], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sum_multi_nt(size_t nele, int nins, const float *ins[], float *out)
{
    nele = nele / 16;
    asm volatile ("" : "+r"(nele) :: "memory");
    for (size_t i = 0; i < nele; i++) {
        auto res = _mm512_load_ps(&ins[0][i * 16]);
        for (int in = 1; in < nins; in++)
            res += _mm512_load_ps(&ins[in][i * 16]);
        _mm512_stream_ps(&out[i * 16], res);
    }
}

} // namespace avx512
#endif

}
