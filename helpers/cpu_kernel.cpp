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

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_single(float *out, unsigned nsteps, unsigned nrep,
                        ChnParamFixed param)
{
    double phase = 0;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i++) {
            out[i] = float(param.amp) * sinpif_pi(float(phase));
            phase += param.freq;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                const ChnParamFixed *params, unsigned nparams)
{
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i++) {
            float v = 0;
            for (unsigned c = 0; c < nparams; c++) {
                auto phase = phases[c];
                auto param = params[c];
                v += float(param.amp) * sinpif_pi(float(phase));
                phase += param.freq;
                phases[c] = phase;
            }
            out[i] = v;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                  const ChnParamFixed *params, unsigned nparams)
{
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        auto phase = phases[0];
        auto param = params[0];
        for (unsigned i = 0; i < nsteps; i++) {
            out[i] = float(param.amp) * sinpif_pi(float(phase));
            phase += param.freq;
        }
        phases[0] = phase;

        for (unsigned c = 1; c < nparams; c++) {
            auto phase = phases[c];
            auto param = params[c];
            for (unsigned i = 0; i < nsteps; i++) {
                out[i] += float(param.amp) * sinpif_pi(float(phase));
                phase += param.freq;
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

} // namespace scalar

#if NACS_CPU_AARCH64

namespace asimd {

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

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_single(float *out, unsigned nsteps, unsigned nrep,
                        ChnParamFixed param)
{
    float32x4_t inc{0, 1, 2, 3};
    double phase = 0;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 4) {
            auto phase_v = float(phase) + inc * float(param.freq);
            vst1q_f32(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 4;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                const ChnParamFixed *params, unsigned nparams)
{
    float32x4_t inc{0, 1, 2, 3};
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 4) {
            float32x4_t v{0, 0, 0, 0};
            for (unsigned c = 0; c < nparams; c++) {
                auto phase = phases[c];
                auto param = params[c];
                auto phase_v = float(phase) + inc * float(param.freq);
                v += float(param.amp) * sinpif_pi(phase_v);
                phase += param.freq * 4;
                phases[c] = phase;
            }
            vst1q_f32(&out[i], v);
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                  const ChnParamFixed *params, unsigned nparams)
{
    float32x4_t inc{0, 1, 2, 3};
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        auto phase = phases[0];
        auto param = params[0];
        for (unsigned i = 0; i < nsteps; i += 4) {
            auto phase_v = float(phase) + inc * float(param.freq);
            vst1q_f32(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 4;
        }
        phases[0] = phase;

        for (unsigned c = 1; c < nparams; c++) {
            auto phase = phases[c];
            auto param = params[c];
            for (unsigned i = 0; i < nsteps; i += 4) {
                auto phase_v = float(phase) + inc * float(param.freq);
                vst1q_f32(&out[i], vld1q_f32(&out[i]) +
                          float(param.amp) * sinpif_pi(phase_v));
                phase += param.freq * 4;
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

} // namespace asimd

namespace sve {

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::sin_range(float *out, double start, double step, unsigned nsteps)
{
    auto svelen_2 = svcntd();
    auto svelen = svelen_2 * 2;
    auto ptrue = svptrue_b32();
    auto vstep = svdup_f64_x(ptrue, step);
    auto vstart = svdup_f64_x(ptrue, start);
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, (size_t)nsteps)));
         i += svelen) {
        auto idx1 = svindex_u64(i, 1);
        auto phase_lo64 = svmad_x(ptrue, svcvt_f64_x(ptrue, idx1), vstep, vstart);
        auto phase_lo = svcvt_f32_x(ptrue, phase_lo64);
        auto idx2 = svindex_u64(i + svelen_2, 1);
        auto phase_hi64 = svmad_x(ptrue, svcvt_f64_x(ptrue, idx2), vstep, vstart);
        auto phase_hi = svcvt_f32_x(ptrue, phase_hi64);
        svst1(pg, &out[i], sinpif_pi(svuzp1(phase_lo, phase_hi)));
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::calc_dry(size_t nrep, size_t ncalc, float t, float freq, float amp)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = svdup_f32(t);
    auto fp = svdup_f32(freq);
    auto ap = svdup_f32(amp);
    for (size_t j = 0; j < nrep; j++) {
        for (size_t i = 0; i < ncalc; i += svelen) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res = svmul_x(ptrue, ap, sinpif_pi(svmul_x(ptrue, tp, fp)));
            asm volatile ("" :: "w"(res) : "memory");
        }
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::calc_fill(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = svdup_f32(t);
    auto fp = svdup_f32(freq);
    auto ap = svdup_f32(amp);
    for (size_t j = 0; j < nrep; j++) {
        svbool_t pg;
        for (size_t i = 0;
             svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, ncalc)));
             i += svelen) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res = svmul_x(ptrue, ap, sinpif_pi(svmul_x(ptrue, tp, fp)));
            svst1(pg, &buff[i], res);
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::calc_fill_nt(size_t nrep, size_t ncalc, float *buff, float t, float freq, float amp)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto tp = svdup_f32(t);
    auto fp = svdup_f32(freq);
    auto ap = svdup_f32(amp);
    for (size_t j = 0; j < nrep; j++) {
        svbool_t pg;
        for (size_t i = 0;
             svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, ncalc)));
             i += svelen) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res = svmul_x(ptrue, ap, sinpif_pi(svmul_x(ptrue, tp, fp)));
            svstnt1(pg, &buff[i], res);
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::fill(size_t nrep, size_t ncalc, int *buff, int v)
{
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto svelen = svcntw();
    auto vp = svdup_s32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+w"(vp) :: "memory");
        svbool_t pg;
        for (size_t i = 0;
             svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, ncalc)));
             i += svelen) {
            svst1(pg, &buff[i], vp);
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::fill_nt(size_t nrep, size_t ncalc, int *buff, int v)
{
    asm volatile ("" : "+r"(nrep) :: "memory");
    asm volatile ("" : "+r"(ncalc) :: "memory");
    auto svelen = svcntw();
    auto vp = svdup_s32(v);
    for (size_t j = 0; j < nrep; j++) {
        asm volatile ("" : "+w"(vp) :: "memory");
        svbool_t pg;
        for (size_t i = 0;
             svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, ncalc)));
             i += svelen) {
            svstnt1(pg, &buff[i], vp);
        }
        asm volatile ("" :: "r"(buff) : "memory");
    }
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::copy(size_t nele, const int *in, int *out)
{
    auto svelen = svcntw();
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        svst1(pg, &out[i], svld1(pg, &in[i]));
    }
    asm volatile ("" :: "r"(out) : "memory");
}
NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::copy_nt(size_t nele, const int *in, int *out)
{
    auto svelen = svcntw();
    asm volatile ("" : "+r"(nele), "+r"(in), "+r"(out) :: "memory");
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        svstnt1(pg, &out[i], svld1(pg, &in[i]));
    }
    asm volatile ("" :: "r"(out) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::fill1(size_t nele, int *buff, int v)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    auto svelen = svcntw();
    auto vp = svdup_s32(v);
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        svst1(pg, &buff[i], vp);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::read1(size_t nele, const int *buff)
{
    asm volatile ("" : "+r"(nele) :: "memory");
    auto svelen = svcntw();
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        auto v = svld1(pg, &buff[i]);
        asm volatile ("" : "+w"(v) :: "memory");
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::sum(size_t nele, const float *buff1, const float *buff2)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" :: "r"(buff1) : "memory");
    asm volatile ("" :: "r"(buff2) : "memory");
    for (size_t i = 0; i < nele; i += svelen) {
        auto res = svadd_x(ptrue, svld1(ptrue, &buff1[i]),
                           svld1(ptrue, &buff2[i]));
        asm volatile ("" :: "w"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::calc_multi_fill(size_t nele, int nchn, float *buff, float t, float freq, float amp)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = svdup_f32(t);
    auto fp = svdup_f32(freq);
    auto ap = svdup_f32(amp);
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        auto res = svdup_f32(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res0 = svmul_x(ptrue, ap, sinpif_pi(svmul_x(ptrue, tp, fp)));
            res = svadd_x(ptrue, res, res0);
        }
        svst1(pg, &buff[i], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::calc_multi_fill_nt(size_t nele, int nchn, float *buff,
                                float t, float freq, float amp)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = svdup_f32(t);
    auto fp = svdup_f32(freq);
    auto ap = svdup_f32(amp);
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        auto res = svdup_f32(0);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res0 = svmul_x(ptrue, ap, sinpif_pi(svmul_x(ptrue, tp, fp)));
            res = svadd_x(ptrue, res, res0);
        }
        svstnt1(pg, &buff[i], res);
    }
    asm volatile ("" :: "r"(buff) : "memory");
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::read_calc_multi(size_t nele, int nchn, const float *buff,
                             float t, float freq, float amp)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    asm volatile ("" :: "r"(buff) : "memory");
    auto tp = svdup_f32(t);
    auto fp = svdup_f32(freq);
    auto ap = svdup_f32(amp);
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        auto res = svld1(pg, &buff[i]);
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res0 = svmul_x(ptrue, ap, sinpif_pi(svmul_x(ptrue, tp, fp)));
            res = svadd_x(ptrue, res, res0);
        }
        asm volatile ("" :: "w"(res) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::calc_multi(size_t nele, int nchn, int nins, const float *ins[],
                        float *out, float t, float freq, float amp)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = svdup_f32(t);
    auto fp = svdup_f32(freq);
    auto ap = svdup_f32(amp);
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        auto res = svld1(pg, &ins[0][i]);
        for (int in = 1; in < nins; in++)
            res = svadd_x(ptrue, res, svld1(pg, &ins[in][i]));
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res0 = svmul_x(ptrue, ap, sinpif_pi(svmul_x(ptrue, tp, fp)));
            res = svadd_x(ptrue, res, res0);
        }
        svst1(pg, &out[i], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::calc_multi_nt(size_t nele, int nchn, int nins, const float *ins[],
                           float *out, float t, float freq, float amp)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nele) :: "memory");
    asm volatile ("" : "+r"(nchn) :: "memory");
    auto tp = svdup_f32(t);
    auto fp = svdup_f32(freq);
    auto ap = svdup_f32(amp);
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        auto res = svld1(pg, &ins[0][i]);
        for (int in = 1; in < nins; in++)
            res = svadd_x(ptrue, res, svld1(pg, &ins[in][i]));
        for (int c = 0; c < nchn; c++) {
            asm volatile ("" : "+w"(tp), "+w"(fp), "+w"(ap) :: "memory");
            auto res0 = svmul_x(ptrue, ap, sinpif_pi(svmul_x(ptrue, tp, fp)));
            res = svadd_x(ptrue, res, res0);
        }
        svstnt1(pg, &out[i], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::sum_multi(size_t nele, int nins, const float *ins[], float *out)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nele) :: "memory");
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        auto res = svld1(pg, &ins[0][i]);
        for (int in = 1; in < nins; in++)
            res = svadd_x(ptrue, res, svld1(pg, &ins[in][i]));
        svst1(pg, &out[i], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::sum_multi_nt(size_t nele, int nins, const float *ins[], float *out)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    asm volatile ("" : "+r"(nele) :: "memory");
    svbool_t pg;
    for (size_t i = 0;
         svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, nele))); i += svelen) {
        auto res = svld1(pg, &ins[0][i]);
        for (int in = 1; in < nins; in++)
            res = svadd_x(ptrue, res, svld1(pg, &ins[in][i]));
        svstnt1(pg, &out[i], res);
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::sin_single(float *out, unsigned nsteps, unsigned nrep,
                        ChnParamFixed param)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    auto inc = svcvt_f32_x(ptrue, svindex_u32(0, 1));
    double phase = 0;
    for (unsigned j = 0; j < nrep; j++) {
        svbool_t pg;
        for (size_t i = 0;
             svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, (size_t)nsteps)));
             i += svelen) {
            auto vphase = svdup_f32_x(ptrue, float(phase));
            auto vfreq = svdup_f32_x(ptrue, float(param.freq));
            auto vamp = svdup_f32_x(ptrue, float(param.amp));
            vphase = svmad_x(ptrue, inc, vfreq, vphase);
            svst1(pg, &out[i], svmul_x(ptrue, vamp, sinpif_pi(vphase)));
            phase += param.freq * double(svelen);
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                const ChnParamFixed *params, unsigned nparams)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    auto inc = svcvt_f32_x(ptrue, svindex_u32(0, 1));
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        svbool_t pg;
        for (size_t i = 0;
             svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, (size_t)nsteps)));
             i += svelen) {
            auto v = svdup_f32_x(ptrue, 0);
            for (unsigned c = 0; c < nparams; c++) {
                auto phase = phases[c];
                auto param = params[c];
                auto vphase = svdup_f32_x(ptrue, float(phase));
                auto vfreq = svdup_f32_x(ptrue, float(param.freq));
                auto vamp = svdup_f32_x(ptrue, float(param.amp));
                vphase = svmad_x(ptrue, inc, vfreq, vphase);
                v = svmad_x(ptrue, vamp, sinpif_pi(vphase), v);
                phase += param.freq * double(svelen);
                phases[c] = phase;
            }
            svst1(pg, &out[i], v);
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("+sve"),flatten))
void Kernel::sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                  const ChnParamFixed *params, unsigned nparams)
{
    auto svelen = svcntw();
    auto ptrue = svptrue_b32();
    auto inc = svcvt_f32_x(ptrue, svindex_u32(0, 1));
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        auto phase = phases[0];
        auto param = params[0];
        auto vfreq = svdup_f32_x(ptrue, float(param.freq));
        auto vamp = svdup_f32_x(ptrue, float(param.amp));
        svbool_t pg;
        for (size_t i = 0;
             svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, (size_t)nsteps)));
             i += svelen) {
            auto vphase = svdup_f32_x(ptrue, float(phase));
            vphase = svmad_x(ptrue, inc, vfreq, vphase);
            svst1(pg, &out[i], svmul_x(ptrue, vamp, sinpif_pi(vphase)));
            phase += param.freq * double(svelen);
        }
        phases[0] = phase;

        for (unsigned c = 1; c < nparams; c++) {
            auto phase = phases[c];
            auto param = params[c];
            auto vfreq = svdup_f32_x(ptrue, float(param.freq));
            auto vamp = svdup_f32_x(ptrue, float(param.amp));
            for (size_t i = 0;
                 svptest_first(svptrue_b32(), (pg = svwhilelt_b32(i, (size_t)nsteps)));
                 i += svelen) {
                auto vphase = svdup_f32_x(ptrue, float(phase));
                vphase = svmad_x(ptrue, inc, vfreq, vphase);
                auto v = svmad_x(ptrue, vamp, sinpif_pi(vphase),
                                 svld1(pg, &out[i]));
                svst1(pg, &out[i], v);
                phase += param.freq * double(svelen);
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

} // namespace sve
#endif

#if NACS_CPU_X86 || NACS_CPU_X86_64

namespace sse2 {

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
