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

#if NACS_CPU_AARCH64

NACS_EXPORT() bool hassve()
{
    return host.test_feature(NaCs::AArch64::Feature::sve);
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
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
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
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
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
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        phases[0] = phase;

        for (unsigned c = 1; c < nparams; c++) {
            auto phase = phases[c];
            auto param = params[c];
            for (unsigned i = 0; i < nsteps; i++) {
                out[i] += float(param.amp) * sinpif_pi(float(phase));
                phase += param.freq;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

static constexpr float new_amp_coeff[] = {
    0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375,
    0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375};
static constexpr float old_amp_coeff[] = {
    1.0, 0.9375, 0.875, 0.8125, 0.75, 0.6875, 0.625, 0.5625,
    0.5, 0.4375, 0.375, 0.3125, 0.25, 0.1875, 0.125, 0.0625};
static constexpr float new_freq_coeff[] = {
    0.0, 0.03125, 0.125, 0.28125, 0.5, 0.78125, 1.125, 1.53125,
    2.0, 2.53125, 3.125, 3.78125, 4.5, 5.28125, 6.125, 7.03125};
static constexpr float old_freq_coeff[] = {
    0.0, 0.96875, 1.875, 2.71875, 3.5, 4.21875, 4.875, 5.46875,
    6.0, 6.46875, 6.875, 7.21875, 7.5, 7.71875, 7.875, 7.96875};

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_single(float *out, unsigned nsteps, unsigned nrep,
                             ChnParamMod amp_param, ChnParamMod freq_param)
{
    double phase = 0;
    double amp = ramp_func(0, amp_param);
    double freq = ramp_func(0, freq_param);
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 16) {
            double new_amp = ramp_func(double(t), amp_param);
            double new_freq = ramp_func(double(t), freq_param);
            t += 16;
            for (unsigned si = 0; si < 16; si++) {
                float a = float(new_amp) * new_amp_coeff[si] +
                    float(amp) * old_amp_coeff[si];
                float p = float(phase) +
                    float(new_freq) * new_freq_coeff[si] +
                    float(freq) * old_freq_coeff[si];
                out[i + si] = a * sinpif_pi(p);
            }
            phase = phase + (freq + new_freq) * 8;
            freq = new_freq;
            amp = new_amp;
            if (phase > 32) {
                phase -= 64;
            }
            else if (phase < -32) {
                phase += 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                     const ChnParamMod *amp_params,
                                     const ChnParamMod *freq_params,
                                     unsigned nparams)
{
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(0, amp_params[c]);
        freqs[c] = ramp_func(0, freq_params[c]);
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 16) {
            double new_freqs[nparams];
            double new_amps[nparams];
            for (unsigned c = 0; c < nparams; c++) {
                new_amps[c] = ramp_func(double(t), amp_params[c]);
                new_freqs[c] = ramp_func(double(t), freq_params[c]);
            }
            t += 16;
            for (unsigned si = 0; si < 16; si++) {
                float v = 0;
                for (unsigned c = 0; c < nparams; c++) {
                    float a = float(new_amps[c]) * new_amp_coeff[si] +
                        float(amps[c]) * old_amp_coeff[si];
                    float p = float(phases[c]) +
                        float(new_freqs[c]) * new_freq_coeff[si] +
                        float(freqs[c]) * old_freq_coeff[si];
                    v += a * sinpif_pi(p);
                }
                out[i + si] = v;
            }
            for (unsigned c = 0; c < nparams; c++) {
                auto phase = phases[c] + (freqs[c] + new_freqs[c]) * 8;
                amps[c] = new_amps[c];
                freqs[c] = new_freqs[c];
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
                phases[c] = phase;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_chnblk_loop(float *out, unsigned nsteps, unsigned nrep,
                                        const ChnParamMod *amp_params,
                                        const ChnParamMod *freq_params,
                                        unsigned nparams)
{
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(0, amp_params[c]);
        freqs[c] = ramp_func(0, freq_params[c]);
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto chn_func = [&] (unsigned c) {
                auto new_amp = ramp_func(double(t), amp_params[c]);
                auto new_freq = ramp_func(double(t), freq_params[c]);
                auto amp = amps[c];
                auto freq = freqs[c];
                auto phase = phases[c];
                for (unsigned si = 0; si < 16; si++) {
                    float a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    float p = float(phase) + float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    if (c == 0) {
                        out[i + si] = a * sinpif_pi(p);
                    }
                    else {
                        out[i + si] += a * sinpif_pi(p);
                    }
                }
                phase += (freq + new_freq) * 8;
                amps[c] = new_amp;
                freqs[c] = new_freq;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
                phases[c] = phase;
            };

            // Encourage the compiler to specialize for c=0
            chn_func(0);
            for (unsigned c = 1; c < nparams; c++)
                chn_func(c);
            t += 16;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                       const ChnParamMod *amp_params,
                                       const ChnParamMod *freq_params,
                                       unsigned nparams)
{
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(0, amp_params[c]);
        freqs[c] = ramp_func(0, freq_params[c]);
    }
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto new_amp = ramp_func(double(t), amp_param);
                auto new_freq = ramp_func(double(t), freq_param);
                t += 16;
                for (unsigned si = 0; si < 16; si++) {
                    float a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    float p = float(phase) + float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    if (c == 0) {
                        out[i + si] = a * sinpif_pi(p);
                    }
                    else {
                        out[i + si] += a * sinpif_pi(p);
                    }
                }
                phase += (freq + new_freq) * 8;
                amp = new_amp;
                freq = new_freq;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_single_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                  ChnParamMod amp_param, ChnParamMod freq_param)
{
    double phase = 0;
    double amp = ramp_func(0, amp_param);
    double freq = ramp_func(0, freq_param);
    uint64_t t = 16;
    double ampbuf[nsteps / 16];
    double freqbuf[nsteps / 16];
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps / 16; i++) {
            ampbuf[i] = ramp_func(double(t), amp_param);
            freqbuf[i] = ramp_func(double(t), freq_param);
            t += 16;
        }
        for (unsigned i = 0; i < nsteps; i += 16) {
            double new_amp = ampbuf[i / 16];
            double new_freq = freqbuf[i / 16];
            for (unsigned si = 0; si < 16; si++) {
                float a = float(new_amp) * new_amp_coeff[si] +
                    float(amp) * old_amp_coeff[si];
                float p = float(phase) +
                    float(new_freq) * new_freq_coeff[si] +
                    float(freq) * old_freq_coeff[si];
                out[i + si] = a * sinpif_pi(p);
            }
            phase = phase + (freq + new_freq) * 8;
            freq = new_freq;
            amp = new_amp;
            if (phase > 32) {
                phase -= 64;
            }
            else if (phase < -32) {
                phase += 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                            const ChnParamMod *amp_params,
                                            const ChnParamMod *freq_params,
                                            unsigned nparams)
{
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(0, amp_params[c]);
        freqs[c] = ramp_func(0, freq_params[c]);
    }
    double ampbuf[nsteps / 16];
    double freqbuf[nsteps / 16];
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps / 16; i++) {
                ampbuf[i] = ramp_func(double(t), amp_param);
                freqbuf[i] = ramp_func(double(t), freq_param);
                t += 16;
            }
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto new_amp = ampbuf[i / 16];
                auto new_freq = freqbuf[i / 16];
                for (unsigned si = 0; si < 16; si++) {
                    float a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    float p = float(phase) + float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    if (c == 0) {
                        out[i + si] = a * sinpif_pi(p);
                    }
                    else {
                        out[i + si] += a * sinpif_pi(p);
                    }
                }
                phase += (freq + new_freq) * 8;
                amp = new_amp;
                freq = new_freq;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_single_pbuf_full(float *out, unsigned nsteps, unsigned nrep,
                                       ChnParamMod amp_param, ChnParamMod freq_param)
{
    double phase = 0;
    double amp = ramp_func(0, amp_param);
    double freq = ramp_func(0, freq_param);
    uint64_t t = 16;
    float ampbuf[nsteps];
    float phasebuf[nsteps];
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto new_amp = ramp_func(double(t), amp_param);
            auto new_freq = ramp_func(double(t), freq_param);
            t += 16;
            for (unsigned si = 0; si < 16; si++) {
                ampbuf[i + si] = float(new_amp) * new_amp_coeff[si] +
                    float(amp) * old_amp_coeff[si];
                phasebuf[i + si] = float(phase) +
                    float(new_freq) * new_freq_coeff[si] +
                    float(freq) * old_freq_coeff[si];
            }
            phase = phase + (freq + new_freq) * 8;
            freq = new_freq;
            amp = new_amp;
            if (phase > 32) {
                phase -= 64;
            }
            else if (phase < -32) {
                phase += 64;
            }
        }
        for (unsigned i = 0; i < nsteps; i++) {
            out[i] = ampbuf[i] * sinpif_pi(phasebuf[i]);
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf_full(
    float *out, unsigned nsteps, unsigned nrep,
    const ChnParamMod *amp_params, const ChnParamMod *freq_params, unsigned nparams)
{
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(0, amp_params[c]);
        freqs[c] = ramp_func(0, freq_params[c]);
    }
    float ampbuf[nsteps];
    float phasebuf[nsteps];
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto new_amp = ramp_func(double(t), amp_param);
                auto new_freq = ramp_func(double(t), freq_param);
                t += 16;
                for (unsigned si = 0; si < 16; si++) {
                    ampbuf[i + si] = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    phasebuf[i + si] = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                }
                phase = phase + (freq + new_freq) * 8;
                freq = new_freq;
                amp = new_amp;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
            for (unsigned i = 0; i < nsteps; i++) {
                if (c == 0) {
                    out[i] = ampbuf[i] * sinpif_pi(phasebuf[i]);
                }
                else {
                    out[i] += ampbuf[i] * sinpif_pi(phasebuf[i]);
                }
            }
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
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
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
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
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
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
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
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
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

static constexpr float32x4_t new_amp_coeff[] = {
    {0.0, 0.0625, 0.125, 0.1875}, {0.25, 0.3125, 0.375, 0.4375},
    {0.5, 0.5625, 0.625, 0.6875}, {0.75, 0.8125, 0.875, 0.9375}};
static constexpr float32x4_t old_amp_coeff[] = {
    {1.0, 0.9375, 0.875, 0.8125}, {0.75, 0.6875, 0.625, 0.5625},
    {0.5, 0.4375, 0.375, 0.3125}, {0.25, 0.1875, 0.125, 0.0625}};
static constexpr float32x4_t new_freq_coeff[] = {
    {0.0, 0.03125, 0.125, 0.28125}, {0.5, 0.78125, 1.125, 1.53125},
    {2.0, 2.53125, 3.125, 3.78125}, {4.5, 5.28125, 6.125, 7.03125}};
static constexpr float32x4_t old_freq_coeff[] = {
    {0.0, 0.96875, 1.875, 2.71875}, {3.5, 4.21875, 4.875, 5.46875},
    {6.0, 6.46875, 6.875, 7.21875}, {7.5, 7.71875, 7.875, 7.96875}};

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_single(float *out, unsigned nsteps, unsigned nrep,
                             ChnParamMod amp_param, ChnParamMod freq_param)
{
    float64x2_t pinc{0, 16};

    double phase = 0;
    auto amp = ramp_func(float64x2_t{0, 0}, amp_param)[0];
    auto freq = ramp_func(float64x2_t{0, 0}, freq_param)[0];
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 32) {
            auto new_amps = ramp_func(double(t) + pinc, amp_param);
            auto new_freqs = ramp_func(double(t) + pinc, freq_param);
            t += 32;
            for (unsigned pi = 0; pi < 2; pi++) {
                auto new_amp = new_amps[pi];
                auto new_freq = new_freqs[pi];
                for (unsigned si = 0; si < 4; si++) {
                    auto a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    vst1q_f32(&out[i + pi * 16 + si * 4], a * sinpif_pi(p));
                }
                phase = phase + (freq + new_freq) * 8;
                freq = new_freq;
                amp = new_amp;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                     const ChnParamMod *amp_params,
                                     const ChnParamMod *freq_params,
                                     unsigned nparams)
{
    float64x2_t pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(float64x2_t{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(float64x2_t{0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 32) {
            float64x2_t new_freqs[nparams];
            float64x2_t new_amps[nparams];
            for (unsigned c = 0; c < nparams; c++) {
                new_amps[c] = ramp_func(double(t) + pinc, amp_params[c]);
                new_freqs[c] = ramp_func(double(t) + pinc, freq_params[c]);
            }
            t += 32;
            for (unsigned pi = 0; pi < 2; pi++) {
                for (unsigned si = 0; si < 4; si++) {
                    float32x4_t v = {0, 0, 0, 0};
                    for (unsigned c = 0; c < nparams; c++) {
                        auto a = float(new_amps[c][pi]) * new_amp_coeff[si] +
                            float(amps[c]) * old_amp_coeff[si];
                        auto p = float(phases[c]) +
                            float(new_freqs[c][pi]) * new_freq_coeff[si] +
                            float(freqs[c]) * old_freq_coeff[si];
                        v += a * sinpif_pi(p);
                    }
                    vst1q_f32(&out[i + pi * 16 + si * 4], v);
                }
                for (unsigned c = 0; c < nparams; c++) {
                    auto phase = phases[c] + (freqs[c] + new_freqs[c][pi]) * 8;
                    amps[c] = new_amps[c][pi];
                    freqs[c] = new_freqs[c][pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    phases[c] = phase;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_chnblk_loop(float *out, unsigned nsteps, unsigned nrep,
                                        const ChnParamMod *amp_params,
                                        const ChnParamMod *freq_params,
                                        unsigned nparams)
{
    float64x2_t pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(float64x2_t{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(float64x2_t{0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 32) {
            auto chn_func = [&] (unsigned c) {
                auto new_amp = ramp_func(double(t) + pinc, amp_params[c]);
                auto new_freq = ramp_func(double(t) + pinc, freq_params[c]);
                auto amp = amps[c];
                auto freq = freqs[c];
                auto phase = phases[c];
                for (unsigned pi = 0; pi < 2; pi++) {
                    for (unsigned si = 0; si < 4; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        auto v = a * sinpif_pi(p);
                        if (c != 0)
                            v += vld1q_f32(&out[i + pi * 16 + si * 4]);
                        vst1q_f32(&out[i + pi * 16 + si * 4], v);
                    }
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
                amps[c] = amp;
                freqs[c] = freq;
                phases[c] = phase;
            };

            // Encourage the compiler to specialize for c=0
            chn_func(0);
            for (unsigned c = 1; c < nparams; c++)
                chn_func(c);
            t += 32;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                       const ChnParamMod *amp_params,
                                       const ChnParamMod *freq_params,
                                       unsigned nparams)
{
    float64x2_t pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(float64x2_t{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(float64x2_t{0, 0}, freq_params[c])[0];
    }
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 32) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 32;
                for (unsigned pi = 0; pi < 2; pi++) {
                    for (unsigned si = 0; si < 4; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        auto v = a * sinpif_pi(p);
                        if (c != 0)
                            v += vld1q_f32(&out[i + pi * 16 + si * 4]);
                        vst1q_f32(&out[i + pi * 16 + si * 4], v);
                    }
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_single_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                  ChnParamMod amp_param, ChnParamMod freq_param)
{
    float64x2_t pinc{0, 16};

    double phase = 0;
    auto amp = ramp_func(float64x2_t{0, 0}, amp_param)[0];
    auto freq = ramp_func(float64x2_t{0, 0}, freq_param)[0];
    uint64_t t = 16;
    double ampbuf[nsteps / 16] __attribute__((aligned(16)));
    double freqbuf[nsteps / 16] __attribute__((aligned(16)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps / 32; i++) {
            vst1q_f64(&ampbuf[i * 2], ramp_func(double(t) + pinc, amp_param));
            vst1q_f64(&freqbuf[i * 2], ramp_func(double(t) + pinc, freq_param));
            t += 32;
        }
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto new_amp = ampbuf[i / 16];
            auto new_freq = freqbuf[i / 16];
            for (unsigned si = 0; si < 4; si++) {
                auto a = float(new_amp) * new_amp_coeff[si] +
                    float(amp) * old_amp_coeff[si];
                auto p = float(phase) +
                    float(new_freq) * new_freq_coeff[si] +
                    float(freq) * old_freq_coeff[si];
                vst1q_f32(&out[i + si * 4], a * sinpif_pi(p));
            }
            phase = phase + (freq + new_freq) * 8;
            freq = new_freq;
            amp = new_amp;
            if (phase > 32) {
                phase -= 64;
            }
            else if (phase < -32) {
                phase += 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                            const ChnParamMod *amp_params,
                                            const ChnParamMod *freq_params,
                                            unsigned nparams)
{
    float64x2_t pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(float64x2_t{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(float64x2_t{0, 0}, freq_params[c])[0];
    }
    double ampbuf[nsteps / 16] __attribute__((aligned(16)));
    double freqbuf[nsteps / 16] __attribute__((aligned(16)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps / 32; i++) {
                vst1q_f64(&ampbuf[i * 2], ramp_func(double(t) + pinc, amp_param));
                vst1q_f64(&freqbuf[i * 2], ramp_func(double(t) + pinc, freq_param));
                t += 32;
            }
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto new_amp = ampbuf[i / 16];
                auto new_freq = freqbuf[i / 16];
                for (unsigned si = 0; si < 4; si++) {
                    auto a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    auto v = a * sinpif_pi(p);
                    if (c != 0)
                        v += vld1q_f32(&out[i + si * 4]);
                    vst1q_f32(&out[i + si * 4], v);
                }
                phase += (freq + new_freq) * 8;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
                amp = new_amp;
                freq = new_freq;
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_single_pbuf_full(float *out, unsigned nsteps, unsigned nrep,
                                       ChnParamMod amp_param, ChnParamMod freq_param)
{
    float64x2_t pinc{0, 16};

    double phase = 0;
    auto amp = ramp_func(float64x2_t{0, 0}, amp_param)[0];
    auto freq = ramp_func(float64x2_t{0, 0}, freq_param)[0];
    uint64_t t = 16;
    float ampbuf[nsteps] __attribute__((aligned(16)));
    float phasebuf[nsteps] __attribute__((aligned(16)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 32) {
            auto new_amp = ramp_func(double(t) + pinc, amp_param);
            auto new_freq = ramp_func(double(t) + pinc, freq_param);
            t += 32;
            for (unsigned pi = 0; pi < 2; pi++) {
                for (unsigned si = 0; si < 4; si++) {
                    auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq[pi]) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    vst1q_f32(&ampbuf[i + pi * 16 + si * 4], a);
                    vst1q_f32(&phasebuf[i + pi * 16 + si * 4], p);
                }
                phase = phase + (freq + new_freq[pi]) * 8;
                freq = new_freq[pi];
                amp = new_amp[pi];
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        for (unsigned i = 0; i < nsteps; i += 4) {
            vst1q_f32(&out[i], vld1q_f32(&ampbuf[i]) *
                      sinpif_pi(vld1q_f32(&phasebuf[i])));
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf_full(
    float *out, unsigned nsteps, unsigned nrep,
    const ChnParamMod *amp_params, const ChnParamMod *freq_params, unsigned nparams)
{
    float64x2_t pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(float64x2_t{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(float64x2_t{0, 0}, freq_params[c])[0];
    }
    float ampbuf[nsteps] __attribute__((aligned(16)));
    float phasebuf[nsteps] __attribute__((aligned(16)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 32) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 32;
                for (unsigned pi = 0; pi < 2; pi++) {
                    for (unsigned si = 0; si < 4; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        vst1q_f32(&ampbuf[i + pi * 16 + si * 4], a);
                        vst1q_f32(&phasebuf[i + pi * 16 + si * 4], p);
                    }
                    phase = phase + (freq + new_freq[pi]) * 8;
                    freq = new_freq[pi];
                    amp = new_amp[pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
            for (unsigned i = 0; i < nsteps; i += 4) {
                auto v = vld1q_f32(&ampbuf[i]) * sinpif_pi(vld1q_f32(&phasebuf[i]));
                if (c != 0)
                    v += vld1q_f32(&out[i]);
                vst1q_f32(&out[i], v);
            }
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
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
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
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
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
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
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
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
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
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

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_single(float *out, unsigned nsteps, unsigned nrep,
                        ChnParamFixed param)
{
    float inc_buf[] = {0, 1, 2, 3};
    auto inc = _mm_load_ps(inc_buf);
    double phase = 0;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 4) {
            auto phase_v = float(phase) + inc * float(param.freq);
            _mm_store_ps(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 4;
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                const ChnParamFixed *params, unsigned nparams)
{
    float inc_buf[] = {0, 1, 2, 3};
    auto inc = _mm_load_ps(inc_buf);
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 4) {
            auto v = _mm_set1_ps(0);
            for (unsigned c = 0; c < nparams; c++) {
                auto phase = phases[c];
                auto param = params[c];
                auto phase_v = float(phase) + inc * float(param.freq);
                v += float(param.amp) * sinpif_pi(phase_v);
                phase += param.freq * 4;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
                phases[c] = phase;
            }
            _mm_store_ps(&out[i], v);
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                  const ChnParamFixed *params, unsigned nparams)
{
    float inc_buf[] = {0, 1, 2, 3};
    auto inc = _mm_load_ps(inc_buf);
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        auto phase = phases[0];
        auto param = params[0];
        for (unsigned i = 0; i < nsteps; i += 4) {
            auto phase_v = float(phase) + inc * float(param.freq);
            _mm_store_ps(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 4;
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        phases[0] = phase;

        for (unsigned c = 1; c < nparams; c++) {
            auto phase = phases[c];
            auto param = params[c];
            for (unsigned i = 0; i < nsteps; i += 4) {
                auto phase_v = float(phase) + inc * float(param.freq);
                _mm_store_ps(&out[i], _mm_load_ps(&out[i]) +
                             float(param.amp) * sinpif_pi(phase_v));
                phase += param.freq * 4;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

static constexpr __m128 new_amp_coeff[] = {
    {0.0, 0.0625, 0.125, 0.1875}, {0.25, 0.3125, 0.375, 0.4375},
    {0.5, 0.5625, 0.625, 0.6875}, {0.75, 0.8125, 0.875, 0.9375}};
static constexpr __m128 old_amp_coeff[] = {
    {1.0, 0.9375, 0.875, 0.8125}, {0.75, 0.6875, 0.625, 0.5625},
    {0.5, 0.4375, 0.375, 0.3125}, {0.25, 0.1875, 0.125, 0.0625}};
static constexpr __m128 new_freq_coeff[] = {
    {0.0, 0.03125, 0.125, 0.28125}, {0.5, 0.78125, 1.125, 1.53125},
    {2.0, 2.53125, 3.125, 3.78125}, {4.5, 5.28125, 6.125, 7.03125}};
static constexpr __m128 old_freq_coeff[] = {
    {0.0, 0.96875, 1.875, 2.71875}, {3.5, 4.21875, 4.875, 5.46875},
    {6.0, 6.46875, 6.875, 7.21875}, {7.5, 7.71875, 7.875, 7.96875}};

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_ramp_single(float *out, unsigned nsteps, unsigned nrep,
                             ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m128d pinc{0, 16};

    double phase = 0;
    auto amp = ramp_func(__m128d{0, 0}, amp_param)[0];
    auto freq = ramp_func(__m128d{0, 0}, freq_param)[0];
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 32) {
            auto new_amps = ramp_func(double(t) + pinc, amp_param);
            auto new_freqs = ramp_func(double(t) + pinc, freq_param);
            t += 32;
            for (unsigned pi = 0; pi < 2; pi++) {
                auto new_amp = new_amps[pi];
                auto new_freq = new_freqs[pi];
                for (unsigned si = 0; si < 4; si++) {
                    auto a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    _mm_store_ps(&out[i + pi * 16 + si * 4], a * sinpif_pi(p));
                }
                phase = phase + (freq + new_freq) * 8;
                freq = new_freq;
                amp = new_amp;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_ramp_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                     const ChnParamMod *amp_params,
                                     const ChnParamMod *freq_params,
                                     unsigned nparams)
{
    __m128d pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m128d{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m128d{0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 32) {
            __m128d new_freqs[nparams];
            __m128d new_amps[nparams];
            for (unsigned c = 0; c < nparams; c++) {
                new_amps[c] = ramp_func(double(t) + pinc, amp_params[c]);
                new_freqs[c] = ramp_func(double(t) + pinc, freq_params[c]);
            }
            t += 32;
            for (unsigned pi = 0; pi < 2; pi++) {
                for (unsigned si = 0; si < 4; si++) {
                    __m128 v = {0, 0, 0, 0};
                    for (unsigned c = 0; c < nparams; c++) {
                        auto a = float(new_amps[c][pi]) * new_amp_coeff[si] +
                            float(amps[c]) * old_amp_coeff[si];
                        auto p = float(phases[c]) +
                            float(new_freqs[c][pi]) * new_freq_coeff[si] +
                            float(freqs[c]) * old_freq_coeff[si];
                        v += a * sinpif_pi(p);
                    }
                    _mm_store_ps(&out[i + pi * 16 + si * 4], v);
                }
                for (unsigned c = 0; c < nparams; c++) {
                    auto phase = phases[c] + (freqs[c] + new_freqs[c][pi]) * 8;
                    amps[c] = new_amps[c][pi];
                    freqs[c] = new_freqs[c][pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    phases[c] = phase;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_ramp_multi_chnblk_loop(float *out, unsigned nsteps, unsigned nrep,
                                        const ChnParamMod *amp_params,
                                        const ChnParamMod *freq_params,
                                        unsigned nparams)
{
    __m128d pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m128d{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m128d{0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 32) {
            auto chn_func = [&] (unsigned c) {
                auto new_amp = ramp_func(double(t) + pinc, amp_params[c]);
                auto new_freq = ramp_func(double(t) + pinc, freq_params[c]);
                auto amp = amps[c];
                auto freq = freqs[c];
                auto phase = phases[c];
                for (unsigned pi = 0; pi < 2; pi++) {
                    for (unsigned si = 0; si < 4; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        auto v = a * sinpif_pi(p);
                        if (c != 0)
                            v += _mm_load_ps(&out[i + pi * 16 + si * 4]);
                        _mm_store_ps(&out[i + pi * 16 + si * 4], v);
                    }
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
                amps[c] = amp;
                freqs[c] = freq;
                phases[c] = phase;
            };

            // Encourage the compiler to specialize for c=0
            chn_func(0);
            for (unsigned c = 1; c < nparams; c++)
                chn_func(c);
            t += 32;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_ramp_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                       const ChnParamMod *amp_params,
                                       const ChnParamMod *freq_params,
                                       unsigned nparams)
{
    __m128d pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m128d{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m128d{0, 0}, freq_params[c])[0];
    }
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 32) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 32;
                for (unsigned pi = 0; pi < 2; pi++) {
                    for (unsigned si = 0; si < 4; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        auto v = a * sinpif_pi(p);
                        if (c != 0)
                            v += _mm_load_ps(&out[i + pi * 16 + si * 4]);
                        _mm_store_ps(&out[i + pi * 16 + si * 4], v);
                    }
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_ramp_single_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                  ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m128d pinc{0, 16};

    double phase = 0;
    auto amp = ramp_func(__m128d{0, 0}, amp_param)[0];
    auto freq = ramp_func(__m128d{0, 0}, freq_param)[0];
    uint64_t t = 16;
    double ampbuf[nsteps / 16] __attribute__((aligned(16)));
    double freqbuf[nsteps / 16] __attribute__((aligned(16)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps / 32; i++) {
            _mm_store_pd(&ampbuf[i * 2], ramp_func(double(t) + pinc, amp_param));
            _mm_store_pd(&freqbuf[i * 2], ramp_func(double(t) + pinc, freq_param));
            t += 32;
        }
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto new_amp = ampbuf[i / 16];
            auto new_freq = freqbuf[i / 16];
            for (unsigned si = 0; si < 4; si++) {
                auto a = float(new_amp) * new_amp_coeff[si] +
                    float(amp) * old_amp_coeff[si];
                auto p = float(phase) +
                    float(new_freq) * new_freq_coeff[si] +
                    float(freq) * old_freq_coeff[si];
                _mm_store_ps(&out[i + si * 4], a * sinpif_pi(p));
            }
            phase = phase + (freq + new_freq) * 8;
            freq = new_freq;
            amp = new_amp;
            if (phase > 32) {
                phase -= 64;
            }
            else if (phase < -32) {
                phase += 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                            const ChnParamMod *amp_params,
                                            const ChnParamMod *freq_params,
                                            unsigned nparams)
{
    __m128d pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m128d{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m128d{0, 0}, freq_params[c])[0];
    }
    double ampbuf[nsteps / 16] __attribute__((aligned(16)));
    double freqbuf[nsteps / 16] __attribute__((aligned(16)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps / 32; i++) {
                _mm_store_pd(&ampbuf[i * 2], ramp_func(double(t) + pinc, amp_param));
                _mm_store_pd(&freqbuf[i * 2], ramp_func(double(t) + pinc, freq_param));
                t += 32;
            }
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto new_amp = ampbuf[i / 16];
                auto new_freq = freqbuf[i / 16];
                for (unsigned si = 0; si < 4; si++) {
                    auto a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    auto v = a * sinpif_pi(p);
                    if (c != 0)
                        v += _mm_load_ps(&out[i + si * 4]);
                    _mm_store_ps(&out[i + si * 4], v);
                }
                phase += (freq + new_freq) * 8;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
                amp = new_amp;
                freq = new_freq;
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_ramp_single_pbuf_full(float *out, unsigned nsteps, unsigned nrep,
                                       ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m128d pinc{0, 16};

    double phase = 0;
    auto amp = ramp_func(__m128d{0, 0}, amp_param)[0];
    auto freq = ramp_func(__m128d{0, 0}, freq_param)[0];
    uint64_t t = 16;
    float ampbuf[nsteps] __attribute__((aligned(16)));
    float phasebuf[nsteps] __attribute__((aligned(16)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 32) {
            auto new_amp = ramp_func(double(t) + pinc, amp_param);
            auto new_freq = ramp_func(double(t) + pinc, freq_param);
            t += 32;
            for (unsigned pi = 0; pi < 2; pi++) {
                for (unsigned si = 0; si < 4; si++) {
                    auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq[pi]) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    _mm_store_ps(&ampbuf[i + pi * 16 + si * 4], a);
                    _mm_store_ps(&phasebuf[i + pi * 16 + si * 4], p);
                }
                phase = phase + (freq + new_freq[pi]) * 8;
                freq = new_freq[pi];
                amp = new_amp[pi];
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        for (unsigned i = 0; i < nsteps; i += 4) {
            _mm_store_ps(&out[i], _mm_load_ps(&ampbuf[i]) *
                         sinpif_pi(_mm_load_ps(&phasebuf[i])));
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("sse2"),flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf_full(
    float *out, unsigned nsteps, unsigned nrep,
    const ChnParamMod *amp_params, const ChnParamMod *freq_params, unsigned nparams)
{
    __m128d pinc{0, 16};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m128d{0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m128d{0, 0}, freq_params[c])[0];
    }
    float ampbuf[nsteps] __attribute__((aligned(16)));
    float phasebuf[nsteps] __attribute__((aligned(16)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 32) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 32;
                for (unsigned pi = 0; pi < 2; pi++) {
                    for (unsigned si = 0; si < 4; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        _mm_store_ps(&ampbuf[i + pi * 16 + si * 4], a);
                        _mm_store_ps(&phasebuf[i + pi * 16 + si * 4], p);
                    }
                    phase = phase + (freq + new_freq[pi]) * 8;
                    freq = new_freq[pi];
                    amp = new_amp[pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
            for (unsigned i = 0; i < nsteps; i += 4) {
                auto v = _mm_load_ps(&ampbuf[i]) * sinpif_pi(_mm_load_ps(&phasebuf[i]));
                if (c != 0)
                    v += _mm_load_ps(&out[i]);
                _mm_store_ps(&out[i], v);
            }
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
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

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_single(float *out, unsigned nsteps, unsigned nrep,
                        ChnParamFixed param)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7};
    auto inc = _mm256_load_ps(inc_buf);
    double phase = 0;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 8) {
            auto phase_v = float(phase) + inc * float(param.freq);
            _mm256_store_ps(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 8;
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                const ChnParamFixed *params, unsigned nparams)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7};
    auto inc = _mm256_load_ps(inc_buf);
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 8) {
            auto v = _mm256_set1_ps(0);
            for (unsigned c = 0; c < nparams; c++) {
                auto phase = phases[c];
                auto param = params[c];
                auto phase_v = float(phase) + inc * float(param.freq);
                v += float(param.amp) * sinpif_pi(phase_v);
                phase += param.freq * 8;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
                phases[c] = phase;
            }
            _mm256_store_ps(&out[i], v);
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                  const ChnParamFixed *params, unsigned nparams)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7};
    auto inc = _mm256_load_ps(inc_buf);
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        auto phase = phases[0];
        auto param = params[0];
        for (unsigned i = 0; i < nsteps; i += 8) {
            auto phase_v = float(phase) + inc * float(param.freq);
            _mm256_store_ps(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 8;
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        phases[0] = phase;

        for (unsigned c = 1; c < nparams; c++) {
            auto phase = phases[c];
            auto param = params[c];
            for (unsigned i = 0; i < nsteps; i += 8) {
                auto phase_v = float(phase) + inc * float(param.freq);
                _mm256_store_ps(&out[i], _mm256_load_ps(&out[i]) +
                                float(param.amp) * sinpif_pi(phase_v));
                phase += param.freq * 8;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

static constexpr __m256 new_amp_coeff[] = {
    {0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375},
    {0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375}};
static constexpr __m256 old_amp_coeff[] = {
    {1.0, 0.9375, 0.875, 0.8125, 0.75, 0.6875, 0.625, 0.5625},
    {0.5, 0.4375, 0.375, 0.3125, 0.25, 0.1875, 0.125, 0.0625}};
static constexpr __m256 new_freq_coeff[] = {
    {0.0, 0.03125, 0.125, 0.28125, 0.5, 0.78125, 1.125, 1.53125},
    {2.0, 2.53125, 3.125, 3.78125, 4.5, 5.28125, 6.125, 7.03125}};
static constexpr __m256 old_freq_coeff[] = {
    {0.0, 0.96875, 1.875, 2.71875, 3.5, 4.21875, 4.875, 5.46875},
    {6.0, 6.46875, 6.875, 7.21875, 7.5, 7.71875, 7.875, 7.96875}};

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_ramp_single(float *out, unsigned nsteps, unsigned nrep,
                             ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m256d pinc{0, 16, 32, 48};

    double phase = 0;
    auto amp = ramp_func(__m256d{0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m256d{0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 64) {
            auto new_amps = ramp_func(double(t) + pinc, amp_param);
            auto new_freqs = ramp_func(double(t) + pinc, freq_param);
            t += 64;
            for (unsigned pi = 0; pi < 4; pi++) {
                auto new_amp = new_amps[pi];
                auto new_freq = new_freqs[pi];
                for (unsigned si = 0; si < 2; si++) {
                    auto a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    _mm256_store_ps(&out[i + pi * 16 + si * 8], a * sinpif_pi(p));
                }
                phase = phase + (freq + new_freq) * 8;
                freq = new_freq;
                amp = new_amp;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_ramp_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                     const ChnParamMod *amp_params,
                                     const ChnParamMod *freq_params,
                                     unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 64) {
            __m256d new_freqs[nparams];
            __m256d new_amps[nparams];
            for (unsigned c = 0; c < nparams; c++) {
                new_amps[c] = ramp_func(double(t) + pinc, amp_params[c]);
                new_freqs[c] = ramp_func(double(t) + pinc, freq_params[c]);
            }
            t += 64;
            for (unsigned pi = 0; pi < 4; pi++) {
                for (unsigned si = 0; si < 2; si++) {
                    __m256 v = {0, 0, 0, 0, 0, 0, 0, 0};
                    for (unsigned c = 0; c < nparams; c++) {
                        auto a = float(new_amps[c][pi]) * new_amp_coeff[si] +
                            float(amps[c]) * old_amp_coeff[si];
                        auto p = float(phases[c]) +
                            float(new_freqs[c][pi]) * new_freq_coeff[si] +
                            float(freqs[c]) * old_freq_coeff[si];
                        v += a * sinpif_pi(p);
                    }
                    _mm256_store_ps(&out[i + pi * 16 + si * 8], v);
                }
                for (unsigned c = 0; c < nparams; c++) {
                    auto phase = phases[c] + (freqs[c] + new_freqs[c][pi]) * 8;
                    amps[c] = new_amps[c][pi];
                    freqs[c] = new_freqs[c][pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    phases[c] = phase;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_ramp_multi_chnblk_loop(float *out, unsigned nsteps, unsigned nrep,
                                        const ChnParamMod *amp_params,
                                        const ChnParamMod *freq_params,
                                        unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 64) {
            auto chn_func = [&] (unsigned c) __attribute__((target("avx"))) {
                auto new_amp = ramp_func(double(t) + pinc, amp_params[c]);
                auto new_freq = ramp_func(double(t) + pinc, freq_params[c]);
                auto amp = amps[c];
                auto freq = freqs[c];
                auto phase = phases[c];
                for (unsigned pi = 0; pi < 4; pi++) {
                    for (unsigned si = 0; si < 2; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        auto v = a * sinpif_pi(p);
                        if (c != 0)
                            v += _mm256_load_ps(&out[i + pi * 16 + si * 8]);
                        _mm256_store_ps(&out[i + pi * 16 + si * 8], v);
                    }
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
                amps[c] = amp;
                freqs[c] = freq;
                phases[c] = phase;
            };

            // Encourage the compiler to specialize for c=0
            chn_func(0);
            for (unsigned c = 1; c < nparams; c++)
                chn_func(c);
            t += 64;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_ramp_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                       const ChnParamMod *amp_params,
                                       const ChnParamMod *freq_params,
                                       unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 64) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 64;
                for (unsigned pi = 0; pi < 4; pi++) {
                    for (unsigned si = 0; si < 2; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        auto v = a * sinpif_pi(p);
                        if (c != 0)
                            v += _mm256_load_ps(&out[i + pi * 16 + si * 8]);
                        _mm256_store_ps(&out[i + pi * 16 + si * 8], v);
                    }
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_ramp_single_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                  ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m256d pinc{0, 16, 32, 48};

    double phase = 0;
    auto amp = ramp_func(__m256d{0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m256d{0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    double ampbuf[nsteps / 16] __attribute__((aligned(32)));
    double freqbuf[nsteps / 16] __attribute__((aligned(32)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps / 64; i++) {
            _mm256_store_pd(&ampbuf[i * 4], ramp_func(double(t) + pinc, amp_param));
            _mm256_store_pd(&freqbuf[i * 4], ramp_func(double(t) + pinc, freq_param));
            t += 64;
        }
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto new_amp = ampbuf[i / 16];
            auto new_freq = freqbuf[i / 16];
            for (unsigned si = 0; si < 2; si++) {
                auto a = float(new_amp) * new_amp_coeff[si] +
                    float(amp) * old_amp_coeff[si];
                auto p = float(phase) +
                    float(new_freq) * new_freq_coeff[si] +
                    float(freq) * old_freq_coeff[si];
                _mm256_store_ps(&out[i + si * 8], a * sinpif_pi(p));
            }
            phase = phase + (freq + new_freq) * 8;
            freq = new_freq;
            amp = new_amp;
            if (phase > 32) {
                phase -= 64;
            }
            else if (phase < -32) {
                phase += 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                            const ChnParamMod *amp_params,
                                            const ChnParamMod *freq_params,
                                            unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    double ampbuf[nsteps / 16] __attribute__((aligned(32)));
    double freqbuf[nsteps / 16] __attribute__((aligned(32)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps / 64; i++) {
                _mm256_store_pd(&ampbuf[i * 4],
                                ramp_func(double(t) + pinc, amp_param));
                _mm256_store_pd(&freqbuf[i * 4],
                                ramp_func(double(t) + pinc, freq_param));
                t += 64;
            }
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto new_amp = ampbuf[i / 16];
                auto new_freq = freqbuf[i / 16];
                for (unsigned si = 0; si < 2; si++) {
                    auto a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    auto v = a * sinpif_pi(p);
                    if (c != 0)
                        v += _mm256_load_ps(&out[i + si * 8]);
                    _mm256_store_ps(&out[i + si * 8], v);
                }
                phase += (freq + new_freq) * 8;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
                amp = new_amp;
                freq = new_freq;
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_ramp_single_pbuf_full(float *out, unsigned nsteps, unsigned nrep,
                                       ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m256d pinc{0, 16, 32, 48};

    double phase = 0;
    auto amp = ramp_func(__m256d{0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m256d{0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    float ampbuf[nsteps] __attribute__((aligned(32)));
    float phasebuf[nsteps] __attribute__((aligned(32)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 64) {
            auto new_amp = ramp_func(double(t) + pinc, amp_param);
            auto new_freq = ramp_func(double(t) + pinc, freq_param);
            t += 64;
            for (unsigned pi = 0; pi < 4; pi++) {
                for (unsigned si = 0; si < 2; si++) {
                    auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq[pi]) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    _mm256_store_ps(&ampbuf[i + pi * 16 + si * 8], a);
                    _mm256_store_ps(&phasebuf[i + pi * 16 + si * 8], p);
                }
                phase = phase + (freq + new_freq[pi]) * 8;
                freq = new_freq[pi];
                amp = new_amp[pi];
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        for (unsigned i = 0; i < nsteps; i += 8) {
            _mm256_store_ps(&out[i], _mm256_load_ps(&ampbuf[i]) *
                            sinpif_pi(_mm256_load_ps(&phasebuf[i])));
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx"),flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf_full(
    float *out, unsigned nsteps, unsigned nrep,
    const ChnParamMod *amp_params, const ChnParamMod *freq_params, unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    float ampbuf[nsteps] __attribute__((aligned(32)));
    float phasebuf[nsteps] __attribute__((aligned(32)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 64) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 64;
                for (unsigned pi = 0; pi < 4; pi++) {
                    for (unsigned si = 0; si < 2; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        _mm256_store_ps(&ampbuf[i + pi * 16 + si * 8], a);
                        _mm256_store_ps(&phasebuf[i + pi * 16 + si * 8], p);
                    }
                    phase = phase + (freq + new_freq[pi]) * 8;
                    freq = new_freq[pi];
                    amp = new_amp[pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
            for (unsigned i = 0; i < nsteps; i += 8) {
                auto v = _mm256_load_ps(&ampbuf[i]) *
                    sinpif_pi(_mm256_load_ps(&phasebuf[i]));
                if (c != 0)
                    v += _mm256_load_ps(&out[i]);
                _mm256_store_ps(&out[i], v);
            }
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
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

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_single(float *out, unsigned nsteps, unsigned nrep,
                        ChnParamFixed param)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7};
    auto inc = _mm256_load_ps(inc_buf);
    double phase = 0;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 8) {
            auto phase_v = float(phase) + inc * float(param.freq);
            _mm256_store_ps(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 8;
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                const ChnParamFixed *params, unsigned nparams)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7};
    auto inc = _mm256_load_ps(inc_buf);
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 8) {
            auto v = _mm256_set1_ps(0);
            for (unsigned c = 0; c < nparams; c++) {
                auto phase = phases[c];
                auto param = params[c];
                auto phase_v = float(phase) + inc * float(param.freq);
                v += float(param.amp) * sinpif_pi(phase_v);
                phase += param.freq * 8;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
                phases[c] = phase;
            }
            _mm256_store_ps(&out[i], v);
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                  const ChnParamFixed *params, unsigned nparams)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7};
    auto inc = _mm256_load_ps(inc_buf);
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        auto phase = phases[0];
        auto param = params[0];
        for (unsigned i = 0; i < nsteps; i += 8) {
            auto phase_v = float(phase) + inc * float(param.freq);
            _mm256_store_ps(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 8;
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        phases[0] = phase;

        for (unsigned c = 1; c < nparams; c++) {
            auto phase = phases[c];
            auto param = params[c];
            for (unsigned i = 0; i < nsteps; i += 8) {
                auto phase_v = float(phase) + inc * float(param.freq);
                _mm256_store_ps(&out[i], _mm256_load_ps(&out[i]) +
                                float(param.amp) * sinpif_pi(phase_v));
                phase += param.freq * 8;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

static constexpr __m256 new_amp_coeff[] = {
    {0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375},
    {0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375}};
static constexpr __m256 old_amp_coeff[] = {
    {1.0, 0.9375, 0.875, 0.8125, 0.75, 0.6875, 0.625, 0.5625},
    {0.5, 0.4375, 0.375, 0.3125, 0.25, 0.1875, 0.125, 0.0625}};
static constexpr __m256 new_freq_coeff[] = {
    {0.0, 0.03125, 0.125, 0.28125, 0.5, 0.78125, 1.125, 1.53125},
    {2.0, 2.53125, 3.125, 3.78125, 4.5, 5.28125, 6.125, 7.03125}};
static constexpr __m256 old_freq_coeff[] = {
    {0.0, 0.96875, 1.875, 2.71875, 3.5, 4.21875, 4.875, 5.46875},
    {6.0, 6.46875, 6.875, 7.21875, 7.5, 7.71875, 7.875, 7.96875}};

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_ramp_single(float *out, unsigned nsteps, unsigned nrep,
                             ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m256d pinc{0, 16, 32, 48};

    double phase = 0;
    auto amp = ramp_func(__m256d{0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m256d{0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 64) {
            auto new_amps = ramp_func(double(t) + pinc, amp_param);
            auto new_freqs = ramp_func(double(t) + pinc, freq_param);
            t += 64;
            for (unsigned pi = 0; pi < 4; pi++) {
                auto new_amp = new_amps[pi];
                auto new_freq = new_freqs[pi];
                for (unsigned si = 0; si < 2; si++) {
                    auto a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    _mm256_store_ps(&out[i + pi * 16 + si * 8], a * sinpif_pi(p));
                }
                phase = phase + (freq + new_freq) * 8;
                freq = new_freq;
                amp = new_amp;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_ramp_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                     const ChnParamMod *amp_params,
                                     const ChnParamMod *freq_params,
                                     unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 64) {
            __m256d new_freqs[nparams];
            __m256d new_amps[nparams];
            for (unsigned c = 0; c < nparams; c++) {
                new_amps[c] = ramp_func(double(t) + pinc, amp_params[c]);
                new_freqs[c] = ramp_func(double(t) + pinc, freq_params[c]);
            }
            t += 64;
            for (unsigned pi = 0; pi < 4; pi++) {
                for (unsigned si = 0; si < 2; si++) {
                    __m256 v = {0, 0, 0, 0, 0, 0, 0, 0};
                    for (unsigned c = 0; c < nparams; c++) {
                        auto a = float(new_amps[c][pi]) * new_amp_coeff[si] +
                            float(amps[c]) * old_amp_coeff[si];
                        auto p = float(phases[c]) +
                            float(new_freqs[c][pi]) * new_freq_coeff[si] +
                            float(freqs[c]) * old_freq_coeff[si];
                        v += a * sinpif_pi(p);
                    }
                    _mm256_store_ps(&out[i + pi * 16 + si * 8], v);
                }
                for (unsigned c = 0; c < nparams; c++) {
                    auto phase = phases[c] + (freqs[c] + new_freqs[c][pi]) * 8;
                    amps[c] = new_amps[c][pi];
                    freqs[c] = new_freqs[c][pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    phases[c] = phase;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_ramp_multi_chnblk_loop(float *out, unsigned nsteps, unsigned nrep,
                                        const ChnParamMod *amp_params,
                                        const ChnParamMod *freq_params,
                                        unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 64) {
            auto chn_func = [&] (unsigned c) __attribute__((target("avx2,fma"))) {
                auto new_amp = ramp_func(double(t) + pinc, amp_params[c]);
                auto new_freq = ramp_func(double(t) + pinc, freq_params[c]);
                auto amp = amps[c];
                auto freq = freqs[c];
                auto phase = phases[c];
                for (unsigned pi = 0; pi < 4; pi++) {
                    for (unsigned si = 0; si < 2; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        auto v = a * sinpif_pi(p);
                        if (c != 0)
                            v += _mm256_load_ps(&out[i + pi * 16 + si * 8]);
                        _mm256_store_ps(&out[i + pi * 16 + si * 8], v);
                    }
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
                amps[c] = amp;
                freqs[c] = freq;
                phases[c] = phase;
            };

            // Encourage the compiler to specialize for c=0
            chn_func(0);
            for (unsigned c = 1; c < nparams; c++)
                chn_func(c);
            t += 64;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_ramp_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                       const ChnParamMod *amp_params,
                                       const ChnParamMod *freq_params,
                                       unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx2,fma"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 64) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 64;
                for (unsigned pi = 0; pi < 4; pi++) {
                    for (unsigned si = 0; si < 2; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        auto v = a * sinpif_pi(p);
                        if (c != 0)
                            v += _mm256_load_ps(&out[i + pi * 16 + si * 8]);
                        _mm256_store_ps(&out[i + pi * 16 + si * 8], v);
                    }
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_ramp_single_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                  ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m256d pinc{0, 16, 32, 48};

    double phase = 0;
    auto amp = ramp_func(__m256d{0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m256d{0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    double ampbuf[nsteps / 16] __attribute__((aligned(32)));
    double freqbuf[nsteps / 16] __attribute__((aligned(32)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps / 64; i++) {
            _mm256_store_pd(&ampbuf[i * 4], ramp_func(double(t) + pinc, amp_param));
            _mm256_store_pd(&freqbuf[i * 4], ramp_func(double(t) + pinc, freq_param));
            t += 64;
        }
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto new_amp = ampbuf[i / 16];
            auto new_freq = freqbuf[i / 16];
            for (unsigned si = 0; si < 2; si++) {
                auto a = float(new_amp) * new_amp_coeff[si] +
                    float(amp) * old_amp_coeff[si];
                auto p = float(phase) +
                    float(new_freq) * new_freq_coeff[si] +
                    float(freq) * old_freq_coeff[si];
                _mm256_store_ps(&out[i + si * 8], a * sinpif_pi(p));
            }
            phase = phase + (freq + new_freq) * 8;
            freq = new_freq;
            amp = new_amp;
            if (phase > 32) {
                phase -= 64;
            }
            else if (phase < -32) {
                phase += 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                            const ChnParamMod *amp_params,
                                            const ChnParamMod *freq_params,
                                            unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    double ampbuf[nsteps / 16] __attribute__((aligned(32)));
    double freqbuf[nsteps / 16] __attribute__((aligned(32)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx2,fma"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps / 64; i++) {
                _mm256_store_pd(&ampbuf[i * 4],
                                ramp_func(double(t) + pinc, amp_param));
                _mm256_store_pd(&freqbuf[i * 4],
                                ramp_func(double(t) + pinc, freq_param));
                t += 64;
            }
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto new_amp = ampbuf[i / 16];
                auto new_freq = freqbuf[i / 16];
                for (unsigned si = 0; si < 2; si++) {
                    auto a = float(new_amp) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    auto v = a * sinpif_pi(p);
                    if (c != 0)
                        v += _mm256_load_ps(&out[i + si * 8]);
                    _mm256_store_ps(&out[i + si * 8], v);
                }
                phase += (freq + new_freq) * 8;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
                amp = new_amp;
                freq = new_freq;
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_ramp_single_pbuf_full(float *out, unsigned nsteps, unsigned nrep,
                                       ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m256d pinc{0, 16, 32, 48};

    double phase = 0;
    auto amp = ramp_func(__m256d{0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m256d{0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    float ampbuf[nsteps] __attribute__((aligned(32)));
    float phasebuf[nsteps] __attribute__((aligned(32)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 64) {
            auto new_amp = ramp_func(double(t) + pinc, amp_param);
            auto new_freq = ramp_func(double(t) + pinc, freq_param);
            t += 64;
            for (unsigned pi = 0; pi < 4; pi++) {
                for (unsigned si = 0; si < 2; si++) {
                    auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                        float(amp) * old_amp_coeff[si];
                    auto p = float(phase) +
                        float(new_freq[pi]) * new_freq_coeff[si] +
                        float(freq) * old_freq_coeff[si];
                    _mm256_store_ps(&ampbuf[i + pi * 16 + si * 8], a);
                    _mm256_store_ps(&phasebuf[i + pi * 16 + si * 8], p);
                }
                phase = phase + (freq + new_freq[pi]) * 8;
                freq = new_freq[pi];
                amp = new_amp[pi];
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        for (unsigned i = 0; i < nsteps; i += 8) {
            _mm256_store_ps(&out[i], _mm256_load_ps(&ampbuf[i]) *
                            sinpif_pi(_mm256_load_ps(&phasebuf[i])));
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx2,fma"),flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf_full(
    float *out, unsigned nsteps, unsigned nrep,
    const ChnParamMod *amp_params, const ChnParamMod *freq_params, unsigned nparams)
{
    __m256d pinc{0, 16, 32, 48};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m256d{0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m256d{0, 0, 0, 0}, freq_params[c])[0];
    }
    float ampbuf[nsteps] __attribute__((aligned(32)));
    float phasebuf[nsteps] __attribute__((aligned(32)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx2,fma"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 64) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 64;
                for (unsigned pi = 0; pi < 4; pi++) {
                    for (unsigned si = 0; si < 2; si++) {
                        auto a = float(new_amp[pi]) * new_amp_coeff[si] +
                            float(amp) * old_amp_coeff[si];
                        auto p = float(phase) +
                            float(new_freq[pi]) * new_freq_coeff[si] +
                            float(freq) * old_freq_coeff[si];
                        _mm256_store_ps(&ampbuf[i + pi * 16 + si * 8], a);
                        _mm256_store_ps(&phasebuf[i + pi * 16 + si * 8], p);
                    }
                    phase = phase + (freq + new_freq[pi]) * 8;
                    freq = new_freq[pi];
                    amp = new_amp[pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
            for (unsigned i = 0; i < nsteps; i += 8) {
                auto v = _mm256_load_ps(&ampbuf[i]) *
                    sinpif_pi(_mm256_load_ps(&phasebuf[i]));
                if (c != 0)
                    v += _mm256_load_ps(&out[i]);
                _mm256_store_ps(&out[i], v);
            }
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
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

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_single(float *out, unsigned nsteps, unsigned nrep,
                        ChnParamFixed param)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    auto inc = _mm512_load_ps(inc_buf);
    double phase = 0;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto phase_v = float(phase) + inc * float(param.freq);
            _mm512_store_ps(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 16;
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                const ChnParamFixed *params, unsigned nparams)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    auto inc = _mm512_load_ps(inc_buf);
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto v = _mm512_set1_ps(0);
            for (unsigned c = 0; c < nparams; c++) {
                auto phase = phases[c];
                auto param = params[c];
                auto phase_v = float(phase) + inc * float(param.freq);
                v += float(param.amp) * sinpif_pi(phase_v);
                phase += param.freq * 16;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
                phases[c] = phase;
            }
            _mm512_store_ps(&out[i], v);
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                  const ChnParamFixed *params, unsigned nparams)
{
    float inc_buf[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    auto inc = _mm512_load_ps(inc_buf);
    double phases[nparams] = {};
    for (unsigned j = 0; j < nrep; j++) {
        auto phase = phases[0];
        auto param = params[0];
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto phase_v = float(phase) + inc * float(param.freq);
            _mm512_store_ps(&out[i], float(param.amp) * sinpif_pi(phase_v));
            phase += param.freq * 16;
            if (phase > 32) {
                // assume positive frequency
                phase -= 64;
            }
        }
        phases[0] = phase;

        for (unsigned c = 1; c < nparams; c++) {
            auto phase = phases[c];
            auto param = params[c];
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto phase_v = float(phase) + inc * float(param.freq);
                _mm512_store_ps(&out[i], _mm512_load_ps(&out[i]) +
                                float(param.amp) * sinpif_pi(phase_v));
                phase += param.freq * 16;
                if (phase > 32) {
                    // assume positive frequency
                    phase -= 64;
                }
            }
            phases[c] = phase;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

static constexpr __m512 new_amp_coeff ={
    0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375,
    0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375};
static constexpr __m512 old_amp_coeff = {
    1.0, 0.9375, 0.875, 0.8125, 0.75, 0.6875, 0.625, 0.5625,
    0.5, 0.4375, 0.375, 0.3125, 0.25, 0.1875, 0.125, 0.0625};
static constexpr __m512 new_freq_coeff = {
    0.0, 0.03125, 0.125, 0.28125, 0.5, 0.78125, 1.125, 1.53125,
    2.0, 2.53125, 3.125, 3.78125, 4.5, 5.28125, 6.125, 7.03125};
static constexpr __m512 old_freq_coeff = {
    0.0, 0.96875, 1.875, 2.71875, 3.5, 4.21875, 4.875, 5.46875,
    6.0, 6.46875, 6.875, 7.21875, 7.5, 7.71875, 7.875, 7.96875};

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_ramp_single(float *out, unsigned nsteps, unsigned nrep,
                             ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m512d pinc{0, 16, 32, 48, 64, 80, 96, 112};

    double phase = 0;
    auto amp = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 128) {
            auto new_amps = ramp_func(double(t) + pinc, amp_param);
            auto new_freqs = ramp_func(double(t) + pinc, freq_param);
            t += 128;
            for (unsigned pi = 0; pi < 8; pi++) {
                auto new_amp = new_amps[pi];
                auto new_freq = new_freqs[pi];
                auto a = float(new_amp) * new_amp_coeff +
                    float(amp) * old_amp_coeff;
                auto p = float(phase) +
                    float(new_freq) * new_freq_coeff +
                    float(freq) * old_freq_coeff;
                _mm512_store_ps(&out[i + pi * 16], a * sinpif_pi(p));
                phase = phase + (freq + new_freq) * 8;
                freq = new_freq;
                amp = new_amp;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_ramp_multi_chn_loop(float *out, unsigned nsteps, unsigned nrep,
                                     const ChnParamMod *amp_params,
                                     const ChnParamMod *freq_params,
                                     unsigned nparams)
{
    __m512d pinc{0, 16, 32, 48, 64, 80, 96, 112};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 128) {
            __m512d new_freqs[nparams];
            __m512d new_amps[nparams];
            for (unsigned c = 0; c < nparams; c++) {
                new_amps[c] = ramp_func(double(t) + pinc, amp_params[c]);
                new_freqs[c] = ramp_func(double(t) + pinc, freq_params[c]);
            }
            t += 128;
            for (unsigned pi = 0; pi < 8; pi++) {
                __m512 v = {0, 0, 0, 0, 0, 0, 0, 0};
                for (unsigned c = 0; c < nparams; c++) {
                    auto a = float(new_amps[c][pi]) * new_amp_coeff +
                        float(amps[c]) * old_amp_coeff;
                    auto p = float(phases[c]) +
                        float(new_freqs[c][pi]) * new_freq_coeff +
                        float(freqs[c]) * old_freq_coeff;
                    v += a * sinpif_pi(p);
                }
                _mm512_store_ps(&out[i + pi * 16], v);
                for (unsigned c = 0; c < nparams; c++) {
                    auto phase = phases[c] + (freqs[c] + new_freqs[c][pi]) * 8;
                    amps[c] = new_amps[c][pi];
                    freqs[c] = new_freqs[c][pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    phases[c] = phase;
                }
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_ramp_multi_chnblk_loop(float *out, unsigned nsteps, unsigned nrep,
                                        const ChnParamMod *amp_params,
                                        const ChnParamMod *freq_params,
                                        unsigned nparams)
{
    __m512d pinc{0, 16, 32, 48, 64, 80, 96, 112};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 128) {
            auto chn_func = [&] (unsigned c)
                __attribute__((target("avx512f,avx512dq"))) {
                auto new_amp = ramp_func(double(t) + pinc, amp_params[c]);
                auto new_freq = ramp_func(double(t) + pinc, freq_params[c]);
                auto amp = amps[c];
                auto freq = freqs[c];
                auto phase = phases[c];
                for (unsigned pi = 0; pi < 8; pi++) {
                    auto a = float(new_amp[pi]) * new_amp_coeff +
                        float(amp) * old_amp_coeff;
                    auto p = float(phase) +
                        float(new_freq[pi]) * new_freq_coeff +
                        float(freq) * old_freq_coeff;
                    auto v = a * sinpif_pi(p);
                    if (c != 0)
                        v += _mm512_load_ps(&out[i + pi * 16]);
                    _mm512_store_ps(&out[i + pi * 16], v);
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
                amps[c] = amp;
                freqs[c] = freq;
                phases[c] = phase;
            };

            // Encourage the compiler to specialize for c=0
            chn_func(0);
            for (unsigned c = 1; c < nparams; c++)
                chn_func(c);
            t += 128;
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_ramp_multi_block_loop(float *out, unsigned nsteps, unsigned nrep,
                                       const ChnParamMod *amp_params,
                                       const ChnParamMod *freq_params,
                                       unsigned nparams)
{
    __m512d pinc{0, 16, 32, 48, 64, 80, 96, 112};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, freq_params[c])[0];
    }
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx512f,avx512dq"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 128) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 128;
                for (unsigned pi = 0; pi < 8; pi++) {
                    auto a = float(new_amp[pi]) * new_amp_coeff +
                        float(amp) * old_amp_coeff;
                    auto p = float(phase) +
                        float(new_freq[pi]) * new_freq_coeff +
                        float(freq) * old_freq_coeff;
                    auto v = a * sinpif_pi(p);
                    if (c != 0)
                        v += _mm512_load_ps(&out[i + pi * 16]);
                    _mm512_store_ps(&out[i + pi * 16], v);
                    phase += (freq + new_freq[pi]) * 8;
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                    amp = new_amp[pi];
                    freq = new_freq[pi];
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_ramp_single_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                  ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m512d pinc{0, 16, 32, 48, 64, 80, 96, 112};

    double phase = 0;
    auto amp = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    double ampbuf[nsteps / 16] __attribute__((aligned(64)));
    double freqbuf[nsteps / 16] __attribute__((aligned(64)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps / 128; i++) {
            _mm512_store_pd(&ampbuf[i * 8], ramp_func(double(t) + pinc, amp_param));
            _mm512_store_pd(&freqbuf[i * 8], ramp_func(double(t) + pinc, freq_param));
            t += 128;
        }
        for (unsigned i = 0; i < nsteps; i += 16) {
            auto new_amp = ampbuf[i / 16];
            auto new_freq = freqbuf[i / 16];
            auto a = float(new_amp) * new_amp_coeff +
                float(amp) * old_amp_coeff;
            auto p = float(phase) +
                float(new_freq) * new_freq_coeff +
                float(freq) * old_freq_coeff;
            _mm512_store_ps(&out[i], a * sinpif_pi(p));
            phase = phase + (freq + new_freq) * 8;
            freq = new_freq;
            amp = new_amp;
            if (phase > 32) {
                phase -= 64;
            }
            else if (phase < -32) {
                phase += 64;
            }
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf(float *out, unsigned nsteps, unsigned nrep,
                                            const ChnParamMod *amp_params,
                                            const ChnParamMod *freq_params,
                                            unsigned nparams)
{
    __m512d pinc{0, 16, 32, 48, 64, 80, 96, 112};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, freq_params[c])[0];
    }
    double ampbuf[nsteps / 16] __attribute__((aligned(64)));
    double freqbuf[nsteps / 16] __attribute__((aligned(64)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx512f,avx512dq"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps / 128; i++) {
                _mm512_store_pd(&ampbuf[i * 8],
                                ramp_func(double(t) + pinc, amp_param));
                _mm512_store_pd(&freqbuf[i * 8],
                                ramp_func(double(t) + pinc, freq_param));
                t += 128;
            }
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto new_amp = ampbuf[i / 16];
                auto new_freq = freqbuf[i / 16];
                auto a = float(new_amp) * new_amp_coeff +
                    float(amp) * old_amp_coeff;
                auto p = float(phase) +
                    float(new_freq) * new_freq_coeff +
                    float(freq) * old_freq_coeff;
                auto v = a * sinpif_pi(p);
                if (c != 0)
                    v += _mm512_load_ps(&out[i]);
                _mm512_store_ps(&out[i], v);
                phase += (freq + new_freq) * 8;
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
                amp = new_amp;
                freq = new_freq;
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_ramp_single_pbuf_full(float *out, unsigned nsteps, unsigned nrep,
                                       ChnParamMod amp_param, ChnParamMod freq_param)
{
    __m512d pinc{0, 16, 32, 48, 64, 80, 96, 112};

    double phase = 0;
    auto amp = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, amp_param)[0];
    auto freq = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, freq_param)[0];
    uint64_t t = 16;
    float ampbuf[nsteps] __attribute__((aligned(64)));
    float phasebuf[nsteps] __attribute__((aligned(64)));
    for (unsigned j = 0; j < nrep; j++) {
        for (unsigned i = 0; i < nsteps; i += 128) {
            auto new_amp = ramp_func(double(t) + pinc, amp_param);
            auto new_freq = ramp_func(double(t) + pinc, freq_param);
            t += 128;
            for (unsigned pi = 0; pi < 8; pi++) {
                auto a = float(new_amp[pi]) * new_amp_coeff +
                    float(amp) * old_amp_coeff;
                auto p = float(phase) +
                    float(new_freq[pi]) * new_freq_coeff +
                    float(freq) * old_freq_coeff;
                _mm512_store_ps(&ampbuf[i + pi * 16], a);
                _mm512_store_ps(&phasebuf[i + pi * 16], p);
                phase = phase + (freq + new_freq[pi]) * 8;
                freq = new_freq[pi];
                amp = new_amp[pi];
                if (phase > 32) {
                    phase -= 64;
                }
                else if (phase < -32) {
                    phase += 64;
                }
            }
        }
        for (unsigned i = 0; i < nsteps; i += 16) {
            _mm512_store_ps(&out[i], _mm512_load_ps(&ampbuf[i]) *
                            sinpif_pi(_mm512_load_ps(&phasebuf[i])));
        }
        asm volatile ("" :: "r"(out) : "memory");
    }
}

NACS_EXPORT() NACS_NOINLINE __attribute__((target("avx512f,avx512dq"),flatten))
void Kernel::sin_ramp_multi_block_loop_pbuf_full(
    float *out, unsigned nsteps, unsigned nrep,
    const ChnParamMod *amp_params, const ChnParamMod *freq_params, unsigned nparams)
{
    __m512d pinc{0, 16, 32, 48, 64, 80, 96, 112};
    double phases[nparams] = {};
    double freqs[nparams];
    double amps[nparams];
    for (unsigned c = 0; c < nparams; c++) {
        amps[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, amp_params[c])[0];
        freqs[c] = ramp_func(__m512d{0, 0, 0, 0, 0, 0, 0, 0}, freq_params[c])[0];
    }
    float ampbuf[nsteps] __attribute__((aligned(64)));
    float phasebuf[nsteps] __attribute__((aligned(64)));
    uint64_t _t = 16;
    for (unsigned j = 0; j < nrep; j++) {
        auto chn_func = [&] (unsigned c) __attribute__((target("avx512f,avx512dq"))) {
            auto t = _t;
            auto amp_param = amp_params[c];
            auto freq_param = freq_params[c];
            auto amp = amps[c];
            auto freq = freqs[c];
            auto phase = phases[c];
            for (unsigned i = 0; i < nsteps; i += 128) {
                auto new_amp = ramp_func(double(t) + pinc, amp_param);
                auto new_freq = ramp_func(double(t) + pinc, freq_param);
                t += 128;
                for (unsigned pi = 0; pi < 8; pi++) {
                    auto a = float(new_amp[pi]) * new_amp_coeff +
                        float(amp) * old_amp_coeff;
                    auto p = float(phase) +
                        float(new_freq[pi]) * new_freq_coeff +
                        float(freq) * old_freq_coeff;
                    _mm512_store_ps(&ampbuf[i + pi * 16], a);
                    _mm512_store_ps(&phasebuf[i + pi * 16], p);
                    phase = phase + (freq + new_freq[pi]) * 8;
                    freq = new_freq[pi];
                    amp = new_amp[pi];
                    if (phase > 32) {
                        phase -= 64;
                    }
                    else if (phase < -32) {
                        phase += 64;
                    }
                }
            }
            amps[c] = amp;
            freqs[c] = freq;
            phases[c] = phase;
            for (unsigned i = 0; i < nsteps; i += 16) {
                auto v = _mm512_load_ps(&ampbuf[i]) *
                    sinpif_pi(_mm512_load_ps(&phasebuf[i]));
                if (c != 0)
                    v += _mm512_load_ps(&out[i]);
                _mm512_store_ps(&out[i], v);
            }
        };
        // Encourage the compiler to specialize for c=0
        chn_func(0);
        for (unsigned c = 1; c < nparams; c++)
            chn_func(c);
        _t += nsteps;
        asm volatile ("" :: "r"(out) : "memory");
    }
}

} // namespace avx512
#endif

}
