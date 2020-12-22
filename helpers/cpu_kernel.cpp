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
    calc_fill_nt(nrep, ncalc, buff, t, freq, amp);
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
    calc_fill_nt(nrep, ncalc, buff, t, freq, amp);
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

} // namespace avx512
#endif

}
