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
#include <nacs-utils/processor.h>

namespace CPUKernel {

NACS_EXPORT(ds_helper) extern const NaCs::CPUInfo &host;
NACS_EXPORT(ds_helper) bool hasavx();
NACS_EXPORT(ds_helper) bool hasavx2();
NACS_EXPORT(ds_helper) bool hasavx512();

namespace scalar {

struct Kernel {
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
};

} // namespace scalar

#if NACS_CPU_AARCH64

namespace asimd {

struct Kernel {
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
};

} // namespace asimd
#endif

#if NACS_CPU_X86 || NACS_CPU_X86_64

namespace sse2 {

struct Kernel {
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
};

} // namespace sse2

namespace avx {

struct Kernel {
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
};

} // namespace avx

namespace avx2 {

struct Kernel {
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
};

} // namespace avx2

namespace avx512 {

struct Kernel {
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
};

} // namespace avx512
#endif

}

#endif
