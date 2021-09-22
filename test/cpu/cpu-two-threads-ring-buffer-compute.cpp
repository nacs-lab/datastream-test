/*************************************************************************
 *   Copyright (c) 2016 - 2020 Yichao Yu <yyc1992@gmail.com>             *
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
#include "helpers/print.h"
#include "helpers/test.h"
#include "helpers/threads.h"

#include <nacs-utils/mem.h>
#include <nacs-utils/number.h>
#include <nacs-utils/thread.h>

#include <atomic>
#include <iostream>

#include <assert.h>
#include <stdint.h>

using namespace NaCs;
using namespace CPUKernel;

static int worker_cpu;
static int nchn_write;
static int nchn_read;

static int localbuff_sz = 0;
static int localbuff_blksz = 0;

#if NACS_CPU_X86 || NACS_CPU_X86_64
struct Backoff {
    NACS_INLINE void wake()
    {
    }
    NACS_INLINE void pause()
    {
        for (uint32_t i = 0; i < count; i++)
            CPU::pause();
        count += count / 4;
    }
    NACS_INLINE void reset()
    {
        count = start_count;
    }
private:
    static constexpr uint32_t start_count = 8;
    uint32_t count = start_count;
};
#else
struct Backoff {
    NACS_INLINE void wake()
    {
        CPU::wake();
    }
    NACS_INLINE void pause()
    {
        CPU::pause();
    }
    NACS_INLINE void reset()
    {
    }
};
#endif

template<typename T>
class BlockRing {
    struct BlkData {
        std::atomic<bool> written{false};
    } __attribute__((aligned(64)));
    template<bool is_write>
    class ReturnData {
        using ptr_t = std::conditional_t<is_write,T*,const T*>;
    public:
        operator bool() const
        {
            return m_ptr != nullptr;
        }
        ptr_t get() const
        {
            return m_ptr;
        }
        operator ptr_t() const
        {
            return get();
        }
        void done()
        {
            m_written->store(is_write, std::memory_order_release);
        }
    private:
        ReturnData()
            : m_ptr(nullptr),
              m_written(nullptr)
        {
        }
        ReturnData(ptr_t ptr, std::atomic<bool> *written)
            : m_ptr(ptr),
              m_written(written)
        {
        }
        ptr_t m_ptr;
        std::atomic<bool> *m_written;
        friend class BlockRing;
    };
public:
    using WrData = ReturnData<true>;
    using RdData = ReturnData<false>;
    BlockRing(T *buff, uint32_t size, uint32_t blksz)
        : m_buff(*buff),
          m_blksz(blksz),
          m_nblk(size / blksz),
          m_meta(*(new BlkData[m_nblk]))
    {
    }
    ~BlockRing()
    {
        // `m_meta` becomes a dangling reference after this and must not be accessed anymore.
        delete[] &m_meta;
    }
    NACS_INLINE WrData next_write()
    {
        auto idx = m_wridx;
        auto &meta = (&m_meta)[idx];
        if (meta.written.load(std::memory_order_acquire))
            return WrData();
        auto ptr = &(&m_buff)[m_blksz * idx];
        idx++;
        if (unlikely(idx >= m_nblk))
            idx = 0;
        m_wridx = idx;
        return WrData(ptr, &meta.written);
    }
    NACS_INLINE RdData next_read()
    {
        auto idx = m_rdidx;
        auto &meta = (&m_meta)[idx];
        if (!meta.written.load(std::memory_order_acquire))
            return RdData();
        auto ptr = &(&m_buff)[m_blksz * idx];
        idx++;
        if (unlikely(idx >= m_nblk))
            idx = 0;
        m_rdidx = idx;
        return RdData(ptr, &meta.written);
    }
    NACS_INLINE uint32_t blksz() const
    {
        return m_blksz;
    }
private:
    // Use references to guarantee the constness of these field to the compiler.
    T &m_buff;
    const uint32_t m_blksz;
    const uint32_t m_nblk;
    BlkData &m_meta;
    uint32_t m_wridx __attribute__((aligned(64))) {0};
    uint32_t m_rdidx __attribute__((aligned(64))) {0};
};

struct BlockCounters {
    uint64_t rw{0};
    uint64_t _syncs{0};
    uint64_t stall{0};
    // std::vector<uint32_t> syncs;
    // std::vector<uint64_t> sync_ts;
    void record(uint32_t nsync)
    {
        // sync_ts.push_back(__builtin_ia32_rdtsc());
        // size_t syncs_sz = syncs.size();
        // if ((syncs_sz & 1) == 0) {
        //     // last element is not zero count
        //     syncs.push_back(0);
        //     syncs_sz += 1;
        // }
        // if (nsync) {
        //     stall++;
        //     syncs.push_back(nsync);
        // }
        // else {
        //     syncs[syncs_sz - 1]++;
        // }
        if (nsync) {
            stall++;
            _syncs += nsync;
        }
        rw++;
    }
    uint64_t sync() const
    {
        // uint64_t res = 0;
        // for (size_t i = 1; i < syncs.size(); i += 2)
        //     res += syncs[i];
        // return res;
        return _syncs;
    }
    // void print_syncs(uint64_t offset=0)
    // {
    //     {
    //         bool first = true;
    //         for (size_t i = 1; i < syncs.size(); i++) {
    //             if ((i & 1) == 0 && syncs[i] == 0)
    //                 continue;
    //             if (first) {
    //                 first = false;
    //             }
    //             else {
    //                 std::cout << " - ";
    //             }
    //             std::cout << syncs[i];
    //             if ((i & 1) == 0) {
    //                 std::cout << "x0";
    //             }
    //         }
    //         std::cout << std::endl;
    //     }
    //     {
    //         bool first = true;
    //         for (auto t: sync_ts) {
    //             if (first) {
    //                 first = false;
    //             }
    //             else {
    //                 std::cout << ", ";
    //             }
    //             std::cout << t - offset;
    //         }
    //         std::cout << std::endl;
    //     }
    // }
};

template<typename Kernel, bool use_localbuff>
static void write_block(BlockRing<float> &ring, size_t nrep, uint32_t nele,
                        float t, float freq, float amp, BlockCounters &counter)
{
    auto blksz = ring.blksz();
    auto nchn = nchn_write;
    uint64_t ntotal = uint64_t(nrep) * nele / blksz;
    float *localbuff = (float*)alignTo((uintptr_t)alloca(localbuff_sz * 4 + 64), 64);
    for (uint64_t i = 0; i < ntotal; i++) {
        Backoff backoff;
        auto wrdata = ring.next_write();
        uint32_t nsync = 0;
        int localbuff_fill = 0;
        while (!wrdata) {
            if (use_localbuff && localbuff_fill < localbuff_sz) {
                Kernel::calc_multi_fill(localbuff_blksz, nchn, &localbuff[localbuff_fill],
                                        t, freq, amp);
                localbuff_fill += localbuff_blksz;
            }
            else {
                nsync++;
                backoff.pause();
            }
            wrdata = ring.next_write();
        }

        if (use_localbuff) {
            Kernel::copy(localbuff_fill, (int*)localbuff, (int*)wrdata.get());
            Kernel::calc_multi_fill(blksz - localbuff_fill, nchn,
                                    wrdata.get() + localbuff_fill, t, freq, amp);
        }
        else {
            Kernel::calc_multi_fill(blksz, nchn, wrdata.get(), t, freq, amp);
        }
        wrdata.done();
        backoff.wake();

        counter.record(nsync);
    }
}

template<typename Kernel, bool use_localbuff>
static void read_block(BlockRing<float> &ring, size_t nrep, uint32_t nele,
                       float t, float freq, float amp, BlockCounters &counter)
{
    auto blksz = ring.blksz();
    auto nchn = nchn_read;
    uint64_t ntotal = uint64_t(nrep) * nele / blksz;
    float *localbuff = (float*)alignTo((uintptr_t)alloca(localbuff_sz * 4 + 64), 64);
    for (uint64_t i = 0; i < ntotal; i++) {
        Backoff backoff;
        auto rddata = ring.next_read();
        uint32_t nsync = 0;
        int localbuff_fill = 0;
        while (!rddata) {
            if (use_localbuff && localbuff_fill < localbuff_sz) {
                Kernel::calc_multi_fill(localbuff_blksz, nchn, &localbuff[localbuff_fill],
                                        t, freq, amp);
                localbuff_fill += localbuff_blksz;
            }
            else {
                nsync++;
                backoff.pause();
            }
            rddata = ring.next_read();
        }

        if (use_localbuff) {
            Kernel::sum(localbuff_fill, localbuff, rddata.get());
            Kernel::read_calc_multi(blksz - localbuff_fill, nchn,
                                    rddata.get() + localbuff_fill, t, freq, amp);
        }
        else {
            Kernel::read_calc_multi(blksz, nchn, rddata.get(), t, freq, amp);
        }
        rddata.done();
        backoff.wake();

        counter.record(nsync);
    }
}

template<typename Kernel, bool localbuff>
static void test_block(size_t nrep, uint32_t nele, uint32_t block_size)
{
    if (nrep < 128)
        nrep = 128;
    float t = Gen::rand_single(0, 1);
    float freq = Gen::rand_single(0, 1);
    float amp = Gen::rand_single(0, 1000);
    auto buff = (float*)mapAnonPage(alignTo(nele * 4, page_size), Prot::RW);
    BlockRing<float> ring(buff, nele, block_size);
    std::map<std::string,double> read_perf;
    std::atomic<bool> done{false};
    BlockCounters rd_counter;
    BlockCounters wr_counter;
    // rd_counter.syncs.reserve(nrep * nele / ring.blksz());
    // wr_counter.syncs.reserve(nrep * nele / ring.blksz());
    // rd_counter.sync_ts.reserve(nrep * nele / ring.blksz());
    // wr_counter.sync_ts.reserve(nrep * nele / ring.blksz());
    ::Thread::start(std::vector<int>{worker_cpu}, [=, &ring, &read_perf, &done,
                                                 &rd_counter] (int) {
        Test::Timer timer;
        timer.enable_cache();
        timer.restart();
        read_block<Kernel,localbuff>(ring, nrep, nele, t, freq, amp, rd_counter);
        read_perf = timer.get_res(nrep, nele);
        read_perf["pipe_rw"] = (double)rd_counter.rw / (double)nrep / (double)nele;
        read_perf["pipe_sync"] = (double)rd_counter.sync() / (double)nrep / (double)nele;
        read_perf["pipe_stall"] = (double)rd_counter.stall / (double)nrep / (double)nele;
        done.store(true, std::memory_order_release);
        CPU::wake();
    });

    Test::Timer timer;
    timer.enable_cache();
    timer.restart();
    write_block<Kernel,localbuff>(ring, nrep, nele, t, freq, amp, wr_counter);
    auto write_perf = timer.get_res(nrep, nele);
    write_perf["pipe_rw"] = (double)wr_counter.rw / (double)nrep / (double)nele;
    write_perf["pipe_sync"] = (double)wr_counter.sync() / (double)nrep / (double)nele;
    write_perf["pipe_stall"] = (double)wr_counter.stall / (double)nrep / (double)nele;
    while (!done.load(std::memory_order_acquire))
        CPU::pause();

    // uint64_t toffset = min(rd_counter.sync_ts[0], wr_counter.sync_ts[0]);
    // rd_counter.print_syncs(toffset);
    // wr_counter.print_syncs(toffset);

    std::map<std::string,decltype(read_perf)> total_perf{
        std::make_pair("read", read_perf),
        std::make_pair("write", write_perf),
    };
    Print::json(std::cout, total_perf);
    std::cout << std::endl;
}

template<bool localbuff>
static void runtests(uint32_t size, uint32_t block_size)
{
    assert((size % block_size) == 0);
    size_t base_size = size_t(8ull * 1024 * 1024 * 1024 / size / max(nchn_read, nchn_write));
#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cout << "AVX512:" << std::endl;
        test_block<avx512::Kernel,localbuff>(32 * base_size, size / 4, block_size / 4);
        return;
    }
    if (CPUKernel::hasavx2()) {
        std::cout << "AVX2:" << std::endl;
        test_block<avx2::Kernel,localbuff>(16 * base_size, size / 4, block_size / 4);
        return;
    }
    if (CPUKernel::hasavx()) {
        std::cout << "AVX:" << std::endl;
        test_block<avx::Kernel,localbuff>(12 * base_size, size / 4, block_size / 4);
        return;
    }
    std::cout << "SSE2:" << std::endl;
    test_block<sse2::Kernel,localbuff>(6 * base_size, size / 4, block_size / 4);
    return;
#endif

#if NACS_CPU_AARCH64
    std::cout << "ASIMD:" << std::endl;
    test_block<asimd::Kernel,localbuff>(6 * base_size, size / 4, block_size / 4);
    return;
#endif

    std::cout << "Scalar:" << std::endl;
    test_block<scalar::Kernel,localbuff>(base_size, size / 4, block_size / 4);
}

static inline long parse_int(const char *s)
{
    if (!*s) {
        fprintf(stderr, "Error parsing integer '%s'\n", s);
        exit(1);
    }
    char *ep;
    long res = strtol(s, &ep, 10);
    if (s == ep || *ep != 0) {
        fprintf(stderr, "Error parsing integer '%s'\n", s);
        exit(1);
    }
    return res;
}

int main(int argc, char **argv)
{
    if (argc < 8) {
        fprintf(stderr, "Needs at least seven arguments\n");
        exit(1);
    }
    uint32_t size = (uint32_t)parse_int(argv[1]);
    uint32_t block_size = (uint32_t)parse_int(argv[2]);
    nchn_write = (int)parse_int(argv[3]);
    nchn_read = (int)parse_int(argv[4]);
    localbuff_sz = (int)parse_int(argv[5]) / 4;
    localbuff_blksz = (int)parse_int(argv[6]) / 4;
    int cpu0 = (int)parse_int(argv[7]);
    ::Thread::pin(cpu0);
    if (argc >= 4) {
        worker_cpu = (int)parse_int(argv[8]);
    }
    else if (cpu0 == 1) {
        worker_cpu = 0;
    }
    else {
        worker_cpu = 1;
    }
    if (localbuff_sz > 0) {
        runtests<true>(size, block_size);
    }
    else {
        runtests<false>(size, block_size);
    }
    return 0;
}
