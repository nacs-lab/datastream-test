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

template<typename Kernel>
static void write_block(BlockRing<int> &ring, size_t nrep, uint32_t nele, int v,
                        BlockCounters &counter)
{
    auto blksz = ring.blksz();
    uint64_t ntotal = uint64_t(nrep) * nele / blksz;
    for (uint64_t i = 0; i < ntotal; i++) {
        Backoff backoff;
        auto wrdata = ring.next_write();
        uint32_t nsync = 0;
        while (!wrdata) {
            nsync++;
            backoff.pause();
            wrdata = ring.next_write();
        }

        if (!Test::empty)
            Kernel::fill1(blksz, wrdata.get(), v);
        wrdata.done();
        backoff.wake();

        counter.record(nsync);
    }
}

template<typename Kernel>
static void read_block(BlockRing<int> &ring, size_t nrep, uint32_t nele,
                       BlockCounters &counter)
{
    auto blksz = ring.blksz();
    uint64_t ntotal = uint64_t(nrep) * nele / blksz;
    for (uint64_t i = 0; i < ntotal; i++) {
        Backoff backoff;
        auto rddata = ring.next_read();
        uint32_t nsync = 0;
        while (!rddata) {
            nsync++;
            backoff.pause();
            rddata = ring.next_read();
        }

        if (!Test::empty)
            Kernel::read1(blksz, rddata.get());
        rddata.done();
        backoff.wake();

        counter.record(nsync);
    }
}

template<typename Kernel>
static void test_block(size_t nrep, uint32_t nele, uint32_t block_size)
{
    if (nrep < 128)
        nrep = 128;
    int v = (int)Gen::rand_single(10, 1000000);
    auto buff = (int*)mapAnonPage(alignTo(nele * 4, page_size), Prot::RW);
    BlockRing<int> ring(buff, nele, block_size);
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
        read_block<Kernel>(ring, nrep, nele, rd_counter);
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
    write_block<Kernel>(ring, nrep, nele, v, wr_counter);
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

struct PipeCounters {
    uint64_t rw{0};
    uint64_t sync{0};
    uint64_t stall{0};
};

template<typename Kernel>
static void write_pipe(DataPipe<int> &pipe, size_t nrep, uint32_t nele, int v,
                       PipeCounters &counter)
{
    uint64_t ntotal = uint64_t(nrep) * nele;
    assert(nele > 0);
    assert(nele % 16 == 0);
    for (uint64_t i = 0; i < ntotal; ) {
        size_t sz;
        int *ptr;
        uint32_t nsync = 0;
        Backoff backoff;
        while (true) {
            ptr = pipe.get_write_ptr(&sz);
            if (sz <= 15 && sz > 0) {
                pipe.sync_writer();
            }
            sz &= ~(size_t)15;
            if (sz != 0)
                break;
            nsync++;
            backoff.pause();
        }
        if (!Test::empty)
            Kernel::fill1(sz, ptr, v);
        if (nsync) {
            counter.stall++;
            counter.sync += nsync;
        }
        counter.rw++;
        // This can produce more data then n
        pipe.wrote_size(sz);
        backoff.wake();
        i += sz;
    }
}

template<typename Kernel>
static void read_pipe(DataPipe<int> &pipe, size_t nrep, uint32_t nele,
                      PipeCounters &counter)
{
    uint64_t ntotal = uint64_t(nrep) * nele;
    assert(nele > 0);
    assert(nele % 16 == 0);
    for (uint64_t i = 0; i < ntotal; ) {
        size_t sz;
        const int *ptr;
        uint32_t nsync = 0;
        Backoff backoff;
        while (true) {
            ptr = pipe.get_read_ptr(&sz);
            if (sz <= 15 && sz > 0) {
                pipe.sync_reader();
            }
            sz &= ~(size_t)15;
            if (sz != 0)
                break;
            nsync++;
            backoff.pause();
        }
        if (!Test::empty)
            Kernel::read1(sz, ptr);
        if (nsync) {
            counter.stall++;
            counter.sync += nsync;
        }
        counter.rw++;
        pipe.read_size(sz);
        backoff.wake();
        i += sz;
    }
}

template<typename Kernel>
static void test_pipe(size_t nrep, uint32_t nele, uint32_t block_size)
{
    if (nrep < 128)
        nrep = 128;
    int v = (int)Gen::rand_single(10, 1000000);
    auto buff = (int*)mapAnonPage(alignTo(nele * 4, page_size), Prot::RW);
    DataPipe<int> pipe(buff, nele, block_size);
    std::map<std::string,double> read_perf;
    std::atomic<bool> done{false};
    ::Thread::start(std::vector<int>{worker_cpu}, [=, &pipe, &read_perf, &done] (int) {
        Test::Timer timer;
        timer.enable_cache();
        timer.restart();
        PipeCounters counter;
        read_pipe<Kernel>(pipe, nrep, nele, counter);
        read_perf = timer.get_res(nrep, nele);
        read_perf["pipe_rw"] = (double)counter.rw / (double)nrep / (double)nele;
        read_perf["pipe_sync"] = (double)counter.sync / (double)nrep / (double)nele;
        read_perf["pipe_stall"] = (double)counter.stall / (double)nrep / (double)nele;
        done.store(true, std::memory_order_release);
        CPU::wake();
    });

    Test::Timer timer;
    timer.enable_cache();
    timer.restart();
    PipeCounters counter;
    write_pipe<Kernel>(pipe, nrep, nele, v, counter);
    auto write_perf = timer.get_res(nrep, nele);
    write_perf["pipe_rw"] = (double)counter.rw / (double)nrep / (double)nele;
    write_perf["pipe_sync"] = (double)counter.sync / (double)nrep / (double)nele;
    write_perf["pipe_stall"] = (double)counter.stall / (double)nrep / (double)nele;
    while (!done.load(std::memory_order_acquire))
        CPU::pause();

    std::map<std::string,decltype(read_perf)> total_perf{
        std::make_pair("read", read_perf),
        std::make_pair("write", write_perf),
    };
    Print::json(std::cout, total_perf);
    std::cout << std::endl;
}

static void runtests(uint32_t size, uint32_t block_size)
{
    assert((size % block_size) == 0);
#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cout << "AVX512:" << std::endl;
        test_pipe<avx512::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                                  size / 4, block_size / 4);
        test_block<avx512::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                                   size / 4, block_size / 4);
        return;
    }
    if (CPUKernel::hasavx2()) {
        std::cout << "AVX2:" << std::endl;
        test_pipe<avx2::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                                size / 4, block_size / 4);
        test_block<avx2::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                                 size / 4, block_size / 4);
        return;
    }
    if (CPUKernel::hasavx()) {
        std::cout << "AVX:" << std::endl;
        test_pipe<avx::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                               size / 4, block_size / 4);
        test_block<avx::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                                size / 4, block_size / 4);
        return;
    }
    std::cout << "SSE2:" << std::endl;
    test_pipe<sse2::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                            size / 4, block_size / 4);
    test_block<sse2::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                             size / 4, block_size / 4);
    return;
#endif

#if NACS_CPU_AARCH64
    std::cout << "ASIMD:" << std::endl;
    test_pipe<asimd::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                             size / 4, block_size / 4);
    test_block<asimd::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                              size / 4, block_size / 4);
    return;
#endif

    std::cout << "Scalar:" << std::endl;
    test_pipe<scalar::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                              size / 4, block_size / 4);
    test_block<scalar::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size),
                               size / 4, block_size / 4);
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
    if (argc < 4) {
        fprintf(stderr, "Needs at least three arguments\n");
        exit(1);
    }
    int cpu0 = (int)parse_int(argv[3]);
    ::Thread::pin(cpu0);
    if (argc >= 4) {
        worker_cpu = (int)parse_int(argv[4]);
    }
    else if (cpu0 == 1) {
        worker_cpu = 0;
    }
    else {
        worker_cpu = 1;
    }
    runtests((uint32_t)parse_int(argv[1]), (uint32_t)parse_int(argv[2]));
    return 0;
}
