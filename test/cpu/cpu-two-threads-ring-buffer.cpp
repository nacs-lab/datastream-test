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
    static constexpr uint32_t start_count = 16;
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
    template<bool doneval>
    struct ReturnData {
        operator bool() const
        {
            return m_ptr != nullptr;
        }
        T *get() const
        {
            return m_ptr;
        }
        operator T*() const
        {
            return get();
        }
        void done()
        {
            m_written->store(doneval, std::memory_order_release);
        }
    private:
        ReturnData()
            : m_ptr(nullptr),
              m_written(nullptr)
        {
        }
        ReturnData(T *ptr, std::atomic<bool> *written)
            : m_ptr(ptr),
              m_written(written)
        {
        }
        T *m_ptr;
        std::atomic<bool> *m_written;
        friend class BlockRing;
    };
public:
    using WrData = struct ReturnData<true>;
    using RdData = struct ReturnData<false>;
    BlockRing(T *buff, size_t size, size_t blksz)
        : m_buff(buff),
          m_blksz(blksz),
          m_meta(size / blksz)
    {
    }
    NACS_INLINE WrData next_write()
    {
        auto idx = m_wridx;
        auto &meta = m_meta[idx];
        if (meta.written.load(std::memory_order_acquire))
            return WrData();
        auto ptr = &m_buff[m_blksz * idx];
        idx++;
        if (unlikely(idx >= m_meta.size()))
            idx = 0;
        m_wridx = idx;
        return WrData(ptr, &meta.written);
    }
    NACS_INLINE RdData next_read()
    {
        auto idx = m_rdidx;
        auto &meta = m_meta[idx];
        if (!meta.written.load(std::memory_order_acquire))
            return RdData();
        auto ptr = &m_buff[m_blksz * idx];
        idx++;
        if (unlikely(idx >= m_meta.size()))
            idx = 0;
        m_rdidx = idx;
        return RdData(ptr, &meta.written);
    }
    NACS_INLINE size_t blksz() const
    {
        return m_blksz;
    }
private:
    T *const m_buff;
    const size_t m_blksz;
    std::vector<BlkData> m_meta;
    size_t m_wridx __attribute__((aligned(64))) {0};
    size_t m_rdidx __attribute__((aligned(64))) {0};
};

struct BlockCounters {
    uint64_t rw{0};
    uint64_t sync{0};
    uint64_t stall{0};
};

template<typename Kernel>
static void write_block(BlockRing<int> &ring, size_t nrep, size_t nele, int v,
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

        Kernel::fill1(blksz, wrdata.get(), v);
        wrdata.done();
        backoff.wake();

        if (nsync) {
            counter.stall++;
            counter.sync += nsync;
        }
        counter.rw++;
    }
}

template<typename Kernel>
static void read_block(BlockRing<int> &ring, size_t nrep, size_t nele,
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

        Kernel::read1(blksz, rddata.get());
        rddata.done();
        backoff.wake();

        if (nsync) {
            counter.stall++;
            counter.sync += nsync;
        }
        counter.rw++;
    }
}

template<typename Kernel>
static void test_block(size_t nrep, size_t nele)
{
    if (nrep < 128)
        nrep = 128;
    Thread::pin(0);
    int v = (int)Gen::rand_single(10, 1000000);
    auto buff = (int*)mapAnonPage(alignTo(nele * 4, page_size), Prot::RW);
    BlockRing<int> ring(buff, nele, nele / 512);
    std::map<std::string,double> read_perf;
    std::atomic<bool> done{false};
    Thread::start(std::vector<int>{1}, [=, &ring, &read_perf, &done] (int) {
        Test::Timer timer;
        timer.enable_cache();
        timer.restart();
        BlockCounters counter;
        read_block<Kernel>(ring, nrep, nele, counter);
        read_perf = timer.get_res(nrep, nele);
        read_perf["ring_rw"] = (double)counter.rw / (double)nrep / (double)nele;
        read_perf["ring_sync"] = (double)counter.sync / (double)nrep / (double)nele;
        read_perf["ring_stall"] = (double)counter.stall / (double)nrep / (double)nele;
        done.store(true, std::memory_order_release);
        CPU::wake();
    });

    Test::Timer timer;
    timer.enable_cache();
    timer.restart();
    BlockCounters counter;
    write_block<Kernel>(ring, nrep, nele, v, counter);
    auto write_perf = timer.get_res(nrep, nele);
    write_perf["ring_rw"] = (double)counter.rw / (double)nrep / (double)nele;
    write_perf["ring_sync"] = (double)counter.sync / (double)nrep / (double)nele;
    write_perf["ring_stall"] = (double)counter.stall / (double)nrep / (double)nele;
    while (!done.load(std::memory_order_acquire))
        CPU::pause();

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
static void write_pipe(DataPipe<int> &pipe, size_t nrep, size_t nele, int v,
                       PipeCounters &counter)
{
    uint64_t ntotal = uint64_t(nrep) * nele;
    assert(nele > 0);
    assert(nele % 16 == 0);
    size_t nmin = nele / 16;
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
            if (sz > nmin * 2)
                sz /= 2;
            sz &= ~(size_t)15;
            if (sz != 0)
                break;
            nsync++;
            backoff.pause();
        }
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
static void read_pipe(DataPipe<int> &pipe, size_t nrep, size_t nele,
                      PipeCounters &counter)
{
    uint64_t ntotal = uint64_t(nrep) * nele;
    assert(nele > 0);
    assert(nele % 16 == 0);
    size_t nmin = nele / 16;
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
            if (sz > nmin * 2)
                sz /= 2;
            sz &= ~(size_t)15;
            if (sz != 0)
                break;
            nsync++;
            backoff.pause();
        }
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
static void test_pipe(size_t nrep, size_t nele)
{
    if (nrep < 128)
        nrep = 128;
    Thread::pin(0);
    int v = (int)Gen::rand_single(10, 1000000);
    auto buff = (int*)mapAnonPage(alignTo(nele * 4, page_size), Prot::RW);
    DataPipe<int> pipe(buff, nele);
    std::map<std::string,double> read_perf;
    std::atomic<bool> done{false};
    Thread::start(std::vector<int>{1}, [=, &pipe, &read_perf, &done] (int) {
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

static void runtests(long size)
{
#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cout << "AVX512:" << std::endl;
        test_pipe<avx512::Kernel>(size_t(32 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
        test_block<avx512::Kernel>(size_t(32 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
        return;
    }
    if (CPUKernel::hasavx2()) {
        std::cout << "AVX2:" << std::endl;
        test_pipe<avx2::Kernel>(size_t(16 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
        test_block<avx2::Kernel>(size_t(16 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
        return;
    }
    if (CPUKernel::hasavx()) {
        std::cout << "AVX:" << std::endl;
        test_pipe<avx::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
        test_block<avx::Kernel>(size_t(12 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
        return;
    }
    std::cout << "SSE2:" << std::endl;
    test_pipe<sse2::Kernel>(size_t(6 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
    test_block<sse2::Kernel>(size_t(6 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
    return;
#endif

#if NACS_CPU_AARCH64
    std::cout << "ASIMD:" << std::endl;
    test_pipe<asimd::Kernel>(size_t(6 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
    test_block<asimd::Kernel>(size_t(6 * 16ull * 1024 * 1024 * 1024 / size), size / 4);
    return;
#endif

    std::cout << "Scalar:" << std::endl;
    test_pipe<scalar::Kernel>(size_t(16ull * 1024 * 1024 * 1024 / size), size / 4);
    test_block<scalar::Kernel>(size_t(16ull * 1024 * 1024 * 1024 / size), size / 4);
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
    if (argc < 2) {
        fprintf(stderr, "Needs at least one argument\n");
        exit(1);
    }
    runtests(parse_int(argv[1]));
    return 0;
}
