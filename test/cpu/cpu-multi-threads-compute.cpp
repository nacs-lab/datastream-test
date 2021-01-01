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

#include <yaml-cpp/yaml.h>

#include <atomic>
#include <chrono>
#include <fstream>
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
    template<bool is_write,bool dummy>
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
            if constexpr (!dummy) {
                m_written->store(is_write, std::memory_order_release);
            }
        }
        ReturnData()
            : m_ptr(nullptr),
              m_written(nullptr)
        {
        }
    private:
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
    template<bool dummy>
    using WrData = ReturnData<true,dummy>;
    template<bool dummy>
    using RdData = ReturnData<false,dummy>;
    BlockRing(T *buff, uint32_t size, uint32_t blksz)
        : m_buff(*buff),
          m_blksz(blksz),
          m_nblk(size / blksz),
          m_meta(*(new BlkData[m_nblk]))
    {
    }
    // BlockRing(BlockRing &&o)
    //     : m_buff(o.m_buff),
    //       m_blksz(o.m_blksz),
    //       m_nblk(o.m_nblk),
    //       m_meta(*(new BlkData[o.m_nblk])),
    //       m_wridx(m_wridx),
    //       m_rdidx(m_rdidx)
    // {
    //     for (uint32_t i = 0; i < m_nblk; i++) {
    //         m_meta[i] = std::move(o.m_meta[i]);
    //     }
    // }
    ~BlockRing()
    {
        // `m_meta` becomes a dangling reference after this and must not be accessed anymore.
        delete[] &m_meta;
    }
    template<bool dummy=false>
    NACS_INLINE WrData<dummy> next_write()
    {
        if constexpr (dummy) {
            auto idx = m_wridx;
            auto ptr = &(&m_buff)[m_blksz * idx];
            idx++;
            if (unlikely(idx >= m_nblk))
                idx = 0;
            m_wridx = idx;
            return WrData<dummy>(ptr, nullptr);
        }
        else {
            auto idx = m_wridx;
            auto &meta = (&m_meta)[idx];
            if (meta.written.load(std::memory_order_acquire))
                return WrData<dummy>();
            auto ptr = &(&m_buff)[m_blksz * idx];
            idx++;
            if (unlikely(idx >= m_nblk))
                idx = 0;
            m_wridx = idx;
            return WrData<dummy>(ptr, &meta.written);
        }
    }
    template<bool dummy=false>
    NACS_INLINE RdData<dummy> next_read()
    {
        if constexpr (dummy) {
            auto idx = m_rdidx;
            auto ptr = &(&m_buff)[m_blksz * idx];
            idx++;
            if (unlikely(idx >= m_nblk))
                idx = 0;
            m_rdidx = idx;
            return RdData<dummy>(ptr, nullptr);
        }
        else {
            auto idx = m_rdidx;
            auto &meta = (&m_meta)[idx];
            if (!meta.written.load(std::memory_order_acquire))
                return RdData<dummy>();
            auto ptr = &(&m_buff)[m_blksz * idx];
            idx++;
            if (unlikely(idx >= m_nblk))
                idx = 0;
            m_rdidx = idx;
            return RdData<dummy>(ptr, &meta.written);
        }
    }
    NACS_INLINE uint32_t blksz() const
    {
        return m_blksz;
    }
    NACS_INLINE uint32_t nblk() const
    {
        return m_nblk;
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

struct DummyType {
    template<typename... Args>
    DummyType(Args&&... args)
    {
    }
};

template<typename Kernel, bool use_localbuff, bool has_in, bool is_final>
class Worker {
    static constexpr uint32_t localbuff_sz = use_localbuff ? 16 * 1024 / 4 : 0;
    static constexpr uint32_t localbuff_blksz = localbuff_sz / 4;
    static BlockRing<float> *const *copy_ins_array(BlockRing<float> *const *ins, int nins)
    {
        // Make sure when there's no input channels, we do not dereference `ins`
        if constexpr (has_in) {
            auto res = new BlockRing<float>*[nins];
            memcpy(res, ins, nins * sizeof(BlockRing<float>*));
            return res;
        }
        else {
            return nullptr;
        }
    }
public:
    Worker(BlockRing<float> &out, BlockRing<float> *const *ins, int nins, int nchn,
           std::atomic<bool> &done)
        : m_out(out),
          m_ins(copy_ins_array(ins, nins)),
          m_nins(nins),
          m_nchn(nchn),
          m_done(done)
    {
        if constexpr (use_localbuff) {
            m_localbuff_fill = 0;
        }
    }
    ~Worker()
    {
        if constexpr (has_in) {
            delete[] m_ins;
        }
    }

    NACS_INLINE void run(float t, float freq, float amp)
    {
        auto blksz = m_out.blksz();
        while (true) {
            auto wrdata = m_out.next_write<is_final>();
            while (!wrdata) {
                if (unlikely(m_done.load(std::memory_order_relaxed)))
                    return;
                backoff(t, freq, amp);
                wrdata = m_out.next_write<is_final>();
            }

            fill_out(wrdata.get(), blksz, t, freq, amp);

            wrdata.done();
            reset();
            if constexpr (is_final) {
                m_blocktime[m_blockcnt++] = Test::cycleclock();
                if (m_blockcnt >= m_blocktime.size()) {
                    m_done.store(true, std::memory_order_relaxed);
                    return;
                }
            }
        }
    }

    void set_maxrun(size_t total)
    {
        if constexpr (is_final) {
            m_blockcnt = 0;
            m_blocktime.resize(total);
        }
    }

    __attribute__((flatten))
    static void threadfun(BlockRing<float> &out, BlockRing<float> *const *ins, int nins,
                          int nchn, std::atomic<bool> &done, float t, float freq, float amp,
                          size_t ntotalblk, std::vector<uint64_t> *blocktime=nullptr)
    {
        Worker worker(out, ins, nins, nchn, done);
        worker.set_maxrun(ntotalblk);
        worker.run(t, freq, amp);
        if constexpr (is_final) {
            if (blocktime) {
                *blocktime = std::move(worker.m_blocktime);
            }
        }
    }

private:
    NACS_INLINE void fill_out_in(float *out, uint32_t blksz, float t, float freq, float amp)
    {
        const float *ins[m_nins + 1];
        decltype(m_ins[0]->next_read()) rddatas[m_nins];
        for (int i = 0; i < m_nins; i++) {
            rddatas[i] = m_ins[i]->next_read();
            while (!rddatas[i]) {
                if (unlikely(m_done.load(std::memory_order_relaxed)))
                    return;
                backoff(t, freq, amp);
                rddatas[i] = m_ins[i]->next_read();
            }
            ins[i] = rddatas[i].get();
        }

        if constexpr (use_localbuff) {
            if (m_localbuff_fill) {
                ins[m_nins] = m_localbuff;
                if (is_final) {
                    Kernel::sum_multi_nt(m_localbuff_fill, m_nins + 1, ins, out);
                }
                else {
                    Kernel::sum_multi(m_localbuff_fill, m_nins + 1, ins, out);
                }
            }
            if (is_final) {
                Kernel::calc_multi_nt(blksz - m_localbuff_fill, m_nchn, m_nins, ins,
                                      out + m_localbuff_fill, t, freq, amp);
            }
            else {
                Kernel::calc_multi(blksz - m_localbuff_fill, m_nchn, m_nins, ins,
                                   out + m_localbuff_fill, t, freq, amp);
            }
        }
        else {
            if (is_final) {
                Kernel::calc_multi_nt(blksz, m_nchn, m_nins, ins, out, t, freq, amp);
            }
            else {
                Kernel::calc_multi(blksz, m_nchn, m_nins, ins, out, t, freq, amp);
            }
        }

        for (auto &rddata: rddatas) {
            rddata.done();
        }
    }
    NACS_INLINE void fill_out_noin(float *out, uint32_t blksz, float t, float freq, float amp)
    {
        if constexpr (use_localbuff) {
            if (is_final) {
                Kernel::copy_nt(m_localbuff_fill, (int*)m_localbuff, (int*)out);
                Kernel::calc_multi_fill_nt(blksz - m_localbuff_fill, m_nchn,
                                           out + m_localbuff_fill, t, freq, amp);
            }
            else {
                Kernel::copy(m_localbuff_fill, (int*)m_localbuff, (int*)out);
                Kernel::calc_multi_fill(blksz - m_localbuff_fill, m_nchn,
                                        out + m_localbuff_fill, t, freq, amp);
            }
        }
        else {
            if (is_final) {
                Kernel::calc_multi_fill_nt(blksz, m_nchn, out, t, freq, amp);
            }
            else {
                Kernel::calc_multi_fill(blksz, m_nchn, out, t, freq, amp);
            }
        }
    }
    NACS_INLINE void fill_out(float *out, uint32_t blksz, float t, float freq, float amp)
    {
        if constexpr (!has_in) {
            fill_out_noin(out, blksz, t, freq, amp);
        }
        else {
            fill_out_in(out, blksz, t, freq, amp);
        }
    }

    NACS_INLINE void backoff(float t, float freq, float amp)
    {
        if constexpr (use_localbuff) {
            if (m_localbuff_fill < localbuff_sz) {
                Kernel::calc_multi_fill(localbuff_blksz, m_nchn,
                                        &m_localbuff[m_localbuff_fill], t, freq, amp);
                m_localbuff_fill += localbuff_blksz;
                return;
            }
        }
        m_backoff.pause();
    }
    NACS_INLINE void reset()
    {
        if constexpr (use_localbuff) {
            m_localbuff_fill = 0;
        }
        m_backoff.wake();
        m_backoff.reset();
    }

    BlockRing<float> &m_out;
    std::conditional_t<has_in,BlockRing<float> *const *const,DummyType> m_ins;
    std::conditional_t<has_in,int,DummyType> m_nins;
    int m_nchn;
    Backoff m_backoff;
    std::conditional_t<use_localbuff,uint32_t,DummyType> m_localbuff_fill;
    std::atomic<bool> &m_done;
    float m_localbuff[localbuff_sz] __attribute__((aligned(64)));
    std::conditional_t<is_final,std::vector<uint64_t>,DummyType> m_blocktime;
    std::conditional_t<is_final,size_t,DummyType> m_blockcnt;
};

template<typename Kernel, bool use_localbuff, bool is_final>
static NACS_INLINE void _threadfun1(BlockRing<float> &out, BlockRing<float> *const *ins,
                                    int nins, int nchn, std::atomic<bool> &done,
                                    float t, float freq, float amp, size_t ntotalblk,
                                    std::vector<uint64_t> *blocktime=nullptr)
{
    if (nins) {
        Worker<Kernel,use_localbuff,true,is_final>::threadfun(out, ins, nins, nchn,
                                                              done, t, freq, amp, ntotalblk,
                                                              blocktime);
    }
    else {
        Worker<Kernel,use_localbuff,false,is_final>::threadfun(out, nullptr, 0, nchn,
                                                               done, t, freq, amp, ntotalblk,
                                                               blocktime);
    }
}

template<typename Kernel, bool use_localbuff>
static NACS_INLINE void _threadfun2(BlockRing<float> &out, BlockRing<float> *const *ins,
                                    int nins, int nchn, std::atomic<bool> &done,
                                    float t, float freq, float amp, size_t ntotalblk,
                                    bool is_final, std::vector<uint64_t> *blocktime=nullptr)
{
    if (is_final) {
        _threadfun1<Kernel,use_localbuff,true>(out, ins, nins, nchn, done,
                                               t, freq, amp, ntotalblk, blocktime);
    }
    else {
        _threadfun1<Kernel,use_localbuff,false>(out, ins, nins, nchn, done,
                                                t, freq, amp, ntotalblk, blocktime);
    }
}

template<typename Kernel>
static void NACS_INLINE threadfun(BlockRing<float> &out, BlockRing<float> *const *ins,
                                  int nins, int nchn, std::atomic<bool> &done,
                                  float t, float freq, float amp, size_t ntotalblk,
                                  bool use_localbuff, bool is_final,
                                  std::vector<uint64_t> *blocktime=nullptr)
{
    if (use_localbuff) {
        _threadfun2<Kernel,true>(out, ins, nins, nchn, done,
                                 t, freq, amp, ntotalblk, is_final, blocktime);
    }
    else {
        _threadfun2<Kernel,false>(out, ins, nins, nchn, done,
                                  t, freq, amp, ntotalblk, is_final, blocktime);
    }
}

struct WorkerConfig {
    std::vector<int> ins;
    int cpu;
    int nchn{1};
    bool use_localbuff{false};
    bool is_final{false};
};

struct Config {
    std::vector<WorkerConfig> workers;
    size_t ntotalblk;
    uint32_t blksz;
    uint32_t nblk;
    uint32_t nblk_final;
    float amp;
    float t;
    float freq;
    static Config loadYAML(const char *fname)
    {
        Config conf;
        auto file = YAML::LoadFile(fname);
        auto required_key = [&] (auto name) {
            if (auto node = file[name])
                return node;
            throw std::runtime_error(std::string("Required key '") + name + "' missing.");
        };
        conf.ntotalblk = required_key("total_block").as<size_t>();
        conf.blksz = required_key("block_size").as<uint32_t>();
        conf.nblk = required_key("num_block").as<uint32_t>();
        conf.nblk_final = required_key("num_block_final").as<uint32_t>();
        auto workers = required_key("workers");
        if (!workers.IsSequence())
            throw std::runtime_error("Invalid workers config.");
        int def_nchn = 1;
        if (auto node = file["num_channel"])
            def_nchn = node.as<int>();
        bool def_localbuff = false;
        if (auto node = file["localbuff"])
            def_localbuff = node.as<bool>();
        std::map<int,int> cpu_map;
        for (const auto &worker: workers) {
            if (!worker.IsMap())
                throw std::runtime_error("Invalid workers config.");
            auto required_key = [&] (auto name) {
                if (auto node = worker[name])
                    return node;
                throw std::runtime_error(std::string("Required worker key '") + name
                                         + "' missing.");
            };
            int worker_id = (int)conf.workers.size();
            auto &wc = conf.workers.emplace_back();
            wc.cpu = required_key("cpu").as<int>();
            if (cpu_map.count(wc.cpu))
                throw std::runtime_error("Multiple workers on the same CPU");
            cpu_map[wc.cpu] = worker_id;
            wc.nchn = def_nchn;
            wc.use_localbuff = def_localbuff;
            if (auto node = worker["num_channel"])
                wc.nchn = node.as<int>();
            if (auto node = worker["localbuff"])
                wc.use_localbuff = node.as<bool>();
            if (auto node = worker["final"])
                wc.is_final = node.as<bool>();
            if (auto ins_node = worker["input_cpu"]) {
                if (!ins_node.IsSequence())
                    throw std::runtime_error("Invalid worker input config.");
                wc.ins = ins_node.as<std::vector<int>>();
            }
        }

        // Map from CPU id to worker ID
        for (auto &wc: conf.workers) {
            for (auto &in: wc.ins) {
                if (!cpu_map.count(in))
                    throw std::runtime_error("Invalid input CPU.");
                in = cpu_map[in];
            }
        }

        conf.t = Gen::rand_single(0, 1);
        conf.freq = Gen::rand_single(0, 1);
        conf.amp = Gen::rand_single(0, 1000);

        return conf;
    }
};

struct WorkerParam {
    BlockRing<float> *out;
    std::vector<BlockRing<float>*> ins;
    int nchn;
    std::atomic<bool> *done;
    float t;
    float freq;
    float amp;
    size_t ntotalblk;
    bool use_localbuff;
    bool is_final;
    template<typename Kernel>
    void NACS_INLINE run(std::vector<uint64_t> *blocktime=nullptr)
    {
        threadfun<Kernel>(*out, &ins[0], (int)ins.size(), nchn, *done, t, freq, amp, ntotalblk,
                          use_localbuff, is_final, blocktime);
    }
};

struct MMapDeleter {
    template<typename T>
    void operator() (T *ptr)
    {
        unmapPage(ptr, size);
    }
    size_t size;
};

struct TestRes {
    std::vector<uint64_t> blocktime;
    void summary(std::ostream &stm)
    {
        stm << (double)blocktime.back() / (double)blocktime.size() << std::endl;
    }
    void save(std::ostream &stm)
    {
        stm.write((const char*)&blocktime[0], blocktime.size() * sizeof(uint64_t));
    }
};

template<typename Kernel>
static TestRes test_threads(const Config &config)
{
    auto nworkers = config.workers.size();
    std::vector<bool> consumed(nworkers, false);
    std::vector<std::unique_ptr<float,MMapDeleter>> buffs(nworkers);
    std::vector<std::unique_ptr<BlockRing<float>>> rings(nworkers);
    std::vector<int> cpus(nworkers);
    if ((config.blksz * 4) % page_size != 0)
        throw std::runtime_error("Block size should be multiple of page size.");
    if (config.ntotalblk < config.nblk_final || config.ntotalblk < config.nblk)
        throw std::runtime_error("Total blocks should fill all the buffers.");
    bool has_final = false;
    for (size_t i = 0; i < nworkers; i++) {
        auto &wc = config.workers[i];
        for (auto in: wc.ins) {
            if (consumed[in])
                throw std::runtime_error("More than one consumer for worker");
            consumed[in] = true;
        }
        if (wc.nchn < 1)
            throw std::runtime_error("No work for worker");
        if (wc.is_final) {
            if (consumed[i])
                throw std::runtime_error("Consumer exists for final worker");
            consumed[i] = true;
            if (has_final)
                throw std::runtime_error("Multiple final output");
            has_final = true;
            auto nele = config.blksz * config.nblk_final;
            auto size = 4 * nele;
            buffs[i] = std::unique_ptr<float,MMapDeleter>((float*)mapAnonPage(size, Prot::RW),
                                                          MMapDeleter{size});
            rings[i] = std::make_unique<BlockRing<float>>(buffs[i].get(), nele, config.blksz);
        }
        else {
            auto nele = config.blksz * config.nblk;
            auto size = 4 * nele;
            buffs[i] = std::unique_ptr<float,MMapDeleter>((float*)mapAnonPage(size, Prot::RW),
                                                          MMapDeleter{size});
            rings[i] = std::make_unique<BlockRing<float>>(buffs[i].get(), nele, config.blksz);
        }
        if (wc.cpu < 0)
            throw std::runtime_error("Invalid CPU ID.");
        cpus[i] = wc.cpu;
    }
    for (auto c: consumed) {
        if (!c) {
            throw std::runtime_error("Unused output from worker.");
        }
    }
    std::vector<WorkerParam> workers(nworkers);
    std::atomic<bool> done{false};
    for (size_t i = 0; i < nworkers; i++) {
        auto &wc = config.workers[i];
        auto &w = workers[i];
        w.out = rings[i].get();
        for (auto in: wc.ins)
            w.ins.push_back(rings[in].get());
        w.nchn = wc.nchn;
        w.done = &done;
        w.t = config.t;
        w.freq = config.freq;
        w.amp = config.amp;
        w.ntotalblk = config.ntotalblk;
        w.use_localbuff = wc.use_localbuff;
        w.is_final = wc.is_final;
    }
    TestRes res;
    std::atomic<int> done_count{0};
    auto t0 = Test::cycleclock();
    Thread::start(cpus, [&] (int i) {
        auto &w = workers[i];
        w.run<Kernel>(w.is_final ? &res.blocktime : nullptr);
        done_count.fetch_add(1, std::memory_order_release);
    });
    while (done_count.load(std::memory_order_acquire) != nworkers) {
        using namespace std::chrono_literals;
        std::this_thread::sleep_for(10ms);
    }
    for (auto &t: res.blocktime)
        t -= t0;
    return res;
}

static TestRes runtests(const Config &config)
{
#if NACS_CPU_X86 || NACS_CPU_X86_64
    if (CPUKernel::hasavx512()) {
        std::cout << "AVX512:" << std::endl;
        return test_threads<avx512::Kernel>(config);
    }
    if (CPUKernel::hasavx2()) {
        std::cout << "AVX2:" << std::endl;
        return test_threads<avx2::Kernel>(config);
    }
    if (CPUKernel::hasavx()) {
        std::cout << "AVX:" << std::endl;
        return test_threads<avx::Kernel>(config);
    }
    std::cout << "SSE2:" << std::endl;
    return test_threads<sse2::Kernel>(config);
#endif

#if NACS_CPU_AARCH64
    std::cout << "ASIMD:" << std::endl;
    return test_threads<asimd::Kernel>(config);
#endif

    std::cout << "Scalar:" << std::endl;
    return test_threads<scalar::Kernel>(config);
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
    auto res = runtests(Config::loadYAML(argv[1]));
    res.summary(std::cout);
    if (argc >= 3) {
        std::ofstream stm(argv[2]);
        res.save(stm);
    }
    return 0;
}
