/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef HELPERS_CHANNEL_WORKER_H
#define HELPERS_CHANNEL_WORKER_H

#include "cpu_kernel.h"

#include <nacs-utils/utils.h>
#include <nacs-utils/processor.h>

namespace CPUKernel {

// # Sampling rates
//
// The parameters of the sine wave (amplitude and frequency) is evaluated
// once per 16 output samples. This number is selected based on the vector size for avx512
// but should be used independent of the kernel implementation so that the ramp resolution
// is independent of the hardware used. This gives a parameter update rate of 39.0625 MHz
// for a 625 MHz sampling rate, or a time resolution of 25.6 ns which should be enough for now.

// # Value scaling and rounding
// We round the frequencies to integer values so that the user can select exact frequencies
// more easily when dealing with interference between multiple tones.
// The phase values are similarly rounded to integers to achieve accurate phase accumulations.

// The integer frequency value is computed by scaling and rounding the return values
// from the ramp functions is such that 16 output samples at a frequency number `f`
// corresponds to a phase accumulation of `2 * f`.
// (Similarly, a linear frequency ramp from `f1` to `f2` within 16 output samples
//  corresponds to a phase accumulation of `f1 + f2`.)
// The integer phase value will be scaled by a constant (passed in as template parameter
// in this test implementation) to get the real phase for sine function evaluation.
// The phase wrap around value (i.e. a value corresponds to `2pi * n`
// when the phase number is reset) is also passed in as the template parameter.

// ## Phase wrapping
// The maximum valid frequency value corresponds to `pi` phase difference on each sample.
// We would like to do phase computation every 16 output samples
// (to match when we need to potentially update the amplitude and frequency values,
//  so that we don't introduce unnecessary branches),
// which means that every time we would accumulate up to `16pi` (8 cycles) phase.
// We would like to do phase reset in a way that doesn't require division
// and using addition or subtraction instead.
// If we pick `16pi` as the maximum phase and subtract `32pi` of phase
// everytime we run over the limit, we can make sure that we'll never have any phase runaway
// (i.e. continuously increasing phase over time)
// while also keeping the phase always within `[-16pi, 32pi]`.
// This corresponds to losing up to 4 bits (16 cycles) of precision
// or 20 bits left in the single percision format which should still be good enough
// since our output percision is 16bit.

// # Stream / channel command
// ## Branch elimination/reduction
// For the real generation, we cannot hard code everything directly
// since the timing and ordering of the pulses are not necessarily known at compile time.
// Some dispatch at runtime is definitely necessary. We can, however,
// reduce the total number of dispatch by doing most of the dispatch in one/few go.
// The dispatches/branches we need are,
//
// 1. Which channel to compute
// 2. Whether the channel has ramp
// 3. If the channel ramp is on frequency or amplitude (or both)
// 4. The type of ramp function (scalar or vector) and the call of the function itself.
// 5. (For vector) whether there are cached value for the parameter
// 6. Whether/when we need to finish the current command and/or forward to the next command.
// 7. How many samples to compute
//
// Out of these, the ones that we basically can't avoid are
//
// 1. Sample iteration loop (how many samples to compute)
// 2. Channel loop (which channel(s) to compute).
//    This can in principle be skipped when there's only one active channel.
// 3. Channel ramp dispatch
//
// Out of the ones above, 2 and 3 can in principle be avoided but doing so requires
// hardcoding the ramp functions so we can't really do that ahead of time.
// The best we could do at this point is to do all the branches/dispatches we need
// in one of these three necessary branches.
//
// 1. For the sample iteration loop.
//    In addition to/instead of setting the iteration limit to the buffer size
//    or when we want to pass the data to the consumer,
//    we can also set the limit based on when we need to process the next command
//    or finish the current ramp. This limit can also be pre-computed
//    when computing the command stream.
// 2. For the channel loop, I don't think there's much we can do
//    since this is mostly a constant loop within each block of data generation.
// 3. For each channel during each iteration, we need to figure out what ramp, if any,
//    we are doing for each channel. Since we need to be able to handle arbitrary ramp
//    this is equivalent to calling an arbitrary function (though we may inline the function
//    and essentially turn the function call into a switch) and we should be able to combine
//    most of what we need to branch/dispatch on into this function call
//    (effectively doing jump threading) by specialing the code that computes the sine functions
//    on the ramp functions. This includes
//
//    1. Ramp vs constant
//    2. Amp ramp vs freq ramp vs ramp on both
//    3. Vector ramp vs scalar ramp
//
//    Since we can't necessarily match amp and freq ramps
//    and we don't really want to have N^2 specializations by creating a specialization
//    for each combination of amp and freq ramp functions,
//    the function to ramp both amplitude and frequency at the same time will most likely
//    need to have an additional dispatch for either the amp or freq ramp function.
//    (though we can most likely still specialize it for scalar vs vector ramp)
//    The check for cached ramping parameter for vector ramp functions will most likely
//    still need to be in the ramp function since it is different on each call to the function
//    and may not be synchronized between channels.
//    (Note that we could potentially eliminate this check by letting each function
//     compute all the samples for the ramp parameters it computed
//     but this would require passing the result of the sine wave in memory therefore
//     increasing memory access, which seems to negatively impact performance
//     even if it fits fully into L1 cache. It will also require additional loop end checking
//     in the specialized ramp function (since we are now running a loop over multiple samples
//     and may need to stop if the current ramp ends early), offsetting the branch elimination
//     advantage. This may be useful for GPU version using local memory though.)
//
// This specialization means that the "command" seen by each stream of data
// will include raw commands from multiple channels and these stream command will
// not have any overlap between them. In another word, the stream commands
// will basically be a list of pre-processed markers marking
// the start and end time for each raw commands and due to command merging
// the end of a raw command may correspond to the start of a new stream command
// (e.g. the old command being one that ramps freq and amp at the same time whereas
//  the new one only ramps amp).

// From this point on, we'll use "raw command" to refer to the sequence level pulses
// from the user. These applies to a single parameter (amplitude, frequency or phase)
// on a single output channel. "Channel command" will mean the processed version of this
// specific to a single channel that is passed to each stream. Compared to the raw command,
// in addition to having different structure layout to store the information,
// the ramp functions are replaced by the specialized sine wave computation functions
// mentioned above. Finally, "stream command" will be used to refer to all the command
// that are passed to a stream. Among other things, the stream command will
// have the channel ID on top of the channel command.
// It'll also include command to add/remove streams, which would obviously be required
// if/when we support dynamically moving channels between streams,
// but would also be required at the beginning of a new sequence,
// even if we only do static assignment of channel to streams.

// ## Command format
// Among other formation translation details, the merging of amplitude and frequency ramps
// on a channel into a single channel command means that the time when a new channel
// command start may not be the same as the time offset we need to convert sequence time
// to pulse (i.e. ramp function input) time.
// Since the time offset is also in general not known at compile time
// we can't store this info in the code either and have to put it in the channel command.

// We expect the channel command (in fact, all commands) to be issued relatively rarely
// compared to the sampling rate and for most, if not all, of the commands
// to be generated before the sequence starts so the memory bandwidth on the command stream
// shouldn't be very important. It's more important to have a command format that's
// as uniform (less branches) and simple to decode as possible.
// (Unfortunately, due to the requirement for setting some end parameters
//  we can't specialize all branches in the sine wave computation function themselves
//  and can't avoid branches in decoding the commands.)

// The information we need to store in the channel command includes
//
// * Function (index or pointer)
// * Length (time)
// * Ramp time offset (amplitude and/or frequency), required if there's ramp
// * End values (amplitude and/or frequency or phase),
//   required if the ramp ends by the end of the command.
//   (We need this in case the ramp end doesn't align with the end of a 16-sample block)
//    in which case we would no have evaluated the end value
//    if we only evaluate the ramp function on the grid.
// * And we also need some metadata to record the type of command
//
// For other stream commands, like adding or removing channels, or triggering,
// we'll use a structure of the same size.

// # Channel numbering
// Within each stream, we use an ID number to identify the channel.
// For maximum performance during sine wave generation, we'd like to store the current information
// of the channels in a compat array.
// For maximum performance during command processing on the stream,
// we will use the number to index directly into this array.
// This means that as we add and remove channels from the stream,
// the channel ID used on the stream, and therefore the ID in the incoming commands,
// will change, so the code generating the commands must be aware of the mapping
// between the real channel and the ID used on a stream.
// Because of this, we might as well let that code fully manage the channels.

// It might also be possible to let the managing code to also manage the state buffer
// completely as well, i.e. instead of passing commands to the channels, simply
// pass the resulting buffer storing the states of the channels to the workers
// and each (group of) command would correspond to passing a new buffer to the worker.
// This could further simplify the code on the worker, making it easier
// to implement in LLVM IR (so that it can be included in the generated worker function)
// and to support GPU (also easier codegen and less need to do pointer chasing/complex branches).
// It should also make it easier to move channels around
// since the workers would be completely agnostic to channel movement.
// However, doing so fully would require waiting for the computation to finish
// before the managing thread can generate the state buffer and it might be more efficient
// to just let the worker thread to do the channel state management work
// to minimize thread communication (and if the worker thread needs to wait for
// the managing thread before continuing the main thread should have also
// issued some other work for the worker to do in the mean time, meaning that
// we'll need to create more logical streams to keep the workers busy
// and that will increase memory access...).
// The trade-off on GPU might be different,
// but I suspect we still need to do some limited state processing in the GPU code...

enum class CmdType {
    // Normal channel command taking care of amplitude and frequency
    Channel,
    // Set phase
    Phase,
    // Add channel (initial parameter taken from a global array)
    Add,
    // Remove channel (final parameter stored to a global array)
    Remove,
    // Sync channel (current parameter stored to a global array)
    Sync,
    // Move channel (override channel at the target)
    Move,
};

struct StreamCmd {
    uint32_t chn;
    CmdType type;
    uint8_t has_amp_end: 1;
    uint8_t has_freq_end: 1;
    union {
        uint32_t func;
        uint32_t global_idx; // for add and remove command
    };
    uint64_t end_time;
    uint64_t amp_start_time;
    uint64_t freq_start_time;
    union {
        // Phase setting is so infrequent so I really don't want to
        // increase the size of every command by that much just for this...
        int64_t phase;
        struct {
            float amp;
            uint32_t freq;
        };
    } end;
};

struct ChnState {
    // By the end of a 16-output-sample block, the numbers stored in here should always
    // be the ones that should be used at the end of the block.
    // Each implementation can choose to use these slots differently at other time.
    float amp;
    uint32_t freq;
    int64_t phase;

    uint64_t cmd_time;
    ChnCmd *cmd;
};

namespace scalar {

template<typename Params>
struct ChnFunc {
    struct Cache {
        uint32_t old_freq;
        uint32_t old_amp;
    };
    static NACS_INLINE float kernel_cacf(ChnState &state, uint8_t subidx)
    {
        NaCs::assume(subidx < 16);
        auto phase = state.phase;
        auto freq = state.freq;
        if (subidx == 0) {
            phase = phase + 2 * (int64_t)freq;
            state.phase = phase > Params::max_phase ? phase - Params::max_phase * 2 : phase;
        }
        // `phase` is the phase value at the end of the 16 sample block
        // phase number for the current sample is `phase - freq * (2 - subidx / 8)`
        float phasef = (float)phase - (float)freq * (2.0f - subidx * 0.125f);
        return state.amp * Kernel::sinpif_pi(phasef * Params::phase_scale);
    }
    template<typename RampF>
    static NACS_INLINE float kernel_carf(ChnState &state, RampF &&rampf,
                                         Cache &cache, uint8_t subidx)
    {
        NaCs::assume(subidx < 16);
        auto phase = state.phase;
        auto new_freq = state.freq;
        if (subidx == 0) {
            cache.old_freq = new_freq;
            state.freq = new_freq = rampf();
            phase = phase + cache.old_freq + new_freq;
            state.phase = phase > Params::max_phase ? phase - Params::max_phase * 2 : phase;
        }
        auto old_freq = cache.old_freq;
        // `phase` is the phase value at the end of the 16 sample block.
        // The frequency at any a particular subidx is
        // `old_freq * (1 - subidx / 16) + new_freq * subidx / 16`.
        auto subidx2_256 = float(subidx * subidx) / 256;
        // phase number for the current sample is
        // `original_phase + f0 * (subidx / 8 - subidx^2/16^2) + f1 * subidx^2 / 16^2`
        // or
        // `phase + f0 * (subidx / 8 - subidx^2/16^2 - 1) + f1 * (subidx^2 / 16^2 - 1)`
        float phasef = (float)phase + (float)old_freq * (subidx * 0.125f - subidx2_256 - 1.0f);
        phasef += (float)new_freq * (subidx2_256 - 1.0f);
        return state.amp * Kernel::sinpif_pi(phasef * Params::phase_scale);
    }
    template<typename RampA>
    static NACS_INLINE float kernel_racf(ChnState &state, RampA &&rampa,
                                         Cache &cache, uint8_t subidx)
    {
        NaCs::assume(subidx < 16);
        auto phase = state.phase;
        auto freq = state.freq;
        auto new_amp = state.amp;
        if (subidx == 0) {
            cache.old_amp = new_amp;
            state.amp = new_amp = rampa();
            phase = phase + 2 * (int64_t)freq;
            state.phase = phase > Params::max_phase ? phase - Params::max_phase * 2 : phase;
        }
        float phasef = (float)phase - (float)freq * (2.0f - subidx * 0.125f);
        float amp = (cache.old_amp * (16 - subidx) + new_amp * subidx) / 16;
        return amp * Kernel::sinpif_pi(phasef * Params::phase_scale);
    }
    template<typename RampA, typename RampF>
    static NACS_INLINE float kernel_rarf(ChnState &state, RampA &&rampa, RampF &&rampf,
                                         Cache &cache, uint8_t subidx)
    {
        NaCs::assume(subidx < 16);
        auto phase = state.phase;
        auto new_freq = state.freq;
        auto new_amp = state.amp;
        if (subidx == 0) {
            cache.old_amp = new_amp;
            state.amp = new_amp = rampa();

            cache.old_freq = new_freq;
            state.freq = new_freq = rampf();
            phase = phase + cache.old_freq + new_freq;
            state.phase = phase > Params::max_phase ? phase - Params::max_phase * 2 : phase;
        }
        auto old_freq = cache.old_freq;
        auto subidx2_256 = float(subidx * subidx) / 256;
        float phasef = (float)phase + (float)old_freq * (subidx * 0.125f - subidx2_256 - 1.0f);
        phasef += (float)new_freq * (subidx2_256 - 1.0f);
        float amp = (cache.old_amp * (16 - subidx) + new_amp * subidx) / 16;
        return amp * Kernel::sinpif_pi(phasef * Params::phase_scale);
    }
    struct FuncCaCf {
        static NACS_INLINE float calc(ChnState &state, Cache&, uint8_t subidx, uint64_t)
        {
            return kernel_cacf(state, subidx);
        }
    };
    template<typename RampF>
    struct FuncCaRf {
        static NACS_INLINE float calc(ChnState &state, Cache &cache, uint8_t subidx, uint64_t idx)
        {
            return kernel_carf(state, [&] __attribute__((always_inline)) {
                return RampF::calc(idx - state.cmd_time);
            }, cache, subidx);
        }
    };
    template<typename RampA>
    struct FuncRaCf {
        static NACS_INLINE float calc(ChnState &state, Cache &cache, uint8_t subidx, uint64_t idx)
        {
            return kernel_racf(state, [&] __attribute__((always_inline)) {
                return RampA::calc(idx - state.cmd_time);
            }, cache, subidx);
        }
    };
    template<typename RampA, typename RampF>
    struct FuncRaRf {
        static NACS_INLINE float calc(ChnState &state, Cache &cache, uint8_t subidx, uint64_t idx)
        {
            return kernel_racf(state, [&] __attribute__((always_inline)) {
                return RampA::calc(idx - state.cmd_time);
            }, [&] __attribute__((always_inline)) {
                return RampF::calc(idx - state.cmd_time); // TODO
            }, cache, subidx);
        }
    };
};

} // namespace scalar

}

#endif
