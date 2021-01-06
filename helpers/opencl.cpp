/*************************************************************************
 *   Copyright (c) 2020 - 2020 Yichao Yu <yyc1992@gmail.com>             *
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

#include "opencl.h"

#include <iostream>
#include <regex>

namespace OCL {

NACS_EXPORT() const char *strerror(int err)
{
    switch (err) {
#define handle_error(x) case x: return #x
        handle_error(CL_SUCCESS);
        handle_error(CL_DEVICE_NOT_FOUND);
        handle_error(CL_DEVICE_NOT_AVAILABLE);
        handle_error(CL_COMPILER_NOT_AVAILABLE);
        handle_error(CL_MEM_OBJECT_ALLOCATION_FAILURE);
        handle_error(CL_OUT_OF_RESOURCES);
        handle_error(CL_OUT_OF_HOST_MEMORY);
        handle_error(CL_PROFILING_INFO_NOT_AVAILABLE);
        handle_error(CL_MEM_COPY_OVERLAP);
        handle_error(CL_IMAGE_FORMAT_MISMATCH);
        handle_error(CL_IMAGE_FORMAT_NOT_SUPPORTED);
        handle_error(CL_BUILD_PROGRAM_FAILURE);
        handle_error(CL_MAP_FAILURE);
#ifdef CL_VERSION_1_1
        handle_error(CL_MISALIGNED_SUB_BUFFER_OFFSET);
        handle_error(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST);
#endif
#ifdef CL_VERSION_1_2
        handle_error(CL_COMPILE_PROGRAM_FAILURE);
        handle_error(CL_LINKER_NOT_AVAILABLE);
        handle_error(CL_LINK_PROGRAM_FAILURE);
        handle_error(CL_DEVICE_PARTITION_FAILED);
        handle_error(CL_KERNEL_ARG_INFO_NOT_AVAILABLE);
#endif
        handle_error(CL_INVALID_VALUE);
        handle_error(CL_INVALID_DEVICE_TYPE);
        handle_error(CL_INVALID_PLATFORM);
        handle_error(CL_INVALID_DEVICE);
        handle_error(CL_INVALID_CONTEXT);
        handle_error(CL_INVALID_QUEUE_PROPERTIES);
        handle_error(CL_INVALID_COMMAND_QUEUE);
        handle_error(CL_INVALID_HOST_PTR);
        handle_error(CL_INVALID_MEM_OBJECT);
        handle_error(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR);
        handle_error(CL_INVALID_IMAGE_SIZE);
        handle_error(CL_INVALID_SAMPLER);
        handle_error(CL_INVALID_BINARY);
        handle_error(CL_INVALID_BUILD_OPTIONS);
        handle_error(CL_INVALID_PROGRAM);
        handle_error(CL_INVALID_PROGRAM_EXECUTABLE);
        handle_error(CL_INVALID_KERNEL_NAME);
        handle_error(CL_INVALID_KERNEL_DEFINITION);
        handle_error(CL_INVALID_KERNEL);
        handle_error(CL_INVALID_ARG_INDEX);
        handle_error(CL_INVALID_ARG_VALUE);
        handle_error(CL_INVALID_ARG_SIZE);
        handle_error(CL_INVALID_KERNEL_ARGS);
        handle_error(CL_INVALID_WORK_DIMENSION);
        handle_error(CL_INVALID_WORK_GROUP_SIZE);
        handle_error(CL_INVALID_WORK_ITEM_SIZE);
        handle_error(CL_INVALID_GLOBAL_OFFSET);
        handle_error(CL_INVALID_EVENT_WAIT_LIST);
        handle_error(CL_INVALID_EVENT);
        handle_error(CL_INVALID_OPERATION);
        handle_error(CL_INVALID_GL_OBJECT);
        handle_error(CL_INVALID_BUFFER_SIZE);
        handle_error(CL_INVALID_MIP_LEVEL);
        handle_error(CL_INVALID_GLOBAL_WORK_SIZE);
#ifdef CL_VERSION_1_1
        handle_error(CL_INVALID_PROPERTY);
#endif
#ifdef CL_VERSION_1_2
        handle_error(CL_INVALID_IMAGE_DESCRIPTOR);
        handle_error(CL_INVALID_COMPILER_OPTIONS);
        handle_error(CL_INVALID_LINKER_OPTIONS);
        handle_error(CL_INVALID_DEVICE_PARTITION_COUNT);
#endif
#ifdef CL_VERSION_2_0
        handle_error(CL_INVALID_PIPE_SIZE);
        handle_error(CL_INVALID_DEVICE_QUEUE);
#endif
#ifdef CL_VERSION_2_2
        handle_error(CL_INVALID_SPEC_ID);
        handle_error(CL_MAX_SIZE_RESTRICTION_EXCEEDED);
#endif
#undef handle_error
    default: return "UNKNOWN";
    }
}

NACS_EXPORT() bool catch_error(std::function<void()> func)
{
    try {
        func();
        return true;
    }
    catch (cl::BuildError &err) {
        std::cerr << err.what() << " " << strerror(err.err()) << std::endl;
        for (const auto &log: err.getBuildLog()) {
            std::cerr << "  " << log.first.getInfo<CL_DEVICE_NAME>() << ": "
                      << log.second << std::endl;
        }
    }
    catch (cl::Error &err) {
        std::cerr << err.what() << " " << strerror(err.err()) << std::endl;
    }
    return false;
}

const std::regex oclver_regex("OpenCL ([0-9]+)\\.");

NACS_EXPORT() std::vector<cl::Device> all_ocl2_devices(bool includecpu)
{
    auto type = CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_ACCELERATOR;
    if (includecpu)
        type |= CL_DEVICE_TYPE_CPU;
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    std::vector<cl::Device> devices;
    std::smatch match;
    for (auto &p: platforms) {
        std::vector<cl::Device> devs;
        p.getDevices(type, &devs);
        for (auto &dev: devs) {
            auto ver = dev.getInfo<CL_DEVICE_VERSION>();
            if (std::regex_search(ver, match, oclver_regex)) {
                if (std::stoi(match[1].str()) >= 2) {
                    devices.push_back(dev);
                }
            }
        }
    }
    return devices;
}

}
