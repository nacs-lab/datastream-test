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

#ifndef HELPERS_CPU_KERNEL_H
#define HELPERS_CPU_KERNEL_H

#define CL_HPP_TARGET_OPENCL_VERSION 200
// not sure if we want this for the real program but it's good enough for the test.
#define CL_HPP_ENABLE_EXCEPTIONS

#include <CL/cl2.hpp>

#include <yaml-cpp/yaml.h>

#include <nacs-utils/utils.h>

#include <functional>
#include <vector>

namespace OCL {

NACS_EXPORT(ds_helper) const char *strerror(int err);
NACS_EXPORT(ds_helper) bool catch_error(std::function<void()> func);
NACS_EXPORT(ds_helper) std::vector<cl::Device> all_ocl2_devices(bool includecpu=true);
NACS_EXPORT(ds_helper) void get_device_ids(const cl::Device &dev, YAML::Node &out);
NACS_EXPORT(ds_helper) YAML::Node get_device_ids(const cl::Device &dev);

}

#endif
