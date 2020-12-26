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

#include <ostream>
#include <type_traits>

namespace Print {

namespace {

template<typename T>
using remove_cv_ref = std::remove_cv_t<std::remove_reference_t<T>>;

template<typename T>
using is_str_map = std::integral_constant<bool,
                                          std::is_base_of<
                                              std::string,
                                              remove_cv_ref<typename T::key_type>>::value &&
                                          !std::is_void<typename T::mapped_type>::value>;

}

template<typename T>
static inline std::enable_if_t<std::is_arithmetic<std::remove_reference_t<T>>::value>
json(std::ostream &stm, const T &v)
{
    stm << v;
}

template<typename T>
static inline std::enable_if_t<is_str_map<std::remove_reference_t<T>>::value>
json(std::ostream &stm, const T &d)
{
    bool first = true;
    stm << "{";
    for (const auto &kv: d) {
        if (first) {
            first = false;
            stm << " \"";
        }
        else {
            stm << ", \"";
        }
        // TODO better quote
        stm << kv.first << "\": ";
        json(stm, kv.second);
    }
    if (first) {
        stm << "}";
    }
    else {
        stm << " }";
    }
}

}
