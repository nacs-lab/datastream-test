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

#include <fstream>
#include <iterator>
#include <vector>

#include <stdint.h>
#include <string.h>

// Given the number range, this is slightly more efficient than leb128
void encode(std::vector<uint8_t> &out, uint64_t v)
{
    auto b2 = uint16_t(v & 0x7fff);
    v = v >> 15;
    if (v)
        b2 = b2 | 0x8000;
    out.push_back(uint8_t(b2 & 0xff));
    out.push_back(uint8_t(b2 >> 8));
    while (v) {
        auto b = uint8_t(v & 0x7f);
        v = v >> 7;
        if (v)
            b = b | 0x80;
        out.push_back(b);
    }
}

// void encode_leb128(std::vector<uint8_t> &out, uint64_t v)
// {
//     while (true) {
//         auto b = uint8_t(v & 0x7f);
//         v = v >> 7;
//         if (!v) {
//             out.push_back(b);
//             return;
//         }
//         out.push_back(b | 0x80);
//     }
// }

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "Needs two arguments\n");
        exit(1);
    }
    std::ifstream istm(argv[1]);
    std::vector<uint64_t> data;
    uint64_t buff;
    while (istm.read((char*)&buff, sizeof(uint64_t)))
        data.push_back(buff);
    istm.close();
    auto v0 = data[0];
    auto ndata = data.size();
    uint64_t vmin = data[1] - v0;
    for (size_t i = 0; i < ndata - 1; i++) {
        auto v1 = data[i + 1];
        auto vd = v1 - v0;
        if (vd < vmin)
            vmin = vd;
        data[i + 1] = vd;
        v0 = v1;
    }
    for (size_t i = 1; i < ndata; i++)
        data[i] -= vmin;
    data.insert(data.begin() + 1, vmin);
    std::vector<uint8_t> out;
    for (auto d: data)
        encode(out, d);
    std::ofstream ostm(argv[2]);
    ostm.write((const char*)&out[0], out.size());
    return 0;
}
