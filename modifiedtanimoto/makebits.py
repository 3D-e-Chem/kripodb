# Copyright 2013 Netherlands eScience Center
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import intbitset


def read_header(line):
    (format_name, format_version, fp_size, label) = line.split(' ')
    fp_size = int(fp_size)
    return (format_name, format_version, fp_size, label)


def read_bitset(line, fp_size):
    row = line.split(' ', fp_size + 3)
    fid = row.pop(0)
    nr_onbits = int(row.pop())
    # ignore 0, is seperator between metadata and data
    row.pop()
    bits = [int(d) for d in row]
    bitset = intbitset(bits)
    if len(bitset) != nr_onbits:
        raise Exception('On bit checksum incorrect for {}'.format(fid))
    return (fid, bitset)


def read_file(file):
    header = file.readline()
    (format_name, format_version, fp_size, label) = read_header(header)
    bitsets = {}
    for line in file:
        (fid, bitset) = read_bitset(line, fp_size)
        bitsets[fid] = bitset
    return (bitsets, fp_size)
