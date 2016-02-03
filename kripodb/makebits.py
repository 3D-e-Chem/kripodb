# Copyright 2016 Netherlands eScience Center
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

from intbitset import intbitset


def read_header(line):
    (format_name, format_version, fp_size, label) = line.strip().split(' ')
    fp_size = int(fp_size)
    return format_name, format_version, fp_size, label


def read_fp_size(header):
    return header[2]


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
    return fid, bitset


def read_file(infile):
    header = infile.readline()
    (format_name, format_version, fp_size, label) = read_header(header)
    bitsets = {}
    for line in infile:
        (fid, bitset) = read_bitset(line, fp_size)
        bitsets[fid] = bitset
    return bitsets, fp_size


def iter_file(infile):
    """
    Yields header first then tuples of idenfitier and intbitset object

    :param infile:
    :return:
    """
    header = read_header(infile.readline())
    fp_size = read_fp_size(header)
    yield header
    for line in infile:
        (fid, bitset) = read_bitset(line, fp_size)
        yield fid, bitset


def write_header(fp_size):
    return "MAKEBITS 1.0 {} BigGrid\n".format(fp_size)


def write_bitset(fid, bitset):
    bits = bitset.extract_finite_list()
    bits.extend([0, len(bitset)])
    return fid + " " + " ".join([str(d) for d in bits]) + "\n"


def write_file(fp_size, bitsets, fn):
    fn.write(write_header(fp_size))
    for fid, bitset in bitsets.iteritems():
        fn.write(write_bitset(fid, bitset))
