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


def calc_mean_onbit_density(bitsets, number_of_bits):
    all_nr_onbits = [len(v) for v in bitsets.itervalues()]
    mean_onbit = sum(all_nr_onbits) / float(len(bitsets))
    density = mean_onbit / number_of_bits
    return float(density)


def corrections(mean_onbit_density):
    p0 = mean_onbit_density
    corr_st = (2 - p0) / 3
    corr_sto = (1 + p0) / 3
    return corr_st, corr_sto


def distance(bitset1, bitset2, number_of_bits, corr_st, corr_sto):
    a = len(bitset1)
    b = len(bitset2)
    c = len(bitset1 & bitset2)
    n = number_of_bits
    st = float(c) / (a + b - c)
    st0 = (n - a - b - + c) / float(n - c)
    smt = corr_st * st + corr_sto * st0
    return smt


def distances(bitsets1, bitsets2, number_of_bits, corr_st, corr_sto, cutoff, full_matrix=False):
    for (label1, bs1) in bitsets1.iteritems():
        for (label2, bs2) in bitsets2.iteritems():
            if label1 == label2:
                # always skip self
                continue
            if not full_matrix and label1 > label2:
                continue

            d = distance(bs1, bs2, number_of_bits, corr_st, corr_sto)

            if d >= cutoff:
                yield label1, label2, d
