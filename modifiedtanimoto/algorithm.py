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


def mean_onbit_density(bitsets, fp_size):
    all_nr_onbits = [len(d) for d in bitsets]
    mean_onbit = sum(all_nr_onbits) / float(len(bitsets))
    mean_onbit_density = mean_onbit / fp_size
    return float(mean_onbit_density)


def corrections(mean_onbit_density):
    p0 = mean_onbit_density
    corr_st = (2 - p0) / 3
    corr_sto = (1 + p0) / 3
    return (corr_st, corr_sto)


def distance(bitset1, bitset2, fp_size, corr_st, corr_sto):
    a = len(bitset1)
    b = len(bitset2)
    c = len(bitset1 & bitset2)
    n = fp_size
    st = float(c) / (a + b - c)
    st0 = (n - a - b - + c) / float(n - c)
    smt = corr_st * st + corr_sto * st0
    return smt
