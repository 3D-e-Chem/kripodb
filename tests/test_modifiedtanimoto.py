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

from nose.tools import assert_almost_equal
from intbitset import intbitset
from kripodb import modifiedtanimoto


class TestAlgorithm(object):
    number_of_bits = None
    corr_st = None
    corr_sto = None

    def setup(self):
        self.number_of_bits = 100
        self.corr_st = 0.663333333333
        self.corr_sto = 0.336666666667

    def test_calc_mean_onbit_density(self):
        bitsets = {
            'a': intbitset([1, 2, 3]),
            'b': intbitset([1, 2, 4, 5, 8]),
            'c': intbitset([1, 2, 4, 8])
        }

        result = modifiedtanimoto.calc_mean_onbit_density(bitsets, self.number_of_bits)

        expected = 0.04
        assert result == expected

    def test_corrections(self):
        corr_st, corr_sto = modifiedtanimoto.corrections(0.01)

        assert_almost_equal(corr_st, 0.663333333333)
        assert_almost_equal(corr_sto, 0.336666666667)

    def test_distance(self):
        bitset1 = intbitset([1, 2, 3])
        bitset2 = intbitset([1, 2, 4, 8])

        result = modifiedtanimoto.distance(bitset1, bitset2,
                                    self.number_of_bits,
                                    self.corr_st, self.corr_sto)

        expected = 0.5779523809525572
        assert_almost_equal(result, expected)

    def test_distances(self):
        bitsets = {
            'a': intbitset([1, 2, 3]),
            'b': intbitset([1, 2, 4, 5, 8]),
            'c': intbitset([1, 2, 4, 8])
        }

        iterator = modifiedtanimoto.distances(bitsets, bitsets,
                                       self.number_of_bits,
                                       self.corr_st, self.corr_sto,
                                       0.55)
        result = [r for r in iterator]

        expected = [
            ('a', 'c', 0.5779523809525572),
            ('b', 'c', 0.8357708333333689)]
        # pair a-c is below cutoff with distance of 0.53

        assert len(result) == len(expected)
        for i, r in enumerate(result):
            assert r[0] == expected[i][0]
            assert r[1] == expected[i][1]
            assert_almost_equal(r[2], expected[i][2])

    def test_distances_fullmatrix(self):
        bitsets = {
            'a': intbitset([1, 2, 3]),
            'b': intbitset([1, 2, 4, 5, 8]),
            'c': intbitset([1, 2, 4, 8])
        }

        iterator = modifiedtanimoto.distances(bitsets, bitsets,
                                       self.number_of_bits,
                                       self.corr_st, self.corr_sto,
                                       0.55, True)
        result = [r for r in iterator]

        expected = [
            ('a', 'c', 0.5779523809525572),
            ('c', 'a', 0.5779523809525572),
            ('c', 'b', 0.8357708333333689),
            ('b', 'c', 0.8357708333333689)]
        # pair a-c is below cutoff with distance of 0.53

        print result

        assert len(result) == len(expected)
        for i, r in enumerate(result):
            assert r[0] == expected[i][0]
            assert r[1] == expected[i][1]
            assert_almost_equal(r[2], expected[i][2])
