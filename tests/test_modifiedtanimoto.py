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

from __future__ import absolute_import

import pytest
from pyroaring import BitMap

from kripodb import modifiedtanimoto


def assert_similarities(result, expected):
    result = sorted(result)
    expected = sorted(expected)
    assert len(result) == len(expected)
    for i, r in enumerate(result):
        assert r[0] == expected[i][0]
        assert r[1] == expected[i][1]
        pytest.approx(r[2], expected[i][2])


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
            'a': BitMap([1, 2, 3]),
            'b': BitMap([1, 2, 4, 5, 8]),
            'c': BitMap([1, 2, 4, 8])
        }

        result = modifiedtanimoto.calc_mean_onbit_density(bitsets.values(), self.number_of_bits)

        expected = 0.04
        assert result == expected

    def test_corrections(self):
        corr_st, corr_sto = modifiedtanimoto.corrections(0.01)

        pytest.approx(corr_st, 0.663333333333)
        pytest.approx(corr_sto, 0.336666666667)

    def test_similarity(self):
        bitset1 = BitMap([1, 2, 3])
        bitset2 = BitMap([1, 2, 4, 8])

        result = modifiedtanimoto.similarity(bitset1, bitset2,
                                             self.number_of_bits,
                                             self.corr_st, self.corr_sto)

        expected = 0.5779523809525572
        pytest.approx(result, expected)

    def test_similarities_ignore_upper_triangle(self):
        bitsets = {
            'a': BitMap([1, 2, 3]),
            'b': BitMap([1, 2, 4, 5, 8]),
            'c': BitMap([1, 2, 4, 8])
        }

        iterator = modifiedtanimoto.similarities(bitsets, bitsets,
                                                 self.number_of_bits,
                                                 self.corr_st, self.corr_sto,
                                                 0.55, True)
        result = [r for r in iterator]

        expected = [
            ('a', 'c', 0.5779523809525572),
            ('b', 'c', 0.8357708333333689)]
        # pair a-c is below cutoff with similarity of 0.53
        assert_similarities(result, expected)

    def test_similarities(self):
        bitsets = {
            'a': BitMap([1, 2, 3]),
            'b': BitMap([1, 2, 4, 5, 8]),
            'c': BitMap([1, 2, 4, 8])
        }

        iterator = modifiedtanimoto.similarities(bitsets, bitsets,
                                                 self.number_of_bits,
                                                 self.corr_st, self.corr_sto,
                                                 0.55, False)
        result = [r for r in iterator]

        expected = [
            ('a', 'c', 0.5779523809525572),
            ('c', 'a', 0.5779523809525572),
            ('c', 'b', 0.8357708333333689),
            ('b', 'c', 0.8357708333333689)]
        # pair a-c is below cutoff with similarity of 0.53
        assert_similarities(result, expected)


@pytest.mark.parametrize("bitset2,expected_score", (
    (BitMap((1, 2, 3, 4)), 1.0),
    (BitMap((5, 6, 7, 8)), 0.33),
    (BitMap((3, 4, 5, 6)), 0.55),
    (BitMap((1, 2, 5, 6)), 0.55),
    (BitMap((2, 3, 4, 5)), 0.73),
    (BitMap((1, 5, 6, 7)), 0.424),
    (BitMap((1, 2, 3, 4, 5)), 0.86),
    (BitMap((1, 2, 3)), 0.83),
))
def test_similarity_numberofbits400(bitset2, expected_score):
    number_of_bits = 400
    corr_st, corr_sto = modifiedtanimoto.corrections(0.01)
    bitset1 = BitMap([1, 2, 3, 4])

    result = modifiedtanimoto.similarity(bitset1, bitset2,
                                         number_of_bits,
                                         corr_st, corr_sto)

    assert result == pytest.approx(expected_score, rel=1e-2)
