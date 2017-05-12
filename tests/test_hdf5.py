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
from numpy.testing import assert_array_almost_equal, assert_almost_equal

from kripodb.hdf5 import SimilarityMatrix
from .utils import SimilarityMatrixInMemory


@pytest.fixture
def matrix():
    sim_matrix = SimilarityMatrix('data/similarities.h5')
    yield sim_matrix
    sim_matrix.close();


@pytest.fixture
def empty_matrix():
    with SimilarityMatrixInMemory() as sim_matrix:
        yield sim_matrix


@pytest.fixture
def example_matrix():
    with SimilarityMatrixInMemory() as sim_matrix:
        labels = {'a': 0, 'b': 1, 'c': 2, 'd': 3}
        similarities = [
            ('a', 'b', 0.9),
            ('a', 'c', 0.6),
            ('b', 'c', 0.6),
            ('d', 'c', 0.7)
        ]
        sim_matrix.update(similarities, labels)
        yield sim_matrix


class TestSimilarityMatrix(object):
    def test_find_1(self, matrix):
        result = list(matrix.find('2n6i_4FU_frag1', 0.98))

        expected = [('2n6i_4FU_frag2', 1.0), ('2n6i_4FU_frag6', 1.0)]
        assert_array_almost_equal([r[1] for r in result], [r[1] for r in expected], 3)
        assert [r[0] for r in result] == [r[0] for r in result]

    def test_iter_first2(self, matrix):
        myiter = iter(matrix)

        result = [next(myiter), next(myiter)]

        expected = [('2mlm_2W7_frag1', '2mlm_2W7_frag2', 0.5877), ('2mlm_2W7_frag2', '3wvm_STE_frag1', 0.4633)]
        assert_array_almost_equal([r[2] for r in result], [r[2] for r in expected], 3)
        assert [(r[0], r[1],) for r in result] == [(r[0], r[1],) for r in result]

    def test_iter_last(self, matrix):
        myiter = iter(matrix)

        result = None
        for row in myiter:
            result = row

        expected = ('3wyl_3KB_frag20', '3wyl_3KB_frag21', 0.999496452277409)
        assert_almost_equal(result[2], expected[2], 5)
        assert result[:2] == expected[:2]

    def test_keep(self, example_matrix, empty_matrix):
        in_matrix = example_matrix
        out_matrix = empty_matrix
        frags2keep = {'a', 'b'}
        in_matrix.keep(out_matrix, frags2keep)

        expected_labels = {'a', 'b', 'c'}
        assert set(out_matrix.labels.label2ids().keys()) == expected_labels
        expected_similarities = {
            ('a', 'b', 0.9),
            ('a', 'c', 0.6),
            ('b', 'c', 0.6)
        }
        assert set(out_matrix) == expected_similarities

    def test_skip(self, example_matrix, empty_matrix):
        in_matrix = example_matrix
        out_matrix = empty_matrix
        frags2skip = {'b'}
        in_matrix.skip(out_matrix, frags2skip)

        expected_labels = {'a', 'c', 'd'}
        assert set(out_matrix.labels.label2ids().keys()) == expected_labels
        expected_similarities = {
            ('a', 'c', 0.6),
            ('d', 'c', 0.7),
        }
        assert set(out_matrix) == expected_similarities


class TestPairsTable(object):
    def test_count(self, example_matrix):
        counts = list(example_matrix.count(100000))

        expected = [(0.6,  2),
                    (0.7,  1),
                    (0.9,  1)]
        assert_array_almost_equal(counts, expected, 6)

    def test_count_rawscore(self, example_matrix):
            counts = list(example_matrix.count(100000, True))
            expected = [(39321,  2),
                        (45874,  1),
                        (58981,  1)]
            assert_array_almost_equal(counts, expected, 6)

    def test_count_multiframe(self, example_matrix):
            counts = list(example_matrix.count(2))

            expected = [(0.6,  2),
                        (0.7,  1),
                        (0.9,  1)]
            assert_array_almost_equal(counts, expected, 6)
