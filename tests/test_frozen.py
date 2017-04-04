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
import numpy as np
from numpy.testing import assert_array_almost_equal
import pandas as pd
import pandas.util.testing as pdt

from .utils import FrozenSimilarityMatrixInMemory, SimilarityMatrixInMemory


@pytest.fixture
def similarity_matrix():
    pair_matrix_inmem = SimilarityMatrixInMemory()
    pair_matrix = pair_matrix_inmem.matrix
    labels = {'a': 0, 'b': 1, 'c': 2, 'd': 3}
    similarities = [
        ('a', 'b', 0.9),
        ('a', 'c', 0.5),
        ('b', 'c', 0.6),
        ('d', 'c', 0.7)
    ]
    pair_matrix.update(similarities, labels)
    yield pair_matrix
    pair_matrix_inmem.close()


@pytest.fixture
def frozen_similarity_matrix():
    matrix_inmem = FrozenSimilarityMatrixInMemory()
    matrix = matrix_inmem.matrix
    yield matrix
    matrix_inmem.close()


class TestFrozenSimilarityMatrix(object):
    def test_from_pairs_defaults(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        result = frozen_similarity_matrix.to_pandas()
        labels = ['a', 'b', 'c', 'd']
        expected = pd.DataFrame([
            [0.0, 0.9, 0.5, 0.0],
            [0.9, 0.0, 0.6, 0.0],
            [0.5, 0.6, 0.0, 0.7],
            [0.0, 0.0, 0.7, 0.0]
        ], index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_from_pairs_multiframe(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 1, None, False)

        result = frozen_similarity_matrix.to_pandas()
        labels = ['a', 'b', 'c', 'd']
        expected = pd.DataFrame([
            [0.0, 0.9, 0.5, 0.0],
            [0.9, 0.0, 0.6, 0.0],
            [0.5, 0.6, 0.0, 0.7],
            [0.0, 0.0, 0.7, 0.0]
        ], index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_from_pairs_limited(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 1, 1, False)

        result = frozen_similarity_matrix.to_pandas()
        labels = ['a', 'b', 'c', 'd']
        expected = pd.DataFrame([
            [0.0, 0.9, 0.0, 0.0],
            [0.9, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0]
        ], index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_from_pairs_singlesided(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10, None, True)

        result = frozen_similarity_matrix.to_pandas()
        print(result)
        labels = ['a', 'b', 'c', 'd']
        expected = pd.DataFrame([
            [0.0, 0.9, 0.5, 0.0],
            [0.0, 0.0, 0.6, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.7, 0.0]
        ], index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_find_defaults(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        hits = frozen_similarity_matrix.find('c', 0.55)
        expected = [('d', 0.7), ('b', 0.6)]
        assert hits == expected

    def test_find_limit(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        hits = frozen_similarity_matrix.find('c', 0.55 , 1)
        expected = [('d', 0.7)]
        assert hits == expected

    def test_find_cutoffhigh_nohits(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        hits = frozen_similarity_matrix.find('c', 0.9)
        expected = []
        assert hits == expected

    def test_find_badkey_keyerror(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        with pytest.raises(KeyError):
            frozen_similarity_matrix.find('f', 0.45)

    def test_find_singlesided(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10, None, True)
        print(frozen_similarity_matrix.scores.read())
        hits = frozen_similarity_matrix.find('c', 0.0)
        expected = []
        assert hits == expected

    def test_from_array(self, similarity_matrix, frozen_similarity_matrix):
        labels = ['a', 'b', 'c', 'd']
        data = [
            [0.0, 0.9, 0.5, 0.0],
            [0.9, 0.0, 0.6, 0.0],
            [0.5, 0.6, 0.0, 0.7],
            [0.0, 0.0, 0.7, 0.0],
        ]
        frozen_similarity_matrix.from_array(np.array(data), labels)

        result = frozen_similarity_matrix.to_pandas()
        expected = pd.DataFrame(data, index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_getitem(self, similarity_matrix, frozen_similarity_matrix):
        labels = ['a', 'b', 'c', 'd']
        data = [
            [0.0, 0.9, 0.5, 0.0],
            [0.9, 0.0, 0.6, 0.0],
            [0.5, 0.6, 0.0, 0.7],
            [0.0, 0.0, 0.7, 0.0],
        ]
        frozen_similarity_matrix.from_array(np.array(data), labels)

        result = []
        for label in labels:
            result.append(frozen_similarity_matrix[label])
        expected = [
            [(u'b', 0.9), (u'c', 0.5), (u'd', 0.0)],
            [(u'a', 0.9), (u'c', 0.6), (u'd', 0.0)],
            [(u'a', 0.5), (u'b', 0.6), (u'd', 0.7)],
            [(u'a', 0.0), (u'b', 0.0), (u'c', 0.7)],
        ]
        assert result == expected

    def test_to_pairs(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)
        with SimilarityMatrixInMemory() as thawed_matrix:

            frozen_similarity_matrix.to_pairs(thawed_matrix)

            # compare scores
            assert_array_almost_equal(
                [d[2] for d in thawed_matrix],
                [d[2] for d in similarity_matrix],
                5
            )

    def test_count(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        counts = list(frozen_similarity_matrix.count())
        expected = [(0.5, 1),
                    (0.6, 1),
                    (0.7, 1),
                    (0.9, 1)]
        assert_array_almost_equal(counts, expected, 6)

    def test_count_lower_triangle(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        counts = list(frozen_similarity_matrix.count(lower_triangle=True))
        expected = [(0.5, 1),
                    (0.6, 1),
                    (0.7, 1),
                    (0.9, 1)]
        assert_array_almost_equal(counts, expected, 6)

    def test_count_cmp_triangle(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        upper_triangle = list(frozen_similarity_matrix.count(lower_triangle=False))
        lower_triangle = list(frozen_similarity_matrix.count(lower_triangle=True))
        assert_array_almost_equal(upper_triangle, lower_triangle, 6)

    def test_count_raw_score(self, similarity_matrix, frozen_similarity_matrix):
        frozen_similarity_matrix.from_pairs(similarity_matrix, 10)

        counts = list(frozen_similarity_matrix.count(raw_score=True))
        expected = [(32767, 1),
                    (39321, 1),
                    (45874, 1),
                    (58981, 1)]
        assert_array_almost_equal(counts, expected, 6)
