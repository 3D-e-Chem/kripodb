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

from nose.tools import eq_
from numpy.testing import assert_array_almost_equal, assert_almost_equal

from kripodb.hdf5 import SimilarityMatrix
from .utils import SimilarityMatrixInMemory


class TestSimilarityMatrix(object):
    matrix = None

    def setUp(self):
        self.matrix = SimilarityMatrix('data/similarities.h5')

    def tearDown(self):
        self.matrix.close()

    def test_find_1(self):
        result = list(self.matrix.find('2n6i_4FU_frag1', 0.98))

        expected = [('2n6i_4FU_frag2', 1.0), ('2n6i_4FU_frag6', 1.0)]
        assert_array_almost_equal([r[1] for r in result], [r[1] for r in expected], 3)
        eq_([r[0] for r in result], [r[0] for r in result])

    def test_iter_first2(self):
        myiter = iter(self.matrix)

        result = [next(myiter), next(myiter)]

        expected = [('2mlm_2W7_frag1', '2mlm_2W7_frag2', 0.5877), ('2mlm_2W7_frag2', '3wvm_STE_frag1', 0.4633)]
        assert_array_almost_equal([r[2] for r in result], [r[2] for r in expected], 3)
        eq_([(r[0], r[1],) for r in result], [(r[0], r[1],) for r in result])

    def test_iter_last(self):
        myiter = iter(self.matrix)

        result = None
        for row in myiter:
            result = row

        expected = ('3wyl_3KB_frag20', '3wyl_3KB_frag21', 0.999496452277409)
        assert_almost_equal(result[2], expected[2], 5)
        eq_(result[:2], expected[:2])


class TestPairsTable(object):
    def test_count(self):
        with SimilarityMatrixInMemory() as matrix:
            labels = {'a': 0, 'b': 1, 'c': 2, 'd': 3}
            similarities = [
                ('a', 'b', 0.9),
                ('a', 'c', 0.6),
                ('b', 'c', 0.6),
                ('d', 'c', 0.7)
            ]
            matrix.update(similarities, labels)

            counts = list(matrix.count(100000))

            expected = [(0.6,  2),
                        (0.7,  1),
                        (0.9,  1)]
            assert_array_almost_equal(counts, expected, 6)

    def test_count_multiframe(self):
        with SimilarityMatrixInMemory() as matrix:
            labels = {'a': 0, 'b': 1, 'c': 2, 'd': 3}
            similarities = [
                ('a', 'b', 0.9),
                ('a', 'c', 0.6),
                ('b', 'c', 0.6),
                ('d', 'c', 0.7)
            ]
            matrix.update(similarities, labels)

            counts = list(matrix.count(2))

            expected = [(0.6,  2),
                        (0.7,  1),
                        (0.9,  1)]
            assert_array_almost_equal(counts, expected, 6)
