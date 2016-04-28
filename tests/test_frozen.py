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
import os
from nose.tools import eq_, assert_raises
import pandas as pd
import pandas.util.testing as pdt
from tests.test_pairs import tmpname
from kripodb.frozen import FrozenDistanceMatrix
from kripodb.hdf5 import DistanceMatrix


class TestFrozenDistanceMatrix(object):
    pair_matrix_fn = None
    pair_matrix = None
    matrix_fn = None
    matrix = None

    def setUp(self):
        self.pair_matrix_fn = tmpname()
        self.pair_matrix = DistanceMatrix(self.pair_matrix_fn, 'a', driver='H5FD_CORE', driver_core_backing_store=0)
        labels = {'a': 0, 'b': 1, 'c': 2, 'd': 3}
        distances = [
            ('a', 'b', 0.9),
            ('a', 'c', 0.5),
            ('b', 'c', 0.6),
            ('d', 'c', 0.7)
        ]
        self.pair_matrix.update(distances, labels)
        self.matrix_fn = tmpname()
        self.matrix = FrozenDistanceMatrix(self.matrix_fn, 'a', driver='H5FD_CORE', driver_core_backing_store=0)

    def tearDown(self):
        self.pair_matrix.close()
        if os.path.isfile(self.pair_matrix_fn):
            os.remove(self.pair_matrix_fn)
        self.matrix.close()
        if os.path.isfile(self.matrix_fn):
            os.remove(self.matrix_fn)

    def test_from_pairs_defaults(self):
        self.matrix.from_pairs(self.pair_matrix, 10)

        result = self.matrix.to_pandas()
        labels = ['a', 'b', 'c', 'd']
        expected = pd.DataFrame([
            [0.0, 0.9, 0.5, 0.0],
            [0.9, 0.0, 0.6, 0.0],
            [0.5, 0.6, 0.0, 0.7],
            [0.0, 0.0, 0.7, 0.0]
        ], index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_from_pairs_multiframe(self):
        self.matrix.from_pairs(self.pair_matrix, 1, None, False)

        result = self.matrix.to_pandas()
        labels = ['a', 'b', 'c', 'd']
        expected = pd.DataFrame([
            [0.0, 0.9, 0.5, 0.0],
            [0.9, 0.0, 0.6, 0.0],
            [0.5, 0.6, 0.0, 0.7],
            [0.0, 0.0, 0.7, 0.0]
        ], index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_from_pairs_limited(self):
        self.matrix.from_pairs(self.pair_matrix, 1, 1, False)

        result = self.matrix.to_pandas()
        labels = ['a', 'b', 'c', 'd']
        expected = pd.DataFrame([
            [0.0, 0.9, 0.0, 0.0],
            [0.9, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0]
        ], index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_from_pairs_singlesided(self):
        self.matrix.from_pairs(self.pair_matrix, 10, None, True)

        result = self.matrix.to_pandas()
        labels = ['a', 'b', 'c', 'd']
        expected = pd.DataFrame([
            [0.0, 0.9, 0.5, 0.0],
            [0.0, 0.0, 0.6, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.7, 0.0]
        ], index=labels, columns=labels)
        pdt.assert_almost_equal(result, expected)

    def test_find_defaults(self):
        self.matrix.from_pairs(self.pair_matrix, 10)

        hits = self.matrix.find('c', 0.55)
        expected = [('d', 0.7), ('b', 0.6)]
        eq_(hits, expected)

    def test_find_limit(self):
        self.matrix.from_pairs(self.pair_matrix, 10)

        hits = self.matrix.find('c', 0.55 , 1)
        expected = [('d', 0.7)]
        eq_(hits, expected)

    def test_find_cutoffhigh_nohits(self):
        self.matrix.from_pairs(self.pair_matrix, 10)

        hits = self.matrix.find('c', 0.9)
        expected = []
        eq_(hits, expected)

    def test_find_badkey_keyerror(self):
        self.matrix.from_pairs(self.pair_matrix, 10)

        with assert_raises(KeyError):
            self.matrix.find('f', 0.45)
