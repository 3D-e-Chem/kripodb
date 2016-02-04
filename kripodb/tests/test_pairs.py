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

from collections import Mapping
import StringIO

import kripodb.hdf5
import tables

import os
import tempfile

from intbitset import intbitset
from kripodb.hdf5 import PairsTable
from mock import Mock
from nose.tools import eq_, assert_raises

import kripodb.pairs as pairs


def tmpname():
    tmpf = tempfile.NamedTemporaryFile()
    out_file = tmpf.name
    tmpf.close()
    return out_file


class MockedIntbitsetDict(Mapping):
    def __init__(self, thedict, number_of_bits):
        self.dict = thedict
        self.number_of_bits = number_of_bits

    def __len__(self):
        return len(self.dict)

    def __getitem__(self, key):
        return self.dict[key]

    def __iter__(self):
        return self.dict.__iter__()

    def materialize(self):
        return self.dict


class Testpairs(object):
    def setUp(self):
        self.number_of_bits = 8
        self.bitsets = MockedIntbitsetDict({
            'a': intbitset([1, 2, 3]),
            'b': intbitset([1, 2, 4, 5, 8]),
            'c': intbitset([1, 2, 4, 8])
        }, 8)
        self.pairs = [
            ('a', 'c', 0.7523),
            ('b', 'c', 0.8342)]
        self.id2label = {
            1: 'a',
            2: 'b',
            3: 'c'
        }
        self.label2id = {
            'a': 1,
            'b': 2,
            'c': 3
        }
        self.compact_pairs = [(1, 3, 75), (2, 3, 83)]
        self.precision = 100
        self.h5filename = tmpname()

    def tearDown(self):
        if os.path.isfile(self.h5filename):
            os.remove(self.h5filename)

    def test_dump_pairs_tsv(self):
        out = StringIO.StringIO()

        pairs.dump_pairs_tsv(self.pairs, out)
        result = out.getvalue()

        expected = "a\tc\t0.7523\nb\tc\t0.8342\n"
        eq_(result, expected)

    def fill_matrix(self):
        kripodb.pairs.dump_pairs_hdf5(self.pairs,
                                      self.label2id,
                                      self.precision,
                                      2,
                                      self.h5filename)

    def empty_matrix(self):
        kripodb.pairs.dump_pairs_hdf5([],
                                      {},
                                      self.precision,
                                      0,
                                      self.h5filename)

    def test_similar_run(self):
        self.fill_matrix()
        out = StringIO.StringIO()

        pairs.similar_run('a', self.h5filename, 0.55, out)

        result = out.getvalue()
        expected = "a\tc\t0.75\n"
        eq_(result, expected)

    def test_similar_run_nohits(self):
        self.fill_matrix()
        out = StringIO.StringIO()

        hits = pairs.similar_run('a', self.h5filename, 0.99, out)

        result = out.getvalue()
        expected = ""
        eq_(result, expected)

    def test_dump_pairs_ashdf5(self):
        expectedrows = 2
        kripodb.pairs.dump_pairs_hdf5(self.pairs,
                                      self.label2id,
                                      self.precision,
                                      expectedrows,
                                      self.h5filename)

        h5file = tables.open_file(self.h5filename)
        mypairs = []
        for row in h5file.root.pairs:
            mypairs.append((row[0], row[1], row[2]))
        eq_(mypairs, self.compact_pairs)
        h5file.close()

    def test_dump_pairs_badformat(self):
        with assert_raises(LookupError) as cm:
            pairs.dump_pairs(self.bitsets,
                             self.bitsets,
                             'bikes',
                             'StringIO',
                             None,
                             self.number_of_bits,
                             0.4,
                             0.05,
                             self.label2id,
                             self.precision,
                             False
                             )

        eq_(cm.exception.message, 'Invalid output format')

    def test_dump_pairs_badffn(self):
        with assert_raises(Exception) as cm:
            pairs.dump_pairs(self.bitsets,
                             self.bitsets,
                             'hdf5_compact',
                             '-',
                             None,
                             self.number_of_bits,
                             0.4,
                             0.05,
                             self.label2id,
                             self.precision,
                             False
                             )

        eq_(cm.exception.message, "hdf5 formats can't be outputted to stdout")

    def test_dump_pairs_astsv(self):
        out = StringIO.StringIO()

        pairs.dump_pairs(self.bitsets,
                         self.bitsets,
                         'tsv',
                         'StringIO',
                         out,
                         self.number_of_bits,
                         0.4,
                         0.05,
                         self.label2id,
                         self.precision,
                         False
                         )
        result = out.getvalue()

        expected = "a\tc\t0.135555555556\n"
        eq_(result, expected)

    def test_dump_pairs_astsv_nomem(self):
        out = StringIO.StringIO()

        pairs.dump_pairs(self.bitsets,
                         self.bitsets,
                         'tsv',
                         'StringIO',
                         out,
                         self.number_of_bits,
                         0.4,
                         0.05,
                         self.label2id,
                         self.precision,
                         True
                         )
        result = out.getvalue()

        expected = "a\tc\t0.135555555556\n"
        eq_(result, expected)

    def test_dump_pairs_ashdf5(self):
        pairs.dump_pairs(self.bitsets,
                         self.bitsets,
                         'hdf5',
                         self.h5filename,
                         None,
                         self.number_of_bits,
                         0.4,
                         0.05,
                         self.label2id,
                         self.precision,
                         True
                         )

        h5file = tables.open_file(self.h5filename)
        mypairs = []
        for row in h5file.root.pairs:
            mypairs.append((row[0], row[1], row[2]))
        expected = [(1, 3, 13)]
        eq_(mypairs, expected)
        h5file.close()

    def test_distance2query(self):
        out = StringIO.StringIO()

        pairs.distance2query(self.bitsets,
                             'a',
                             out,
                             0.4,
                             0.05,
                             True
                             )
        result = out.getvalue()

        expected = "a\tc\t0.135555555556\n"
        eq_(result, expected)

    def test_total_number_of_pairs(self):
        self.fill_matrix()

        result = pairs.total_number_of_pairs([self.h5filename])

        eq_(result, 2)

    def test_labels_consistency_check_ok(self):
        self.fill_matrix()

        pairs.labels_consistency_check([self.h5filename])
        assert True
