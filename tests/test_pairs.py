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
from collections import Mapping
import os
import tempfile

import tables
from six import StringIO
from intbitset import intbitset
from nose.tools import eq_, assert_raises

import kripodb.hdf5
import kripodb.pairs as pairs
from kripodb.hdf5 import DistanceMatrix


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
        self.compact_pairs = [(1, 3, 49301), (2, 3, 54669)]
        self.h5filename = tmpname()

    def tearDown(self):
        if os.path.isfile(self.h5filename):
            os.remove(self.h5filename)

    def test_dump_pairs_tsv(self):
        out = StringIO()

        pairs.dump_pairs_tsv(self.pairs, out)
        result = out.getvalue()

        expected = "a\tc\t0.7523\nb\tc\t0.8342\n"
        eq_(result, expected)

    def fill_matrix(self):
        kripodb.pairs.dump_pairs_hdf5(self.pairs,
                                      self.label2id,
                                      2,
                                      self.h5filename)

    def empty_matrix(self):
        kripodb.pairs.dump_pairs_hdf5([],
                                      {},
                                      0,
                                      self.h5filename)

    def test_similar_run(self):
        self.fill_matrix()
        out = StringIO()

        pairs.similar_run('a', self.h5filename, 0.55, out)

        result = out.getvalue()
        expected = "a\tc\t0.7523\n"
        eq_(result, expected)

    def test_similar_run_nohits(self):
        self.fill_matrix()
        out = StringIO()

        pairs.similar_run('a', self.h5filename, 0.99, out)

        result = out.getvalue()
        expected = ""
        eq_(result, expected)

    def test_dump_pairs_ashdf5(self):
        expectedrows = 2
        kripodb.pairs.dump_pairs_hdf5(self.pairs,
                                      self.label2id,
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
                             False
                             )

        eq_(cm.exception.args, ('Invalid output format',))

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
                             False
                             )

        eq_(cm.exception.args, ("hdf5 formats can't be outputted to stdout",))

    def test_dump_pairs_astsv(self):
        out = StringIO()

        pairs.dump_pairs(self.bitsets,
                         self.bitsets,
                         'tsv',
                         'StringIO',
                         out,
                         self.number_of_bits,
                         0.4,
                         0.05,
                         self.label2id,
                         False
                         )
        result = out.getvalue()

        expected = "a\tc\t0.13556\n"
        eq_(result, expected)

    def test_dump_pairs_astsv_nomem(self):
        out = StringIO()

        pairs.dump_pairs(self.bitsets,
                         self.bitsets,
                         'tsv',
                         'StringIO',
                         out,
                         self.number_of_bits,
                         0.4,
                         0.05,
                         self.label2id,
                         True
                         )
        result = out.getvalue()

        expected = "a\tc\t0.13556\n"
        eq_(result, expected)

    def test_distance2query(self):
        out = StringIO()

        pairs.distance2query(self.bitsets,
                             'a',
                             out,
                             0.4,
                             0.05,
                             True
                             )
        result = out.getvalue()

        expected = "a\tc\t0.13556\n"
        eq_(result, expected)

    def test_total_number_of_pairs(self):
        self.fill_matrix()

        result = pairs.total_number_of_pairs([self.h5filename])

        eq_(result, 2)


def test_merge():
    infiles = [tmpname(), tmpname(), tmpname()]

    outfile = tmpname()
    try:
        # fill infiles
        inmatrix1 = DistanceMatrix(infiles[0], 'w', 1, 2**16-1, 2)
        inmatrix1.update([('a', 'b', 0.2)], {'a': 1, 'b': 2, 'c': 3})
        inmatrix1.close()

        # matrix with same labels -> copy pairs table by dump/append, ignores labels tables
        inmatrix2 = DistanceMatrix(infiles[1], 'w', 2, 2**16-1, 3)
        inmatrix2.update([('a', 'c', 0.6)], {'a': 1, 'b': 2, 'c': 3})
        inmatrix2.close()

        # matrix generated with different labels -> copy pairs table by iterate/update, adds missing labels
        inmatrix3 = DistanceMatrix(infiles[2], 'w', 2, 2**16-1, 3)
        inmatrix3.update([('b', 'e', 0.4), ('e', 'f', 0.8)], {'b': 1, 'e': 2, 'f': 3})
        inmatrix3.close()

        pairs.merge(infiles, outfile)

        # compare it
        outmatrix = DistanceMatrix(outfile)
        result = list(outmatrix)
        outmatrix.close()
        expected = [('a', 'b', 0.2), ('a', 'c', 0.6), ('b', 'e', 0.4), ('e', 'f', 0.8)]
        eq_(result, expected)
    finally:
        for infile in infiles:
            if os.path.isfile(infile):
                os.remove(infile)
        if os.path.isfile(outfile):
            os.remove(outfile)
