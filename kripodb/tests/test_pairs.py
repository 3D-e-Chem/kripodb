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

import tables

import os
import tempfile

from intbitset import intbitset
from kripodb.db import FragmentsDb
from mock import Mock, call, patch
from nose.tools import eq_, assert_raises
from tables.table import Table

import kripodb.pairs as pairs


def mypairs(query):
    return_values = {
        '(a == 1) & (score >= 55)': [(1, 3, 75)],
        '(b == 3) & (score >= 55)': [(2, 3, 85)]
    }
    return return_values.get(query, [])


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


class MockedLookup(object):
    def __init__(self, thedict):
        self.lookup = {}
        for k, v in thedict.iteritems():
            self.lookup['frag_id == {}'.format(k)] = iter([(k, v)])
            self.lookup['label == "{}"'.format(v)] = iter([(k, v)])
        print self.lookup

    def where(self, query):
        return self.lookup[query]


class Testpairs(object):
    def setup(self):
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
        self.lookup = MockedLookup(self.id2label)

    def test_dump_pairs_tsv(self):
        out = StringIO.StringIO()

        pairs.dump_pairs_tsv(self.pairs, out)
        result = out.getvalue()

        expected = "a\tc\t0.7523\nb\tc\t0.8342\n"
        eq_(result, expected)

    def mock_pairsdb(self):
        pairsdb = Mock(Table)
        pairsdb.attrs = {'score_precision': self.precision}
        pairsdb.where = Mock(return_value=[])
        return pairsdb

    def test_similar_run(self):
        h5file = tmpname()
        try:
            pairs.dump_pairs_hdf5_compact(self.pairs,
                                          self.label2id,
                                          self.precision,
                                          2,
                                          h5file)

            fragments = Mock(FragmentsDb)
            fragments.label2id.return_value = self.label2id
            fragments.id2label.return_value = self.id2label
            out = StringIO.StringIO()

            pairs.similar_run('a', h5file, 0.55, out)

            result = out.getvalue()
            expected = "a\tc\t0.75\n"
            eq_(result, expected)
        finally:
            os.remove(h5file)

    def test_similar_nohits(self):
        pairsdb = self.mock_pairsdb()

        hits = pairs.similar(1, pairsdb, self.lookup, 0.55)

        eq_(hits, [])
        pairsdb.where.assert_has_calls([
            call('(a == 1) & (score >= 55)'),
            call('(b == 1) & (score >= 55)')
        ])

    def test_similar_hitsleft(self):
        pairsdb = Mock(Table)
        pairsdb.attrs = {'score_precision': self.precision}
        pairsdb.where.side_effect = mypairs

        hits = pairs.similar(1, pairsdb, self.lookup, 0.55)

        expected = [('a', 'c', 0.75)]
        eq_(hits, expected)

    def test_similar_hitsright(self):
        pairsdb = Mock(Table)
        pairsdb.attrs = {'score_precision': self.precision}
        pairsdb.where.side_effect = mypairs

        hits = pairs.similar(3, pairsdb, self.lookup, 0.55)

        expected = [('c', 'b', 0.85)]
        eq_(hits, expected)

    def test_dump_pairs_ashdf5_compact(self):
        expectedrows = 2
        out_file = tmpname()

        try:
            pairs.dump_pairs_hdf5_compact(self.pairs,
                                          self.label2id,
                                          self.precision,
                                          expectedrows,
                                          out_file)

            h5file = tables.open_file(out_file)
            mypairs = []
            for row in h5file.root.pairs:
                mypairs.append((row[0], row[1], row[2]))
            eq_(mypairs, self.compact_pairs)
            h5file.close()
        finally:
            os.remove(out_file)

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

    def test_dump_pairs_ashdf5_compact(self):
        out_file = tmpname()

        try:
            pairs.dump_pairs(self.bitsets,
                             self.bitsets,
                             'hdf5_compact',
                             out_file,
                             None,
                             self.number_of_bits,
                             0.4,
                             0.05,
                             self.label2id,
                             self.precision,
                             True
                             )

            h5file = tables.open_file(out_file)
            mypairs = []
            for row in h5file.root.pairs:
                mypairs.append((row[0], row[1], row[2]))
            expected = [(1, 3, 13)]
            eq_(mypairs, expected)
            h5file.close()
        finally:
            os.remove(out_file)

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