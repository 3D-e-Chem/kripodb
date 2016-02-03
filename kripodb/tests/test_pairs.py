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

import StringIO

import tables

import os
import tempfile

from intbitset import intbitset
from mock import Mock, call
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


class Testpairs(object):
    def setup(self):
        self.bitsets = {
            'a': intbitset([1, 2, 3]),
            'b': intbitset([1, 2, 4, 5, 8]),
            'c': intbitset([1, 2, 4, 8])
        }
        self.pairs = [
            ('a', 'c', 0.7523),
            ('b', 'c', 0.8342)]
        self.number_of_bits = 8
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

    def test_dump_pairs_tsv(self):
        out = StringIO.StringIO()

        pairs.dump_pairs_tsv(self.pairs, out)
        result = out.getvalue()

        expected = "a\tc\t0.7523\nb\tc\t0.8342\n"
        eq_(result, expected)

    def test_similar_nohits(self):
        pairsdb = Mock(Table)
        pairsdb.attrs = {'score_precision': self.precision}
        pairsdb.where = Mock(return_value=[])

        hits = pairs.similar(1, pairsdb, self.id2label, 0.55)

        eq_(hits, [])
        pairsdb.where.assert_has_calls([
            call('(a == 1) & (score >= 55)'),
            call('(b == 1) & (score >= 55)')
        ])

    def test_similar_hitsleft(self):
        pairsdb = Mock(Table)
        pairsdb.attrs = {'score_precision': self.precision}
        pairsdb.where.side_effect = mypairs

        hits = pairs.similar(1, pairsdb, self.id2label, 0.55)

        expected = [('a', 0.75, 'c')]
        eq_(hits, expected)

    def test_similar_hitsright(self):
        pairsdb = Mock(Table)
        pairsdb.attrs = {'score_precision': self.precision}
        pairsdb.where.side_effect = mypairs

        hits = pairs.similar(3, pairsdb, self.id2label, 0.55)

        expected = [('c', 0.85, 'b')]
        eq_(hits, expected)

    def test_dump_pairs_hdf5_compact(self):

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


    def test_dump_pairs_tsv(self):
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

    def test_dump_pairs_tsv_nomem(self):
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