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

import tables
from six import StringIO
from intbitset import intbitset
import pytest

import kripodb.hdf5
import kripodb.pairs as pairs
from kripodb.hdf5 import SimilarityMatrix
from .utils import tmpname


class MockedIntbitsetDict(Mapping):
    dict = {}
    number_of_bits = 0

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


@pytest.fixture
def h5filename():
    fn = tmpname()
    yield fn
    if os.path.isfile(fn):
        os.remove(fn)


@pytest.fixture
def number_of_bits():
    return 8


@pytest.fixture
def bitsets():
    return MockedIntbitsetDict({
        'a': intbitset([1, 2, 3]),
        'b': intbitset([1, 2, 4, 5, 8]),
        'c': intbitset([1, 2, 4, 8])
    }, 8)


@pytest.fixture
def sample_pairs():
    return [
        ('a', 'c', 0.7523),
        ('b', 'c', 0.8342)]


@pytest.fixture
def id2label():
    return {
        1: 'a',
        2: 'b',
        3: 'c'
    }


@pytest.fixture
def label2id():
    return {
        'a': 1,
        'b': 2,
        'c': 3
    }


@pytest.fixture
def compact_pairs():
    return [(1, 3, 49301), (2, 3, 54669)]


class Testpairs(object):
    def test_dump_pairs_tsv(self, sample_pairs):
        out = StringIO()

        pairs.dump_pairs_tsv(sample_pairs, out)
        result = out.getvalue()

        expected = "a\tc\t0.7523\nb\tc\t0.8342\n"
        assert result == expected

    def fill_matrix(self, sample_pairs, label2id, h5filename):
        kripodb.pairs.dump_pairs_hdf5(sample_pairs,
                                      label2id,
                                      2,
                                      h5filename)

    def test_similar_run(self, sample_pairs, label2id, h5filename):
        self.fill_matrix(sample_pairs, label2id, h5filename)
        out = StringIO()

        pairs.similar_run('a', h5filename, 0.55, out)

        result = out.getvalue()
        expected = "a\tc\t0.7523\n"
        assert result == expected

    def test_similar_run_nohits(self, sample_pairs, label2id, h5filename):
        self.fill_matrix(sample_pairs, label2id, h5filename)
        out = StringIO()

        pairs.similar_run('a', h5filename, 0.99, out)

        result = out.getvalue()
        expected = ""
        assert result == expected

    def test_dump_pairs_ashdf5(self, sample_pairs, label2id, h5filename, compact_pairs):
        expectedrows = 2
        kripodb.pairs.dump_pairs_hdf5(sample_pairs,
                                      label2id,
                                      expectedrows,
                                      h5filename)

        h5file = tables.open_file(h5filename)
        mypairs = []
        for row in h5file.root.pairs:
            mypairs.append((row[0], row[1], row[2]))
        assert mypairs == compact_pairs
        h5file.close()

    def test_dump_pairs_badformat(self, bitsets, number_of_bits, label2id):
        with pytest.raises(LookupError) as cm:
            pairs.dump_pairs(bitsets,
                             bitsets,
                             'bikes',
                             'StringIO',
                             None,
                             number_of_bits,
                             0.4,
                             0.05,
                             label2id,
                             True
                             )

        assert cm.value.args == ('Invalid output format',)

    def test_dump_pairs_badffn(self, bitsets, number_of_bits, label2id):
        with pytest.raises(Exception) as cm:
            pairs.dump_pairs(bitsets,
                             bitsets,
                             'hdf5_compact',
                             '-',
                             None,
                             number_of_bits,
                             0.4,
                             0.05,
                             label2id,
                             True
                             )

        assert cm.value.args == ("hdf5 formats can't be outputted to stdout",)

    def test_dump_pairs_astsv(self, bitsets, number_of_bits, label2id):
        out = StringIO()

        pairs.dump_pairs(bitsets,
                         bitsets,
                         'tsv',
                         'StringIO',
                         out,
                         number_of_bits,
                         0.4,
                         0.05,
                         label2id,
                         False,
                         True,
                         )
        result = out.getvalue()

        expected = "a\tc\t0.13556\n"
        assert result == expected

    def test_dump_pairs_astsv_nomem(self, bitsets, number_of_bits, label2id):
        out = StringIO()

        pairs.dump_pairs(bitsets,
                         bitsets,
                         'tsv',
                         'StringIO',
                         out,
                         number_of_bits,
                         0.4,
                         0.05,
                         label2id,
                         True,
                         True,
                         )
        result = out.getvalue()

        expected = "a\tc\t0.13556\n"
        assert result == expected

    def test_similarity2query(self, bitsets):
        out = StringIO()

        pairs.similarity2query(bitsets,
                             'a',
                             out,
                             0.4,
                             0.05,
                             True
                             )
        result = out.getvalue()

        expected = "a\tc\t0.13556\n"
        assert result == expected

    def test_total_number_of_pairs(self, sample_pairs, label2id, h5filename):
        self.fill_matrix(sample_pairs, label2id, h5filename)

        result = pairs.total_number_of_pairs([h5filename])

        assert result == 2


def test_merge():
    infiles = [tmpname(), tmpname(), tmpname()]

    outfile = tmpname()
    try:
        # fill infiles
        inmatrix1 = SimilarityMatrix(infiles[0], 'w', 1, 2**16-1, 2)
        inmatrix1.update([('a', 'b', 0.2)], {'a': 1, 'b': 2, 'c': 3})
        inmatrix1.close()

        # matrix with same labels -> copy pairs table by dump/append, ignores labels tables
        inmatrix2 = SimilarityMatrix(infiles[1], 'w', 2, 2**16-1, 3)
        inmatrix2.update([('a', 'c', 0.6)], {'a': 1, 'b': 2, 'c': 3})
        inmatrix2.close()

        # matrix generated with different labels -> copy pairs table by iterate/update, adds missing labels
        inmatrix3 = SimilarityMatrix(infiles[2], 'w', 2, 2**16-1, 3)
        inmatrix3.update([('b', 'e', 0.4), ('e', 'f', 0.8)], {'b': 1, 'e': 2, 'f': 3})
        inmatrix3.close()

        pairs.merge(infiles, outfile)

        # compare it
        outmatrix = SimilarityMatrix(outfile)
        result = list(outmatrix)
        outmatrix.close()
        expected = [('a', 'b', 0.2), ('a', 'c', 0.6), ('b', 'e', 0.4), ('e', 'f', 0.8)]
        assert result == expected
    finally:
        for infile in infiles:
            if os.path.isfile(infile):
                os.remove(infile)
        if os.path.isfile(outfile):
            os.remove(outfile)
