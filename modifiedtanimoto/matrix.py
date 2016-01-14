# Copyright 2013 Netherlands eScience Center
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

import sys
import logging
import tables

from algorithm import distance, corrections
from dbm import load_bitsets
from id2label import read_id2label, swap_label2id


def dump_matrix(bs_file1, bs_file2,
                out_format, out_file,
                size,
                mean_onbit_density,
                cutoff,
                id2label_file):

    logging.warn('Reading {}'.format(bs_file1))
    bitsets1 = load_bitsets(bs_file1)
    if bs_file1 == bs_file2:
        bitsets2 = bitsets1
    else:
        logging.warn('Reading {}'.format(bs_file2))
        bitsets2 = load_bitsets(bs_file2)

    out = sys.stdout
    if out_format.beginswith('tsv') and out_file != '-':
        out = open(out_file)

    (corr_st, corr_sto) = corrections(mean_onbit_density)
    if out_format == 'tsv':
        dump_matrix_tsv(bitsets1, bitsets2, size, corr_st, corr_sto, cutoff, out)
    elif out_format == 'tsv_numbered':
        label2id = swap_label2id(read_id2label(id2label_file))
        dump_matrix_tsv_numbered(bitsets1, bitsets2, size, corr_st, corr_sto, cutoff, label2id, out)
    elif out_format == 'hdf5':
        dump_matrix_hdf5(bitsets1, bitsets2, size, corr_st, corr_sto, cutoff, out_file)
    else:
        raise LookupError('Invalid output format')


def dump_matrix_tsv(bitsets1, bitsets2, size, corr_st, corr_sto, cutoff, out):
    """

    Pro:
    * when stored in sqlite can be used outside of Python
    Con:
    * big


    Parameters
    ----------
    bitsets1
    bitsets2
    size
    corr_st
    corr_sto
    cutoff
    out

    Returns
    -------

    """
    for (label1, bs1) in bitsets1.iteritems():
        for (label2, bs2) in bitsets2.iteritems():
            if label1 >= label2:
                continue

            d = distance(bs1, bs2, size, corr_st, corr_sto)
            if d >= cutoff:
                out.write('{}\t{}\t{}\n'.format(label1, label2, d))


def dump_matrix_tsv_numbered(bitsets1, bitsets2, size, corr_st, corr_sto, cutoff, label2id, out):
    """

    Pro:
    * more compact, because label string is replaced with a integer
    * when stored in sqlite can be used outside of Python
    Con:
    * Requires a lookup tabel

    Parameters
    ----------
    bitsets1: dictionary of bitset identifier as key and a intbitset object as value
    bitsets2: dictionary of bitset identifier as key and a intbitset object as value
    size: size of bitset
    corr_st: correction value for modified tanimoto
    corr_sto: correction value for modified tanimoto
    cutoff: distance greater than cutoff are allowed
    id2label_file: name of file with id to label lookup in tsv format
    out: name of output file

    Returns
    -------

    """

    for (label1, bs1) in bitsets1.iteritems():
        for (label2, bs2) in bitsets2.iteritems():
            if label1 >= label2:
                continue

            d = distance(bs1, bs2, size, corr_st, corr_sto)
            if d >= cutoff:
                id1 = label2id[label1]
                id2 = label2id[label2]
                out.write('{}\t{}\t{}\n'.format(id1, id2, d))


class Pair(tables.IsDescription):
    a = tables.StringCol(15)
    b = tables.StringCol(15)
    score = tables.Float32Col()


def dump_matrix_hdf5(bitsets1, bitsets2, size, corr_st, corr_sto, cutoff, out_file):
    """

    Pro:
    * small
    Con:
    * requires hdf5 library to access


    Parameters
    ----------
    bitsets1
    bitsets2
    size
    corr_st
    corr_sto
    cutoff
    out_file

    Returns
    -------

    """
    filters = tables.Filters(complevel=6, complib='blosc')
    h5file = tables.open_file(out_file, mode='w', filters=filters)
    group = h5file.create_group('/', 'matrix', 'Distance matrix')
    table = h5file.create_table(group, 'pairs', Pair, 'Distance matrix pairs')
    hit = table.row
    for (label1, bs1) in bitsets1.iteritems():
        for (label2, bs2) in bitsets2.iteritems():
            if label1 >= label2:
                continue

            d = distance(bs1, bs2, size, corr_st, corr_sto)
            if d >= cutoff:
                hit['a'] = label1
                hit['b'] = label2
                hit['score'] = d
                hit.append()
        table.flush()
    table.cols.a.create_index(filters=filters)
    table.cols.b.create_index(filters=filters)
    h5file.close()
