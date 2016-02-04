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

import logging
import tables

from modifiedtanimoto import distances, corrections
from kripodb.db import FingerprintsDb, FragmentsDb


def dump_pairs(bitsets1,
               bitsets2,
               out_format,
               out_file,
               out,
               number_of_bits,
               mean_onbit_density,
               cutoff,
               label2id,
               precision,
               nomemory):
    """Dump pairs of bitset collection

    :param bitsets1: dictionary of bitset identifier as key
        and a intbitset object as value
    :param bitsets2: dictionary of bitset identifier as key
        and a intbitset object as value
    :param out_format:
    :param out_file:
    :param out:
    :param number_of_bits: Maximum number of bits in bitset
    :param mean_onbit_density:
    :param cutoff:
    :param label2id: dict to translate label to id (string to int)
    :param precision:
    :param nomemory: If true bitset2 is not loaded into memory
    :return:
    """
    if out_file == '-' and out_format.startswith('hdf5'):
        raise Exception("hdf5 formats can't be outputted to stdout")

    if not nomemory:
        # load whole dict in memory so it can be reused for each bitset1
        # deserialization of bitsets2 is only done one time
        bitsets2 = bitsets2.materialize()

    expectedrows = len(bitsets1) * len(bitsets2) * cutoff * 0.025

    (corr_st, corr_sto) = corrections(mean_onbit_density)

    logging.warn('Generating pairs')

    distances_iter = distances(bitsets1, bitsets2,
                               number_of_bits, corr_st, corr_sto,
                               cutoff)

    if out_format == 'tsv':
        dump_pairs_tsv(distances_iter, out)
    elif out_format == 'hdf5_compact':
        dump_pairs_hdf5_compact(distances_iter,
                                label2id,
                                precision,
                                expectedrows,
                                out_file)
    else:
        raise LookupError('Invalid output format')


def dump_pairs_tsv(distances_iter, out):
    """Dump pairs as

    Pro:
    * when stored in sqlite can be used outside of Python
    Con:
    * big, unless output is compressed

    :param distances_iter:
    :param out:
    :return:

    """
    for label1, label2, distance in distances_iter:
        out.write('{}\t{}\t{}\n'.format(label1, label2, distance))


class PairCompact(tables.IsDescription):
    a = tables.UInt32Col()
    b = tables.UInt32Col()
    score = tables.UInt16Col()


def dump_pairs_hdf5_compact(distances_iter,
                            label2id,
                            precision,
                            expectedrows,
                            out_file):
    """

    Pro:
    * very small, 9 bytes for each pair
    * index on pair ids
    Con:
    * requires hdf5 library to access
    * Requires a lookup table

    :param distances_iter:
    :param label2id: dict to translate label to id (string to int)
    :param precision:
    :param expectedrows:
    :param out_file:
    :return:
    """
    filters = tables.Filters(complevel=6, complib='blosc')
    h5file = tables.open_file(out_file, mode='w', filters=filters)
    table = h5file.create_table('/',
                                'pairs',
                                PairCompact,
                                'Distance pairs',
                                expectedrows=expectedrows)
    table.attrs['score_precision'] = precision
    hit = table.row
    for label1, label2, distance in distances_iter:
        hit['a'] = label2id[label1]
        hit['b'] = label2id[label2]
        hit['score'] = int(distance * precision)
        hit.append()
    table.cols.a.create_index(filters=filters)
    table.cols.b.create_index(filters=filters)
    h5file.close()


class Id2Label(tables.IsDescription):
    frag_id = tables.UInt32Col()
    label = tables.StringCol(16)


def add_lookup2h5(h5file, filters, label2id):
    expectedrows = len(label2id)
    table = h5file.create_table('/',
                                'labels',
                                Id2Label,
                                'Labels lookup',
                                expectedrows=expectedrows)

    for label, frag_id in label2id.iteritems():
        table.row['frag_id'] = frag_id
        table.row['label'] = label
        table.row.append()

    table.cols.frag_id.create_index(filters=filters)
    table.cols.label.create_index(filters=filters)
    table.flush()
    return table


def distance2query(bitsets2, query, out, mean_onbit_density, cutoff, memory):
    number_of_bits = bitsets2.number_of_bits
    if query in bitsets2:
        # exact match
        query_bitset = bitsets2[query]
        bitsets1 = {
            query: query_bitset
        }
    else:
        # all bitsets which have a key that starts with query
        bitsets1 = {k: v for k, v in bitsets2.iteritems_startswith(query)}

        if memory:
            # load whole dict in memory so it can be reused for each bitset1
            # deserialization of bitset2 is only done one time
            bitsets2 = bitsets2.materialize()

    (corr_st, corr_sto) = corrections(mean_onbit_density)

    distances_iter = distances(bitsets1, bitsets2,
                               number_of_bits, corr_st, corr_sto,
                               cutoff, True)
    sorted_distances = sorted(distances_iter, key=lambda row: row[2], reverse=True)
    dump_pairs_tsv(sorted_distances, out)


def similar_run(query, pairsdbfn, fragments, cutoff, out, memory):
    id2label = fragments.id2label()
    if memory:
        # when there are many hits,
        # it is more efficient to fetch all with a single query
        # instead of many queries
        id2label = id2label.materialize()
    frag_id = fragments.label2id()[query]
    h5file = tables.open_file(pairsdbfn)
    pairs = h5file.root.pairs

    hits = similar(frag_id, pairs, id2label, cutoff)
    dump_pairs_tsv(hits, out)

    h5file.close()


def similar(frag_id, pairsdb, id2label, cutoff):
    hits = []

    query = id2label[frag_id]
    precision = float(pairsdb.attrs['score_precision'])
    scutoff = int(cutoff * precision)

    query1 = '(a == {}) & (score >= {})'.format(frag_id, scutoff)
    for row in pairsdb.where(query1):
        score = row[2] / precision
        label = id2label[row[1]]
        hits.append((query, label, score))

    query2 = '(b == {}) & (score >= {})'.format(frag_id, scutoff)
    for row in pairsdb.where(query2):
        score = row[2] / precision
        label = id2label[row[0]]
        hits.append((query, label, score))

    # most similar first
    sorted_hits = sorted(hits, reverse=True, key=lambda r: (r[1], r[2]))

    return sorted_hits
