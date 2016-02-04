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

from kripodb.hdf5 import DistanceMatrix
from modifiedtanimoto import distances, corrections


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
    elif out_format == 'hdf5':
        dump_pairs_hdf5(distances_iter,
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


def dump_pairs_hdf5(distances_iter,
                    label2id,
                    precision,
                    expectedrows,
                    out_file):
    """

    Pro:
    * very small, 10 bytes for each pair + compression
    Con:
    * requires hdf5 library to access

    :param distances_iter:
    :param label2id: dict to translate label to id (string to int)
    :param precision:
    :param expectedrows:
    :param out_file:
    :return:
    """
    matrix = DistanceMatrix(out_file, 'w')

    pairs = matrix.pairs(expectedrows, precision)
    pairs.update(distances_iter, label2id)
    pairs.add_indexes()

    labels = matrix.labels(len(label2id))
    labels.update(label2id)

    matrix.close()


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


def similar_run(query, pairsdbfn, cutoff, out):
    matrix = DistanceMatrix(pairsdbfn)
    pairs = matrix.pairs()
    labels = matrix.labels()

    hits = similar(query, pairs, labels, cutoff)
    dump_pairs_tsv(hits, out)

    matrix.close()


def similar(query, pairsdb, labels, cutoff):
    frag_id = labels.by_label(query)
    raw_hits = pairsdb.find(frag_id, cutoff)

    # replace ids with labels + add query column
    hits = []
    for hit_id, score in raw_hits.iteritems():
        hit = (query, labels.by_id(hit_id), score)
        hits.append(hit)

    # highest score/most similar first
    sorted_hits = sorted(hits, reverse=True, key=lambda r: r[2])

    return sorted_hits


def total_number_of_pairs(fingerprintfilenames):
    sizes = []
    for filename in fingerprintfilenames:
        matrix = DistanceMatrix(filename)
        pairs = matrix.pairs()
        sizes.append(len(pairs))
        matrix.close()
    return sum(sizes)


def labels_consistency_check(fingerprintfilenames):
    nr_labels = set()
    for filename in fingerprintfilenames:
        matrix = DistanceMatrix(filename)
        labels = matrix.labels()
        nr_labels.add(len(labels))
        # TODO implement more checks
        matrix.close()
    assert len(nr_labels) == 1


def merge(ins, out):
    expectedrows = total_number_of_pairs(ins)
    labels_consistency_check(ins)
    out_matrix = DistanceMatrix(out, 'w')

    # copy labels
    first_in_matrix = DistanceMatrix(ins[0])
    first_in_labels = first_in_matrix.labels()
    out_labels = out_matrix.labels(len(first_in_labels))
    out_labels.append(first_in_labels)

    # copy precision
    first_in_pairs = first_in_matrix.pairs()
    out_pairs = out_matrix.pairs(expectedrows, first_in_pairs.score_precision)

    first_in_matrix.close()

    # copy pairs
    for filename in ins:
        matrix = DistanceMatrix(filename)
        pairs = matrix.pairs()
        out_pairs.append(pairs)
        matrix.close()

    out_pairs.add_indexes()

    out_matrix.close()
