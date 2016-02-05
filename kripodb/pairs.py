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
"""Module handling generation and retrieval of distance matrix"""

import logging

from .hdf5 import DistanceMatrix
from .modifiedtanimoto import distances, corrections


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
    """Dump pairs in tab delimited file

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
    """Dump pairs in hdf5 file

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
    """Calculate distance of query against all fingerprints in bitsets2 and write to tab delimited file.

    Args:
        bitsets2 (kripodb.db.IntbitsetDict):
        query (str): Query identifier or beginning of it
        out (File): File object to write output to
        mean_onbit_density (flaot): Mean on bit density
        cutoff (float): Cutoff, distance scores below cutoff are discarded.
        memory (Optional[bool]): When true will load bitset2 into memory, when false it doesn't

    """
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
    """Find similar fragments to query based on distance matrix and write to tab delimited file.

    Args:
        query (str): Query fragment identifier
        pairsdbfn (str): Filename of distance matrix file
        cutoff (float): Cutoff, distance scores below cutoff are discarded.
        out (File): File object to write output to

    """
    matrix = DistanceMatrix(pairsdbfn)
    pairs = matrix.pairs()
    labels = matrix.labels()

    hits = similar(query, pairs, labels, cutoff)
    dump_pairs_tsv(hits, out)

    matrix.close()


def similar(query, pairsdb, labels, cutoff):
    """Find similar fragments to query based on distance matrix.

    Args:
        query (str): Query fragment identifier
        pairsdb (kripodb.db.PairsTable): Pairs table
        labels (kripodb.db.LabelsLookup): Labels lookup table
        cutoff (float): Cutoff, distance scores below cutoff are discarded.

    Returns:
        List[(str, str, float)]: List of (query fragment identifier, hit fragment identifier, distance score) sorted on distance score

    """
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
    """Concatenate distance matrix files into a single one.

    Args:
        ins (list[str]): List of input distance matrix filenames
        out (str):  Output distance matrix filenames

    Raises:
        AssertionError: When nr of labels of input files is not the same

    """
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
