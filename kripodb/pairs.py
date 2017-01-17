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
"""Module handling generation and retrieval of similarity of fingerprint pairs"""

from __future__ import absolute_import

import tables

import logging
from kripodb.frozen import FrozenSimilarityMatrix

from .hdf5 import SimilarityMatrix
from .modifiedtanimoto import similarities, corrections
from .webservice.client import WebserviceClient


def dump_pairs(bitsets1,
               bitsets2,
               out_format,
               out_file,
               out,
               number_of_bits,
               mean_onbit_density,
               cutoff,
               label2id,
               nomemory,
               ignore_upper_triangle=False):
    """Dump pairs of bitset collection.

    A pairs are rows of the bitset identifier of both bitsets with a similarity score.

    Args:
        bitsets1 (Dict{str, intbitset.intbitset}): First dict of fingerprints
            with fingerprint label as key and intbitset as value
        bitsets2 (Dict{str, intbitset.intbitset}): Second dict of fingerprints
            with fingerprint label as key and intbitset as value
        out_format: 'tsv' or 'hdf5'
        out_file: Filename of output file where 'hdf5' format is written to.
        out (File): File object where 'tsv' format is written to.
        number_of_bits (int): Number of bits for all bitsets
        mean_onbit_density (float): Mean on bit density
        cutoff (float): Cutoff, similarity scores below cutoff are discarded.
        label2id: dict to translate label to id (string to int)
        nomemory: If true bitset2 is not loaded into memory
        ignore_upper_triangle: When true returns similarity where label1 > label2,
            when false returns all similarities

    """
    if out_file == '-' and out_format.startswith('hdf5'):
        raise Exception("hdf5 formats can't be outputted to stdout")

    if not nomemory:
        # load whole dict in memory so it can be reused for each bitset1
        # deserialization of bitsets2 is only done one time
        bitsets2 = bitsets2.materialize()

    expectedrows = len(bitsets1) * len(bitsets2) * cutoff * 0.025

    (corr_st, corr_sto) = corrections(mean_onbit_density)

    logging.warning('Generating pairs')

    similarities_iter = similarities(bitsets1, bitsets2,
                                     number_of_bits, corr_st, corr_sto,
                                     cutoff,
                                     ignore_upper_triangle)

    if out_format == 'tsv':
        dump_pairs_tsv(similarities_iter, out)
    elif out_format == 'hdf5':
        dump_pairs_hdf5(similarities_iter,
                        label2id,
                        expectedrows,
                        out_file)
    else:
        raise LookupError('Invalid output format')


def dump_pairs_tsv(similarities_iter, out):
    """Dump pairs in tab delimited file

    Pro:
    * when stored in sqlite can be used outside of Python
    Con:
    * big, unless output is compressed

    Args:
        similarities_iter (Iterator): Iterator with tuple with fingerprint 1 label, fingerprint 2 label, similarity as members
        out (File): Writeable file
    """
    for label1, label2, similarity in similarities_iter:
        out.write('{0}\t{1}\t{2:.5}\n'.format(label1, label2, similarity))


def dump_pairs_hdf5(similarities_iter,
                    label2id,
                    expectedrows,
                    out_file):
    """Dump pairs in hdf5 file

    Pro:
    * very small, 10 bytes for each pair + compression
    Con:
    * requires hdf5 library to access

    Args:
        similarities_iter (Iterator): Iterator with tuple with fingerprint 1 label, fingerprint 2 label, similarity as members
        label2id (dict): dict to translate label to id (string to int)
        expectedrows:
        out_file:

    """
    matrix = SimilarityMatrix(out_file, 'w',
                              expectedpairrows=expectedrows,
                              expectedlabelrows=len(label2id))

    matrix.update(similarities_iter, label2id)

    matrix.close()


def similarity2query(bitsets2, query, out, mean_onbit_density, cutoff, memory):
    """Calculate similarity of query against all fingerprints in bitsets2 and write to tab delimited file.

    Args:
        bitsets2 (kripodb.db.IntbitsetDict):
        query (str): Query identifier or beginning of it
        out (File): File object to write output to
        mean_onbit_density (flaot): Mean on bit density
        cutoff (float): Cutoff, similarity scores below cutoff are discarded.
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

    similarities_iter = similarities(bitsets1, bitsets2,
                               number_of_bits, corr_st, corr_sto,
                               cutoff, True)
    sorted_similarities = sorted(similarities_iter, key=lambda row: row[2], reverse=True)
    dump_pairs_tsv(sorted_similarities, out)


def similar_run(query, pairsdbfn, cutoff, out):
    """Find similar fragments to query based on similarity matrix and write to tab delimited file.

    Args:
        query (str): Query fragment identifier
        pairsdbfn (str): Filename of similarity matrix file or url of kripodb webservice
        cutoff (float): Cutoff, similarity scores below cutoff are discarded.
        out (File): File object to write output to

    """
    if pairsdbfn.startswith('http'):
        client = WebserviceClient(pairsdbfn)
        hits = client.similar_fragments(query, cutoff)
        hits = [(h['query_frag_id'], h['hit_frag_id'], h['score']) for h in hits]
        dump_pairs_tsv(hits, out)
    else:
        matrix = open_similarity_matrix(pairsdbfn)
        hits = similar(query, matrix, cutoff)
        dump_pairs_tsv(hits, out)
        matrix.close()


def open_similarity_matrix(fn):
    """Open read-only similarity matrix file.

    Args:
        fn (str): Filename of similarity matrix

    Returns:
        SimilarityMatrix | FrozenSimilarityMatrix: A read-only similarity matrix object

    """
    # peek in file to detect format
    f = tables.open_file(fn, 'r')
    is_frozen = 'scores' in f.root
    f.close()
    if is_frozen:
        matrix = FrozenSimilarityMatrix(fn)
    else:
        matrix = SimilarityMatrix(fn, cache_labels=True)
    return matrix


def similar(query, similarity_matrix, cutoff, limit=None):
    """Find similar fragments to query based on similarity matrix.

    Args:
        query (str): Query fragment identifier
        similarity_matrix (kripodb.db.SimilarityMatrix): Similarity matrix
        cutoff (float): Cutoff, similarity scores below cutoff are discarded.
        limit (int): Maximum number of hits. Default is None for no limit.

    Yields:
        Tuple[(str, str, float)]: List of (query fragment identifier, hit fragment identifier, similarity score) sorted on similarity score

    """
    raw_hits = similarity_matrix.find(query, cutoff, limit)
    # add query column
    for hit_id, score in raw_hits:
        yield query, hit_id, score


def total_number_of_pairs(fingerprint_filenames):
    """Count number of pairs in similarity matrix files

    Args:
        fingerprint_filenames (list[str]): List of file names of similarity matrices

    Returns:
        int: Total number of pairs
    """
    sizes = []
    for filename in fingerprint_filenames:
        matrix = SimilarityMatrix(filename)
        pairs = matrix.pairs
        sizes.append(len(pairs))
        matrix.close()
    return sum(sizes)


def merge(ins, out):
    """Concatenate similarity matrix files into a single one.

    Args:
        ins (list[str]): List of input similarity matrix filenames
        out (str):  Output similarity matrix filenames

    Raises:
        AssertionError: When nr of labels of input files is not the same

    """
    expectedrows = total_number_of_pairs(ins)
    out_matrix = SimilarityMatrix(out, 'w', expectedpairrows=expectedrows)

    # copy pairs
    for in_filename in ins:
        in_matrix = SimilarityMatrix(in_filename)
        out_matrix.append(in_matrix)
        in_matrix.close()

    out_matrix.close()
