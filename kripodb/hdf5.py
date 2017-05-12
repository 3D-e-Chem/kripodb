# Copyright 2016 Netherlands eScience Center
#
# Licensed under the Apache License, Version 2.0 (the 'License");
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
"""Similarity matrix using hdf5 as storage backend."""
from __future__ import absolute_import

from math import log10, ceil, floor

import numpy as np
from progressbar import ProgressBar
import tables
import six


class SimilarityMatrix(object):
    """Similarity matrix

    Args:
        filename (str): File name of hdf5 file to write or read similarity matrix from
        mode (str): Can be 'r' for reading or 'w' for writing
        expectedpairrows (int): Expected number of pairs to be added.
            Required when similarity matrix is opened in write mode, helps optimize storage
        expectedlabelrows (int): Expected number of labels to be added.
            Required when similarity matrix is opened in write mode, helps optimize storage
        cache_labels (bool): Cache labels, speed up label lookups

    Attributes:
        h5file (tables.File): Object representing an open hdf5 file
        pairs (PairsTable): HDF5 Table that contains pairs
        labels (LabelsLookup): Table to look up label of fragment by id or id of fragment by label
    """
    filters = tables.Filters(complevel=6, complib='blosc')

    def __init__(self, filename, mode='r', expectedpairrows=None, expectedlabelrows=None, cache_labels=False, **kwargs):
        self.h5file = tables.open_file(filename, mode, filters=self.filters, **kwargs)
        self.pairs = PairsTable(self.h5file, expectedpairrows)
        self.labels = LabelsLookup(self.h5file, expectedlabelrows)
        self.cache_i2l = {}
        self.cache_l2i = {}
        if cache_labels:
            self._build_label_cache()

    def _build_label_cache(self):
        if not self.cache_l2i:
            self.cache_l2i = self.labels.label2ids()
            self.cache_i2l = {v: k for k, v in six.iteritems(self.cache_l2i)}

    def close(self):
        """Closes the hdf5file"""
        self.h5file.close()

    def append(self, other):
        """Append data from other similarity matrix to me

        Args:
            other (SimilarityMatrix): Other similarity matrix
        """
        if len(self.labels) == 0:
            # copy labels when self has no labels
            self.labels.append(other.labels)

        if self.labels == other.labels:
            # same id 2 labels mapping, can safely copy pairs
            self.pairs.append(other.pairs)
        else:
            # different id 2 labels mapping, must remap labels of other into labels of self
            # for labels missing in self, generate new id and add to self
            self.labels.merge(other.labels)
            # for other pairs map a,b to labels of other and map other labels to self id using self labels
            self.pairs.update(other, self.labels.label2ids())

    def __iter__(self):
        self._build_label_cache()
        for pair in self.pairs:
            yield self.cache_i2l[pair['a']], self.cache_i2l[pair['b']], pair['score']

    def update(self, similarities_iter, label2id):
        """Store pairs of fragment identifier with their similarity score and label 2 id lookup

        Args:
            similarities_iter (iterator): Iterator which yields (label1, label2, similarity_score)
            label2id (dict): Dictionary with fragment label as key and fragment identifier as value.

        """
        self.pairs.update(similarities_iter, label2id)
        self.labels.update(label2id)

    def find(self, query, cutoff, limit=None):
        """Find similar fragments to query.

        Args:
            query (str): Query fragment identifier
            cutoff (float): Cutoff, similarity scores below cutoff are discarded.
            limit (int): Maximum number of hits. Default is None for no limit.

        Yields:
            (str, float): Hit fragment idenfier and similarity score
        """
        if self.cache_l2i:
            frag_id = self.cache_l2i[query]
            for hit_frag_id, score in self.pairs.find(frag_id, cutoff, limit):
                yield self.cache_i2l[hit_frag_id], score
        else:
            frag_id = self.labels.by_label(query)
            for hit_frag_id, score in self.pairs.find(frag_id, cutoff, limit):
                yield self.labels.by_id(hit_frag_id), score

    def count(self, frame_size, raw_score=False, lower_triangle=False):
        """Count occurrences of each score

        Args:
            frame_size (int): Size of matrix loaded each time. Larger requires more memory and smaller is slower.
            raw_score (bool): Return raw int16 score or fraction score
            lower_triangle (bool): Dummy argument to force same interface for thawed and frozen matrix

        Returns:
            (str, int): Score and number of occurrences
        """
        return self.pairs.count(frame_size, raw_score)

    def keep(self, other, keep):
        """Copy content of self to other and only keep given fragment labels and the labels they pair with

        Args:
            other (SimilarityMatrix): Writable matrix to fill
            keep (set[str]): Fragment labels to keep
        """
        frag_ids2keep = self.labels.by_labels(keep)
        all_frag_ids2keep = self.pairs.keep(other.pairs, frag_ids2keep)
        self.labels.keep(other.labels, all_frag_ids2keep)

    def skip(self, other, skip):
        """Copy content of self to other and skip all given fragment labels

        Args:
            other (SimilarityMatrix): Writable matrix to fill
            skip (set[str]): Fragment labels to skip
        """
        frag_ids2skip = self.labels.by_labels(skip)
        self.pairs.skip(other.pairs, frag_ids2skip)
        self.labels.skip(other.labels, frag_ids2skip)


class AbstractSimpleTable(object):
    """Abstract wrapper around a HDF5 table

    Args:
        table (tables.Table): HDF5 table
        append_chunk_size (int): Size of chunk to append in one go.
            Defaults to 1e8, which when table description is 10bytes will require 2Gb during append.

    Attributes
        table (tables.Table): HDF5 table
        append_chunk_size (int): Number of rows to read from other table during append.

    """

    def __init__(self, table, append_chunk_size=int(1e8)):
        self.table = table
        self.append_chunk_size = append_chunk_size

    def append(self, other):
        """Append rows of other table to self

        Args:
            other: Table of same type as self

        """
        limit = len(other.table)
        for start in range(0, limit, self.append_chunk_size):
            stop = self.append_chunk_size + start
            self.table.append(other.table.read(start=start, stop=stop))

    def __len__(self):
        """

        Returns:
            int: Number of rows in table
        """
        return len(self.table)

    def __iter__(self):
        return self.table.__iter__()


class SimilarityPair(tables.IsDescription):
    """Table description for similarity pair"""
    a = tables.UInt32Col()
    b = tables.UInt32Col()
    score = tables.UInt16Col()


class PairsTable(AbstractSimpleTable):
    """Tabel to store similarity score of a pair of fragment fingerprints

    When table does not exist in h5file it is created.

    Args:
        h5file (tables.File): Object representing an open hdf5 file
        expectedrows (int): Expected number of pairs to be added.
            Required when similarity matrix is opened in write mode, helps optimize storage

    Attributes:
        score_precision (int): Similarity score is a fraction,
            the score is converted to an int by multiplying it with the precision
        full_matrix (bool): Matrix is filled above and below diagonal.
    """
    table_name = 'pairs'
    filters = tables.Filters(complevel=6, complib='blosc')

    def __init__(self, h5file, expectedrows=0):
        if self.table_name in h5file.root:
            table = h5file.root.__getattr__(self.table_name)
        else:
            table = h5file.create_table('/',
                                        self.table_name,
                                        SimilarityPair,
                                        'Similarity pairs',
                                        expectedrows=expectedrows)

        super(PairsTable, self).__init__(table)
        self.score_precision = 2 ** 16 - 1

    @property
    def full_matrix(self):
        try:
            return self.table.attrs['full_matrix']
        except KeyError:
            return False

    @full_matrix.setter
    def full_matrix(self, value):
        self.table.attrs['full_matrix'] = value

    def update(self, similarities_iter, label2id):
        """Store pairs of fragment identifier with their similarity score

        Args:
            similarities_iter (Iterator): Iterator which yields (label1, label2, similarity_score)
            label2id (Dict): Lookup with fragment label as key and fragment identifier as value

        """
        hit = self.table.row
        for label1, label2, similarity in similarities_iter:
            hit['a'] = label2id[label1]
            hit['b'] = label2id[label2]
            hit['score'] = int(similarity * self.score_precision)
            hit.append()
        self.table.flush()

    def find(self, frag_id, cutoff, limit):
        """Find fragment hits which has a similarity score with frag_id above cutoff.

        Args:
            frag_id (int): query fragment identifier
            cutoff (float): Cutoff, similarity scores below cutoff are discarded.
            limit (int): Maximum number of hits. Default is None for no limit.

        Returns:
            List[Tuple]: Where first tuple value is hit fragment identifier and second value is similarity score

        """
        precision = float(self.score_precision)
        precision10 = float(10**(floor(log10(precision))))
        scutoff = int(cutoff * precision)

        hits = {}
        query1 = '(a == {0}) & (score >= {1})'.format(frag_id, scutoff)
        for row in self.table.where(query1):
            hit_id, score = row[1], row[2]
            score = ceil(precision10 * score / precision) / precision10
            hits[hit_id] = score

        if not self.full_matrix:
            query2 = '(b == {0}) & (score >= {1})'.format(frag_id, scutoff)
            for row in self.table.where(query2):
                hit_id, score = row[0], row[2]
                score = ceil(precision10 * score / precision) / precision10
                hits[hit_id] = score

        # highest score==most similar first
        sorted_hits = sorted(six.iteritems(hits), reverse=True, key=lambda r: r[1])

        if limit is not None:
            sorted_hits = sorted_hits[:limit]

        return sorted_hits

    def append(self, other):
        """Append rows of other table to self

        Args:
            other: Table of same type as self

        """
        super(PairsTable, self).append(other)

        if self.score_precision is None:
            self.score_precision = other.score_precision

    def __iter__(self):
        precision = float(self.score_precision)
        precision10 = float(10**(floor(log10(precision))))
        for pair in super(PairsTable, self).__iter__():
            score = ceil(precision10 * pair['score'] / precision) / precision10
            yield {'a': pair['a'], 'b': pair['b'], 'score': score}

    def count(self, frame_size, raw_score=False):
        """Count occurrences of each score

        Args:
            frame_size (int): Size of matrix loaded each time. Larger requires more memory and smaller is slower.
            raw_score (bool): Return raw int16 score or fraction score

        Returns:
            Tuple[(str, int)]: Score and number of occurrences
        """
        # Count occurrences of each score by reading pairs table in frames
        nr_rows = len(self.table)
        nr_bins = self.score_precision + 1
        counts = np.zeros(shape=nr_bins, dtype=np.int64)
        bar = ProgressBar()
        for start in bar(six.moves.range(0, nr_rows, frame_size)):
            stop = frame_size + start
            frame = self.table.read(start=start, stop=stop)
            frame_counts = np.bincount(frame['score'], minlength=nr_bins)
            counts += frame_counts

        if raw_score:
            for raw_score in counts.nonzero()[0]:
                yield (raw_score, counts[raw_score])
        else:
            # Convert int score into fraction
            precision = float(self.score_precision)
            precision10 = float(10 ** (ceil(log10(precision))))
            for raw_score in counts.nonzero()[0]:
                score = ceil(precision10 * raw_score / precision) / precision10
                count = counts[raw_score]
                yield (score, count)

    def keep(self, other, keep):
        """Copy pairs from self to other and keep given fragment identifiers and the identifiers they pair with.

        Args:
            other (PairsTable): Pairs table to fill
            keep (set[int]): Fragment identifiers to keep

        Returns:
            set[int]: Fragment identifiers that have been copied to other
        """
        all_frags2keep = set(keep)
        hit = other.table.row
        for row in self.table:
            if row['a'] in keep and row['b'] in keep:
                hit['a'] = row['a']
                hit['b'] = row['b']
                hit['score'] = row['score']
                hit.append()
            elif row['a'] in keep:
                hit['a'] = row['a']
                hit['b'] = row['b']
                hit['score'] = row['score']
                hit.append()
                all_frags2keep.add(row['b'])
            elif row['a'] in keep:
                hit['a'] = row['a']
                hit['b'] = row['b']
                hit['score'] = row['score']
                hit.append()
                all_frags2keep.add(row['a'])
        other.table.flush()
        return all_frags2keep

    def skip(self, other, skip):
        """Copy content from self to other and skip given fragment identifiers

        Args:
            other (PairsTable): Pairs table to fill
            skip (set[int]): Fragment identifiers to skip
        """
        hit = other.table.row
        for row in self.table:
            if row['a'] not in skip and row['b'] not in skip:
                hit['a'] = row['a']
                hit['b'] = row['b']
                hit['score'] = row['score']
                hit.append()
        other.table.flush()


class Id2Label(tables.IsDescription):
    """Table description of id 2 label table."""
    frag_id = tables.UInt32Col()
    label = tables.StringCol(16)


class LabelsLookup(AbstractSimpleTable):
    """Table to look up label of fragment by id or id of fragment by label

    When table does not exist in h5file it is created.

    Args:
        h5file (tables.File): Object representing an open hdf5 file
        expectedrows (int): Expected number of pairs to be added.
            Required when similarity matrix is opened in write mode, helps optimize storage
    """
    table_name = 'labels'
    filters = tables.Filters(complevel=6, complib='blosc')

    def __init__(self, h5file, expectedrows=0):
        if self.table_name in h5file.root:
            table = h5file.root.__getattr__(self.table_name)
        else:
            table = h5file.create_table('/',
                                        self.table_name,
                                        Id2Label,
                                        'Labels lookup',
                                        expectedrows=expectedrows)
            table.cols.frag_id.create_index(filters=self.filters)
            table.cols.label.create_index(filters=self.filters)

        super(LabelsLookup, self).__init__(table)

    def by_id(self, frag_id):
        """Look up label of fragment by id

        Args:
            frag_id (int): Fragment identifier

        Raises:
            IndexError: When id of fragment is not found

        Returns:
            str: Label of fragment
        """
        query = 'frag_id == {}'.format(frag_id)
        result = self.table.where(query)
        first_row = next(result)
        label_col = first_row[1].decode()
        return label_col

    def by_label(self, label):
        """Look up id of fragment by label

        Args:
            label (str): Fragment label

        Raises:
            IndexError: When label of fragment is not found

        Returns:
            int: Fragment identifier
        """
        query = 'label == b"{}"'.format(label)
        result = self.table.where(query)
        first_row = next(result)
        id_col = first_row[0]
        return id_col

    def by_labels(self, labels):
        """Look up ids of fragments by label

        Args:
            labels (set[str]): Set of fragment labels

        Raises:
            IndexError: When label of fragment is not found

        Returns:
            set[int]: Set of fragment identifiers
        """
        ids = set()
        for frag_label, frag_id in six.iteritems(self.label2ids()):
            if frag_label in labels:
                ids.add(frag_id)
        return ids

    def label2ids(self):
        """Return whole table as a dictionary

        Returns:
            dict: Dictionary with label as key and frag_id as value.

        """
        return {r['label'].decode(): r['frag_id'] for r in self.table}

    def update(self, label2id):
        """Update labels lookup by adding labels in label2id.

        Args:
            label2id (dict): Dictionary with fragment label as key and fragment identifier as value.

        """
        for label, frag_id in six.iteritems(label2id):
            self.table.row['frag_id'] = frag_id
            self.table.row['label'] = label
            self.table.row.append()
        self.table.flush()

    def __eq__(self, other):
        # same length
        if len(self.table) != len(other.table):
            return False

        # same columns
        if self.table.coldescrs != other.table.coldescrs:
            return False

        # same content
        mydict = self.label2ids()
        otherdict = other.label2ids()
        return mydict == otherdict

    def merge(self, label2id):
        """Merge label2id dict into self

        When label does not exists an id is generated and the label/id is added.
        When label does exist the id of the label in self is kept.

        Args:
            label2id (dict]): Dictionary with fragment label as key and fragment identifier as value.

        Returns:
            dict: Dictionary of label/id which where in label2id, but missing in self
        """
        id_offset = max([r['frag_id'] for r in self.table]) + 1

        mylabels = frozenset([r['label'].decode() for r in self.table])
        missing_label2ids = {}
        for row in label2id:
            if row['label'] not in mylabels:
                missing_label2ids[row['label']] = id_offset
                id_offset += 1

        self.update(missing_label2ids)
        return missing_label2ids

    def __iter__(self):
        for r in self.table.__iter__():
            yield {'frag_id': r['frag_id'], 'label': r['label'].decode()}

    def _copy(self, other, condition):
        hit = other.table.row
        for row in self.table:
            if condition(row['frag_id']):
                hit['frag_id'] = row['frag_id']
                hit['label'] = row['label']
                hit.append()
        other.table.flush()

    def keep(self, other, keep):
        """Copy content of self to other and only keep given fragment identifiers

        Args:
            other (LabelsLookup): Labels table to fill
            keep (set[int]): Fragment identifiers to keep
        """
        self._copy(other, lambda d: d in keep)

    def skip(self, other, skip):
        """Copy content of self to other and skip given fragment identifiers

        Args:
            other (LabelsLookup): Labels table to fill
            skip (set[int]): Fragment identifiers to skip
        """
        self._copy(other, lambda d: d not in skip)
