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
"""Distance matrix using hdf5 as storage backend."""
from __future__ import absolute_import
from math import log10, ceil, floor

import tables
import six


class DistanceMatrix(object):
    """Distance matrix

    Args:
        filename (str): File name of hdf5 file to write or read distance matrix from
        mode (str): Can be 'r' for reading or 'w' for writing
        expectedpairrows (int): Expected number of pairs to be added.
            Required when distance matrix is opened in write mode, helps optimize storage
        precision (int): Distance score is a fraction,
            the score is converted to an int by multiplying it with the precision
        expectedlabelrows (int): Expected number of labels to be added.
            Required when distance matrix is opened in write mode, helps optimize storage
        cache_labels (bool): Cache labels, speed up label lookups

    Attributes:
        h5file (tables.File): Object representing an open hdf5 file
        pairs (PairsTable): HDF5 Table that contains pairs
        labels (LabelsLookup): Table to look up label of fragment by id or id of fragment by label
    """
    filters = tables.Filters(complevel=6, complib='blosc')

    def __init__(self, filename, mode='r', expectedpairrows=None, precision=None, expectedlabelrows=None, cache_labels=False):
        self.h5file = tables.open_file(filename, mode, filters=self.filters)
        self.pairs = PairsTable(self.h5file, expectedpairrows, precision)
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
        """Append data from other distance matrix to me

        Args:
            other (DistanceMatrix): Other distance matrix
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

    def update(self, distances_iter, label2id):
        """Store pairs of fragment identifier with their distance score and label 2 id lookup

        Args:
            distances_iter (Iterator): Iterator which yields (label1, label2, distance_score)
            label2id (Dict): Dictionary with fragment label as key and fragment identifier as value.

        """
        self.pairs.update(distances_iter, label2id)
        self.pairs.add_indexes()
        self.labels.update(label2id)

    def find(self, query, cutoff, limit=None):
        """Find similar fragments to query.

        Args:
            query (str): Query fragment identifier
            cutoff (float): Cutoff, distance scores below cutoff are discarded.
            limit (int): Maximum number of hits. Default is None for no limit.

        Yields:
            Tuple[(str, float)]: Hit fragment idenfier and distance score
        """
        if self.cache_l2i:
            frag_id = self.cache_l2i[query]
            for hit_frag_id, score in self.pairs.find(frag_id, cutoff, limit):
                yield self.cache_i2l[hit_frag_id], score
        else:
            frag_id = self.labels.by_label(query)
            for hit_frag_id, score in self.pairs.find(frag_id, cutoff, limit):
                yield self.labels.by_id(hit_frag_id), score


class AbstractSimpleTable(object):
    """Abstract wrapper around a HDF5 table

    Args:
        table (tables.Table): HDF5 table

    Attributes
        table (tables.Table): HDF5 table

    """

    def __init__(self, table):
        self.table = table

    def append(self, other):
        """Append rows of other table to self

        Args:
            other: Table of same type as self

        """
        self.table.append(other.table.read())

    def __len__(self):
        """

        Returns:
            int: Number of rows in table
        """
        return len(self.table)

    def __iter__(self):
        return self.table.__iter__()


class DistancePair(tables.IsDescription):
    """Table description for distance pair"""
    a = tables.UInt32Col()
    b = tables.UInt32Col()
    score = tables.UInt16Col()


class PairsTable(AbstractSimpleTable):
    """Tabel to store distance score of a pair of fragment fingerprints

    When table does not exist in h5file it is created.

    Args:
        h5file (tables.File): Object representing an open hdf5 file
        expectedrows (int): Expected number of pairs to be added.
            Required when distance matrix is opened in write mode, helps optimize storage
        precision (int): Distance score is a fraction,
                the score is converted to an int by multiplying it with the precision

    Attributes:
        score_precision (int): Distance score is a fraction,
            the score is converted to an int by multiplying it with the precision
        full_matrix (bool): Matrix is filled above and below diagonal.
    """
    table_name = 'pairs'
    filters = tables.Filters(complevel=6, complib='blosc')

    def __init__(self, h5file, expectedrows=0, precision=None):
        if self.table_name in h5file.root:
            table = h5file.root.__getattr__(self.table_name)
        else:
            table = h5file.create_table('/',
                                        self.table_name,
                                        DistancePair,
                                        'Distance pairs',
                                        expectedrows=expectedrows)

        self.table = table
        if precision is not None:
            self.score_precision = precision

    @property
    def score_precision(self):
        try:
            return self.table.attrs['score_precision']
        except KeyError:
            return None

    @score_precision.setter
    def score_precision(self, value):
        self.table.attrs['score_precision'] = value

    @property
    def full_matrix(self):
        try:
            return self.table.attrs['full_matrix']
        except KeyError:
            return False

    @full_matrix.setter
    def full_matrix(self, value):
        self.table.attrs['full_matrix'] = value

    def add_indexes(self):
        """Add indexes on identifier columns

        Best done after calling update().
        """
        self.table.cols.a.create_index(filters=self.filters)
        if not self.full_matrix:
            self.table.cols.b.create_index(filters=self.filters)

    def update(self, distances_iter, label2id):
        """Store pairs of fragment identifier with their distance score

        Args:
            distances_iter (Iterator): Iterator which yields (label1, label2, distance_score)
            label2id (Dict): Lookup with fragment label as key and fragment identifier as value

        """
        hit = self.table.row
        for label1, label2, distance in distances_iter:
            hit['a'] = label2id[label1]
            hit['b'] = label2id[label2]
            hit['score'] = int(distance * self.score_precision)
            hit.append()
        self.table.flush()

    def find(self, frag_id, cutoff, limit):
        """Find fragment hits which has a distance score with frag_id above cutoff.

        Args:
            frag_id (int): query fragment identifier
            cutoff (float): Cutoff, distance scores below cutoff are discarded.
            limit (int): Maximum number of hits. Default is None for no limit.

        Returns:
            List[Tuple]: Where first tuple value is hit fragment identifier and second value is distance score

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
            Required when distance matrix is opened in write mode, helps optimize storage
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

        self.table = table

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

    def label2ids(self):
        """Return whole table as a dictionary

        Returns:
            Dict: Dictionary with label as key and frag_id as value.

        """
        return {r['label'].decode(): r['frag_id'] for r in self.table}

    def update(self, label2id):
        """Update labels lookup by adding labels in label2id.

        Args:
            label2id (Dict): Dictionary with fragment label as key and fragment identifier as value.

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
            label2id (Dict): Dictionary with fragment label as key and fragment identifier as value.

        Returns:
            Dict: Dictionary of label/id which where in label2id, but missing in self
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
