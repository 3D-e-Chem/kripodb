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
from math import log10, ceil

import tables


class DistanceMatrix(object):
    filters = tables.Filters(complevel=6, complib='blosc')

    def __init__(self, filename, mode='r'):
        self.h5file = tables.open_file(filename, mode, filters=self.filters)

    def pairs(self, expectedrows=0, precision=None):
        return PairsTable(self.h5file, expectedrows, precision)

    def labels(self, expectedrows=0):
        return LabelsLookup(self.h5file, expectedrows)

    def close(self):
        self.h5file.close()


class AbstractSimpleTable(object):
    def append(self, other):
        col_names = [col_name for col_name in self.table.colpathnames]
        nrows = 0
        for srcRow in other.table.iterrows():
            for col_name in col_names:
                self.table.row[col_name] = srcRow[col_name]
            self.table.row.append()
            nrows += 1
        self.table.flush()
        return nrows

    def __len__(self):
        return len(self.table)


class DistancePair(tables.IsDescription):
    a = tables.UInt32Col()
    b = tables.UInt32Col()
    score = tables.UInt16Col()


class PairsTable(AbstractSimpleTable):
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
        return self.table.attrs['score_precision']

    @score_precision.setter
    def score_precision(self, value):
        self.table.attrs['score_precision'] = value

    def add_indexes(self):
        self.table.cols.a.create_index(filters=self.filters)
        self.table.cols.b.create_index(filters=self.filters)

    def update(self, distances_iter, label2id):
        hit = self.table.row
        for label1, label2, distance in distances_iter:
            hit['a'] = label2id[label1]
            hit['b'] = label2id[label2]
            hit['score'] = int(distance * self.score_precision)
            hit.append()
        self.table.flush()

    def find(self, frag_id, cutoff):
        precision = float(self.score_precision)
        ndigits = int(ceil(log10(precision)))
        scutoff = int(cutoff * precision)

        hits = {}
        query1 = '(a == {}) & (score >= {})'.format(frag_id, scutoff)
        for row in self.table.where(query1):
            hit_id, score = row[1], row[2]
            score = round(score / precision, ndigits)
            hits[hit_id] = score

        query2 = '(b == {}) & (score >= {})'.format(frag_id, scutoff)
        for row in self.table.where(query2):
            query_id, score = row[0], row[2]
            score = round(score / precision, ndigits)
            hits[hit_id] = score

        return hits


class Id2Label(tables.IsDescription):
    frag_id = tables.UInt32Col()
    label = tables.StringCol(16)


class LabelsLookup(AbstractSimpleTable):
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
        return self.table.where('frag_id == {}'.format(frag_id)).next()[1]

    def by_label(self, label):
        return self.table.where('label == "{}"'.format(label)).next()[0]

    def update(self, label2id):
        for label, frag_id in label2id.iteritems():
            self.table.row['frag_id'] = frag_id
            self.table.row['label'] = label
            self.table.row.append()
        self.table.flush()


