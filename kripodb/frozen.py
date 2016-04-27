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
"""Distance matrix using pytables carray"""
from __future__ import absolute_import
from math import log10, ceil, floor
import numpy as np
from tables import parameters
import tables


class FrozenDistanceMatrix(object):
    filters = tables.Filters(complevel=6, complib='blosc')

    def __init__(self, filename, mode='r'):
        self.h5file = tables.open_file(filename, mode, filters=self.filters)
        self.score_precision = 2**16-1
        if 'labels' in self.h5file.root:
            self.labels = self.h5file.root.labels
        else:
            self.labels = None
        if 'scores' in self.h5file.root:
            self.scores = self.h5file.root.scores
        else:
            self.scores = None
        self.cache_i2l = {}
        self.cache_l2i = {}
        if self.labels is not None:
            self.build_label_cache()

    def close(self):
        """Closes the hdf5file"""
        self.h5file.close()

    def find(self, query, cutoff, limit=None):
        """Find similar fragments to query.

        Args:
            query (str): Query fragment identifier
            cutoff (float): Cutoff, distance scores below cutoff are discarded.
            limit (int): Maximum number of hits. Default is None for no limit.

        Returns:
            Tuple[(str, float)]: Hit fragment idenfier and distance score
        """
        precision = float(self.score_precision)
        precision10 = float(10**(floor(log10(precision))))
        scutoff = int(cutoff * precision)
        query_id = self.cache_l2i[query]
        subjects = self.table.root.scores[query_id, ...]
        filled_subjects_ids = subjects.nonzero()[0]
        filled_subjects = [(i, subjects[i]) for i in filled_subjects_ids]
        hits = [(self.cache_i2l[k], ceil(precision10 * v / precision) / precision10) for k, v in filled_subjects if v >= scutoff]
        sorted_hits = sorted(hits,  key=lambda r: r[1], reverse=True)
        if limit is not None:
            sorted_hits = sorted_hits[:limit]
        return sorted_hits

    def build_label_cache(self):
        self.cache_i2l = dict(enumerate(self.labels))
        self.cache_l2i = {v: k for k, v in self.cache_i2l.items()}

    def from_pairs(self, distance_matrix, frame_size, memory, limit=None):
        parameters.CHUNK_CACHE_SIZE = memory * 1024 ** 3
        parameters.CHUNK_CACHE_NELMTS = 2 ** 14

        nr_frags = len(distance_matrix.labels)

        print('Filling labels...', end='')

        id2labels = {v: k for k, v in distance_matrix.labels.label2ids().items()}
        id2nid = {v: k for k, v in enumerate(id2labels)}
        labels2nid = [None] * nr_frags
        for myid in id2nid:
            labels2nid[id2nid[myid]] = np.string_(id2labels[myid])
        self.labels = self.h5file.create_carray('/', 'labels', obj=labels2nid, filters=self.filters)
        self.h5file.flush()
        self.build_label_cache()

        print('Done')
        print('Filling matrix')
        self.scores = self.h5file.create_carray('/', 'scores', atom=tables.UInt16Atom(),
                                                shape=(nr_frags, nr_frags), chunkshape=(1, nr_frags),
                                                filters=self.filters)
        if limit is None:
            limit = len(distance_matrix.pairs)

        self._ingest_pairs(distance_matrix.pairs.table, id2nid, frame_size, limit)
        self.h5file.flush()

    def _ingest_pairs(self, pairs, id2nid, frame_size, limit):
        nr_cells = 0
        for start in range(0, limit, frame_size):
            stop = frame_size + start
            raw_frame = pairs.read(start=start, stop=stop)
            frame = self._translate_frame(raw_frame, id2nid)
            nr_cells += frame.shape[0] * 2
            self._ingest_pairs_frame(frame)
            print('Filled {0} cells in matrix from {1}:{2}'.format(nr_cells, start, stop))
            del frame

    def _translate_frame(self, frame, id2nid):
        framet = []
        for pair in frame:
            framet.append((id2nid[pair[0]], id2nid[pair[1]], pair[2]))
        framenp = np.array(framet, dtype=[('a', '<u4'), ('b', '<u4'), ('score', '<u2')])
        return framenp

    def _ingest_pairs_frame(self, frame):
        print('Sorting a>b')
        frame.sort(order=('a', 'b'))
        print('Inserting a>b')
        scores = self.scores
        for row in frame:
            scores[row[0], row[1]] = row[2]
        print('Sorting b>a')
        # sort on b, to write to same chunk longer
        frame.sort(order=('b', 'a'))
        print('Inserting b>a')
        # loop in reverse order, to have cache hits at start
        for row in frame[::-1]:
            scores[row[1], row[0]] = row[2]
