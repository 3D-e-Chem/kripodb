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
from __future__ import absolute_import, print_function
from math import log10, ceil, floor
try:
    from time import process_time
except ImportError:
    from time import clock as process_time

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import six
import tables


class FrozenDistanceMatrix(object):
    """Frozen distances matrix

    Compression is used so None/0 take no space.

    Warning! Can not be enlarged.

    Args:
        filename (str): File name of hdf5 file to write or read distance matrix from
        mode (str): Can be 'r' for reading or 'w' for writing
        **kwargs: Passed though to tables.open_file()

    Attributes:
        h5file (tables.File): Object representing an open hdf5 file
        scores (tables.CArray): HDF5 Table that contains matrix
        labels (tables.CArray): Table to look up label of fragment by id or id of fragment by label

    """
    filters = tables.Filters(complevel=6, complib='blosc')

    def __init__(self, filename, mode='r', **kwargs):
        self.h5file = tables.open_file(filename, mode, filters=self.filters, **kwargs)
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
        subjects = self.h5file.root.scores[query_id, ...]
        filled_subjects_ids = subjects.nonzero()[0]
        filled_subjects = [(i, subjects[i]) for i in filled_subjects_ids]
        hits = [(self.cache_i2l[k], ceil(precision10 * v / precision) / precision10) for k, v in filled_subjects if v >= scutoff]
        sorted_hits = sorted(hits,  key=lambda r: r[1], reverse=True)
        if limit is not None:
            sorted_hits = sorted_hits[:limit]
        return sorted_hits

    def build_label_cache(self):
        self.cache_i2l = {k: v.decode() for k, v in enumerate(self.labels)}
        self.cache_l2i = {v: k for k, v in self.cache_i2l.items()}

    def from_pairs(self, distance_matrix, frame_size, limit=None, single_sided=False):
        """Fills self with matrix which is stored in pairs.

        Also known as COOrdinate format, the 'ijv' or 'triplet' format.

        Args:
            distance_matrix (kripodb.hdf5.DistanceMatrix):
            frame_size (int): Number of pairs to append in a single go
            limit (int|None): Number of pairs to add, None for no limit, default is None.
            single_sided (bool): If false add stored direction and reverse direction. Default is False.


        time kripodb distances freeze --limit 200000 -f 100000 data/feb2016/01-01_to_13-13.out.h5 percell.h5
        47.2s
        time kripodb distances freeze --limit 200000 -f 100000 data/feb2016/01-01_to_13-13.out.h5 coo.h5
        0.2m - 2m6s
        .4m - 2m19s
        .8m - 2m33s
        1.6m - 2m48s
        3.2m - 3m4s
        6.4m  - 3m50s
        12.8m - 4m59s
        25.6m - 7m27s
        """
        nr_frags = len(distance_matrix.labels)

        six.print_('Filling labels ... ', end='')

        id2labels = {v: k for k, v in distance_matrix.labels.label2ids().items()}
        id2nid = {v: k for k, v in enumerate(id2labels)}
        labels2nid = [None] * nr_frags
        for myid in id2nid:
            labels2nid[id2nid[myid]] = np.string_(id2labels[myid])
        self.labels = self.h5file.create_carray('/', 'labels', obj=labels2nid, filters=self.filters)
        self.h5file.flush()
        self.build_label_cache()

        six.print_('Done')
        six.print_('Filling matrix')
        self.scores = self.h5file.create_carray('/', 'scores', atom=tables.UInt16Atom(),
                                                shape=(nr_frags, nr_frags), chunkshape=(1, nr_frags),
                                                filters=self.filters)
        if limit is None:
            limit = len(distance_matrix.pairs)

        self._ingest_pairs(distance_matrix.pairs.table, id2nid, frame_size, limit, single_sided)
        self.h5file.flush()

    def _ingest_pairs(self, pairs, id2nid, frame_size, limit, single_sided):
        i = 0
        for start in range(0, limit, frame_size):
            stop = frame_size + start
            t1 = process_time()
            six.print_('Fetching pairs {0}:{1} of {2} ... '.format(start, stop, limit), end='', flush=True)
            raw_frame = pairs.read(start=start, stop=stop)
            t2 = process_time()
            six.print_('{0}s, Parsing '.format(int(t2 - t1)), end='', flush=True)
            frame = self._translate_frame(raw_frame, id2nid, single_sided)
            t3 = process_time()
            six.print_(' {0}s, Writing ... '.format(int(t3 - t2)), end='', flush=True)
            # alternate direction, to make use of cached chunks of prev frame
            direction = 1
            if i % 2 == 0:
                direction = -1
            self._ingest_pairs_frame(frame, direction)
            del frame
            t4 = process_time()
            six.print_('{0}s, Done'.format(int(t4 - t3)), flush=True)
            i += 1

    def _translate_frame(self, raw_frame, id2nid, single_sided):
        frame = np.array([], dtype=[('a', '<u4'), ('b', '<u4'), ('score', '<u2')])
        six.print_('.', end='', flush=True)
        if single_sided:
            frame.resize((len(raw_frame),))
        else:
            frame.resize((len(raw_frame) * 2,))
        i = 0
        six.print_('.', end='', flush=True)
        for pair in raw_frame:
            a = id2nid[pair[0]]
            b = id2nid[pair[1]]
            frame[i] = (a, b, pair[2])
            i += 1
        six.print_('.', end='', flush=True)
        a = frame['a']
        b = frame['b']
        data = frame['score']
        nr_frags = len(id2nid)
        smat = coo_matrix((data, (b, a)), shape=(nr_frags, nr_frags)).tocsc()
        if not single_sided:
            smat += smat.transpose()
        return smat

    def _ingest_pairs_frame(self, frame, direction):
        scores = self.scores
        for current_col_idx in six.moves.range(frame.shape[0]):
            new_col = frame.getcol(current_col_idx)
            if not new_col.nnz:
                # new col has only zeros, skipping
                continue
            current_col = scores[current_col_idx, ...]
            # write whole column, so chunk compression + shuffle only performed once
            scores[current_col_idx, ...] = current_col + new_col.toarray()[:, 0]
            # TODO use [x,y] = v when nnz is low fraction of nr_frags

    def to_pandas(self):
        """Pandas dataframe with labelled colums and rows.

        Warning! Only use on matrices that fit in memory

        Returns:
            pd.DataFrame

        """
        precision = float(self.score_precision)
        decimals = int(log10(precision))
        labels = [v.decode() for v in self.labels]
        df = pd.DataFrame(self.scores.read(), index=labels, columns=labels)
        df /= precision
        df = df.round(decimals)
        return df

