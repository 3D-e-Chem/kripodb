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
"""Similarity matrix using pytables carray"""
from __future__ import absolute_import, print_function
from math import log10, ceil, floor
try:
    # for Python >3.3
    from time import process_time
except ImportError:
    from time import clock as process_time
import numpy as np
import pandas as pd
from progressbar import ProgressBar
from scipy.sparse import coo_matrix
import six
import tables


class FrozenSimilarityMatrix(object):
    """Frozen similarities matrix

    Can retrieve whole column of a specific row fairly quickly.
    Store as compressed dense matrix.
    Due to compression the zeros use up little space.

    Warning! Can not be enlarged.

    Compared find performance FrozenSimilarityMatrix with SimilarityMatrix::

    >>> from kripodb.db import FragmentsDb
    >>> db = FragmentsDb('data/feb2016/Kripo20151223.sqlite')
    >>> ids = [v[0] for v in db.cursor.execute('SELECT frag_id FROM fragments ORDER BY RANDOM() LIMIT 20')]
    >>> from kripodb.frozen import FrozenSimilarityMatrix
    >>> fdm = FrozenSimilarityMatrix('01-01_to_13-13.out.frozen.blosczlib.h5')
    >>> from kripodb.hdf5 import SimilarityMatrix
    >>> dm = SimilarityMatrix('data/feb2016/01-01_to_13-13.out.h5', cache_labels=True)
    >>> %timeit list(dm.find(ids[0], 0.45, None))
    ... 1 loop, best of 3: 1.96 s per loop
    >>>  %timeit list(fdm.find(ids[0], 0.45, None))
    ... The slowest run took 6.21 times longer than the fastest. This could mean that an intermediate result is being cached.
    ... 10 loops, best of 3: 19.3 ms per loop
    >>> ids = [v[0] for v in db.cursor.execute('SELECT frag_id FROM fragments ORDER BY RANDOM() LIMIT 20')]
    >>> %timeit -n1 [list(fdm.find(v, 0.45, None)) for v in ids]
    ... 1 loop, best of 3: 677 ms per loop
    >>> %timeit -n1 [list(dm.find(v, 0.45, None)) for v in ids]
    ... 1 loop, best of 3: 29.7 s per loop

    Args:
        filename (str): File name of hdf5 file to write or read similarity matrix from
        mode (str): Can be 'r' for reading or 'w' for writing
        **kwargs: Passed though to tables.open_file()

    Attributes:
        h5file (tables.File): Object representing an open hdf5 file
        scores (tables.CArray): HDF5 Table that contains matrix
        labels (tables.CArray): Table to look up label of fragment by id or id of fragment by label

    """
    filters = tables.Filters(complevel=6, complib='blosc', shuffle=True)

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
            cutoff (float): Cutoff, similarity scores below cutoff are discarded.
            limit (int): Maximum number of hits. Default is None for no limit.

        Returns:
            list[tuple[str,float]]: Hit fragment identifier and similarity score
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

    def __getitem__(self, item):
        """Get all similarities of fragment.

        Self is excluded.

        Args:
            item (STR): Label of a fragment

        Returns:
            list[tuple[str, float]]: list of (fragment_label, score)

        """
        precision = float(self.score_precision)
        precision10 = float(10**(floor(log10(precision))))
        query_id = self.cache_l2i[item]
        subjects = self.h5file.root.scores[query_id, ...]
        hits = [(self.cache_i2l[k], ceil(precision10 * v / precision) / precision10) for k, v in enumerate(subjects) if k != query_id]
        return hits

    def build_label_cache(self):
        self.cache_i2l = {k: v.decode() for k, v in enumerate(self.labels)}
        self.cache_l2i = {v: k for k, v in self.cache_i2l.items()}

    def from_pairs(self, similarity_matrix, frame_size, limit=None, single_sided=False):
        """Fills self with matrix which is stored in pairs.

        Also known as COOrdinate format, the 'ijv' or 'triplet' format.

        Args:
            similarity_matrix (kripodb.hdf5.SimilarityMatrix):
            frame_size (int): Number of pairs to append in a single go
            limit (int|None): Number of pairs to add, None for no limit, default is None.
            single_sided (bool): If false add stored direction and reverse direction. Default is False.


        time kripodb similarities freeze --limit 200000 -f 100000 data/feb2016/01-01_to_13-13.out.h5 percell.h5
        47.2s
        time kripodb similarities freeze --limit 200000 -f 100000 data/feb2016/01-01_to_13-13.out.h5 coo.h5
        0.2m - 2m6s
        .4m - 2m19s
        .8m - 2m33s
        1.6m - 2m48s
        3.2m - 3m4s
        6.4m  - 3m50s
        12.8m - 4m59s
        25.6m - 7m27s
        """
        nr_frags = len(similarity_matrix.labels)

        six.print_('Filling labels ... ', end='')

        id2labels = {v: k for k, v in similarity_matrix.labels.label2ids().items()}
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
            limit = len(similarity_matrix.pairs)

        self._ingest_pairs(similarity_matrix.pairs.table, id2nid, frame_size, limit, single_sided)
        self.h5file.flush()

    def _ingest_pairs(self, pairs, oid2nid, frame_size, limit, single_sided):
        oid2nid_v = np.vectorize(oid2nid.get)
        # whole pairs set does not fit in memory, so split it in frames with `frame_size` number of pairs.
        for start in range(0, limit, frame_size):
            stop = frame_size + start
            t1 = process_time()
            six.print_('Fetching pairs {0}:{1} of {2} ... '.format(start, stop, limit), end='', flush=True)
            raw_frame = pairs.read(start=start, stop=stop)
            t2 = process_time()
            six.print_('{0}s, Parsing ... '.format(int(t2 - t1)), flush=True)
            frame = self._translate_frame(raw_frame, oid2nid_v, single_sided)
            t3 = process_time()
            six.print_('Writing ... '.format(int(t3 - t2)), flush=True)
            # alternate direction, to make use of cached chunks of prev frame
            self._ingest_pairs_frame(frame)
            del frame
            t4 = process_time()
            six.print_('{0}s, Done with {1}:{2} in {3}s'.format(int(t4 - t3), start, stop, int(t4 - t1)), flush=True)

    def _translate_frame(self, raw_frame, oid2nid, single_sided):
        bar = ProgressBar(max_value=4)
        bar.update(0)
        a = oid2nid(raw_frame['a'])
        bar.update(1)
        b = oid2nid(raw_frame['b'])
        bar.update(2)
        data = raw_frame['score']
        nr_frags = len(self.labels)
        smat = coo_matrix((data, (b, a)), shape=(nr_frags, nr_frags)).tocsc()
        bar.update(3)
        if not single_sided:
            smat += smat.transpose()
        bar.update(4)
        return smat

    def _ingest_pairs_frame(self, frame):
        scores = self.scores
        # loop each
        bar = ProgressBar()
        for row_idx in bar(six.moves.range(frame.shape[0])):
            new_row = frame.getcol(row_idx)
            if not new_row.nnz:
                # new col has only zeros, skipping
                continue
            current_row = scores[row_idx, ...]
            # write whole column, so chunk compression + shuffle only performed once per row_idx
            scores[row_idx, ...] = current_row + new_row.toarray()[:, 0]

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

    def from_array(self, data, labels):
        """Fill matrix from 2 dimensional array

        Args:
            data (np.array): 2 dimensional square array with scores
            labels (list): List of labels for each column and row index
        """
        labels = [np.string_(d) for d in labels]
        self.labels = self.h5file.create_carray('/', 'labels', obj=labels, filters=self.filters)
        self.h5file.flush()
        self.build_label_cache()

        nr_frags = len(labels)
        self.scores = self.h5file.create_carray('/', 'scores', atom=tables.UInt16Atom(),
                                                shape=(nr_frags, nr_frags), chunkshape=(1, nr_frags),
                                                filters=self.filters)
        self.scores[0:nr_frags, 0:nr_frags] = (data * self.score_precision).astype('uint16')

    def to_pairs(self, pairs):
        """Copies labels and scores from self to pairs matrix.

        Args:
            pairs (SimilarityMatrix):

        """
        six.print_('copy labels', flush=True)
        self.build_label_cache()
        pairs.labels.update(self.cache_l2i)

        six.print_('copy matrix to pairs', flush=True)
        limit = self.scores.shape[0]
        bar = ProgressBar()
        for query_id in bar(six.moves.range(0, limit)):
            subjects = self.scores[query_id, ...]
            filled_subjects_ids = subjects.nonzero()[0]
            filled_subjects = [(query_id, i, subjects[i]) for i in filled_subjects_ids if query_id < i]
            if filled_subjects:
                pairs.pairs.table.append(filled_subjects)

    def count(self, frame_size=None, raw_score=False, lower_triangle=False):
        """Count occurrences of each score

        Only scores are counted of the upper triangle or lower triangle.
        Zero scores are skipped.

        Args:
            frame_size (int): Dummy argument to force same interface for thawed and frozen matrix
            raw_score (bool): When true return raw int16 score else fraction score
            lower_triangle (bool): When true return scores from lower triangle else return scores from upper triangle

        Returns:
            Tuple[(str, int)]: Score and number of occurrences
        """
        nr_rows = self.scores.shape[0]
        nr_bins = self.score_precision + 1
        counts = np.zeros(shape=nr_bins, dtype=np.int64)
        bar = ProgressBar()
        for query_id in bar(six.moves.range(0, nr_rows)):
            if lower_triangle:
                subjects = self.scores[query_id, query_id + 1:]
            else:
                subjects = self.scores[query_id, :query_id + 1]
            frame_counts = np.bincount(subjects[subjects.nonzero()], minlength=nr_bins)
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
