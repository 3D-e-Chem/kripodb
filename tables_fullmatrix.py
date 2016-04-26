"""

Create h5py node for each query with arrays

"""
import numpy as np
from kripodb.hdf5 import DistanceMatrix
from tables import parameters
import tables

parameters.CHUNK_CACHE_SIZE = 8 * 1024 ** 3
parameters.CHUNK_CACHE_NELMTS = 2 ** 14

dm = DistanceMatrix('/media/sf_data/3d-e-chem/feb2016/01-01_to_13-13.out.h5')
nr_frags = len(dm.labels)

f = tables.open_file("mytestfile.h5", mode="w")

id2labels = {v: k for k, v in dm.labels.label2ids().items()}
id2nid = {v: k for k, v in enumerate(id2labels)}

print('Labels lookup')

labels2nid = [None] * nr_frags
for myid in id2nid:
    labels2nid[id2nid[myid]] = np.string_(id2labels[myid])
f.create_carray('/', 'labels', obj=labels2nid, filters=dm.filters)
f.flush()

print('Filling')


def fill_matrix(frame, scores, id2nid, id2labels):
    i = 0
    last_query = None
    s2q = {}
    for row in frame:
        query = row['a']
        querynid = id2nid[query]
        subject = row['b']
        subjectnid = id2nid[subject]
        score = row['score']
        scores[querynid, subjectnid] = score
        if last_query != query and last_query:
            last_query_nid = id2nid[last_query]
            print('Writing complementary {0} pairs of {1}: {2}  ... '.format(len(s2q),
                                                                             last_query_nid,
                                                                             id2labels[last_query]), end='')
            for k, v in s2q.items():
                scores[k, last_query_nid] = v
                i += 1

            print('Done')
            s2q = {}
        i += 1
        s2q[subjectnid] = score
        last_query = query

    last_query_nid = id2nid[last_query]
    print('Writing complementary {0} pairs of {1}: {2}  ... '.format(len(s2q), last_query_nid, id2labels[last_query]),
          end='')
    for k, v in s2q.items():
        scores[k, last_query_nid] = v
        i += 1
    print('Done')
    return i

scores = f.create_carray('/', 'scores', atom=tables.UInt16Atom(),
                         shape=(nr_frags, nr_frags), chunkshape=(1, nr_frags),
                         filters=dm.filters)
nr_cells = 0
frame_size = 250000
end = len(dm.pairs)
end = 1000000
for start in range(0, end, frame_size):
    frame = dm.pairs.table.read(start=start, stop=frame_size + start)
    nr_cells += fill_matrix(frame, scores, id2nid, id2labels)
    print('Filled {0} cells in matrix of {1}x{1}'.format(nr_cells, nr_frags))

f.flush()
f.close()

dm.close()

# stop=10000
# No cache - 21.4s
# Cache 1G - 20.6s
# Cache 4G - 12.8s
# Single flush 4G - 9.5s
# stop=100000
# 4G multiflush - 3m0.852s
# 4G - 2m24.9s
# stop=200000
# 4G - 5m48s
# tables
# 9.1s
# stop=100000 no cache- 1m12s
# stop=100000 4G - 53s
# stop=100000 2G - 1m1.6s
# stop=100000 5G - 48s
# stop=100000 6G - 49s
# stop=200000 6G - 1m28s - 108M
# stop=300000 6G - 2m11.8s - 116M
# stop=400000 6G - 2m49s - 123M
# func - 49s
# jit func - 55s

# Query pattern
# i2l = dict(enumerate(f.root.labels))
# l2i = {v: k for k, v in i2l.items()}
# ql = b'1wnt_NAP_frag4'
# qi = l2i[ql]
# rh = f.root.scores[qi, ...]
# p = rh.nonzero()[0]
# z = {i2l[i]: rh[i] for i in p}
# h = sorted([(k, v) for k, v in z.items() if v > 40000], key=lambda r: r[1], reverse=True)[:10]
#
