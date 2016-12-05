import os
import tempfile

from kripodb.hdf5 import SimilarityMatrix
from kripodb.frozen import FrozenSimilarityMatrix


def tmpname():
    tmpf = tempfile.NamedTemporaryFile()
    out_file = tmpf.name
    tmpf.close()
    return out_file


class SimilarityMatrixInMemory(object):
    def __init__(self):
        self.matrix_fn = tmpname()
        self.matrix = SimilarityMatrix(self.matrix_fn, 'a', driver='H5FD_CORE', driver_core_backing_store=0)

    def __enter__(self):
        return self.matrix

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.matrix.close()
        if os.path.isfile(self.matrix_fn):
            os.remove(self.matrix_fn)


class FrozenSimilarityMatrixInMemory(object):
    def __init__(self):
        self.matrix_fn = tmpname()
        self.matrix = FrozenSimilarityMatrix(self.matrix_fn, 'a', driver='H5FD_CORE', driver_core_backing_store=0)

    def __enter__(self):
        return self.matrix

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.matrix.close()
        if os.path.isfile(self.matrix_fn):
            os.remove(self.matrix_fn)
