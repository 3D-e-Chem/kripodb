import os

from .utils import tmpname
from kripodb.pharmacophores import PharmacophoresDb


class PharmacophoresDbInMemory(object):
    def __init__(self):
        self.db_fn = tmpname()
        self.db = PharmacophoresDb(self.db_fn, 'a', driver='H5FD_CORE', driver_core_backing_store=0)

    def __enter__(self):
        return self.db

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.db.close()
        if os.path.isfile(self.db_fn):
            os.remove(self.db_fn)


class TestPharmacophoresDb(object):
    def test_len_empty(self):
        with PharmacophoresDbInMemory() as db:
            assert len(db.points) == 0

    def test_add_fragment(self):
        with PharmacophoresDbInMemory() as db:
            db.points.add_fragment('frag1', [0], [('POSC', 3.4, 5.6, 7.8)])
            db.points.table.flush()
            assert len(db.points) == 1
