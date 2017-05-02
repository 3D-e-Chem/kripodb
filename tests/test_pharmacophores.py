import os
from six import BytesIO, StringIO

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

from .utils import tmpname
from kripodb.pharmacophores import PharmacophoresDb, read_pphore_sdfile, as_phar, read_fragtxtfile_as_file


@pytest.fixture
def example1_sdfile():
    body = b'''
 OpenBabel07051617103D

 14  0  0  0  0  0  0  0  0  0999 V2000
   22.5699   -6.3076   36.8593 O   0  6  0  0  0
   23.7871   -3.9004   36.3395 P   0  7  0  0  0
   23.7871   -3.9004   36.3395 P   0  7  0  0  0
   23.1923   -6.6223   36.0325 P   0  7  0  0  0
   23.1923   -6.6223   36.0325 P   0  7  0  0  0
   18.7201   -9.8937   40.4312 As  0  7  0  0  0
   18.3503   -9.5392   39.3836 Ne  0  0  0  0  0
   20.4171  -10.6185   40.3362 Ne  0  0  0  0  0
   19.7922   -6.5243   39.4320 He  0  0  0  0  0
   17.9406   -2.4043   34.9401 O   0  6  0  0  0
   14.6641   -7.4275   36.2138 As  0  7  0  0  0
   15.4420   -8.2931   36.1398 Ne  0  0  0  0  0
   14.4007   -6.8416   35.2404 Ne  0  0  0  0  0
   22.3608   -5.1679   39.2345 Rn  0  0  0  0  0
M  CHG  8   1  -2   2  -3   3  -3   4  -3   5  -3   6  -3  10  -2  11  -3
M  END
$$$$
'''
    return BytesIO(body)


@pytest.fixture
def example1_points():
    return [
        ('HACC', 22.5699, -6.3076, 36.8593),
        ('POSC', 23.7871, -3.9004, 36.3395),
        ('POSC', 23.7871, -3.9004, 36.3395),
        ('POSC', 23.1923, -6.6223, 36.0325),
        ('POSC', 23.1923, -6.6223, 36.0325),
        ('HDON', 18.7201, -9.8937, 40.4312),
        ('NEGC', 18.3503, -9.5392, 39.3836),
        ('NEGC', 20.4171, -10.6185, 40.3362),
        ('LIPO', 19.7922, -6.5243, 39.432),
        ('HACC', 17.9406, -2.4043, 34.9401),
        ('HDON', 14.6641, -7.4275, 36.2138),
        ('NEGC', 15.442, -8.2931, 36.1398),
        ('NEGC', 14.4007, -6.8416, 35.2404),
        ('AROM', 22.3608, -5.1679, 39.2345),
    ]


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


def assert_points(result, expected):
    npresult = np.array(result)
    npexpected = np.array(expected)
    assert_array_almost_equal(
        npresult[:, (1, 2, 3)].astype('float'),
        npexpected[:, (1, 2, 3)].astype('float'),
        4,
        'Positions not equal'
    )
    assert npresult[:, 0].tolist() == npexpected[:, 0].tolist()


def test_read_pphore_sdfile(example1_sdfile, example1_points):
    result = read_pphore_sdfile(example1_sdfile)

    assert_points(result, example1_points)


@pytest.fixture
def filled_PharmacophorePointsTable(example1_points):
    with PharmacophoresDbInMemory() as db:
        db.points.add_fragment('frag1', [0], example1_points)
        db.points.add_fragment('frag2', [1, 2], example1_points)
        db.points.add_fragment('frag3', [0, 2, 13], example1_points)
        db.points.table.flush()
        yield db.points


class TestPharmacophorePointsTable(object):
    def test_len_empty(self):
        with PharmacophoresDbInMemory() as db:
            assert len(db.points) == 0

    def test_len(self, filled_PharmacophorePointsTable):
        assert len(filled_PharmacophorePointsTable) == 6

    def test_getitem_present(self, filled_PharmacophorePointsTable, example1_points):
        result = filled_PharmacophorePointsTable['frag1']

        expected = [example1_points[0]]
        assert_points(result, expected)

    def test_getitem_absent(self, filled_PharmacophorePointsTable):
        with pytest.raises(KeyError) as excinfo:
            filled_PharmacophorePointsTable['frag999']
        assert excinfo.value.args[0] == 'frag999'

    def test_contains_present(self, filled_PharmacophorePointsTable, example1_points):
        assert 'frag1' in filled_PharmacophorePointsTable

    def test_contains_absent(self, filled_PharmacophorePointsTable):
        assert 'frag999' not in filled_PharmacophorePointsTable

    def test_add_fragment_duplicate(self, filled_PharmacophorePointsTable, example1_points):
        with pytest.raises(ValueError) as excinfo:
            filled_PharmacophorePointsTable.add_fragment('frag1', [0], example1_points)
        assert str(excinfo.value) == "Duplicate key 'frag1' found"
        assert len(filled_PharmacophorePointsTable) == 6


@pytest.fixture
def example1_phar():
    return '''frag1
HACC 22.5699 -6.3076 36.8593 0 0 0 0 0
$$$$
'''


@pytest.fixture
def example3_phar():
    return '''frag3
HACC 22.5699 -6.3076 36.8593 0 0 0 0 0
POSC 23.7871 -3.9004 36.3395 0 0 0 0 0
AROM 22.3608 -5.1679 39.2345 0 0 0 0 0
$$$$
'''


def test_as_phar(filled_PharmacophorePointsTable, example3_phar):
    frag_id = 'frag3'
    points = filled_PharmacophorePointsTable[frag_id]

    result = as_phar(frag_id, points)

    assert result == example3_phar


@pytest.fixture
def example4fragtxtfile():
    return StringIO('''frag1 1 2 3
frag2 
frag3 2 3
frag4 1
''')


def test_read_fragtxtfile_as_file(example4fragtxtfile):
    result = read_fragtxtfile_as_file(example4fragtxtfile)

    expected = {
        'frag1': [0, 1, 2],
        'frag2': [],
        'frag3': [1, 2],
        'frag4': [0]
    }
    assert result == expected
