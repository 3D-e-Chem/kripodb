import StringIO
from intbitset import intbitset
from nose.tools import assert_raises
from modifiedtanimoto import makebits as makebits


def test_read_header():
    line = 'MAKEBITS 1.0 574331 BigGrid'

    result = makebits.read_header(line)

    expected = ('MAKEBITS', '1.0', 574331, 'BigGrid')
    assert result == expected


def test_read_bitset():
    line = '3frb_TOP_frag24 1 2 3 4 6 10 11 12 15 0 9'

    (fid, bitset) = makebits.read_bitset(line, 100)

    expected = intbitset([1, 2, 3, 4, 6, 10, 11, 12, 15])
    assert fid == '3frb_TOP_frag24'
    assert bitset == expected


def test_read_bitset_toolong():
    line = '3frb_TOP_frag24 1 2 3 4 6 10 11 12 15 0 5'
    with assert_raises(Exception) as e:
        makebits.read_bitset(line, 100)
    expected = 'On bit checksum incorrect for 3frb_TOP_frag24'
    assert e.exception.message == expected


def test_read_file():
    infile = StringIO.StringIO('MAKEBITS 1.0 574331 BigGrid\n3frb_TOP_frag24 1 2 3 4 6 10 11 12 15 0 9\n')

    (bitsets, fp_size) = makebits.read_file(infile)

    expected = {'3frb_TOP_frag24': intbitset([1, 2, 3, 4, 6, 10, 11, 12, 15])}
    assert fp_size == 574331
    assert bitsets == expected
