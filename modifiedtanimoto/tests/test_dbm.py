from intbitset import intbitset
import modifiedtanimoto.dbm as dbm


def test_fastdump():
    bitsets = {
        'a': intbitset([2, 5, 8, 23])
    }
    result = {}
    dbm.fastdump(bitsets, result)

    expected = bitsets['a'].fastdump()
    assert result['a'] == expected


def test_fastload():
    bs = intbitset([2, 5, 8, 23])
    db = {'a': bs.fastdump()}

    bitsets = {}
    dbm.fastload(db, bitsets)

    assert bitsets['a'] == bs