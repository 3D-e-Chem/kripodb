from nose.tools import assert_almost_equal
from intbitset import intbitset
import modifiedtanimoto.algorithm as algorithm


def test_mean_onbit_density():
    bitsets = [
        intbitset([1, 2, 3]),
        intbitset([1, 2, 4, 5, 8]),
        intbitset([1, 2, 4, 8])
    ]
    fp_size = 8

    result = algorithm.mean_onbit_density(bitsets, fp_size)

    expected = 0.5
    assert result == expected


def test_corrections():
    corr_st, corr_sto = algorithm.corrections(0.01)

    assert_almost_equal(corr_st, 0.663333333333)
    assert_almost_equal(corr_sto, 0.336666666667)


def test_distance():
    bitset1 = intbitset([1, 2, 3])
    bitset2 = intbitset([1, 2, 4, 8])
    fp_size = 8
    corr_st = 0.663333333333
    corr_sto = 0.336666666667

    result = algorithm.distance(bitset1, bitset2, fp_size, corr_st, corr_sto)

    expected = 0.20922222222203335
    assert_almost_equal(result, expected)