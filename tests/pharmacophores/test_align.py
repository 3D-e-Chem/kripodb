# Copyright 2016 Netherlands eScience Center
#
# Licensed under the Apache License, Version 2.0 (the "License");
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
from __future__ import absolute_import

import pytest
from numpy.testing import assert_array_almost_equal

from kripodb.pharmacophores.align import Aligner, NoOverlapFound, align, distances


@pytest.fixture
def reference_pharmacophore():
    return [
        ('LIPO', 26.0369, 11.9800, 4.3352),
        ('LIPO', 25.4947, 12.4949, 2.8223),
        ('HDON', 19.5809, 17.5262, 6.2020),
        ('NEGC', 20.3908, 17.2075, 6.9782),
        ('NEGC', 19.9137, 17.6862, 5.0959),
        ('NEGC', 18.4865, 18.9091, 6.7297),
        ('NEGC', 19.0379, 25.9383, 15.2001),
        ('HDON', 23.7426, 26.2927, 14.6188),
        ('LIPO', 16.9648, 16.8184, 12.5531),
        ('LIPO', 16.1188, 17.5788, 13.7996),
        ('HACC', 19.2926, 17.1766, 15.7563),
        ('HDON', 21.8512, 17.5877, 17.7648),
        ('LIPO', 21.6134, 15.9862, 15.1898),
        ('LIPO', 26.7691, 19.8439, 16.4741),
        ('LIPO', 29.1815, 19.1760, 16.1647),
        ('LIPO', 25.7991, 21.2142, 14.4885),
        ('LIPO', 27.2749, 20.4918, 14.1040),
        ('LIPO', 24.5309, 21.1357, 13.3778),
        ('LIPO', 23.9291, 23.3071, 12.7534),
        ('LIPO', 28.3741, 16.1960, 19.2624),
        ('LIPO', 25.9780, 15.8915, 15.4798),
        ('LIPO', 24.0567, 16.5469, 12.4554),
        ('LIPO', 24.9876, 13.9507, 11.4379),
        ('LIPO', 23.5648, 15.1316, 15.5744),
        ('AROM', 22.3755, 20.4466, 6.0574),
        ('LIPO', 23.8808, 20.8064, 10.9438),
        ('LIPO', 26.5379, 14.6626, 7.0921),
        ('LIPO', 25.3444, 15.7184, 4.3437),
        ('HDON', 23.2381, 11.3095, 3.8981),
        ('HDON', 23.5066, 12.4505, 3.4648),
        ('LIPO', 15.6875, 21.3788, 11.7234),
        ('AROM', 24.7701, 15.9888, 2.9569),
        ('LIPO', 24.4291, 19.9667, 5.6702),
        ('LIPO', 21.8130, 19.6184, 3.9196),
        ('LIPO', 21.1677, 14.4833, 8.2061),
        ('LIPO', 22.5933, 12.8494, 6.5634),
        ('LIPO', 19.3598, 20.0103, 15.4982),
        ('LIPO', 29.4735, 17.1238, 15.8539),
        ('LIPO', 28.4358, 15.3412, 17.9075),
        ('LIPO', 22.6927, 21.2814, 8.5115),
        ('LIPO', 20.6163, 21.6970, 6.7951),
    ]


@pytest.fixture
def probe_pharmacophore():
    return [
        ('HACC', -1.7076, 2.2682, 22.7126),
        ('HDON', -0.0317, 4.6294, 22.4973),
        ('HDON', 3.9657, -4.2182, 19.4535),
        ('AROM', -4.7420, 5.8751, 25.9774),
        ('HDON', 3.7079, 4.2267, 24.4837),
        ('NEGC', 3.8882, 3.0747, 24.4667),
        ('HDON', 1.1487, 5.4662, 26.7621),
        ('NEGC', 0.5654, 6.4209, 26.4331),
        ('NEGC', 0.5855, 4.6737, 27.4061),
        ('NEGC', 2.1046, 4.7328, 25.3711),
        ('LIPO', -0.9244, -0.6887, 25.2947),
        ('LIPO', 1.3499, 0.1450, 23.3679),
        ('LIPO', -0.2041, 0.2003, 24.0542),
        ('LIPO', -2.4841, -1.2942, 25.0740),
        ('LIPO', -5.5423, 0.1490, 24.4589),
        ('LIPO', -2.2450, 5.3313, 29.7474),
        ('LIPO', -2.8484, 1.4969, 28.8643),
        ('LIPO', -0.7917, 3.4229, 28.1349),
        ('LIPO', 3.7043, 2.5349, 25.3078),
        ('LIPO', 1.4568, 1.3568, 25.5041),
        ('LIPO', 5.7106, 2.4540, 25.7827),
        ('LIPO', 5.4877, -3.2515, 15.6444),
        ('LIPO', -0.1181, 1.2310, 19.4872),
        ('LIPO', -0.7460, 1.0889, 21.8024),
        ('LIPO', 6.2402, 3.1961, 22.3875),
        ('LIPO', -4.3778, 0.0278, 30.1472),
        ('HACC', 3.3976, -3.8262, 16.5781),
        ('HDON', 2.9732, -3.4432, 13.1614),
        ('LIPO', 1.9186, -2.1014, 15.8518),
        ('LIPO', 1.0415, -0.2159, 14.6303),
        ('LIPO', 2.9680, 1.9694, 17.1866),
        ('LIPO', 4.0759, -3.1885, 20.6854),
        ('LIPO', 3.7108, -1.5415, 20.7330),
        ('LIPO', 5.8368, 0.5870, 12.9128),
        ('LIPO', 5.5328, 2.5229, 16.8755),
        ('LIPO', 6.1775, 0.2833, 16.1111),
        ('LIPO', -4.4496, -1.8999, 26.3564),
    ]


def test_distances():
    points = [
        ('HACC', -1.7076, 2.2682, 22.7126),
        ('HDON', -0.0317, 4.6294, 22.4973),
        ('HDON', 3.9657, -4.2182, 19.4535),
    ]
    result = distances(points)

    expected = [[0.0, 2.9034910607749422, 9.213112973365735],
                [2.9034910607749422, 0.0, 10.174672032060789],
                [9.213112973365735, 10.174672032060789, 0.0]]

    assert_array_almost_equal(result, expected)


@pytest.fixture
def aligner(reference_pharmacophore, probe_pharmacophore):
    return Aligner(reference_pharmacophore, probe_pharmacophore)


class TestAligner_SomeContent(object):
    def test_construct(self, aligner, reference_pharmacophore, probe_pharmacophore):
        assert aligner.reference == reference_pharmacophore
        assert aligner.probe == probe_pharmacophore

    def test_pairs(self, aligner):
        aligner.calculate_possible_pairs()

        assert len(aligner.pairs) > 0
        assert len(aligner.nodes) > 0

    def test_cliques(self, aligner):
        aligner.calculate_possible_pairs()

        cliques = aligner.cliques()
        assert len(cliques) > 0

    def test_cliqued_pharmacophores(self, aligner):
        aligner.calculate_possible_pairs()
        aligner.cliques()
        (reference, probe) = aligner.cliqued_pharmacophores()

        assert len(reference) > 0
        assert len(probe) > 0

    def test_transformation(self, aligner):
        rmsd, matrix = aligner.transformation()

        assert rmsd > 0
        assert len(matrix) == 4


def test_aligner_transformation_SingleSamePoint():
    ref = probe = [('LIPO', 1., 1., 1.)]
    # TODO should return matrix which does nothing
    aligner = Aligner(ref, probe)
    with pytest.raises(NoOverlapFound):
        aligner.transformation()


def test_align(aligner, probe_pharmacophore):
    rmsd, matrix = aligner.transformation()
    aligned_probe = align(probe_pharmacophore, matrix)

    assert len(aligned_probe) == len(probe_pharmacophore)
