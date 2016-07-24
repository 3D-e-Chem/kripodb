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

import numpy as np
from kripodb.frozen import FrozenSimilarityMatrix
from six import StringIO
from nose.tools import eq_

from kripodb.script.cclustera import enrich_fragments, dense_dump_iter
from tests.test_pairs import tmpname


def test_enrich_fragments():
    data = {
        "2v3v_LCP_frag1": {
            "Path": [],
            "Coordinates": [0.2864670506986739, -0.24369131639863245, -0.23773760870587482],
            "Categories": [
                "LCP",
                "2v3v",
                "frag1"
            ],
            "Properties": []
        }
    }
    mapping = StringIO('''Entry	Gene names  (primary )	Protein families	Cross-reference (PDB)
P81186	napA	Prokaryotic molybdopterin-containing oxidoreductase family, NasA/NapA/NarB subfamily	2JIM;2JIO;2JIP;2JIQ;2JIR;2NAP;2V3V;2V45;
''')

    enrich_fragments(data, mapping)

    expected = {
        "2v3v_LCP_frag1": {
            "Path": [],
            "Coordinates": [0.2864670506986739, -0.24369131639863245, -0.23773760870587482],
            "Categories": [
                "2v3v",
                "LCP",
                "frag1",
                "P81186",
                "gene_napA",
                "Prokaryotic molybdopterin-containing oxidoreductase family",
                "NasA/NapA/NarB subfamily",
            ],
            "Properties": []
        }
    }
    eq_(data, expected)


def test_dense_dump_allfrags():
    matrix_fn = tmpname()
    matrix = FrozenSimilarityMatrix(matrix_fn, 'a', driver='H5FD_CORE', driver_core_backing_store=0)
    labels = ['a', 'b', 'c', 'd']
    data = [
        [0.0, 0.9, 0.5, 0.0],
        [0.9, 0.0, 0.6, 0.0],
        [0.5, 0.6, 0.0, 0.7],
        [0.0, 0.0, 0.7, 0.0],
    ]
    matrix.from_array(np.array(data), labels)

    result = list(dense_dump_iter(matrix, frag1only=False))
    expected = [
        (u'a', u'b', 0.9),
        (u'a', u'c', 0.5),
        (u'a', u'd', 0.0),
        (u'b', u'c', 0.6),
        (u'b', u'd', 0.0),
        (u'c', u'd', 0.7),
    ]
    eq_(result, expected)

    matrix.close()


def test_dense_dump_frag1only():
    matrix_fn = tmpname()
    matrix = FrozenSimilarityMatrix(matrix_fn, 'a', driver='H5FD_CORE', driver_core_backing_store=0)
    labels = ['a_frag1', 'a_frag2', 'c_frag1', 'c_frag2']
    data = [
        [0.0, 0.9, 0.5, 0.0],
        [0.9, 0.0, 0.6, 0.0],
        [0.5, 0.6, 0.0, 0.7],
        [0.0, 0.0, 0.7, 0.0],
    ]
    matrix.from_array(np.array(data), labels)

    result = list(dense_dump_iter(matrix, frag1only=True))
    expected = [(u'a_frag1', u'c_frag1', 0.5)]
    eq_(result, expected)

    matrix.close()