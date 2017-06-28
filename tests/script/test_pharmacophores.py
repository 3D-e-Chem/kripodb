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

import json
import sys
from six import StringIO

from numpy.testing import assert_array_almost_equal
import pytest

import kripodb.script as script
import kripodb.script.pharmacophores as pscript


def test_align_subcommand_defaults():
    parser = script.make_parser()

    args = parser.parse_args([
        'pharmacophores',
        'align',
        '--pharmacophoresdb', 'phardb',
        'fp1',
        'fp2'
    ])

    assert args.func == pscript.align_run
    assert args.output == sys.stdout

    fargs = vars(args)
    del (fargs['func'])
    del (fargs['output'])
    expected = {
        'subcommand': 'pharmacophores',
        'pharmacophoresdb': 'phardb',
        'reference': 'fp1',
        'probe': 'fp2',
        'break_num_cliques': 3000,
        'cutoff': 1.0,
    }
    assert fargs == expected


def test_align_run():
    out = StringIO()

    pscript.align_run('data/pharmacophores.h5', '3jat_G2P_frag1', '2n2k_MTN_frag1', 1.0, 3000, out)

    expected_rmsd = 0.3575553403252601
    expected_matrix = [
        [0.5725482868782117, -0.8178244756237104, -0.05789288612280747, 299.4765555664346],
        [0.6968909249848213, 0.5226488480675373, -0.4911020467148112, 429.4803697809889],
        [0.4318929240756618, 0.2408346085687844, 0.8691761578925732, 222.43723988101792],
        [0.0, 0.0, 0.0, 1.0]
    ]

    result = json.loads(out.getvalue())

    assert result[0] == pytest.approx(expected_rmsd)
    assert_array_almost_equal(result[1], expected_matrix)
