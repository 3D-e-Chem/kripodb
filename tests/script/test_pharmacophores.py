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

import pytest
from numpy.testing import assert_array_almost_equal

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

    (rmsd, matrix) = json.loads(out.getvalue())

    assert rmsd == pytest.approx(0.19767105)
    expected = [[5.253691e-01, -8.392010e-01, 1.404603e-01, 3.040672e+02],
                [-8.467798e-01, -4.994923e-01, 1.829521e-01, 4.662223e+02],
                [-8.337473e-02, -2.150564e-01, -9.730362e-01, 2.896292e+02],
                [0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00]]
    assert_array_almost_equal(matrix, expected, 4)
