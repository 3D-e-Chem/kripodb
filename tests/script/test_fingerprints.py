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

from six import StringIO

import kripodb.script as script
import kripodb.script.fingerprints


def test_pairs_subcommand_defaults():
    parser = script.make_parser()

    args = parser.parse_args(['fingerprints', 'similarities', '--fragmentsdbfn', 'fragdb', 'fp1', 'fp2', 'outfn'])

    assert args.func == kripodb.script.fingerprints.pairs_run

    fargs = vars(args)
    del(fargs['func'])
    expected = {
        'subcommand': 'fingerprints',
        'out_format': 'hdf5',
        'cutoff': 0.45,
        'out_file': 'outfn',
        'fragmentsdbfn': 'fragdb',
        'mean_onbit_density': 0.01,
        'nomemory': False,
        'fingerprintsfn2': 'fp2',
        'fingerprintsfn1': 'fp1',
        'ignore_upper_triangle': False,
    }
    assert fargs == expected


def test_meanbitdensity():
    out = StringIO()

    kripodb.script.fingerprints.meanbitdensity_run('data/fingerprints.sqlite', out)

    assert out.getvalue() == '0.0077683\n'


