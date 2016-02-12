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
import os
from cStringIO import StringIO

from nose.tools import eq_

from kripodb.hdf5 import DistanceMatrix
import kripodb.script as script
from tests.test_pairs import tmpname


def test_pairs_subcommand_defaults():
    parser = script.make_parser()

    args = parser.parse_args(['pairs', '--fragmentsdbfn', 'fragdb', 'fp1', 'fp2', 'outfn'])

    eq_(args.func, script.pairs_run)

    fargs = vars(args)
    del(fargs['func'])
    expected = {
        'out_format': 'hdf5',
        'cutoff': 0.45,
        'out_file': 'outfn',
        'fragmentsdbfn': 'fragdb',
        'mean_onbit_density': 0.01,
        'precision': 65535,
        'nomemory': False,
        'fingerprintsfn2': 'fp2',
        'fingerprintsfn1': 'fp1'
    }
    eq_(fargs, expected)


def test_meanbitdensity():
    out = StringIO()

    script.meanbitdensity_run('data/fingerprints.sqlite', out)

    eq_(out.getvalue(), '0.01228\n')


def test_distmatrix_import_run():
    output_fn = tmpname()

    tsv = '''frag_id1	frag_id2	score
2mlm_2W7_frag1	2mlm_2W7_frag2	0.5877164873731594
2mlm_2W7_frag2	3wvm_STE_frag1	0.4633096818493935
'''
    inputfile = StringIO(tsv)

    try:
        script.distmatrix_import_run(inputfile=inputfile,
                                     distmatrixfn=output_fn,
                                     fragmentsdb='data/fragments.sqlite',
                                     precision=65535,
                                     nrrows=2)

        distmatrix = DistanceMatrix(output_fn)
        rows = [r for r in distmatrix]
        distmatrix.close()
        expected = [('2mlm_2W7_frag1', '2mlm_2W7_frag2', 0.5877164873731594), ('2mlm_2W7_frag2', '3wvm_STE_frag1', 0.4633096818493935)]
        eq_(rows, expected)
    finally:
        os.remove(output_fn)


def test_distmatrix_export_run():
    outputfile = StringIO()
    script.distmatrix_export_run('data/distances.h5', outputfile)

    # go back to start of file
    outputfile.seek(0)
    output = outputfile.getvalue()
    expected = 'frag_id1\tfrag_id2\tscore\n' \
               '2mlm_2W7_frag1\t2mlm_2W7_frag2\t0.5877164873731594\n' \
               '2mlm_2W7_frag2\t3wvm_STE_frag1\t0.4633096818493935\n'
    assert output.startswith(expected)
