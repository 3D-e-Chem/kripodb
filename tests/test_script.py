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
import os
from six import StringIO

from nose.tools import eq_
from numpy.testing import assert_array_almost_equal

from kripodb.hdf5 import SimilarityMatrix
import kripodb.script as script
from tests.test_pairs import tmpname


def test_pairs_subcommand_defaults():
    parser = script.make_parser()

    args = parser.parse_args(['fingerprints', 'similarities', '--fragmentsdbfn', 'fragdb', 'fp1', 'fp2', 'outfn'])

    eq_(args.func, script.pairs_run)

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
    eq_(fargs, expected)


def test_meanbitdensity():
    out = StringIO()

    script.meanbitdensity_run('data/fingerprints.sqlite', out)

    eq_(out.getvalue(), '0.0077683\n')


def test_simmatrix_import_run():
    output_fn = tmpname()

    tsv = '''frag_id1	frag_id2	score
2mlm_2W7_frag1	2mlm_2W7_frag2	0.5877164873731594
2mlm_2W7_frag2	3wvm_STE_frag1	0.4633096818493935
'''
    inputfile = StringIO(tsv)

    try:
        script.simmatrix_import_run(inputfile=inputfile,
                                    format='tsv',
                                    simmatrixfn=output_fn,
                                    fragmentsdb='data/fragments.sqlite',
                                    nrrows=2)

        simmatrix = SimilarityMatrix(output_fn)
        result = [r for r in simmatrix]
        simmatrix.close()
        expected = [('2mlm_2W7_frag1', '2mlm_2W7_frag2xx', 0.5877), ('2mlm_2W7_frag2', '3wvm_STE_frag1', 0.4633)]
        assert_array_almost_equal([r[2] for r in result], [r[2] for r in expected], 3)
        eq_([(r[0], r[1],) for r in result], [(r[0], r[1],) for r in result])
    finally:
        if os.path.exists(output_fn):
            os.remove(output_fn)


def test_simmatrix_import_run_ignore_upper_triangle():
    output_fn = tmpname()

    tsv = '''frag_id1	frag_id2	score
2mlm_2W7_frag1	2mlm_2W7_frag1	1.0000000000000000
2mlm_2W7_frag2	2mlm_2W7_frag2	1.0000000000000000
2mlm_2W7_frag1	2mlm_2W7_frag2	0.5877164873731594
2mlm_2W7_frag2	3wvm_STE_frag1	0.4633096818493935
2mlm_2W7_frag2	2mlm_2W7_frag1	0.5877164873731594
3wvm_STE_frag1	2mlm_2W7_frag2	0.4633096818493935
'''
    inputfile = StringIO(tsv)

    try:
        script.simmatrix_import_run(inputfile=inputfile,
                                    format='tsv',
                                    simmatrixfn=output_fn,
                                    fragmentsdb='data/fragments.sqlite',
                                    nrrows=2,
                                    ignore_upper_triangle=True)

        simmatrix = SimilarityMatrix(output_fn)
        result = [r for r in simmatrix]
        simmatrix.close()
        print(result)
        expected = [('2mlm_2W7_frag1', '2mlm_2W7_frag2xx', 0.5877), ('2mlm_2W7_frag2', '3wvm_STE_frag1', 0.4633)]
        assert_array_almost_equal([r[2] for r in result], [r[2] for r in expected], 3)
        eq_([(r[0], r[1],) for r in result], [(r[0], r[1],) for r in result])
    finally:
        if os.path.exists(output_fn):
            os.remove(output_fn)


def test_simmatrix_export_run():
    outputfile = StringIO()
    script.simmatrix_export_run('data/similarities.h5', outputfile)

    # go back to start of file
    outputfile.seek(0)
    output = outputfile.getvalue()
    expected = 'frag_id1\tfrag_id2\tscore\n' \
               '2mlm_2W7_frag1\t2mlm_2W7_frag2\t0.5878\n' \
               '2mlm_2W7_frag2\t3wvm_STE_frag1\t0.4634\n'
    assert output.startswith(expected)


def test_read_fpneighpairs_file():
    text = StringIO('''Compounds similar to 2xry_FAD_frag4:
2xry_FAD_frag4   1.0000
3cvv_FAD_frag3   0.5600
Compounds similar to 1wnt_NAP_frag1:
1wnt_NAP_frag1   1.0000
1wnt_NAP_frag3   0.8730
''')

    result = list(script.read_fpneighpairs_file(text))

    expected = [('2xry_FAD_frag4', '3cvv_FAD_frag3', 0.56), ('1wnt_NAP_frag1', '1wnt_NAP_frag3', 0.873)]
    eq_(result, expected)


def test_simmatrix_importfpneigh_run():
    output_fn = tmpname()

    tsv = '''Compounds similar to 2mlm_2W7_frag1:
2mlm_2W7_frag1   1.0000
2mlm_2W7_frag2   0.5877
Compounds similar to 2mlm_2W7_frag2:
2mlm_2W7_frag2   1.0000
3wvm_STE_frag1   0.4633
'''
    inputfile = StringIO(tsv)

    try:
        script.simmatrix_importfpneigh_run(inputfile=inputfile,
                                           simmatrixfn=output_fn,
                                           fragmentsdb='data/fragments.sqlite',
                                           nrrows=3)

        simmatrix = SimilarityMatrix(output_fn)
        rows = [r for r in simmatrix]
        simmatrix.close()
        expected = [(u'2mlm_2W7_frag1', u'2mlm_2W7_frag2', 0.5877), (u'2mlm_2W7_frag2', u'3wvm_STE_frag1', 0.4633)]
        eq_(rows, expected)
    finally:
        os.remove(output_fn)


def test_simmatrix_importfpneigh_run_ignore_upper_triangle():
    output_fn = tmpname()

    tsv = '''Compounds similar to 2mlm_2W7_frag1:
2mlm_2W7_frag1   1.0000
2mlm_2W7_frag2   0.5877
Compounds similar to 2mlm_2W7_frag2:
2mlm_2W7_frag2   1.0000
2mlm_2W7_frag1   0.5877
3wvm_STE_frag1   0.4633
'''
    inputfile = StringIO(tsv)

    try:
        script.simmatrix_importfpneigh_run(inputfile=inputfile,
                                           simmatrixfn=output_fn,
                                           fragmentsdb='data/fragments.sqlite',
                                           nrrows=3,
                                           ignore_upper_triangle=True)

        simmatrix = SimilarityMatrix(output_fn)
        rows = [r for r in simmatrix]
        simmatrix.close()
        expected = [(u'2mlm_2W7_frag1', u'2mlm_2W7_frag2', 0.5877), (u'2mlm_2W7_frag2', u'3wvm_STE_frag1', 0.4633)]
        eq_(rows, expected)
    finally:
        os.remove(output_fn)


def test_fpneigh2tsv_run():
    fpneigh_in = '''Compounds similar to 2mlm_2W7_frag1:
2mlm_2W7_frag1   1.0000
2mlm_2W7_frag2   0.5877
Compounds similar to 2mlm_2W7_frag2:
2mlm_2W7_frag2   1.0000
3wvm_STE_frag1   0.4633
'''
    inputfile = StringIO(fpneigh_in)

    outputfile = StringIO()

    script.fpneigh2tsv_run(inputfile, outputfile)

    expected = '''frag_id1\tfrag_id2\tscore
2mlm_2W7_frag1\t2mlm_2W7_frag2\t0.5877
2mlm_2W7_frag2\t3wvm_STE_frag1\t0.4633
'''
    eq_(outputfile.getvalue(), expected)
