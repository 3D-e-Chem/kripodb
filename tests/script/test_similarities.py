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

from numpy.testing import assert_array_almost_equal
import pytest
from six import StringIO

from kripodb.db import FragmentsDb
from kripodb.hdf5 import SimilarityMatrix
import kripodb.script.similarities as script
from ..utils import tmpname


def test_simmatrix_import_run():
    output_fn = tmpname()

    tsv = '''frag_id1	frag_id2	score
2mlm_2W7_frag1	2mlm_2W7_frag2	0.5877164873731594
2mlm_2W7_frag2	3wvm_STE_frag1	0.4633096818493935
'''
    inputfile = StringIO(tsv)

    try:
        script.simmatrix_import_run(inputfile=inputfile,
                                    inputformat='tsv',
                                    simmatrixfn=output_fn,
                                    fragmentsdb='data/fragments.sqlite',
                                    nrrows=2)

        simmatrix = SimilarityMatrix(output_fn)
        result = [r for r in simmatrix]
        simmatrix.close()
        expected = [('2mlm_2W7_frag1', '2mlm_2W7_frag2xx', 0.5877), ('2mlm_2W7_frag2', '3wvm_STE_frag1', 0.4633)]
        assert_array_almost_equal([r[2] for r in result], [r[2] for r in expected], 3)
        assert [(r[0], r[1],) for r in result] == [(r[0], r[1],) for r in result]
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
                                    inputformat='tsv',
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
        assert [(r[0], r[1],) for r in result] == [(r[0], r[1],) for r in result]
    finally:
        if os.path.exists(output_fn):
            os.remove(output_fn)


def test_simmatrix_export_run():
    outputfile = StringIO()
    script.simmatrix_export_run('data/similarities.h5', outputfile, False, False, None)

    # go back to start of file
    outputfile.seek(0)
    output = outputfile.getvalue()
    expected = 'frag_id1\tfrag_id2\tscore\n' \
               '2mlm_2W7_frag1\t2mlm_2W7_frag2\t0.5878\n' \
               '2mlm_2W7_frag2\t3wvm_STE_frag1\t0.4634\n'
    assert output.startswith(expected)


def test_simmatrix_export_run_noheader():
    outputfile = StringIO()
    script.simmatrix_export_run('data/similarities.h5', outputfile, True, False, None)

    # go back to start of file
    outputfile.seek(0)
    output = outputfile.getvalue()
    expected = '2mlm_2W7_frag1\t2mlm_2W7_frag2\t0.5878\n' \
               '2mlm_2W7_frag2\t3wvm_STE_frag1\t0.4634\n'
    assert output.startswith(expected)


def test_simmatrix_export_run_frag1():
    outputfile = StringIO()
    script.simmatrix_export_run('data/similarities.h5', outputfile, False, True, None)

    # go back to start of file
    outputfile.seek(0)
    output = outputfile.getvalue()
    expected = 'frag_id1\tfrag_id2\tscore\n' \
               '2mlr_PX4_frag1\t2n2k_MTN_frag1\t0.4661\n' \
               '2mm3_CHO_frag1\t3wt8_RET_frag1\t0.4689\n'
    assert output.startswith(expected)


def test_simmatrix_export_run_pdb():
    pdbio = StringIO('2mlm\n')
    outputfile = StringIO()
    script.simmatrix_export_run('data/similarities.h5', outputfile, False, False, pdbio)

    # go back to start of file
    outputfile.seek(0)
    output = outputfile.getvalue()
    expected = 'frag_id1\tfrag_id2\tscore\n' \
               '2mlm_2W7_frag1\t2mlm_2W7_frag2\t0.5878\n' \
               '2mlm_2W7_frag2\t2mlm_2W7_frag4\t0.4609\n'
    assert output.startswith(expected)


def test_simmatrix_export_run_frag1_and_pdb():
    pdbio = StringIO('2mm3\n3wt8\n')
    outputfile = StringIO()
    script.simmatrix_export_run('data/similarities.h5', outputfile, False, True, pdbio)

    # go back to start of file
    outputfile.seek(0)
    output = outputfile.getvalue()
    expected = 'frag_id1\tfrag_id2\tscore\n' \
               '2mm3_CHO_frag1\t3wt8_RET_frag1\t0.4689\n' \
               '2mm3_CHO_frag1\t2mm3_GCH_frag1\t0.571\n'
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
    assert result == expected


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
        assert rows == expected
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
        assert rows == expected
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
    assert outputfile.getvalue() == expected


@pytest.fixture
def output_simmatrix_fn():
    fn = tmpname()
    yield fn
    if os.path.exists(fn):
        os.remove(fn)


@pytest.fixture
def keep_fragments_db():
    fn = tmpname()
    with FragmentsDb(fn) as db:
        shelve = {
            '2mm3-CHO-frag1': {
                'atomCodes': 'C1,C2,C3,O3,C4,C5,C6,C7,O7,C8,C9,C10,C11,C12...',
                'hashcode': '47df1bbad2834aca',
                'ligID': '2mm3-A-CHO-202-B',
                'numRgroups': '0'
            }
        }
        db.add_fragments_from_shelve(shelve)
    yield fn
    if os.path.exists(fn):
        os.remove(fn)


def test_simmatrix_filter_keep(keep_fragments_db, output_simmatrix_fn):
    script.simmatrix_filter('data/similarities.h5', output_simmatrix_fn, keep_fragments_db, None)

    result = SimilarityMatrix(output_simmatrix_fn)
    nr_pairs = len(result.pairs)
    result.close()
    assert nr_pairs == 93


def test_simmatrix_filter_skip(output_simmatrix_fn):
    labels2skip = StringIO('2mm3_CHO_frag1\n')

    script.simmatrix_filter('data/similarities.h5', output_simmatrix_fn, None, labels2skip)

    result = SimilarityMatrix(output_simmatrix_fn)
    nr_pairs = len(result.pairs)
    result.close()
    assert nr_pairs == 11857
