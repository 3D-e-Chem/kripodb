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
import numpy as np
from six import StringIO, BytesIO
import pytest
from mock import patch

import kripodb.dive as dive
from .utils import FrozenSimilarityMatrixInMemory


@pytest.fixture
def mock_fetch_response():
    mresponse = BytesIO()
    mresponse.write(b'structureId,source\n')
    mresponse.write(b'"2n2k","Homo sapiens"\n')
    mresponse.seek(0)
    return mresponse


@patch('kripodb.pdb.urlopen')
def test_dive_export(mocked_urlopen, mock_fetch_response):
    mocked_urlopen.return_value = mock_fetch_response

    fragmentsdb = 'data/fragments.sqlite'
    uniprot_annot = StringIO(
        'Entry\tGene names  (primary )\tProtein families\tCross-reference (PDB)' + '\n' +
        'P0CG48\tUBC\tUbiquitin family\t1C3T;2N2K' + '\n'
    )
    pdbtag = StringIO('2n2k' + '\n')
    pdbtag.name = 'mytag'

    propnames = StringIO()
    props = StringIO()

    dive.dive_export(fragmentsdb, uniprot_annot, [pdbtag], propnames, props)

    assert '["pdb", "het", "fragment", "title", "smiles", "weight", "uniprot", "protein", "organism", "gene", "pdbtag", "family0", "family1", "family2", "family3", "family4"]' == propnames.getvalue()
    props_lines = props.getvalue().split('\n')
    result = list(filter(lambda d: d.startswith('2n2k_MTN_frag1'), props_lines))[0]
    expected = '2n2k_MTN_frag1 pdb:2n2k het:MTN fragment:1 "title:Ensemble structure of the closed state of Lys63-linked diubiquitin in the absence of a ligand" smiles:CC1(C)C=C(C[S-])C(C)(C)[NH+]1O 170.17 uniprot:P0CG48 "protein:Polyubiquitin-C" "organism:Homo sapiens" "gene:UBC" pdbtag:mytag "family0:Ubiquitin family"'
    assert result == expected


def test_dense_dump_allfrags():
    with FrozenSimilarityMatrixInMemory() as matrix:
        labels = ['a', 'b', 'c', 'd']
        data = [
            [0.0, 0.9, 0.5, 0.0],
            [0.9, 0.0, 0.6, 0.0],
            [0.5, 0.6, 0.0, 0.7],
            [0.0, 0.0, 0.7, 0.0],
        ]
        matrix.from_array(np.array(data), labels)

        result = list(dive.dense_dump_iter(matrix, frag1only=False))
        expected = [
            (u'a', u'b', 0.9),
            (u'a', u'c', 0.5),
            (u'b', u'c', 0.6),
            (u'c', u'd', 0.7),
        ]
        assert result == expected


def test_dense_dump_frag1only():
    with FrozenSimilarityMatrixInMemory() as matrix:
        labels = ['a_frag1', 'a_frag2', 'c_frag1', 'c_frag2']
        data = [
            [0.0, 0.9, 0.5, 0.0],
            [0.9, 0.0, 0.6, 0.0],
            [0.5, 0.6, 0.0, 0.7],
            [0.0, 0.0, 0.7, 0.0],
        ]
        matrix.from_array(np.array(data), labels)

        result = list(dive.dense_dump_iter(matrix, frag1only=True))
        expected = [(u'a_frag1', u'c_frag1', 0.5)]
        assert result == expected


def test_dive_sphere():
    inputfile = 'data/fragments.sqlite'
    outputfile = StringIO()
    onlyfrag1 = False

    dive.dive_sphere(inputfile, outputfile, onlyfrag1)

    # Dont test layout, only test a single item has correct props
    sphere = json.loads(outputfile.getvalue())
    example = sphere['3wtj_TH4_frag1']
    assert example['Categories'], ['3wtj' == 'TH4']
    assert len(example['Coordinates']) == 3


def test_dive_sphere_frag1():
    inputfile = 'data/fragments.sqlite'
    outputfile = StringIO()
    onlyfrag1 = True

    dive.dive_sphere(inputfile, outputfile, onlyfrag1)

    sphere = json.loads(outputfile.getvalue())
    # Dont test layout, only test a single item has correct props
    example = sphere['3wtj_TH4_frag1']
    assert example['Categories'], ['3wtj' == 'TH4']
    assert len(example['Coordinates']) == 3

    with pytest.raises(KeyError):
        sphere['3wtj_TH4_frag2']


def test_dump_props():
    props = {
        '2v3v_LCP_frag1': {
            'pdb': '2v3v',
            'het': 'LCP',
            'fragment': 1,
            'title': 'A NEW CATALYTIC MECHANISM OF PERIPLASMIC NITRATE REDUCTASE FROM DESULFOVIBRIO DESULFURICANS ATCC 27774 FROM CRYSTALLOGRAPHIC AND EPR DATA AND BASED ON DETAILED ANALYSIS OF THE SIXTH LIGAND',
            'smiles': 'Nc1nc2N[C@@H]3O[C@H](CO[P@@](O)(=O)O[P@@](O)(=O)OC[C@H]4O[C@H]([C@H](O)[C@@H]4O)n4cnc5c4nc(N)[nH]c5=O)C(S)=C(S)[C@@H]3Nc2c(=O)[nH]1',
            'weight': 740.56,
            'protein': 'Periplasmic nitrate reductase',
            'organism': 'Desulfovibrio desulfuricans',
            'uniprot': 'P81186',
            'gene': 'napA',
            'families': [
                'Prokaryotic molybdopterin-containing oxidoreductase family',
                'NasA/NapA/NarB subfamily'
            ]
        }
    }
    props_file = StringIO()

    dive.dump_props(props, props_file)

    expected = '2v3v_LCP_frag1 pdb:2v3v het:LCP fragment:1 "title:A NEW CATALYTIC MECHANISM OF PERIPLASMIC NITRATE REDUCTASE FROM DESULFOVIBRIO DESULFURICANS ATCC 27774 FROM CRYSTALLOGRAPHIC AND EPR DATA AND BASED ON DETAILED ANALYSIS OF THE SIXTH LIGAND" smiles:Nc1nc2N[C@@H]3O[C@H](CO[P@@](O)(=O)O[P@@](O)(=O)OC[C@H]4O[C@H]([C@H](O)[C@@H]4O)n4cnc5c4nc(N)[nH]c5=O)C(S)=C(S)[C@@H]3Nc2c(=O)[nH]1 740.56 uniprot:P81186 "protein:Periplasmic nitrate reductase" "organism:Desulfovibrio desulfuricans" "gene:napA"  "family0:Prokaryotic molybdopterin-containing oxidoreductase family" "family1:NasA/NapA/NarB subfamily"\n'
    assert props_file.getvalue() == expected
