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
from six import StringIO
from nose.tools import eq_, assert_raises

import kripodb.script.dive as dive
from .utils import FrozenSimilarityMatrixInMemory


def uniprot_mapping():
    return  StringIO('''Entry	Gene names  (primary )	Protein families	Cross-reference (PDB)
P81186	napA	Prokaryotic molybdopterin-containing oxidoreductase family, NasA/NapA/NarB subfamily	2JIM;2JIO;2JIP;2JIQ;2JIR;2NAP;2V3V;2V45;
''')


def test_add_uniprot():
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
    mapping = uniprot_mapping()

    dive.add_uniprot(data, mapping)

    expected = {
        "2v3v_LCP_frag1": {
            "Path": [],
            "Coordinates": [0.2864670506986739, -0.24369131639863245, -0.23773760870587482],
            "Categories": [
                "2v3v",
                "LCP",
                "frag1",
                "P81186",
                "gene:napA",
                "Prokaryotic molybdopterin-containing oxidoreductase family",
                "NasA/NapA/NarB subfamily",
            ],
            "Properties": []
        }
    }
    eq_(data, expected)


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
        eq_(result, expected)


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
        eq_(result, expected)


def test_dive_sphere():
    inputfile = 'data/fragments.sqlite'
    outputfile = StringIO()
    onlyfrag1 = False

    dive.dive_sphere(inputfile, outputfile, onlyfrag1)

    # Dont test layout, only test a single item has correct props
    sphere = json.loads(outputfile.getvalue())
    example = sphere['3wtj_TH4_frag1']
    eq_(example['Categories'], ['3wtj', 'TH4'])
    eq_(len(example['Coordinates']), 3)


def test_dive_sphere_frag1():
    inputfile = 'data/fragments.sqlite'
    outputfile = StringIO()
    onlyfrag1 = True

    dive.dive_sphere(inputfile, outputfile, onlyfrag1)

    sphere = json.loads(outputfile.getvalue())
    # Dont test layout, only test a single item has correct props
    example = sphere['3wtj_TH4_frag1']
    eq_(example['Categories'], ['3wtj', 'TH4'])
    eq_(len(example['Coordinates']), 3)

    with assert_raises(KeyError):
        sphere['3wtj_TH4_frag2']


def test_dive_merge_uniprot_nodata():
    mapping = uniprot_mapping()
    data = {}

    dive.dive_merge_uniprot(mapping, data)

    expected = {}
    eq_(data, expected)


def test_dive_merge_uniprot():
    mapping = uniprot_mapping()
    data = {
        '2v3v_LCP_frag1': {
            'pdb': '2v3v',
            'uniprot': 'P81186'
        }
    }

    dive.dive_merge_uniprot(mapping, data)

    expected = {
        '2v3v_LCP_frag1': {
            'pdb': '2v3v',
            'uniprot': 'P81186',
            'gene': 'napA',
            'families': [
                'Prokaryotic molybdopterin-containing oxidoreductase family',
                'NasA/NapA/NarB subfamily'
            ]
        }
    }
    eq_(data, expected)


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

    expected = '2v3v_LCP_frag1 pdb:2v3v het:LCP fragment:1 "title:A NEW CATALYTIC MECHANISM OF PERIPLASMIC NITRATE REDUCTASE FROM DESULFOVIBRIO DESULFURICANS ATCC 27774 FROM CRYSTALLOGRAPHIC AND EPR DATA AND BASED ON DETAILED ANALYSIS OF THE SIXTH LIGAND" smiles:Nc1nc2N[C@@H]3O[C@H](CO[P@@](O)(=O)O[P@@](O)(=O)OC[C@H]4O[C@H]([C@H](O)[C@@H]4O)n4cnc5c4nc(N)[nH]c5=O)C(S)=C(S)[C@@H]3Nc2c(=O)[nH]1 740.56 uniprot:P81186 "protein:Periplasmic nitrate reductase" "organism:Desulfovibrio desulfuricans" gene:napA "family0:Prokaryotic molybdopterin-containing oxidoreductase family" "family1:NasA/NapA/NarB subfamily"\n'
    eq_(props_file.getvalue(), expected)