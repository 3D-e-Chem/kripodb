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
import socket

import pytest
import requests
import signal
from pandas.util.testing import assert_frame_equal
import pandas as pd

from kripodb.canned import similarities, fragments_by_pdb_codes, fragments_by_id
from kripodb.webservice.client import IncompleteFragments
from kripodb.webservice.server import serve_app


def test_similarities():
    queries = pd.Series(['3j7u_NDP_frag24'])

    result = similarities(queries, 'data/similarities.h5', 0.85)

    expected = [
        {'query_frag_id': '3j7u_NDP_frag24', 'hit_frag_id': '3j7u_NDP_frag23', 'score': 0.8991},
    ]
    assert_frame_equal(result, pd.DataFrame(expected))


def test_similarities_limitof1():
    queries = pd.Series(['3j7u_NDP_frag24'])

    result = similarities(queries, 'data/similarities.h5', 0.55, 1)

    expected = [
        {'query_frag_id': '3j7u_NDP_frag24', 'hit_frag_id': '3j7u_NDP_frag23', 'score': 0.8991},
    ]
    assert_frame_equal(result, pd.DataFrame(expected))


def test_fragments_by_pdb_codes():
    pdb_codes = pd.Series(['2n2k'])

    result = fragments_by_pdb_codes(pdb_codes, 'data/fragments.sqlite')

    # ignoring molecules
    result.drop('mol', axis=1, inplace=True, errors='ignore')
    expected = [{
        'nr_r_groups': 0, 'smiles': 'CC1(C)C=C(C[S-])C(C)(C)[NH+]1O', 'pdb_code': '2n2k',
        'atom_codes': 'O1,N1,C1,C2,C3,C4,S1,C5,C6,C7,C8,C9', 'het_code': 'MTN', 'hash_code': 'd491952cd7c9dc30',
        'frag_nr': 1, 'frag_id': '2n2k_MTN_frag1', 'rowid': 175992, 'het_chain': 'A', 'het_seq_nr': 101,
        'prot_chain': 'A', 'uniprot_acc': 'P0CG48', 'uniprot_name': 'Polyubiquitin-C', 'prot_name': 'ubiquitin',
        'ec_number': None,
        'pdb_title': 'Ensemble structure of the closed state of Lys63-linked diubiquitin in the absence of a ligand',
    }, {
        'nr_r_groups': 1, 'smiles': '[*]C[S-]', 'pdb_code': '2n2k', 'atom_codes': 'C4,S1', 'het_code': 'MTN',
        'hash_code': '8b8dc32f7e8a9db3', 'frag_nr': 2, 'frag_id': '2n2k_MTN_frag2', 'rowid': 175950, 'het_chain': 'A',
        'het_seq_nr': 101,
        'prot_chain': 'A', 'uniprot_acc': 'P0CG48', 'uniprot_name': 'Polyubiquitin-C', 'prot_name': 'ubiquitin',
        'ec_number': None,
        'pdb_title': 'Ensemble structure of the closed state of Lys63-linked diubiquitin in the absence of a ligand',
    }, {
        'nr_r_groups': 1, 'smiles': '[*]C1=CC(C)(C)[NH+](O)C1(C)C', 'pdb_code': '2n2k',
        'atom_codes': 'O1,N1,C1,C2,C3,C5,C6,C7,C8,C9', 'het_code': 'MTN', 'hash_code': '17c58abf7bdf33ba',
        'frag_nr': 3, 'frag_id': '2n2k_MTN_frag3', 'rowid': 175971, 'het_chain': 'A', 'het_seq_nr': 101,
        'prot_chain': 'A', 'uniprot_acc': 'P0CG48', 'uniprot_name': 'Polyubiquitin-C', 'prot_name': 'ubiquitin',
        'ec_number': None,
        'pdb_title': 'Ensemble structure of the closed state of Lys63-linked diubiquitin in the absence of a ligand',
    }]
    assert_frame_equal(result, pd.DataFrame(expected))


def test_fragments_by_pdb_codes_with_prefix():
    pdb_codes = pd.Series(['3wxj'])

    result = fragments_by_pdb_codes(pdb_codes, 'data/fragments.sqlite', 'prefix_')

    # ignoring molecules
    result.drop('prefix_mol', axis=1, inplace=True, errors='ignore')

    expected = [{
        'prefix_nr_r_groups': 0, 'prefix_smiles': 'O=P([O-])([O-])OCC(O)CO', 'prefix_pdb_code': '3wxj',
        'prefix_atom_codes': 'O1,C1,C2,O2,C3,O1P,O4P,O2P,O3P,P', 'prefix_het_code': 'G3P', 'prefix_hash_code': 'ee9013689ff298d4',
        'prefix_frag_nr': 1, 'prefix_frag_id': '3wxj_G3P_frag1', 'prefix_rowid': 352104, 'prefix_het_chain': 'B', 'prefix_het_seq_nr': 601,
        'prefix_prot_chain': 'B', 'prefix_uniprot_acc': 'D3KVM3', 'prefix_uniprot_name': None, 'prefix_prot_name': 'Glycerol kinase',
        'prefix_ec_number': '2.7.1.30',
        'prefix_pdb_title': 'Crystal structure of trypanosoma brucei gambiense glycerol kinase in complex with glycerol 3-phosphate',
    }]
    assert_frame_equal(result, pd.DataFrame(expected))


@pytest.fixture
def open_port():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("", 0))
    s.listen(1)
    port = s.getsockname()[1]
    s.close()
    return port


@pytest.fixture
def server_url(open_port):
    pid = os.fork()
    if pid == 0:
        serve_app(matrix='data/similarities.h5', db='data/fragments.sqlite', internal_port=open_port)
    else:
        while True:
            try:
                requests.get('http://localhost:{0}'.format(open_port))
                break
            except requests.ConnectionError:
                continue
    yield 'http://localhost:{0}'.format(open_port)
    os.kill(pid, signal.SIGTERM)


@pytest.mark.skip(reason='Server is still booting when this test is run, so it fails for the wrong reason')
def test_fragments_by_pdb_codes__usingwebservice_withbadids(server_url):
    pdb_codes = pd.Series(['0000'])

    with pytest.raises(IncompleteFragments) as e:
        fragments_by_pdb_codes(pdb_codes, server_url)

    assert e.value.fragments == pd.DataFrame()
    assert e.value.absent_identifiers == ['0000']


def test_fragments_by_id():
    frag_ids = pd.Series(['2n2k_MTN_frag1'])

    result = fragments_by_id(frag_ids, 'data/fragments.sqlite')

    # ignoring molecules
    result.drop('mol', axis=1, inplace=True)
    expected = [{
        'nr_r_groups': 0, 'smiles': 'CC1(C)C=C(C[S-])C(C)(C)[NH+]1O', 'pdb_code': '2n2k',
        'atom_codes': 'O1,N1,C1,C2,C3,C4,S1,C5,C6,C7,C8,C9', 'het_code': 'MTN', 'hash_code': 'd491952cd7c9dc30',
        'frag_nr': 1, 'frag_id': '2n2k_MTN_frag1', 'rowid': 175992, 'het_chain': 'A', 'het_seq_nr': 101,
        'prot_chain': 'A', 'uniprot_acc': 'P0CG48', 'uniprot_name': 'Polyubiquitin-C', 'prot_name': 'ubiquitin',
        'ec_number': None,
        'pdb_title': 'Ensemble structure of the closed state of Lys63-linked diubiquitin in the absence of a ligand',
    }]
    assert_frame_equal(result, pd.DataFrame(expected))


def test_fragments_by_id_with_prefix():
    frag_ids = pd.Series(['2n2k_MTN_frag1'])

    result = fragments_by_id(frag_ids, 'data/fragments.sqlite', 'prefix_')

    # ignoring molecules
    result.drop('prefix_mol', axis=1, inplace=True)
    expected = [{
        'prefix_nr_r_groups': 0, 'prefix_smiles': 'CC1(C)C=C(C[S-])C(C)(C)[NH+]1O', 'prefix_pdb_code': '2n2k',
        'prefix_atom_codes': 'O1,N1,C1,C2,C3,C4,S1,C5,C6,C7,C8,C9', 'prefix_het_code': 'MTN',
        'prefix_hash_code': 'd491952cd7c9dc30', 'prefix_frag_nr': 1, 'prefix_frag_id': '2n2k_MTN_frag1',
        'prefix_rowid': 175992, 'prefix_het_chain': 'A', 'prefix_het_seq_nr': 101, 'prefix_prot_chain': 'A',
        'prefix_uniprot_acc': 'P0CG48', 'prefix_uniprot_name': 'Polyubiquitin-C', 'prefix_prot_name': 'ubiquitin',
        'prefix_ec_number': None,
        'prefix_pdb_title': 'Ensemble structure of the closed state of Lys63-linked diubiquitin in the absence of a ligand',
    }]
    assert_frame_equal(result, pd.DataFrame(expected))