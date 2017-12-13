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

import pytest
from rdkit.Chem.AllChem import MolFromSmiles

from kripodb.pharmacophores import PharmacophoresDb
from kripodb.webservice import server
from kripodb.pairs import open_similarity_matrix
from kripodb.version import __version__
from kripodb.webservice.server import KripodbJSONEncoder


class TestKripodbJSONEncoder(object):
    def test_mol(self):
        encoder = KripodbJSONEncoder()
        obj = MolFromSmiles('[*]C1OC(COP(=O)([O-])OP(=O)([O-])OCC2OC(N3C=CCC(C(N)=O)=C3)C(O)C2O)C(O)C1[*]')
        result = encoder.default(obj)
        assert result.endswith('END\n')


@pytest.fixture
def similarity_matrix():
    matrix = open_similarity_matrix('data/similarities.frozen.h5')
    yield matrix
    matrix.close()


@pytest.fixture
def fragsdb_filename():
    return 'data/fragments.sqlite'


@pytest.fixture
def pharmacophores_db():
    db = PharmacophoresDb('data/pharmacophores.h5')
    yield db
    db.close()


@pytest.fixture
def app(similarity_matrix, fragsdb_filename, pharmacophores_db):
    return server.wsgi_app(similarity_matrix, fragsdb_filename, pharmacophores_db)


@pytest.fixture
def expected_fragments_info():
    return [
        {'smiles': '[*]C1OC(COP(=O)([O-])OP(=O)([O-])OCC2OC(N3C=CCC(C(N)=O)=C3)C(O)C2O)C(O)C1[*]',
         'pdb_code': '3j7u',
         'pdb_title': 'Catalase structure determined by electron crystallography of thin 3D crystals',
         'atom_codes': 'PA,O1A,O2A,O5B,C5B,C4B,O4B,C3B,O3B,C2B,C1B,O3,PN,O1N,O2N,O5D,C5D,C4D,O4D,C3D,O3D,C2D,O2D,C1D,N1N,C2N,C3N,C7N,O7N,N7N,C4N,C5N,C6N',
         'uniprot_acc': 'P00432',
         'prot_chain': 'A', 'het_seq_nr': 602, 'het_code': 'NDP', 'prot_name': 'Catalase',
         'ec_number': '1.11.1.6', 'frag_nr': 24, 'frag_id': '3j7u_NDP_frag24', 'rowid': 7059,
         'uniprot_name': 'Catalase', 'nr_r_groups': 2, 'het_chain': 'A', 'hash_code': '6ef5a609fb192dba'}
    ]


@pytest.fixture
def expected_fragments_info_with_mol(expected_fragments_info):
    expected_fragments_info[0][
        'mol'] = '3j7u_NDP_frag24\n     RDKit          3D\n\n 35 37  0  0  0  0  0  0  0  0999 V2000\n  -15.1410  -11.1250  -79.4200 P   0  0  0  0  0  0  0  0  0  0  0  0\n  -14.6900  -10.9960  -80.8600 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -16.5040  -11.6890  -79.0770 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -14.9990   -9.6870  -78.7060 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -15.1870   -8.4550  -79.4050 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -14.6700   -7.3160  -78.5260 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -13.2400   -7.2390  -78.5880 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -15.2130   -5.9510  -78.9460 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -16.1600   -5.4570  -77.9880 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -14.0000   -5.0420  -79.0650 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -14.1790   -3.8250  -78.3260 R   0  0  0  0  0  1  0  0  0  0  0  0\n  -12.8370   -5.8690  -78.5180 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -11.5470   -5.6210  -79.2410 R   0  0  0  0  0  1  0  0  0  0  0  0\n  -14.0270  -11.9960  -78.6490 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -14.1810  -13.5930  -78.4870 P   0  0  0  0  0  0  0  0  0  0  0  0\n  -14.5480  -14.2030  -79.8230 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -15.0330  -13.8500  -77.2690 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -12.6800  -14.0730  -78.1770 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -12.1840  -14.2350  -76.8490 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -11.1340  -13.1670  -76.6050 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -11.6880  -11.8550  -76.6770 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -10.5070  -13.2750  -75.2350 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -9.4070  -14.1780  -75.3000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -10.0970  -11.8400  -74.9280 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -8.6920  -11.6460  -75.1050 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -10.8280  -10.9760  -75.9460 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -11.5890   -9.8540  -75.3660 N   0  0  0  0  0  0  0  0  0  0  0  0\n  -12.7860  -10.0630  -74.7850 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -13.5340   -9.0090  -74.2510 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -14.8620   -9.2740  -73.5990 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -15.1890  -10.4300  -73.3940 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -15.6600   -8.2650  -73.2400 N   0  0  0  0  0  0  0  0  0  0  0  0\n  -13.0230   -7.5870  -74.3390 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -11.7130   -7.4960  -74.9740 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -11.0640   -8.6200  -75.4710 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  1  3  1  0\n  1  4  1  0\n  1 14  1  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  6  8  1  0\n  7 12  1  0\n  8  9  1  0\n  8 10  1  0\n 10 11  1  0\n 10 12  1  0\n 12 13  1  0\n 14 15  1  0\n 15 16  2  0\n 15 17  1  0\n 15 18  1  0\n 18 19  1  0\n 19 20  1  0\n 20 21  1  0\n 20 22  1  0\n 21 26  1  0\n 22 23  1  0\n 22 24  1  0\n 24 25  1  0\n 24 26  1  0\n 26 27  1  0\n 27 28  1  0\n 27 35  1  0\n 28 29  2  0\n 29 30  1  0\n 29 33  1  0\n 30 31  2  0\n 30 32  1  0\n 33 34  1  0\n 34 35  2  0\nM  CHG  2   3  -1  17  -1\nM  END\n'
    return expected_fragments_info


def test_get_similar_fragments(app):
    fragment_id = '3j7u_NDP_frag24'
    cutoff = 0.85

    with app.app.test_request_context():
        result = server.get_similar_fragments(fragment_id, cutoff, 1000)
        expected = [
            {'query_frag_id': '3j7u_NDP_frag24', 'hit_frag_id': '3j7u_NDP_frag23', 'score': 0.8991},
        ]
    assert result == expected


def test_get_similar_fragments_notfound(app):
    fragment_id = 'foo-bar'
    cutoff = 0.85

    with app.app.test_request_context():
        response = server.get_similar_fragments(fragment_id, cutoff, 1000)
        assert response.status_code == 404
        body = response_json(response)
        assert fragment_id in body['detail']
        assert fragment_id == body['identifier']


def test_get_fragments__fragid(app, expected_fragments_info):
    fragment_id = '3j7u_NDP_frag24'

    with app.app.test_request_context():
        result = server.get_fragments(fragment_ids=[fragment_id])
        # Remove RDKit mol as equals is to strict when comparing Molecule objects
        del result[0]['mol']
        assert result == expected_fragments_info


def response_json(response):
    return response.body


def test_get_fragments__fragid_notfound(app):
    fragment_id = 'foo-bar'
    with app.app.test_request_context():
        response = server.get_fragments(fragment_ids=[fragment_id])
        assert response.status_code == 404
        body = response_json(response)
        assert fragment_id in body['detail']
        assert [fragment_id] == body['absent_identifiers']


def test_fragments_by_id_withsomenotfound(app, expected_fragments_info_with_mol):
    absent_fragment_id = 'foo-bar'
    present_fragment_id = '3j7u_NDP_frag24'
    with app.app.test_request_context():
        response = server.get_fragments(fragment_ids=[present_fragment_id, absent_fragment_id])
        assert response.status_code == 404
        body = response_json(response)
        assert body['fragments'] == expected_fragments_info_with_mol
        assert absent_fragment_id in body['detail']
        assert [absent_fragment_id] == body['absent_identifiers']


def test_get_fragments__pdbcode(app):
    pdb_code = '3j7u'

    with app.app.test_request_context():
        result = server.get_fragments(pdb_codes=[pdb_code])
        assert len(result) == 32


def test_get_fragments__pdbcode_notfound(app):
    pdb_code = 'foo-bar'
    with app.app.test_request_context():
        response = server.get_fragments(pdb_codes=[pdb_code])
        assert response.status_code == 404
        body = response_json(response)
        assert pdb_code in body['detail']
        assert [pdb_code] == body['absent_identifiers']


def test_get_version():
    result = server.get_version()

    expected = {'version': __version__}
    assert result == expected


def test_wsgi_app(similarity_matrix, fragsdb_filename, pharmacophores_db, app):
    assert app.app.config['similarities'] == similarity_matrix
    assert app.app.config['fragments'] == fragsdb_filename
    assert app.app.config['pharmacophores'] == pharmacophores_db
    assert app.app.json_encoder == KripodbJSONEncoder


def test_get_fragment_svg(app):
    fragment_id = '3j7u_NDP_frag24'
    with app.app.test_request_context():
        result = server.get_fragment_svg(fragment_id, 400, 150)
        assert b'NH' in result.data
        assert b'400' in result.data
        assert b'150' in result.data
        assert result.mimetype == 'image/svg+xml'


def test_get_fragment_svg_notfound(app):
    fragment_id = 'foo-bar'
    with app.app.test_request_context():
        response = server.get_fragment_svg(fragment_id, 400, 150)
        assert response.status_code == 404
        body = response_json(response)
        assert fragment_id in body['detail']
        assert fragment_id == body['identifier']


def test_get_fragment_phar(app):
    fragment_id = '3j7u_NDP_frag24'
    with app.app.test_request_context():
        result = server.get_fragment_phar(fragment_id)
        assert b'LIPO' in result.data
        assert result.mimetype == 'text/plain'


def test_get_fragment_phar_notfound(app):
    fragment_id = 'foo-bar'
    with app.app.test_request_context():
        response = server.get_fragment_phar(fragment_id)
        assert response.status_code == 404
        body = response_json(response)
        assert fragment_id in body['detail']
        assert fragment_id == body['identifier']
