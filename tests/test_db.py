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
import sqlite3

from intbitset import intbitset
from mock import call, Mock
import pytest
from rdkit.Chem import MolFromSmiles, MolToSmiles
import six

import kripodb.db as db


def test_adapt_intbitset():
    bs = intbitset([1, 3, 5, 8])

    result = db.adapt_intbitset(bs)

    expected = bs.fastdump()
    assert result == expected


def test_convert_intbitset():
    bs = intbitset([1, 3, 5, 8])

    result = db.convert_intbitset(bs.fastdump())

    assert result == bs


class TestFastInserter(object):
    @pytest.fixture
    def cursor(self):
        return Mock()

    @pytest.fixture
    def fast_inserter(self, cursor):
        return db.FastInserter(cursor)

    def test_with(self, cursor, fast_inserter):
        with fast_inserter:
            cursor.execute.assert_has_calls([call('PRAGMA journal_mode=WAL'),
                                             call('PRAGMA synchronous=OFF')])

        cursor.execute.assert_has_calls([call('PRAGMA journal_mode=DELETE'),
                                         call('PRAGMA synchronous=FULL')])


@pytest.fixture
def fragmentsdb():
    return db.FragmentsDb(':memory:')


class TestFragmentsDBEmpty(object):
    def test_id2label(self, fragmentsdb):
        assert fragmentsdb.id2label() == {}

    def test_label2id(self, fragmentsdb):
        assert fragmentsdb.label2id() == {}

    def test_getitem_keyerror(self, fragmentsdb):
        key = 'id1'
        with pytest.raises(KeyError):
            fragmentsdb[key]

    def test_by_pdb_code(self, fragmentsdb):
        pdb_code = '1kvm'
        with pytest.raises(LookupError):
            fragmentsdb.by_pdb_code(pdb_code)

    def test_add_fragments_from_shelve_weirdid(self, fragmentsdb):
        result = fragmentsdb.add_fragment_from_shelve('1muu-GDX', {})
        assert result is None

    def test_add_fragments_from_shelve_weirdid2(self, fragmentsdb):
        result = fragmentsdb.add_fragment_from_shelve('1muu-GDX-B', {})
        assert result is None

    def test_len(self, fragmentsdb):
        assert len(fragmentsdb) == 0


@pytest.fixture
def myshelve():
    return {
        '1muu-GDX-frag7': {
            'atomCodes': 'C5D,O5D,PA,O1A,O2A,O3A,PB,O2B,O3B,O1B,C1*,O5*,C5*,C6*,O6A,O6B,C2*,O2*,C3*,O3*,C4*,O4*',
            'hashcode': '0d6ced7ce686f4da',
            'ligID': '1muu-A-GDX-1005-B',
            'numRgroups': '1'
        }
    }


@pytest.fixture
def filled_fragmentsdb(fragmentsdb, myshelve):
    fragmentsdb.add_fragments_from_shelve(myshelve)

    mol = MolFromSmiles('[*]COP(=O)([O-])OP(=O)([O-])OC1OC(C(=O)[O-])C(O)C(O)C1O')
    mol.SetProp('_Name', '1muu_GDX_frag7')
    fragmentsdb.add_molecule(mol)
    pdbs = [{
        'chainId': 'A',
        'structureId': '1muu',
        'structureTitle': '2.0 A crystal structure of GDP-mannose dehydrogenase',
        'ecNo': '1.1.1.132',
        'uniprotAcc': 'P11759',
        'compound': 'GDP-mannose 6-dehydrogenase',
        'uniprotRecommendedName': 'GDP-mannose 6-dehydrogenase',
    }, {
        # pdbs which has no fragment should be skipped
        'chainId': 'A',
        'structureId': '2n2k',
        'structureTitle': 'Ensemble structure of the closed state of Lys63-linked diubiquitin in the absence of a ligand',
        'ecNo': None,
        'uniprotAcc': 'P0CG48',
        'compound': 'ubiquitin',
        'uniprotRecommendedName': 'Polyubiquitin-C',
    }]
    fragmentsdb.add_pdbs(pdbs)
    return fragmentsdb


@pytest.fixture
def expected_fragment():
    return {
        'nr_r_groups': 1,
        'smiles': '[*]COP(=O)([O-])OP(=O)([O-])OC1OC(C(=O)[O-])C(O)C(O)C1O',
        'pdb_code': '1muu',
        'atom_codes': 'C5D,O5D,PA,O1A,O2A,O3A,PB,O2B,O3B,O1B,C1*,O5*,C5*,C6*,O6A,O6B,C2*,O2*,C3*,O3*,C4*,O4*',
        'het_code': 'GDX',
        'hash_code': '0d6ced7ce686f4da',
        'frag_nr': 7,
        'frag_id': '1muu_GDX_frag7',
        'rowid': 1,
        'het_seq_nr': 1005,
        'het_chain': 'B',
        'prot_chain': 'A',
        'pdb_title': '2.0 A crystal structure of GDP-mannose dehydrogenase',
        'prot_name': 'GDP-mannose 6-dehydrogenase',
        'ec_number': '1.1.1.132',
        'uniprot_acc': 'P11759',
        'uniprot_name': 'GDP-mannose 6-dehydrogenase',
    }


class TestFragmentsDBFilled(object):
    def test_getitem(self, filled_fragmentsdb, expected_fragment):
        fragment = filled_fragmentsdb['1muu_GDX_frag7']

        assert MolToSmiles(fragment['mol']) == '[*]COP(=O)([O-])OP(=O)([O-])OC1OC(C(=O)[O-])C(O)C(O)C1O'
        del fragment['mol']

        assert fragment == expected_fragment

    def test_id2label(self, filled_fragmentsdb):
        assert filled_fragmentsdb.id2label() == {1: '1muu_GDX_frag7'}

    def test_label2id(self, filled_fragmentsdb):
        assert filled_fragmentsdb.label2id() == {'1muu_GDX_frag7': 1}

    def test_by_pdb_code(self, filled_fragmentsdb, expected_fragment):
        pdb_code = '1muu'

        fragments = filled_fragmentsdb.by_pdb_code(pdb_code)

        del fragments[0]['mol']
        assert fragments == [expected_fragment]

    def test_len(self, filled_fragmentsdb):
        assert len(filled_fragmentsdb) == 1

    def test_duplicate(self, filled_fragmentsdb, myshelve):
        with pytest.raises(sqlite3.IntegrityError):
            filled_fragmentsdb.add_fragments_from_shelve(myshelve)


class TestHetSeqNr(object):
    def test_isnumber(self, filled_fragmentsdb):
        fragment = filled_fragmentsdb['1muu_GDX_frag7']

        assert fragment['het_seq_nr'] == 1005

    def test_nan(self, fragmentsdb):
        myshelve = {
            '1hoo-GNP-frag1': {
                'atomCodes': 'PG,O1G,O2G,O3G,N3B,PB,O1B,O2B,O3A,PA,O1A,O2A,O5*,C5*,C4*,O4*,C3*,O3*,C2*,O2*,C1*,N9,C8,N7,C5,C6,O6,N1,C2,N2,N3,C4',
                'hashcode': 'be4ce041f2a35721',
                'ligID': '1hoo-A-GNP-432B-A',
                'numRgroups': '0'
            }
        }
        fragmentsdb.add_fragments_from_shelve(myshelve)
        fragmentsdb.add_pdbs([{
            'chainId': 'A',
            'structureId': '1hoo',
            'structureTitle': 'STRUCTURE OF GUANINE NUCLEOTIDE (GPPCP) COMPLEX OF ADENYLOSUCCINATE SYNTHETASE FROM E. COLI AT PH 6.5 AND 25 DEGREES CELSIUS',
            'ecNo': '6.3.4.4',
            'uniprotAcc': 'P0A7D4',
            'uniprotRecommendedName': 'Adenylosuccinate synthetase',
            'compound': 'ADENYLOSUCCINATE SYNTHETAS',
        }])

        fragment = fragmentsdb['1hoo_GNP_frag1']

        assert fragment['het_seq_nr'] == 432


@pytest.fixture
def fingerprintsdb():
    fdb = db.FingerprintsDb(':memory:')
    yield fdb
    fdb.close()


@pytest.fixture
def bitsets(fingerprintsdb):
    return fingerprintsdb.as_dict(100)


class TestIntbitsetDictEmpty(object):
    def test_default_number_of_bits(self, fingerprintsdb):
        bitsets = db.IntbitsetDict(fingerprintsdb)

        assert bitsets.number_of_bits is None

    def test_get_number_of_bits(self, bitsets):
        assert bitsets.number_of_bits == 100

    def test_set_number_of_bits(self, bitsets):
        bitsets.number_of_bits = 200

        assert bitsets.number_of_bits == 200

    def test_delete_number_of_bits(self, bitsets):
        del bitsets.number_of_bits

        assert bitsets.number_of_bits is None

    def test_len_empty(self, bitsets):
        assert len(bitsets) == 0

    def test_contains_false(self, bitsets):
        assert 'id1' not in bitsets

    def test_update(self, bitsets):
        bs = intbitset([1, 3, 5, 8])
        other = {'id1': bs}

        bitsets.update(other)

        result = {k: v for k, v in six.iteritems(bitsets)}
        assert result == other

    def test_getitem_keyerror(self, bitsets):
        with pytest.raises(KeyError) as e:
            bitsets['id1']
        assert e.value.args == ('id1',)


@pytest.fixture
def sample_intbitset():
    return intbitset([1, 3, 5, 8])


@pytest.fixture
def filled_bitsets(bitsets, sample_intbitset):
    bid = 'id1'
    bitsets[bid] = sample_intbitset
    return bitsets


class TestIntbitsetDictFilled(object):

    def test_getitem(self, filled_bitsets, sample_intbitset):
        result = filled_bitsets['id1']

        assert result == sample_intbitset

    def test_len_filled(self, filled_bitsets):
        assert len(filled_bitsets) == 1

    def test_contains_true(self, filled_bitsets):
        assert 'id1' in filled_bitsets

    def test_del(self, filled_bitsets):
        del filled_bitsets['id1']

        assert len(filled_bitsets) == 0

    def test_keys(self, filled_bitsets):
        result = list(filled_bitsets.keys())

        expected = ['id1']
        assert result == expected

    def test_iteritems(self, filled_bitsets, sample_intbitset):
        result = {k: v for k, v in six.iteritems(filled_bitsets)}

        expected = {'id1': sample_intbitset}
        assert result == expected

    def test_iteritems_startswith(self, filled_bitsets, sample_intbitset):
        filled_bitsets['someid'] = sample_intbitset

        result = {k: v for k, v in filled_bitsets.iteritems_startswith('id')}

        expected = {'id1': sample_intbitset}
        assert result == expected
        assert 'someid' not in result

    def test_itervalues(self, filled_bitsets, sample_intbitset):
        result = [v for v in six.itervalues(filled_bitsets)]

        expected = [sample_intbitset]
        assert result == expected

    def test_materialize(self, filled_bitsets, sample_intbitset):
        result = filled_bitsets.materialize()

        expected = {'id1': sample_intbitset}
        assert result == expected
