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

from intbitset import intbitset

import logging
from nose.tools import eq_, raises, assert_raises
from mock import call, Mock
from rdkit.Chem import MolFromSmiles, MolToSmiles

import kripodb.db as db


def test_adapt_intbitset():
    bs = intbitset([1, 3, 5, 8])

    result = db.adapt_intbitset(bs)

    expected = bs.fastdump()
    eq_(result, expected)


def test_convert_intbitset():
    bs = intbitset([1, 3, 5, 8])

    result = db.convert_intbitset(bs.fastdump())

    eq_(result, bs)


class TestFastInserter(object):
    def setUp(self):
        self.cursor = Mock()
        self.unit = db.FastInserter(self.cursor)

    def test_with(self):
        with self.unit:
            self.cursor.execute.assert_has_calls([call('PRAGMA journal_mode=WAL'),
                                                  call('PRAGMA synchronous=OFF')])

        self.cursor.execute.assert_has_calls([call('PRAGMA journal_mode=DELETE'),
                                              call('PRAGMA synchronous=FULL')])


class TestFragmentsDBEmpty(object):
    def setUp(self):
        self.fdb = db.FragmentsDb(':memory:')

    def test_id2label(self):
        eq_(self.fdb.id2label(), {})

    def test_label2id(self):
        eq_(self.fdb.label2id(), {})

    @raises(KeyError)
    def test_getitem_keyerror(self):
        key = 'id1'
        self.fdb[key]

    def test_by_pdb_code(self):
        pdb_code = '1kvm'

        fragments = self.fdb.by_pdb_code(pdb_code)

        eq_(fragments, [])

    def test_add_fragments_from_shelve_weirdid(self):
        result = self.fdb.add_fragment_from_shelve('1muu-GDX', {})
        eq_(result, None)

    def test_add_fragments_from_shelve_weirdid2(self):
        result = self.fdb.add_fragment_from_shelve('1muu-GDX-B', {})
        eq_(result, None)

    def test_len(self):
        eq_(len(self.fdb), 0)


class TestFragmentsDBFilled(object):
    def setUp(self):
        self.fdb = db.FragmentsDb(':memory:')
        self.myshelve ={
            '1muu-GDX-frag7': {
                'atomCodes': 'C5D,O5D,PA,O1A,O2A,O3A,PB,O2B,O3B,O1B,C1*,O5*,C5*,C6*,O6A,O6B,C2*,O2*,C3*,O3*,C4*,O4*',
                'hashcode': '0d6ced7ce686f4da',
                'ligID': '1muu-A-GDX-1005-B',
                'numRgroups': '1'
            }
        }
        self.fdb.add_fragments_from_shelve(self.myshelve)

        self.mol = MolFromSmiles('[*]COP(=O)([O-])OP(=O)([O-])OC1OC(C(=O)[O-])C(O)C(O)C1O')
        self.mol.SetProp('_Name', '1muu_GDX_frag7')
        self.fdb.add_molecule(self.mol)
        self.expected_fragment = {
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
        self.pdbs = [{
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
        self.fdb.add_pdbs(self.pdbs)

    def test_getitem(self):
        fragment = self.fdb['1muu_GDX_frag7']

        eq_(MolToSmiles(fragment['mol']), '[*]COP(=O)([O-])OP(=O)([O-])OC1OC(C(=O)[O-])C(O)C(O)C1O')
        del fragment['mol']

        eq_(fragment, self.expected_fragment)

    def test_id2label(self):
        eq_(self.fdb.id2label(), {1: '1muu_GDX_frag7'})

    def test_label2id(self):
        eq_(self.fdb.label2id(), {'1muu_GDX_frag7': 1})

    def test_by_pdb_code(self):
        pdb_code = '1muu'

        fragments = self.fdb.by_pdb_code(pdb_code)

        del fragments[0]['mol']
        eq_(fragments, [self.expected_fragment])

    def test_len(self):
        eq_(len(self.fdb), 1)


class TestHetSeqNr(object):
    def test_isnumber(self):
        fdb = db.FragmentsDb(':memory:')
        myshelve = {
            '1muu-GDX-frag7': {
                'atomCodes': 'C5D,O5D,PA,O1A,O2A,O3A,PB,O2B,O3B,O1B,C1*,O5*,C5*,C6*,O6A,O6B,C2*,O2*,C3*,O3*,C4*,O4*',
                'hashcode': '0d6ced7ce686f4da',
                'ligID': '1muu-A-GDX-1005-B',
                'numRgroups': '1'
            }
        }
        fdb.add_fragments_from_shelve(myshelve)
        fdb.add_pdbs([{
            'structureTitle': '2.0 A crystal structure of GDP-mannose dehydrogenase',
            'ecNo': '1.1.1.132',
            'uniprotAcc': 'P11759',
            'compound': 'GDP-mannose 6-dehydrogenase',
            'uniprotRecommendedName': 'GDP-mannose 6-dehydrogenase',
            'chainId': 'A',
            'structureId': '1muu'
        }])

        fragment = fdb['1muu_GDX_frag7']

        eq_(fragment['het_seq_nr'], 1005)

    def test_nan(self):
        fdb = db.FragmentsDb(':memory:')
        myshelve = {
            '1hoo-GNP-frag1': {
                'atomCodes': 'PG,O1G,O2G,O3G,N3B,PB,O1B,O2B,O3A,PA,O1A,O2A,O5*,C5*,C4*,O4*,C3*,O3*,C2*,O2*,C1*,N9,C8,N7,C5,C6,O6,N1,C2,N2,N3,C4',
                'hashcode': 'be4ce041f2a35721',
                'ligID': '1hoo-A-GNP-432B-A',
                'numRgroups': '0'
            }
        }
        fdb.add_fragments_from_shelve(myshelve)
        fdb.add_pdbs([{
            'chainId': 'A',
            'structureId': '1hoo',
            'structureTitle': 'STRUCTURE OF GUANINE NUCLEOTIDE (GPPCP) COMPLEX OF ADENYLOSUCCINATE SYNTHETASE FROM E. COLI AT PH 6.5 AND 25 DEGREES CELSIUS',
            'ecNo': '6.3.4.4',
            'uniprotAcc': 'P0A7D4',
            'uniprotRecommendedName': 'Adenylosuccinate synthetase',
            'compound': 'ADENYLOSUCCINATE SYNTHETAS',
        }])

        fragment = fdb['1hoo_GNP_frag1']

        eq_(fragment['het_seq_nr'], 432)


class TestIntbitsetDictEmpty(object):
    def setUp(self):
        self.fdb = db.FingerprintsDb(':memory:')
        self.bitsets = self.fdb.as_dict(100)
        self.bid = 'id1'
        self.bs = intbitset([1, 3, 5, 8])

    def tearDown(self):
        self.fdb.close()

    def test_default_number_of_bits(self):
        fdb = db.FingerprintsDb(':memory:')
        bitsets = db.IntbitsetDict(fdb)

        eq_(bitsets.number_of_bits, None)

    def test_get_number_of_bits(self):
        eq_(self.bitsets.number_of_bits, 100)

    def test_set_number_of_bits(self):
        self.bitsets.number_of_bits = 200

        eq_(self.bitsets.number_of_bits, 200)

    def test_delete_number_of_bits(self):
        del self.bitsets.number_of_bits

        eq_(self.bitsets.number_of_bits, None)

    def test_len_empty(self):
        eq_(len(self.bitsets), 0)

    def test_contains_false(self):
        assert 'id1' not in self.bitsets

    def test_update(self):
        other = {self.bid: self.bs}

        self.bitsets.update(other)

        result = {k: v for k, v in self.bitsets.iteritems()}
        eq_(result, other)

    def test_getitem_keyerror(self):
        with assert_raises(KeyError) as e:
            self.bitsets['id1']
        eq_(e.exception.message, 'id1')



class TestIntbitsetDictFilled(object):
    def setUp(self):
        self.fdb = db.FingerprintsDb(':memory:')
        self.bitsets = db.IntbitsetDict(self.fdb, 100)
        self.bid = 'id1'
        self.bs = intbitset([1, 3, 5, 8])
        self.bitsets[self.bid] = self.bs

    def test_getitem(self):
        result = self.bitsets['id1']

        eq_(result, self.bs)

    def test_len_filled(self):
        eq_(len(self.bitsets), 1)

    def test_contains_true(self):
        assert 'id1' in self.bitsets

    def test_del(self):
        del self.bitsets['id1']

        eq_(len(self.bitsets), 0)

    def test_keys(self):
        result = self.bitsets.keys()

        expected = ['id1']
        eq_(result, expected)

    def test_iteritems(self):
        result = {k: v for k, v in self.bitsets.iteritems()}

        expected = {self.bid: self.bs}
        eq_(result, expected)

    def test_iteritems_startswith(self):
        self.bitsets['someid'] = self.bs

        result = {k: v for k, v in self.bitsets.iteritems_startswith('id')}

        expected = {self.bid: self.bs}
        eq_(result, expected)
        assert 'someid' not in result

    def test_itervalues(self):
        result = [v for v in self.bitsets.itervalues()]

        expected = [self.bs]
        eq_(result, expected)

    def test_materialize(self):
        result = self.bitsets.materialize()

        expected = {self.bid: self.bs}
        eq_(result, expected)
