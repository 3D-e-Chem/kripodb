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

from collections import MutableMapping
import sqlite3
import logging
import zlib
from intbitset import intbitset
from rdkit.Chem import MolToMolBlock, MolFromMolBlock, Mol, MolToSmiles

ATTR_NUMBER_OF_BITS = 'number_of_bits'


def adapt_intbitset(ibs):
    return ibs.fastdump()


def convert_intbitset(s):
    ibs = intbitset()
    ibs.fastload(s)
    return ibs


def adapt_molblockgz(mol):
    return zlib.compress(MolToMolBlock(mol))


def convert_molblockgz(molgz):
    return MolFromMolBlock(zlib.decompress(molgz))


sqlite3.register_adapter(intbitset, adapt_intbitset)
sqlite3.register_converter('intbitset', convert_intbitset)
sqlite3.register_adapter(Mol, adapt_molblockgz)
sqlite3.register_converter('molblockgz', convert_molblockgz)


class FastInserter(object):
    def __init__(self, cursor):
        """

        Use with to make inserting faster, but less sage

        Parameters
        ----------
        cursor sqlite3 cursor

        Returns
        -------

        """
        self.cursor = cursor

    def __enter__(self):
        # increase insert speed, this is less safe
        self.cursor.execute('PRAGMA journal_mode=WAL')
        self.cursor.execute('PRAGMA synchronous=OFF')

    def __exit__(self, exc_type, exc_val, exc_tb):
        # switch back to default journal, so db file can be read-only and is safe again
        self.cursor.execute('PRAGMA journal_mode=DELETE')
        self.cursor.execute('PRAGMA synchronous=FULL')


class SqliteDb(object):
    def __init__(self, filename):
        """Connect to fragment db.

        Database is created if it does not exist.

        :param filename: Sqlite filename
        :return:
        """
        self.filename = filename
        self.connection = sqlite3.connect(filename, detect_types=sqlite3.PARSE_DECLTYPES)
        self.connection.text_factory = str
        self.connection.row_factory = sqlite3.Row

        self.cursor = self.connection.cursor()

        self.create_tables()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def commit(self):
        self.connection.commit()

    def close(self):
        self.connection.close()


class FragmentsDb(SqliteDb):
    """Fragments database"""

    def create_tables(self):
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS fragments (
            frag_id TEXT PRIMARY KEY,
            frag_nr INT,
            pdb_code TEXT,
            het_code TEXT,
            atomCodes TEXT,
            hashcode TEXT,
            ligID TEXT,
            numRgroups INT,
            ContactReses TEXT,
            SCOP_family TEXT,
            SCOP_fold TEXT,
            SCOP_species TEXT,
            SCOP_super TEXT
        )''')
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS molecules (
            frag_id TEXT PRIMARY KEY,
            smiles TEXT,
            molfile molblockgz
        )''')

    def add_molecules(self, mols):
        with FastInserter(self.cursor):
            for mol in mols:
                self.add_molecule(mol)

    def add_fragments_from_shelve(self, myshelve):
        with FastInserter(self.cursor):
            for k, v in myshelve.iteritems():
                self.add_fragment_from_shelve(k, v)

    def add_molecule(self, mol):
        sql = '''INSERT OR REPLACE INTO molecules (frag_id, smiles, molfile) VALUES (?, ?, ?)'''

        if mol is None:
            return

        self.cursor.execute(sql, (
            mol.GetProp('_Name'),
            MolToSmiles(mol),
            mol,
        ))

        self.connection.commit()

    def add_fragment_from_shelve(self, frag_id, fragment):
        sql = '''INSERT OR REPLACE INTO fragments (
            frag_id,
            pdb_code,
            het_code,
            frag_nr,
            atomCodes,
            hashcode,
            ligID,
            numRgroups,
            ContactReses,
            SCOP_family,
            SCOP_fold,
            SCOP_species,
            SCOP_super
        ) VALUES (
            :frag_id,
            :pdb_code,
            :het_code,
            :frag_nr,
            :atomCodes,
            :hashcode,
            :ligID,
            :numRgroups,
            :ContactReses,
            :SCOP_family,
            :SCOP_fold,
            :SCOP_species,
            :SCOP_super
        )'''

        splitted_frag_id = frag_id.split('-')
        try:
            frag_nr = int(splitted_frag_id[2].replace('frag', ''))
        except IndexError:
            frag_nr = None
        except ValueError as e:
            logging.warn('Weird id {}, skipping'.format(frag_id))
            return

        row = {
            'frag_id': frag_id.replace('-', '_'),
            'pdb_code': splitted_frag_id[0],
            'het_code': splitted_frag_id[1],
            'frag_nr': frag_nr,
            'hashcode': None,
            'atomCodes': None,
            'ligID': None,
            'numRgroups': None,
            'ContactReses': None,
            'SCOP_family': None,
            'SCOP_fold': None,
            'SCOP_species': None,
            'SCOP_super': None
        }
        for k, v in fragment.iteritems():
            row[k] = v

        if row['numRgroups'] is not None:
            row['numRgroups'] = int(row['numRgroups'])

        self.cursor.execute(sql, row)

    def __getitem__(self, key):
        sql = '''SELECT m.rowid, * FROM fragments LEFT JOIN molecules m USING (frag_id) WHERE frag_id=?'''
        self.cursor.execute(sql, (key,))
        row = self.cursor.fetchone()

        if row is None:
            raise KeyError("'{}' not found".format(key))

        fragment = {}
        for idx, v in enumerate(row.keys()):
            fragment[v] = row[idx]
        return fragment

    def _to_dict(self, sql):
        lookup = {}
        for k, v in self.cursor.execute(sql):
            lookup[k] = v
        return lookup

    def id2label(self):
        sql = '''SELECT rowid, frag_id FROM molecules'''
        return self._to_dict(sql)

    def label2id(self):
        sql = '''SELECT frag_id, rowid FROM molecules'''
        return self._to_dict(sql)


class FingerprintsDb(SqliteDb):

    def create_tables(self):
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS bitsets (
            frag_id TEXT PRIMARY KEY,
            bitset intbitset
        )''')
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS attributes (
            key TEXT PRIMARY KEY,
            value TEXT
        )''')

    def as_dict(self):
        return IntbitsetDict(self)


class IntbitsetDict(MutableMapping):
    """
    Dictionary of intbitset with sqlite3 backend.

    Convert dbm to sqlite:
    ```
    from kripodb.dbm import IntbitsetDictDbm
    from kripodb.db import FragmentsDb
    dbm = IntbitsetDictDbm('data/fingerprint12.fp.db')
    frags = FragmentsDb('data/fragments12.db')
    bitsets = frags.bitsets()
    bitsets.update(dbm)
    ```

    """

    def __init__(self, db, number_of_bits=None):
        self.db = db
        self.cursor = db.cursor
        if number_of_bits is not None:
            self.number_of_bits = number_of_bits

    def __delitem__(self, key):
        sql = '''DELETE FROM bitsets WHERE frag_id=?'''

        self.cursor.execute(sql, (key,))
        self.db.commit()

    def __iter__(self):
        sql = '''SELECT frag_id FROM bitsets'''
        for row in self.cursor.execute(sql):
            yield row[0]

    def __getitem__(self, key):
        sql = '''SELECT bitset FROM bitsets WHERE frag_id=?'''
        self.cursor.execute(sql, (key,))
        row = self.cursor.fetchone()

        if row is None:
            raise KeyError("'{}' not found".format(key))

        return row[0]

    def __setitem__(self, key, value):
        sql = '''INSERT OR REPLACE INTO bitsets (frag_id, bitset) VALUES (?, ?)'''

        self.cursor.execute(sql, (key, value))

        self.db.commit()

    def __len__(self):
        self.cursor.execute('SELECT count(*) FROM bitsets')
        row = self.cursor.fetchone()
        return row[0]

    def iteritems(self):
        sql = '''SELECT frag_id, bitset FROM bitsets'''
        for row in self.cursor.execute(sql):
            yield row

    def itervalues(self):
        sql = '''SELECT bitset FROM bitsets'''
        for row in self.cursor.execute(sql):
            yield row[0]

    def __contains__(self, key):
        sql = '''SELECT count(*) FROM bitsets WHERE frag_id=?'''
        self.cursor.execute(sql, (key,))
        row = self.cursor.fetchone()
        return row[0] == 1

    def update(*args, **kwds):
        self = args[0]

        with FastInserter(self.cursor):
            MutableMapping.update(*args, **kwds)
            # make table and index stored contiguously
            self.cursor.execute('VACUUM')

    def iteritems_startswith(self, prefix):
        """item iterator over keys with prefix

        :param prefix: Prefix of key
        :return:
        """
        sql = '''SELECT frag_id, bitset FROM bitsets WHERE frag_id LIKE ?'''
        for row in self.cursor.execute(sql, (prefix + '%',)):
            yield row

    @property
    def number_of_bits(self):
        """
        Number of bits the bitsets consist of
        :return: int
        """
        self.cursor.execute('SELECT value FROM attributes WHERE key=?', (ATTR_NUMBER_OF_BITS,))
        row = self.cursor.fetchone()
        if row is None:
            return None
        return int(row[0])

    @number_of_bits.setter
    def number_of_bits(self, value):
        sql = 'INSERT OR REPLACE INTO attributes (key, value) VALUES (?, ?)'
        self.cursor.execute(sql, (ATTR_NUMBER_OF_BITS, str(value)))
        self.db.commit()

    @number_of_bits.deleter
    def number_of_bits(self):
        sql = 'DELETE FROM attributes WHERE key=?'
        self.cursor.execute(sql, (ATTR_NUMBER_OF_BITS,))
        self.db.commit()
