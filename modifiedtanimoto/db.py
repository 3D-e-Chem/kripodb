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


class FastInserter(object):
    def __init__(self, cursor):
        self.cursor = cursor

    def __enter__(self):
        # increase insert speed, this is less safe
        self.cursor.execute('PRAGMA journal_mode=WAL')
        self.cursor.execute('PRAGMA synchronous=OFF')

    def __exit__(self, exc_type, exc_val, exc_tb):
        # switch back to default journal, so db file can be read-only and is safe again
        self.cursor.execute('PRAGMA journal_mode=DELETE')
        self.cursor.execute('PRAGMA synchronous=FULL')


class FragmentsDb(object):
    """Fragments database"""
    def __init__(self, filename):
        """Connect to fragment db.

        Database is created if it does not exist.

        :param filename: Sqlite filename
        :return:
        """
        self.filename = filename

        sqlite3.register_adapter(intbitset, adapt_intbitset)
        sqlite3.register_converter('intbitset', convert_intbitset)
        sqlite3.register_adapter(Mol, adapt_molblockgz)
        sqlite3.register_converter('molblockgz', convert_molblockgz)
        self.connection = sqlite3.connect(filename, detect_types=sqlite3.PARSE_DECLTYPES)
        self.connection.text_factory = str

        self.cursor = self.connection.cursor()

        self.create_tables()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def create_tables(self):
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS fragments (
            frag_id TEXT PRIMARY KEY,
            pdb_code TEXT,
            het_code TEXT,
            atom_codes TEXT
            hash_code TEXT,
            ligand_id TEXT,
            nr_r_groups INT,
            smile TEXT,
            molfile molblockgz
        )''');
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS sites (
            site_id TEXT PRIMARY KEY,
            pdb_code TEXT,
            het_code TEXT,
            contact_residues TEXT,
            scop_family TEXT,
            scop_fold TEXT,
            scop_species TEXT
            ligand_id TEXT
        )''')
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS bitsets (
            frag_id TEXT PRIMARY KEY,
            bitset intbitset
        )''')
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS attributes (
            key TEXT PRIMARY KEY,
            value TEXT
        )''')

    def commit(self):
        self.connection.commit()

    def close(self):
        self.connection.close()

    def bitsets(self):
        return IntbitsetDict(self)

    def add_molecules(self, mols):
        with FastInserter(self.cursor):
            for mol in mols:
                self.add_molecule(mol)

    def add_fragments_from_shelve(self, myshelve):
        with FastInserter(self.cursor):
            for k, v in myshelve.iteritems():
                self.add_fragment_from_shelve(k, v)

    def add_molecule(self, mol):
        sql = '''UPDATE fragments SET smile=?, molfile=? WHERE frag_id=?'''

        self.cursor.execute(sql, (
            MolToSmiles(mol),
            mol,
            mol.GetProp('_Name')
        ))

        self.db.commit()

    def add_fragment_from_shelve(self, frag_id, fragment):
        sql = '''INSERT OR REPLACE INTO bitsets (frag_id, bitset) VALUES (?, ?)'''


class IntbitsetDict(MutableMapping):
    """
    Dictionary of intbitset with sqlite3 backend.

    Convert dbm to sqlite:
    ```
    from modifiedtanimoto.dbm import IntbitsetDictDbm
    from modifiedtanimoto.db import FragmentsDb
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
