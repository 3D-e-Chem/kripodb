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
        # sqlite3 defaults to unicode as text_factory, unicode can't be used for byte string
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
            numRgroups INT
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

        self.cursor.execute('CREATE INDEX IF NOT EXISTS fragments_pdb_code_i ON fragments (pdb_code)')

    def add_molecule(self, mol):
        sql = '''INSERT OR REPLACE INTO molecules (frag_id, smiles, molfile) VALUES (?, ?, ?)'''

        if mol is None:
            logging.warn('Empty molecule, skipping')
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
            numRgroups
        ) VALUES (
            :frag_id,
            :pdb_code,
            :het_code,
            :frag_nr,
            :atomCodes,
            :hashcode,
            :ligID,
            :numRgroups
        )'''

        splitted_frag_id = frag_id.split('-')
        if len(splitted_frag_id) != 3:
            logging.warn('Weird id {}, skipping'.format(frag_id))
            return

        try:
            frag_nr = int(splitted_frag_id[2].replace('frag', ''))
        except ValueError:
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
            'numRgroups': None
        }
        for k, v in fragment.iteritems():
            row[k] = v

        row['numRgroups'] = int(row['numRgroups'])

        self.cursor.execute(sql, row)

    def __getitem__(self, key):
        sql = '''SELECT f.rowid, * FROM fragments f LEFT JOIN molecules USING (frag_id) WHERE frag_id=?'''
        self.cursor.execute(sql, (key,))
        row = self.cursor.fetchone()

        if row is None:
            raise KeyError("'{}' not found".format(key))

        return self._row2fragment(row)

    def _row2fragment(self, row):
        fragment = {}
        for idx, v in enumerate(row.keys()):
            fragment[v] = row[idx]
        return fragment

    def by_pdb_code(self, pdb_code):
        fragments = []
        sql = '''SELECT m.rowid, * FROM fragments JOIN molecules m USING (frag_id) WHERE pdb_code=? ORDER BY frag_id'''
        for row in self.cursor.execute(sql, (pdb_code,)):
            fragments.append(self._row2fragment(row))

        return fragments

    def id2label(self):
        return SqliteDict(self.connection, 'fragments', 'rowid', 'frag_id')

    def label2id(self):
        return SqliteDict(self.connection, 'fragments', 'frag_id', 'rowid')

    def __len__(self):
        self.cursor.execute('SELECT count(*) FROM fragments')
        row = self.cursor.fetchone()
        return row[0]


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

    def as_dict(self, number_of_bits=None):
        return IntbitsetDict(self, number_of_bits)


class SqliteDict(MutableMapping):
    def __init__(self, connection, table_name, key_column, value_column):
        self.connection = connection
        self.cursor = connection.cursor()
        kwargs = {
            'key_column': key_column,
            'table_name': table_name,
            'value_column': value_column
        }
        self.sqls = {
            'iter': 'SELECT {key_column} FROM {table_name}'.format(**kwargs),
            'getitem': 'SELECT {value_column} FROM {table_name} WHERE {key_column}=?'.format(**kwargs),
            'delitem': 'DELETE FROM {table_name} WHERE {key_column}=?'.format(**kwargs),
            'setitem': '''INSERT OR REPLACE INTO {table_name}
                          ({key_column}, {value_column}) VALUES (?, ?)'''.format(**kwargs),
            'len': 'SELECT count(*) FROM {table_name}'.format(**kwargs),
            'iteritems': 'SELECT {key_column}, {value_column} FROM {table_name}'.format(**kwargs),
            'itervalues': 'SELECT {value_column} FROM {table_name}'.format(**kwargs),
            'contains': 'SELECT count(*) FROM {table_name} WHERE {key_column}=?'.format(**kwargs),
            'iteritems_startswith': '''SELECT {key_column}, {value_column} FROM {table_name}
                                    WHERE {key_column} LIKE ?'''.format(**kwargs),
        }

    def __iter__(self):
        sql = self.sqls['iter']
        for row in self.cursor.execute(sql):
            yield row[0]

    def __getitem__(self, key):
        sql = self.sqls['getitem']
        self.cursor.execute(sql, (key,))
        row = self.cursor.fetchone()

        if row is None:
            raise KeyError(key)

        return row[0]

    def __delitem__(self, key):
        sql = self.sqls['delitem']
        self.cursor.execute(sql, (key,))
        self.connection.commit()

    def __setitem__(self, key, value):
        sql = '''INSERT OR REPLACE INTO bitsets (frag_id, bitset) VALUES (?, ?)'''
        self.cursor.execute(sql, (key, value))
        self.connection.commit()

    def __len__(self):
        sql = self.sqls['len']
        self.cursor.execute(sql)
        row = self.cursor.fetchone()
        return row[0]

    def iteritems(self):
        sql = self.sqls['iteritems']
        for row in self.cursor.execute(sql):
            yield row

    def itervalues(self):
        sql = self.sqls['itervalues']
        for row in self.cursor.execute(sql):
            yield row[0]

    def __contains__(self, key):
        sql = self.sqls['contains']
        self.cursor.execute(sql, (key,))
        row = self.cursor.fetchone()
        return row[0] == 1

    def iteritems_startswith(self, prefix):
        """item iterator over keys with prefix

        :param prefix: Prefix of key
        :return:
        """
        sql = self.sqls['iteritems_startswith']
        for row in self.cursor.execute(sql, (prefix + '%',)):
            yield row

    def materialize(self):
        """Fetches all kev/value pairs from the sqlite database.

        Useful when dictionary is iterated multiple times and the cost of fetching is to high.

        Returns
        -------
        Dictionary with all kev/value pairs

        """
        return {k: v for k, v in self.iteritems()}


class IntbitsetDict(SqliteDict):
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
        super(IntbitsetDict, self).__init__(db.connection, 'bitsets', 'frag_id', 'bitset')
        if number_of_bits is not None:
            self.number_of_bits = number_of_bits

    def update(*args, **kwds):
        self = args[0]

        with FastInserter(self.cursor):
            MutableMapping.update(*args, **kwds)
            # make table and index stored contiguously
            self.cursor.execute('VACUUM')

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
        self.connection.commit()

    @number_of_bits.deleter
    def number_of_bits(self):
        sql = 'DELETE FROM attributes WHERE key=?'
        self.cursor.execute(sql, (ATTR_NUMBER_OF_BITS,))
        self.connection.commit()
