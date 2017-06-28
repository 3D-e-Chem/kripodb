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
"""Fragments and fingerprints sqlite based data storage.

Registers `intbitset` and `molblockgz` data types in sqlite.
"""

from __future__ import absolute_import
from collections import MutableMapping
import sqlite3
import logging
import zlib
import re

from intbitset import intbitset
from rdkit.Chem import MolToMolBlock, MolFromMolBlock, MolToSmiles
from rdkit.Chem.rdchem import Mol
import six

ATTR_NUMBER_OF_BITS = 'number_of_bits'


def adapt_intbitset(ibs):
    """Convert intbitset to fast dumped intbitset

    Args:
        ibs (intbitset): bitset

    Examples:
        Serialize intbitset

        >>> adapt_intbitset(intbitset([1, 2, 3, 4]))
        'x\x9c\x93c@\x05\x00\x01\xf0\x00\x1f'

    Returns:
        str: Fast dumped intbitset
    """
    return ibs.fastdump()


def convert_intbitset(s):
    """Convert fast dumped intbitset to intbitset

    Args:
        s (str): Fast dumped intbitset

    Examples:
        Deserialize intbitset

        >>> ibs = convert_intbitset('x\x9c\x93c@\x05\x00\x01\xf0\x00\x1f')
        intbitset([1, 2, 3, 4])

    Returns:
        intbitset: bitset
    """
    ibs = intbitset()
    ibs.fastload(s)
    return ibs


def adapt_molblockgz(mol):
    """Convert RDKit molecule to compressed molblock

    Args:
        mol (rdkit.Chem.Mol): molecule

    Returns:
        str: Compressed molblock
    """
    molblock = MolToMolBlock(mol).encode()
    return zlib.compress(molblock)


def convert_molblockgz(molgz):
    """Convert compressed molblock to RDKit molecule

    Args:
        molgz: (str) zlib compressed molblock

    Returns:
        rdkit.Chem.Mol: molecule
    """
    return MolFromMolBlock(zlib.decompress(molgz))


sqlite3.register_adapter(intbitset, adapt_intbitset)
sqlite3.register_converter('intbitset', convert_intbitset)
sqlite3.register_adapter(Mol, adapt_molblockgz)
sqlite3.register_converter('molblockgz', convert_molblockgz)


class FastInserter(object):
    """Use with to make inserting faster, but less safe

    By setting journal mode to WAL and turn synchronous off.

    Args:
        cursor (sqlite3.Cursor): Sqlite cursor

    Examples:

        >>> with FastInserter(cursor):
                cursor.executemany('INSERT INTO table VALUES (?), rows))

    """

    def __init__(self, cursor):
        self.cursor = cursor

    def __enter__(self):
        # increase insert speed, this is less safe
        self.cursor.connection.commit()
        self.cursor.execute('PRAGMA journal_mode=WAL')
        self.cursor.execute('PRAGMA synchronous=OFF')

    def __exit__(self, exc_type, exc_val, exc_tb):
        # switch back to default journal, so db file can be read-only and is safe again
        self.cursor.connection.commit()
        self.cursor.execute('PRAGMA journal_mode=DELETE')
        self.cursor.execute('PRAGMA synchronous=FULL')


class SqliteDb(object):
    """Wrapper around a sqlite database connection

    Database is created if it does not exist.

    Args:
        filename (str):  Sqlite filename

    Attributes:
        connection (sqlite3.Connection): Sqlite connection
        cursor (sqlite3.Cursor): Sqlite cursor
    """

    def __init__(self, filename):
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
        """Commit pending changes"""
        self.connection.commit()

    def close(self):
        """Close database"""
        self.connection.close()


def _row2fragment(row):
    fragment = {}
    for idx, v in enumerate(row.keys()):
        fragment[v] = row[idx]
    return fragment


class FragmentsDb(SqliteDb):
    """Fragments database"""
    select_sql = '''SELECT f.rowid, * FROM fragments f
                    JOIN pdbs USING (pdb_code, prot_chain)
                    LEFT JOIN molecules USING (frag_id)'''

    def create_tables(self):
        """Create tables if they don't exist"""
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS fragments (
            frag_id TEXT PRIMARY KEY,
            frag_nr INT NOT NULL,
            pdb_code TEXT NOT NULL,
            prot_chain TEXT NOT NULL,
            het_chain TEXT NOT NULL,
            het_code TEXT NOT NULL,
            het_seq_nr INT,
            atom_codes TEXT,
            hash_code TEXT,
            nr_r_groups INT
        )''')
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS molecules (
            frag_id TEXT PRIMARY KEY,
            smiles TEXT,
            mol molblockgz
        )''')
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS pdbs (
            pdb_code TEXT NOT NULL,
            prot_chain TEXT NOT NULL,
            pdb_title TEXT,
            prot_name TEXT,
            uniprot_acc TEXT,
            uniprot_name TEXT,
            ec_number TEXT,
            PRIMARY KEY (pdb_code, prot_chain)
        )''')

    def add_molecules(self, mols):
        """Adds molecules to to molecules table.

        Args:
            mols (list[rdkit.Chem.Mol]): List of molecules
        """
        with FastInserter(self.cursor):
            for mol in mols:
                self.add_molecule(mol)

    def add_pdbs(self, pdbs):
        """Adds pdb meta data to to pdbs table.

        Args:
            pdbs (Iterable[Dict]): List of pdb meta data
        """

        rows = self.cursor.execute('SELECT pdb_code || prot_chain FROM fragments')
        pdbs_in_fragments = frozenset([r[0] for r in rows])

        with FastInserter(self.cursor):
            for pdb in pdbs:
                if pdb['structureId'].lower() + pdb['chainId'] in pdbs_in_fragments:
                    self.add_pdb(pdb)

    def add_fragments_from_shelve(self, myshelve, skipdups=False):
        """Adds fragments from shelve to fragments table.

        Also creates index on pdb_code column.

        Args:
            myshelve (Dict[Fragment]): Dictionary with fragment identifier as key and fragment as value.
            skipdups (bool): Skip duplicates, instead of dieing one first duplicate

        """
        with FastInserter(self.cursor):
            for k, v in six.iteritems(myshelve):
                self.add_fragment_from_shelve(k, v, skipdups)

        self.cursor.execute('CREATE INDEX IF NOT EXISTS fragments_pdb_code_i ON fragments (pdb_code)')

    def add_molecule(self, mol):
        """Adds molecule to molecules table

        Args:
            mol (rdkit.Chem.AllChem.Mol): the rdkit molecule

        """
        sql = '''INSERT OR REPLACE INTO molecules (frag_id, smiles, mol) VALUES (?, ?, ?)'''

        if mol is None:
            logging.warning('Empty molecule, skipping')
            return

        self.cursor.execute(sql, (
            mol.GetProp('_Name'),
            MolToSmiles(mol),
            mol,
        ))

        self.connection.commit()

    def add_fragment_from_shelve(self, frag_id, fragment, skipdups=False):
        sql = '''INSERT INTO fragments (
            frag_id,
            pdb_code,
            prot_chain,
            het_code,
            frag_nr,
            atom_codes,
            hash_code,
            het_chain,
            het_seq_nr,
            nr_r_groups
        ) VALUES (
            :frag_id,
            :pdb_code,
            :prot_chain,
            :het_code,
            :frag_nr,
            :atom_codes,
            :hash_code,
            :het_chain,
            :het_seq_nr,
            :nr_r_groups
        )'''

        splitted_frag_id = frag_id.split('-')
        if len(splitted_frag_id) != 3:
            logging.warning('Weird id {}, skipping'.format(frag_id))
            return

        try:
            frag_nr = int(splitted_frag_id[2].replace('frag', ''))
        except ValueError:
            logging.warning('Weird id {}, skipping'.format(frag_id))
            return

        lig_id = fragment['ligID'].split('-')
        het_seq_nr = int(re.sub('[A-Z]$', '', lig_id[3]))

        frag_id = frag_id.replace('-', '_')
        row = {
            'frag_id': frag_id,
            'pdb_code': splitted_frag_id[0],
            'prot_chain': lig_id[1],
            'het_code': splitted_frag_id[1],
            'het_seq_nr': het_seq_nr,
            'het_chain': lig_id[4],
            'frag_nr': frag_nr,
            'hash_code': fragment['hashcode'],
            'atom_codes': fragment['atomCodes'],
            'nr_r_groups': int(fragment['numRgroups']),
        }

        try:
            self.cursor.execute(sql, row)
        except sqlite3.IntegrityError as e:
            logging.warning('Duplicate ID: {}, skipping'.format(frag_id))
            if not skipdups:
                raise e

    def add_pdb(self, pdb):
        sql = '''INSERT OR REPLACE INTO pdbs (
            pdb_code,
            prot_chain,
            pdb_title,
            prot_name,
            uniprot_acc,
            uniprot_name,
            ec_number
        ) VALUES (
            :pdb_code,
            :prot_chain,
            :pdb_title,
            :prot_name,
            :uniprot_acc,
            :uniprot_name,
            :ec_number
        )'''
        pdb2col = {
            'structureId': 'pdb_code',
            'chainId': 'prot_chain',
            'structureTitle': 'pdb_title',
            'compound': 'prot_name',
            'uniprotAcc': 'uniprot_acc',
            'uniprotRecommendedName': 'uniprot_name',
            'ecNo': 'ec_number',
        }
        row = {pdb2col[k]: v for k, v in six.iteritems(pdb)}
        row['pdb_code'] = row['pdb_code'].lower()
        self.cursor.execute(sql, row)

    def __getitem__(self, key):
        """Retrieve fragment based on it's identifier.

        Args:
            key (str): Fragment identifier

        Returns:
            Fragment

        """
        sql = self.select_sql + 'WHERE frag_id=?'
        self.cursor.execute(sql, (key,))
        row = self.cursor.fetchone()

        if row is None:
            raise KeyError(key)

        return _row2fragment(row)

    def by_pdb_code(self, pdb_code):
        """Retrieve fragments which are part of a PDB structure.

        Args:
            pdb_code (str): PDB code

        Returns:
            List[Fragment]: List of fragments

        Raises:
            LookupError: When pdb_code could not be found

        """
        fragments = []
        sql = self.select_sql + 'WHERE pdb_code=? ORDER BY frag_id'
        for row in self.cursor.execute(sql, (pdb_code,)):
            fragments.append(_row2fragment(row))

        if len(fragments) == 0:
            raise LookupError(pdb_code)

        return fragments

    def id2label(self):
        """Lookup table of fragments from an number to a label.

        Returns:
            SqliteDict

        """
        return SqliteDict(self.connection, 'fragments', 'rowid', 'frag_id')

    def label2id(self):
        """Lookup table of fragments from an label to a number.

        Returns:
            SqliteDict

        """
        return SqliteDict(self.connection, 'fragments', 'frag_id', 'rowid')

    def __len__(self):
        self.cursor.execute('SELECT count(*) FROM fragments')
        row = self.cursor.fetchone()
        return row[0]


class FingerprintsDb(SqliteDb):
    """Fingerprints database"""

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
        """Returns a dict-like object to query and alter fingerprints db

        Args:
            number_of_bits (Optional[int]): Number of bits that all fingerprints have

        Returns:
            IntbitsetDict
        """
        return IntbitsetDict(self, number_of_bits)


class SqliteDict(MutableMapping):
    """Dict-like object of 2 columns of a sqlite table.

    Can be used to query and alter the table.

    Args:
        connection (sqlite3.Connection): Sqlite connection
        table_name (str): Table name
        key_column (str): Column name used as key
        value_column (str): Column name used as value

    Attributes:
        connection (sqlite3.Connection): Sqlite connection
        cursor (sqlite3.Cursor): Sqlite cursor

    """
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
        sql = self.sqls['setitem']
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

    def items(self):
        sql = self.sqls['iteritems']
        for row in self.cursor.execute(sql):
            yield row

    def values(self):
        sql = self.sqls['itervalues']
        for row in self.cursor.execute(sql):
            yield row[0]

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

        Args:
            prefix (str): Prefix of key

        Examples:
           All items with key starting with letter 'a' are returned.

            >>> for frag_id, fragment in fragments.iteritems_startswith('a'):
                    # do something with frag_id and fragment

        Returns:
            List[Tuple[key, value]]

        """
        sql = self.sqls['iteritems_startswith']
        for row in self.cursor.execute(sql, (prefix + '%',)):
            yield row

    def materialize(self):
        """Fetches all kev/value pairs from the sqlite database.

        Useful when dictionary is iterated multiple times and the cost of fetching is to high.

        Returns:
            Dict: Dictionary with all kev/value pairs
        """
        return {k: v for k, v in six.iteritems(self)}


class IntbitsetDict(SqliteDict):
    """Dictionary of intbitset with sqlite3 backend.

    Args:
        db (FingerprintsDb): Fingerprints db
        number_of_bits (int): Number of bits

    Attributes:
        number_of_bits (int): Number of bits the bitsets consist of

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
