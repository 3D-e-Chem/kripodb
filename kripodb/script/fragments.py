import argparse
import logging
import shelve

from rdkit.Chem.rdmolfiles import SDMolSupplier

from ..db import FragmentsDb
from ..hdf5 import SimilarityMatrix
from ..pdb import PdbReport


def make_fragments_parser(subparsers):
    """Creates a parser for fragments sub commands

    Args:
        subparsers (argparse.ArgumentParser): Parser to which to add sub commands to
    """
    sc = subparsers.add_parser('fragments', help='Fragments').add_subparsers()
    shelve2fragmentsdb_sc(sc)
    sdf2fragmentsdb_sc(sc)
    pdb2fragmentsdb_sc(sc)
    fragmentsdb_filter_sc(sc)


def shelve2fragmentsdb_sc(subparsers):
    sc = subparsers.add_parser('shelve', help='Add fragments from shelve to sqlite')
    sc.add_argument('--skipdups', action='store_true', help='Skip duplicates, instead of dieing one first duplicate')
    sc.add_argument('shelvefn', type=str)
    sc.add_argument('fragmentsdb',
                    default='fragments.db',
                    help='Name of fragments db file (default: %(default)s)')
    sc.set_defaults(func=shelve2fragmentsdb_run)


def shelve2fragmentsdb_run(shelvefn, fragmentsdb, skipdups):
    myshelve = shelve.open(shelvefn, 'r')
    frags = FragmentsDb(fragmentsdb)
    frags.add_fragments_from_shelve(myshelve, skipdups)


def sdf2fragmentsdb_sc(subparsers):
    sc = subparsers.add_parser('sdf', help='Add fragments sdf to sqlite')
    sc.add_argument('sdffns', help='SDF filename', nargs='+')
    sc.add_argument('fragmentsdb',
                    default='fragments.db',
                    help='Name of fragments db file (default: %(default)s)')

    sc.set_defaults(func=sdf2fragmentsdb_run)


def sdf2fragmentsdb_run(sdffns, fragmentsdb):
    frags = FragmentsDb(fragmentsdb)
    for sdffn in sdffns:
        logging.warning('Parsing {}'.format(sdffn))
        suppl = SDMolSupplier(sdffn)
        frags.add_molecules(suppl)


def pdb2fragmentsdb_sc(subparsers):
    sc = subparsers.add_parser('pdb', help='Add pdb metadata from RCSB PDB website to fragment sqlite db')
    sc.add_argument('fragmentsdb',
                    default='fragments.db',
                    help='Name of fragments db file (default: %(default)s)')

    sc.set_defaults(func=pdb2fragmentsdb_run)


def pdb2fragmentsdb_run(fragmentsdb):
    pdb_report = PdbReport()
    pdbs = pdb_report.fetch()
    frags = FragmentsDb(fragmentsdb)
    frags.add_pdbs(pdbs)


def fragmentsdb_filter_sc(subparsers):
    sc = subparsers.add_parser('filter', help='Filter fragments database')
    sc.add_argument('input', type=str,
                    help='Name of fragments db input file')
    sc.add_argument('output', type=str,
                    help='Name of fragments db output file, will overwrite file if it exists')
    sc.add_argument('--pdbs', type=argparse.FileType('r'),
                    help='Keep fragments from any of the supplied pdb codes, one pdb code per line, use - for stdin')
    sc.add_argument('--matrix', type=str, help='Keep fragments which are in similarity matrix file')
    sc.set_defaults(func=fragmentsdb_filter)


def fragmentsdb_filter(input, output, pdbs, matrix):
    if matrix:
        fragmentsdb_filter_matrix(input, output, matrix)
    else:
        fragmentsdb_filter_pdbs(input, output, pdbs)


def fragmentsdb_filter_matrix(input, output, matrix):
    output_db = FragmentsDb(output)

    # mount input into output db
    print('Reading: ' + input)
    output_db.cursor.execute('ATTACH DATABASE ? AS orig', (input,))

    # create temp table with pdbs
    output_db.cursor.execute('CREATE TEMPORARY TABLE filter (frag_id TEXT PRIMARY KEY)')
    sql = 'INSERT OR REPLACE INTO filter (frag_id) VALUES (?)'
    print('Matrix labels')
    simmatrix = SimilarityMatrix(matrix)
    for frag_id in simmatrix.labels.label2ids().keys():
        output_db.cursor.execute(sql, (frag_id,))
    simmatrix.close()

    # insert select
    output_db.cursor.execute('INSERT INTO fragments SELECT * FROM orig.fragments JOIN filter USING (frag_id)')
    output_db.cursor.execute(
        'INSERT INTO pdbs SELECT * FROM orig.pdbs WHERE pdb_code IN (SELECT pdb_code FROM fragments)')
    output_db.cursor.execute(
        'INSERT INTO molecules SELECT * FROM orig.molecules WHERE frag_ID IN (SELECT frag_id FROM fragments)')

    # drop temp table with pdbs
    output_db.cursor.execute('DROP TABLE filter')

    # vacuum
    output_db.cursor.execute('VACUUM')

    print('Wrote: ' + output)


def fragmentsdb_filter_pdbs(input, output, pdbs):
    output_db = FragmentsDb(output)

    # mount input into output db
    print('Reading: ' + input)
    output_db.cursor.execute('ATTACH DATABASE ? AS orig', (input,))

    # create temp table with pdbs
    output_db.cursor.execute('CREATE TEMPORARY TABLE filter (pdb_code TEXT PRIMARY KEY)')
    sql = 'INSERT OR REPLACE INTO filter (pdb_code) VALUES (?)'
    for pdb in pdbs:
        output_db.cursor.execute(sql, (pdb.rstrip().lower(),))

    # insert select
    output_db.cursor.execute('INSERT INTO pdbs SELECT * FROM orig.pdbs JOIN filter USING (pdb_code)')
    output_db.cursor.execute('INSERT INTO fragments SELECT * FROM orig.fragments JOIN filter USING (pdb_code)')
    output_db.cursor.execute('INSERT INTO molecules SELECT * FROM orig.molecules WHERE frag_ID IN (SELECT frag_id FROM fragments)')

    # drop temp table with pdbs
    output_db.cursor.execute('DROP TABLE filter')

    # vacuum
    output_db.cursor.execute('VACUUM')

    print('Wrote: ' + output)


