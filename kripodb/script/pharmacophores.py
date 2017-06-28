import argparse
import json

from kripodb.db import FragmentsDb
from ..pharmacophores import PharmacophoresDb
from ..pharmacophores.align import Aligner


def dir2db_run(startdir, pharmacophoresdb, nrrows):
    with PharmacophoresDb(pharmacophoresdb, 'a', expectedrows=nrrows) as db:
        db.add_dir(startdir)


def add_sc(sc):
    parser = sc.add_parser('add', help='Add pharmacophores from directory to database')
    parser.add_argument('startdir', help='Directory to start finding *.pphores.sd.gz and *.pphores.txt files in')
    parser.add_argument('pharmacophoresdb', help='Name of pharmacophore db file')
    parser.add_argument('--nrrows',
                        type=int,
                        default=2 ** 16,
                        help='''Number of expected pharmacophores,
                        only used when database is created
                        (default: %(default)s)''')
    parser.set_defaults(func=dir2db_run)


def get_run(pharmacophoresdb, query, output):
    with PharmacophoresDb(pharmacophoresdb) as db:
       db.write_phar(output, query)


def get_sc(sc):
    parser = sc.add_parser('get', help='Retrieve pharmacophore of a fragment')
    parser.add_argument('pharmacophoresdb', help='Name of pharmacophore db file')
    parser.add_argument('query', type=str, help='Query fragment identifier')
    parser.add_argument('--output', type=argparse.FileType('w'), default='-')
    parser.set_defaults(func=get_run)


def filter_run(inputfn, fragmentsdb, outputfn):
    frags = FragmentsDb(fragmentsdb)
    fragids2keep = set([f.encode() for f in frags.id2label().values()])
    with PharmacophoresDb(inputfn) as dbin:
        expectedrows = len(dbin.points)
        with PharmacophoresDb(outputfn, 'w', expectedrows=expectedrows) as dbout:
            col_names = [colName for colName in dbin.points.table.colpathnames]
            rowout = dbout.points.table.row
            for rowin in dbin.points.table.iterrows():
                if rowin['frag_id'] in fragids2keep:
                    for col_name in col_names:
                        rowout[col_name] = rowin[col_name]
                    rowout.append()
            dbout.points.table.flush()


def filter_sc(sc):
    parser = sc.add_parser('filter', help='Filter pharmacophores')
    parser.add_argument('inputfn', help='Name of input pharmacophore db file')
    parser.add_argument('--fragmentsdb',
                        default='fragments.db',
                        help='Name of fragments db file, fragments present in db are passed '
                             '(default: %(default)s)')
    parser.add_argument('outputfn', help='Name of output pharmacophore db file')
    parser.set_defaults(func=filter_run)


def align_sc(sc):
    parser = sc.add_parser('align', help='Align pharmacophore to another pharmacophore')
    parser.add_argument('pharmacophoresdb', help='Name of pharmacophore db file')
    parser.add_argument('reference', type=str, help='Reference fragment identifier')
    parser.add_argument('probe', type=str, help='Probe fragment identifier')
    parser.add_argument('--break_num_cliques',
                        type=int,
                        default=3000,
                        help='''Break when {break_num_cliques} cliques found
                        (default: %(default)s)''')
    parser.add_argument('--cutoff',
                        type=float,
                        default=1.0,
                        help='Tolerance threshold for considering two distances to be equivalent (default: %(default)s)'
                        )
    parser.set_defaults(func=align_run)

def align_run(pharmacophoresdb, reference, probe, cutoff, break_num_cliques):
    with PharmacophoresDb(pharmacophoresdb) as db:
        reference_points = db[reference]
        probe_points = db[probe]
        aligner = Aligner(reference_points, probe_points)
        transform = aligner.transformation(cutoff, break_num_cliques)
        print(json.dumps(transform))

def make_pharmacophores_parser(subparsers):
    """Creates a parser for pharmacophores sub commands

    Args:
        subparsers (argparse.ArgumentParser): Parser to which sub commands are added

    """
    sc = subparsers.add_parser('pharmacophores', help='Pharmacophores').add_subparsers()
    add_sc(sc)
    get_sc(sc)
    filter_sc(sc)
    align_sc(sc)
