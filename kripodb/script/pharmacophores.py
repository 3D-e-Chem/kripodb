import argparse

from ..db import FragmentsDb
from ..pharmacophores import PharmacophoresDb, read_pphore_sdfile, as_phar


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
    parser.add_argument('--query', type=str, help='Query fragment identifier', default=None)
    parser.add_argument('--output', type=argparse.FileType('w'), default='-', help="Phar formatted text file")
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


def merge_sc(sc):
    parser = sc.add_parser('merge', help='Merge pharmacophore database files into new one')
    parser.add_argument('ins', nargs='+', help='Input pharmacophore database files')
    parser.add_argument('out', help='Output pharmacophore database file')
    parser.set_defaults(func=merge_pharmacophore_dbs)


def merge_pharmacophore_dbs(ins, out):
    nr_rows = 0
    for in_fn in ins:
        with PharmacophoresDb(in_fn) as in_db:
            nr_rows += len(in_db)

    with PharmacophoresDb(out, 'a', expectedrows=nr_rows) as db:

        for in_fn in ins:
            with PharmacophoresDb(in_fn) as in_db:
                db.append(in_db)


def phar2db_sc(sc):
    parser = sc.add_parser('import', help='Convert phar formatted file to pharmacophore database file')
    parser.add_argument('infile', type=argparse.FileType('r'), help='Input phar formatted file')
    parser.add_argument('outfile', help='Output pharmacophore database file')
    parser.add_argument('--nrrows',
                        type=int,
                        default=2 ** 16,
                        help='''Number of expected pharmacophores,
                            only used when database is created
                            (default: %(default)s)''')
    parser.set_defaults(func=phar2db)


def phar2db(infile, outfile, nrrows):
    with PharmacophoresDb(outfile, 'a', expectedrows=nrrows) as out_db:
        out_db.read_phar(infile)


def sd2phar_sc(sc):
    parser = sc.add_parser('sd2phar', help='Convert sd formatted pharmacophore file to phar formatted file')
    parser.add_argument('infile', type=argparse.FileType('rb'), help='Input sd formatted file')
    parser.add_argument('outfile', type=argparse.FileType('w'), help='Output phar formatted file')
    parser.add_argument('--frag_id', type=str, help='Fragment identifier', default='frag')
    parser.set_defaults(func=sd2phar)


def sd2phar(infile, outfile, frag_id):
    points = read_pphore_sdfile(infile)
    phar = as_phar(frag_id, points)
    outfile.write(phar)


def make_pharmacophores_parser(subparsers):
    """Creates a parser for pharmacophores sub commands

    Args:
        subparsers (argparse.ArgumentParser): Parser to which sub commands are added

    """
    sc = subparsers.add_parser('pharmacophores', help='Pharmacophores').add_subparsers()
    add_sc(sc)
    get_sc(sc)
    filter_sc(sc)
    merge_sc(sc)
    phar2db_sc(sc)
    sd2phar_sc(sc)
