import argparse

from ..pharmacophores import PharmacophoresDb


def dir2db_run(startdir, pharmacophoresdb, nrrows):
    with PharmacophoresDb(pharmacophoresdb, 'a', expectedrows=nrrows) as db:
        db.add_dir(startdir)


def add_sc(sc):
    parser = sc.add_parser('add', help='Add pharmacophores from directory to database')
    parser.add_argument('startdir', help='Directory to start finding *.pphores.sd.gz and *.pphores.txt files in')
    parser.add_argument('pharmacophoresdb',
                        default='pharmacophores.db',
                        help='Name of pharmacophore db file (default: %(default)s)')
    parser.add_argument('--nrrows',
                        type=int,
                        default=2 ** 16,
                        help='''Number of expected pharmacophores,
                        only used when database is created
                        (default: %(default)s)''')
    parser.set_defaults(func=dir2db_run)


def make_pharmacophores_parser(subparsers):
    """Creates a parser for pharmacophores sub commands
    
    Args:
        subparsers (argparse.ArgumentParser): Parser to which to add sub commands to 

    """
    sc = subparsers.add_parser('pharmacophores', help='Pharmacophores').add_subparsers()
    add_sc(sc)
    get_sc(sc)


def get_run(pharmacophoresdb, query):
    with PharmacophoresDb(pharmacophoresdb) as db:
        print(db[query])


def get_sc(sc):
    parser = sc.add_parser('get', help='Retrieve pharmacophore of a fragment')
    parser.add_argument('pharmacophoresdb',
                        default='pharmacophores.db',
                        help='Name of pharmacophore db file (default: %(default)s)')
    parser.add_argument('query', type=str, help='Query fragment identifier')
    parser.set_defaults(func=get_run)
