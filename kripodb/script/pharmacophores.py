import argparse

from ..pharmacophores import PharmacophoresDb


def dir2db_run(startdir, pharmacophoresdb):
    db = PharmacophoresDb(pharmacophoresdb)
    db.add_dir(startdir)


def add_sc(sc):
    parser = sc.addparser('add', help='Add pharmacophores from directory to database')
    parser.add_argument('startdir', help='Directory to start finding *.pphores.sd.gz and *.pphores.txt files in')
    parser.add_argument('pharmacophoresdb',
                        default='pharmacophores.db',
                        help='Name of pharmacophore db file (default: %(default)s)')
    parser.set_defaults(func=dir2db_run)


def make_pharmacophores_parser(subparsers):
    """Creates a parser for pharmacophores sub commands
    
    Args:
        subparsers (argparse.ArgumentParser): Parser to which to add sub commands to 

    """
    sc= subparsers.add_parser('pharmacophores', help='Pharmacophores').add_subparsers()
    add_sc(sc)

