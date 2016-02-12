# Copyright 2016 Netherlands eScience Center
#
# Licensed under the Apache License, Version 2.0 (the 'License");
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
import argparse
import csv
import gzip
import logging
import shelve
import sys
import tarfile

from rdkit.Chem.rdmolfiles import SDMolSupplier

from . import makebits
from . import pairs
from .db import FragmentsDb, FingerprintsDb
from .hdf5 import DistanceMatrix
from .pdb import PdbReport
from .modifiedtanimoto import calc_mean_onbit_density
from .version import __version__


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version=__version__)
    subparsers = parser.add_subparsers()

    makebits2fingerprintsdb_sc(subparsers)

    fingerprintsdb2makebits_sc(subparsers)

    meanbitdensity_sc(subparsers)

    distance2query_sc(subparsers)

    similar_sc(subparsers)

    pairs_sc(subparsers)

    shelve2fragmentsdb_sc(subparsers)

    sdf2fragmentsdb_sc(subparsers)

    pdb2fragmentsdb_sc(subparsers)

    merge_pairs_sc(subparsers)

    distmatrix_export_sc(subparsers)

    distmatrix_import_sc(subparsers)

    return parser


def pairs_sc(subparsers):
    sc_help = '''Calculate modified tanimoto distance between fingerprints'''
    sc_description = '''

    Output formats:
    * tsv, tab seperated id1,id2, distance
    * hdf5, hdf5 file contstructed with pytables with a, b and score, but but a and b have been replaced
      by numbers and distance has been converted to scaled int
    '''
    out_formats = ['tsv', 'hdf5']
    sc = subparsers.add_parser('pairs',
                               help=sc_help,
                               description=sc_description)
    sc.add_argument("fingerprintsfn1",
                    help="Name of reference fingerprints db file")
    sc.add_argument("fingerprintsfn2",
                    help="Name of query fingerprints db file")
    sc.add_argument("out_file",
                    help="Name of output file (use - for stdout)")
    sc.add_argument("--out_format",
                    choices=out_formats,
                    default='hdf5',
                    help="Format of output")
    sc.add_argument("--fragmentsdbfn",
                    help='Name of fragments db file (only required for hdf5 format)')
    sc.add_argument("--mean_onbit_density",
                    type=float,
                    default=0.01)
    sc.add_argument("--cutoff",
                    type=float,
                    default=0.45,
                    help="Set Tanimoto cutoff")
    ph = '''Distance precision for compact formats,
    distance range from 0..<precision>'''
    sc.add_argument("--precision",
                    type=int,
                    default=65535,
                    help=ph)
    sc.add_argument("--nomemory",
                    action='store_true',
                    help='Do not store query fingerprints in memory')
    sc.set_defaults(func=pairs_run)


def pairs_run(fingerprintsfn1, fingerprintsfn2,
              out_format, out_file,
              mean_onbit_density,
              cutoff,
              fragmentsdbfn,
              precision, nomemory):

    if 'hdf5' in out_format and fragmentsdbfn is None:
        raise Exception('Hdf5 format requires fragments db')

    label2id = {}
    if fragmentsdbfn is not None:
        label2id = FragmentsDb(fragmentsdbfn).label2id().materialize()

    bitsets1 = FingerprintsDb(fingerprintsfn1).as_dict()
    bitsets2 = FingerprintsDb(fingerprintsfn2).as_dict()

    if bitsets1.number_of_bits != bitsets2.number_of_bits:
        raise Exception('Number of bits is not the same')

    out = sys.stdout
    if out_file != '-' and out_format.startswith('tsv'):
        if out_file.endswith('gz'):
            out = gzip.open(out_file, 'w')
        else:
            out = open(out_file, 'w')

    pairs.dump_pairs(bitsets1,
                     bitsets2,
                     out_format,
                     out_file,
                     out,
                     bitsets1.number_of_bits,
                     mean_onbit_density,
                     cutoff,
                     label2id,
                     precision,
                     nomemory)


def makebits2fingerprintsdb_sc(subparsers):
    sc = subparsers.add_parser('makebits2fingerprintsdb', help='Add Makebits file to fingerprints db')
    sc.add_argument('infiles', nargs='+', type=argparse.FileType('r'), metavar='infile',
                    help='Name of makebits formatted fingerprint file (.tar.gz or not packed or - for stdin)')
    sc.add_argument('outfile', help='Name of fingerprints db file', default='fingerprints.db')
    sc.set_defaults(func=makebits2fingerprintsdb)


def makebits2fingerprintsdb_single(infile, bitsets):
    gen = makebits.iter_file(infile)
    header = next(gen)
    number_of_bits = makebits.read_fp_size(header)
    bitsets.number_of_bits = number_of_bits
    bitsets.update(gen)


def makebits2fingerprintsdb(infiles, outfile):
    bitsets = FingerprintsDb(outfile).as_dict()
    for infile in infiles:
        if infile.name.endswith('tar.gz'):
            with tarfile.open(fileobj=infile) as tar:
                for tarinfo in tar:
                    if tarinfo.isfile():
                        f = tar.extractfile(tarinfo)
                        makebits2fingerprintsdb_single(f, bitsets)
                        f.close()
        else:
            makebits2fingerprintsdb_single(infile, bitsets)


def fingerprintsdb2makebits_sc(subparsers):
    sc = subparsers.add_parser('fingerprintsdb2makebits',
                               help='Dump bitsets in fingerprints db to makebits file')

    sc.add_argument('infile',
                    default='fingerprints.db',
                    help='Name of fingerprints db file')
    sc.add_argument('outfile',
                    type=argparse.FileType('w'),
                    help='Name of makebits formatted fingerprint file (or - for stdout)')
    sc.set_defaults(func=fingerprintsdb2makebits)


def fingerprintsdb2makebits(infile, outfile):
    bitsets = FingerprintsDb(infile).as_dict()
    makebits.write_file(bitsets.number_of_bits, bitsets, outfile)


def distance2query_sc(subparsers):
    sc_help = 'Find the fragments closests to query based on fingerprints'
    sc = subparsers.add_parser('distance2query', help=sc_help)
    sc.add_argument("fingerprintsdb",
                    default='fingerprints.db',
                    help="Name of fingerprints db file")
    sc.add_argument("query", type=str, help='Query identifier or beginning of it')
    sc.add_argument("out", type=argparse.FileType('w'), help='Output file tabdelimited (query, hit, score)')
    sc.add_argument("--mean_onbit_density",
                    type=float,
                    default=0.01)
    sc.add_argument("--cutoff",
                    type=float,
                    default=0.55,
                    help="Set Tanimoto cutoff")
    sc.add_argument("--memory",
                    action='store_true',
                    help='Store bitsets in memory')
    sc.set_defaults(func=pairs.distance2query)


def distance2query_run(fingerprintsdb, query, out, mean_onbit_density, cutoff, memory):
    bitsets = FingerprintsDb(fingerprintsdb).as_dict()
    pairs.distance2query(bitsets, query, out, mean_onbit_density, cutoff, memory)


def similar_sc(subparsers):
    sc_help = 'Find the fragments closets to query based on distance matrix'
    sc = subparsers.add_parser('similar', help=sc_help)
    sc.add_argument("pairsdbfn", type=str, help='Compact hdf5 distance matrix file')
    sc.add_argument("query", type=str, help='Query fragment identifier')
    sc.add_argument("--out", type=argparse.FileType('w'), default='-',
                    help='Output file tab delimited (query, hit, distance score)')
    sc.add_argument("--cutoff",
                    type=float,
                    default=0.55,
                    help="Distance cutoff")
    sc.set_defaults(func=pairs.similar_run)


def meanbitdensity_sc(subparsers):
    sc = subparsers.add_parser('meanbitdensity', help='Compute mean bit density of fingerprints')
    sc.add_argument("fingerprintsdb",
                    default='fingerprints.db',
                    help="Name of fingerprints db file")
    sc.add_argument("--out", type=argparse.FileType('w'),
                    default='-',
                    help='Output file, default is stdout')
    sc.set_defaults(func=meanbitdensity_run)


def meanbitdensity_run(fingerprintsdb, out):
    bitsets = FingerprintsDb(fingerprintsdb).as_dict()
    density = calc_mean_onbit_density(bitsets.values(), bitsets.number_of_bits)
    out.write("{0:.5f}\n".format(density))


def shelve2fragmentsdb_sc(subparsers):
    sc = subparsers.add_parser('shelve2fragmentsdb', help='Add fragments from shelve to sqlite')
    sc.add_argument('shelvefn', type=str)
    sc.add_argument('fragmentsdb',
                    default='fragments.db',
                    help="Name of fragments db file")
    sc.set_defaults(func=shelve2fragmentsdb_run)


def shelve2fragmentsdb_run(shelvefn, fragmentsdb):
    myshelve = shelve.open(shelvefn, 'r')
    frags = FragmentsDb(fragmentsdb)
    frags.add_fragments_from_shelve(myshelve)


def sdf2fragmentsdb_sc(subparsers):
    sc = subparsers.add_parser('sdf2fragmentsdb', help='Add fragments sdf to sqlite')
    sc.add_argument('sdffns', help='SDF filename', nargs='+')
    sc.add_argument("fragmentsdb",
                    default='fragments.db',
                    help="Name of fragments db file")

    sc.set_defaults(func=sdf2fragmentsdb_run)


def sdf2fragmentsdb_run(sdffns, fragmentsdb):
    frags = FragmentsDb(fragmentsdb)
    for sdffn in sdffns:
        logging.warn('Parsing {}'.format(sdffn))
        suppl = SDMolSupplier(sdffn)
        frags.add_molecules(suppl)


def pdb2fragmentsdb_sc(subparsers):
    sc = subparsers.add_parser('pdb2fragmentsdb', help='Add pdb metadata from RCSB PDB website to fragment sqlite db')
    sc.add_argument("fragmentsdb",
                    default='fragments.db',
                    help="Name of fragments db file")

    sc.set_defaults(func=pdb2fragmentsdb_run)


def pdb2fragmentsdb_run(fragmentsdb):
    pdb_report = PdbReport()
    pdbs = pdb_report.fetch()
    frags = FragmentsDb(fragmentsdb)
    frags.add_pdbs(pdbs)


def merge_pairs_sc(subparsers):
    sc = subparsers.add_parser('mergepairs', help='Combine pairs files into a new file')
    sc.add_argument('ins', help='Input pair file in hdf5_compact format', nargs='+')
    sc.add_argument('out', help='Output pair file in hdf5_compact format')
    sc.set_defaults(func=pairs.merge)


def distmatrix_export_sc(subparsers):
    sc = subparsers.add_parser('distmatrix_export', help='Export distance matrix to tab delimited file')
    sc.add_argument("distmatrixfn", type=str, help='Compact hdf5 distance matrix filename')
    sc.add_argument("outputfile", type=argparse.FileType('w'),
                    help='Tab delimited output file, use - for stdout')
    sc.set_defaults(func=distmatrix_export_run)


def distmatrix_export_run(distmatrixfn, outputfile):
    """Export distance matrix to tab delimited file

    Args:
        distmatrixfn (str): Compact hdf5 distance matrix filename
        outputfile (file): Tab delimited output file

    """
    distmatrix = DistanceMatrix(distmatrixfn)
    writer = csv.writer(outputfile, delimiter="\t", lineterminator='\n')
    writer.writerow(['frag_id1', 'frag_id2', 'score'])
    writer.writerows(distmatrix)
    distmatrix.close()


def distmatrix_import_sc(subparsers):
    sc = subparsers.add_parser('distmatrix_import', help='Import distance matrix from tab delimited file')
    sc.add_argument("inputfile", type=argparse.FileType('r'),
                    help='Input file, use - for stdin')
    sc.add_argument("fragmentsdb",
                    default='fragments.db',
                    help="Name of fragments db file")
    sc.add_argument("distmatrixfn", type=str, help='Compact hdf5 distance matrix file, will overwrite file if it exists')
    ph = '''Distance precision for compact formats,
    distance range from 0..<precision>'''
    sc.add_argument("--precision",
                    type=int,
                    default=65535,
                    help=ph)
    # Have to ask, because inputfile can be stdin so can't do 2 passes through file
    sc.add_argument("--nrrows",
                    type=int,
                    default=2**16,
                    help='Number of rows in inputfile')
    sc.set_defaults(func=distmatrix_import_run)


def distmatrix_import_run(inputfile, fragmentsdb, distmatrixfn, precision, nrrows):
    frags = FragmentsDb(fragmentsdb)
    label2id = frags.label2id()
    distmatrix = DistanceMatrix(distmatrixfn, 'w',
                                precision=precision,
                                expectedlabelrows=len(label2id),
                                expectedpairrows=nrrows)

    reader = csv.reader(inputfile, delimiter="\t")
    # ignore header
    reader.next()

    # distmatrix wants score as float instead of str
    def csv_iter(rows):
        for row in rows:
            row[2] = float(row[2])
            yield row

    distmatrix.update(csv_iter(reader), label2id)
    distmatrix.close()


def main(argv=sys.argv[1:]):
    """Main script function.

    Calls run method of selected sub commandos.

    Args:
        argv (list[str]): List of command line arguments

    """
    parser = make_parser()
    args = parser.parse_args(argv)
    fargs = vars(args)
    func = args.func
    del(fargs['func'])
    func(**fargs)
