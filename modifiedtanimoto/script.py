# Copyright 2013 Netherlands eScience Center
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
import gzip
import sys
import tarfile
from modifiedtanimoto.db import FragmentsDb
from . import makebits
from . import pairs
from algorithm import calc_mean_onbit_density


def make_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    makebits2fragmentsdb_sc(subparsers)

    fragmentsdb2makebits_sc(subparsers)

    id2label_sc(subparsers)

    meanbitdensity_sc(subparsers)

    distance2query_sc(subparsers)

    pairs_sc(subparsers)

    return parser


def pairs_sc(subparsers):
    sc_help = '''Generate pairs from 2 bitset dicts
    with modified tanimoto distance'''
    sc_description = '''

    Output formats:
    * tsv, tab seperated id1,id2, distance
    * tsv_compact, same format as tsv, but id1 and id2 have been replaced
      by numbers and distance has been converted to scaled int
    * hdf5, hdf5 file contstructed with pytables with id1, id2 and distance
    * hdf5_compact, same format as hdf5, but same compacting tsv_distance

    '''
    out_formats = ['tsv', 'tsv_compact', 'hdf5', 'hdf5_compact']
    sc = subparsers.add_parser('pairs',
                               help=sc_help,
                               description=sc_description)
    sc.add_argument("fragmentsfn1",
                    help="Name of reference fragments db file")
    sc.add_argument("fragmentsfn2",
                    help="Name of query fragments db file")
    sc.add_argument("out_file",
                    help="Name of output file (use - for stdout)")
    sc.add_argument("--out_format",
                    choices=out_formats,
                    default='tsv',
                    help="Format of output")
    sc.add_argument("--id2label",
                    dest="id2label_file",
                    help="Id to label lookup tsv file")
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
                    default=100,
                    help=ph)
    sc.add_argument("--memory",
                    action='store_true',
                    help='Store query bitsets in memory')
    sc.set_defaults(func=pairs_run)


def pairs_run(fragmentsfn1, fragmentsfn2,
              out_format, out_file,
              mean_onbit_density,
              cutoff,
              id2label_file,
              precision, memory):

    bitsets1 = FragmentsDb(fragmentsfn1).bitsets()
    bitsets2 = FragmentsDb(fragmentsfn2).bitsets()

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
                     id2label_file,
                     precision,
                     memory
                     )


def makebits2fragmentsdb_sc(subparsers):
    sc = subparsers.add_parser('makebits2fragmentsdb', help='Add Makebits file to fragments db')
    sc.add_argument('infiles', nargs='+', type=argparse.FileType('r'), metavar='infile',
                    help='Name of makebits formatted fingerprint file (.tar.gz or not packed)')
    sc.add_argument('outfile', help='Name of fragments db file', default='fragments.db')
    sc.set_defaults(func=makebits2fragmentsdb)


def makebits2fragmentsdb_single(infile, bitsets):
    gen = makebits.iter_file(infile)
    header = next(gen)
    number_of_bits = makebits.read_fp_size(header)
    bitsets.number_of_bits = number_of_bits
    bitsets.update(gen)


def makebits2fragmentsdb(infiles, outfile):
    bitsets = FragmentsDb(outfile).bitsets()
    for infile in infiles:
        if infile.name.endswith('tar.gz'):
            with tarfile.open(fileobj=infile) as tar:
                for tarinfo in tar:
                    if tarinfo.isfile():
                        f = tar.extractfile(tarinfo)
                        makebits2fragmentsdb_single(f, bitsets)
                        f.close()
        else:
            makebits2fragmentsdb_single(infile, bitsets)


def fragmentsdb2makebits_sc(subparsers):
    sc = subparsers.add_parser('fragmentsdb2makebits',
                               help='Dump bitsets in fragments db to makebits file')

    sc.add_argument('infile',
                    default='fragments.db',
                    help='Name of fragments db file')
    sc.add_argument('outfile',
                    type=argparse.FileType('w'),
                    help='Name of makebits formatted fingerprint file (or - for stdout)')
    sc.set_defaults(func=fragmentsdb2makebits)


def fragmentsdb2makebits(infile, outfile):
    bitsets = FragmentsDb(infile).bitsets()
    makebits.write_file(bitsets.number_of_bits, bitsets, outfile)


def bitsets2id2label(infile):
    bitsets = FragmentsDb(infile).bitsets()
    for bsid, bslabel in enumerate(bitsets):
        print("{}\t{}".format(bsid, bslabel))


def id2label_sc(subparsers):
    sc_help = 'Create id2file from fragments db file'
    sc_desc = '''Write bitset id and label to stdout'''
    sc = subparsers.add_parser('id2label',
                               help=sc_help,
                               description=sc_desc)
    sc.add_argument('infile',
                    default='fragments.db',
                    help='Name of fragments db file')
    sc.set_defaults(func=bitsets2id2label)


def distance2query_sc(subparsers):
    sc_help = 'Find the fragments closests to query'
    sc = subparsers.add_parser('distance2query', help=sc_help)
    sc.add_argument("fragmentsdb",
                    default='fragments.db',
                    help="Name of fragments db file")
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


def meanbitdensity_sc(subparsers):
    sc = subparsers.add_parser('meanbitdensity', help='Compute mean bit density of bitsets')
    sc.add_argument("fragmentsdb",
                    default='fragments.db',
                    help="Name of fragments db file")
    sc.set_defaults(func=meanbitdensity_run)


def meanbitdensity_run(fragmentsdb):
    bitsets = FragmentsDb(fragmentsdb).bitsets()
    print(calc_mean_onbit_density(bitsets, bitsets.number_of_bits))


def main(argv=sys.argv[1:]):
    parser = make_parser()
    args = parser.parse_args(argv)
    fargs = vars(args)
    func = args.func
    del(fargs['func'])
    func(**fargs)
