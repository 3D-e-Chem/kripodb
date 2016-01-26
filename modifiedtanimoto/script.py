# Copyright 2013 Netherlands eScience Center
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
import sys
import argparse
import logging
import tarfile
from dbm import IntbitsetDictDbm, combine_intbitsetdbms
from . import makebits
from . import pairs

DEFAULT_NUMBER_OF_BITS = 574331


def make_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    pairs_sc(subparsers)

    makebits2intbitsetdbm_sc(subparsers)

    id2label_sc(subparsers)

    combine_intbitsetdbms_sc(subparsers)

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
    sc.add_argument("bs_file1",
                    help="Name of reference fingerprint file")
    sc.add_argument("bs_file2",
                    help="Name of query fingerprint file")
    sc.add_argument("out_file",
                    help="Name of output file (use - for stdout)")
    sc.add_argument("--out_format",
                    choices=out_formats,
                    default='tsv',
                    help="Format of output")
    sc.add_argument("--id2label",
                    dest="id2label_file",
                    help="Id to label lookup tsv file")
    sc.add_argument("--number_of_bits",
                    type=int,
                    default=DEFAULT_NUMBER_OF_BITS)
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
    sc.set_defaults(func=pairs.dump_pairs)


def makebits2intbitsetdbm_sc(subparsers):
    sc_help = 'Convert Makebits file to intbitset dbm'
    sc_desc = '''Creates a dbm with anydbm package where
    key is the fingerprint identifier and the value is a
    fastdump-ed intbitset serialization.

    Creates files in working directory.

    The value can be converted back into a intbitset object by
    ```
    from intbitset import intbitset
    db = open('dbfilename')
    bs = intbitset()
    bs.fastload(db['some_id'])
    ```
    '''
    fc = argparse.RawDescriptionHelpFormatter
    sc = subparsers.add_parser('makebits2intbitsetdbm',
                               formatter_class=fc,
                               help=sc_help,
                               description=sc_desc)
    sc.add_argument("infile",
                    help="Name of makebits formatted fingerprint tar.gz file")
    sc.set_defaults(func=makebits2intbitsetdbm_run)


def makebits2intbitsetdbm(infile, outfile):
    gen = makebits.iter_file(infile)
    header = next(gen)
    fp_size = makebits.read_fp_size(header)
    logging.warn('Fingerprint size is {}'.format(fp_size))
    with IntbitsetDictDbm(outfile, fp_size) as bitsets:
        for fid, bitset in gen:
            bitsets[fid] = bitset


def makebits2intbitsetdbm_run(infile):
    if infile.endswith('tar.gz'):
        with tarfile.open(infile) as tar:
            for tarinfo in tar:
                if tarinfo.isfile():
                    outfile = tarinfo.name + '.db'
                    msg = 'Reading {} and writing {}'.format(infile, outfile)
                    logging.warn()
                    f = tar.extractfile(tarinfo)
                    makebits2intbitsetdbm(f, outfile)
                    f.close()
    elif infile.endswith('.fp'):
        f = open(infile)
        outfile = infile + '.db'
        logging.warn('Reading {} and writing {}'.format(infile, outfile))
        makebits2intbitsetdbm(f, outfile)
    else:
        msg = 'Unable to open {}, format not supported'.format(infile)
        raise Exception(msg)


def bitsets2id2label(infile):
    with IntbitsetDictDbm(infile, 0) as bitsets:
        for bsid, bslabel in enumerate(bitsets):
            print("{}\t{}".format(bsid, bslabel))


def id2label_sc(subparsers):
    sc_help = 'Create id2file from fingerprint/bitset file'
    sc_desc = '''Write bitset id and label to stdout'''
    sc = subparsers.add_parser('id2label',
                               help=sc_help,
                               description=sc_desc)
    sc.add_argument("infile",
                    help="Name of makebits formatted fingerprint tar.gz file")
    sc.set_defaults(func=bitsets2id2label)


def combine_intbitsetdbms_sc(subparsers):
    sc_help = 'Combine multiple intbitset dbm files into a single file'
    sc = subparsers.add_parser('combineintbitsetdbms', help=sc_help)
    sc.add_argument("infiles", type=str, nargs='+')
    sc.add_argument("outfile", type=str)
    sc.set_defaults(func=combine_intbitsetdbms)


def main(argv=sys.argv[1:]):
    parser = make_parser()
    args = parser.parse_args(argv)
    fargs = vars(args)
    func = args.func
    del(fargs['func'])
    func(**fargs)
