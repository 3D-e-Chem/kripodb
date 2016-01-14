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
from . import dbm
from . import makebits
from . import matrix


def make_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    matrix_sc(subparsers)

    makebits2intbitsetdbm_sc(subparsers)

    return parser


def matrix_sc(subparsers):
    sc_help = 'Generate pairs from 2 bitset collections with modified tanimoto distance'
    sc_description = '''

    '''

    sc = subparsers.add_parser('pairs', help=sc_help, description=sc_description)
    sc.add_argument("--infile",
                    dest="bs_file1",
                    help="Name of reference fingerprint file")
    sc.add_argument("--queryfile",
                    dest="bs_file2",
                    help="Name of query fingerprint file")
    sc.add_argument("--out_format",
                    choices=['tsv', 'tsv_numbered', 'hdf5'],
                    help="Format of output")
    sc.add_argument("--out_file",
                    help="Name of output file (use - for stdout)")
    sc.add_argument("--id2label",
                    dest="id2label_file",
                    help="Id to label lookup tsv file")
    sc.add_argument("--size",
                    type=int,
                    default=574331)
    sc.add_argument("--mean_onbit_density",
                    type=float,
                    default=0.01)
    sc.add_argument("--cutoff",
                    type=float,
                    default=0.45,
                    help="Set Tanimoto cutoff")
    sc.set_defaults(func=matrix.dump_matrix)


def makebits2intbitsetdbm_sc(subparsers):
    sc_help = 'Convert Makebits file to intbitset dbm'
    sc_desc = '''Creates a dbm with anydbm package where
    key is the fingerprint identifier and the value is a
    fastdump-ed intbitset serialization.

    The value can be converted back into a intbitset object by
    ```
    from intbitset import intbitset
    db = open('dbfilename')
    bs = intbitset()
    bs.fastload(db['some_id'])
    ```
    '''
    sc = subparsers.add_parser('makebits2intbitsetdbm',
                               formatter_class=argparse.RawDescriptionHelpFormatter,
                               help=sc_help,
                               description=sc_desc)
    sc.add_argument("infile",
                    help="Name of makebits formatted fingerprint tar.gz file")
    sc.set_defaults(func=makebits2intbitsetdbm)


def makebits2intbitsetdbm(infile):
    if infile.endswith('tar.gz'):
        with tarfile.open(infile) as tar:
            for tarinfo in tar:
                if tarinfo.isfile():
                    f = tar.extractfile(tarinfo)
                    logging.warn('Reading {}'.format(tarinfo.name))
                    (bitsets, fp_size) = makebits.read_file(f)
                    f.close()
                    logging.warn('Fingerprint size is {}'.format(fp_size))
                    outfile = tarinfo.name + '.db'
                    logging.warn('Writing {}'.format(outfile))
                    dbm.dump_bitsets(bitsets, outfile)
    else:
        msg = 'Unable to open {}, format not supported'.format(infile)
        raise Exception(msg)


def main(argv=sys.argv[1:]):
    parser = make_parser()
    args = parser.parse_args(argv)
    fargs = vars(args)
    func = args.func
    del(fargs['func'])
    func(**fargs)
