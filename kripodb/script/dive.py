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
from __future__ import absolute_import, print_function

import argparse

from ..dive import dive_export, dive_sphere, dense_dump


def make_dive_parsers(subparsers):
    sc = subparsers.add_parser('dive', help='DiVE visualization utils').add_subparsers()
    fragments_sphere_sc(sc)
    dense_dump_sc(sc)
    dive_export_sc(sc)


def fragments_sphere_sc(subparsers):
    sc = subparsers.add_parser('fragments', help='Export fragments as DiVE formatted sphere')
    sc.add_argument('inputfile', type=str,
                    help='Name of fragments db input file')
    sc.add_argument('outputfile', type=argparse.FileType('w'),
                    help='Name of fragments dive output file, use - for stdout')
    sc.add_argument('--onlyfrag1', action='store_true',
                    help='Only *_frag1 (default: %(default)s)')
    sc.set_defaults(func=dive_sphere)


def dive_export_sc(sc):
    sc = sc.add_parser('export', help='Writes props for DiVE visualization')
    sc.add_argument('fragmentsdb', type=str,
                    help='Name of fragments db input file')
    uniprot_annot_help = '''Uniprot download accession 2 gene symbol, family mapping.
    Fetch "http://www.uniprot.org/uniprot/?query=database:pdb&format=tab&columns=id,genes(PREFERRED),families,database(PDB)"
    '''
    sc.add_argument('uniprot_annot', type=argparse.FileType('r'), help=uniprot_annot_help)
    sc.add_argument('--propnames',
                    type=argparse.FileType('w'),
                    help='Name of prop names file',
                    default='kripo.propnames.txt')
    sc.add_argument('--props',
                    type=argparse.FileType('w'),
                    help='Name of props file',
                    default='kripo.props.txt')
    sc.add_argument('--pdbtags', type=argparse.FileType('r'), action='append', help='Tag pdb in file by filename')
    sc.set_defaults(func=dive_export)


def dense_dump_sc(sc):
    """Dump dense matrix with zeros"""
    sc = sc.add_parser('dump', help='Dump dense matrix with zeros')
    sc.add_argument('inputfile', type=str,
                    help='Name of dense similarity matrix')
    sc.add_argument('outputfile', type=argparse.FileType('w'),
                    help='Name of output file, use - for stdout')
    sc.add_argument('--frag1only', action='store_true', help='Only *frag1 (default: %(default)s)')
    sc.set_defaults(func=dense_dump)



