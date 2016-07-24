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
from __future__ import absolute_import

import argparse
import csv
import json
import math
from six.moves.urllib.request import urlopen

from progressbar import ProgressBar

from ..frozen import FrozenSimilarityMatrix
from ..db import FragmentsDb


def make_cclustera_parsers(subparsers):
    sc = subparsers.add_parser('cclustera', help='CClustera visualization utils').add_subparsers()
    fragments_sphere_sc(sc)
    cclustera_enrich_sc(sc)
    dense_dump_sc(sc)


def fragments_sphere_sc(subparsers):
    sc = subparsers.add_parser('fragments', help='Export fragments as cclustera sphere')
    sc.add_argument('inputfile', type=str,
                    help='Name of fragments db input file')
    sc.add_argument('outputfile', type=argparse.FileType('w'),
                    help='Name of fragments cclustera output file, use - for stdout')
    sc.add_argument('--onlyfrag1', action='store_true',
                    help='Only *_frag1 (default: %(default)s)')
    sc.set_defaults(func=cclustera_sphere)


def cclustera_sphere(inputfile, outputfile, onlyfrag1):
    frags_db = FragmentsDb(inputfile)
    nodes = {}

    # distribute fragments evenly on sphere using Fibonacci sphere algorithm
    # from http://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    samples = len(frags_db)

    sql = 'SELECT frag_id, pdb_code, het_code FROM fragments'
    if onlyfrag1:
        sql += ' WHERE frag_id LIKE "%_frag1"'
        frags_db.cursor.execute('SELECT count(*) FROM fragments WHERE frag_id LIKE "%_frag1"')
        samples = frags_db.cursor.fetchone()[0]

    rnd = 1.
    offset = 2. / samples
    increment = math.pi * (3. - math.sqrt(5.));

    frag_ids = frags_db.cursor.execute(sql)
    for i, frag in enumerate(frag_ids):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y, 2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        node_info = {
            'Path': [],
            'Coordinates': [x, y, z],
            'Categories': [frag[1], frag[2]],
            'Properties': []
        }
        nodes[frag[0]] = node_info

    json.dump(nodes, outputfile)


def cclustera_enrich_sc(sc):
    sc = sc.add_parser('enrich', help='Enrich cclustera data file')
    sc.add_argument('inputfile', type=argparse.FileType('r'),
                    help='Name of input file, user - for stdin')
    sc.add_argument('outputfile', type=argparse.FileType('w'),
                    help='Name of output file, use - for stdout')
    sc.set_defaults(func=cclustera_enrich)


def cclustera_enrich(inputfile, outputfile):
    # From uniprot download accession 2 gene symbol, family mapping
    # uniprot_url = 'http://www.uniprot.org/uniprot/?query=database:pdb&format=tab&columns=id,genes(PREFERRED),families,database(PDB)'
    # mapping = urlopen(uniprot_url)
    mapping = open('unprot_gene_family_pdb.csv', 'r')
    data = json.load(inputfile)

    enrich_fragments(data, mapping)

    json.dump(data, outputfile)


def enrich_fragments(data, mapping):
    """Adds Uniprot mappings to categories field of each fragment.

    Args:
        data (dict): Fragments in CClustera format
        mapping (File): Tab separated file with Uniprot mappings

    Returns:
        data
    """
    pdb2uniprot_accs = {}
    uniprot_acc2gene = {}
    uniprot_acc2family = {}
    reader = csv.reader(mapping, delimiter='\t')
    next(reader)
    for row in reader:
        if row[1]:
            uniprot_acc2gene[row[0]] = 'gene_' + row[1]
        if row[2]:
            uniprot_acc2family[row[0]] = row[2].split(', ')
        if row[3]:
            for pdb in row[3].split(';'):
                # Kripo uses lowercase pdb code, while rest of world uses uppercase
                pdb2uniprot_accs[pdb.lower()] = row[0]
    for frag_id in data:
        (pdb_code, het_code, frag_nr) = frag_id.split('_')
        cats = [pdb_code, het_code, frag_nr]
        if pdb_code in pdb2uniprot_accs:
            uniprot_acc = pdb2uniprot_accs[pdb_code]
            cats.append(uniprot_acc)
            if uniprot_acc in uniprot_acc2gene:
                cats.append(uniprot_acc2gene[uniprot_acc])
            if uniprot_acc in uniprot_acc2family:
                for fam in uniprot_acc2family[uniprot_acc]:
                    cats.append(fam)

        data[frag_id]['Categories'] = list(cats)


def dense_dump_sc(sc):
    """Dump dense matrix with zeros"""
    sc = sc.add_parser('dump', help='Dump dense matrix with zeros')
    sc.add_argument('inputfile', type=str,
                    help='Name of dense similarity matrix')
    sc.add_argument('outputfile', type=argparse.FileType('w'),
                    help='Name of output file, use - for stdout')
    sc.add_argument('--frag1only', action='store_true', help='Only *frag1 (default: %(default)s)')
    sc.set_defaults(func=dense_dump)


def dense_dump(inputfile, outputfile, frag1only):
    matrix = FrozenSimilarityMatrix(inputfile)
    writer = csv.writer(outputfile, delimiter='\t', lineterminator='\n')
    writer.writerow(['frag_id1', 'frag_id2', 'score'])
    writer.writerows(dense_dump_iter(matrix, frag1only))
    matrix.close()


def dense_dump_iter(matrix, frag1only):
    """Iterate dense matrix with zeros

    Args:
        matrix (FrozenSimilarityMatrix): Dense similarity matrix
        frag1only (bool): True to iterate over \*frag1 only

    Yields:
        (str, str, float): Fragment label pair and score
    """
    completed_frags = set()
    bar = ProgressBar()
    labels = [v.decode() for v in matrix.labels]
    for row_label in bar(labels):
        if frag1only and not row_label.endswith('frag1'):
            continue
        completed_frags.add(row_label)
        cols = matrix[row_label]
        for (col_label, score) in cols:
            if frag1only and not col_label.endswith('frag1'):
                continue
            if col_label in completed_frags:
                continue
            yield (row_label, col_label, score)
