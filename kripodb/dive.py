import csv
import json
import logging
import math
from os.path import basename

from progressbar import ProgressBar
from rdkit.Chem.Descriptors import HeavyAtomMolWt
import six

from .db import FragmentsDb
from .frozen import FrozenSimilarityMatrix
from .pdb import PdbReport


def dive_sphere(inputfile, outputfile, onlyfrag1):
    """Export fragments as DiVE formatted sphere

    Args:
        inputfile (str): fragments db input file
        outputfile (file): fragments dive output file
        onlyfrag1 (bool): Only \*_frag1

    """
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


def dive_export(fragmentsdb, uniprot_annot, pdbtags, propnames, props):
    """Writes metdata props for DiVE visualization

    Args:
        fragmentsdb (str): Filename fo fragments db file
        uniprot_annot (file): Readable file object with uniprot gene and family mapping as tsv
        pdbtags (list): List of readable file objects to tag pdb by filename
        propnames (file): Writable file object to write prop names to
        props (file): Writeable file object to write props to
    """
    db = FragmentsDb(fragmentsdb)

    data = {}
    dive_get_fragments(db, data)
    dive_merge_uniprot(uniprot_annot, data)
    dive_merge_pdb(data)
    dive_merge_pdb_tag(pdbtags, data)

    dump_propnames(propnames, pdbtags is not None)
    dump_props(data, props)


def dive_get_fragments(db, data):
    # TODO add organism column to pdb data table
    sql = '''SELECT
            frag_id,
            pdb_code as pdb,
            het_code as het,
            frag_nr as fragment,
            pdb_title as title,
            uniprot_acc as uniprot,
            uniprot_name as protein,
            smiles,
            mol
          FROM
            fragments
            JOIN pdbs USING (pdb_code)
            LEFT JOIN molecules USING (frag_id)
            '''
    for row in db.cursor.execute(sql):
        cols = row.keys()
        frag_id = row[0]
        data[frag_id] = {}
        mol = row[-1]
        if mol:
            data[frag_id]['weight'] = HeavyAtomMolWt(mol)
            # TODO add other Lipinski parameters aswell http://www.rdkit.org/Python_Docs/rdkit.Chem.Lipinski-module.html
        for col in cols[1:-1]:
            data[frag_id][col] = row[col]


def dive_merge_uniprot(uniprot_annot_fn, data):
    pdb2uniprot_accs = {}
    uniprot_acc2gene = {}
    uniprot_acc2family = {}
    logging.warning('Loading uniprot')
    reader = csv.reader(uniprot_annot_fn, delimiter='\t')
    next(reader)
    for row in reader:
        if row[1]:
            uniprot_acc2gene[row[0]] = row[1]
        if row[2]:
            uniprot_acc2family[row[0]] = row[2].split(', ')
        if row[3]:
            for pdb in row[3].split(';'):
                # Kripo uses lowercase pdb code, while rest of world uses uppercase
                pdb2uniprot_accs[pdb.lower()] = row[0]

    for frag_id in data:
        record = data[frag_id]
        pdb_code = record['pdb']
        if pdb_code in pdb2uniprot_accs:
            uniprot_acc = pdb2uniprot_accs[pdb_code]
            if uniprot_acc != record['uniprot']:
                record['uniprot'] = uniprot_acc
            if uniprot_acc in uniprot_acc2gene:
                record['gene'] = uniprot_acc2gene[uniprot_acc]
            if uniprot_acc in uniprot_acc2family:
                record['families'] = uniprot_acc2family[uniprot_acc]


def dive_merge_pdb(data):
    logging.warning('Loading pdb from internet')
    pdb_report = PdbReport(fields=['source'])
    pdb2organism = {pdb['structureId'].lower(): pdb['source'] for pdb in pdb_report.fetch() if pdb['source']}
    for frag_id in data:
        record = data[frag_id]
        pdb_code = record['pdb']
        if pdb_code in pdb2organism:
            organism = pdb2organism[pdb_code]
            record['organism'] = organism


def dive_merge_pdb_tag(pdbtags, data):
    logging.warning('Loading pdb tags')
    tags = {}
    for pdbtagfile in pdbtags:
        tagname = basename(pdbtagfile.name)
        for line in pdbtagfile:
            tags[line.strip().lower()] = tagname
    for frag_id in data:
        record = data[frag_id]
        pdb_code = record['pdb']
        if pdb_code in tags:
            record['pdbtag'] = tags[pdb_code]


def dump_propnames(propnamesfn, has_pdbtag):
    propnames = [
        'pdb',
        'het',
        'fragment',
        'title',
        'smiles',
        'weight',
        'uniprot',
        'protein',
        'organism',
        'gene',
    ]
    if has_pdbtag:
        propnames.append('pdbtag')
    propnames.extend([
        'family0',
        'family1',
        'family2',
        'family3',
        'family4',
    ])
    json.dump(propnames, propnamesfn)


def dump_props(props, propsfn):
    for frag_id, v in six.iteritems(props):
        propsfn.write(frag_id)
        propsfn.write(' ')
        fields = [
            'pdb:' + v['pdb'],
            'het:' + v['het'],
            'fragment:' + str(v['fragment']),
            '"title:' + v['title'] + '"',
        ]
        if 'smiles' in v and v['smiles']:
            fields.append('smiles:' + v['smiles'])
        else:
            fields.append('')
        if 'weight' in v:
            fields.append('{0:.2f}'.format(v['weight']))
        else:
            fields.append('')
        if v['uniprot']:
            fields.append('uniprot:' + v['uniprot'].split('#')[0])
        else:
            fields.append('')
        if 'protein' in v and v['protein']:
            fields.append('"protein:' + v['protein'] + '"')
        else:
            fields.append('')
        if 'organism' in v:
            fields.append('"organism:' + v['organism'] + '"')
        else:
            fields.append('')
        if 'gene' in v:
            fields.append('"gene:' + v['gene'] + '"')
        else:
            fields.append('')
        if 'pdbtag' in v:
            fields.append('pdbtag:' + v['pdbtag'])
        else:
            fields.append('')
        if 'families' in v:
            for idx, fam in enumerate(v['families']):
                fields.append('"family' + str(idx) + ':' + fam + '"')

        propsfn.write(' '.join(fields))
        propsfn.write("\n")


def dense_dump(inputfile, outputfile, frag1only):
    """Dump dense matrix with zeros included

    Args:
        inputfile (str): Filename of dense similarity matrix
        outputfile (file): Writeable file object
        frag1only (bool): Only dump frag1 fragments

    Returns:

    """
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
            if not score:
                continue
            yield (row_label, col_label, score)
