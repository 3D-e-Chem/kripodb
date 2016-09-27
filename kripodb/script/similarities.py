import argparse
import csv

import six
from tables import parameters
from .. import pairs
from ..db import FragmentsDb
from ..frozen import FrozenSimilarityMatrix
from ..hdf5 import SimilarityMatrix
from ..webservice.server import serve_app


def make_similarities_parser(subparsers):
    """Creates a parser for similarities sub commands

    Args:
        subparsers (argparse.ArgumentParser): Parser to which to add sub commands to
    """
    dm_sc = subparsers.add_parser('similarities', help='Similarity matrix').add_subparsers()
    similar_sc(dm_sc)
    merge_pairs_sc(dm_sc)
    simmatrix_export_sc(dm_sc)
    simmatrix_import_sc(dm_sc)
    simmatrix_filter_sc(dm_sc)
    similarity_freeze_sc(dm_sc)
    fpneigh2tsv_sc(dm_sc)
    serve_sc(dm_sc)


def similar_sc(subparsers):
    sc_help = 'Find the fragments closets to query based on similarity matrix'
    sc = subparsers.add_parser('similar', help=sc_help)
    sc.add_argument('pairsdbfn', type=str, help='hdf5 similarity matrix file or base url of kripodb webservice')
    sc.add_argument('query', type=str, help='Query fragment identifier')
    sc.add_argument('--out', type=argparse.FileType('w'), default='-',
                    help='Output file tab delimited (query, hit, similarity score)')
    sc.add_argument('--cutoff',
                    type=float,
                    default=0.55,
                    help='Similarity cutoff (default: %(default)s)')
    sc.set_defaults(func=pairs.similar_run)


def merge_pairs_sc(subparsers):
    sc = subparsers.add_parser('merge', help='Combine pairs files into a new file')
    sc.add_argument('ins', help='Input pair file in hdf5_compact format', nargs='+')
    sc.add_argument('out', help='Output pair file in hdf5_compact format')
    sc.set_defaults(func=pairs.merge)


def simmatrix_export_sc(subparsers):
    sc = subparsers.add_parser('export', help='Export similarity matrix to tab delimited file')
    sc.add_argument('simmatrixfn', type=str, help='Compact hdf5 similarity matrix filename')
    sc.add_argument('outputfile', type=argparse.FileType('w'),
                    help='Tab delimited output file, use - for stdout')
    sc.set_defaults(func=simmatrix_export_run)


def simmatrix_export_run(simmatrixfn, outputfile):
    """Export similarity matrix to tab delimited file

    Args:
        simmatrixfn (str): Compact hdf5 similarity matrix filename
        outputfile (file): Tab delimited output file

    """
    simmatrix = SimilarityMatrix(simmatrixfn)
    writer = csv.writer(outputfile, delimiter="\t", lineterminator='\n')
    writer.writerow(['frag_id1', 'frag_id2', 'score'])
    writer.writerows(simmatrix)
    simmatrix.close()


def simmatrix_import_sc(subparsers):
    sc = subparsers.add_parser('import',
                               help='Import similarity matrix from tab delimited file',
                               description='''When input has been split into chunks,
                                           use `--ignore_upper_triangle` flag for similarities between same chunk.
                                           This prevents storing pair a->b also as b->a.''')
    sc.add_argument('inputfile', type=argparse.FileType('r'),
                    help='Input file, use - for stdin')
    sc.add_argument('fragmentsdb',
                    default='fragments.db',
                    help='Name of fragments db file (default: %(default)s)')
    sc.add_argument('simmatrixfn', type=str, help='Compact hdf5 similarity matrix file, will overwrite file if it exists')
    sc.add_argument('--format',
                    choices=['tsv', 'fpneigh'],
                    default='fpneigh',
                    help='tab delimited (tsv) or fpneigh formatted input (default: %(default)s)')
    # Have to ask, because inputfile can be stdin so can't do 2 passes through file
    sc.add_argument('--nrrows',
                    type=int,
                    default=2**16,
                    help='Number of rows in inputfile (default: %(default)s)')
    sc.add_argument('--ignore_upper_triangle',
                    action='store_true',
                    help='Ignore upper triangle (default: %(default)s)')
    sc.set_defaults(func=simmatrix_import_run)


def simmatrix_import_run(inputfile, fragmentsdb, simmatrixfn, format, nrrows, ignore_upper_triangle=False):
    if format == 'tsv':
        simmatrix_import_tsv(inputfile, fragmentsdb, simmatrixfn, nrrows, ignore_upper_triangle)
    elif format == 'fpneigh':
        simmatrix_importfpneigh_run(inputfile, fragmentsdb, simmatrixfn, nrrows, ignore_upper_triangle)


def simmatrix_import_tsv(inputfile, fragmentsdb, simmatrixfn, nrrows, ignore_upper_triangle=False):
    frags = FragmentsDb(fragmentsdb)
    label2id = frags.label2id().materialize()
    simmatrix = SimilarityMatrix(simmatrixfn, 'w',
                                expectedlabelrows=len(label2id),
                                expectedpairrows=nrrows)

    reader = csv.reader(inputfile, delimiter="\t")
    # ignore header
    next(reader)

    # simmatrix wants score as float instead of str
    def csv_iter(rows):
        for row in rows:
            if row[0] == row[1]:
                continue
            if ignore_upper_triangle and row[0] > row[1]:
                continue
            row[2] = float(row[2])
            yield row

    simmatrix.update(csv_iter(reader), label2id)
    simmatrix.close()


def simmatrix_importfpneigh_run(inputfile, fragmentsdb, simmatrixfn, nrrows, ignore_upper_triangle=False):
    frags = FragmentsDb(fragmentsdb)
    label2id = frags.label2id().materialize()
    simmatrix = SimilarityMatrix(simmatrixfn, 'w',
                                expectedlabelrows=len(label2id),
                                expectedpairrows=nrrows)

    simmatrix.update(read_fpneighpairs_file(inputfile, ignore_upper_triangle), label2id)
    simmatrix.close()


def simmatrix_filter_sc(subparsers):
    sc = subparsers.add_parser('filter', help='Filter similarity matrix')
    sc.add_argument('input', type=str,
                    help='Input hdf5 similarity matrix file')
    sc.add_argument('output', type=str,
                    help='Output hdf5 similarity matrix file, will overwrite file if it exists')
    sc.add_argument('--fragmentsdb',
                    default='fragments.db',
                    help='Name of fragments db file (default: %(default)s)')
    sc.set_defaults(func=simmatrix_filter)


def simmatrix_filter(input, output, fragmentsdb):
    simmatrix_in = SimilarityMatrix(input)
    frags = FragmentsDb(fragmentsdb)
    print('Counting')
    expectedlabelrows = len(frags)
    labelsin = len(simmatrix_in.labels)
    expectedpairrows = int(len(simmatrix_in.pairs) * (float(expectedlabelrows) / labelsin))

    simmatrix_out = SimilarityMatrix(output,
                                    'w',
                                    expectedlabelrows=expectedlabelrows,
                                    expectedpairrows=expectedpairrows,
                                    )

    print('Building frag_id keep list')
    frag_labels2keep = set(frags.id2label().values())
    frag_ids2keep = set()
    for frag_label, frag_id in six.iteritems(simmatrix_in.labels.label2ids()):
        if frag_label in frag_labels2keep:
            frag_ids2keep.add(frag_id)

    print('Copying subset of pairs table')
    all_frags2keep = set(frag_ids2keep)
    hit = simmatrix_out.pairs.table.row
    for row in simmatrix_in.pairs.table:
        if row[0] in frag_ids2keep and row[1] in frag_ids2keep:
            hit['a'] = row[0]
            hit['b'] = row[1]
            hit['score'] = row[2]
            hit.append()
            hit['b'] = row[0]
            hit['a'] = row[1]
            hit['score'] = row[2]
            hit.append()
        elif row[0] in frag_ids2keep:
            hit['a'] = row[0]
            hit['b'] = row[1]
            hit['score'] = row[2]
            hit.append()
            all_frags2keep.add(row[1])
        elif row[1] in frag_ids2keep:
            hit['a'] = row[1]
            hit['b'] = row[0]
            hit['score'] = row[2]
            hit.append()
            all_frags2keep.add(row[0])

    print('Adding indices')
    simmatrix_out.pairs.add_indexes()

    print('Copying subset of labels table')
    hit = simmatrix_out.labels.table.row
    for row in simmatrix_in.labels.table:
        if row[0] in all_frags2keep:
            hit['frag_id'] = row[0]
            hit['label'] = row[1]
            hit.append()

    simmatrix_in.close()
    simmatrix_out.close()


def similarity_freeze_sc(subparsers):
    sc = subparsers.add_parser('freeze', help='Optimize similarity matrix for reading')
    sc.add_argument('in_fn', type=str, help='Input pairs file')
    sc.add_argument('out_fn', type=str, help='Output array file, file is overwritten')
    sc.add_argument('-f', '--frame_size', type=int, default=10**8, help='Size of frame (default: %(default)s)')
    sc.add_argument('-m', '--memory', type=int, default=1, help='Memory cache in Gigabytes (default: %(default)s)')
    sc.add_argument('-l', '--limit', type=int, help='Number of pairs to copy, None for no limit (default: %(default)s)')
    sc.add_argument('-s', '--single_sided', action='store_true', help='Store half matrix (default: %(default)s)')
    sc.set_defaults(func=similarity_freeze_run)


def similarity_freeze_run(in_fn, out_fn, frame_size, memory, limit, single_sided):
    dm = SimilarityMatrix(in_fn, 'r')
    parameters.CHUNK_CACHE_SIZE = memory * 1024 ** 3
    parameters.CHUNK_CACHE_NELMTS = 2 ** 14
    dfm = FrozenSimilarityMatrix(out_fn, 'w')
    dfm.from_pairs(dm, frame_size, limit, single_sided)
    dm.close()
    dfm.close()


def read_fpneighpairs_file(inputfile, ignore_upper_triangle=False):
    """Read fpneigh formatted similarity matrix file.

    Args:
        inputfile (file): File object to read
        ignore_upper_triangle (bool): Ignore upper triangle of input

    Yields:
        Tuple((Str,Str,Float)): List of (query fragment identifier, hit fragment identifier, similarity score)

    """
    current_query = None
    reader = csv.reader(inputfile, delimiter=' ', skipinitialspace=True)
    for row in reader:
        if len(row) == 2 and current_query != row[0]:
            if ignore_upper_triangle and current_query > row[0]:
                continue
            yield (current_query, row[0], float(row[1]))
        elif len(row) == 4:
            current_query = row[3][:-1]


def fpneigh2tsv_sc(subparsers):
    sc = subparsers.add_parser('fpneigh2tsv', help='Convert fpneigh formatted file to tab delimited file')
    sc.add_argument('inputfile', type=argparse.FileType('r'),
                    help='Input file, use - for stdin')
    sc.add_argument('outputfile', type=argparse.FileType('w'),
                    help='Tab delimited output file, use - for stdout')
    sc.set_defaults(func=fpneigh2tsv_run)


def fpneigh2tsv_run(inputfile, outputfile):
    reader = read_fpneighpairs_file(inputfile)
    writer = csv.writer(outputfile, delimiter="\t", lineterminator='\n')
    writer.writerow(['frag_id1', 'frag_id2', 'score'])
    writer.writerows(reader)


def serve_sc(subparsers):
    sc = subparsers.add_parser('serve', help='Serve similarity matrix as webservice')
    sc.add_argument('matrix', type=str, help='Filename of similarity matrix hdf5 file')
    sc.add_argument('--internal_port',
                    type=int,
                    default=8084,
                    help='TCP port on which to listen (default: %(default)s)')
    sc.add_argument('--external_url',
                    type=str,
                    default='http://localhost:8084/kripo',
                    help='URL which should be used in Swagger spec (default: %(default)s)')

    sc.set_defaults(func=serve_app)