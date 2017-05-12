import argparse
import csv

from tables import parameters
from .. import pairs
from ..db import FragmentsDb
from ..frozen import FrozenSimilarityMatrix
from ..hdf5 import SimilarityMatrix


def make_similarities_parser(subparsers):
    """Creates a parser for similarities sub commands

    Args:
        subparsers (argparse.ArgumentParser): Parser to which to add sub commands to
    """
    sc = subparsers.add_parser('similarities', help='Similarity matrix').add_subparsers()
    similar_sc(sc)
    merge_pairs_sc(sc)
    simmatrix_export_sc(sc)
    simmatrix_import_sc(sc)
    simmatrix_filter_sc(sc)
    similarity_freeze_sc(sc)
    similarity_thaw_sc(sc)
    fpneigh2tsv_sc(sc)
    histogram_sc(sc)


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
    sc.add_argument('--no_header', action='store_true', help='Output no header (default: %(default)s)')
    sc.add_argument('--frag1', action='store_true', help='Only output *frag1 fragments (default: %(default)s)')
    pdbhelp = 'Only output fragments which are from pdb code in file, one pdb code per line (default: %(default)s)'
    sc.add_argument('--pdb', type=argparse.FileType('r'), help=pdbhelp)
    sc.set_defaults(func=simmatrix_export_run)


def load_pdb_filter_file(pdbs_file):
    pdbs = set()
    for line in pdbs_file:
        pdbs.add(line.strip().lower())
    return pdbs


def pdb_filter(rows, pdbs):
    for row in rows:
        if row[0][:4] in pdbs and row[1][:4] in pdbs:
            yield row


def frag1_filter(rows):
    for row in rows:
        if row[0].endswith('frag1') and row[1].endswith('frag1'):
            yield row


def simmatrix_export_run(simmatrixfn, outputfile, no_header, frag1, pdb):
    """Export similarity matrix to tab delimited file

    Args:
        simmatrixfn (str): (Compact) hdf5 similarity matrix filename
        outputfile (file): Tab delimited output file
        no_header (bool): Output no header
        frag1 (bool): Only output \*frag1
        pdb (str): Filename with pdb codes inside

    """
    simmatrix = pairs.open_similarity_matrix(simmatrixfn)
    if pdb:
        pdbs = load_pdb_filter_file(pdb)
    else:
        pdbs = None
    writer = csv.writer(outputfile, delimiter="\t", lineterminator='\n')

    with_header = not no_header
    if with_header:
        writer.writerow(['frag_id1', 'frag_id2', 'score'])

    if frag1 and pdb:
        writer.writerows(pdb_filter(frag1_filter(simmatrix), pdbs))
    elif frag1:
        writer.writerows(frag1_filter(simmatrix))
    elif pdb:
        writer.writerows(pdb_filter(simmatrix, pdbs))
    else:
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
    sc.add_argument('--inputformat',
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


def simmatrix_import_run(inputfile, fragmentsdb, simmatrixfn, inputformat, nrrows, ignore_upper_triangle=False):
    if inputformat == 'tsv':
        simmatrix_import_tsv(inputfile, fragmentsdb, simmatrixfn, nrrows, ignore_upper_triangle)
    elif inputformat == 'fpneigh':
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
    group = sc.add_mutually_exclusive_group()
    group.add_argument('--fragmentsdb',
                    help='Name of fragments db file, '
                         'fragments in it will be kept as well as their pair counter parts.')
    group.add_argument('--skip', type=argparse.FileType('r'), help='File with fragment identifiers on each line to skip')
    sc.set_defaults(func=simmatrix_filter)


def simmatrix_filter(input, output, fragmentsdb, skip):
    simmatrix_in = SimilarityMatrix(input)
    if fragmentsdb:
        frags = FragmentsDb(fragmentsdb)
        expectedlabelrows = len(frags)
        labelsin = len(simmatrix_in.labels)
        expectedpairrows = int(len(simmatrix_in.pairs) * (float(expectedlabelrows) / labelsin))

        simmatrix_out = SimilarityMatrix(output,
                                         'w',
                                         expectedlabelrows=expectedlabelrows,
                                         expectedpairrows=expectedpairrows,
                                         )

        frag_labels2keep = set(frags.id2label().values())
        simmatrix_in.keep(simmatrix_out, frag_labels2keep)
    if skip:
        labels2skip = set()
        for line in skip:
            labels2skip.add(line.strip())
        labelsin = len(simmatrix_in.labels)
        expectedlabelrows = labelsin - len(labels2skip)
        expectedpairrows = int(len(simmatrix_in.pairs) * (float(expectedlabelrows) / labelsin))

        simmatrix_out = SimilarityMatrix(output,
                                         'w',
                                         expectedlabelrows=expectedlabelrows,
                                         expectedpairrows=expectedpairrows,
                                         )

        simmatrix_in.skip(simmatrix_out, labels2skip)

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


def similarity_thaw_sc(subparsers):
    sc = subparsers.add_parser('thaw', help='Optimize similarity matrix for writing')
    sc.add_argument('in_fn', type=str, help='Input packed frozen matrix file')
    sc.add_argument('out_fn', type=str, help='Output pairs file, file is overwritten')
    sc.add_argument('--nonzero_fraction',
                    type=float,
                    default=0.012,
                    help='Fraction of pairs which have score above threshold (default: %(default)s)')
    sc.set_defaults(func=similarity_thaw_run)


def similarity_thaw_run(in_fn, out_fn, nonzero_fraction):
    fsm = FrozenSimilarityMatrix(in_fn, 'r')
    nr_scores = int(fsm.scores.shape[0] * fsm.scores.shape[1] * nonzero_fraction)
    nr_labels = fsm.labels.shape[0]
    sm = SimilarityMatrix(out_fn, 'w', expectedpairrows=nr_scores, expectedlabelrows=nr_labels)
    fsm.to_pairs(sm)
    sm.close()
    fsm.close()


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


def histogram_sc(subparsers):
    sc = subparsers.add_parser('histogram', help='Distribution of similarity scores')
    sc.add_argument('inputfile', type=str, help='Filename of similarity matrix hdf5 file')
    sc.add_argument('outputfile', type=argparse.FileType('w'),
                    help='Tab delimited output file, use - for stdout')
    sc.add_argument('-f', '--frame_size', type=int, default=10**8, help='Size of frame (default: %(default)s)')
    sc.add_argument('-r', '--raw_score',
                    action='store_true',
                    help='Return raw score (16 bit integer) instead of fraction score')
    sc.add_argument('-l', '--lower_triangle',
                    action='store_true',
                    help='Return scores from lower triangle else return scores from upper triangle')
    sc.set_defaults(func=histogram)


def histogram(inputfile, outputfile, frame_size, raw_score, lower_triangle):
    matrix = pairs.open_similarity_matrix(inputfile)
    counts = matrix.count(frame_size=frame_size, raw_score=raw_score, lower_triangle=lower_triangle)
    writer = csv.writer(outputfile, delimiter="\t", lineterminator='\n')
    writer.writerow(['score', 'count'])
    writer.writerows(counts)
    matrix.close()
