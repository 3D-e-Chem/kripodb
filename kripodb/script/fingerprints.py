import argparse
import gzip
import sys
import tarfile

from .. import pairs, makebits
from ..db import FragmentsDb, FingerprintsDb
from ..modifiedtanimoto import calc_mean_onbit_density


def make_fingerprints_parser(subparsers):
    """Creates a parser for fingerprints sub commands

    Args:
        subparsers (argparse.ArgumentParser): Parser to which to add sub commands to
    """
    fp_sc = subparsers.add_parser('fingerprints', help='Fingerprints').add_subparsers()
    makebits2fingerprintsdb_sc(fp_sc)
    fingerprintsdb2makebits_sc(fp_sc)
    meanbitdensity_sc(fp_sc)
    similarity2query_sc(fp_sc)
    pairs_sc(fp_sc)


def pairs_sc(subparsers):
    sc_help = '''Calculate modified tanimoto similarity between fingerprints'''
    sc_description = '''

    Output formats:
    * tsv, tab separated id1,id2, similarity
    * hdf5, hdf5 file constructed with pytables with a, b and score, but but a and b have been replaced
      by numbers and similarity has been converted to scaled int

    When input has been split into chunks,
    use `--ignore_upper_triangle` flag for computing similarities between same chunk.
    This prevents storing pair a->b also as b->a.
    '''
    out_formats = ['tsv', 'hdf5']
    sc = subparsers.add_parser('similarities',
                               help=sc_help,
                               description=sc_description)
    sc.add_argument('fingerprintsfn1',
                    help='Name of reference fingerprints db file')
    sc.add_argument('fingerprintsfn2',
                    help='Name of query fingerprints db file')
    sc.add_argument('out_file',
                    help='Name of output file (use - for stdout)')
    sc.add_argument('--out_format',
                    choices=out_formats,
                    default='hdf5',
                    help='Format of output (default: %(default)s)')
    sc.add_argument('--fragmentsdbfn',
                    help='Name of fragments db file (only required for hdf5 format)')
    sc.add_argument('--mean_onbit_density',
                    help='Mean on bit density (default: %(default)s)',
                    type=float,
                    default=0.01)
    sc.add_argument('--cutoff',
                    type=float,
                    default=0.45,
                    help='Set Tanimoto cutoff (default: %(default)s)')
    sc.add_argument('--nomemory',
                    action='store_true',
                    help='Do not store query fingerprints in memory (default: %(default)s)')
    sc.add_argument('--ignore_upper_triangle',
                    action='store_true',
                    help='Ignore upper triangle (default: %(default)s)')
    sc.set_defaults(func=pairs_run)


def pairs_run(fingerprintsfn1, fingerprintsfn2,
              out_format, out_file,
              mean_onbit_density,
              cutoff,
              fragmentsdbfn,
              nomemory,
              ignore_upper_triangle):

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
                     nomemory,
                     ignore_upper_triangle)


def makebits2fingerprintsdb_sc(subparsers):
    sc = subparsers.add_parser('import', help='Add Makebits file to fingerprints db')
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
    sc = subparsers.add_parser('export',
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


def similarity2query_sc(subparsers):
    sc_help = 'Find the fragments closests to query based on fingerprints'
    sc = subparsers.add_parser('similar', help=sc_help)
    sc.add_argument('fingerprintsdb',
                    default='fingerprints.db',
                    help='Name of fingerprints db file')
    sc.add_argument('query', type=str, help='Query identifier or beginning of it')
    sc.add_argument('out', type=argparse.FileType('w'), help='Output file tabdelimited (query, hit, score)')
    sc.add_argument('--mean_onbit_density',
                    help='Mean on bit density (default: %(default)s)',
                    type=float,
                    default=0.01)
    sc.add_argument('--cutoff',
                    type=float,
                    default=0.55,
                    help='Set Tanimoto cutoff (default: %(default)s)')
    sc.add_argument('--memory',
                    action='store_true',
                    help='Store bitsets in memory (default: %(default)s)')
    sc.set_defaults(func=pairs.similarity2query)


def similarity2query_run(fingerprintsdb, query, out, mean_onbit_density, cutoff, memory):
    bitsets = FingerprintsDb(fingerprintsdb).as_dict()
    pairs.similarity2query(bitsets, query, out, mean_onbit_density, cutoff, memory)


def meanbitdensity_sc(subparsers):
    sc = subparsers.add_parser('meanbitdensity', help='Compute mean bit density of fingerprints')
    sc.add_argument('fingerprintsdb',
                    default='fingerprints.db',
                    help='Name of fingerprints db file (default: %(default)s)')
    sc.add_argument('--out', type=argparse.FileType('w'),
                    default='-',
                    help='Output file, default is stdout (default: %(default)s)')
    sc.set_defaults(func=meanbitdensity_run)


def meanbitdensity_run(fingerprintsdb, out):
    bitsets = FingerprintsDb(fingerprintsdb).as_dict()
    density = calc_mean_onbit_density(bitsets.values(), bitsets.number_of_bits)
    out.write("{0:.5}\n".format(density))