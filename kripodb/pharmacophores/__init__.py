from os import path, walk

import tables
import gzip
from rdkit.Chem import ForwardSDMolSupplier

FEATURE_TYPES = [{
    'key': 'LIPO',
    'label': 'Hydrophobe',
    'color': 'ff33cc',
    'element': 'He',
}, {
    'key': 'POSC',
    'label': 'Positive charge',
    'color': 'ff9933',
    'element': 'P'
}, {
    'key': 'NEGC',
    'label': 'Negative charge',
    'color': '376092',
    'element': 'Ne'
}, {
    'key': 'HDON',
    'label': 'H-bond donor',
    'color': '00ff00',
    'element': 'As'
}, {
    'key': 'HACC',
    'label': 'H-bond acceptor',
    'color': 'bfbfbf',
    'element': 'O'
}, {
    'key': 'AROM',
    'label': 'Aromatic',
    'color': '00ffff',
    'element': 'Rn'
}]
"""
Types of pharmacophore feature types. List of dictionaries with the following keys:
    * key, short identifier of type
    * label, human readable label
    * color, hex rrggbb color
    * element, Element used in kripo pharmacophore sdfile for this type
"""
FEATURE_TYPE_KEYS = [r['key'] for r in FEATURE_TYPES]
FEATURE_TYPE_ATOM2KEY = {r['element']: r['key'] for r in FEATURE_TYPES}

PYTABLE_FILTERS = tables.Filters(complevel=6, complib='blosc')


class PharmacophoreRow(tables.IsDescription):
    """Table description for similarity pair

    Attributes:
        frag_id (str): Fragment identifier
    """
    frag_id = tables.StringCol(16)
    type = tables.EnumCol(FEATURE_TYPE_KEYS, FEATURE_TYPE_KEYS[0], base='uint8')
    x = tables.Float32Col()
    y = tables.Float32Col()
    z = tables.Float32Col()


class PharmacophoresDb(object):
    """Database for pharmacophores of fragments aka sub-pockets.

    Args:
        filename (str): File name of hdf5 file to write or read pharmacophores to/from
        mode (str): Can be 'r' for reading or 'w' for writing or 'a' for appending
        expectedrows (int): Expected number of pharmacophores.
            Required when hdf5 file is created, helps optimize compression
        **kwargs: Passed to tables.open_file

    Pharmacophore points of a fragment can be retrieved using::

        points = db['frag_id1']

    `points` is a list of points, each point is a tuple with following columns feature type key, x, y and z coordinate.
    The feature type key is defined in FEATURE_TYPES.

    Attributes:
        h5file (tables.File): Object representing an open hdf5 file
        points (PharmacophorePointsTable): HDF5 table that contains pharmacophore points

    """

    def __init__(self, filename, mode='r', expectedrows=0, **kwargs):
        self.h5file = tables.open_file(filename, mode, filters=PYTABLE_FILTERS, **kwargs)
        self.points = PharmacophorePointsTable(self.h5file, expectedrows)

    def close(self):
        """Closes the hdf5file

        Instead of calling close() explicitly, use context manager::

            with PharmacophoresDb('data/pharmacophores.h5') as db:
                points = db['frag_id1']
        """
        self.h5file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def add_dir(self, startdir):
        """Find \*_pphore.sd.gz \*_pphores.txt file pairs recursively in start directory and add them.

        Args:
            startdir (str): Path to a start directory
        """
        self.points.add_dir(startdir)

    def __getitem__(self, item):
        return self.points[item]

    def write_phar(self, outfile, frag_id):
        """Write pharmacophore of frag_id as phar format to outfile

        Args:
            outfile (file): File object to write to
            frag_id (str): Fragment identifier
        """
        outfile.write(as_phar(frag_id, self[frag_id]))


def read_pphore_gzipped_sdfile(sdfile):
    """Read a gzipped sdfile which contains pharmacophore points as atoms

    Args:
        sdfile (string): Path to filename

    Returns:
        list: List of Pharmacophore points
    """
    with gzip.open(sdfile) as gzfile:
        return read_pphore_sdfile(gzfile)


def read_pphore_sdfile(sdfile):
    """Read a gzipped sdfile which contains pharmacophore points as atoms

    Args:
        sdfile (file): File object with sdfile contents

    Returns:
        list: List of pharmacophore points
    """
    mols = list(ForwardSDMolSupplier(sdfile))
    mol = mols[0]
    conf = mol.GetConformer(0)
    points = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        point = (
            FEATURE_TYPE_ATOM2KEY[atom.GetSymbol()],
            float(pos.x),
            float(pos.y),
            float(pos.z),
        )
        points.append(point)
    return points


def read_fragtxtfile(fragtxtfile):
    """Read a fragment text file

    Args:
        fragtxtfile: Filename of fragment text file

    Returns:
        dict: Dictionary where key is fragment identifier and value is a list of pharmacophore point indexes.
    """
    with open(fragtxtfile) as f:
        return read_fragtxtfile_as_file(f)


def read_fragtxtfile_as_file(fileobject):
    """Read a fragment text file object which contains the pharmacophore point indexes for each fragment identifier.

    File format is a fragment on each line,
    the line is space separated with fragment_identifier followed by the pharmacophore point indexes.

    Args:
        fileobject (file): File object to read

    Returns:
        dict: Dictionary where key is fragment identifier and value is a list of pharmacophore point indexes.

    """
    frags = {}
    for line in fileobject:
        point_ids = line.strip().split(' ')
        frag_id = point_ids.pop(0)
        point_ids = [int(r) - 1 for r in point_ids]
        frags[frag_id] = point_ids
    return frags


def as_phar(frag_id, points):
    """Return pharmacophore in \*.phar format.

    See `align-it <http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/align-it/1.0.4/align-it.html#format>`_ for format description.

    Args:
        frag_id (str): Fragment identifier
        points (list): List of points where each point is (key,x,y,z)

    Returns:
        str: Pharmacophore is \*.phar format

    """
    lines = [frag_id]
    for point in points:
        # Kripo only has code and x/y/z position,
        # it has no alpha and no normal position
        # so use zeros
        line = '{0} {1:.4f} {2:.4f} {3:.4f} 0 0 0 0 0'.format(*point)
        lines.append(line)
    lines.append('$$$$')
    return '\n'.join(lines) + '\n'


class PharmacophorePointsTable(object):
    """Wrapper around pytables table to store pharmacohpore points

    Args:
        h5file (tables.File):  Pytables hdf5 file object which contains the pharmacophores table
        expectedrows (int): Expected number of pharmacophores.
            Required when hdf5 file is created, helps optimize compression

    Pharmacophore points of a fragment can be retrieved using::

        points = table['frag_id1']

    `points` is a list of points, each point is a tuple with following columns feature type key, x, y and z coordinate.
    The feature type key is defined in FEATURE_TYPES.

    Number of pharmacophore points can be requested using::

        nr_points = len(table)

    To check whether fragment identifier is contained use::

        'frag_id1' in table

    Attributes:
        table (tables.Table): Pytables table with rows of type PharmacophoreRow
    """
    table_name = 'pharmacophores'

    def __init__(self, h5file, expectedrows=0):
        if self.table_name in h5file.root:
            table = h5file.root.__getattr__(self.table_name)
        else:
            table = h5file.create_table('/',
                                        self.table_name,
                                        PharmacophoreRow,
                                        'Pharmacophore points of Kripo sub-pockets',
                                        expectedrows=expectedrows)
            table.cols.frag_id.create_index(filters=PYTABLE_FILTERS)
        self.table = table

    def add_dir(self, startdir):
        """Find \*_pphore.sd.gz \*_pphores.txt file pairs recursively in start directory and add them.

        Args:
            startdir (str): Path to a start directory
        """
        for root, _dirs, files in walk(startdir):
            sdfiles = [path.join(root, file) for file in files if file.endswith('_pphore.sd.gz')]
            for sdfile in sdfiles:
                fragtxtfile = sdfile.replace('_pphore.sd.gz', '_pphores.txt')
                self.add_pocket(sdfile, fragtxtfile)

    def add_pocket(self, sdfile, fragtxtfile):
        points = read_pphore_gzipped_sdfile(sdfile)
        for frag_id, point_ids in read_fragtxtfile(fragtxtfile):
            self.add_fragment(frag_id, point_ids, points)
        self.table.flush()

    def add_fragment(self, frag_id, points_ids, points):
        if frag_id in self:
            raise ValueError("Duplicate key '{0}' found".format(frag_id))
        row = self.table.row
        types = self.table.get_enum('type')
        for i in points_ids:
            row['frag_id'] = frag_id
            point = points[i]
            row['type'] = types[point[0]]
            row['x'] = point[1]
            row['y'] = point[2]
            row['z'] = point[3]
            row.append()

    def __contains__(self, item):
        query = 'frag_id == z'
        binds = {'z': item}
        nr_hits = len(self.table.get_where_list(query, binds))
        return nr_hits > 0

    def __len__(self):
        return len(self.table)

    def __getitem__(self, key):
        points = []
        types = self.table.get_enum('type')
        query = 'frag_id == z'
        binds = {'z': key}
        for row in self.table.where(query, binds):
            points.append((
                types(row['type']),
                row['x'],
                row['y'],
                row['z'],
            ))
        if len(points) == 0:
            raise KeyError(key)
        return points
