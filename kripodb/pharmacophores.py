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
FEATURE_TYPE_KEYS = [r['key'] for r in FEATURE_TYPES]
FEATURE_TYPE_ATOM2KEY = {r['element']: r['key'] for r in FEATURE_TYPES}

PYTABLE_FILTERS = tables.Filters(complevel=6, complib='blosc')


class PharmacophoreRow(tables.IsDescription):
    """Table description for similarity pair"""
    frag_id = tables.StringCol(16)
    type = tables.EnumCol(FEATURE_TYPE_KEYS, FEATURE_TYPE_KEYS[0], base='uint8')
    x = tables.Float32Col()
    y = tables.Float32Col()
    z = tables.Float32Col()


class PharmacophoresDb(object):
    def __init__(self, filename, mode='r', expectedrows=0, **kwargs):
        self.h5file = tables.open_file(filename, mode, filters=PYTABLE_FILTERS, **kwargs)
        self.points = PharmacophorePointsTable(self.h5file, expectedrows)

    def close(self):
        self.h5file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def add_dir(self, startdir):
        self.points.add_dir(startdir)

    def __getitem__(self, item):
        return self.points[item]


def read_pphore_gzipped_sdfile(sdfile):
    """Read a gzipped sdfile which contains pharmacophore points as atoms

    Args:
        sdfile (string): Path to filename

    Returns: List of Pharmacophore points

    """
    with gzip.open(sdfile) as gzfile:
        return read_pphore_sdfile(gzfile)


def read_pphore_sdfile(sdfile):
    """Read a gzipped sdfile which contains pharmacophore points as atoms

    Args:
        sdfile (file): File object with sdfile contents

    Returns: List of pharmacophore points

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


class PharmacophorePointsTable(object):
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
        for root, dirs, files in walk(startdir):
            sdfiles = [path.join(root, file) for file in files if file.endswith('_pphore.sd.gz')]
            for sdfile in sdfiles:
                fragtxtfile = sdfile.replace('_pphore.sd.gz', '_pphores.txt')
                self.add_pocket(sdfile, fragtxtfile)

    def add_pocket(self, sdfile, fragtxtfile):
        points = read_pphore_gzipped_sdfile(sdfile)
        with open(fragtxtfile) as f:
            for line in f:
                point_ids = line.split(' ')
                frag_id = point_ids.pop(0)
                point_ids = [int(r) - 1 for r in point_ids]
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
