from os import walk

from rdkit.Chem.rdmolfiles import SDMolSupplier

from db import SqliteDb

symbol2point_type = {
    'He': 'HDON',
    # TODO add other types and correct He
}


def _mol2points(point_ids, mol):
    conf = mol.GetConformer(0)
    points = []
    for point_id in point_ids:
        atom_idx = point_id + 1
        pos = conf.GetAtomPosition(atom_idx)
        symbol = mol.GetAtomWithIdx(atom_idx).GetSymbol()
        point = (
            symbol2point_type[symbol],
            pos.x,
            pos.y,
            pos.z
        )
        points.append(point)
    return points


def _row2pharmacophore(row):
    mol = row['points']
    point_ids = row['point_ids'].split(',')
    return {
        'frag_id': row['frag_id'],
        'points': _mol2points(point_ids, mol)
    }


class PharmacophoresDb(SqliteDb):
    def create_tables(self):
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS ligand_pharmacophores (
            phar_id TEXT PRIMARY KEY,
            points molblockgz
        )''')
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS fragments_pharmacophores (
            frag_id TEXT PRIMARY KEY,
            phar_id TEXT NOT NULL,
            point_ids TEXT NOT NULL,
            FOREIGN KEY(phar_id) REFERENCES ligand_pharmacophores(phar_id)
        )''')
        self.cursor.execute('''CREATE INDEX IF NOT EXISTS fp_fk ON ligand_pharmacophores(phar_id)''')

    def add_dir(self, startdir):
        for root, dirs, files in walk(startdir):
            sdfiles = [file for file in files if file.endswith('.pphores.sd.gz')]
            for sdfile in sdfiles:
                fragtxtfile = sdfile.replace('.sd.gz', '.txt')
                self.add_ligand(sdfile, fragtxtfile)

    def add_ligand_pharmacophores(self, sdfile, fragtxtfile):
        phar_id = self.add_ligand_pharmacophore(sdfile)
        with open(fragtxtfile) as f:
            for line in f:
                # TODO add test with example file to parse
                point_ids = line.split(' ')
                frag_id = point_ids.pop(0)
                self.add_fragment(phar_id, frag_id, point_ids)

    def __getitem__(self, key):
        sql = 'SELECT * FROM fragments_pharmacophores JOIN ligand_pharmacophores USING (phar_id) WHERE frag_id=?'
        self.cursor.execute(sql, (key,))
        row = self.cursor.fetchone()

        if row is None:
            raise KeyError(key)

        return _row2pharmacophore(row)

    def __len__(self):
        self.cursor.execute('SELECT count(*) FROM fragments_pharmacophores')
        row = self.cursor.fetchone()
        return row[0]

    def add_ligand_pharmacophore(self, sdfile):
        suppl = SDMolSupplier(sdfile)
        mol = suppl.next()
        sql = '''INSERT OR REPLACE INTO ligand_pharmacophores (phar_id, points) VALUES (?, ?)'''
        phar_id =mol.GetProp('_Name')  # TODO find way to get frag_id similiar to the ones in frag dbs
        self.cursor.execute(sql, (phar_id, mol,))
        self.connection.commit()
        return phar_id

    def add_fragment(self, phar_id, frag_id, point_ids):
        sql = '''INSERT OR REPLACE INTO fragments_pharmacophores (frag_id, phar_id, points) VALUES (?, ?, ?)'''
        self.cursor.execute(sql, (frag_id, phar_id, ','.join(point_ids),))
        self.connection.commit()

