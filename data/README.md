# Example data set

* fragments.sqlite - Fragments sqlite database containing a small number of fragments with their smiles string and molblock.
* fingerprints.sqlite - Fingerprints sqlite database with fingerprint stored as [fastdumped intbitset](http://intbitset.readthedocs.org/en/latest/index.html#intbitset.intbitset.fastdump)
* similarities.h5 - HDF5 file with similarities matrix of fingerprints using modified tanimoto similarity index 

## Creating tiny data set

1. Create fingerprints db with 1000 fingerprints
```
gunzip -c fingerprint01.fp.gz | head -1001 | kripodb fingerprints import - fingerprints.sqlite
```

2. Shrink fragments db to only contain fragments which have a fingerprint
```
cat | sqlite3 fragments.sqlite <<EOF
ATTACH DATABASE 'fingerprints.sqlite' AS fp;
DELETE FROM molecules WHERE frag_id NOT IN (SELECT frag_id FROM fp.bitsets);
DELETE FROM fragments WHERE frag_id NOT IN (SELECT frag_id FROM fp.bitsets);
DELETE FROM pdbs WHERE pdb_code || prot_chain NOT IN (SELECT pdb_code || prot_chain FROM fragments);
VACUUM;
EOF

```

3. Create similarity matrix

```
kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite fingerprints.sqlite fingerprints.sqlite similarities.h5
```

