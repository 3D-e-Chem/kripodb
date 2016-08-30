# Kripo DB

[![Build Status](https://travis-ci.org/3D-e-Chem/kripodb.svg?branch=master)](https://travis-ci.org/3D-e-Chem/kripodb)
[![Codacy Grade Badge ](https://api.codacy.com/project/badge/Grade/4878758675a0402bb75019672fa6e45c)](https://www.codacy.com/app/3D-e-Chem/kripodb?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/kripodb&amp;utm_campaign=Badge_Grade)
[![Codacy Coverage Badge](https://api.codacy.com/project/badge/Coverage/4878758675a0402bb75019672fa6e45c)](https://www.codacy.com/app/3D-e-Chem/kripodb?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/kripodb&amp;utm_campaign=Badge_Coverage)
[![DOI](https://zenodo.org/badge/19641/3D-e-Chem/kripodb.svg)](https://zenodo.org/badge/latestdoi/19641/3D-e-Chem/kripodb)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/3dechem/kripodb/)
[![Documentation Status](https://readthedocs.org/projects/kripodb/badge/?version=latest)](http://kripodb.readthedocs.io/en/latest/?badge=latest)

Library to interact with Kripo fragment, fingerprint and similarity data files.

KRIPO stands for Key Representation of Interaction in POckets, see [reference](http://dx.doi.org/10.1186/1758-2946-6-S1-O26) for more information.

# Glossary

* Pocket, binding site of the ligand in the protein of a crystal structure
* Fragment, part of the ligand
* Subpocket, part of the protein pocket which binds with the fragment
* Fingerprint, fingerprint of structure-based pharmacophore of subpocket
* Similarity matrix, similarities between all fingerprint pairs calculated using the modified tanimoto similarity index
* Kripo fragment identifier, used as identifier for fragment, subpocket and fingerprint

# Install

Requirements:

* rdkit, http://rdkit.org, to read SDF files and generate smile strings from molecules
* libhdf5 headers, to read/write similarity matrix in hdf5 format

```
pip install -U setuptools
pip install numpy
python setup.py install
```

# Usage

To see available commands
```
kripodb --help
```

## Create all

Commands to create all data files
```
kripodb fragments shelve fragments.shelve fragments.sqlite
kripodb fragments sdf fragment??.sdf fragments.sqlite
kripodb fragments pdb fragments.sqlite
kripodb fingerprints import 01.fp 01.fp.db
kripodb fingerprints import 02.fp 02.fp.db
kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite --ignore_upper_triangle 01.fp.db 01.fp.db sim_01_01.h5
kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite --ignore_upper_triangle 02.fp.db 02.fp.db sim_02_02.h5
kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite 01.fp.db 02.fp.db sim_01_02.h5
kripodb similarities merge sim_*_*.h5  sim_all.h5
kripodb similarities freeze sim_all.h5 sim_all.frozen.h5
# Make froze similarity matrix smaller, by using slower compression
ptrepack --complevel 6 --complib blosc:zlib sim_all.frozen.h5 sim_all.packedfrozen.h5
rm sim_all.frozen.h5
kripodb similarities serve sim_all.packedfrozen.h5
```

## Search for most similar fragments

Command to find fragments most similar to `3kxm_K74_frag1` fragment.
```
kripodb similar sim_all.h5 3kxm_K74_frag1 --cutoff 0.45
```

## Create similarity matrix from text files

Input files `sim_??_??.txt.gz` looks like:
```
Compounds similar to 2xry_FAD_frag4:
2xry_FAD_frag4   1.0000
3cvv_FAD_frag3   0.5600
```

To create a single similarity matrix from multiple text files:
```
gunzip -c sim_01_01.txt.gz | kripodb similarities import --ignore_upper_triangle - fragments.sqlite sim_01_01.h5
gunzip -c sim_01_02.txt.gz | kripodb similarities import - fragments.sqlite sim_01_02.h5
gunzip -c sim_02_02.txt.gz | kripodb similarities import --ignore_upper_triangle - fragments.sqlite sim_02_02.h5
kripodb similarities merge sim_??_??.h5 sim_all.h5
```

The `--ignore_upper_triangle` flag is used to prevent scores corruption when freezing similarity matrix.

# Data sets

## Example

An example data set included in the [data/](data/) directory of this repo. See [data/README.md](data/README.md) for more information.

## GPCR

All fragments based on GPCR proteins compared with all proteins in PDB.

* kripo.gpcrandhits.sqlite - Fragments sqlite database
* kripo.gpcr.h5 - HDF5 file with similarity matrix

The data set has been published at [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.50835.svg)](http://dx.doi.org/10.5281/zenodo.50835)

## Protein Data Bank

All fragments form all proteins-ligand complexes in PDB compared with all.

* Fragments sqlite database - Download from http://3d-e-chem.vu-compmedchem.nl/kripodb/fragments.sqlite
* Similarity matrix - Can be queried on webservice at http://3d-e-chem.vu-compmedchem.nl/kripodb. For build instructions see http://kripodb.readthedocs.io/en/latest/data-update.html
* Fragment fingerprints - See http://kripodb.readthedocs.io/en/latest/data-update.html for instructions how to convert to a similarity matrix

A data set with PDB entries till 23 December 2015 has been published at [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.55254.svg)](http://dx.doi.org/10.5281/zenodo.55254)

# Knime

The [Knime-KripoDB-example.zip](https://github.com/3D-e-Chem/knime-kripodb/blob/master/examples/Knime-KripoDB-example.zip) file is an example workflow showing how to use KripoDB python package inside Knime (http://www.knime.org).
It can be run by importing it into Knime.
Make sure the Python used by Knime is the same as the Python with kripodb package installed.

The https://github.com/3D-e-Chem/knime-kripodb repo adds KripoDB code templates to Knime.

# Development of KripoDB

Install the development deps with:
```
pip install -r requirements.txt
```

# Docker

## Create image

```
docker build -t 3dechem/kripodb .
```

## Run container

Show the kripodb help with
```
docker run --rm 3dechem/kripodb kripodb --help
```

To calculate the mean bit density of the fingerprints in the `fingerprints.sqlite` file in the current working directory use following command.
```
docker run --rm -u $UID -v $PWD:/data 3dechem/kripodb kripodb meanbitdensity /data/fingerprints.sqlite
```

# Web service

The Kripo data files can be queried using a web service.

Start webservice with:
```
kripodb serve data/similarities.h5 data/fragments.sqlite
```
It will print the urls for the swagger spec and UI.

Note! The webservice returns a limited amount of results. To get all results use local files.

On http://3d-e-chem.vu-compmedchem.nl/kripodb/ui/ there is a KripoDB webservice with the full PDB fragment all vs all matrix.

# Documentation

API and data update pipeline documentation can be found at http://kripodb.readthedocs.io/en/latest/.

# Reference

KRIPO – a structure-based pharmacophores approach explains polypharmacological effects;
Tina Ritschel, Tom JJ Schirris, and Frans GM Russel; J Cheminform. 2014; 6(Suppl 1): O26;
Published online 2014 Mar 11; http://dx.doi.org/10.1186/1758-2946-6-S1-O26
