# Kripo DB

[![Build Status](https://travis-ci.org/3D-e-Chem/kripodb.svg?branch=master)](https://travis-ci.org/3D-e-Chem/kripodb)
[![Build status](https://ci.appveyor.com/api/projects/status/diign2fenvai0dst?svg=true)](https://ci.appveyor.com/project/3D-e-Chem/kripodb)
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

* RDKit, http://www.rdkit.org/docs/Install.html, to read SDF files and generate smile strings from molecules
* pip, version 8.0.0 or greater, for wheel support
* git, to clone kripodb repository during installation

```
pip install -U setuptools
pip install numpy
pip install git+https://github.com/3D-e-Chem/kripodb.git
```

# Usage

To see available commands
```
kripodb --help
```

## Create all

Commands to create all data files see [update documentation](docs/data-update.rst).

## Search for most similar fragments

Command to find fragments most similar to `3kxm_K74_frag1` fragment.
```
kripodb similar sim_all.h5 3kxm_K74_frag1 --cutoff 0.45
```

## Create similarity matrix from text files

Commands to create similarity matrix see [update documentation](docs/data-update.rst).

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

# KNIME

The [Knime-KripoDB-example.zip](https://github.com/3D-e-Chem/knime-kripodb/blob/master/examples/Knime-KripoDB-example.zip) file is an example workflow showing how to use KripoDB python package inside KNIME (http://www.knime.org).
It can be run by importing it into KNIME.
Make sure the Python used by KNIME is the same as the Python with kripodb package installed.

The https://github.com/3D-e-Chem/knime-kripodb repo adds KripoDB code templates to KNIME.

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
kripodb serve data/similarities.h5 data/fragments.sqlite data/pharmacophores.h5
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
