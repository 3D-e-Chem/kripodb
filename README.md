# Kripo DB

[![Build Status](https://travis-ci.org/3D-e-Chem/kripodb.svg?branch=master)](https://travis-ci.org/3D-e-Chem/kripodb)
[![Codacy Badge](https://api.codacy.com/project/badge/grade/4878758675a0402bb75019672fa6e45c)](https://www.codacy.com/app/NLeSC/kripodb)
[![Codacy Badge](https://api.codacy.com/project/badge/coverage/4878758675a0402bb75019672fa6e45c)](https://www.codacy.com/app/NLeSC/kripodb)

Library to interact with Kripo fragment, fingerprint and similarity data files.

KRIPO stands for Key Representation of Interaction in POckets, see [reference](http://dx.doi.org/10.1186/1758-2946-6-S1-O26) for more information.

# Install

Requirements:

* rdkit, http://rdkit.org, to read SDF files and generate smile strings from molecules

```
pip install -U setuptools
pip install numpy
python setup.py develop
```

# Usage

To see available commands
```
kripodb --help
```

## Create all

Commands to create all data files
```
kripodb shelve2fragmentsdb fragments.shelve fragments.sqlite
kripodb sdf2fragmentsdb fragment??.sdf fragments.sqlite
kripodb makebits2fingerprintsdb 01.fp 01.fp.db
kripodb makebits2fingerprintsdb 02.fp 02.fp.db
kripodb pairs --fragmentsdbfn fragments.sqlite 01.fp.db 01.fp.db dist_01_01.h5
kripodb pairs --fragmentsdbfn fragments.sqlite 01.fp.db 02.fp.db dist_01_02.h5
kripodb pairs --fragmentsdbfn fragments.sqlite 02.fp.db 01.fp.db dist_02_01.h5
kripodb pairs --fragmentsdbfn fragments.sqlite 02.fp.db 02.fp.db dist_02_02.h5
kripodb mergepairs dist_*_*.h5  dist_all.h5
```

## Search for most similar fragments

Command to find fragments most similar to `3kxm_K74_frag1` fragment.
```
kripodb similar dist_all.h5 3kxm_K74_frag1 --cutoff 0.45
```

# Docker

## Create image

```
docker build -t nlesc/kripodb .
```

## Run container

To calculate the mean bit density of the fingerprints in the `fingerprints.sqlite` file in the current working directory use following command.
```
docker run --rm -u $UID -v $PWD:/data nlesc/kripodb kripodb meanbitdensity /data/fingerprints.sqlite
```

# Reference

KRIPO – a structure-based pharmacophores approach explains polypharmacological effects;
Tina Ritschel, Tom JJ Schirris, and Frans GM Russel; J Cheminform. 2014; 6(Suppl 1): O26;
Published online 2014 Mar 11; http://dx.doi.org/10.1186/1758-2946-6-S1-O26