# Kripo DB

[![Build Status](https://travis-ci.org/3D-e-Chem/kripodb.svg?branch=master)](https://travis-ci.org/3D-e-Chem/kripodb)
[![Codacy Badge](https://api.codacy.com/project/badge/grade/4878758675a0402bb75019672fa6e45c)](https://www.codacy.com/app/NLeSC/kripodb)
[![Codacy Badge](https://api.codacy.com/project/badge/coverage/4878758675a0402bb75019672fa6e45c)](https://www.codacy.com/app/NLeSC/kripodb)

Library to interact with Kripo fragment, fingerprint and similarity data files.

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
