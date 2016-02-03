# Kripo DB

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
