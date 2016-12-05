# Change log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
Formatted as described on http://keepachangelog.com/.

## Unreleased

### Added

- Published documentation on http://kripodb.readthedocs.io
- Documented update pipeline
- Documented command line interface in Sphinx
- Retrieve fragments from webservice based on fragment id and pdb code

### Changed

- merge similarity pairs files in chunks instead of loading whole source file in memory
- canned fragments_by_* methods can use local file or webservice
- error when duplicate fragment insert is performed
- Renamed `kripodb similarities serve` to `kripodb serve`, as it now also serves the fragments

### Removed

- no longer create indices for similarity pairs file, querying is done on dense matrix

## [2.0.0] - 2016-07-14

### Changed

- Renamed distance to similarity (#21)
- Flag to ignore upper triangle when calculating distances, instead of always ignore (#20)

## [1.4.2] - 2016-06-03

### Changed

- Lower webservice cutoff to 0.45 (#18)

## [1.4.1] - 2016-05-31

### Added

- Webservice online at http://3d-e-chem.vu-compmedchem.nl/kripodb/ui/
- Ignore_upper triangle option in distance import sub command

## [1.4.0] - 2016-05-03

### Changed

- Using nested sub-commands instead of long sub-command. For example `kripodb distmatrix_import` now is `kripodb distances import`

### Added

- Faster distance matrix storage format
- Python3 support (#12)
- Automated build to docker hub.

### Removed

- CLI argument `--precision`

## [1.3.0] - 2016-04-23

### Added

- webservice server/client for distance matrix (#16). The CLI and canned commands can now take a local file or a url.

### Fixed

- het_seq_nr contains non-numbers (#15)

## [1.2.5] - 2016-03-24

### Fixed

- fpneigh2tsv not available as sub command

## [1.2.4] - 2016-03-24

### Added

- Sub command to convert fpneight distance file to tsv.

## [1.2.3] - 2016-03-01

### Changed

- Converting distances matrix will load id2label lookup into memory to speed up conversion

## [1.2.2] - 2016-02-22

### Added

- Added sub command to read fpneigh formatted distance matrix file (#14)

## [1.2.1] - 2016-02-12

### Added

- Added sub commands to read/write distance matrix in tab delimited format (#13)
- Created repo for Knime example and plugin at https://github.com/3D-e-Chem/knime-kripodb (#8)

## [1.2.0] - 2016-02-11

### Added

- Prefix to canned fragments lookups (#11)
- PDB meta data to fragments db (#6)
- Limit to distance matrix searches (#5)

### Changed

- Merging of distance matrix files more robust (#10)
- Tanimoto coefficient is rounded up (#7)

## [1.0.0] - 2016-02-05

### Added

- Convert fragments shelve to sqlite
- Convert SDF molecules file to sqlite
- Convert Makebits formated file to sqlite
- Create distance matrix using modified tanimoto coefficient in hdf5 format
