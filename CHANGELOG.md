# Change log

## Unreleased

### Changed

* Flag to ignore upper triangle when calculating distances, instead of always ignore (#20)

## 1.4.2 - 3 June 2016

### Changed

* Lower webservice cutoff to 0.45 (#18)

## 1.4.1 - 31 May 2016

### Added

* Webservice online at http://3d-e-chem.vu-compmedchem.nl/kripodb/ui/
* Ignore_upper triangle option in distance import sub command

## 1.4.0 - 3 May 2016

### Changed

* Using nested sub-commands instead of long sub-command. For example `kripodb distmatrix_import` now is `kripodb distances import`

### Added

* Faster distance matrix storage format
* Python3 support (#12)
* Automated build to docker hub.

### Removed

* CLI argument `--precision`

## 1.3.0 - 23 Apr 2016

### Added

* webservice server/client for distance matrix (#16). The CLI and canned commands can now take a local file or a url.

### Fixed

* het_seq_nr contains non-numbers (#15)

## 1.2.5 - 24 Mar 2016

### Fixed

* fpneigh2tsv not available as sub command

## 1.2.4 - 24 Mar 2016

### Added

* Sub command to convert fpneight distance file to tsv.

## 1.2.3 - 1 Mar 2016

### Changed

* Converting distances matrix will load id2label lookup into memory to speed up conversion

## 1.2.2 - 22 Feb 2016

### Added

- Added sub command to read fpneigh formatted distance matrix file (#14)

## 1.2.1 - 12 Feb 2016

### Added

- Added sub commands to read/write distance matrix in tab delimited format (#13)
- Created repo for Knime example and plugin at https://github.com/3D-e-Chem/knime-kripodb (#8)

## 1.2.0 - 11 Feb 2016

### Added

- Prefix to canned fragments lookups (#11)
- PDB meta data to fragments db (#6)
- Limit to distance matrix searches (#5)

### Changed

- Merging of distance matrix files more robust (#10)
- Tanimoto coefficient is rounded up (#7)

## 1.0.0 - 5 Feb 2016

### Added

- Convert fragments shelve to sqlite
- Convert SDF molecules file to sqlite
- Convert Makebits formated file to sqlite
- Create distance matrix using modified tanimoto coefficient in hdf5 format
