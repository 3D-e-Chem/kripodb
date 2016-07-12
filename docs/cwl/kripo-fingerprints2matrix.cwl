#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
label: Compute distance matrix of Kripo fingerprints
doc: |
    Workflow that does:
    kripodb fingerprints import fingerprints_01.txt fingerprints_01.db
    kripodb fingerprints distances --fragmentsdbfn fragments.sqlite --ignore_upper_triangle fingerprints_01.db fingerprints_01.db dist_01_01.h5
    kripodb distances freeze dist_01_01.h5 dist_01_01.frozen.h5
    ptrepack --complevel 6 --complib blosc:zlib dist_01_01.frozen.h5 dist_01_01.packedfrozen.h5
inputs:
  fingerprinttxt: File
  fragmentsdb: File
  distmatrixpackedfrozen: string
outputs:
  distmatrixpackedfrozen_file:
    type: File
    outputSource: ptrepack/destfile
steps:
  import-fingerprint:
    run: kripodb-fingerprints-import.cwl
    in:
      fingerprinttxt: fingerprinttxt
      fingerprintdb_name:
        default: fingerprints.sqlite
    out:
      - fingerprintdb
  distance-generate:
    run: kripodb-fingerprints-distances.cwl
    in:
      fragmentsdb: fragmentsdb
      fingerprintdb1: import-fingerprint/fingerprintdb
      fingerprintdb2: import-fingerprint/fingerprintdb
      sparsematrix_name:
        default: sparse_matrix.h5
    out:
      - sparsematrix
  distance-freeze:
    run: kripodb-distances-freeze.cwl
    in:
      sparsematrix: distance-generate/sparsematrix
      frozenmatrix_name:
        default: "frozen_matrix.h5"
    out:
      - frozenmatrix
  ptrepack:
    run: ptrepack.cwl
    in:
      sourcefile: distance-freeze/frozenmatrix
      destfile_name: distmatrixpackedfrozen
      complib:
        default: blosc:zlib
      complevel:
        default: 6
    out:
      - destfile
