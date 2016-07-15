#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
label: Compute similarity matrix of Kripo fingerprints
doc: |
    Workflow that does:
    kripodb fingerprints import fingerprints_01.txt fingerprints_01.db
    kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite --ignore_upper_triangle fingerprints_01.db fingerprints_01.db sim_01_01.h5
    kripodb similarities freeze sim_01_01.h5 sim_01_01.frozen.h5
    ptrepack --complevel 6 --complib blosc:zlib sim_01_01.frozen.h5 sim_01_01.packedfrozen.h5
inputs:
  fingerprinttxt: File
  fragmentsdb: File
  sim_matrixpackedfrozen: string
outputs:
  sim_matrixpackedfrozen_file:
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
  similarity-generate:
    run: kripodb-fingerprints-similarities.cwl
    in:
      fragmentsdb: fragmentsdb
      fingerprintdb1: import-fingerprint/fingerprintdb
      fingerprintdb2: import-fingerprint/fingerprintdb
      sparsematrix_name:
        default: sparse_matrix.h5
    out:
      - sparsematrix
  similarity-freeze:
    run: kripodb-similarities-freeze.cwl
    in:
      sparsematrix: similarity-generate/sparsematrix
      frozenmatrix_name:
        default: "frozen_matrix.h5"
    out:
      - frozenmatrix
  ptrepack:
    run: ptrepack.cwl
    in:
      sourcefile: similarity-freeze/frozenmatrix
      destfile_name: sim_matrixpackedfrozen
      complib:
        default: blosc:zlib
      complevel:
        default: 6
    out:
      - destfile
