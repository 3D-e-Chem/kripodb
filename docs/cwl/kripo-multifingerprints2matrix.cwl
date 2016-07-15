#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
label: Compute similarity matrix of multiple Kripo fingerprint files
doc: |
    Workflow that does:
    kripodb fingerprints import 01.fp 01.fp.db
    kripodb fingerprints import 02.fp 02.fp.db
    kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite --ignore_upper_triangle 01.fp.db 01.fp.db sim_01_01.h5
    kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite --ignore_upper_triangle 02.fp.db 02.fp.db sim_02_02.h5
    kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite 01.fp.db 02.fp.db sim_01_02.h5
    kripodb similarities merge sim_*_*.h5  sim_all.h5
    kripodb similarities freeze sim_all.h5 sim_all.frozen.h5
    ptrepack --complevel 6 --complib blosc:zlib sim_all.frozen.h5 sim_all.packedfrozen.h5
inputs:
  fingerprinttxt:
    type:
      - type: array
        items: File
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
      destfile_name: distmatrixpackedfrozen
      complib:
        default: blosc:zlib
      complevel:
        default: 6
    out:
      - destfile
