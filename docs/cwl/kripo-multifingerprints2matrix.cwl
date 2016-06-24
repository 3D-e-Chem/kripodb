#!/usr/bin/env cwl-runner
cwlVersion: cwl:draft-3
class: Workflow
label: Compute distance matrix of multiple Kripo fingerprint files
description: |
    Workflow that does:
    kripodb fingerprints import 01.fp 01.fp.db
    kripodb fingerprints import 02.fp 02.fp.db
    kripodb fingerprints distances --fragmentsdbfn fragments.sqlite --ignore_upper_triangle 01.fp.db 01.fp.db dist_01_01.h5
    kripodb fingerprints distances --fragmentsdbfn fragments.sqlite --ignore_upper_triangle 02.fp.db 02.fp.db dist_02_02.h5
    kripodb fingerprints distances --fragmentsdbfn fragments.sqlite 01.fp.db 02.fp.db dist_01_02.h5
    kripodb distances merge dist_*_*.h5  dist_all.h5
    kripodb distances freeze dist_all.h5 dist_all.frozen.h5
    ptrepack --complevel 6 --complib blosc:zlib dist_all.frozen.h5 dist_all.packedfrozen.h5
inputs:
  - id: fingerprinttxt
    type:
      - type: array
        items: File
  - id: fragmentsdb
    type: File
  - id: distmatrixpackedfrozen
    type: string
outputs:
  - id: distmatrixpackedfrozen_file
    type: File
    source: "#ptrepack/destfile"
steps:
  - id: import-fingerprint
    run: kripodb-fingerprints-import.cwl
    inputs:
      - id: fingerprinttxt
        source: "#fingerprinttxt"
      - id: fingerprintdb_name
        default: fingerprints.sqlite
    outputs:
      - id: fingerprintdb
  - id: distance-generate
    run: kripodb-fingerprints-distances.cwl
    inputs:
      - id: fragmentsdb
        source: "#fragmentsdb"
      - id: fingerprintdb1
        source: "#import-fingerprint/fingerprintdb"
      - id: fingerprintdb2
        source: "#import-fingerprint/fingerprintdb"
      - id: sparsematrix_name
        default: sparse_matrix.h5
    outputs:
      - id: sparsematrix
  - id: distance-freeze
    run: kripodb-distances-freeze.cwl
    inputs:
      - id: sparsematrix
        source: "#distance-generate/sparsematrix"
      - id: frozenmatrix_name
        default: "frozen_matrix.h5"
    outputs:
      - id: frozenmatrix
  - id: ptrepack
    run: ptrepack.cwl
    inputs:
      - id: sourcefile
        source: "#distance-freeze/frozenmatrix"
      - id: destfile_name
        source: "#distmatrixpackedfrozen"
      - id: complib
        default: blosc:zlib
      - id: complevel
        default: 6
    outputs:
      - id: destfile
