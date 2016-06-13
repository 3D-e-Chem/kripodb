# Workflow that does:
# kripodb fingerprints import fingerprints_01.txt fingerprints_01.db
# kripodb fingerprints distances --fragmentsdbfn fragments.sqlite --ignore_upper_triangle fingerprints_01.db fingerprints_01.db dist_01_01.h5
# kripodb distances freeze dist_01_01.h5 dist_01_01.frozen.h5
# ptrepack --complevel 0 --complib blosc distances.blosc6.h5

cwlVersion: cwl:draft-3
class: Workflow
inputs:
  - id: fingerprinttxt
    type: File
  - id: fragmentsdb
    type: File
  - id: frozenmatrix
    type: str
outputs:
  - id: frozenmatrixout
    type: File
    source: "#kripodb-distances-freeze/frozenmatrix"
steps:
  - id: import-fingerprint
    run: kripodb-fingerprints-import.cwl
    inputs:
      - id: fingerprinttxt
        source: "#fingerprinttxt"
      - id: fingerprintdb
        source: "#fingerpintdb"
  - id: distance-generate
    run: kripodb-fingerprints-distances.cwl
    inputs:
      - id: fragmentsdb
        source: "#fingerpintdb"
      - id: fingerprintdb1
        source: "#fingerpintdb"
      - id: fingerprintdb2
        source: "#fingerpintdb"
