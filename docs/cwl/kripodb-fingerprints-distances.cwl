#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"
class: CommandLineTool
label: kripodb-fingerprints-distances
description: Calculate distance matrix from 2 fingerprint files
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - fingerprints
  - distances
inputs:
  - id: fragmentsdbfn
    type: File
    description: Name of fragments db file (only required for hdf5 format)
    inputBinding:
      prefix: --fragmentsdbfn
  - id: ignore_upper_triangle
    type: boolean
    description: Ignore upper triangle
    inputBinding:
      prefix: --ignore_upper_triangle
  - id: fingerprintdb1
    type: File
    inputBinding:
      position: 1
  - id: fingerprintdb2
    type: File
    inputBinding:
      position: 2
  - id: destfile
    type: string
    inputBinding:
      position: 3
outputs:
  - id: dest
    type: File
    outputBinding:
      glob: $(inputs.destfile)
