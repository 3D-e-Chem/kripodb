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
  - id: fragmentsdb
    type: File
    description: Name of fragments db file (only required for hdf5 format)
    inputBinding:
      prefix: --fragmentsdbfn
  - id: ignore_upper_triangle
    type: boolean
    description: Ignore upper triangle
    default: false
    inputBinding:
      prefix: --ignore_upper_triangle
  - id: nomemory
    type: boolean
    description: Do not store query fingerprints in memory
    default: false
    inputBinding:
      prefix: --nomemory
  - id: cutoff
    type: float
    description: Set Tanimoto cutoff
    default: 0.45
    inputBinding:
      prefix: --cutoff
  - id: mean_onbit_density
    type: float
    description: Mean on bit density
    default: 0.01
    inputBinding:
      prefix: --mean_onbit_density
  - id: out_format
    type:
      type: enum
      name: out formats
      symbols:
        - tsv
        - hdf5
    description: Format of output
    default: hdf5
    inputBinding:
      prefix: --out_format
  - id: fingerprintdb1
    type: File
    inputBinding:
      position: 1
  - id: fingerprintdb2
    type: File
    inputBinding:
      position: 2
  - id: sparsematrix_name
    type: string
    inputBinding:
      position: 3
outputs:
  - id: sparsematrix
    type: File
    outputBinding:
      glob: $(inputs.sparsematrix_name)
