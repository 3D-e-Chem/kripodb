#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: kripodb-fingerprints-distances
doc: Calculate distance matrix from 2 fingerprint files
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - fingerprints
  - distances
inputs:
  fragmentsdb:
    type: File
    doc: Name of fragments db file (only required for hdf5 format)
    inputBinding:
      prefix: --fragmentsdbfn
  ignore_upper_triangle:
    type: boolean
    doc: Ignore upper triangle
    default: false
    inputBinding:
      prefix: --ignore_upper_triangle
  nomemory:
    type: boolean
    doc: Do not store query fingerprints in memory
    default: false
    inputBinding:
      prefix: --nomemory
  cutoff:
    type: float
    doc: Set Tanimoto cutoff
    default: 0.45
    inputBinding:
      prefix: --cutoff
  mean_onbit_density:
    type: float
    doc: Mean on bit density
    default: 0.01
    inputBinding:
      prefix: --mean_onbit_density
  out_format:
    type:
      type: enum
      name: out formats
      symbols:
        - tsv
        - hdf5
    doc: Format of output
    default: hdf5
    inputBinding:
      prefix: --out_format
  fingerprintdb1:
    type: File
    inputBinding:
      position: 1
  fingerprintdb2:
    type: File
    inputBinding:
      position: 2
  sparsematrix_name:
    type: string
    inputBinding:
      position: 3
outputs:
  sparsematrix:
    type: File
    outputBinding:
      glob: $(inputs.sparsematrix_name)
