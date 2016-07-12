#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: kripodb-fingerprints-import
doc: Convert Kripo fingerprint text file to Kripo internal fingerprint format
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - fingerprints
  - import
inputs:
  fingerprinttxt:
    type: File
    inputBinding:
      position: 1
  fingerprintdb_name:
    type: string
    inputBinding:
      position: 2
outputs:
  fingerprintdb:
    type: File
    outputBinding:
      glob: $(inputs.fingerprintdb_name)
