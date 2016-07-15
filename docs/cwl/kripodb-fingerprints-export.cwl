#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: kripodb-fingerprints-export
doc: Convert Kripo internal fingerprint format to Kripo fingerprint text file
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - fingerprints
  - export
inputs:
  fingerprintdb:
    type: File
    inputBinding:
      position: 1
  fingerprinttxt_name:
    type: string
    inputBinding:
      position: 2
  start:
    type: int
    inputBinding:
      prefix: --start
  stop:
    type: int
    inputBinding:
      prefix: --stop
outputs:
  fingerprinttxt:
    type: File
    outputBinding:
      glob: $(inputs.fingerprinttxt_name)
