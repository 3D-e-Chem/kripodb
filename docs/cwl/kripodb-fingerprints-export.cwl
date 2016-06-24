#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"
class: CommandLineTool
label: kripodb-fingerprints-export
description: Convert Kripo internal fingerprint format to Kripo fingerprint text file
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - fingerprints
  - export
inputs:
  - id: fingerprintdb
    type: File
    inputBinding:
      position: 1
  - id: fingerprinttxt_name
    type: string
    inputBinding:
      position: 2
outputs:
  - id: fingerprinttxt
    type: File
    outputBinding:
      glob: $(inputs.fingerprinttxt_name)
