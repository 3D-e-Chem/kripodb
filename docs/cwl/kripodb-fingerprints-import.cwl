#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"
class: CommandLineTool
label: kripodb-fingerprints-import
description: Convert Kripo fingerprint text file to Kripo internal fingerprint format
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - fingerprints
  - import
inputs:
  - id: fingerprinttxt
    type: File
    inputBinding:
      position: 1
  - id: destfile
    type: string
    inputBinding:
      position: 2
outputs:
  - id: dest
    type: File
    outputBinding:
      glob: $(inputs.destfile)
