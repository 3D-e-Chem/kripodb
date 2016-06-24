#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"
class: CommandLineTool
label: kripodb-distances-merge
description: Merge distance pairs files into a new file
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - distances
  - merge
inputs:
  - id: ins
    description: Input distance pairs files
    type: array
      - type: File
    inputBinding:
      position: 1
  - id: out_name
    description: Output distance pairs file
    type: string
    inputBinding:
      position: 2
outputs:
  - id: out
    type: File
    outputBinding:
      glob: $(inputs.out_name)
