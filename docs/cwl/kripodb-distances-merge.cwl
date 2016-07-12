#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: kripodb-distances-merge
doc: Merge distance pairs files into a new file
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - distances
  - merge
inputs:
  ins:
    doc: Input distance pairs files
    type:
      type: array
      items: File
    inputBinding:
      position: 1
  out_name:
    doc: Output distance pairs file
    type: string
    inputBinding:
      position: 2
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out_name)
