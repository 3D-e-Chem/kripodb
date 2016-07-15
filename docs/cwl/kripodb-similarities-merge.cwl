#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: kripodb-similarities-merge
doc: Merge similarity pairs files into a new file
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - similarities
  - merge
inputs:
  ins:
    doc: Input similarity pairs files
    type:
      type: array
      items: File
    inputBinding:
      position: 1
  out_name:
    doc: Output similarity pairs file
    type: string
    inputBinding:
      position: 2
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out_name)
