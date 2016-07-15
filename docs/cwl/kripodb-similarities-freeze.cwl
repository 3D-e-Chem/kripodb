#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: kripodb-similarities-freeze
doc: Convert kripo sparse similarity matrix to a dense compressed version
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - similarities
  - freeze
inputs:
  frame_size:
    type: int
    default: 100000000
    doc: Size of frame
    inputBinding:
      prefix: --frame_size
  memory_cache:
    type: int
    default: 1
    doc: Memory cache in Gigabytes
    inputBinding:
      prefix: --memory
  limit:
    type: ["null", int]
    doc: Number of pairs to copy, None for no limit
    inputBinding:
      prefix: --limit
  single_sided:
    type: boolean
    doc: Store half matrix
    default: false
    inputBinding:
      prefix: --single_sided
  sparsematrix:
    doc: Input pairs file
    type: File
    inputBinding:
      position: 1
  frozenmatrix_name:
    doc: Output array file, file is overwritten
    type: string
    inputBinding:
      position: 2
outputs:
  frozenmatrix:
    type: File
    outputBinding:
      glob: $(inputs.frozenmatrix_name)
