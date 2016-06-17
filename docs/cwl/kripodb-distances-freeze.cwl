#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"
class: CommandLineTool
label: kripodb-distances-freeze
description: Convert kripo sparse distance matrix to a dense compressed version
requirements:
  - class: DockerRequirement
    dockerPull: 3dechem/kripodb
baseCommand: kripodb
arguments:
  - distances
  - freeze
inputs:
  - id: frame_size
    type: int
    default: 100000000
    description: Size of frame
    inputBinding:
      prefix: --frame_size
  - id: memory_cache
    type: int
    default: 1
    description: Memory cache in Gigabytes
    inputBinding:
      prefix: --memory
  - id: limit
    type: int
    description: Number of pairs to copy, None for no limit
    inputBinding:
      prefix: --limit
  - id: single_sided
    type: boolean
    description: Store half matrix
    inputBinding:
      prefix: --single_sided
  - id: sourcefile
    type: File
    inputBinding:
      position: 1
  - id: frozenmatrix
    type: string
    inputBinding:
      position: 2
outputs:
  - id: dest
    type: File
    outputBinding:
      glob: $(inputs.destfile)
