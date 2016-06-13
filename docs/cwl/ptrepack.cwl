#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"
class: CommandLineTool
label: ptrepack
description: Recompress PyTable formatted file
baseCommand: ptrepack
inputs:
  - id: complib
    type:
      name: complibs
      type: enum
      symbols:
        - zlib
        - lzo
        - bzip2
        - blosc
        - "blosc:blosclz"
        - "blosc:lz4"
        - "blosc:lz4hc"
        - "blosc:snappy"
        - "blosc:zlib"
    description: compression library
    default: zlib
    inputBinding:
      prefix: --complib
  - id: complevel
    type: int
    default: 0
    description: Compression level
    inputBinding:
      prefix: --complevel
  - id: sourcefile
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
