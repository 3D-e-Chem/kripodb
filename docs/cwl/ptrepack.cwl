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
    description: Source file
    inputBinding:
      position: 1
  - id: destfile_name
    type: string
    description: Destination file
    inputBinding:
      position: 2
outputs:
  - id: destfile
    type: File
    outputBinding:
      glob: $(inputs.destfile_name)
