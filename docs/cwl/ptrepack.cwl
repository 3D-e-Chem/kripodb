#!/usr/bin/env cwl-runner
cwlVersion: "cwl:draft-3"
class: CommandLineTool
label: ptrepack
doc: Recompress PyTable formatted file
baseCommand: ptrepack
requirements:
  - class: DockerRequirement
    dockerPull: rcsa/python2-hdf5
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
    doc: compression library
    default: zlib
    inputBinding:
      prefix: --complib
  - id: complevel
    type: int
    default: 0
    doc: Compression level
    inputBinding:
      prefix: --complevel
  - id: sourcefile
    type: File
    doc: Source file
    inputBinding:
      position: 1
  - id: destfile_name
    type: string
    doc: Destination file
    inputBinding:
      position: 2
outputs:
  - id: destfile
    type: File
    outputBinding:
      glob: $(inputs.destfile_name)
