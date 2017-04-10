#!/usr/bin/env bash

kripodb similarities freeze -f 400000000 similarities.h5 similarities.frozen.h5
ptrepack --complevel 6 --complib blosc:zlib similarities.frozen.h5 similarities.packedfrozen.h5 && rm similarities.frozen.h5