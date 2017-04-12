#!/usr/bin/env bash
kripodb similarities merge similarities.new__*[0-9].h5 similarities.new__existing.h5 && \
rm similarities.new__*[0-9].h5

# Compact the fingerprint file (makebits ascii format)
gzip out.fp
mv out.fp.gz out.$(date +%Y%V).fp.gz

# Add new similarities to existing similarities file
kripodb similarities merge ../current/similarities.h5 similarities.new__existing.h5 similarities.new__new.h5 similarities.h5 && \
rm similarities.new__existing.h5 similarities.new__new.h5