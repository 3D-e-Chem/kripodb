#!/bin/sh

kripodb similarities export --no_header --frag1 ../current/similarities.h5 similarities.frag1.txt

LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.frag1.txt -output largevis3.similarities.frag1.txt
LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.frag1.txt -output largevis2.similarities.frag1.txt

rm -f similarities.frag1.txt
