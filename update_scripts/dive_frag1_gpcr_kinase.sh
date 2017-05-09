#!/bin/sh

# Requires pdb.gpcr.txt + pdb.gpcr.kinase.txt + pdb.kinase.txt in current working directory

kripodb similarities export --no_header --frag1 --pdb pdb.gpcr.txt ../current/similarities.h5 similarities.frag1.gpcr.txt
kripodb similarities export --no_header --frag1 --pdb pdb.gpcr.kinase.txt ../current/similarities.h5 similarities.frag1.gpcr.kinase.txt
kripodb similarities export --no_header --frag1 --pdb pdb.kinase.txt ../current/similarities.h5 similarities.frag1.kinase.txt

LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.frag1.gpcr.txt -output largevis3.similarities.frag1.gpcr.txt
LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.frag1.gpcr.txt -output largevis2.similarities.frag1.gpcr.txt

LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.frag1.gpcr.kinase.txt -output largevis3.similarities.frag1.gpcr.kinase.txt
LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.frag1.gpcr.kinase.txt -output largevis2.similarities.frag1.gpcr.kinase.txt

LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.frag1.kinase.txt -output largevis3.similarities.frag1.kinase.txt
LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.frag1.kinase.txt -output largevis2.similarities.frag1.kinase.txt
