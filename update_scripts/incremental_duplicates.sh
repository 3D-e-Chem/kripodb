#!/usr/bin/env bash

# find duplicate ids
for x in ../current/*fp.gz
do
export x
gunzip -c $x | perl -p -e 's/^(.*?) .*$/$1,$ENV{"x"}/' >> ids.csv
done
perl -p -e 's/^(.*?) .*$/$1,out.fp/' out.fp >> ids.csv

python $SCRIPTS/fix_dup_ids.py