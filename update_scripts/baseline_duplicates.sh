#!/usr/bin/env bash

# find duplicate ids
for x in *.fp.gz
do
export x
gunzip -c $x | perl -p -e 's/^(.*?) .*$/$1,$ENV{"x"}/' >> ids.csv
done

python $SCRIPTS/fix_dup_ids.py
