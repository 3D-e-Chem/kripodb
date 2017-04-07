#!/usr/bin/env python

from __future__ import print_function
import os

# Find duplicates
ids = {}
dups = {}
with open('ids.csv') as f:
    for l in f:
        (k, v) = l.strip().split(',')
        if k == 'MAKEBITS':
            continue
        if k in ids:
            if v in dups:
                dups[v].add(k)
            else:
                dups[v] = set([k])
        else:
            ids[k] = v

if len(dups) == 0:
    print('No duplicates found')

# Remove duplicates
for fn in dups:
    print('Removing ' + ','.join(dups[fn]) + ' from ' + fn)
    fn_dudup = fn + '.dedup'
    with open(fn) as fin:
        with open(fn_dudup, 'w') as fout:
            for l in fin:
                k = l.strip().split(' ')[0]
                if k not in dups[fn]:
                    fout.write(l)

    os.rename(fn, fn + '.dups')
    os.rename(fn_dudup, fn)
