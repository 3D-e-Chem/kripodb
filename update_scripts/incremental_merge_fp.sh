#!/usr/bin/env bash

merged_fp=kripo_fingerprint_${1}_fp.gz
echo $merged_fp

echo 'MAKEBITS 1.0 574331 BigGrid' | gzip -c > $merged_fp
for x in out.*.fp.gz
do
gnuzip -c $x | tail -n +2 | gzip -c >> $merged_fp
done

rm out.*.fp.gz
