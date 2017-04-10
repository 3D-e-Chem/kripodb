#!/usr/bin/env bash

for x in $(ls *.fp)
do
srun -n 1 gzip $x &
done
wait