#!/usr/bin/env bash

for z in $(ls *fp.gz)
do
srun -Jgunzip$(basename $z .fp.gz) -n 1 /bin/sh -c "gunzip -c $z > $(basename $z .gz)" &
done
wait