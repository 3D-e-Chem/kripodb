#!/usr/bin/env bash

# Compute similarities and write as h5 files
nrrows=10000000
for x in $(ls *.fp)
do
for y in $(ls *.fp)
do
if [ "$x" = "$y" ]
then
srun -J $x$y -n 1 /bin/sh -c "fpneigh -m Mod_Tanimoto=0.01 -d 0.45 -q $x $y | kripodb similarities import --nrrows $nrrows --ignore_upper_triangle - fragments.sqlite similarities.$(basename $x .fp)__$(basename $y .fp).h5" &
elif [[ $x < $y ]]
then
srun -J $x$y -n 1 /bin/sh -c "fpneigh -m Mod_Tanimoto=0.01 -d 0.45 -q $x $y | kripodb similarities import --nrrows $nrrows - fragments.sqlite similarities.$(basename $x .fp)__$(basename $y .fp).h5" &
fi
done
done
wait
