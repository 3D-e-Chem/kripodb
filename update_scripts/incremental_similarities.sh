#!/usr/bin/env bash

# Compute similarities against itself
# to prevent duplicates similarities of a chunk against itself should ignore the upper triangle.
nrrows=10000000
srun -J new__new -n 1 /bin/sh -c "fpneigh -m Mod_Tanimoto=0.01 -d 0.45 -q out.fp out.fp > similarities.new_new.txt && tail similarities.new_new.txt && kripodb similarities import --nrrows $nrrows --ignore_upper_triangle similarities.new_new.txt fragments.sqlite similarities.new__new.h5 && rm similarities.new_new.txt" &

# Compute similarities against existing fingerprint chunks
for x in `ls ../current/*fp.gz`
do
srun -J new__$x -n 1 /bin/sh -c "gunzip -c $x > $(basename $x .gz) && fpneigh -m Mod_Tanimoto=0.01 -d 0.45 -q out.fp $(basename $x .gz) > similarities.new__$(basename $x .fp.gz).txt && tail similarities.new__$(basename $x .fp.gz).txt && kripodb similarities import --nrrows $nrrows similarities.new__$(basename $x .fp.gz).txt fragments.sqlite similarities.new__$(basename $x .fp.gz).h5 && rm similarities.new__$(basename $x .fp.gz).txt $(basename $x .gz)" &
done
wait
