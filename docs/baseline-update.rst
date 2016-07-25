Baseline update
===============

.. contents::

The Kripo data set is generated from scratch every year or when algorithms change.

1. Create staging directory
---------------------------

Create a new directory::

  mkdir staging
  cd ..

2. Create sub-pocket pharmacophore fingerprints
-----------------------------------------------

Use directory listing of new pdb files as input::

  ls $PDBS_ADDED_DIR | pdblist2fps_final_local.py

.. todo:: Too slow when run on single cpu.
    Chunkify input, run in parallel and merge results

.. _create-fragment-information:

3. Create fragment information
------------------------------

1. Fragment shelve
^^^^^^^^^^^^^^^^^^

Where the fragment came from is stored in a Python shelve file.
It can be generated from the pharmacophore files using::

  compiledDatabase.py

2. Fragment sdf
^^^^^^^^^^^^^^^

The data generated thus far contains the molblocks of the ligands and atom nrs of each fragment.
The fragment molblocks can be generated into a fragment sdf file with::

  fragid2sd.py > fragments.sd

4. Add new fragment information to fragment sqlite db
-----------------------------------------------------

The following commands add the fragment shelve and sdf to the fragments database::

    cp ../current/fragments.sqlite .
    kripodb fragments shelve fragments.shelve fragments.sqlite
    kripodb fragments sdf fragments.sd fragments.sqlite

5. Populate PDB metadata in fragments database
----------------------------------------------
The following command will updated the PDB metadata to fragments database::

    kripodb fragments pdb fragments.sqlite

6. Calculate similarity scores between fingerprints
---------------------------------------------------

The similarities between fingerprints can be calculated with::

    nrrows = 10000000
    for x in $(ls *.fp)
    do
    for y in $(ls *.fp)
    do
    if [ "$x" = "$y" ]
    then
    fpneigh -m Mod_Tanimoto -d 0.45 -q $x $y | kripodb similarities import --nrrows $nrrows --ignore_upper_triangle - fragments.sqlite similarities.$(basename $x .fp)_$(basename $y .fp).h5
    elif [[ $x < $y ]]
    then
    fpneigh -m Mod_Tanimoto -d 0.45 -q $x $y | kripodb similarities import --nrrows $nrrows - fragments.sqlite similarities.$(basename $x .fp)_$(basename $y .fp).h5
    fi
    done
    done

    # Compact the fingerprint file (makebits ascii format)
    for x in $(ls *.fp)
    do
    gzip $x
    done

    # Merge
    kripodb similarities merge similarities.*.h5 similarities.h5

To prevent duplicates similarities of a chunk against itself should ignore the upper triangle.

.. todo:: Don't fpneigh run sequentially but submit to batch queue system and run in parallel

7. Convert pairs file into dense similarity matrix
--------------------------------------------------

.. tip:: Converting the pairs file into a dense matrix goes quicker with more memory.

The following commands converts the pairs into a compressed dense matrix::

    kripodb similarities freeze similarities.h5 similarities.frozen.h5
    ptrepack --complevel 6 --complib blosc:zlib similarities.frozen.h5 similarities.packedfrozen.h5
    rm similarities.h5 similarities.frozen.h5

The output of this step is ready to be served as a webservice using the `kripodb serve` command.

8. Switch staging to current
----------------------------

The webserver and webservice are configure to look in the `current` directory for files.

The staging can be made current with the following commands::

    mv current old
    mv staging current
