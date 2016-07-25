Incremental update
==================

.. contents::

The Kripo data set can be incrementally updated with new PDB entries.

1. Create staging directory
---------------------------

Create a new directory::

  mkdir staging
  cd ..

2. Create sub-pocket pharmacophore fingerprints
-----------------------------------------------

Use directory listing of new pdb files as input::

  ls $PDBS_ADDED_DIR | pdblist2fps_final_local.py

.. note:: The process can run out of memory, so rerun can be required

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

The similarities between the new and existing fingerprints and between new fingerprints themselves can be calculated with::

    current_chunks=$(ls ../current/*fp.gz |wc -l)
    all_chunks=$(($current_chunks + 1))
    sbatch -n $all_chunks <<
    #!/bin/sh

    # Compute similarities against itself
    nrrows=10000000
    srun -J new__new -n 1 /bin/sh -c "fpneigh -m Mod_Tanimoto=0.01 -d 0.45 -q out.fp out.fp | kripodb similarities import --nrrows $nrrows --ignore_upper_triangle - fragments.sqlite similarities.new__new.h5" &

    # Compute similarities against existing fingerprint chunks
    for x in `ls ../current/*fp.gz`
    do
    srun -J new__$x -n 1 /bin/sh -c "gunzip -c $x | fpneigh -m Mod_Tanimoto=0.01 -d 0.45 -q out.fp | kripodb similarities import --nrrows $nrrows - fragments.sqlite similarities.new__$(basename $x .fp.gz).h5" &
    done
    wait
    EOF

    # Compact the fingerprint file (makebits ascii format)
    gzip out.fp
    mv out.fp.gz out.$(date +%Y%U).fp.gz

    # Add new similarities to existing similarities file
    kripodb similarities merge ../current/similarities.h5 similarities.new_existing.h5 similarities.new_new.h5 similarities.h5

To prevent duplicates similarities of a chunk against itself should ignore the upper triangle.

.. todo:: Don't fpneigh run sequentially but submit to batch queue system and run in parallel

7. Convert pairs file into dense similarity matrix
--------------------------------------------------

.. note:: Converting the pairs file into a dense matrix goes quicker with more memory.

The following commands converts the pairs into a compressed dense matrix::

    kripodb similarities freeze similarities.h5 similarities.frozen.h5
    ptrepack --complevel 6 --complib blosc:zlib similarities.frozen.h5 similarities.packedfrozen.h5
    rm similarities.h5 similarities.frozen.h5

The output of this step is ready to be served as a webservice using the `kripodb serve` command.

8. Switch staging to current
----------------------------

The webserver and webservice are configure to look in the `current` directory for files.

The current and new pharmacophores need to be combined::

    mv staging/FRAGMENT_PPHORES staging/FRAGMENT_PPHORES.new
    rsync -a current/FRAGMENT_PPHORES staging/FRAGMENT_PPHORES
    rm -r staging/FRAGMENT_PPHORES.new

.. todo:: rsync of current/FRAGMENT_PPHORES to destination, maybe too slow due large number of files.
    Switch to move old pharmacohores and rsync new pharmacophores into it when needed.

The current and new fingerprints need to be combined::

    cp -n current/*.fp.gz staging/

The staging can be made current with the following commands::

    mv current old
    mv staging current

