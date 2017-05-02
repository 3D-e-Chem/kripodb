Incremental update
==================

.. contents::

The Kripo data set can be incrementally updated with new PDB entries.

1. Create staging directory
---------------------------

Setup path with update scripts using::

    export SCRIPTS=$PWD/../kripodb/update_scripts

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


3. Pharmacophores
^^^^^^^^^^^^^^^^^

The raw pharmacophores are stored in the FRAGMENT_PPHORES sub-directory.
Each pocket has a \*_pphore.sd.gz file which contains the pharmacophore points of the whole pocket and
a \*_pphores.txt file which contains the indexes of pharmacophore points for each sub pocket or fragment.
The raw pharmacophores of the update can be added to the existing pharmacophores datafile with::

    cp ../current/pharmacophores.h5 .
    kripodb pharmacophores add FRAGMENT_PPHORES pharmacophores.h5

4. Add new fragment information to fragment sqlite db
-----------------------------------------------------

The following commands add the fragment shelve and sdf to the fragments database::

    cp ../current/fragments.sqlite .
    kripodb fragments shelve fragments.shelve fragments.sqlite
    kripodb fragments sdf fragments.sd fragments.sqlite

Step 4 and 5 can be submitted to scheduler with::

   jid_db=$(sbatch --parsable -n 1 -J db_append $SCRIPTS/db_append.sh)

5. Populate PDB metadata in fragments database
----------------------------------------------
The following command will updated the PDB metadata to fragments database::

    kripodb fragments pdb fragments.sqlite

6. Check no fragments are duplicated
------------------------------------

The similarity matrix can not handle duplicates. It will result in addition of scores::

    jid_dups=$(sbatch --parsable -n 1 -J check_dups --dependency=afterok:$jid_db $SCRIPTS/incremental_duplicates.sh)

7. Calculate similarity scores between fingerprints
---------------------------------------------------

The similarities between the new and existing fingerprints and between new fingerprints themselves can be calculated with::

    current_chunks=$(ls ../current/*fp.gz |wc -l)
    all_chunks=$(($current_chunks + 1))
    jid_fpneigh=$(sbatch --parsable -n $all_chunks -J fpneigh --dependency=afterok:$jid_dups $SCRIPTS/incremental_similarities.sh)
    jid_merge_matrices=$(sbatch --parsable -n 1 -J merge_matrices --dependency=afterok:$jid_fpneigh $SCRIPTS/incremental_merge_similarities.sh)

8. Convert pairs file into dense similarity matrix
--------------------------------------------------

.. note:: Converting the pairs file into a dense matrix goes quicker with more memory.

    The frame size (-f) should be as big as possible, 100000000 requires 6Gb RAM.

The following commands converts the pairs into a compressed dense matrix::

    jid_compress_matrix=$(sbatch --parsable -n 1 -J compress_matrix --dependency=afterok:$jid_merge_matrices $SCRIPTS/freeze_similarities.sh)

The output of this step is ready used to find similar fragments,
using either the webservice with the `kripodb serve` command or with the `kripodb similarities similar` command directly.

9. Switch staging to current
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

    mv current old && mv staging current

9.1 Merge fingerprint files (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To keep the number of files to a minimum it is advised to merge the fingerprint files from incremental updates of a year.

The incremental fingerprint files are named like `out.<year><week>.fp.gz`, to generate kripo_fingerprints_<year>_fp.gz run::

    sbatch --parsable -n 1 -J merge_fp $SCRIPTS/incremental_merge_fp.sh <year>

