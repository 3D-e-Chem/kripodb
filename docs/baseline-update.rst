Baseline update
===============

.. contents::

The Kripo data set is generated from scratch every year or when algorithms change.

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

3. Pharmacophores
^^^^^^^^^^^^^^^^^

The raw pharmacophores are stored in the FRAGMENT_PPHORES sub-directory.
Each pocket has a \*_pphore.sd.gz file which contains the pharmacophore points of the whole pocket and
a \*_pphores.txt file which contains the indexes of pharmacophore points for each sub pocket or fragment.
The raw pharmacophores need to be added to the pharmacophores datafile with::

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

    jid_dups=$(sbatch --parsable -n 1 -J check_dups --dependency=afterok:$jid_db $SCRIPTS/baseline_duplicates.sh)

7. Calculate similarity scores between fingerprints
---------------------------------------------------

The similarities between fingerprints can be calculated with::

    all_chunks=$(ls *fp.gz |wc -l)
    jid_fpunzip=$(sbatch --parsable -n $all_chunks -J fpunzip --dependency=afterok:$jid_dups $SCRIPTS/baseline_fpunzip.sh)
    nr_chunks=($(ls *.fp.gz|wc -l) * $(ls *.fp.gz|wc -l)/2 - $(ls *.fp.gz|wc -l))
    jid_fpneigh=$(sbatch --parsable -n $nr_chunks -J fpneigh --dependency=afterok:$jid_fpunzip $SCRIPTS/baseline_similarities.sh)
    jid_fpzip=$(sbatch --parsable -n $all_chunks -J fpzip --dependency=afterok:$jid_fpneigh $SCRIPTS/baseline_fpzip.sh)
    jid_merge_matrices=$(sbatch --parsable -n 1 -J merge_matrices --dependency=afterok:$jid_fpneigh $SCRIPTS/baseline_merge_similarities.sh)

To prevent duplicates similarities of a chunk against itself should ignore the upper triangle.

.. todo:: Don't fpneigh run sequentially but submit to batch queue system and run in parallel

8. Convert pairs file into dense similarity matrix
--------------------------------------------------

.. tip:: Converting the pairs file into a dense matrix goes quicker with more memory.

The following commands converts the pairs into a compressed dense matrix::

    jid_compress_matrix=$(sbatch --parsable -n 1 -J compress_matrix --dependency=afterok:$jid_merge_matrices $SCRIPTS/freeze_similarities.sh)

The output of this step is ready to be served as a webservice using the `kripodb serve` command.

9. Switch staging to current
----------------------------

The webserver and webservice are configure to look in the `current` directory for files.

The staging can be made current with the following commands::

    mv current old
    mv staging current
