=============
Weekly update
=============

The Kripo data set is updated weekly. Below is a description of the layout and steps involved.

Layout
======

Directory called `current` which holds current pharmacophores, fingerprints, fragment information and similarity matrix.
Directory called `staging` which is used to compute new items and combine new and old items.

Steps
=====

:ref:`create-staging-directory`

:ref:`create-sub-pocket-pharmacophore-fingerprints`

:ref:`create-fragment-information`

:ref:`calculate-similarity-scores-between-fingerprints`

:ref:`Add new fragment information to fragment sqlite db`

:ref:`Add new similarity scores to similarity pairs file`

:ref:`Convert pairs file into dense similarity matrix`

:ref:`Switch staging to current`

.. note:: Steps 4..7 can be done using the KripDB Python library.

.. _create-staging-directory:

1. Create staging directory
---------------------------

Create a new directory::

  mkdir staging
  cd ..

.. _create-sub-pocket-pharmacophore-fingerprints:

2. Create sub-pocket pharmacophore fingerprints
-----------------------------------------------

Use directory listing of new pdb files as input::

  ls $NEW_PDBS_DIR | pdblist2fps_final_local.py

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

  fragid2sd.py

.. _calculate-similarity-scores-between-fingerprints:

4. Calculate similarity scores between fingerprints
---------------------------------------------------

The similarities between the new and existing fingerprints and between new fingerprints themselves can be calculated with::

    kripodb fingerprints import out.fp out.fp.db
    kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite ../current/out.fp.sqlite out.fp.sqlite similarities.new_existing.h5
    kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite out.fp.sqlite out.fp.sqlite similarities.new_new.h5

.. _Add new fragment information to fragment sqlite db:

5. Add new fragment information to fragment sqlite db
-----------------------------------------------------

The following commands add the fragment shelve and sdf to the fragments database::

    cp ../current/fragments.sqlite .
    kripodb fragments shelve fragments.shelve fragments.sqlite
    kripodb fragments sdf fragments.sd fragments.sqlite

.. _Add new similarity scores to similarity pairs file:

6. Add new similarity scores to similarity pairs file
-----------------------------------------------------

The following command merges the current pairs file with the new pairs files::

    kripodb similarities merge ../staging/similarities.h5 similarities.new_existing.h5 similarities.new_new.h5 similarities.h5

.. _Convert pairs file into dense similarity matrix:

7. Convert pairs file into dense similarity matrix
--------------------------------------------------

.. note:: Converting the pairs file into a dense matrix goes quicker with more memory.

The following commands converts the pairs into a compressed dense matrix::

    kripodb similarities freeze similarities.h5 similarities.frozen.h5
    ptrepack --complevel 6 --complib blosc:zlib similarities.frozen.h5 similarities.packedfrozen.h5
    rm similarities.h5 similarities.frozen.h5

The output of this step is ready to be served as a webservice using the `kripodb serve` command.

.. _Switch staging to current:

8. Switch staging to current
----------------------------

The webserver and webservice are configure to look in the `current` directory for files.

The staging can be made current with the following commands::

    mv current old
    mv staging current
