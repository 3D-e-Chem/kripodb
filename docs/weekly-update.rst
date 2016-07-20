Weekly update
=============

The Kripo data set is updated weekly with new PDB entries.

1. Create staging directory
---------------------------

Create a new directory::

  mkdir staging
  cd ..

2. Create sub-pocket pharmacophore fingerprints
-----------------------------------------------

Use directory listing of new pdb files as input::

  ls $NEW_PDBS_DIR | pdblist2fps_final_local.py

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

4. Add new fragment information to fragment sqlite db
-----------------------------------------------------

The following commands add the fragment shelve and sdf to the fragments database::

    cp ../current/fragments.sqlite .
    kripodb fragments shelve fragments.shelve fragments.sqlite
    kripodb fragments sdf fragments.sd fragments.sqlite

5. Calculate similarity scores between fingerprints
---------------------------------------------------

The similarities between the new and existing fingerprints and between new fingerprints themselves can be calculated with::

    kripodb fingerprints import out.fp out.fp.db
    kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite ../current/out.fp.sqlite out.fp.sqlite similarities.new_existing.h5
    kripodb fingerprints similarities --fragmentsdbfn fragments.sqlite out.fp.sqlite out.fp.sqlite similarities.new_new.h5

6. Add new similarity scores to similarity pairs file
-----------------------------------------------------

The following command merges the current pairs file with the new pairs files::

    kripodb similarities merge ../staging/similarities.h5 similarities.new_existing.h5 similarities.new_new.h5 similarities.h5

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

The staging can be made current with the following commands::

    mv current old
    mv staging current
