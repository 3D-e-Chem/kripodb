Data update
===========

The Kripo data can be updated in 2 ways:

.. toctree::
    :maxdepth: 1

    Baseline update <baseline-update.rst>
    Incremental update <incremental-update.rst>

Steps
-----

Overview of steps involved in updating Kripo:

1. Create staging directory
2. Create sub-pocket pharmacophore fingerprints
3. Create fragment information
4. Add new fragment information to fragment sqlite db
5. Populate PDB metadata in fragments database
6. Check no fragments are duplicated
7. Calculate similarity scores between fingerprints
8. Convert pairs file into dense similarity matrix
9. Switch staging to current

.. note:: Steps 2 through 3 require undisclosed scripts
.. note:: Steps 4 and 6 through 7 can be done using the KripoDB Python library.

.. todo:: Remove Kripo fragment/fingerprints of obsolete PDBs (ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat)

Disk layout
-----------

Directories for Kripo:

* **current/**, directory which holds current dataset
* **staging/**, which is used to compute new items and combine new and old items.
* **old/**, which is used as a backup containing the previous update.

Files and directories for a data set (inside current, staging and old directories):

* **FRAGMENT_PPHORES/**, directory with pharmacophores
* **out.fp.sqlite**, fingerprints file
* **fragments.sqlite**, fragment information database file
* **similarities.h5**, similarities as pairs table
* **similarities.packedfrozen.h5**, similarities as dense matrix

Input directories:

* **$PDBS_ADDED_DIR**, directory containing new PDB files to be processed

Requirements
------------

* Slurm batch scheduler
* KripoDB and it's dependencies installed and in path
* Posix filesystem, NFS of Virtualbox share dp not accept writing of hdf5 or sqlite files
