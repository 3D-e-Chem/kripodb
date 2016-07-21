Update pipelines
================

The Kripo data can be updated using pipelines.
A pipeline is a sequence of steps, where each steps executes CLI commands.

.. toctree::
    :maxdepth: 1

    Baseline update <baseline-update.rst>
    Weekly update <weekly-update.rst>

Steps
-----

Overview of steps involved in updating Kripo:

1. Create staging directory
2. Create sub-pocket pharmacophore fingerprints
3. Create fragment information
4. Add new fragment information to fragment sqlite db
5. Calculate similarity scores between fingerprints
6. Add new similarity scores to similarity pairs file
7. Convert pairs file into dense similarity matrix
8. Switch staging to current

.. note:: Steps 2 through 3 require undisclosed scripts
.. note:: Steps 4 and 6 through 7 can be done using the KripoDB Python library.

Disk layout
-----------

Directories for Kripo:

* **current/**, directory which holds current pharmacophores, fingerprints, fragment information and similarity matrix.
* **staging/**, which is used to compute new items and combine new and old items.
* **old/**, which is used as a backup containing the previous update.

Input directories:

* **$PDBS_ADDED_DIR**, directory containing new PDB files to be processed

.. note: TODO remove Kripo fragment/fingerprints of obsolete PDBs
