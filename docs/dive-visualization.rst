DiVE visualization
==================

DiVE homepage at https://github.com/NLeSC/DiVE

The Kripo similarity matrix can be embedded to 2D or 3D using largevis and then visualized using DiVE.

Steps

1. LargeVis input file from Kripo similarity matrix
2. Perform embedding using LargeVis
3. Generate DiVE metadata datafiles
4. Create DiVE input file

Input datasets

1. only fragment1 or whole unfragmented ligands
2. all fragments
3. only gpcr frag1
4. only kinase frag1
5. only gpcr and kinase frag1

Output datasets

1. 2D
2. 3D

1. LargeVis input file from Kripo similarity matrix
---------------------------------------------------

Dump the similarity matrix to csv of \*frag1 fragments::

    kripodb similarities export --no_header --frag1 similarities.h5 similarities.frag1.txt

Similarities between GPCR pdb entries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the `GPCRDB web service <http://gpcrdb.org/services/reference/#!/structure/structure_list>`_ to fetch a list of PDB codes which contain GPCR proteins::

    curl -X GET --header 'Accept: application/json' 'http://gpcrdb.org/services/structure/' | jq  -r '.[] | .pdb_code' > pdb.gpcr.txt

Dump the similarity matrix to csv::

    kripodb similarities export --no_header --frag1 --pdb pdb.gpcr.txt similarities.h5 similarities.frag1.gpcr.txt

Similarities between GPCR and Kinase pdb entries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the `KLIFS KNIME nodes <https://github.com/3D-e-Chem/knime-klifs>`_ to create a file with of PDB codes of Kinases called `pdb.kinase.txt`.

Dump the similarity matrix to csv::

    cat pdb.gpcr.txt pdb.kinase.txt > pdb.gpcr.kinase.txt
    kripodb similarities export --no_header --frag1 --pdb pdb.gpcr.kinase.txt similarities.h5 similarities.frag1.gpcr.kinase.txt

2. Perform embedding using LargeVis
-----------------------------------

Get or compile LargeVis binaries from https://github.com/lferry007/LargeVis

Compile using miniconda::

    conda install gsl gcc
    cd LargeVis/Linux
    c++ LargeVis.cpp main.cpp -o LargeVis -lm -pthread -lgsl -lgslcblas -Ofast -Wl,-rpath,$CONDA_PREFIX/lib -march=native -ffast-math
    cp LargeVis $CONDA_PREFIX/bin/

Then embed frag1 similarity matrix in 3D with::

    LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.frag1.txt -output largevis.frag1.3d.txt

Then embed frag1 similarity matrix in 2D with::

    LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.frag1.txt -output largevis.frag1.2d.txt

Then embed similarity matrix in 3D with::

    LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.txt -output largevis.3d.txt

Then embed similarity matrix in 2D with::

    LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.txt -output largevis.2d.txt


The `kripo export` in step 1 and the LargeVis command can be submitted to scheduler with::

   sbatch -n 1 $SCRIPTS/dive_frag1.sh
   sbatch -n 1 $SCRIPTS/dive_frag1_gpcr_kinase.sh

3. Generate DiVE metadata datafiles
-----------------------------------

Command to generate properties files::

    wget -O uniprot.txt 'http://www.uniprot.org/uniprot/?query=database:pdb&format=tab&columns=id,genes(PREFERRED),families,database(PDB)'
    kripodb dive export --pdbtags pdb.gpcr.txt --pdbtags pdb.kinase.txt fragments.sqlite uniprot.txt

Will generate in current working directory the following files:

* kripo.props.txt
* kripo.propnames.txt

4. Create DiVE input file
-------------------------

DiVE has a script which can combine the LargeVis coordinates together with metadata. 
Download the MakeVizDataWithProperMetadata.py script from https://github.com/NLeSC/DiVE/blob/master/scripts_prepareData/MakeVizDataWithProperMetadata.py

For more information about the script see https://github.com/NLeSC/DiVE#from-output-of-largevis-to-input-of-dive .

Example command to generate new DiVE input file::

    python MakeVizDataWithProperMetadata.py -coord largevis2.similarities.frag1.gpcr.kinase.txt -metadata kripo.props.txt -np kripo.propnames.txt -json largevis2.similarities.frag1.gpcr.kinase.json -dir .

The generated file (largevis2.similarities.frag1.gpcr.kinase.json) can be uploaded at https://nlesc.github.io/DiVE/ to visualize.