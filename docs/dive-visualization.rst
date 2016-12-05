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

Output datasets

1. 2D
2. 3D

1. LargeVis input file from Kripo similarity matrix
---------------------------------------------------

Dump the similarity matrix to csv::

    kripodb similarities export --no_header similarities.h5 similarities.txt
    perl -ne 'print if /.*frag1\s.*frag1/' < similarities.txt > similarities.frag1.txt

2. Perform embedding using LargeVis
-----------------------------------

Get or compile LargeVis binaries from https://github.com/lferry007/LargeVis

Compile using miniconda::

    conda install gsl gcc
    cd LargeVis/Linux
    c++ LargeVis.cpp main.cpp -o LargeVis -lm -pthread -lgsl -lgslcblas -Ofast -march=native -ffast-math
    LD_LIBRARY_PATH=$CONDA_PREFIX/lib
    cp LargeVis $CONDA_PREFIX/bin/


Then embed frag1 similarity matrix in 3D with::

    LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.frag1.txt -output largevis.frag1.3d.txt

Then embed frag1 similarity matrix in 2D with::

    LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.frag1.txt -output largevis.frag1.2d.txt

Then embed similarity matrix in 3D with::

    LargeVis -fea 0 -outdim 3 -threads $(nproc) -input similarities.txt -output largevis.3d.txt

Then embed similarity matrix in 2D with::

    LargeVis -fea 0 -outdim 2 -threads $(nproc) -input similarities.txt -output largevis.2d.txt

3. Generate DiVE metadata datafiles
-----------------------------------

Command to generate properties files::

    wget -O uniprot.txt 'http://www.uniprot.org/uniprot/?query=database:pdb&format=tab&columns=id,genes(PREFERRED),families,database(PDB)'
    kripodb dive export fragments.sqlite uniprot.txt

Will generate in current working directory the following files:

* kripo.meta.txt
* kripo.props.txt
* kripo.propnames.txt

4. Create DiVE input file
-------------------------

DiVE has a script which can combine the LargeVis coordinates together with metadata. 
Download the MakeVizDataWithProperMetadata.py script from https://github.com/NLeSC/DiVE/blob/master/scripts_prepareData/MakeVizDataWithProperMetadata.py

For more information about the script see https://github.com/NLeSC/DiVE#from-output-of-largevis-to-input-of-dive .

Commands to generate new DiVE input file::

    python MakeVizDataWithProperMetadata.py -coord largevis.frag1.2d.txt -metadata kripo.meta.txt -np -kripo.props.txt -pif kripo.propnames.txt -dir frag1.3d
    python MakeVizDataWithProperMetadata.py -coord largevis.frag1.3d.txt -metadata kripo.meta.txt -np -kripo.props.txt -pif kripo.propnames.txt -dir frag1.2d
    python MakeVizDataWithProperMetadata.py -coord largevis.2d.txt -metadata kripo.meta.txt -np -kripo.props.txt -pif kripo.propnames.txt -dir 2d
    python MakeVizDataWithProperMetadata.py -coord largevis.3d.txt -metadata kripo.meta.txt -np -kripo.props.txt -pif kripo.propnames.txt -dir 3d

The generated file can be uploaded at https://nlesc.github.io/DiVE/ to visualize.