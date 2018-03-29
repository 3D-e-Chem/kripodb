# Using anaconda for RDKit,
# because OS version is too old,
# compiling it takes too long and
# there is no docker image with an up to date RDKit.
FROM continuumio/miniconda3

MAINTAINER Stefan Verhoeven <s.verhoeven@esciencecenter.nl>

RUN conda update conda -y -q && \
conda install -y -q -c conda-forge -c rdkit rdkit pytables coverage pandas numpy scipy gxx_linux-64

ENV PATH /opt/conda/bin:$PATH

ADD . /app

WORKDIR /app

RUN /opt/conda/bin/python setup.py install

CMD ["kripodb", "--help"]
