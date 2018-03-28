# Using anaconda for RDKit,
# because OS version is too old,
# compiling it takes too long and
# there is no docker image with an up to date RDKit.
FROM continuumio/miniconda3

MAINTAINER Stefan Verhoeven <s.verhoeven@esciencecenter.nl>

RUN /opt/conda/bin/conda install -y -q -c conda-forge -c rdkit rdkit pytables coverage pandas numpy scipy && \
/opt/conda/bin/conda clean -y -s -p -t -l -i && apt-get update && apt-get install -y build-essential && rm -rf /var/lib/apt/lists/*

ENV PATH /opt/conda/bin:$PATH

ADD . /app

WORKDIR /app

RUN /opt/conda/bin/python setup.py install

CMD ["kripodb", "--help"]

