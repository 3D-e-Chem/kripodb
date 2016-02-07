# Using anaconda for RDKit,
# because OS version is too old,
# compiling it takes too long and
# there is no docker image with an up to date RDKit.
FROM continuumio/miniconda

MAINTAINER Stefan Verhoeven <s.verhoeven@esciencecenter.nl>

RUN /opt/conda/bin/conda install -y -q -c https://conda.anaconda.org/rdkit rdkit pytables bsddb coverage pandas && \
/opt/conda/bin/conda clean -y -s -p -t -l -i

ENV PATH /opt/conda/bin:$PATH

# To install the intbitset python package gcc is required
# Anaconda install of gcc failed,
# falling back to OS supplied gcc
RUN apt-get update && apt-get install  -y -q gcc && \
pip install intbitset && \
apt-get remove -y gcc && apt-get autoremove -y && \
rm -rf /var/lib/apt/lists/*

ADD . /app

WORKDIR /app

RUN /opt/conda/bin/python setup.py install

CMD ["kripodb", "--help"]

