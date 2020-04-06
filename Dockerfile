FROM continuumio/miniconda:4.7.12

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install r-monocle3
