FROM continuumio/miniconda:4.7.12

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

#RUN conda install -c conda-forge -c bioconda -c defaults -c bioconda/label/broken bioconductor-singlecellexperiment=1.6.0
# RUN conda install  bioconductor-singlecellexperiment=1.6.0 bioconductor-s4vectors=0.24 r-monocle3=0.2.1 r-seurat=3.0.2 r-getopt=1.20.3 r-optparse=1.6.4 r-dplyr=0.8.0.1

RUN conda install r-monocle3=0.2.1 r-seurat=3.0.2 r-getopt=1.20.3 r-optparse=1.6.4 r-dplyr=0.8.0.1 bioconductor-summarizedexperiment=1.16.0



RUN mkdir /monocle
COPY run_monocle3.R /monocle/

