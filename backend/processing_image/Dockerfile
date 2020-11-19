FROM rocker/r-ver:4.0.3

# Most of these are required for SCE
RUN apt-get update && apt-get -y install \
      git \
      python3-pip \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
      libcurl4-openssl-dev \
      libglpk-dev \
      libxml2-dev \
      wget \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Install R dependencies
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('SingleCellExperiment')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('LoomExperiment')"
RUN R -e "install.packages(c('devtools', 'igraph')); library(devtools); devtools::install_github('cellgeni/sceasy')"
RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install python dependencies
RUN pip3 install loompy==3.0.6 scanpy==1.6.0 python-igraph==0.8.3 louvain==0.7.0
RUN pip3 install git+https://github.com/chanzuckerberg/cellxgene.git#egg=cellxgene
RUN pip3 install awscli

ADD code /

RUN chmod a+x /dataset_processing/run.py

CMD ["/dataset_processing/run.py"]
