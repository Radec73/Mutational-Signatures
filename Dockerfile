FROM ubuntu:22.04

MAINTAINER Radovan Cyprich

# Set environment variables
ENV R_VERSION=4.3.2 \
    PYTHON_VERSION=3.9.7

WORKDIR /dp

# Install system dependencies
RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
      apt-transport-https \
      python3-pip \
      python3 \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      cmake \
      build-essential \
      curl \
      unzip \
      wget \
      r-base \
      sudo

RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('remotes', 'data.table', 'dplyr', 'purrr', 'tidyr', 'furrr', 'Rcpp', 'cowplot', 'NMF', 'ggpubr', 'cli', 'reticulate', 'roxygen2'))" && \
    R -e "BiocManager::install('BSgenome')" && \
    R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" && \
    R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')" && \
    R -e "BiocManager::install('ShixiangWang/sigminer@v2.3.0', dependencies = TRUE)" && \
    R -e "library('sigminer'); load(system.file('extdata', 'toy_copynumber_tally_W.RData', package = 'sigminer', mustWork = TRUE)); mat = cn_tally_W[['nmf_matrix']]; print(mat);"

# Install R packages
#COPY installed_packages.txt .
#RUN R -e 'install.packages(as.character(readLines("installed_packages.txt")), repos="https://cloud.r-project.org/")'

# Install Python dependencies
#COPY requirements.txt .
#RUN pip install --upgrade pip --upgrade setuptools --no-cache-dir -r requirements.txt

RUN pip install signatureanalyzer &&\
    pip install SigProfilerAssignment &&\
    pip install SigProfilerExtractor && \
    pip install --upgrade setuptools

COPY . .

# Set the default command to run when container starts
CMD ["python", "superscript.py"]


