FROM maven:3.6.3-openjdk-11

MAINTAINER Matt Eldridge "matthew.eldridge@cruk.cam.ac.uk" \
           version "1.0.0" \
           description "Variant calling pipeline for amplicon sequencing data"

RUN apt update && apt install -y \
    git \
    curl \
    procps \
    gdebi-core \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev

ENV R_VERSION=4.0.4

RUN curl -O https://cdn.rstudio.com/r/debian-10/pkgs/r-${R_VERSION}_1_amd64.deb && \
    gdebi -n r-${R_VERSION}_1_amd64.deb && \
    rm r-${R_VERSION}_1_amd64.deb

ENV PATH ${PATH}:/opt/R/${R_VERSION}/bin

RUN R -e 'install.packages(c("optparse", "svglite", "tidyverse"), repos = "https://cloud.r-project.org")'

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /opt/miniconda3.sh && \
    sh /opt/miniconda3.sh -b -p /opt/conda && \
    rm /opt/miniconda3.sh

ENV PATH ${PATH}:/opt/conda/bin

COPY conda.yml .

RUN conda env update -n root -f conda.yml && conda clean -a

ADD "https://api.github.com/repos/crukci-bioinformatics/ampliconseq/commits?per_page=1" /opt/ampliconseq_latest_commit
RUN git clone https://github.com/crukci-bioinformatics/ampliconseq /opt/ampliconseq && \
    cd /opt/ampliconseq && \
    mvn package

ENV PATH /opt/ampliconseq/target/tools/bin:${PATH}

