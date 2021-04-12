FROM maven:3.6.3-openjdk-11

MAINTAINER Matt Eldridge "matthew.eldridge@cruk.cam.ac.uk" \
           version "1.0.0" \
           description "Variant calling pipeline for amplicon sequencing data"

RUN apt update && apt install -y git curl procps

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /opt/miniconda3.sh && \
    sh /opt/miniconda3.sh -b -p /opt/conda && \
    rm /opt/miniconda3.sh

ENV PATH ${PATH}:/opt/conda/bin

COPY conda.yml .

RUN conda env create -f conda.yml && conda clean -a

ENV PATH /opt/conda/envs/ampliconseq/bin:$PATH

ADD "https://api.github.com/repos/crukci-bioinformatics/ampliconseq/commits?per_page=1" /opt/ampliconseq/.latest_commit

RUN git clone https://github.com/crukci-bioinformatics/ampliconseq /opt/ampliconseq/build && \
    cd /opt/ampliconseq/build && \
    mvn package && \
    cd .. && \
    tar zxf build/target/ampliconseq-*-tools.tar.gz && \
    ln -s ampliconseq-* current

ENV PATH /opt/ampliconseq/current/bin:${PATH}

