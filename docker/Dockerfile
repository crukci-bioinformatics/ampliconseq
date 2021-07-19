FROM maven:3.6.3-openjdk-11

LABEL authors="Matt Eldridge" \
      version="1.0-snapshot" \
      description="Variant calling pipeline for amplicon sequencing data"

RUN apt-get update && \
    apt-get install -y git curl procps && \
    apt-get clean

ARG CONDA_VERSION=py38_4.9.2

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -o miniconda3.sh && \
    mkdir -p /opt && \
    sh miniconda3.sh -b -p /opt/conda && \
    rm miniconda3.sh

COPY conda.yml .

RUN /opt/conda/bin/conda env create -f conda.yml && /opt/conda/bin/conda clean -a

ENV PATH /opt/conda/envs/ampliconseq/bin:$PATH

# ensures that a new build is run if commits have been made since the last build
ADD "https://api.github.com/repos/crukci-bioinformatics/htsjdk-tools/commits?per_page=1" /opt/htsjdk-tools/.latest_commit

ARG HTSJDK_TOOLS_VERSION=master

RUN git clone --branch ${HTSJDK_TOOLS_VERSION} --depth 1 https://github.com/crukci-bioinformatics/htsjdk-tools /opt/htsjdk-tools/build && \
    cd /opt/htsjdk-tools/build && \
    mvn install

# ensures that a new build is run if commits have been made since the last build
ADD "https://api.github.com/repos/crukci-bioinformatics/ampliconseq/commits?per_page=1" /opt/ampliconseq/.latest_commit

ARG AMPLICONSEQ_VERSION=master

RUN git clone --branch ${AMPLICONSEQ_VERSION} --depth 1 https://github.com/crukci-bioinformatics/ampliconseq /opt/ampliconseq/build && \
    cd /opt/ampliconseq/build && \
    mvn package && \
    cd .. && \
    tar zxf build/target/ampliconseq-*-tools.tar.gz && \
    ln -s ampliconseq-* current

ENV PATH /opt/ampliconseq/current/bin:${PATH}

