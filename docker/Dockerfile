FROM maven:3.9.11

LABEL authors="Matt Eldridge" \
      version="1.1-snapshot" \
      description="Variant calling pipeline for amplicon sequencing data"

RUN apt-get update && \
    apt-get install -y git curl procps && \
    apt-get clean

ARG MINIFORGE3_VERSION=25.3.1-0

RUN curl -LO https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE3_VERSION}/Miniforge3-${MINIFORGE3_VERSION}-Linux-x86_64.sh && \
    curl -LO https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE3_VERSION}/Miniforge3-${MINIFORGE3_VERSION}-Linux-x86_64.sh.sha256 && \
    sha256sum -c Miniforge3-${MINIFORGE3_VERSION}-Linux-x86_64.sh.sha256 && \
    rm Miniforge3-${MINIFORGE3_VERSION}-Linux-x86_64.sh.sha256 && \
    mkdir -p /opt && \
    sh Miniforge3-${MINIFORGE3_VERSION}-Linux-x86_64.sh -b -p /opt/miniforge3 && \
    rm Miniforge3-${MINIFORGE3_VERSION}-Linux-x86_64.sh

COPY conda.yml .

RUN /opt/miniforge3/bin/mamba env create -f conda.yml && \
    /opt/miniforge3/bin/mamba clean -a

ENV PATH /opt/miniforge3/envs/ampliconseq/bin:$PATH

# ensures that a new build is run if commits have been made since the last build
ADD "https://api.github.com/repos/crukci-bioinformatics/htsjdk-tools/commits?per_page=1" /opt/htsjdk-tools/.latest_commit

ARG HTSJDK_TOOLS_VERSION=1.1

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

