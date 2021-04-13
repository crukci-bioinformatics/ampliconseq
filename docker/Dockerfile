FROM maven:3.6.3-openjdk-11

LABEL authors="Matt Eldridge" \
      version="1.0.0" \
      description="Variant calling pipeline for amplicon sequencing data"

RUN apt-get update && \
    apt-get install -y git curl procps && \
    apt-get clean

ARG CONDA_VERSION=py38_4.9.2
ARG CONDA_MD5=122c8c9beb51e124ab32a0fa6426c656

RUN curl https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -o miniconda3.sh && \
    echo "${CONDA_MD5}  miniconda3.sh" > miniconda3.md5 && \
    mkdir -p /opt && \
    sh miniconda3.sh -b -p /opt/conda && \
    rm miniconda3.sh miniconda3.md5

COPY conda.yml .

RUN /opt/conda/bin/conda env create -f conda.yml && /opt/conda/bin/conda clean -a

ENV PATH /opt/conda/envs/ampliconseq/bin:$PATH

# ensures that a new build is run if commits have been made since the last build
ADD "https://api.github.com/repos/crukci-bioinformatics/ampliconseq/commits?per_page=1" /opt/ampliconseq/.latest_commit

ARG VERSION=master

RUN git clone --branch ${VERSION} --depth 1 https://github.com/crukci-bioinformatics/ampliconseq /opt/ampliconseq/build && \
    cd /opt/ampliconseq/build && \
    mvn package && \
    cd .. && \
    tar zxf build/target/ampliconseq-*-tools.tar.gz && \
    ln -s ampliconseq-* current

ENV PATH /opt/ampliconseq/current/bin:${PATH}

