# DO NOT CHANGE
FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

WORKDIR /tmp/docker-build/work/

SHELL [ \
    "/usr/bin/env", "bash", \
    "-o", "errexit", \
    "-o", "pipefail", \
    "-o", "nounset", \
    "-o", "verbose", \
    "-o", "errtrace", \
    "-O", "inherit_errexit", \
    "-O", "shift_verbose", \
    "-c" \
]
ENV TZ='Etc/UTC'
ENV LANG='en_US.UTF-8'

ARG DEBIAN_FRONTEND=noninteractive

# Get secure apt key for CRAN
RUN apt-get update --yes && \
    apt-get install --yes dirmngr apt-transport-https gnupg2 && \
    gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' && \
    gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' | tee /etc/apt/trusted.gpg.d/cran_debian_key.asc

# install R requirements
RUN echo "deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/" | tee /etc/apt/sources.list.d/cran.list && \
    apt-get update --yes && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes r-base r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev wget

COPY latch_environment.R /tmp/docker-build/work/latch_environment.R
RUN Rscript /tmp/docker-build/work/latch_environment.R

COPY clusterProfiler_4.12.2.tar.gz /tmp/docker-build/work/clusterProfiler_4.12.2.tar.gz
RUN R -e "install.packages('/tmp/docker-build/work/clusterProfiler_4.12.2.tar.gz', repos = NULL, type='source')"

COPY KEGG.db_2.8.0.tar.gz /tmp/docker-build/work/KEGG.db_2.8.0.tar.gz
RUN R -e "install.packages('/tmp/docker-build/work/KEGG.db_2.8.0.tar.gz', repos = NULL, type='source')"

# Latch SDK
# DO NOT REMOVE
RUN pip install latch==2.52.2
RUN mkdir /opt/latch

# Copy workflow data (use .dockerignore to skip files)
COPY . /root/

# Latch workflow registration metadata
# DO NOT CHANGE
ARG tag

# DO NOT CHANGE
ENV FLYTE_INTERNAL_IMAGE $tag

WORKDIR /root
