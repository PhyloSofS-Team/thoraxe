FROM ubuntu:bionic

# to avoid interaction with tzdata during the installation
ARG DEBIAN_FRONTEND=noninteractive

# to avoid the UnicodeEncodeError (default LANG is C)
ENV LANG C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    python3 \
    python3-pip \
    python3-wheel \
    python3-setuptools \
    python3-dev \
    build-essential

WORKDIR /app

COPY docker_banner.sh docker_banner.sh

RUN cat /app/docker_banner.sh >> ~/.bashrc

RUN git clone https://github.com/PhyloSofS-Team/ProGraphMSA.git prographmsa && \
    chmod a+x ./prographmsa/bin/ProGraphMSA_64 && \
    cp ./prographmsa/bin/ProGraphMSA_64 /bin/ProGraphMSA && \
    rm -fr prographmsa

RUN git clone https://github.com/PhyloSofS-Team/thoraxe.git && \
    python3 -m pip install ./thoraxe

WORKDIR /project

CMD ["/bin/bash"]
