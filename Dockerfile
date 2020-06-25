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
    build-essential \
    # Utilities:
    tree \
    ipython3

WORKDIR /app

COPY docker_banner.sh docker_banner.sh

RUN cat /app/docker_banner.sh >> ~/.bashrc

# install numpy before socket-bio (ThorAxe's dependency)
RUN python3 -m pip install numpy && \
    python3 -m pip install thoraxe

WORKDIR /project

CMD ["/bin/bash"]
