FROM ubuntu:22.04

LABEL maintainer="Afonso Bravo"
LABEL version="0.0.1"
LABEL description="2fast2q"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       python3 python3-pip python3-venv python3-distutils procps curl ca-certificates python3-tk \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir --upgrade pip \
    && python3 -m pip install --no-cache-dir fast2q 

WORKDIR /data
ENTRYPOINT ["/bin/bash"]
CMD []
