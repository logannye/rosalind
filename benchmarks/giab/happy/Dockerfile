FROM ubuntu:18.04@sha256:dca176c9663a7ba4c1f0e710986f5a25e672842963d95b960191e2d9f7185ebe

ARG DEBIAN_FRONTEND=noninteractive
ARG HAPPY_URL=https://github.com/Illumina/hap.py/archive/refs/tags/v0.3.15.tar.gz
ARG HAPPY_SHA256=838912f1ab805224110af3c1035849d5483b5ce4c565ae0f0db7f7284af06676
ARG RTG_URL=https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-nojre.zip
ARG RTG_SHA256=5efc8cb1a44019b97dce8d324e3a0b3dc3328d022fc9db13c80b09996757c4aa
COPY requirements-py27.txt /opt/requirements-py27.txt

RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf build-essential bzip2 ca-certificates cmake curl cython git \
    libboost-all-dev libbz2-dev libncurses5-dev openjdk-8-jre-headless \
    pkg-config python python2.7 python2.7-dev python-distribute python-numpy \
    python-pandas python-pip python-psutil python-pysam python-scipy \
    python-setuptools python-nose unzip zlib1g-dev && \
    pip install --disable-pip-version-check --no-deps --require-hashes -r /opt/requirements-py27.txt && \
    python -c 'import bx, numpy, pandas, pysam, scipy, six; assert six.__version__ == "1.16.0"' && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build
RUN curl --fail --location --retry 3 "$HAPPY_URL" -o happy.tar.gz && \
    echo "$HAPPY_SHA256  happy.tar.gz" | sha256sum --check --strict && \
    tar xzf happy.tar.gz --strip-components=1 && \
    python install.py /opt/hap.py --no-tests && \
    curl --fail --location --retry 3 "$RTG_URL" -o rtg.zip && \
    echo "$RTG_SHA256  rtg.zip" | sha256sum --check --strict && \
    unzip -q rtg.zip -d /opt && \
    mv /opt/rtg-tools-3.12.1 /opt/rtg-tools && \
    printf '%s\n' 'RTG_MEM=2G' 'RTG_JAVA_OPTS="-XX:-UseContainerSupport"' \
      'RTG_TALKBACK=false' 'RTG_USAGE=false' > /opt/rtg-tools/rtg.cfg && \
    /opt/hap.py/bin/hap.py --version && \
    /opt/hap.py/bin/test_haplotypes && \
    /opt/rtg-tools/rtg version && \
    rm -rf /build

# Ubuntu 18.04's Java 8 cgroup-v2 detection fails on current runners. RTG uses
# an explicit 2 GiB Java heap (overridable with RTG_MEM) rather than host-RAM
# percentage discovery; the heap setting is not a whole-process memory cap.
WORKDIR /data
ENTRYPOINT ["/opt/hap.py/bin/hap.py"]
