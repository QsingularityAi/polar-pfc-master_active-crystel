FROM debian:buster

LABEL MAINTAINER="simon.praetorius@tu-dresden.de"
LABEL VERSION="v0.1"

RUN apt-get update -y \
 && apt-get install -y --no-install-recommends \
      build-essential \
      ca-certificates \
      cmake \
      g++-8 \
      gcc-8 \
      git \
      libalglib-dev \
      libboost-dev \
      libeigen3-dev \
      libfftw3-dev \
      zlib1g-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# download external dependencies
RUN mkdir -p /tmp/sources && cd /tmp/sources \
 && git clone https://github.com/jlblancoc/nanoflann.git --single-branch --branch master \
 && git clone https://bitbucket.org/nschaeff/shtns.git --single-branch --branch master

# install shtns
RUN cd /tmp/sources/shtns \
 && ./configure --prefix=/opt/software/shtns --enable-openmp \
 && make \
 && make install

# install nanoflann
RUN mkdir -p /opt/software/nanoflann/include \
 && cp /tmp/sources/nanoflann/include/*.hpp /opt/software/nanoflann/include

# clean source files
RUN rm -rf /tmp/sources

# compile the project code
COPY ./code /app
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Release \
          -DSHTNS_LIB:PATH=/opt/software/shtns/lib \
          -DSHTNS_INC:PATH=/opt/software/shtns/include \
          -DNANOFLANN_INC:PATH=/opt/software/nanoflann/include \
          -DALGLIB_INC:PATH=/usr/include/libalglib /app \
 && make

ENV PATH=/app/build:$PATH
ENTRYPOINT ["/app/build/polar_pfc"]


