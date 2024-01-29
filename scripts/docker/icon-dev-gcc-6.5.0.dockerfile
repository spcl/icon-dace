FROM ubuntu:18.04 as prepare-spack

ARG SSH_PRIVATE_KEY
ARG SPACK_ROOT=/opt/spack

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -yqq update \
 && apt-get -yqq install \
      build-essential \
      ca-certificates \
      curl \
      file \
      gfortran \
      git \
      python3 \
      unzip

RUN mkdir /root/.ssh \
 && echo "${SSH_PRIVATE_KEY}" > /root/.ssh/id_rsa \
 && chmod 600 /root/.ssh/id_rsa \
 && ssh-keyscan gitlab.dkrz.de > /root/.ssh/known_hosts

RUN git clone --depth 1 git@gitlab.dkrz.de:icon-spack/spack.git ${SPACK_ROOT}

RUN (echo "#! /usr/bin/env bash" \
 &&  echo ". $SPACK_ROOT/share/spack/setup-env.sh" \
 &&  echo "exec bash -c \"\$*\"") > /usr/local/bin/docker-shell \
 && chmod +x /usr/local/bin/docker-shell

SHELL ["docker-shell"]

# Create the package cache
RUN spack spec zlib

FROM prepare-spack as builder

RUN apt-get -yqq update \
 && apt-get -yqq install \
      g++-6 \
      gcc-6 \
      gfortran-6 \
 && spack compiler find

RUN mkdir /opt/spack-environment \
 && cd /opt/spack-environment \
 && (echo "spack:" \
 &&  echo "  specs:" \
 &&  echo "    - autoconf" \
 &&  echo "    - automake" \
 &&  echo "    - libtool" \
 &&  echo "    - netcdf-fortran" \
 &&  echo "    - eccodes" \
 &&  echo "    - serialbox" \
 &&  echo "    - netlib-lapack" \
 &&  echo "    - cdo" \
 &&  echo "  packages:" \
 &&  echo "    all:" \
 &&  echo "      providers:" \
 &&  echo "        mpi: [mpich]" \
 &&  echo "      compiler: [gcc@6.5.0]" \
 &&  echo "    hdf5:" \
 &&  echo "      variants: +szip+threadsafe+mpi" \
 &&  echo "    netcdf-c:" \
 &&  echo "      variants: +mpi+parallel-netcdf" \
 &&  echo "  concretization: together" \
 &&  echo "  config:" \
 &&  echo "    install_tree: /opt/software" \
 &&  echo "  view: /opt/view") > spack.yaml \
 && spack env activate . \
 && spack install \
 && spack gc -y

FROM ubuntu:18.04
MAINTAINER Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>

RUN apt-get -yqq update \
 && apt-get -yqq install \
      --no-install-recommends \
      file \
      g++-6 \
      gcc-6 \
      gfortran-6 \
      git \
      ksh \
      less \
      make \
      nano \
      python \
      rsync \
      vim \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/view /opt/view

RUN useradd -ms /bin/bash icon
USER icon
WORKDIR /home/icon

RUN (echo 'export PATH="/opt/view/bin:$PATH"' \
 &&  echo 'export ICON_SW_PREFIX="/opt/view"' \
 &&  echo 'export PYTHONPATH="/opt/view/python/pp_ser:$PYTHONPATH"' \
 &&  echo 'export LANG=C.UTF-8' \
 &&  echo 'export PS1="\[$(tput bold)\]\[$(tput setaf 2)\]\u@\[$(tput setaf 3)\]dev-gcc-6.5.0\[$(tput sgr0)\]:\w$ "') > ~/.bashrc

