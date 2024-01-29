# syntax=docker/dockerfile:1.2
FROM nvcr.io/nvidia/nvhpc:21.3-devel-cuda11.2-ubuntu20.04
ENV SPACK_PREFIX=/root/c2sm-spack
ENV SPACK_SYSTEM_CONFIG_PATH=/root/c2sm-spack/sysconfigs/unknown
ENV SPACK_ROOT=/root/c2sm-spack/spack

RUN apt-get -yqq update \
 && apt-get -yqq install --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        gnupg2 \
        iproute2 \
        locales \
        python3 \
        python3-pip \
        python3-setuptools \
        tcl \
        unzip \
 && locale-gen en_US.UTF-8 \
 && pip3 install boto3 \
 && rm -rf /var/lib/apt/lists/*

# setup spack-c2sm instance
RUN git clone --depth 1 --recurse-submodules --shallow-submodules -b dev_v0.18.1 https://github.com/C2SM/spack-c2sm.git ${SPACK_PREFIX} && \
    cd ${SPACK_PREFIX} && \
    . ./setup-env.sh

# add additional packages from Sergey
COPY packages/icon ${SPACK_PREFIX}/packages/icon
COPY packages/rttov ${SPACK_PREFIX}/packages/rttov
COPY packages/libcdi-pio ${SPACK_PREFIX}/packages/libcdi-pio
COPY packages/yaxt ${SPACK_PREFIX}/packages/yaxt
COPY packages/yac ${SPACK_PREFIX}/packages/yac
COPY packages/scales-ppm ${SPACK_PREFIX}/packages/scales-ppm

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - eccodes@2.19.1%nvhpc +aec~fortran~ipo~memfs~netcdf~openmp~png~pthreads~python+shared~tools build_type=RelWithDebInfo definitions=auto jp2k=openjpeg samples=auto" \
&&   echo "  - mpich@3.4.3%nvhpc ~argobots~cuda+fortran+hwloc+hydra+libxml2+pci~rocm+romio~slurm~two_level_namespace~verbs+wrapperrpath datatype-engine=auto device=ch4 netmod=ofi pmi=pmi" \
&&   echo "  - netcdf-c@4.8.1%nvhpc +mpi~parallel-netcdf+pic+shared" \
&&   echo "  - netcdf-fortran@4.5.4%nvhpc ~doc+pic+shared" \
&&   echo "  - netlib-lapack@3.10.1%nvhpc ~external-blas~ipo+lapacke+shared~xblas build_type=RelWithDebInfo" \
&&   echo "  - hdf5@1.12.2%nvhpc ~cxx~fortran+hl~ipo~java+mpi+shared+szip~threadsafe+tools api=default build_type=RelWithDebInfo" \
&&   echo "  - findutils%gcc" \
&&   echo "  concretizer:" \
&&   echo "    unify: true" \
&&   echo "  config:" \
&&   echo "    install_tree:" \
&&   echo "      root: /opt/software" \
&&   echo "  view: /opt/view" \
&&   echo "  packages:" \
&&   echo "    cuda:"\
&&   echo "      externals:" \
&&   echo "      - spec: cuda@11.2.142%gcc" \
&&   echo "        prefix: /opt/nvidia/hpc_sdk/Linux_x86_64/21.3/cuda/11.2/" \
&&   echo "      buildable: false") > /opt/spack-environment/spack.yaml

RUN mkdir -p $SPACK_ROOT

RUN ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
    /usr/local/bin/docker-shell \
    && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
    /usr/local/bin/interactive-shell \
    && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
    /usr/local/bin/spack-env

SHELL ["docker-shell"]

RUN cd /opt/spack-environment && \
    spack compiler find && \
    spack external find --not-buildable cmake perl python && \
    spack env activate . && \
    spack install --fail-fast

ENTRYPOINT ["/bin/bash", "/usr/local/bin/spack-env"]
CMD ["interactive-shell"]
