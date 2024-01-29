<!--
This file is written using Markdown language, which might make it difficult to
read it in a plain text editor. Please, visit ICON project page on DKRZ GitLab
(https://gitlab.dkrz.de/icon/icon/-/tree/master/config/generic) to see this file
rendered or use a Markdown viewer of your choice
(https://www.google.com/search?q=markdown+viewer).
-->

# Introduction

This directory contains generic
[configure wrappers](../../README.md##configuration-wrappers) for ICON. You can
use them either inside the [docker containers](#docker-containers) or
[natively](#native-building) on your system.

# Docker containers

The easiest way to build and run ICON on your personal machine is to use Docker
images from the
[`iconmodel` repository on Docker Hub](https://hub.docker.com/repositories/iconmodel).

The recommended workflow is to build and run ICON inside the container and
manage and edit the source code in a separate terminal, using the tools
available on your machine. This scenario implies that the directory with ICON
source code is located on the host machine and mounted to the container. You
will also need to mount a so-called `pool` directory containing ICON input
files, e.g. grid files. The contents and the layout of the directory depends on
the experiment you want to run and are not covered in this document.

Run the container in the interactive mode as follows:
```bash
docker run -it -v /path/to/icon-src:/home/icon/icon -v /path/to/pool:/home/icon/pool iconmodel/icon-dev
```
where `/path/to/icon-src` and `/path/to/pool` are paths to ICON source and
`pool` directories on your machine, and `/home/icon/icon` and `/home/icon/pool`
are respective mount points of the directories inside the container.

As a result of the previous command you will get an interactive command prompt
of the container. You can now configure, build and run ICON using the following
commands as a reference:
```console
icon@dev-gcc$ cd ./icon
icon@dev-gcc$ ./config/generic/gcc
icon@dev-gcc$ make -j4
icon@dev-gcc$ ./make_runscripts atm_tracer_Hadley
icon@dev-gcc$ cd ./run
icon@dev-gcc$ ./exp.atm_tracer_Hadley.run
```
> **_NOTE:_** To be able to run ICON inside the container, you might need to
increase the amount of RAM available to Docker (Preferences->Resources->Memory).

# Native building

The generic wrappers in this directory are written with the assumption that the
required [software libraries](#software-libraries) are installed to the same
prefix. The prefix defaults to `/opt/local` on macOS and to `/usr` on other
platforms. The default values can be overridden by setting the environment
variable `ICON_SW_PREFIX`:
```bash
export ICON_SW_PREFIX='/path/to/icon/prerequisites'
```

## Prerequisites

This section provides a list of software required for building and running
ICON. The users can build and install (to the same prefix) the listed packages
manually or use a
[package managers](https://en.wikipedia.org/wiki/Package_manager) available for
their platform. The basic instructions on how to do it on several popular
platforms are provided in section [Tested platforms](#tested-platforms).

### Building tools

- [GNU Make](https://www.gnu.org/software/make) v3.81+
- [CMake](https://cmake.org) v3.18+
- [Python](https://www.python.org) v2.6+ or v3.5+
- [Perl](https://www.perl.org) v5.10+
- Interoperable C and Fortran compilers

### Software libraries

- [MPICH](https://www.mpich.org), [OpenMPI](https://www.open-mpi.org) or any
other [MPI](https://www.mpi-forum.org) implementation that provides compiler
wrappers `mpicc` and `mpif90` for C and Fortran, respectively, as well as the
job launcher `mpiexec`.
    > **_NOTE:_** The job launcher of [OpenMPI](https://www.open-mpi.org) fails
to run more MPI processes than the number of real processor cores available on
the machine by default. That might lead to failures when configuring or running
ICON. The solution to the problem is to run the configure wrapper with an
additional argument `MPI_LAUNCH='mpiexec --oversubscribe'` (alternatively, you
can set the `OMPI_MCA_rmaps_base_oversubscribe` environment variable to `1`).

    > **_NOTE:_** It is not rare that the latest versions (or the default
versions available via the package managers) of
[OpenMPI](https://www.open-mpi.org) and [MPICH](https://www.mpich.org) are
affected with bugs that make the libraries unusable for ICON. A good way to make
sure that the MPI library does not have significant defects is to switch to the
root source directory of ICON and run the following commands (do not forget the
aforementioned extra arguments for the configure wrapper if you are using
[OpenMPI](https://www.open-mpi.org)):
    > ```bash
    > ./config/generic/gcc --enable-yaxt --enable-cdi-pio --enable-coupling
    > make -j4 check-bundled TESTS= XFAIL_TESTS=  # this step speeds up the next one but can be skipped
    > make check-bundled  # avoid running this step in parallel on a weak machine, i.e. omit the -j argument
    > test $? -eq 0 && echo "Everything is fine" || echo "Something went wrong"
    > ```
    > After that, you can clean up the source directory and reconfigure ICON the
way you need.
- [HDF5](https://support.hdfgroup.org/HDF5) with high-level interface (for
<a href="#netcdf-c">NetCDF-С</a>), thread-safety (for <a href="#cdo">CDO</a>),
and szlib filtering support (only C interface required, not a direct dependency
of ICON)
- <a name="netcdf-c"/> [NetCDF-C](https://www.unidata.ucar.edu/software/netcdf/docs)
with NetCDF-4 support
- [NetCDF-Fortran](https://www.unidata.ucar.edu/software/netcdf/docs-fortran)
- [BLAS](http://www.netlib.org/blas)
- [LAPACK](http://www.netlib.org/lapack)
- [ecCodes](https://confluence.ecmwf.int/display/ECC) with JPEG2000 and AEC
support (only C interface required)
- [libfyaml](https://github.com/pantoniou/libfyaml)
- [Libxml2](http://www.xmlsoft.org)

See section [ICON dependencies](../../README.md#icon-dependencies) for more
details.

### Optional tools

- <a name="cdo"/> [CDO](https://code.mpimet.mpg.de/projects/cdo) for pre- and
post-processing, also used by some of the
[generated runscripts](../../README.md#running)
- [rsync](https://rsync.samba.org/) for the
[generated runscipts](../../README.md#running) in the case of
[out-of-source building](../../README.md#out-of-source-configuration-building)

## Tested platforms and tools

This section provides basic instructions on how to install most commonly
required subset of [ICON dependencies](../../README.md#icon-dependencies) on
different operating systems using relevant
[package managers](https://en.wikipedia.org/wiki/Package_manager).

### macOS with [MacPorts](https://www.macports.org)

**Tested on `macOS Ventura 13.5.1`.**

Most of the required software packages are either already available on the
system or installed together with [Xcode](https://developer.apple.com/xcode) and
[Command Line Tools for Xcode](https://developer.apple.com/download/more/),
which are [prerequisites for MacPorts](https://www.macports.org/install.php).
The rest of the required software can be installed by running the following
commands:

```bash
# Install building tools and ICON dependencies:
sudo port -N install       \
  cmake                    \
  gcc12                    \
  mpich-gcc12              \
  hdf5 +hl+threadsafe+szip \
  netcdf                   \
  netcdf-fortran +gcc12    \
  eccodes                  \
  libxml2

# Install libfyaml, which is currently not available via MacPorts:
curl -OL https://github.com/pantoniou/libfyaml/releases/download/v0.8/libfyaml-0.8.tar.gz
tar xvf libfyaml-0.8.tar.gz
cd libfyaml-0.8
./configure --prefix=/opt/local
make -j
sudo make install

# The command above can be reverted as follows:
# curl -OL https://github.com/pantoniou/libfyaml/releases/download/v0.8/libfyaml-0.8.tar.gz
# tar xvf libfyaml-0.8.tar.gz
# cd libfyaml-0.8
# ./configure --prefix=/opt/local
# sudo make uninstall

# Select the compiler and MPI compiler wrappers:
sudo port select --set gcc mp-gcc12
sudo port select --set mpi mpich-gcc12-fortran
hash -r

# Install optional tools:
sudo port -N install cdo +netcdf
```
> **_NOTE:_** You can try replacing `mpich` with `openmpi` in the commands above
if the version of [MPICH](https://www.mpich.org) that is currently available via
[MacPorts](https://www.macports.org) fails the tests described in the
[Software libraries](#software-libraries) section. Note, however, that
[OpenMPI](https://www.open-mpi.org) is known to have problems with running on
macOS. Although the
[list of known issues](https://www.open-mpi.org/faq/?category=osx) is very
dated, some of them are still relevant. In particular,
[this one](https://www.open-mpi.org/faq/?category=osx#startup-errors-with-open-mpi-2.0.x)
(also see [here](https://github.com/open-mpi/ompi/issues/7393)). The solution
here is to run the configure wrapper with one more argument
`BUILD_ENV="export TMPDIR='/tmp';"` and make sure that the `TMPDIR` environment
variable is set to `/tmp` before running ICON.

### Ubuntu with [Apt](https://wiki.debian.org/Apt)

**Tested on `Ubuntu Jammy Jellyfish 22.04.3 LTS`.**

```bash
# Install building tools and ICON dependencies:
sudo apt install -y \
  build-essential   \
  cmake             \
  python3           \
  gcc               \
  gfortran          \
  libopenmpi-dev    \
  libhdf5-dev       \
  libnetcdf-dev     \
  libnetcdff-dev    \
  libeccodes-dev    \
  libblas-dev       \
  liblapack-dev     \
  libfyaml-dev      \
  libxml2-dev

# Select MPI libraries and compiler wrappers
sudo update-alternatives --set mpi /usr/bin/mpicc.openmpi
sudo update-alternatives --set mpirun /usr/bin/mpirun.openmpi

# If the two non-interactive commands above do not work,
# try the interactive analogues:
# sudo update-alternatives --config mpi
# sudo update-alternatives --config mpirun

# Install optional tools:
sudo apt install -y cdo
```

### Arch Linux with [Pacman](https://wiki.archlinux.org/index.php/pacman)

**Tested on `Arch Linux 2023.07.23`.**

```bash
# Install building tools and ICON dependencies:
sudo pacman -S --noconfirm \
  base-devel               \
  git                      \
  cmake                    \
  python                   \
  gcc-fortran              \
  openmpi                  \
  hdf5                     \
  netcdf                   \
  netcdf-fortran           \
  blas                     \
  lapack                   \
  libxml2

# Install libfyaml from the Arch User Repository:
( git clone https://aur.archlinux.org/libfyaml.git && cd libfyaml && makepkg -csi --noconfirm )
# Install ecCodes from the Arch User Repository:
( git clone https://aur.archlinux.org/eccodes.git && cd eccodes && makepkg -csi --noconfirm )

# Install optional tools:
sudo pacman -S --noconfirm rsync
# Install CDO and its extra dependencies from the Arch User Repository:
( git clone https://aur.archlinux.org/udunits.git && cd udunits && makepkg -csi --noconfirm )
( git clone https://aur.archlinux.org/magics++.git && cd magics++ && makepkg -csi --noconfirm )
( git clone https://aur.archlinux.org/cdo.git && cd cdo && makepkg -csi --noconfirm )
```

### macOS/Linux with [Spack](https://spack.io)

**Tested on `Ubuntu Desktop Jammy Jellyfish 22.04.3 LTS`.**

> **_NOTE:_** Spack has its own list of
[prerequisites](https://spack.readthedocs.io/en/latest/getting_started.html#prerequisites).
The default Ubuntu Desktop installation seems to be missing only few of them:
>```bash
>sudo apt install -y build-essential git gfortran
>```

```bash
# Install Spack:
git clone https://github.com/spack/spack.git
. ./spack/share/spack/setup-env.sh

# Find compilers and building tools that are already installed on the system:
spack compiler find
spack external find

# Install ICON dependencies:
spack install                             \
  cmake                                   \
  openmpi                                 \
  netcdf-fortran ^hdf5+hl+szip+threadsafe \
  eccodes                                 \
  netlib-lapack                           \
  libfyaml                                \
  libxml2

# Symlink the dependencies to a single prefix (e.g. to $HOME/icon-sw):
export ICON_SW_PREFIX="$HOME/icon-sw"
spack view symlink -i "$ICON_SW_PREFIX" \
  cmake                                 \
  openmpi                               \
  netcdf-fortran                        \
  eccodes                               \
  netlib-lapack                         \
  libfyaml                              \
  libxml2

# Install optional tools:
spack install cdo
```
