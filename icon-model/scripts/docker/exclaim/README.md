# ICON-EXCLAIM Containerization:  a path for Continuous Integration and portable high-performance ICON simulations

ICON-EXCLAIM is based on a public release of the ICOsahedral Non-hydrostatic (ICON) model for climate simulations and numerical weather prediction.  
It is being rewritten in the EXCLAIM project to use a workflow written in Python (as opposed to Fortran currently) calling components (or "granules" 
which are also written in Python whose computationally intensive kernels are formulated in a Python-subset domain specific language (DSL) called 
GridTools for Python (GT4Py).  

ICON-EXCLAIM will be deployed on the new "Alps" supercomputer at the Swiss National Supercomputing Centre.  This national flagship platform is split 
into virtual clusters.  Software portability and testing of all granules in the simulation system within a Continuous Integration (CI) framework is 
key.  To this end Containers are key:  these allow us to run the software within an environment based on a particular version of a compiler 
and a well-defined software stack.  The container can be moved to other platforms and run within a container run-time system, e.g., Docker or Sarus, 
and should produce identical (or very similar) results.  Alternatively, one can build containers based on older images of the software stack which
are no longer provided on the target platform.  This allows a model like ICON, whose results are very sensitive to the software stack and compiler 
versions, to run in a stable and predictable manner for extended periods of time, which ensures the reproducibility of scientific results over
a multi-year period.

The Alps CI system also requires that all software tested resides in a container.  Since CI for individual granules and for the ICON model as a 
whole, it is a fundamental requirement to containerize ICON.  ICON is a complex software which relies on source code from multiple repositories 
all over the global.  As such it provides a complex test case for the evolving CI system.  In the subsequent sections we discuss the layout of
the ICON software, the approach taken for its containerization, which include a number of additional repositories.

## Software ecosystem

Many software repositories are required for the ICON-EXCLAIM build and containerization.   


-	github.com/C2SM/container	:  This (public) repository, where the Dockerfiles result.  This is public.  

- github.com/C2SM/spack-c2sm : a public repository which contains a version of the Spack package manager tailored for C2SM users of both COSMO and ICON

-	github.com/C2SM/icon-exclaim :  a private repository where the ICON source resides.  This is accessed through SSH keys defined for registered users.

ICON-EXCLAIM relies on many GIT submodules (sometimes called 'externals') which reside in 

The git submodules:  come from various repositories, some public, some private.  In particular, URLs with the prefix '../..' are within the main ICON repository,  gitlab.dkrz.de:, and thus are all private.  The same SSH keys as ICON apply for these private repositories.

- [submodule "externals/mtime"] url = git@github.com:C2SM/libmtime.git
- [submodule "externals/jsbach"] url = git@github.com:C2SM/jsbach/jsbach.git
- [submodule "externals/yac"]	url = https://gitlab.dkrz.de/dkrz-sw/YAC.git
- [submodule "externals/tixi"] url = git@github.com:C2SM/icon-libraries/libtixi.git
- [submodule "externals/yaxt"] url = https://gitlab.dkrz.de/dkrz-sw/yaxt.git 
- [submodule "externals/rte-rrtmgp"] url = https://github.com/earth-system-radiation/rte-rrtmgp.git
- [submodule "externals/cub"] url = https://github.com/NVlabs/cub.git
- [submodule "externals/cdi"]	url = https://gitlab.dkrz.de/mpim-sw/libcdi.git 
- [submodule "externals/sct"]	url =https://gitlab.dkrz.de/dkrz-sw/sct.git
- [submodule "externals/ecrad"]	url = git@github.com:C2SM/libecrad.git
- [submodule "externals/dace_icon"]	url = git@github.com:C2SM/dace-icon-interface.git
- [submodule "externals/emvorado"]url = git@github.com:C2SM/emvorado-for-icon.git
- [submodule "externals/probtest"]url = git@github.com:C2SM/probtest.git
- [submodule "utils/mkexp"] url = https://git.mpimet.mpg.de/public/mkexp
- [submodule "externals/art"]	url = git@github.com:C2SM/art/art.git
- [submodule "externals/ppm"]	url = git@github.com:C2SM/ppm.git

All the private repositories have been mirrored from their original repositories on github.com:C2SM so that any user who can access the icon-exclaim.git repository can also aceess all submodules.

A future step will be to incorporate the Docker files and the documentation into the icon-exclaim.git repository, such that it can create containers without further clones.  This is a requirement of the CSCS-CI system.

## Spack build



## Container stages

ICON has extensive package dependencies, e.g., LAPACK, BLAS, MPI, NetCDF, HDF5, eccodes libraries.  There is distinct advantage of building the dependencies in ***stages*** so that previously built containers can be 'checkpointed':  this will ensure one need not keep rebuilding 'from scratch'.  After creating as many as 5 stages during development, we settled on two (with two variants of the application container):

- icon-dependencies-mpich : this container contains all the above-mentioned ICON dependencies.  Since these rarely change, this container can usually be reused when the ICON container is built.  This container is agnostic to GPUs.  

- icon-mpich : this container depends on the previous one, and wraps the icon application (and its submodules).  We distinguish icon-mpich from a rarely used (and unsupported) icon-serial variant, which is has no communication and must be run with a single process (possibly with multi-threading).

- icon-mpich-gpu : this container also depends on icon-dependencies-mpich and builds the application container to be run on GPUs.

In other words, a runnable container would require first building icon-dependencies-mpich, and then either icon-mpich or icon-mpich-gpu.  

## Build Instructions

The docker image can naturally be built with ```docker```, but also with OCI-compatible builders, such as Buildah  (https://buildah.io/).  Due to the requirement for root access, the former presents security risks on HPC systems and is virtually never available.  On the other hand, one can easily build images on a personal computer and copy them subsequently to the HPC system to be run with a Docker-compatible framework, such as Sarus (https://products.cscs.ch/sarus/).  Since the latter only builds the image, there are no such security risks, and it should build the image faster than a laptop.  

### Docker build

One starts by checking out the master branch of this repository:

```
git clone git@github.com:C2SM/container
cd container/icon-exclaim
```

The container for the dependencies of ICON is built first:

```
cd icon-dependencies-mpich
docker build -t $(cat TAG) --ssh default .
```

In order to build the application container, the user will need access to the icon-exclaim.git repository which is private due to licensing restrictions.  It is possible that the SSH keys needed for the icon-exclaim have to be scanned at this point:

```
eval $(ssh-agent) > /dev/null
ssh-add ~/.ssh/<private_key>
```

At this point it should be possible to build the CPU or GPU icon containers, e.g., for GPU:

```
cd ../icon-mpich-gpu
docker build -t $(cat TAG) --ssh default .
```

The image has to be saved, e.g.,

```
docker save c2sm/nvhpc:21.3-devel-cuda_multi-ubuntu20.04-icon-mpich -o my_icon.tar
```

This can be run on the PC through docker or copied to an HPC platform for execution there.

### Buildah build on CSCS Infrastructure

Buildah is a tool for building OCI-compatible images. It doesn't depend on a daemon such as Docker and therefore doesn't require root privileges. Buildah provides a command-line tool that replicates all the commands found in a Dockerfile. It is the container builder of choice on CSCS infrastructure. One should first consult the Buildah page https://user.cscs.ch/tools/containers/buildah/ for key instructions, in particular to define the appropriate ```$HOME/.config/containers/storage.conf``` file, e.g.,

```
[storage]
  runroot = "/scratch/local/<username>/runroot"
  graphroot = "/scratch/local/<username>/root"
```

One starts by checking out the master branch of this repository:

```
git clone git@github.com:C2SM/container
cd container/icon-exclaim
```

Due to CSCS limitations, the image must be built on a compute node:
```
salloc -N1 --time=04:00:00 -C "gpu&contbuild" -A csstaff
ssh $SLURM_NODELIST
```

The container for the dependencies of ICON is built first:

```
cd icon-dependencies-mpich
buildah bud --format=docker --tag $(cat TAG)
```

At this point one can optionally archive the dependency image, which can then be loaded ("pulled") at a later time during the build of the application container:

```
buildah push $(cat TAG) docker-archive:/scratch/snx3000/username/icon-dependencies-mpich.tar
```

or one can proceed immediately to complete the container for the application itself and archive it for subsequent execution.  It is possible that the SSH keys needed for the icon-exclaim have to be scanned at this point:

```
eval $(ssh-agent) > /dev/null
ssh-add ~/.ssh/<private_key>
```

For the CPU container:

```
cd ../icon-mpich
buildah bud --format=docker --tag $(cat TAG)
buildah push $(cat TAG) docker-archive:/scratch/snx3000/username/icon-mpich.tar
```

There are separate dockerfiles (largely for clarity, since the differences are minor) for the GPU:

```
cd ../icon-mpich-gpu
buildah bud --format=docker --tag $(cat TAG)
buildah push $(cat TAG) docker-archive:/scratch/snx3000/username/icon-mpich-gpu.tar
```

### Container execution on CSCS Infrastructure

Due to the above-mentioned security limitations of Docker, CSCS has built the OCI-compliant container-engine Sarus (https://products.cscs.ch/sarus/) to run Docker containers safely on HPC infrastructure.   First the container --- whether generated with Docker on a personal computer, or with Buildah --- must be loaded.

```
module load sarus
sarus load /scratch/snx3000/username/icon-mpich-gpu.tar icon:latest
```

The typical run scripts  require some subtle modifications in particular in the START and MODEL environment variables.  It is worthwhile to construct a new script, e.g.,  based on the existing one.  For example, exp.mch_bench_r19b07_dev.run might be modified to exp.mch_bench_r19b07_dev_sarus.run

```
> builddir=/scratch/snx3000/wsawyer/ICON/icon-exclaim/build_nvhpc_gpu
> experiments_dir=$builddir/experiments/$EXPNAME
> module load sarus
< export START="srun -n $mpi_total_procs --ntasks-per-node $mpi_procs_pernode --threads-per-core=1 --cpus-per-task $OMP_NUM_THREADS"
< export MODEL="${basedir}/bin/icon"
---
> export START="srun -n $mpi_total_procs --ntasks-per-node $mpi_procs_pernode --threads-per-core=1 --cpus-per-task $OMP_NUM_THREADS sarus run --mpi --mount=type=bind,src=$icon_data_rootFolder,destination=$icon_data_rootFolder --mount=type=bind,src=$thisdir,destination=$thisdir --mount=type=bind,src=$builddir,destination=$builddir --workdir=$experiments_dir load/library/icon:latest"
> export MODEL="bash -c 'MPICH_RDMA_ENABLED_CUDA=1 LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libcuda.so icon'"
```

There are several locations in the standard ICON scripts where file system checks have to be disabled, since the container has its own file system:

```
> ###ls -ld ${EXPDIR}
> ###check_error $? "${EXPDIR} does not exist?"
> ###ls -l ${MODEL}
> ###check_error $? "${MODEL} does not exist?"
```

With these changes, the script can be launched in the standard way:

```
sbatch exp.mch_bench_r19b07_dev_sarus.run
```










