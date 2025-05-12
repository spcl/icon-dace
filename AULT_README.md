```bash
# cd /path/to/dace
# git checkout f2dace/staging
# cd /path/to/icon-dace/
git checkout serde
cd icon-model/
mkdir -p build/verification
# alias spack=/scratch/ybudanaz/spack/bin/spack7
spack load gcc@13
spack load icon@2024.10%gcc
spack load openblas%gcc
export LIBRARY_PATH=$(spack location -i netcdf-fortran%gcc)/lib:${LIBRARY_PATH}
export LIBRARY_PATH=$(spack location -i openblas)/lib:${LIBRARY_PATH}
export LDFLAGS="-L$(spack location -i netcdf-fortran%gcc)/lib -L$(spack location -i openblas)/lib"
export CC=gcc
export CXX=g++
export F90=gfortran
export FC=gfortran
alias fc=gfortran
alias f90=gfortran
../../config/cscs/ault_ben_dace.cpu.gcc # Ignore rsync error I guess
make -j
```