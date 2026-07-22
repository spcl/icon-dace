# NVHPC + netCDF + VT-library build paths for the icon-gpu spack env.
# Source before configuring/building ICON on ault:
#     spack env activate icon-gpu
#     source config/cscs/ault-build-env.sh
# Override VT_DIR to point at the VT velocity/ tree if it is not the sibling
# icon-vt-dace/velocity of this icon-dace checkout.

env_loc() { spack -e icon-gpu location -i "$1"; }
_nvroot="$(env_loc nvhpc)/Linux_x86_64/26.1"
_ncf="$(env_loc netcdf-fortran)"
_ncc="$(env_loc netcdf-c)"
_cmk="$(env_loc cmake)"
_vtdir="${VT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../icon-vt-dace/velocity" && pwd)}"
_vsql="$(spack -e vt-gpu find --format '{prefix}' sqlite | head -1)"
_vzst="$(spack -e vt-gpu find --format '{prefix}' zstd | head -1)"

export PATH="${_nvroot}/comm_libs/mpi/bin:${_nvroot}/compilers/bin:${_cmk}/bin:${PATH}"
export LIBRARY_PATH="${_nvroot}/cuda/lib64:${_nvroot}/compilers/lib:${_ncf}/lib:${_ncc}/lib:${_vsql}/lib:${_vzst}/lib:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="${LIBRARY_PATH}:${_vtdir}"
export CPATH="${_nvroot}/cuda/include:${_ncf}/include:${_ncc}/include:${CPATH:-}"
export DACE_FCFLAGS="-I${_ncf}/include -I${_ncc}/include"
export DACE_LDFLAGS="-L${_vtdir} -Wl,-rpath,${_vtdir} -Wl,-rpath-link,${_vsql}/lib -Wl,-rpath-link,${_vzst}/lib"
export DACE_LIBS="-lvelocity"
