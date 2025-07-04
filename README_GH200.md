# ICON Verification Pipeline

The complete command then looks like this:
```sh
uenv start --view=default icon/25.2:v1@santis
export _RELEASE=TRUE
export _USE_NVHPC=TRUE
export _USE_CUDA_EVENTS=FALSE
export DACE_LDFLAGS="-L/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE"
export DACE_LIBS="-lvelocity_gpu_stage9"
export DACE_FCFLAGS="-L/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity -Wl,-rpath,/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity"
export LD_LIBRARY_PATH="/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity:${LD_LIBRARY_PATH}"
export F90=nvfortran
export CC=nvc
export CXX=nvc++
../../config/cscs/clariden_ben_dace.gpu.gh200.nvidia
make -j16
```