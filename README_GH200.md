# ICON Verification Pipeline

The complete command then looks like this:
```sh
uenv start --view=default icon/25.2:v1@santis
export _RELEASE=TRUE
export _USE_NVHPC=TRUE
export _USE_CUDA_EVENTS=FALSE
export _CUDA_ARCH=native
cd /capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity/
python -m utils.stages.compile_gpu_stage9
cd /capstor/scratch/cscs/ybudanaz/icon-dace/icon-model/build/verification/run

mkdir -p /capstor/scratch/cscs/ybudanaz/icon-dace/icon-model/build/verification
cd /capstor/scratch/cscs/ybudanaz/icon-dace/icon-model/build/verification
export DACE_LDFLAGS="-L/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE"
export DACE_LIBS="-lvelocity_gpu_stage9"
export DACE_FCFLAGS="-DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE -L/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity -Wl,-rpath,/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity"
export LD_LIBRARY_PATH="/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity:${LD_LIBRARY_PATH}"
export F90=nvfortran
export CC=nvc
export CXX=nvc++
../../config/cscs/daint_yakup_dace.gpu.gh200.nvidia
make -j
cp -R ../../grids .
./make_runscripts --all
cd run
```