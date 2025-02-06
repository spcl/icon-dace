# Commands

Reset
```
sudo docker rm icon_cpu_container
git clone --recursive git@spclgitlab.ethz.ch:benweber/icon-dace-integration.git
cd icon-dace-integration
git clone --recursive git@github.com:spcl/dace.git
cd dace
git checkout multi_sdfg
cd ..
sudo docker run --name icon_gpu_container --gpus all -it -v /home/primrose/Work/IconGrounds/icon-dace-integration:/home/icon/icon icon_nvhpc_rev1:latest
apt update -y && apt upgrade -y
#chmod -R 777 /
su icon
cd ~
python3 -m venv .venv
source .venv/bin/activate
cd icon
pip install --upgrade setuptools wheel pip
pip install -e ./dace[testing]
```

sudo docker run --name icon_gpu_container --gpus all -it -v /home/primrose/Work/IconGrounds/icon-dace:/home/icon/icon icon_nvhpc_rev1:latest

Start ICON CPU
```
sudo docker start -ai icon_gpu_container
```

Attach 2nd terminal
```
sudo docker exec -it icon_cpu_container /bin/bash
```

GOTO ICON:
```
su icon
cd /home/icon/icon
source ~/.venv/bin/activate
```

Build example:
```
mkdir -p build/verification/
cd build/verification/
../../config/generic/gcc \
    CC=mpicc \
    FC=mpif90 \
    FCFLAGS="-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -lnetcdff" \
    LDFLAGS="-L/usr/lib/x86_64-linux-gnu/ " \
    LIBS="-leccodes -lnetcdff -lnetcdf -lopenblas" \
    --enable-acm-license \
    --disable-mixed-precision \
    --disable-edmf \
    --disable-les \
    --disable-ocean \
    --disable-jsbach \
    --disable-coupling \
    --disable-aes \
    --disable-rte-rrtmgp \
    --enable-ecrad \
    --disable-mpi \
    --disable-mpi-checks \
    --disable-openmp \
    --disable-loop-exchange \
    --enable-dace-subst=no \
    --enable-explicit-fpp

../../config/generic/gcc     CC=mpicc     FC=mpif90     FCFLAGS="-g -I/usr/include/ -fcheck=all -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function"     LDFLAGS="-L/usr/lib/x86_64-linux-gnu/ "     LIBS="-leccodes -lnetcdff -lnetcdf -lopenblas"     --enable-acm-license     --disable-mixed-precision     --disable-edmf     --disable-les     --disable-ocean     --disable-jsbach     --disable-coupling     --disable-aes     --disable-rte-rrtmgp     --enable-ecrad     --disable-mpi     --disable-mpi-checks     --disable-openmp     --disable-loop-exchange     --enable-dace-subst=verify     --enable-explicit-fpp --disable-upperatmo -without-external-rte-rrtmgp --disable-rte-rrtmgp

make -j16

./make_runscripts exclaim_ape_R02B04_integration_demo_config

# Do not call outside the build pipeline
# python ../../sdfgs/generate.py  edges_to_cells_bilinear_interpolation ../../sdfgs/integrations.yaml .

cd ./run

./exp.exclaim_ape_R02B04_integration_demo_config.run
```

Skip:
vertical_levels
conv_tracer
turb_tracer
generated cells->edges is problematic

```
```
./pipeline.sh
cd icon-scratchpad/icon-model
mkdir -p build/verification/
cd build/verification/
export SDFG_LIB_PATHS="-L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build"
export SDFG_LIBS="-lradiation"
../../config/generic/gcc \
    CC=mpicc \
    FC=mpif90 \
    FCFLAGS="-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE" \
    LDFLAGS="-L/usr/lib/x86_64-linux-gnu/ ${SDFG_LIB_PATHS}" \
    LIBS="-leccodes -lnetcdff -lnetcdf -lopenblas ${SDFG_LIBS}" \
    --enable-acm-license \
    --disable-mixed-precision \
    --disable-edmf \
    --disable-les \
    --disable-ocean \
    --disable-jsbach \
    --disable-coupling \
    --disable-aes \
    --disable-rte-rrtmgp \
    --enable-ecrad \
    --disable-mpi \
    --disable-mpi-checks \
    --disable-openmp \
    --disable-loop-exchange \
    --enable-dace-subst=no \
    --enable-explicit-fpp
make -j16
```