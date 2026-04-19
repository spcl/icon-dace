# SC2026 — ICON integration

## 1. Get the code
```bash
git clone --branch prat/sc2026 git@github.com:spcl/icon-dace.git
cd icon-dace
# if the build pulls in externals/submodules, initialize them:
git submodule update --init --recursive
```

Pull in the VT-generated Fortran wrapper (see the velocity
`SC2026_HOWTO.md` step 3 for how it's produced). Assuming the VT tree
sits alongside this repo as `../icon-vt-dace/`:
```bash
cp ../icon-vt-dace/wrapper.f90 src/atm_dyn_iconam/wrapper.f90
```
`wrapper.f90` provides `velocity_tendencies_gpu` (the ICON→VT dispatch).
Its ABI must match the `.so`. In practice the three precisions currently
produce an identical `wrapper.f90`, so one copy covers all of them — but
regenerate + recopy if you change the VT build flags in a way that could
alter the ABI.

`src/atm_dyn_iconam/serde.f90` (`module vt_serde`, used by the modified
`mo_nh_stepping` / `mo_solve_nonhydro` / `mo_velocity_advection`) is
already tracked in this repo — it's the canonical copy for SC2026. VT's
`velocity/serde.f90` is a lagging mirror; **do not overwrite icon-dace's
`serde.f90` by syncing from VT**.
<!-- TODO: streamline serde.f90 into a full auto-generated module; it's
     currently a mix of auto-generated bulk (W_*, ctor_*) and a hand-written
     API block (vt_generation / dycore_generation / *_tic / do_serialize /
     runtime knobs). -->


## 2. Prereqs
- On the cluster:
  - One-time: pull the ICON uenv image to your local repo (auto-created
    on first use under `/capstor/scratch/cscs/$USER/.uenv-images`):
    ```bash
    uenv image pull icon/25.2:v1@santis
    ```
  - Every session: start the environment with the ICON view:
    ```bash
    uenv start --view=default icon/25.2:v1@santis
    ```
  - A working `spack` install, already activated in the shell (path/version
    depends on the user — `spack --version` should return something sensible).
    No mandatory `spack load` yet; packages get loaded ad-hoc as needed
    during build.

## 3. Build directory
Create and enter the verification build tree:
```bash
mkdir -p build/verification
cd build/verification
```

## 4. Configure
```bash
../../config/cscs/daint_yakup_dace.gpu.gh200.nvidia
```

## 5. Build
```bash
make -j
```

## 6. Grid file
ICON needs the grid file under `build/verification/`:
```bash
curl -LO http://icon-downloads.mpimet.mpg.de/grids/public/edzw/icon_grid_0010_R02B04_G.nc
# → icon_grid_0010_R02B04_G.nc
```

Grids are identified by a `<gridID>_<refinement>` tag that gets spliced into
`icon_grid_<tag>_G.nc`. Common choices (from the DWD catalog):
| `GRID` tag    | File                            |
|---------------|---------------------------------|
| `0010_R02B04` | `icon_grid_0010_R02B04_G.nc` (default) |
| `0008_R02B05` | `icon_grid_0008_R02B05_G.nc`    |
| `0002_R02B06` | `icon_grid_0002_R02B06_G.nc`    |
| `0050_R02B03` | `icon_grid_0050_R02B03_G.nc`    |

If you plan to run a different `GRID` at submit time, drop the corresponding
`.nc` file here too (same `curl -LO` pattern, different filename).

## 7. Runscript
Template `run/exp.sc2026` is tracked; the wrapper `run/make_sc2026_runscript.sh`
regenerates + patches the generated `.run`:
```bash
./run/make_sc2026_runscript.sh
# → run/exp.sc2026.run   (SBATCH header + ulimit + ATM_TIMESTEP guard patched in)
```

## 8. Submit

The submit wrapper expects the VT shared library at:
```
${VT_DIR}/libvelocity_gpu_stage8_solve_nh_integration_release.${VT_PREC}.so
```
Default `VT_DIR=/capstor/scratch/cscs/pmazumde/sc2026-ad-test/icon-vt-dace`
(override via env). Default `VT_PREC=fp64`. So with defaults it looks for:
```
/capstor/scratch/cscs/pmazumde/sc2026-ad-test/icon-vt-dace/libvelocity_gpu_stage8_solve_nh_integration_release.fp64.so
```
Build this via the velocity HOWTO's integration step.

```bash
./run/sbatch_sc2026.sh <ATM_TIMESTEP> [VT_PREC=fp64] [GRID=0010_R02B04]
# e.g.
./run/sbatch_sc2026.sh 2                     # dt=2, fp64, default grid
./run/sbatch_sc2026.sh 8 fp16 0008_R02B05    # dt=8, fp16, different grid
VT_DIR=/other/path ./run/sbatch_sc2026.sh 2 fp32
```

### vt_serde runtime knobs
The wrapper passes the caller's env through (`sbatch --export=ALL,...`), so
these env vars reach `vt_serde_init` inside the job:

| Env var                   | Default | SC2026 value | Effect |
|---------------------------|---------|--------------|--------|
| `NDYN_SUBSTEPS_OVERRIDE`  | 10      | 5–10         | sets `ndyn_substeps_var(jg)` inside the instrumented window |
| `SERDE_GEN_START`         | 0       | 0            | first `physics_generation` where `do_serialize` flips on |
| `SERDE_GEN_END`           | 51      | 51           | `physics_generation` where `do_serialize` flips off |

Defaults already match the paper's (10, 0, 51) configuration, so for the
standard SC2026 experiment you don't need to set anything:
```bash
./run/sbatch_sc2026.sh 2                     # dt=2, fp64, 10 substeps, gens 0..51
```
For the lower-substep variant:
```bash
NDYN_SUBSTEPS_OVERRIDE=5 ./run/sbatch_sc2026.sh 2 fp32
```

Logs: `run/LOG.SAVEME-SER-<PREC>.sc2026_dt<step>_<prec>_<grid>.<jobid>.o`.
