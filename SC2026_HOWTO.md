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

Wrapper signature:
```
./run/sbatch_sc2026.sh <ATM_TIMESTEP> [VT_PREC=fp64] [GRID=0010_R02B04]
```
Submits `run/exp.sc2026.run` via sbatch. Caller env is passed through, so
all vt_serde knobs (below) can be exported on the command line.

### vt_serde runtime knobs

| Env var                   | Default   | Effect |
|---------------------------|-----------|--------|
| `USE_VT_GPU`              | `.true.`  | if false, use ICON's built-in `velocity_tendencies` (CPU/OpenACC) instead of dispatching to VT's `libvelocity.so` |
| `NDYN_SUBSTEPS_OVERRIDE`  | `10`      | `ndyn_substeps_var(jg)` inside the instrumentation window |
| `SERDE_GEN_START`         | `0`       | first `physics_generation` where `do_serialize` flips on |
| `SERDE_GEN_END`           | `51`      | `physics_generation` where `do_serialize` flips off |

Defaults (10 substeps, VT-GPU on) are what we've been using for SC2026
runs so far. The four variants below all pass `SERDE_GEN_END=0` — this
**disables serialization** (no `.data` dumps), which is the right default
for timing runs. To *enable* serialization, remove `SERDE_GEN_END=0`
(default is `51`, meaning dumps fire for `physics_generation` 0..50). See
"**Collecting serialized data**" below.

### 8.1 Vanilla ICON (reference run)
No `libvelocity.so`. Uses ICON's built-in `velocity_tendencies`. Used to
generate the reference `got` files that the VT runs are validated against.
```bash
USE_VT_GPU=0 SERDE_GEN_END=0 ./run/sbatch_sc2026.sh 8 fp64
```
(`fp64`/grid values here just drive `EXPNAME` and log naming — no `.so` is
loaded when `USE_VT_GPU=0`.)

### 8.2 VT fp64
Requires `libvelocity_gpu_stage8_solve_nh_integration_release.fp64.so`
under `${VT_DIR}` (default
`/capstor/scratch/cscs/pmazumde/sc2026-ad-test/icon-vt-dace`).
```bash
USE_VT_GPU=1 SERDE_GEN_END=0 ./run/sbatch_sc2026.sh 8 fp64
```

### 8.3 VT fp32
Requires `...release.fp32.so`.
```bash
USE_VT_GPU=1 SERDE_GEN_END=0 ./run/sbatch_sc2026.sh 8 fp32
```

### 8.4 VT fp16
Requires `...release.fp16.so`.
```bash
USE_VT_GPU=1 SERDE_GEN_END=0 ./run/sbatch_sc2026.sh 8 fp16
```

Each variant lands in a **separate `experiments/<EXPNAME>/`** directory —
`sbatch_sc2026.sh` builds `EXPNAME=sc2026_dt<dt>_{vanilla|gpu<prec>}_<grid>`
so runs don't overwrite each other's serialized dumps.

### Collecting serialized data
To actually write `.data` dumps (for validation against the reference run),
drop `SERDE_GEN_END=0` — the default `SERDE_GEN_END=51` enables serialization
for `physics_generation = 0..50`. Combine with the reference (vanilla) run
so you have a `got` vs `want` pair.

```bash
# reference (vanilla ICON, full 51-gen dumps)
USE_VT_GPU=0 ./run/sbatch_sc2026.sh 8 fp64

# corresponding VT fp32 run (same gens, same grid, same dt → same filenames)
./run/sbatch_sc2026.sh 8 fp32
```

**Where the dumps land:** each file is written relative to the job's CWD,
which the generated `exp.sc2026.run` sets to the experiment directory:
```
build/verification/experiments/${EXPNAME}/
```
Filenames follow `<field>.p<phys>.d<dycore>.vt<vt>.ss<substep>.data`, e.g.
`p_patch.p0.d1.vt1.ss0.data`, `p_prog.t0.p0.d1.vt1.ss0.data`, etc.

### Lower-substep variant
Any of the above can be run with 5 substeps instead of 10:
```bash
NDYN_SUBSTEPS_OVERRIDE=5 SERDE_GEN_END=0 ./run/sbatch_sc2026.sh 8 fp32
```

Logs: `run/LOG.SAVEME-SER-<PREC>.sc2026_dt<step>_<prec>_<grid>.<jobid>.o`.
