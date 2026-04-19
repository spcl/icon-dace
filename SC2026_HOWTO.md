# SC2026 — ICON integration

## 1. Get the code
```bash
git clone --branch prat/sc2026 git@github.com:spcl/icon-dace.git
cd icon-dace
git submodule update --init --recursive
```

Copy the VT-generated wrapper (see velocity `SC2026_HOWTO.md` step 3).
Assumes the VT tree at `../icon-vt-dace/`:
```bash
cp ../icon-vt-dace/wrapper.f90 src/atm_dyn_iconam/wrapper.f90
```
`wrapper.f90` = `velocity_tendencies_gpu` (ICON→VT dispatch). ABI must
match the `.so`. All three precisions currently produce an identical
wrapper, so one copy covers fp64/fp32/fp16 — regenerate only on VT
build-flag changes that alter the ABI.

`src/atm_dyn_iconam/serde.f90` (`module vt_serde`) is tracked here and is
canonical. **Do not overwrite it from VT** — VT's copy lags.
<!-- TODO: streamline serde.f90 into a fully auto-generated module. -->

## 2. Prereqs
- One-time: `uenv image pull icon/25.2:v1@santis`
- Every session: `uenv start --view=default icon/25.2:v1@santis`
- `spack` activated in the shell (`spack --version` should work).
  Packages load ad-hoc during build.

## 3. Build directory
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
```bash
curl -LO http://icon-downloads.mpimet.mpg.de/grids/public/edzw/icon_grid_0010_R02B04_G.nc
```

`GRID=<gridID>_<refinement>` maps to `icon_grid_<GRID>_G.nc`. DWD catalog:

| `GRID`        | File                            |
|---------------|---------------------------------|
| `0010_R02B04` | `icon_grid_0010_R02B04_G.nc` (default) |
| `0008_R02B05` | `icon_grid_0008_R02B05_G.nc`    |
| `0002_R02B06` | `icon_grid_0002_R02B06_G.nc`    |
| `0050_R02B03` | `icon_grid_0050_R02B03_G.nc`    |

Drop the `.nc` for any other `GRID` you plan to submit with.

## 7. Runscript
`run/exp.sc2026` is the template; `make_sc2026_runscript.sh` generates
+ patches the `.run`:
```bash
./run/make_sc2026_runscript.sh
# → run/exp.sc2026.run
```

## 8. Submit

```
./run/sbatch_sc2026.sh <ATM_TIMESTEP> [VT_PREC=fp64] [GRID=0010_R02B04]
```
Caller env passes through, so vt_serde knobs work via `--export=ALL`.

### vt_serde knobs

| Env var                   | Default  | Effect |
|---------------------------|----------|--------|
| `USE_VT_GPU`              | `.true.` | false → ICON's built-in `velocity_tendencies`; true → VT's `libvelocity.so` |
| `NDYN_SUBSTEPS_OVERRIDE`  | `0`      | 0 = leave namelist (5) alone; `N>0` = force `N` |
| `SERDE_GEN_START`         | `0`      | first `physics_generation` where `do_serialize` flips on |
| `SERDE_GEN_END`           | `51`     | `physics_generation` where `do_serialize` flips off |

All five variants below set `SERDE_GEN_END=51` (explicit default) —
`.data` dumps for `physics_generation 0..50` are produced. Set
`SERDE_GEN_END=0` for pure-timing runs with no dumps.

### 8.1–8.5 Variants

VT variants need `libvelocity_gpu_stage8_solve_nh_integration_release.<VT_PREC>.so`
under `${VT_DIR}` (default `/capstor/scratch/cscs/pmazumde/sc2026-ad-test/icon-vt-dace`).

| # | Purpose | Command (explicit GRID shown on 8.1; same positional arg for the rest) |
|---|---|---|
| 8.1 | Vanilla, 5 substeps (got reference)  | `USE_VT_GPU=0 NDYN_SUBSTEPS_OVERRIDE=5  SERDE_GEN_END=51 ./run/sbatch_sc2026.sh 8 fp64 0010_R02B04` |
| 8.2 | Vanilla, 10 substeps (temporal-refined reference) | `USE_VT_GPU=0 NDYN_SUBSTEPS_OVERRIDE=10 SERDE_GEN_END=51 ./run/sbatch_sc2026.sh 8 fp64 0010_R02B04` |
| 8.3 | VT fp64                              | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  SERDE_GEN_END=51 ./run/sbatch_sc2026.sh 8 fp64 0010_R02B04` |
| 8.4 | VT fp32                              | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  SERDE_GEN_END=51 ./run/sbatch_sc2026.sh 8 fp32 0010_R02B04` |
| 8.5 | VT fp16                              | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  SERDE_GEN_END=51 ./run/sbatch_sc2026.sh 8 fp16 0010_R02B04` |

Swap the last positional for a different grid tag (`0008_R02B05`,
`0002_R02B06`, `0050_R02B03`) after dropping the matching `.nc` into
`build/verification/`. For 8.1/8.2 (`USE_VT_GPU=0`) the precision/grid
args only shape `EXPNAME` and log filenames — no `.so` is loaded.

Each variant writes to `experiments/<EXPNAME>/` where
`EXPNAME=sc2026_dt<dt>_{vanilla|gpu<prec>}_<grid>`.

### Serialized data
Files land in `build/verification/experiments/${EXPNAME}/` as
`<field>.p<phys>.d<dycore>.vt<vt>.ss<substep>.data`
(e.g. `p_prog.t0.p0.d1.vt1.ss0.data`). Pair 8.1 (vanilla) with 8.3–8.5
(VT) for `got`/`want` comparison.

Logs: `run/LOG.SAVEME-SER-<PREC>.<EXPNAME>.<jobid>.o`.
