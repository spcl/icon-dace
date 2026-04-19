# SC2026 — ICON integration

This HOWTO is the **data-generation half** of the SC2026 reproduction.
Its output = serialized `.data` files for all 5 variants × 4 grids. For
analysis (SNR comparisons, table population) return to VT's
`SC2026_HOWTO.md §5`.

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

`sbatch_all_sc2026.sh <GRID>` fires all 5 variants (§8.1–8.5) for one grid:

```bash
SERDE_GEN_END=51 ./run/sbatch_all_sc2026.sh 0050_R02B03   # 320 km
SERDE_GEN_END=51 ./run/sbatch_all_sc2026.sh 0010_R02B04   # 160 km
SERDE_GEN_END=11 ./run/sbatch_all_sc2026.sh 0008_R02B05   #  80 km — trim
SERDE_GEN_END=6  ./run/sbatch_all_sc2026.sh 0002_R02B06   #  40 km — trim hard
```

Drop the matching `icon_grid_<GRID>_G.nc` into `build/verification/`
before submitting.

`SERDE_GEN_END` sets the serialization window (physics_generation
`0..N-1` are dumped). Each refinement ≈4× more cells = 4× fatter `.data`
files, so the finer grids get trimmed windows above. Paper-table numbers
only need the first physics step; bump to `51` if you want the full
50-step curve for SNR-evolution plots.

VT variants need `libvelocity_gpu_stage8_solve_nh_integration_release.<VT_PREC>.so`
under `${VT_DIR}` (default `/capstor/scratch/cscs/pmazumde/sc2026-ad-test/icon-vt-dace`).

### 8.1–8.5 Variants (what sbatch_all_sc2026.sh dispatches)

| # | Purpose                                            | Underlying command (fixed `ATM_TIMESTEP=8`)                                    |
|---|----------------------------------------------------|--------------------------------------------------------------------------------|
| 8.1 | Vanilla FP64, ss5 (reference)                    | `USE_VT_GPU=0 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 fp64 <GRID>` |
| 8.2 | Vanilla FP64, ss10 (temporally-refined reference)| `USE_VT_GPU=0 NDYN_SUBSTEPS_OVERRIDE=10 ./run/sbatch_sc2026.sh 8 fp64 <GRID>` |
| 8.3 | VT FP64                                          | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 fp64 <GRID>` |
| 8.4 | VT FP32                                          | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 fp32 <GRID>` |
| 8.5 | VT FP16                                          | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 fp16 <GRID>` |

All set `SERDE_GEN_END=51`, so `.data` dumps cover `physics_generation 0..50`.

### vt_serde knobs (for single-variant re-runs via `sbatch_sc2026.sh`)

| Env var                   | Default  | Effect |
|---------------------------|----------|--------|
| `USE_VT_GPU`              | `.true.` | false → ICON's built-in `velocity_tendencies`; true → VT's `libvelocity.so` |
| `NDYN_SUBSTEPS_OVERRIDE`  | `0`      | 0 = leave namelist (5) alone; `N>0` = force `N` |
| `SERDE_GEN_START`         | `0`      | first `physics_generation` where `do_serialize` flips on |
| `SERDE_GEN_END`           | `51`     | `physics_generation` where `do_serialize` flips off (0 = no dumps, pure timing) |

## 9. Artifacts

This HOWTO is a data-generation pipeline — analysis lives back in VT's
`SC2026_HOWTO.md §5`. At the end you should have:

### Experiment dirs

One per variant per grid, under `build/verification/experiments/`:

```
experiments/sc2026_dt8_ss5_vanilla_<grid>/     # 8.1 (FP64 reference)
experiments/sc2026_dt8_ss10_vanilla_<grid>/    # 8.2 (FP64 temporally-refined)
experiments/sc2026_dt8_ss5_gpufp64_<grid>/     # 8.3 (VT FP64)
experiments/sc2026_dt8_ss5_gpufp32_<grid>/     # 8.4 (VT FP32)
experiments/sc2026_dt8_ss5_gpufp16_<grid>/     # 8.5 (VT FP16)
```

`EXPNAME = sc2026_dt<dt>_ss<substeps>_{vanilla|gpu<prec>}_<grid>`
(`ss` = effective dyn substeps; `NDYN_SUBSTEPS_OVERRIDE=0` → namelist default 5).

### Serialized `.data` files

Inside each experiment dir, one file per serialized field × physics_generation:

```
<field>.t0.p<phys>.d<dycore>.vt<vt>.ss<substep>.data
```

Example: `p_prog.t0.p0.d1.vt1.ss5.data`. Fields of interest for the SNR
table: `p_prog` (contains `vn`, `w`) and `p_diag` (contains `vt`,
`vn_ie`, `w_concorr_c`). Default `SERDE_GEN_END=51` → `physics_generation
0..50` covered (≈50 dumps per field per run).

### Logs

`run/LOG.SAVEME-SER-<PREC>.<EXPNAME>.<jobid>.o` — stdout + stderr
combined. Look for `D2H <n>` and `Starting velocity_tendencies for
vt_generation <n>` messages to confirm the serialize path fired.

### Next

Head back to VT's `SC2026_HOWTO.md §5` to run the SNR comparisons on
these artifacts.
