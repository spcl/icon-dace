# SC2026 — ICON integration

This HOWTO is the **data-generation half** of the SC2026 reproduction.
Its output = serialized `.data` files for all 6 variants × 4 grids. For
analysis (SNR comparisons, table population) return to VT's
`SC2026_HOWTO.daint.md §3`.

## 1. Get the code
```bash
git clone --branch prat/sc2026 https://github.com/spcl/icon-dace.git
cd icon-dace
git submodule update --init --recursive
```

Copy the VT-generated wrapper (see velocity `SC2026_HOWTO.daint.md` §1).
Assumes the VT tree at `../icon-vt-dace/`:
```bash
cp ../icon-vt-dace/wrapper.f90 src/atm_dyn_iconam/wrapper.f90
```
`wrapper.f90` = `velocity_tendencies_gpu` (ICON→VT dispatch). ABI must
match the `.so`. All four precisions produce an identical wrapper, so one
copy covers fp64/fp32/fp16/bf16 — regenerate only on VT build-flag changes
that alter the ABI. (`bf16` swaps the type behind the `dace::float16` typedef,
which leaves the flattened argument list and every symbol name unchanged.)

`src/atm_dyn_iconam/serde.f90` (`module vt_serde`) is tracked here and is
canonical. **Do not overwrite it from VT** — VT's copy lags.
<!-- TODO: streamline serde.f90 into a fully auto-generated module. -->

## 2. Prereqs
- One-time: `uenv image pull icon/25.2:v1@santis`
- Every session: `uenv start --view=default icon/25.2:v1@santis`
- `spack` activated in the shell (`spack --version` should work).
  Packages load ad-hoc during build.
- `sqlite` and `zstd` installed in that spack (`spack install sqlite zstd`).
  The VT `.so` links against both, and §8 needs them reachable from the
  submitting shell.

## 3. Build directory
```bash
mkdir -p build/verification
cd build/verification
```

## 4. Configure

ICON links against the VT shared library and records its directory as an
RPATH, so that directory must hold a `libvelocity.so` at both link time and
every launch. Point `VT_DIR` at your VT tree and make the symlink:

```bash
export VT_DIR=/abs/path/to/icon-vt-dace          # dir holding the VT .so files
ln -sf "$VT_DIR"/libvelocity_gpu_stage8_solve_nh_integration_release.fp64.so \
       "$VT_DIR"/libvelocity.so

../../config/cscs/daint_yakup_dace.gpu.gh200.nvidia
```

Which precision the symlink points at does not decide what runs — §8 selects
that per job via `LD_PRELOAD`. The link only needs the symbols, which are
identical across precisions.

`VT_DIR` unset falls back to the path baked into the config script, which will
not exist outside the original author's account.

## 5. Build
```bash
make -j
```

## 6. Grid file

`GRID=<gridID>_<refinement>` maps to `icon_grid_<GRID>_G.nc`:

| `GRID`        | File                                   |
|---------------|----------------------------------------|
| `0050_R02B03` | `icon_grid_0050_R02B03_G.nc`           |
| `0010_R02B04` | `icon_grid_0010_R02B04_G.nc` (default) |
| `0008_R02B05` | `icon_grid_0008_R02B05_G.nc`           |
| `0002_R02B06` | `icon_grid_0002_R02B06_G.nc`           |

The public DWD mirror at
`http://icon-downloads.mpimet.mpg.de/grids/public/edzw/` only carries a
subset (missing at least R02B06 at the time of writing). On Daint, all
four live under
`/capstor/store/cscs/userlab/cws01/pool/data/ICON/input/icon/public/grids/`
— copy what you need into `build/verification/`:

```bash
GRIDS=/capstor/store/cscs/userlab/cws01/pool/data/ICON/input/icon/public/grids
shopt -s globstar    # enable '**' recursion
cp "$GRIDS"/**/icon_grid_0050_R02B03_G.nc ./
cp "$GRIDS"/**/icon_grid_0010_R02B04_G.nc ./
cp "$GRIDS"/**/icon_grid_0008_R02B05_G.nc ./
cp "$GRIDS"/**/icon_grid_0002_R02B06_G.nc ./
```

Off-cluster, pull from DWD where available:
```bash
curl -LO http://icon-downloads.mpimet.mpg.de/grids/public/edzw/icon_grid_0010_R02B04_G.nc
```

## 7. Runscript
`run/exp.sc2026` is the template; `make_sc2026_runscript.sh` generates
+ patches the `.run`:
```bash
./run/make_sc2026_runscript.sh
# → run/exp.sc2026.run
```

## 8. Submit

`sbatch_all_sc2026.sh <GRID>` fires all 6 variants (§8.1–8.6) for one grid:

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
under `${VT_DIR}` (default `/capstor/scratch/cscs/pmazumde/sc2026-ad-test/icon-vt-dace`),
so export `VT_DIR` here too unless your VT tree sits at that path.

Submit from a shell that has `sqlite` and `zstd` on `LD_LIBRARY_PATH`:

```bash
spack load sqlite zstd
export LD_LIBRARY_PATH="$(spack location -i sqlite)/lib:$(spack location -i zstd)/lib:$LD_LIBRARY_PATH"
ldd "$VT_DIR"/libvelocity_gpu_stage8_solve_nh_integration_release.fp64.so | grep -c 'not found'   # must print 0
```

`sbatch_sc2026.sh` passes the VT `.so` to the job as `LD_PRELOAD` via
`--export=ALL`, which also puts it in front of SLURM's task prolog. If the
prolog's shell cannot resolve the `.so`'s own dependencies, it exits non-zero
and the job dies within seconds with `TaskProlog failed status=127` — before
ICON starts, and with nothing useful in the log. `spack load` alone does not
populate `LD_LIBRARY_PATH` on this spack configuration, hence the explicit
export. The `ldd` check above is the cheap way to confirm before submitting.

Newer VT builds bake their dependency paths in as RUNPATH and resolve without
this; the check costs nothing either way.

### 8.1–8.6 Variants (what sbatch_all_sc2026.sh dispatches)

| # | Purpose                                            | Underlying command (fixed `ATM_TIMESTEP=8`)                                    |
|---|----------------------------------------------------|--------------------------------------------------------------------------------|
| 8.1 | Vanilla FP64, ss5 (reference)                    | `USE_VT_GPU=0 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 fp64 <GRID>` |
| 8.2 | Vanilla FP64, ss10 (temporally-refined reference)| `USE_VT_GPU=0 NDYN_SUBSTEPS_OVERRIDE=10 ./run/sbatch_sc2026.sh 8 fp64 <GRID>` |
| 8.3 | VT FP64                                          | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 fp64 <GRID>` |
| 8.4 | VT FP32                                          | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 fp32 <GRID>` |
| 8.5 | VT FP16                                          | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 fp16 <GRID>` |
| 8.6 | VT BF16                                          | `USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  ./run/sbatch_sc2026.sh 8 bf16 <GRID>` |

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
`SC2026_HOWTO.daint.md §3`. At the end you should have:

### Experiment dirs

One per variant per grid, under `build/verification/experiments/`:

```
experiments/sc2026_dt8_ss5_vanilla_<grid>/     # 8.1 (FP64 reference)
experiments/sc2026_dt8_ss10_vanilla_<grid>/    # 8.2 (FP64 temporally-refined)
experiments/sc2026_dt8_ss5_gpufp64_<grid>/     # 8.3 (VT FP64)
experiments/sc2026_dt8_ss5_gpufp32_<grid>/     # 8.4 (VT FP32)
experiments/sc2026_dt8_ss5_gpufp16_<grid>/     # 8.5 (VT FP16)
experiments/sc2026_dt8_ss5_gpubf16_<grid>/     # 8.6 (VT BF16)
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

Head back to VT's `SC2026_HOWTO.daint.md §3` to run the SNR comparisons on
these artifacts.
