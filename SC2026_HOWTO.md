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

The VT build bakes its `sqlite`, `zlib` and `libzstd` directories into the
`.so` as RUNPATH, so the job resolves them without `LD_LIBRARY_PATH`. Confirm
with:

```bash
ldd "$VT_DIR"/libvelocity_gpu_stage8_solve_nh_integration_release.fp64.so | grep -c 'not found'   # 0
```

A non-zero count means the `.so` was linked without the `vt-gpu` env active.
`sbatch_sc2026.sh` passes it to the job as `LD_PRELOAD` via `--export=ALL`,
which also puts it in front of SLURM's task prolog: an unresolvable dependency
there kills the job within seconds with `TaskProlog failed status=127`, before
ICON starts and with nothing useful in the log.

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
| `SERDE_GEN_LIST`          | (unset)  | explicit comma/space-separated generations to dump (max 64); overrides the window/stride gate — drives the logspaced long run below |

### Long run for the SNR-evolution figure (Figure 6(a))

Figure 6(a) is not the 50-generation window above: it is a single **long
R02B04 run** with `vn` sampled at logspaced physics generations out to 25 000.
Extend the run with `end_date` and gate the dumps with `SERDE_GEN_LIST`:

```bash
export end_date=2000-01-03T07:33:20Z          # 25000 x 8 s of model time
export SERDE_GEN_START=0 SERDE_GEN_END=25001
export SERDE_GEN_LIST=1,2,3,5,8,14,24,42,71,121,206,352,599,1021,1740,2965,5053,8610,14671,25000

for SS in 5 10; do
  for P in fp64 fp32 fp16; do
    NDYN_SUBSTEPS_OVERRIDE=$SS USE_VT_GPU=1 \
      ./run/sbatch_sc2026.sh 8 "$P" 0010_R02B04
  done
done
```

`end_date - start_date = 200000 s = 25000 × 8 s`, so the run reaches
`physics_generation 25000`. The 20-point list is `round(10**(i/19 · log10 25000))`
for `i = 0..19`. Comparison (VT `SC2026_HOWTO.daint.md §3`) writes `snr_25k.db`,
which `plotbooks/snr_evolution.ipynb` reads (as `vt_snr_25k.db`) to render
Figure 6(a). BF16 is omitted — this run predates it.

## Running on ault (A100)

The daint flow above assumes the `icon/25.2` uenv. ault (CSCS A100 node) has no
uenv; the toolchain comes from a spack environment instead. Only the steps below
differ — §1, §3, §9 are identical.

**§2 prereqs — spack `icon-gpu` env (no uenv).** ICON needs the NVHPC Fortran
compiler and a matching `netcdf-fortran`; build them once from the committed
env spec (`icon-vt-dace/velocity/arch/cscs/ault/spack-icon.yaml`):

```bash
export SPACK_TREE=$SCRATCH/spack-tree
source $SPACK_TREE/spack/share/spack/setup-env.sh
spack env create icon-gpu $SCRATCH/icon-vt-dace/velocity/arch/cscs/ault/spack-icon.yaml
spack -e icon-gpu install
spack env activate icon-gpu     # puts mpif90/mpicc, netcdf, cmake on PATH
```

`icon-gpu` builds `nvhpc@26.1` then reuses it as the fortran compiler for
`netcdf-fortran`, so the env sets `concretizer:unify:when_possible` — a single
unified solve rejects a from-source fortran compiler and fails to concretize.
Use the pinned spack from VT `SC2026_HOWTO.ault.md` (*Spack, one-time*); an
unpinned clone may not carry `nvhpc@26.1`.

**§4 configure — ault config script (sm_80).** Same `VT_DIR` symlink as daint,
then:

```bash
../../config/cscs/ault_ben_dace.gpu.a100.nvidia
```

This config passes `--disable-rpaths`, so `netcdf` and the NVHPC runtime are
*not* baked into `bin/icon` — they must be on `LD_LIBRARY_PATH` at launch (§8).

**§7 runscript — the ault maker.** `make_runscripts` has no ault machine target,
so it falls back to the `default` (CPU) target: no SLURM header, cache-blocking
`nproma`, multiple MPI ranks. `make_sc2026_runscript.ault.sh` rewrites those to
the single-GPU values the daint_gpu target uses (`nproma=0`, `nblocks_c=1`, one
rank, direct binary launch) and injects the ault SLURM header:

```bash
./run/make_sc2026_runscript.ault.sh
```

**§8 submit — three env exports.** ault has no uenv to provide the runtime, and
the grid path is taken from `$grids_folder`. Export all three before submitting;
`--export=ALL` carries them into the job:

```bash
export LD_LIBRARY_PATH=$SPACK_TREE/spack/var/spack/environments/icon-gpu/.spack-env/view/lib:$LD_LIBRARY_PATH
export VT_DIR=$SCRATCH/icon-vt-dace/velocity        # holds the VT .so files
export grids_folder=$SCRATCH/icon-grids             # holds icon_grid_<GRID>_G.nc

./run/sbatch_all_sc2026.sh 0050_R02B03
```

The variants, knobs and outputs (§8.1–8.6, §9) are identical to daint.

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
