# Changelog for JSBACH4+

- [Changelog for JSBACH4+](#changelog-for-jsbach4)
  - [jsbach-4.20p9](#jsbach-420p9)
  - [jsbach-4.20p8](#jsbach-420p8)
  - [jsbach-4.20p7](#jsbach-420p7)
  - [jsbach-4.20p6](#jsbach-420p6)
  - [jsbach-4.20p5](#jsbach-420p5)
  - [jsbach-4.20p4](#jsbach-420p4)
  - [jsbach-4.20p3](#jsbach-420p3)
  - [jsbach-4.20p2](#jsbach-420p2)
  - [jsbach-4.20p1](#jsbach-420p1)
  - [jsbach-4.20](#jsbach-420)

## jsbach-4.20p9

_Release Date: 5 March 2019_

- Bug fixes in `mo_jsb_tile` (missing TRIM around strings)
- Changes in `mo_jsb_grid_class`
  - Bug fix in `new_vgrid` for lbounds in one of the cases
  - Don't allow the case in `new_vgrid` that levels and dz are given but neither lbounds nor ubounds
  - Remove unnecessary TRIM's

## jsbach-4.20p8

_Release Date: 18 February 2019_

- Bug fixes for lakes
  - `fract_lice` was not propagated from the lake tile to the box tile and therefore returned as zero to the atmosphere
  - Set lake fractions on fractional land/sea grid boxes (coast) to zero

    Needs a corresponding fix in the atmosphere in `mo_echam_phy_init`
- Changes for carbon, mostly in preparation for forest management
  - Added dsl macros in preparation for forest management modules
  - Task order changes required for forest management (leads to small changes in results)
  - Removed per-canopy variables from aggregate
  - Commented variables that are not calculated
  - Adapted calculation of isteps for ECHAM in mo_jsb_time_iface.f90

    ... to avoid integer overflow for long simulation times, i.e. large difference between start and end date (>70y)
  - Subroutines for cpool sums and carbon conservation test

    C conservation test for single operations on the tile (still todo: C conservation test on the box for timestep)
  - Minor corrections in driver and forcing module
  - Small change of wind damage interface required for cohorts + minor corrections
  - Set max_C_content_woods as in jsbach3 at svn revision 9089

    (only) changes c wood results when close to previous max
  - Implemented allometric dependance of MaxLai on biomass
    Ported from Kim Naudt's work in JSBACH3.
    - new switch l_forestRegrowth
    - does not change results when switched off!
    - if switched on, the MaxLai of pfts with ForestFlag is derived from the available biomass
  - Bug fix in mo_pheno_init.f90
- Added daily carbon conservation test on box tile
  - Carbon conservation test for each new day in `jsbach_start_timestep`
  - Added calls to `jsbach_start_timestep` and `jsbach_finish_timestep` to driver
  - Corrected unit of fire atm flux `cflux_fire_all_2_atm`
- Changes to output behaviour
  - Append domain (model) ID to stream names only if more than one domain (model)
  - Default output behaviour for ECHAM is now false and a variable is only written out if `loutput=.TRUE.` is explicitely set for that variable (either in the code or via namelist)
- Cleanup to get rid of compiler warnings and other formatting changes
- Various (bug) fixes
  - Bug fix in `mo_hsm_class.f90`
  - Correct inconsistent definition of tasks in `mo_sse_interface.f90`
  - Added `l_aggregate_all=.TRUE.` where necessary to the `jsb_add_var` calls
  - Corrected box value of gross assimilation (`* fract_fpc_max`) and some minor changes in the ASSIMI_ memory output descriptions
  - Pre-populate some variables to avoid problem with WHERE statements with Cray compiler
- Move adapters to ICON/ECHAM

  Needs new adapter file in ICON resp. ECHAM
  - Remove directory `src/adapters`
  - Necessary code adaptations

- Various changes in pre-processing script
  - Make sure that fractional coastal points don't contain lakes
  - Make sure that glacier fraction is either 0 or 1 and adapt land sea mask to be consistent with that (not for coupled)
  - Suppress small lake fractions <0.1 from input
  - Derive vegetation fractions from LUHv2 land cover data
  - Remove partial ocean points around Caspian Sea in FR_LAND_TOPO from extpar file
  - Remove option `old_r2b4_amip`
  - Update output_path

## jsbach-4.20p7

_Release Date: 26 June 2018_

- Re-insert `l_compat401` switch for backward compatibility
- Use 10m wind from atmosphere in the interface (instead of lowest level wind)

  Note: currently, the 10m wind passed from the atmosphere should still be 0.8 * lowest level wind until issue with lookup table overflow is resolved in icon atmosphere
- Various (technical) code fixes
  - Fix pointless TRIM of already minimal-size module names
  - Fix missing communicator argument in p_bcast calls
  - Fix pointless polymorphism in `mo_jsb_tile_class.f90`
  - Fix implicit (re-)allocation on assignment
  
    For performance reasons it is better not to use automatic reallocation of left-hand side (e.g. realloc-lhs compiler option for intel)
  - Use correct character length for version string

## jsbach-4.20p6

_Release Date: 20 June 2018_

- Earlier reading of lctlib so that it can be used in init_usecase to set up tile structure and processes
- Use correct communicator in `mo_jsb_io_netcdf_iface.f90`
- Replaced `c_litter_green_ag` in the albedo calculations with `c_ag_sum_1`

  This is to be consistent with `c_bg_sum` used in the albedo calculations.

## jsbach-4.20p5

_Release Date: 20 June 2018_
(replaced by jsbach-4.20p6 in order to get a few fixes in)

- Timers
  - Implement constructor function for `t_jsb_process_task` abstract type and add timers for each process task
  - Namelist switches are now `debug_level` and `timer_level`, both can take values between 0 and 3 and are independent of each other
  - Move main JSBACH timer from ICON/ECHAM to JSBACH
- Use proper MPI communicator in `p_bcast` calls (created problem on CRAY with asynchronous output)
- Carbon module
  - Add the summation of `faPAR` from `faPAR_cl`
  - Changed calculation of fuel from canopy to potential vegetation area and adapted fire_litter_threshold
  - Correction of fire bug due to fuel and typo correction of subroutine name
- Some re-factoring in `update_radiation_par`
- Deletion of `l_compat401` switch, changes for accompanying tuning of model and pre-processing script for input data
  - Delete `l_compat401` since standard experiments will now use new snow model and freezing of soil water
  - Change canopy albedo in lctlib
  - Set default for switch `use_alb_canopy` to false in order to actually use canopy albedo values from lctlib (instead of from ini file)
  - Change albedo parameters for glaciers in `mo_rad_constants.f90`
  - Use remapdis for ICON grids above R2B5 instead of remapycon in pre-processing script `jsbach4_ini_files_from_gauss.sh`
  - Add function to modify soil and root depth as part of tuning the model in pre-processing script

    Soil water holding capacity is modified by changing root depth and soil depth.
    They are decreased by 30% in the range 30S - 30N to warm the tropical ontinents
    and increased by 0.2m for the other land grid points to cool the extra-tropical
    continents.
  - Use revision for input data; put data for revision r0002 into new directory in pool
- Further smaller fixes and code cleanup


## jsbach-4.20p4

_Release Date: 17 April 2018_

- Reading of carbon pools from ini file
  - New function `jsbach_init_after_restart` which now reads the carbon pools from the ini file **after** the restart files have been read, but **before** the output is initialized. In ICON, the output of the first step now contains the correct pools that have just been read in. Note that the run script has to take care of re-setting `read_cpools` to `.FALSE.` for the jobs following the first one!
  - Some minor changes (comment for `read_cpools` config parameter and changes in script `cpool_file_from_restart_file.sh` for generation of cpool file from restart file)
- Re-structuring of aggregator

  Leads to a few percent of speed gain in aggregation
- Change in dsl4jsb.py pre-preprocessing script

  If the content of the source file hasn't changed compared to the target file, reset the modification time of the source file to the one of the target file, so that make doesn't think the source file is newer and therefore tries to process the file (but `dsl4jsb.py` doesn't create a new target since the content hasn't changed ... this can lead to an infinite loop in make)
- Some changes to `jsbach4_ini_files_from_gauss.sh` script for generation of input files
  - Set URL of svn repo to new address `svn.mpimet.mpg.de`
  - Allow generation of files for coupled setup for non-r2b4 resolutions (might still need some work)

## jsbach-4.20p3

_Release Date: 10 April 2018_

- Give back new saturated specific humidity over lakes (water and ice) to the atmosphere

  With the corresponding fix in the ICON atmosphere which puts `qsat_lwtr` and `qsat_lice` into the `bb_btm` array, this closes a water leakage bug.
- Some bug fixes in the carbon module
- Some other bug fixes
- Some changes to `jsbach4_ini_files_from_gauss.sh` script for generation of input files
  - updated module versions (cdo etc)
  - output files for debugging are in T255 (icon)
  - updated files for coupled setup
  - changed default cdo remap scheme from "con" to "ycon"

## jsbach-4.20p2

_Release Date: 29 March 2018_

- Bug fixes and changes for carbon/fuel/disturbance modules
- Reading of carbon pools from external file at start or restart
- Variable name changes
- Other bug fixes
- Implement (again) setup for terraplanet simulations
- `dsl4jsb.py`: only process source file if it has changed

  If a source file is processed, the md5 hash value of the file is now written into the header of the target file. The next time `dsl4jsb.py` is called on the source file, a new target file is only generated if the hash value of the source file and the hash value stored in the target file are different.
- Re-structure creation and handling of vertical axes in order to work with new features in ICON infrastructure
  - Works now with asynchronous I/O.
  - In preparation for use with CDI/PIO.
  - Initialization of JSBACH had to be separated into several subroutines because these have to be called at different stages in ICON.

## jsbach-4.20p1

_Release Date: 13 February 2018_

## jsbach-4.20

_Release Date: 16 December 2017_
