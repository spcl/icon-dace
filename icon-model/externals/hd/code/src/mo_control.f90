! mo_control.f90 - Control variables for model housekeeping
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_control

  ! Control variables for model housekeeping.
  !
  ! U. Schlese, DKRZ, December 1994
  ! A. Rhodin, MPI, January 1999,
  !      Subroutine m_control renamed to alloc_mods and moved from
  !      module mo_control to module m_alloc_mods.
  ! L. Kornblueh, MPI, June 1999,
  !      added nproca and nprocb for driving the parallel decomposition
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! I. Kirchner, MPI, December 2000, time control
  ! I. Kirchner, MPI, March 2001, revision
  ! M. Esch, MPI, September 2002, add switch for mixed layer ocean
  ! M. Esch, MPI, November  2003, add switch for scenario runs
  ! M. Puetz, IBM, July 2009, add switches for UNITRANS support
  ! L. Kornblueh, MPI, April 2010, more switches for parallel I/O
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  ! ------------------------------------------------------------------

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PUBLIC

  SAVE    ! to assure that the scope of all variables is for the whole program run

  REAL(dp), POINTER :: vct(:)=>NULL() ! vertical coefficients table.

  INTEGER :: nproma   = -1  !   working dimension for grid-point computations
!  INTEGER :: nproca   = 1   !   number of processors in set A
!  INTEGER :: nprocb   = 1   !   number of processors in set B
! Ha Ho-Hagemann {
  INTEGER :: nproca    !   number of processors in set A
  INTEGER :: nprocb    !   number of processors in set B
! Ha Ho-Hagemann }
  INTEGER :: nprocio  = 0   !   number of processors for I/O server
  LOGICAL :: lunitrans = .FALSE.          ! if .true. use UNITRANS modules for transposes
  LOGICAL :: lunitrans_debug = .FALSE.    ! if .true. use UNITRANS debug output is enabled
  LOGICAL :: lunitrans_datatypes = .TRUE. ! if .true. use UNITRANS will use MPI datatypes

  INTEGER :: nm             !   max zonal wave number.
  INTEGER :: nn             !   max meridional wave number for m=0.
  INTEGER :: nk             !   max meridional wave number.
  INTEGER :: ngl            !   number of gaussian latitudes.
  INTEGER :: nlon           !   max number of points on each latitude line.
  INTEGER :: nlev           !   number of vertical levels.
  INTEGER :: nmp1           !   max zonal wave number + 1.
  INTEGER :: nnp1           !   max meridional wave number + 1.
  INTEGER :: nkp1
  INTEGER :: n2mp1          !   2 * (max zonal wave number + 1).
  INTEGER :: n4mp1          !   4 * (max zonal wave number + 1).
  INTEGER :: nlp2           !   max number of points per latitude line + 2.
  INTEGER :: nlevp1         !   *nlev+1.
  INTEGER :: nsp            !   number of spectral coefficients.
  INTEGER :: n2sp           !   2*number of spectral coefficients.
  INTEGER :: nhgl           !   (number of gaussian latitudes)/2.
  INTEGER :: nscan          !   current scan number.

  INTEGER :: nspace1        !   memory manager space for use of root task
  INTEGER :: nspace2        !   memory manager space for use of subtasks

  INTEGER :: nspadd         !   memory manager space increase

  INTEGER :: nvclev         !   number of levels with vertical coefficients.

  LOGICAL :: ldebug    = .FALSE. !   .true. for mass fixer diagnostics
  LOGICAL :: ldailysst = .FALSE. !   .true. for using daily SST and SIC
  LOGICAL :: lamip     = .FALSE. !   .true. for using variable sst
  LOGICAL :: ldiagamip = .FALSE. !   .true. for AMIP diagnostics
  LOGICAL :: lcouple   = .FALSE. !   .true. for a coupled run
  LOGICAL :: lcouple_co2 = .FALSE. !   .true. if coupled run includes co2
  LOGICAL :: lipcc     = .FALSE. !   .true. for run using IPCC parameters
  LOGICAL :: lnwp      = .FALSE. !   .false. for climate mode .true. for NWP mode
  LOGICAL :: lnudge    = .FALSE. !   .true. for Nudging mode
  LOGICAL :: lmidatm   = .TRUE.  !   .true. for middle atmosphere model version
  LOGICAL :: lmlo      = .FALSE. !   .true. for mixed layer ocean
  LOGICAL :: lmeltpond = .TRUE.  !   .true. for calculation of meltponds

  LOGICAL :: lprint_m0 = .FALSE. !   .false. for printing t(m0) and timestep time in stepon 
  LOGICAL :: lnmi      = .FALSE. !   .true. normal mode initialisation
  LOGICAL :: ltdiag    = .FALSE. !   .true. run with additional tendency diagnostics
  LOGICAL :: lcolumn   = .FALSE. !   .true. column model mode
  LOGICAL :: lvctch    = .FALSE. !   .true. if column model has changed vct

  LOGICAL :: l_volc    = .FALSE. !   switch for volcanic forcing

  REAL(dp):: satoverpasstime = 10.5_dp ! satellite overpass time at 10.30 a.m. [s]

  LOGICAL :: lhd      = .FALSE.  !   .true. for hydrologic discharge model

  ! Spectral and grid point initial files.
  INTEGER, PARAMETER :: nisp  = 23   ! unit for initial spectral fields.
  INTEGER, PARAMETER :: nigp  = 24   ! unit for initial grid point fields.

  ! Climate sea surface temperature and sea ice annual cycle file
  INTEGER, PARAMETER :: nist  = 20   ! unit for surf.temp. file
  INTEGER, PARAMETER :: nice  = 96   ! unit for amip ice file

  ! Climate flux correction file
  INTEGER, PARAMETER :: nflu  = 42   ! unit for flux correction file

  ! Climate leaf area index and vegetation ratio annual cycle file
  INTEGER, PARAMETER :: nvltcl   = 90 ! unit for climate leaf area index
  INTEGER, PARAMETER :: nvgratcl = 91 ! unit for climate vegetation ratio

  ! Climate land surface temperature annual cycle file
  INTEGER, PARAMETER :: ntslcl   = 92 ! unit for climate land surface temperature

  INTEGER, PARAMETER :: nsurface = 79 ! unit for boundary atm-surface jsbach

  ! optional files
  INTEGER :: ndiahdf = 10 ! I/O unit for hdiff diagnostics.

  ! Switches
  LOGICAL :: ltimer            = .FALSE. ! to use timer
  LOGICAL :: ldebugio          = .TRUE.  ! to debug IO
  LOGICAL :: ldebugmem         = .FALSE. ! to debug memory
  LOGICAL :: ltctest           = .FALSE. ! to test time control
  LOGICAL :: ldebugs           = .FALSE. ! to enable debug stream
  LOGICAL :: lroot_io          = .TRUE.  ! to disable classical root I/O
  LOGICAL :: lindependent_read = .FALSE. ! to read initial/restart input on each MPI rank separatly 
  LOGICAL :: lcollective_write = .FALSE. ! to write restart in parallel on each MPI rank

  ! output redirection
  INTEGER :: stdout_redir = 0
  INTEGER :: stderr_redir = 0

END MODULE mo_control
