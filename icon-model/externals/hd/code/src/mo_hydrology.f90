! mo_hydrology.f90 - main model routines of HD model
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann and Ha Ho-Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_hydrology

  !
  ! Authors:
  !
  ! S. Hagemann, MPI, October 1999, original source
  ! U. Schlese , MPI, January 2001, cleanup and introduction
  ! L. Kornblueh, MPI, April 2002, cleanup, parallelization,
  !                                and packed in one module
  ! K. Ketelsen, NEC, February 2003, optimization
  ! T. Jahns, DKRZ, November 2009, Speedup
  ! V. Gayler, MPI, Mai 2011, cleanup, correction for lakes,
  !                           new interpolation routines,
  !                           assure water conservation
  !
  ! Note that the coupled 0.5 degree version is also part of MPI-ESM, the Earth System Model of the 
  !    Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !    Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !    version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !    doi: 10.1029/2018MS001400.
  !
  ! H. Ho-Hagemann, HZG, June 2014, introducing OASIS coupling in offline version.
  ! S. Hagemann,MPI, Feb. 2015, further refinements of offline version.
  ! S. Hagemann, HZG, 2018, Refinements for HD 5 Min. version. 
  !

  USE mo_kind,          ONLY: dp
  USE mo_constants,     ONLY: api, a, rhoh2o
  USE mo_control,       ONLY: nlon, ngl, lcouple
  USE mo_exception,     ONLY: message, message_text, finish
  USE mo_grid,          ONLY: domain, read_grid_info
  USE mo_hd_highres_io, ONLY: hd_highres_open, hd_highres_write, hd_highres_close
  USE mo_gaussgrid,     ONLY: gridarea, philat, philon
!          Note that gridarea is used for nremap = 1 (assuming gaussian grid forcing), 
!          water budget calculations for conservations tests and related messsages, 
!          back trafo of discharge to ocean model grid (to calculate gl_disch) 
  USE mo_io,            ONLY: io_open, io_close, io_read, io_write
  USE mo_netcdf,        ONLY: nf_global, nf_double, file_info, io_enddef, &
                              io_inq_dimid, io_inq_varid, io_inq_dimlen, &
                              io_def_dim, io_def_var, io_put_att_text, &
                              io_put_att_int, io_get_var_double, io_put_var_double, &
                              nf_get_att_int, nf_noerr, nf_max_name
  USE mo_memory_g3b,    ONLY: aros, adrain, disch, slm, alake, awfre, &
                              glac, apmecal
  USE mo_mpi,           ONLY: p_io, p_parallel_io, p_parallel, p_bcast
  USE mo_time_control,  ONLY: lstart, get_time_step,    &
                              delta_time, ev_puthd, get_interval_seconds, &
                              io_time_event,  &
                              initial_date, start_date, out_convert_date, inp_convert_date, &
                              current_date
  USE mo_array_utils,   ONLY: dec_monotonic_closest_midpoint, &
                              inc_monotonic_closest_midpoint
  USE mo_coupling,      ONLY: set_grid_dimensions, set_local_partition, &
                              is_coupled_run, fdir_hd, runoff_s, runoff_dr, &
                              lcoupling_atm, lcoupling_oce

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_hydrology, cleanup_hydrology, read_hydrology
  PUBLIC :: hydrology_model, hydrology_restart, redistribute_sinks
  PUBLIC :: hd_open_timeseries, hd_close_timeseries

  INTEGER                ::  highres_file_id   ! netcdf file id
  CHARACTER(nf_max_name) ::  filename          ! file name
  CHARACTER(nf_max_name) ::  varname           ! variable name

  TYPE cart_idx_2d
    INTEGER :: ilon, ilat
  END TYPE cart_idx_2d

  TYPE cart_coord_2d
    REAL(dp) :: longitude, latitude
    INTEGER  :: dir
  END TYPE cart_coord_2d

  TYPE cart_xidx_2d
    INTEGER :: ilon, ilat, extlen
    REAL(dp), ALLOCATABLE :: amod(:), akdiv(:)
  END TYPE cart_xidx_2d

  !
  ! Grid info of HD grid from parameter file (file_para) that replaces hd_domain.inc
  CHARACTER(LEN=80)     :: file_para = 'hdpara.nc'
  TYPE(domain), PUBLIC  :: grid_hd    

  REAL(dp), PARAMETER :: fullcirc = 360.0_dp
!  REAL(dp), PARAMETER :: hd_scal_lon = fullcirc/nl           ! only used for previous interpolation
!  REAL(dp), PARAMETER :: hd_scal_lat = 0.5_dp*fullcirc/nb    ! only used for previous interpolation 

  ! corresponding coordinates on echam grid
  REAL(dp), PARAMETER :: oclorg = 0.0_dp, ocborg = 90.0_dp

  ! HD model time steps
  INTEGER, PUBLIC :: hd_calling_interval     ! calling interval in seconds
  INTEGER         :: riverflow_timestep      ! sub time step used for riverflow
  INTEGER, PUBLIC :: hd_steps_per_day        ! number of hd_model calls per day
  INTEGER         :: riverflow_steps_per_day ! number of riverflow time steps per day
  REAL(dp)        :: div_riverflow_timestep  ! 1/riverflow_timestep

  REAL(dp), ALLOCATABLE :: alf_k(:,:)    ! retention constant k, overflow
  REAL(dp), ALLOCATABLE :: alf_n(:,:)    ! number of reservoirs n, overflow
  REAL(dp), ALLOCATABLE :: arf_k(:,:)    ! retention constant k, riverflow
  REAL(dp), ALLOCATABLE :: arf_n(:,:)    ! number of reservoirs n, riverflow
  REAL(dp), ALLOCATABLE :: agf_k(:,:)    ! retention constant k, baseflow
  REAL(dp), ALLOCATABLE :: agf_n(:,:)    ! number of reservoirs n, baseflow
  INTEGER,  ALLOCATABLE :: fdir(:,:)     ! river direction
  REAL(dp), ALLOCATABLE :: hd_lsm(:,:)   ! hd model land mask
  INTEGER,  ALLOCATABLE :: filnew(:,:)   ! Longitude index of flow destination according to fdir
  INTEGER,  ALLOCATABLE :: fibnew(:,:)   ! Latitude index of flow destination according to fdir

  REAL(dp), ALLOCATABLE :: finfl(:,:)    ! inflow data
  REAL(dp), ALLOCATABLE :: fgmem(:,:,:)  ! intermediate linear baseflow reservoir
  REAL(dp), ALLOCATABLE :: frfmem(:,:,:) ! intermediate reservoirs, inflow cascade
  REAL(dp), ALLOCATABLE :: flfmem(:,:,:) ! intermediate reservoir, linear overflow
  INTEGER, PARAMETER    :: nmemrf = 5    ! number of riverflow reservoir cascades

  REAL(dp), ALLOCATABLE, PUBLIC :: hd_area(:)          ! grid cell area [m2]
  REAL(dp), ALLOCATABLE         :: lon_hd(:,:)         ! longitudes of the HD model grid
  REAL(dp), ALLOCATABLE         :: lat_hd(:,:)         ! latitudes of the HD model grid
  REAL(dp), ALLOCATABLE         :: runoff_hd(:,:)      ! runoff into the ocean 
  REAL(dp), ALLOCATABLE, PUBLIC :: friv(:,:)           ! river flow
  REAL(dp), ALLOCATABLE, PUBLIC :: water_to_ocean(:,:) ! water going to the ocean directly
!
! *** Important arrays for routine hydrology model -> made allocatable 
  REAL(dp), ALLOCATABLE :: input_overlandflow(:,:)
  REAL(dp), ALLOCATABLE :: input_baseflow(:,:)
  REAL(dp), ALLOCATABLE :: overlandflow(:,:)
  REAL(dp), ALLOCATABLE :: baseflow(:,:)
  REAL(dp), ALLOCATABLE :: riverflow(:,:)
  REAL(dp), ALLOCATABLE :: overlandplusbaseflow(:,:)

  TYPE(cart_idx_2d),  ALLOCATABLE :: oclook_cache(:,:)   ! closest ECHAM grid ocean cell
  TYPE(cart_xidx_2d), ALLOCATABLE :: arf_n_kas(:)        !
  TYPE(cart_xidx_2d), ALLOCATABLE :: alf_n_kas(:)        !
  TYPE(cart_xidx_2d), ALLOCATABLE :: agf_n_kas(:)        !

  ! Arrays on the Gaussian grid
  REAL(dp), ALLOCATABLE, TARGET :: gl_aros(:,:)    ! runoff
  REAL(dp), ALLOCATABLE, TARGET :: gl_adrain(:,:)  ! drainage
  REAL(dp), ALLOCATABLE, TARGET :: gl_disch(:,:)   ! discharge
  REAL(dp), ALLOCATABLE, TARGET :: gl_slm(:,:)     ! land sea mask [1,0]
  REAL(dp), ALLOCATABLE, TARGET :: gl_alake(:,:)   ! fractional lake mask
  REAL(dp), ALLOCATABLE, TARGET :: gl_awfre(:,:)   ! P-E over ocean and lakes
  REAL(dp), ALLOCATABLE, TARGET :: gl_apmecal(:,:) ! P-E of glaciers
  REAL(dp), ALLOCATABLE, TARGET :: gl_glac(:,:)    ! glacier mask [1,0]
  REAL(dp), ALLOCATABLE         :: olm(:,:)        ! ocean land mask [1,0] (without lakes)

  ! Arrays for Scrip remapping
  INTEGER,  ALLOCATABLE :: src_address(:)   ! address indices of the source grid
  INTEGER,  ALLOCATABLE :: dst_address(:)   ! address indices of the destination grid
  REAL(dp), ALLOCATABLE :: remap_matrix(:,:)
  INTEGER               :: num_links        ! number of overlaps between source and destination grid
  INTEGER               :: num_weights      ! number of weights (1 for first order conservative remapping)

  ! parameters of the hydrology namelist
  LOGICAL, PUBLIC :: ldebughd                     ! true for debugging
  LOGICAL, PUBLIC :: diag_water_budget            ! true to print water budget diagnostics
  INTEGER, PUBLIC :: nhd_diag                     ! index for hd diagnostics
  LOGICAL, PUBLIC :: lhd_highres = .FALSE.        ! true for additional output on the HD grid 
  LOGICAL, PUBLIC :: lbase = .TRUE.               ! baseflow ON or OFF
  LOGICAL, PUBLIC :: locean = .TRUE.              ! closure of water budget for ocean coupling
  REAL(dp)        :: fllog1, fblog1               ! grid cell for diagnostics (with nhd_diag=99)
  REAL(dp)        :: fllog2, fblog2               ! grid cell for diagnostics (with nhd_diag=99)
  INTEGER, PUBLIC :: nremap                       ! index for type of remapping to the HD model grid
  LOGICAL, PUBLIC :: lhd_rout = .FALSE.           ! true for routing via index arrays 
  REAL(dp),PUBLIC :: fk_rfk                       ! Modification factor of k value for riverflow
  REAL(dp),PUBLIC :: fk_lfk                       ! Modification factor of k value for overland flow
  REAL(dp),PUBLIC :: fk_gfk                       ! Modification factor of k value for baseflow
  INTEGER, PUBLIC :: irf_vel                      ! index for type of discharge dependence of riverflow velocity
  REAL(dp),PUBLIC :: qrf_ref                      ! Reference discharge for discharge dependent riverflow velocity

  INTEGER, SAVE   :: isolog_unit                  ! unit for outflow diagnostics

CONTAINS

  SUBROUTINE init_hydrology (slm, alake)

    REAL(dp), INTENT(in), OPTIONAL :: slm(:,:)
    REAL(dp), INTENT(in), OPTIONAL :: alake(:,:)

    ! Read HD grid info from HD parameter file
    CALL read_grid_info(file_para, grid_hd, 'lon', 'lat', 1, 1)

    WRITE(message_text,*) 'Dimensions: nl, nb', grid_hd%nlon, grid_hd%nlat
    CALL message('init_hydrology', message_text)

    ! read HD namelist
    CALL config_hydrology

    ! set grid dimensions for the coupler
    CALL set_grid_dimensions(grid_hd%nlon, grid_hd%nlat, &
                             grid_hd%origin_lon, grid_hd%origin_lat, &
                             grid_hd%resolution)

    ! set local partition information for the coupler
    ! (MoHa: at this point I assume that HD runs on a single process,
    !        if this is not the case anymore, the partition information
    !        has to be adjusted accordingly)
    CALL set_local_partition(grid_hd%nlon, grid_hd%nlat, 1, 1)

    ! initialize HD model memory
    CALL hd_init_memory(grid_hd%nlon, grid_hd%nlat, nlon, ngl)

    IF (p_parallel_io) THEN

      CALL set_riverflow_timestep

      ! Read parameter fields and restart file for the HD Model
      CALL read_hydrology

      ! open file for 'isolog' timeseries
      IF (nhd_diag > 0) CALL hd_open_timeseries(nhd_diag)

    END IF

    IF (PRESENT(slm) .AND. PRESENT(alake)) THEN
       CALL hydrology_slm_invariants(slm, alake)
    ELSE
       CALL hydrology_slm_invariants
    END IF

       WRITE (message_text,*) ' Boing'
       CALL message('init_hydrology', message_text)

    IF (lhd_highres .AND. p_parallel_io) THEN
       filename = 'hd_higres.nc'
       varname = 'friv'
       CALL hd_highres_open(filename, varname, highres_file_id)
       WRITE (message_text,*) 'River discharge on HD model grid written to: ', TRIM(filename)
       CALL message('init_hydrology', message_text)
    END IF

  END SUBROUTINE init_hydrology

  SUBROUTINE config_hydrology

    !----------------------------------------------------------------------------------
    ! parameters of namelist HYDROLOGY_CTL
    !
    !          ldebughd    additional output for debugging  
    ! diag_water_budget    switch for additional water budget diagnostics  
    !       lhd_highres    switch for outflow diagnostic on HD model grid
    !             lbase    switch for baseflow calculations
    !            locean    closure of water budget for ocean coupling
    !          nhd_diag    region number for outflow diagnostic (formerly isolog)
    !                         0   none
    !                         1   Bothnian Bay/Sea
    !                         2   Torneaelven
    !                         4   St.Lawrence
    !                         5   Paraguay
    !                         6   Oder
    !                         7   Elbe
    !                         8   Oranje
    !                         9   Amudarya
    !                        10   Lena
    !                        99   user defined (fblog1, fllog1, fblog2, fllog2)
    !            fllog1    user defined grid cells for diagnostics (with nhd_diag=99)
    !            fblog1        fllog1, fblog1: longitude, latitude of grid cell 1
    !            fllog2        fllog2, fblog2: longitude, latitude of grid cell 2
    !            fblog2
    !            nremap    Type of Interpolation from input (atmospheric) grid to HD grid
    !                         0   Input = Output
    !                         1   using HDMAP routine by Veronika (default)
    !                         2   0.5 degree to 5 Min.
    !                         3   Input = Output + logitudinal shift by 180 degree.
    !          lhd_rout    Switch for original routing (F) or via index arrays (T)
    !            fk_rfk    Modification factor of k value for riverflow
    !            fk_lfk    Modification factor of k value for overland flow
    !            fk_gfk    Modification factor of k value for baseflow
    !            irf_vel   Index for type of discharge dependence of riverflow velocity
    !                         0    None
    !                         1    velocity ~ 4th squareroot of Q (riverflow)
    !            qrf_ref   Reference discharge for discharge dependent riverflow velocity
    !----------------------------------------------------------------------------------
    USE mo_namelist,         ONLY: open_nml, position_nml, POSITIONED
    USE mo_io_units,         ONLY: nout
 
    ! local variables

    INTEGER                :: read_status, inml, iunit

    INCLUDE 'hydrology_ctl.inc'

    ! set default values of the namelist parmeters

    ldebughd = .FALSE.          ! additional output for debugging
    diag_water_budget = .FALSE. ! prints to diagnose the water budget
    lbase = .TRUE.              ! base flow calculation swiched on
    locean = .FALSE.            ! close water budget for ocean coupling
    nhd_diag = 0
    lhd_highres = .FALSE.
    fllog1 = 0.0_dp
    fblog1 = 0.0_dp
    fllog2 = 0.0_dp
    fblog2 = 0.0_dp
    nremap = 1
!Hag
    lhd_rout = .TRUE.           ! original routing (0) or via Index arrays (1)
    fk_rfk = 1.0_dp
    fk_lfk = 1.0_dp
    fk_gfk = 1.0_dp
    irf_vel = 0
    qrf_ref = 1000.0_dp
     
    ! read namelist hydrology_ctl

    IF (p_parallel_io) THEN
       inml = open_nml('namelist.hdset')
       iunit = position_nml ('HYDROLOGY_CTL', inml, status=read_status)
       SELECT CASE (read_status)
       CASE (POSITIONED)
          READ (iunit, hydrology_ctl)
          CALL message('config_hydrology', 'Namelist HYDROLOGY_CTL: ')
          WRITE(nout, hydrology_ctl)
       END SELECT
       nhd_diag = MIN(nhd_diag, 99)
    END IF
    IF (p_parallel) THEN
       CALL p_bcast(ldebughd, p_io)
       CALL p_bcast(diag_water_budget, p_io)
       CALL p_bcast(locean, p_io)
       CALL p_bcast(lhd_highres, p_io)
       CALL p_bcast(lhd_rout, p_io)
    ENDIF

  END SUBROUTINE config_hydrology

  SUBROUTINE hd_init_memory(nl,nb, nlon,ngl) 

    INTEGER,    INTENT(in)    :: nl      ! Number of HD longitudes
    INTEGER,    INTENT(in)    :: nb      ! Number of HD latitudes
    INTEGER,    INTENT(in)    :: nlon    ! Number of forcing longitudes
    INTEGER,    INTENT(in)    :: ngl     ! Number of forcing latitudes

!!    REAL(dp), ALLOCATABLE,   INTENT(out) :: slm(:,:)           ! integer land sea mask
!!    REAL(dp), ALLOCATABLE,   INTENT(out) :: slf(:,:)           ! fractional land sea mask
!!    REAL(dp), ALLOCATABLE,   INTENT(out) :: lake(:,:)          ! fractional lake mask

    ! Initialize memory for the HD Model
    ! corresponds to offline routine 'hdini.f' by S. Hagemann

    IF (p_parallel_io) THEN

      ALLOCATE (alf_k(nl,nb))          ; alf_k(:,:)    = 0.0_dp
      ALLOCATE (alf_n(nl,nb))          ; alf_n(:,:)    = 0.0_dp
      ALLOCATE (arf_k(nl,nb))          ; arf_k(:,:)    = 0.0_dp
      ALLOCATE (arf_n(nl,nb))          ; arf_n(:,:)    = 0.0_dp
      ALLOCATE (agf_k(nl,nb))          ; agf_k(:,:)    = 0.0_dp
      ALLOCATE (agf_n(nl,nb))          ; agf_n(:,:)    = 0.0_dp
      ALLOCATE (fdir(nl,nb))           ; fdir(:,:)     = 0
      ALLOCATE (hd_lsm(nl,nb))         ; hd_lsm(:,:)   = 0.0_dp
      ALLOCATE (finfl(nl,nb))          ; finfl(:,:)    = 0.0_dp
      ALLOCATE (fgmem(nl,nb,1))        ; fgmem(:,:,:)  = 0.0_dp
      ALLOCATE (frfmem(nl,nb,nmemrf))  ; frfmem(:,:,:) = 0.0_dp
      ALLOCATE (flfmem(nl,nb,1))       ; flfmem(:,:,:) = 0.0_dp
      ALLOCATE (hd_area(nb))           ; hd_area(:)    = 0.0_dp
      ALLOCATE (oclook_cache(nl,nb))   ; oclook_cache  = cart_idx_2d(-1, -1)
      IF (is_coupled_run() .AND. lcoupling_atm) THEN
        ALLOCATE (runoff_s(nl,nb))       ; runoff_s(:,:) = 0.0_dp
        ALLOCATE (runoff_dr(nl,nb))       ; runoff_dr(:,:) = 0.0_dp
      ELSE
        IF (nremap.EQ.0 .OR. nremap.EQ.3) THEN
           ALLOCATE (runoff_s(nl,nb))     ; runoff_s(:,:) = 0.0_dp
           ALLOCATE (runoff_dr(nl,nb))    ; runoff_dr(:,:) = 0.0_dp
        ENDIF
      ENDIF
      IF (is_coupled_run()) THEN 
        ALLOCATE (fdir_hd(nl,nb))        ; fdir_hd(:,:)  = 0
      ENDIF
      IF (lhd_rout) THEN
         ALLOCATE (filnew(nl,nb))      ; filnew(:,:)   = 0
         ALLOCATE (fibnew(nl,nb))      ; fibnew(:,:)   = 0
      ENDIF
    END IF

    ALLOCATE (gl_aros(nlon,ngl))    ; gl_aros(:,:)    = 0.0_dp
    ALLOCATE (gl_adrain(nlon,ngl))  ; gl_adrain(:,:)  = 0.0_dp
    ALLOCATE (gl_disch(nlon,ngl))   ; gl_disch(:,:)   = 0.0_dp
    ALLOCATE (gl_slm(nlon,ngl))     ; gl_slm(:,:)     = 0.0_dp
    ALLOCATE (gl_alake(nlon,ngl))   ; gl_alake(:,:)   = 0.0_dp
    ALLOCATE (gl_awfre(nlon,ngl))   ; gl_awfre(:,:)   = 0.0_dp
    ALLOCATE (gl_apmecal(nlon,ngl)) ; gl_apmecal(:,:) = 0.0_dp
    ALLOCATE (gl_glac(nlon,ngl))    ; gl_glac(:,:)    = 0.0_dp

    ALLOCATE (input_overlandflow(nl, nb))
    ALLOCATE (input_baseflow(nl, nb))
    ALLOCATE (overlandflow(nl, nb))
    ALLOCATE (baseflow(nl,nb))
    ALLOCATE (riverflow(nl,nb))
    ALLOCATE (overlandplusbaseflow(nl, nb))

    IF (locean) THEN
      ALLOCATE (lon_hd(nl,nb))          ! longitudes of the HD model grid
      ALLOCATE (lat_hd(nl,nb))          ! latitudes of the HD model grid
    ENDIF
    ALLOCATE (runoff_hd(nl,nb))       ! runoff into the ocean 
    ALLOCATE (friv(nl,nb))            ! river flow
    ALLOCATE (water_to_ocean(nl,nb))  ! water going to the ocean directly


  END SUBROUTINE hd_init_memory

  SUBROUTINE read_hydrology
    !
    ! reads parameter fields and restart file for the HD model
    !
    !   hdrestart.nc: restart file with reservoir cascade arrays, ...
    !      hdpara.nc:parameter file with land sea mask, runoff directions, ...


    TYPE (FILE_INFO)  :: parafile, hdfile

    INTEGER nvarid, fileid, i, status
    INTEGER hd_steps_per_day_restart    ! hd model timestep of the run the restart file originates from
    INTEGER riverflow_timestep_restart  ! riverflow timestep of the run the restart file originates from 
    INTEGER :: yyyymmdd                 ! Date info of initial date of the whole simulation

    CHARACTER(len=80) :: file_start
    CHARACTER(len= 7) :: varname
    REAL(dp), ALLOCATABLE :: hd_dp_read(:,:)
    REAL(dp) :: factor
    LOGICAL :: lex
    INTEGER :: buf(1)

    ALLOCATE(hd_dp_read(grid_hd%nlon,grid_hd%nlat))

    ! File names

    IF (lstart) THEN
      file_start = 'hdstart.nc'
    ELSE
      file_start = 'hdrestart.nc'
    END IF

    ! Read parameter: Land sea mask, RDF, ...

    INQUIRE (file=file_para, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(file_para),'>'
      CALL message('read_hydrology', message_text)
      CALL finish ('read_hydrology', 'run terminated.')
    ENDIF

    parafile%opened = .FALSE.
    CALL IO_open (file_para, parafile, IO_READ)
    CALL message('', '')
    WRITE (message_text,*) 'Reading hdpara from file ', TRIM(file_para)
    CALL message('read_hydrology', message_text)

    fileID = parafile%file_id
  
    CALL IO_inq_varid (fileID, 'FLAG', nvarid)
    CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
    hd_lsm = REAL(NINT(hd_dp_read), dp)
    CALL IO_inq_varid (fileID, 'FDIR', nvarid)
    CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
    ! convert back to integer directions pointlessly stored as double
    fdir = MIN(NINT(hd_dp_read), 9)
    IF (is_coupled_run()) THEN
      fdir_hd = fdir
    ENDIF
    CALL IO_inq_varid (fileID, 'ALF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, alf_k)
    CALL IO_inq_varid (fileID, 'ALF_N', nvarid)
    CALL IO_get_var_double (fileID, nvarid, alf_n)
    CALL IO_inq_varid (fileID, 'ARF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, arf_k)
    CALL IO_inq_varid (fileID, 'ARF_N', nvarid)
    CALL IO_get_var_double (fileID, nvarid, arf_n)
    CALL IO_inq_varid (fileID, 'AGF_K', nvarid)
    CALL IO_get_var_double (fileID, nvarid, agf_k)
#ifdef HD_5MIN
    CALL IO_inq_varid (fileID, 'AREA', nvarid)
    CALL IO_get_var_double (fileID, nvarid, hd_area)
#endif
    IF (lhd_rout) THEN
      CALL IO_inq_varid (fileID, 'FILNEW', nvarid)
      CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
      filnew = NINT(hd_dp_read)
      CALL IO_inq_varid (fileID, 'FIBNEW', nvarid)
      CALL IO_get_var_double (fileID, nvarid, hd_dp_read)
      fibnew = NINT(hd_dp_read)
    ENDIF

    CALL IO_close(parafile)

    ! Apply modification factors for sensitivity studies
    arf_k = arf_k * fk_rfk
    alf_k = alf_k * fk_lfk
    agf_k = agf_k * fk_gfk

    ! Read restart information: Reservoirs and inflow

    hdfile%opened = .FALSE.
    INQUIRE (file=file_start, exist=lex)
    IF ( .NOT. lex ) THEN
      WRITE (message_text,*) 'Could not open file <', TRIM(file_start), '>'
      CALL message('read_hydrology', message_text)
      CALL finish ('read_hydrology', 'run terminated.')
    ENDIF

    CALL IO_open (file_start, hdfile, IO_READ)
    WRITE (message_text,*) 'Reading hdrestart from file ', TRIM(file_start)
    CALL message('read_hydrology', message_text)
    CALL message('', '')

    fileID = hdfile%file_id

    CALL IO_inq_varid (fileID, 'FLFMEM', nvarid)
    CALL IO_get_var_double (fileID, nvarid, flfmem)

    varname = 'FRFMEM'
    DO i=1, nmemrf
      WRITE(varname(7:7), '(i1)') i
      CALL IO_inq_varid (fileID, varname, nvarid)
      CALL IO_get_var_double (fileID, nvarid, frfmem(:,:,I))
    ENDDO

    CALL IO_inq_varid (fileID, 'FGMEM', nvarid)
    CALL IO_get_var_double (fileID, nvarid, fgmem)
    fgmem(:,:,1) = fgmem(:,:,1) * hd_lsm(:,:) ! third array dimension is 1 

    CALL IO_inq_varid (fileID, 'FINFL', nvarid)
    CALL IO_get_var_double (fileID, nvarid, finfl)

    ! In principle the intermediate reservoirs have the unit [m3]. However,
    ! routine kasglob was designed for daily calls, and the reservoirs are
    ! treated as volume flows with unit [m3/day]. The state of the reservoir
    ! restart variables depend on the hd model time step, the inflow on the
    ! riverflow timestep. If the timesteps are changeing within a simulation
    ! the restart variables need to be adapted.

    status = NF_GET_ATT_INT (fileID, NF_GLOBAL,'hd_steps_per_day', buf)
    hd_steps_per_day_restart = buf(1)
    IF (status /= NF_NOERR)  hd_steps_per_day_restart = 1    ! old restart files
    IF (hd_steps_per_day_restart /= hd_steps_per_day) THEN
       factor = REAL(hd_steps_per_day,dp) / REAL(hd_steps_per_day_restart,dp)
       flfmem = flfmem * factor
       frfmem = frfmem * factor
       fgmem = fgmem * factor
    END IF

    status = NF_GET_ATT_INT (fileID, NF_GLOBAL,'riverflow_timestep', buf)
    riverflow_timestep_restart = buf(1)
    IF (status /= NF_NOERR) THEN
       ! handling of old restart files without riverflow_timestep attribute
       riverflow_timestep_restart = 4
       finfl = finfl / riverflow_timestep_restart
    END IF
    IF (riverflow_timestep_restart /= riverflow_timestep) THEN
        finfl = finfl * REAL(riverflow_timestep_restart,dp) * div_riverflow_timestep
    ENDIF

    status = NF_GET_ATT_INT (fileID, NF_GLOBAL,'initial_date', buf)
    yyyymmdd = buf(1)
    IF (status /= NF_NOERR)  THEN
      initial_date = start_date    ! old restart files or first year
    ELSE
      CALL inp_convert_date(yyyymmdd, 0, initial_date)
    ENDIF

    CALL IO_close(hdfile)
    DEALLOCATE (hd_dp_read)

  END SUBROUTINE read_hydrology

  SUBROUTINE hydrology_restart

    !
    ! **** Routine that writes the restart file for the HD model
    !
    ! ***** Version 1.0 - Oktober 1999
    !            Programmed and developed by Stefan Hagemann, MPI
    !
    !            Remark: Input data of Runoff and Drainage should have the
    !                       unit m/s.
    !
    ! ***** Version 1.1 - January 2001
    !       ECHAM5-Version
    !
    ! S.Legutke MPI M&D, Jan 2002, deallocate variables at end of
    !                              rerun cycle
    !
    ! ****** list of variables
    !
    !  frfmem(nl, nb, nmemrf) = Intermediate content of reservoir cascade
    !                           for the inflows per Gridbox (=5)
    !
    !  flfmem(nl, nb) = Intermediate content of linear reservoir for
    !                           Overland Flow
    !
    !   fgmem = Array of linear baseflow reservoir (Intermediate content)
    !           At Initialization it has the unit [m^3/s]
    !               (daily time step inherently implemented)
    !   finfl = Inflow data array for each gridbox for time step nstep
    !

    TYPE (FILE_INFO)  :: restartfile

    INTEGER :: nvarid, fileID, i, dims(2), xdimid, ydimid, xvarid, yvarid
    INTEGER :: istep
    INTEGER :: yyyymmdd, hhmmss

    CHARACTER(len=80) :: fname, string
    CHARACTER(len=7)  :: varname

    REAL(dp), ALLOCATABLE :: lons(:)
    REAL(dp), ALLOCATABLE :: lats(:)

    IF (p_parallel_io) THEN
      ALLOCATE(lons(grid_hd%nlon))
      ALLOCATE(lats(grid_hd%nlat))

      fname = 'hdrestart.nc'

      istep = get_time_step() + 1               ! get_time_step seems yo yield the step-1

      !    Open HD model restart file

      restartfile%opened = .FALSE.
      CALL IO_open (fname, restartfile, IO_WRITE)
      WRITE (message_text,*) 'Writing hdrestart to file ', TRIM(fname)
      CALL message('hydrology_restart', message_text)

      fileID = restartfile%file_id

      CALL IO_put_att_int (fileID, NF_GLOBAL, 'istep', istep)
      CALL IO_put_att_int (fileID, NF_GLOBAL, 'hd_steps_per_day', hd_steps_per_day)
      CALL IO_put_att_int (fileID, NF_GLOBAL, 'riverflow_timestep', riverflow_timestep)
      CALL out_convert_date (initial_date, yyyymmdd, hhmmss) 
      CALL IO_put_att_int (fileID, NF_GLOBAL, 'initial_date', yyyymmdd)
      CALL out_convert_date (current_date, yyyymmdd, hhmmss) 
      CALL IO_put_att_int (fileID, NF_GLOBAL, 'current_date', yyyymmdd)

      CALL IO_def_dim (fileID, 'lon', grid_hd%nlon, xdimid)
      CALL IO_def_dim (fileID, 'lat', grid_hd%nlat, ydimid)

      dims(1) = xdimid
      CALL IO_def_var (fileID, 'lon', NF_DOUBLE, 1, dims, xvarid)
      dims(1) = ydimid
      CALL IO_def_var (fileID, 'lat', NF_DOUBLE, 1, dims, yvarid)

      string = 'degrees_east'
      CALL IO_put_att_text (fileID, xvarid, 'units', string)
      string = 'Longitude'
      CALL IO_put_att_text (fileID, xvarid, 'long_name', string)
      string = 'degrees_north'
      CALL IO_put_att_text (fileID, yvarid, 'units', string)
      string = 'Latitude'
      CALL IO_put_att_text (fileID, yvarid, 'long_name', string)

      dims(1) = xdimid
      dims(2) = ydimid

      CALL IO_def_var (fileID, 'FLFMEM', NF_DOUBLE, 2, dims, nvarid)
      string = 'Linear overlandflow reservoir'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm3 d s-1'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 710)

      varname = 'FRFMEM'
      DO i=1, nmemrf
        WRITE(varname(7:7), '(i1)') i
        CALL IO_def_var (fileID, varname, NF_DOUBLE, 2, dims, nvarid)
        string = 'Inflow reservoir cascade'
        CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
        string = 'm3 d s-1'
        CALL IO_put_att_text (fileID, nvarid, 'units', string)
        CALL IO_put_att_int (fileID, nvarid, 'code', 710+i)
      ENDDO

      CALL IO_def_var (fileID, 'FGMEM', NF_DOUBLE, 2, dims, nvarid)
      string = 'Linear baseflow reservoir'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm3 d s-1'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 716)

      CALL IO_def_var (fileID, 'FINFL', NF_DOUBLE, 2, dims, nvarid)
      string = 'Inflow for each gridbox'
      CALL IO_put_att_text (fileID, nvarid, 'long_name', string)
      string = 'm3 s-1'
      CALL IO_put_att_text (fileID, nvarid, 'units', string)
      CALL IO_put_att_int (fileID, nvarid, 'code', 717)

      CALL IO_enddef (fileID)

      DO i = 1, grid_hd%nlon
        lons(i) = REAL(i,dp)*grid_hd%resolution+grid_hd%origin_lon-0.5_dp*grid_hd%resolution
        IF (lons(i) >= 180.0_dp) lons(i) = lons(i)-360.0_dp
      END DO
      DO i = 1, grid_hd%nlat
        lats(i) = grid_hd%origin_lat-REAL(i,dp)*grid_hd%resolution+0.5_dp*grid_hd%resolution
      END DO

      CALL IO_put_var_double (fileID, xvarid, lons)
      CALL IO_put_var_double (fileID, yvarid, lats)

      CALL IO_inq_varid (fileID, 'FLFMEM', nvarid)
      CALL IO_put_var_double (fileID, nvarid, flfmem)

      DO i=1, nmemrf
        WRITE(varname(7:7), '(i1)') i
        CALL IO_inq_varid (fileID, varname, nvarid)
        CALL IO_put_var_double (fileID, nvarid, frfmem(:,:,I))
      ENDDO

      CALL IO_inq_varid (fileID, 'FGMEM', nvarid)
      CALL IO_put_var_double (fileID, nvarid, fgmem)

      CALL IO_inq_varid (fileID, 'FINFL', nvarid)
      CALL IO_put_var_double (fileID, nvarid, finfl)

      CALL IO_close (restartfile)

      ! close output streams

      IF (nhd_diag > 0) CALL hd_close_timeseries
      IF (lhd_highres) CALL hd_highres_close(highres_file_id)

      DEALLOCATE (lons)
      DEALLOCATE (lats)
    END IF

  END SUBROUTINE hydrology_restart

  SUBROUTINE cleanup_hydrology

    CALL cleanup_hd_slm_invariants

    IF (p_parallel_io) THEN

      DEALLOCATE  (alf_k)
      DEALLOCATE  (alf_n)
      DEALLOCATE  (arf_k)
      DEALLOCATE  (arf_n)
      DEALLOCATE  (agf_k)
      DEALLOCATE  (agf_n)
      DEALLOCATE  (fdir)
      DEALLOCATE  (hd_lsm)
      DEALLOCATE  (finfl)
      DEALLOCATE  (fgmem)
      DEALLOCATE  (frfmem)
      DEALLOCATE  (flfmem)
      DEALLOCATE  (hd_area)
      DEALLOCATE  (oclook_cache)
      IF (is_coupled_run() .AND. lcoupling_atm) THEN
        DEALLOCATE  (runoff_s)
        DEALLOCATE  (runoff_dr)
      ELSE
        IF (nremap.EQ.0 .OR. nremap.EQ.3) THEN
          DEALLOCATE  (runoff_s)
          DEALLOCATE  (runoff_dr)
        ENDIF
      ENDIF
      IF (is_coupled_run()) DEALLOCATE  (fdir_hd) 
!Hag
      IF (lhd_rout) THEN
        DEALLOCATE  (filnew)
        DEALLOCATE  (fibnew)
      ENDIF
    END IF

    DEALLOCATE (gl_aros)
    DEALLOCATE (gl_adrain)
    DEALLOCATE (gl_disch)
    DEALLOCATE (gl_slm)
    DEALLOCATE (gl_alake)
    DEALLOCATE (gl_awfre)
    DEALLOCATE (gl_apmecal)
    DEALLOCATE (gl_glac)

    DEALLOCATE (input_overlandflow)
    DEALLOCATE (input_baseflow)
    DEALLOCATE (overlandflow)
    DEALLOCATE (baseflow)
    DEALLOCATE (riverflow)
    DEALLOCATE (overlandplusbaseflow)


  END SUBROUTINE cleanup_hydrology

  SUBROUTINE hydrology_model(slm_offline, alake_offline, glac_offline, &
       aros_offline, adrain_offline, disch_offline, awfre_offline, &
       apmecal_offline)

    ! HD Model - Constants and Switches
    !
    ! **** Global/Regional Discharge Simulation as Subroutine for ECHAM5
    !
    !
    ! ***** Version 1.0 - November 1999
    !   Programmed and Developed by Stefan Hagemann, MPI
    !   Program code is based on Offline-Version of the HD model
    !   which is also regionally applicable. (regsim.f)
    !
    !   Anmerkung: Input data of Runoff und Drainage should have the unit m/s.
    !
    ! **** Remarks: Changes with regard to offline version
    !   Runoff array is now passed into routine instead of reading it
    !   in kasglob via echread. In echread, now named hdech and called
    !   before kasglob, only the transformation of the runoff array
    !   to the resolution of 0.5 degree is done if necessary.
    !   Parameter/Variables luinp, area are deleted from kasglob.
    !
    !   Since the input array to be transformed is only passed to echread
    !   but not read in ECHREAD itself, in ECHREAD only
    !   the transformation of the input array to 0.5 degree is done.
    !   Therefor, new calling parameter and routine names are set/given.
    !   old: CALL echread(luinp, ihead, tocode, istep, ique)
    !   new: CALL hdech(code_t42, tocode_0.5grad, ique)
    !
    !
    ! ***** River Direction File (RDF) format:
    !
    !                    7  8  9
    !                     \ | /
    !                      \|/
    !                    4--5--6
    !                      /|\
    !                     / | \
    !                    1  2  3
    !
    !       Remark: Direction 5 = Discharge Trap
    !               Direction -1 = Ocean Point
    !
    ! ****** List of variables
    !
    !  lbase = Baseflow ON or OFF
    ! locean = Closure of Water budget for ocean coupling.
    !
    ! isolog = Logfile output into Iso file (ASCII file , two columns)
    !      0 = no, 1 = Bothnian Bay/Sea, 2 = Torneaelven, 3 = Global...
    !      4 = St.Lawrence, 5 = Paraguay 6 = Odra
    !lhd_que = Log-Output switch  (.false. No Log-Output to STDOUT)
    !
    !  istep = Chosen  time step for Reading of Input
    !     nl = Number of Longitudes
    !     nb = Number of Latitudes
    !
    !     riverflow_timestep = internal sub-timestep for riverflow computation
    !
    ! **** Global Arrays:
    !
    !  inpout_overlandflow, input_baseflow = local input data arrays for time step istep
    !  overlandplusbaseflow =  for time step istep
    ! overlandflow, baseflow, riverflow = local output data arrays for time step istep
    ! finfl = Inflow data array for each gridbox for time step istep
    !  fdir = River direction array
    !  hd_lsm = Land mask array
    ! alf_k = Array of retention constants k  - Overland flow [day]
    ! alf_n = Array of number of reservoirs n - Overland flow
    ! arf_k = Array of retention constants k  - Riverflow [day]
    ! arf_n = Array of number of reservoirs n - Riverflow
    ! agf_k = Array of retention constants k  - Baseflow [day]
    ! agf_n = Array of number of reservoirs n - Baseflow
    !
    !  frfmem(nl, nb, nmemrf) = Intermediate array of reservoir cascade for
    !                           the inflows per Gridbox (new: = nmemrf = 5)
    !
    !  flfmem(nl, nb, nmemlf) = Intermediate array of reservoir for
    !                           Surface Runoffs per Gridbox (new := nmemlf = 1)
    !
    ! fgmem = Array of linear baseflow reservoir (intermediate content)
    !         At initialization it has the unit [m^3/s]
    !  friv = Array of mean riverflow = Mean Inflow per Gridbox
    !
    ! hd_area(jb) = Array of gridbox arreas, Unit = [m^2]
    !
    !
    ! ***** Parameter and arrays of ECHAM grid
    !
    !  nlon = Longitudes of atmosphere grid
    !  ngl  = Latitudes of atmosphere grid
    !
    !  oclorg = Longitudinal origin of global atmosphere grid
    !  ocborg = Latitudinal origin of global atmosphere grid
    !  ocscal = resolution  = Latitudinal width of atmosphere gridbox in degree
    !
    !  aros   = atmospheric runoff array
    !  adrain = atmospheric drainage array
    !
    !  slm   = land sea mask on the Gaussian grid
    !  disch = discharge to the ocean (on the Gaussian grid)
    !  xresi = Residuum (Runoff+Drainage), which results from different
    !          land sea masks of the atmospheric grid and the HD model grid.
    !          In the former HD model versions, a start value was passed to the
    !          HD model that included further residual water terms that
    !          were distributed with the discharge to close the water
    !          balance in the coupled atmosphere ocean system.

    REAL(dp), INTENT(in),    OPTIONAL :: slm_offline(:,:)
    REAL(dp), INTENT(in),    OPTIONAL :: alake_offline(:,:)
    REAL(dp), INTENT(in),    OPTIONAL :: glac_offline(:,:)
    REAL(dp), INTENT(inout), OPTIONAL :: aros_offline(:,:)   ! INTENT(in)
    REAL(dp), INTENT(inout), OPTIONAL :: adrain_offline(:,:) ! INTENT(in)
    REAL(dp), INTENT(out),   OPTIONAL :: disch_offline(:,:)
    REAL(dp), INTENT(inout), OPTIONAL :: awfre_offline(:,:)
    REAL(dp), INTENT(in),    OPTIONAL :: apmecal_offline(:,:)

    REAL(dp) :: xresi
    REAL(dp) :: water_budget           ! budget of all water within the HD model
    REAL(dp) :: conservation_test      ! water budget difference at beginning and end of routine
    REAL(dp) :: sum_dis_to_ocean       ! Sum of discharge going to ocean

! Ha Ho-Hagemann: check the conservation since the second timestep {
    INTEGER :: istep    
! Ha Ho-Hagemann: check the conservation since the second timestep }

    !  Parameter and switches

    INTEGER :: jl, jb, jg, isub, ndum

    ! ECHAM grid characteristics
    ! Origin coordinate (usually upper left corner) & resolution
    ! Grid box centre at Longitude, Northern Border at latitude
    REAL(dp) :: ocscal

    REAL(dp), POINTER :: gl(:,:)
 
    IF (.NOT. is_coupled_run() .OR. .NOT. lcoupling_atm) THEN

       gl_slm = slm_offline
       gl_alake = alake_offline
       gl_glac = glac_offline
       gl_aros = aros_offline
       gl_adrain = adrain_offline
       gl_awfre = awfre_offline
       gl_apmecal = apmecal_offline

    END IF

    ! from now on only work on IO node

    IF (p_parallel_io) THEN

      ! initializations
      xresi = 0.0_dp                ! no longer needed 

      water_to_ocean(:,:) = 0._dp   ! water that is not handled by the HD model
      baseflow(:,:) = 0._dp

      ocscal = fullcirc / REAL(nlon, dp)

      ! initialization of water conservation test

!      WRITE (0,*) 'Water budget change: water_budget check 1:'
!      WRITE (0,*) 'SUM(gl_awfre*SPREAD(gridarea,1,nlon)*(1-gl_slm))=', SUM(gl_awfre*SPREAD(gridarea,1,nlon)*(1-gl_slm))
!      WRITE (0,*) 'SUM(gl_disch*SPREAD(gridarea,1,nlon))=', SUM(gl_disch*SPREAD(gridarea,1,nlon))
!      WRITE (0,*) 'SUM(flfmem)=', SUM(flfmem)
!      WRITE (0,*) 'SUM(fgmem)=', SUM(fgmem)
!      WRITE (0,*) 'SUM(frfmem)=', SUM(frfmem)
!      WRITE (0,*) 'SUM(finfl)=', SUM(finfl)
!      WRITE (0,*) 'xresi=', xresi
!      WRITE (0,*) 'SUM(gl_aros*SPREAD(gridarea,1,nlon))=', SUM(gl_aros*SPREAD(gridarea,1,nlon))
!      WRITE (0,*) 'SUM(gl_adrain*SPREAD(gridarea,1,nlon))=', SUM(gl_adrain*SPREAD(gridarea,1,nlon))
!      WRITE (0,*) 'SUM(gl_apmecal*SPREAD(gridarea,1,nlon))=', SUM(gl_apmecal*SPREAD(gridarea,1,nlon))

!      WRITE (0,*) 'SUM(runoff_s * SPREAD(hd_area,1,nl))=', SUM(runoff_s * SPREAD(hd_area,1,nl))
!      WRITE (0,*) 'SUM(runoff_dr * SPREAD(hd_area,1,nl))=', SUM(runoff_dr * SPREAD(hd_area,1,nl))

      IF (nremap.EQ.1 .OR. nremap.EQ.2) THEN
        water_budget = &
             ! runoff: from JSBACH (land grid cells)
             SUM(gl_aros*SPREAD(gridarea,1,nlon)) &
             ! P-E over glaciers
           + SUM(gl_apmecal*SPREAD(gridarea,1,nlon)) &
             ! drainage: from JSBACH (land grid cells)
           + SUM(gl_adrain*SPREAD(gridarea,1,nlon)) &
             ! fresh water flux (awfre): defined on grid cells with at least a 
             ! small water fraction (SLF<1). It comprises rain and snow, and
             ! evaporation over water (evapw). The values of awfre over lakes 
             ! will be added to the runoff.
           + SUM(gl_awfre*SPREAD(gridarea,1,nlon)*(1-gl_slm))
             ! water in the overland flow reservoirs
      ELSE

        water_budget = 0.0_dp
        DO jb = 1, grid_hd%nlat
          water_budget = &
            water_budget + &
            (SUM(runoff_s(:,jb)) + SUM(runoff_dr(:,jb))) * hd_area(jb) 
        END DO
        ! water_budget = SUM(runoff_s * SPREAD(hd_area,1,nl)) &
        !              + SUM(runoff_dr * SPREAD(hd_area,1,nl))
      ENDIF
      IF (ldebughd) THEN
        WRITE (message_text,*) 'Sum of inputs (surface runoff + drainage): ', water_budget,' m3/s'
        CALL message ('hydrology_model',message_text)
      ENDIF
      water_budget = water_budget &
           + SUM(flfmem) &
             ! water in the baseflow reservoir
           + SUM(fgmem) &
             ! water in the riverflow reservoirs
           + SUM(frfmem) &
             ! inflow data from the hd restart file
           + SUM(finfl)
      conservation_test=water_budget

      IF (ldebughd) THEN
         WRITE (message_text,*) '  Water budget at start: ', water_budget,' m3/s'
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) ' Overlandflow reservoir: ', SUM(flfmem)
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) '     Baseflow reservoir: ', SUM(fgmem)
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) '  River flow reservoirs: ', SUM(frfmem)
         CALL message ('hydrology_model',message_text)
         WRITE (message_text,*) '           River inflow: ', SUM(finfl)
         CALL message ('hydrology_model',message_text)
      END IF

      ! P-E on lakes is added to the surface runoff to be transported by the HD model
      !     Note: gl_aros (from jsbach) is only calculated on land, not over lakes.

      DO jg = 1, ngl
        DO jl = 1, nlon
          IF (gl_alake(jl,jg) .GT. 0.5_dp) THEN
            gl_aros(jl,jg) = gl_aros(jl,jg) + gl_awfre(jl,jg)
            gl_awfre(jl,jg) = 0._dp
          END IF
        END DO
      END DO

      !! Uwe Mikolajewicz, 2009/10/27
      !! Put P-E on glaciers into surface runoff field
      !! disable de facto the old glacier calving by setting input field to 0!!!
      !! requires ability of the HD model to transport negative runoff,
      !! which is given in ECHAM5.
      DO jg = 1, ngl
        DO jl = 1, nlon
          IF (gl_glac(jl,jg) .GT. 0.5_dp)THEN
            gl_aros(jl,jg) = gl_aros(jl,jg) + gl_apmecal(jl,jg)
            gl_apmecal(jl,jg) = 0.0_dp
          END IF
        END DO
      END DO


      ! ----------
      !  1 Runoff
      ! ----------

      IF (ldebughd) THEN
        IF (nremap.EQ.0 .OR. nremap.EQ.3) THEN
          WRITE (message_text,*) 'conservation test 1: No remapping for surface runoff'
        ELSE
          WRITE (message_text,*) 'conservation test 1: ', SUM(gl_aros*SPREAD(gridarea,1,nlon)), &
                                                       ' (SUM(gl_aros*SPREAD(gridarea,1,nlon)))'
        ENDIF
        CALL message ('hydrology_model',message_text)
      END IF

      ! Interpolation of the runoff from the source (e.g. ECHAM) to the HD model grid
      IF (is_coupled_run() .AND. lcoupling_atm) THEN
        ! receive from CCLM
        input_overlandflow (:,:) = runoff_s(:,:)
      ELSE
        IF (nremap.EQ.1) THEN
          CALL hd_remap(nlon, ngl, gl_aros, grid_hd%nlon, grid_hd%nlat, locean, SPREAD(gridarea,1,nlon), &
                      SPREAD(hd_area,1,grid_hd%nlon), input_overlandflow)
        ELSE IF (nremap.EQ.0) THEN
!------ Ha Ho-Hagemann read directly CCLM on the HD's grid {
          input_overlandflow(:,:) = runoff_s(:,:)
!        write(0,*) input_overlandflow
!------ Ha Ho-Hagemann read directly CCLM on the HD's grid }
        ELSE IF (nremap.EQ.2) THEN
!hag    *** Prep. for 5 Min. HD version      
          CALL hd_remap_05degto5min(nlon, ngl, gl_aros, grid_hd%nlon, grid_hd%nlat, input_overlandflow)
        ELSE IF (nremap.EQ.3) THEN
!       *** read directly on HD grid, but with longitudinal shift of 180 degree
          ndum = nlon/2
          DO jl = 1, ndum
            input_overlandflow (jl+ndum,:) = runoff_s(jl,:)
            input_overlandflow (jl,:) = runoff_s(jl+ndum,:)
          ENDDO
        ENDIF
      ENDIF
      !  Attention: Runoff in m/s --> Trafo with  AREA to m^3/s

!!    WRITE(message_text,*) 'hd_area= ', hd_area(1),hd_area(2),hd_area(3),' ... ',hd_area(grid_hd%nlat) 
!!    CALL message('before runoff conversion', message_text)

      DO jl = 1,grid_hd%nlon
         input_overlandflow(jl,:) = input_overlandflow(jl,:) * hd_area(:)
      END DO

      IF (diag_water_budget) THEN
         WRITE (message_text,*) 'DIAGWB: Global Surface runoff input on HD: ', &
                 SUM(input_overlandflow), ' m3/s'
         CALL message ('hydrology_model',message_text)
      ENDIF

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 1: ', SUM(input_overlandflow) + xresi,' (SUM(input_overlandflow) + xresi)'
         CALL message ('hydrology_model',message_text)
      END IF

      overlandflow(:,:) = 0.0_dp

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 2: ', &
           SUM(flfmem) + SUM(input_overlandflow),' (SUM(flfmem) + SUM(input_overlandflow))'
         CALL message ('hydrology_model',message_text)
      END IF

      ! kasglob handles HD land points with positive reservoir numbers, only.
      ! Land points without outflow (fdir=5) have reservoir number zero.
      ! Besides, some HD ocean points have runoff values, due to the missmatch
      ! of the ECHAM and the HD land sea masks. This water will be given to
      ! the ocean directly to close the water balance.

      CALL kasglob(input_overlandflow, overlandflow, alf_k, alf_n, hd_steps_per_day, flfmem, &
           alf_n_kas)
      WHERE (fdir == -1 .OR. fdir == 0 .OR. fdir == 5)
          water_to_ocean(:,:) = water_to_ocean(:,:) + input_overlandflow(:,:)
      END WHERE     

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 2: ', &
              SUM(flfmem) + SUM(overlandflow) + SUM(water_to_ocean), &
	   ' (SUM(flfmem) + SUM(overlandflow) + SUM(water_to_ocean))'
         CALL message ('hydrology_model',message_text)
      END IF

      ! ------------
      !  2 Drainage
      ! ------------

      IF (lbase) THEN

        IF (ldebughd) THEN
          IF (nremap.EQ.0 .OR. nremap.EQ.3) THEN
            WRITE (message_text,*) 'conservation test 3: No remapping for baseflow input'
          ELSE
            WRITE (message_text,*) 'conservation test 3: ', &
                  SUM(gl_adrain*SPREAD(gridarea,1,nlon)), &
	       ' (SUM(gl_adrain*SPREAD(gridarea,1,nlon)))'
          ENDIF
          CALL message ('hydrology_model',message_text)
        END IF


        IF (is_coupled_run() .AND. lcoupling_atm) THEN

          ! receive from CCLM
          input_baseflow (:,:) = runoff_dr (:,:)

        ELSE
          ! Interpolation of the drainage/subsurface runoff from the source
          ! (e.g. ECHAM) to the HD model grid
          IF (nremap.EQ.1) THEN
!vg.mine>>
            CALL hd_remap(nlon, ngl, gl_adrain, grid_hd%nlon, grid_hd%nlat, locean, &
                          SPREAD(gridarea,1,nlon), SPREAD(hd_area,1,grid_hd%nlon), &
                          input_baseflow)
!vg<<
          ELSE IF (nremap.EQ.0) THEN

!----       Ha Ho-Hagemann read directly CCLM on the HD's grid {
            input_baseflow (:,:) = runoff_dr(:,:)
!----       Ha Ho-Hagemann read directly CCLM on the HD's grid }
          ELSE IF (nremap.EQ.2) THEN
!hag        *** Prep. for 5 Min. HD version      
            CALL hd_remap_05degto5min(nlon, ngl, gl_adrain, grid_hd%nlon, grid_hd%nlat, input_baseflow)
          ELSE IF (nremap.EQ.3) THEN
!           *** read directly on HD grid, but with longitudinal shift of 180 degree
            ndum = nlon/2
            DO jl = 1, ndum
              input_baseflow (jl+ndum,:) = runoff_dr(jl,:)
              input_baseflow (jl,:) = runoff_dr(jl+ndum,:)
            ENDDO
          ENDIF
        ENDIF

         ! *** Attention: Drainage in m/s --> Trafo with AREA to m^3/s

         DO jl = 1, grid_hd%nlon
            input_baseflow(jl,:) = input_baseflow(jl,:) * hd_area(:)
         ENDDO

         IF (diag_water_budget) THEN
            WRITE (message_text,*) 'DIAGWB: Global Subsurface runoff input on HD: ', &
                 SUM(input_baseflow), ' m3/s'
            CALL message ('hydrology_model',message_text)
         ENDIF

         IF (ldebughd) THEN
            WRITE (message_text,*) 'conservation test 3: ', &
                 SUM(input_baseflow) + SUM(water_to_ocean) + xresi, &
	      ' (SUM(input_baseflow) + SUM(water_to_ocean) + xresi)'
            CALL message ('hydrology_model',message_text)
            WRITE (message_text,*) 'conservation test 4: ', &
                 SUM(fgmem) + SUM(input_baseflow) + SUM(water_to_ocean), &
	      ' (SUM(fgmem) + SUM(input_baseflow) + SUM(water_to_ocean))'
            CALL message ('hydrology_model',message_text)
         END IF

         ! Some HD ocean points have drainage values, due to the missmatch of the ECHAM
         ! and the HD land sea masks. This water is ignored by kasglob. It is given to
         ! the ocean directly.

         CALL kasglob(input_baseflow, baseflow, agf_k, agf_n, hd_steps_per_day, fgmem, agf_n_kas)
         WHERE (hd_lsm < 0.5)
            water_to_ocean(:,:) = water_to_ocean(:,:) + input_baseflow(:,:)
         END WHERE

         IF (ldebughd) THEN
            WRITE (message_text,*) 'conservation test 4: ', &
                 SUM(fgmem) + SUM(baseflow) + SUM(water_to_ocean), &
	      ' (SUM(fgmem) + SUM(baseflow) + SUM(water_to_ocean))'
            CALL message ('hydrology_model',message_text)
         END IF

      END IF

      ! -------------
      !  3 Riverflow
      ! -------------

      ! initialization of riverflow

      friv(:,:) = 0.0_dp

      ! input for routing: overlandflow (overlandflow) + baseflow (outflow from drainage)

      overlandplusbaseflow(:,:) = overlandflow(:,:) + baseflow(:,:)

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 5: ', &
              SUM(finfl) + SUM(frfmem) + SUM(overlandplusbaseflow) + SUM(water_to_ocean), &
	   ' (SUM(finfl) + SUM(frfmem) + SUM(overlandplusbaseflow) + SUM(water_to_ocean))'
         CALL message ('hydrology_model',message_text)
      END IF

      ! computation of riverflow in internal sub-time steps

      DO isub = 1, riverflow_timestep

         ! computing riverflow with input finfl from preceeding sub-time step

         ! kasglob handles HD land points with positive reservoir numbers, only.
         ! Land points without outflow (fdir=5) have reservoir number zero.
         ! finfl has non-zero values in regions with fdir=5.
!
!        *** Does riverflow velocity depend on discharge volume, i.e. finfl
         IF (irf_vel.eq. 1 ) THEN
           CALL kasglob_qrf(finfl, riverflow, arf_k, arf_n, riverflow_steps_per_day, frfmem, &
              arf_n_kas)
         ELSE
           CALL kasglob(finfl, riverflow, arf_k, arf_n, riverflow_steps_per_day, frfmem, &
              arf_n_kas)
         ENDIF

         WHERE (fdir == -1 .OR. fdir == 0 .OR. fdir == 5)
            water_to_ocean(:,:) = water_to_ocean(:,:) + finfl(:,:)
         END WHERE
!Hag
!        *** Conduct routing in a subroutine for easier exchange of methods
         IF (lhd_rout) THEN
           call routing_via_index(overlandplusbaseflow, riverflow, finfl)
         ELSE
           call routing(overlandplusbaseflow, riverflow, finfl)
         ENDIF  

         DO jb = 1, grid_hd%nlat
           DO jl = 1, grid_hd%nlon
             ! non land point
             IF (fdir(jl, jb) == -1 .OR. &
                 fdir(jl, jb) == 0 .OR. &
                 fdir(jl, jb) == 5) THEN
               water_to_ocean(jl, jb) = water_to_ocean(jl, jb) + finfl(jl, jb)
               finfl(jl, jb) = 0._dp
             END IF
           END DO
         END DO

         friv(:,:) = friv(:,:) + finfl(:,:)

      ENDDO   ! loop over sub-time steps

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 5: ', &
              SUM(frfmem) + SUM(finfl) + SUM(water_to_ocean), &
	   ' (SUM(frfmem) + SUM(finfl) + SUM(water_to_ocean))'
         CALL message ('hydrology_model',message_text)
      END IF

      ! HD outflow diagnostics
      ! friv is defined on HD land (without internal drainage) and water_to_ocean is defined on hd 
      ! ocean (with internal drainage). The outflow diagnostics take into account both of these arrays.

      IF (nhd_diag /= 0 .OR. lhd_highres) CALL hydrology_diags(nhd_diag, water_to_ocean + friv)

      ! Conversion from the HD to the ECHAM Grid

      IF (ldebughd) THEN
         WRITE (message_text,*) 'conservation test 6: ', SUM(water_to_ocean) + xresi, &
                              	                         ' (SUM(water_to_ocean) + xresi)'
         CALL message ('hydrology_model',message_text)
      END IF

      IF (locean) THEN

         CALL hydrology_to_ocean(nlon, ngl, &
              water_to_ocean, fdir, gl_disch, xresi)
         IF (ldebughd) THEN
            WRITE (message_text,*) 'conservation test 6: ', SUM(gl_disch), ' (SUM(gl_disch))'
            CALL message ('hydrology_model',message_text)
         END IF
         ! Convert discharge from m**3/s to m/s for ocean model
         DO jg = 1, ngl
           DO jl = 1, nlon
             gl_disch(jl, jg) = gl_disch(jl, jg)/gridarea(jg)
           END DO
         END DO

         sum_dis_to_ocean = SUM(gl_disch*SPREAD(gridarea,1,nlon))

      ELSE
         sum_dis_to_ocean = SUM(water_to_ocean) + xresi
      ENDIF

         water_budget = &
              ! discharge to the ocean
            sum_dis_to_ocean &
              ! water in the overland flow reservoirs
            + SUM(flfmem) &
              ! water in the baseflow reservoir
            + SUM(fgmem) &
              ! water in the riverflow reservoirs
            + SUM(frfmem) &
              ! inflow data
            + SUM(finfl)

         DO jb = 1, ngl
           ! Rain, snow and evaporation over water from ECHAM/JSBACH. P-E from
           ! lakes was given to the runoff.
           water_budget = &
             water_budget + &
              SUM(gl_awfre(:,jb)*(1-olm(:,jb)))*gridarea(jb)
         END DO


!         WRITE (0,*) 'Water budget change: water_budget check 2:'
!	 WRITE (0,*) 'SUM(gl_awfre*SPREAD(gridarea,1,nlon)*(1-olm))=', SUM(gl_awfre*SPREAD(gridarea,1,nlon)*(1-olm))
!	 WRITE (0,*) 'SUM(gl_disch*SPREAD(gridarea,1,nlon))=', SUM(gl_disch*SPREAD(gridarea,1,nlon))
!	 WRITE (0,*) 'SUM(flfmem)=', SUM(flfmem)
!	 WRITE (0,*) 'SUM(fgmem)=', SUM(fgmem)
!	 WRITE (0,*) 'SUM(frfmem)=', SUM(frfmem)
!	 WRITE (0,*) 'SUM(finfl)=', SUM(finfl)
!	 WRITE (0,*) 'xresi=', xresi
!         WRITE (0,*) 'Water budget change: conservation_test1=', conservation_test,' m3/s', &
!	                       ' water_budget=',water_budget,' m3/s'

         conservation_test = conservation_test - water_budget

         WRITE (message_text,*) 'Water budget change: conservation_test2=', conservation_test,' m3/s'

         CALL message ('hydrology_model',message_text)

!!!!!!!!!!!! due to the crash in the standalone HD, temporially comment these lines below
!---- Ha Ho-Hagemann {
!     IF (is_coupled_run()) THEN
! Ha Ho-Hagemann: check the conservation since the second timestep {
!       istep = get_time_step()
!       if (istep .gt. 1) then
! Ha Ho-Hagemann: check the conservation since the second timestep }
!           IF (ABS(conservation_test/water_budget) > 1.E-7_dp) THEN
!              WRITE (message_text,*) 'Water conservation problem: budget change: ', &
!                   conservation_test,' m3/s'
!              CALL finish ('hydrology_model',message_text)
!           END IF
! Ha Ho-Hagemann: check the conservation since the second timestep {
!       endif
! Ha Ho-Hagemann: check the conservation since the second timestep }
!     ENDIF
!---- Ha Ho-Hagemann }

         IF (diag_water_budget) THEN
            WRITE (message_text,*) 'DIAGWB: Global discharge into the ocean without sink redistribution: ', &
                 sum_dis_to_ocean, ' m3/s'
            CALL message ('hydrology_model',message_text)
            IF (locean) THEN  
                WRITE (message_text,*) '       precip-evap (on ocean): ', &
                    SUM((gl_awfre)*SPREAD(gridarea,1,nlon)*(1._dp-olm)), ' m3/s'
                CALL message ('hydrology_model',message_text)
            ENDIF
         END IF

    END IF

    IF (is_coupled_run() .AND. lcoupling_oce) THEN
       disch_offline =  gl_disch
    END IF

  END SUBROUTINE hydrology_model

#ifndef HD_5MIN

  PURE INTEGER FUNCTION lon_to_idx(lon)
    REAL(dp), INTENT(IN) :: lon
    lon_to_idx = NINT((lon-grid_hd%origin_lon)/grid_hd%resolution + 1._dp)
  END FUNCTION lon_to_idx
  PURE INTEGER FUNCTION lat_to_idx(lat)
    REAL(dp), INTENT(IN) :: lat
    lat_to_idx = NINT(1._dp + (grid_hd%origin_lat-lat)/grid_hd%resolution)
  END FUNCTION lat_to_idx

  SUBROUTINE hydrology_diags(isolog, hd_out)
    INTEGER,  INTENT(in) :: isolog
    REAL(dp), INTENT(in) :: hd_out(grid_hd%nlon,grid_hd%nlat)

    REAL(dp) :: f1, f2
    INTEGER :: jl, jb
    INTEGER :: istep
    INTEGER, SAVE :: counter = 0

    !  Filling  F1, F2 at chosen coordinates with
    !  OUTFLOW per Gridbox --> FINP, not HD_OUT or FINFL
    !  Using grid box NW corner --> (x-x0)/res +1 with x in x0, x0+res, x0+2*res, 
    !  --> NINT in function is ok.

    istep = get_time_step()

    IF (isolog == 1) THEN

      ! *** Log Output for Inflow into Gulf of Bothnia
      ! *** Since INFLOW ==> FINFL bzw. HD_OUT!
      ! *** Bothnian Bay: B=65.5 ,L=21.5 .... (glob = (404, 50))

      jl = lon_to_idx(21.5_dp)
      jb = lat_to_idx(65.5_dp)
      f1 =  hd_out(jl,  jb  ) + hd_out(jl+1,jb  ) + hd_out(jl+2,jb  ) &
          + hd_out(jl+3,jb  ) + hd_out(jl+4,jb  ) + hd_out(jl+5,jb  ) &
          + hd_out(jl+6,jb  ) + hd_out(jl-1,jb+1) + hd_out(jl+5,jb+1) &
          + hd_out(jl-1,jb+2) + hd_out(jl+1,jb+2) + hd_out(jl+2,jb+2) &
          + hd_out(jl+3,jb+2) + hd_out(jl+4,jb+2) + hd_out(jl-2,jb+3) &
          + hd_out(jl-1,jb+4) + hd_out(jl,  jb+4)

      ! *** Bothnian Sea: B=63.5 ,L=19.0 ....

      f2 =  hd_out(jl-5,jb+4 ) + hd_out(jl-4,jb+4 ) + hd_out(jl-7,jb+5 ) &
          + hd_out(jl-8,jb+6 ) + hd_out(jl-2,jb+6 ) + hd_out(jl-8,jb+7 ) &
          + hd_out(jl-2,jb+7 ) + hd_out(jl-1,jb+7 ) + hd_out(jl-9,jb+8 ) &
          + hd_out(jl-1,jb+8 ) + hd_out(jl-8,jb+9 ) + hd_out(jl-1,jb+9 ) &
          + hd_out(jl-6,jb+10) + hd_out(jl-1,jb+10) + hd_out(jl,  jb+10) &
          + hd_out(jl,  jb+11) + hd_out(jl+1,jb+11) + hd_out(jl+2,jb+11)

    ELSE IF (isolog == 2 .OR. isolog == 3) THEN

      ! *** Torneaelven-Outflow = (22.0 E, 65.5 N), (22.5 E, 65.5 N)
      ! ***                       (23.5 E, 65.5 N)
      ! *** regional System:

      jl = lon_to_idx(22.0_dp)
      jb = lat_to_idx(65.5_dp)
      f1 = REAL(istep,dp)
      f2 = hd_out(jl,jb) + hd_out(jl+1,jb) + hd_out(jl+3,jb)

    ELSE IF (isolog == 4) THEN

      ! *** St.Lawrence-Outflow = (-71.5 W, 47.0 N)
      ! *** regional System: Measurement station at (-75.5 W, 45.5 N)

      jl = lon_to_idx(-71.5_dp)
      jb = lat_to_idx( 47.0_dp)
      f1 = hd_out(jl,jb)
      f2 = hd_out(jl-8,jb+3)

    ELSE IF (isolog == 5) THEN

      ! *** Paraguay-Outflow = (-59 W, -27 N)

      jl = lon_to_idx(-59.0_dp)
      jb = lat_to_idx(-27.0_dp)
      f1 = REAL(istep,dp)
      f2 = hd_out(jl,jb)

    ELSE IF (isolog == 6) THEN

      ! *** Oder-Outflow = (14.0 E, 54.5 N), Hohensaaten-Finow (14 E, 53W)

      jl = lon_to_idx(14.0_dp)
      jb = lat_to_idx(54.5_dp)
      f1 = hd_out(jl,jb)
      f2 = hd_out(jl+1,jb+2)+hd_out(jl+2,jb+2)

   ELSE IF (isolog == 7) THEN

      ! *** Elbe-Outflow = (8.5 E, 54.5 N), Neu-Darchau (10.5 E, 53.5 N)

      jl = lon_to_idx( 8.5_dp)
      jb = lat_to_idx(54.5_dp)
      f1 = hd_out(jl,jb)
      f2 = hd_out(jl+4,jb+2)

   ELSE IF (isolog == 8) THEN

      ! *** Oranje-Outflow = (-28.5 S, 16.0 E), Congo (-6.0 S, 12.0 E)

      jl = lon_to_idx( 16.0_dp)
      jb = lat_to_idx(-28.5_dp)
      f1 = hd_out(jl,jb)

      jl = lon_to_idx(12.0_dp)
      jb = lat_to_idx(-6.0_dp)
      f2 = hd_out(jl,jb)

   ELSE IF (isolog == 9) THEN

      ! *** Amudarya-Outflow (47) = (43 N, 59 E), Syrdarya (49) = (46 N, 62 E)

      jl = lon_to_idx(59.0_dp)
      jb = lat_to_idx(43.0_dp)
      f1 = hd_out(jl,jb)

      jl = lon_to_idx(62.0_dp)
      jb = lat_to_idx(46.0_dp)
      f2 = hd_out(jl,jb)

   ELSE IF (isolog == 10) THEN

      ! *** Lena-Outflow (40) = (72 N, 127 E), Ob (46) = (67 N, 71.5 E)

      jl = lon_to_idx(127.0_dp)
      jb = lat_to_idx( 72.0_dp)
      f1 = hd_out(jl,jb)

      jl = lon_to_idx(71.5_dp)
      jb = lat_to_idx(67.0_dp)
      f2 = hd_out(jl,jb)

   ELSE IF (isolog == 99) THEN

      ! *** user defined outflow coordinates (namelist hydrology_ctl)

      jl = lon_to_idx(fllog1)
      jb = lat_to_idx(fblog1)
      f1 = hd_out(jl,jb)

      jl = lon_to_idx(fllog2)
      jb = lat_to_idx(fblog2)
      f2 = hd_out(jl,jb)

   ENDIF

   IF (isolog > 0) CALL hd_write_timeseries (f1, f2)

   IF (lhd_highres) CALL hd_highres_write (highres_file_id, hd_out, counter, grid_hd%nlon, grid_hd%nlat, hd_area)

  END SUBROUTINE hydrology_diags
#else

  PURE INTEGER FUNCTION lon_to_idx(lon)
    REAL(dp), INTENT(IN) :: lon
    lon_to_idx = INT((lon-grid_hd%origin_lon)/grid_hd%resolution + 1.001_dp)
  END FUNCTION lon_to_idx
  PURE INTEGER FUNCTION lat_to_idx(lat)
    REAL(dp), INTENT(IN) :: lat
    lat_to_idx = INT(1.001_dp + (grid_hd%origin_lat-lat)/grid_hd%resolution)
  END FUNCTION lat_to_idx

  SUBROUTINE hydrology_diags(isolog, hd_out)
!
!   *** 5 Min. version of the routine -- HAG -- Dec. 2104
    INTEGER,  INTENT(in) :: isolog
    REAL(dp), INTENT(in) :: hd_out(grid_hd%nlon,grid_hd%nlat)

    REAL(dp) :: f1, f2
    INTEGER :: jl, jb
    INTEGER :: istep
    INTEGER, SAVE :: counter = 0

    !  Filling  F1, F2 at chosen coordinates with
    !  OUTFLOW per Gridbox --> FINP, not HD_OUT or FINFL
    !  Using Grid box centre  --> (x-x0)/res +1 with x in x0+0.5 res, x0+1.5 res ...
    !  INT goes down, NINT can go up or down, but should actually go down --> Use INT in function

    istep = get_time_step()

    IF (isolog == 7) THEN

      ! *** Elbe-Outflow = (9.875 E, 53.5416 N), Neu-Darchau (53.23 N 10.88 E)

      jl = lon_to_idx(  9.875_dp)
      jb = lat_to_idx(53.5416_dp)
      f1 = hd_out(jl,jb)
      f2 = hd_out(jl+12,jb+4)
   ELSE IF (isolog == 99) THEN

      ! *** user defined outflow coordinates (namelist hydrology_ctl)

      jl = lon_to_idx(fllog1)
      jb = lat_to_idx(fblog1)
      f1 = hd_out(jl,jb)

      jl = lon_to_idx(fllog2)
      jb = lat_to_idx(fblog2)
      f2 = hd_out(jl,jb)
   ELSE
      WRITE (message_text,*) '  HD Logoutput Type: ', isolog,  &
                 ' does not exist'
      CALL message('hydrology_diags_5min', message_text)
      CALL finish ('hydrology_diags_5min', 'run terminated.')
   ENDIF

   IF (isolog > 0) CALL hd_write_timeseries (f1, f2)

   IF (lhd_highres) CALL hd_highres_write (highres_file_id, hd_out, counter, grid_hd%nlon, grid_hd%nlat)

  END SUBROUTINE hydrology_diags
#endif

  SUBROUTINE kasglob(finp, ymod, a_k, a_n, steps_per_day, fmem, a_n_kas)
    !
    ! ***** Global Flow Simulation with the conceptual model reservoir cascade
    !   Program was partailly written following the routine lfsim (in gate.for)
    !   and funkas (in modfunct.for).
    !
    ! ***** Programmed and developed by Stefan Hagemann
    !
    ! ***** Version 2.2 - November 1995
    !   Finer resolution of system function for Riverflow
    !   Re-arranging of IN/OUTput-arrays and Passing within a
    !      single reservoir field fmem
    !
    ! ***** Version 3.0 - November 1995
    !   Computation of outflow via Differential Equation
    !   of lin. reservoir cascade, comprising of nn reservoirs with
    !   nn = INT(a_n) = INT(n)
    !
    ! ***** Version 3.1 - Februar 1996
    !   Implementation of possible computation of riverflow with mm
    !   sub (internal) time steps.
    !
    ! ***** Version 5.0 - Oktober 1999
    !   Runoff-Input-data are passed to kasglob instead of reading it within
    !   kasglob itself.
    !   Calling parameters/variables luinp and area are deleted.
    !

    INTEGER,  INTENT(in)           :: steps_per_day  ! number of routine calls per day
    REAL(dp), INTENT(in)           :: a_k(grid_hd%nlon, grid_hd%nlat)    ! array of k-parameter [day]
    REAL(dp), INTENT(in)           :: a_n(grid_hd%nlon, grid_hd%nlat)    ! array of n-parameter
    REAL(dp), INTENT(in)           :: finp(grid_hd%nlon, grid_hd%nlat)   ! input overlandflow/riverflow array
    REAL(dp), INTENT(inout)        :: ymod(grid_hd%nlon, grid_hd%nlat)   ! simulated overlandflow/riverflow
    REAL(dp), INTENT(inout)        :: fmem(:,:,:)    ! intermediate content of reservoir cascade

    TYPE(cart_xidx_2d), INTENT(in) :: a_n_kas(:)     ! kasglob list

    REAL(dp) :: akdiv, fdum, amod, divmm, fmd_sum

    INTEGER :: j, jl, jb
    INTEGER :: nn
    INTEGER :: nx, i, extlen, extelem

    divmm = 1._dp/steps_per_day

    nx = SIZE(a_n_kas)
    ! **** Computing modeled value at each grid point

    DO i = 1, nx
      jb = a_n_kas(i)%ilat
      jl = a_n_kas(i)%ilon
      extlen = a_n_kas(i)%extlen
      DO extelem = 1, extlen
        nn = NINT(a_n(jl,jb))
        amod = a_n_kas(i)%amod(extelem)
#if 0
        ! Test for the amod values defined in create_kasglob_list
        IF (ABS(amod) > ABS(a_k(jl,jb) * a_n(jl,jb)/AINT(a_n(jl,jb)))+EPSILON(1._dp)) THEN
          WRITE (message_text,*) 'amod problem: ', jl, jb, &
               amod, a_k(jl,jb) * a_n(jl,jb)/AINT(a_n(jl,jb)), extelem
          CALL message('hd: kasglob', message_text)
        END IF
#endif

        ! *** Dt=1 day ==> AKDIV = 1./ (AMOD+1.)
        akdiv = a_n_kas(i)%akdiv(extelem)*divmm
        fdum = finp(jl,jb)

        ! *** Nash-Cascade
        ! *** It is [AMOD] = day ==> AMOD(sec) = AMOD(day) * 1 day
        ! *** Remember: In principle,it is: FDUM = FINP * 1 day
        ! ***           ==> FMEM = x * 1 day
        ! ***           ==> FDUM = x * 1 day * 1 / AMOD(sec)
        ! ***                    = x * 1 day / (AMOD(day) * 1 day)
        ! ***                    = x / AMOD(day)
        ! ***           ==> FMEM = x * 1 day - FDUM * 1 day
        ! ***                    = (x - FDUM) * 1 day
        ! *** Outflow FDUM is computed correctly, Intermediate reservoir unit is
        ! *** a volume flow instead of a volume. This is to avoid
        ! *** back and forth multiplication with factor 1 day = 86400 sec

        DO j = 1, nn
          fmd_sum = fmem(jl,jb,j) + fdum
          fdum = fmd_sum * akdiv
          fmem(jl,jb,j) = fmd_sum - fdum
        END DO
        ymod(jl,jb) = fdum
        jl = jl + 1
      ENDDO
    ENDDO

  END SUBROUTINE kasglob

  SUBROUTINE kasglob_qrf(finp, ymod, a_k, a_n, steps_per_day, fmem, a_n_kas)
    !
    ! ***** Global Flow Simulation for discharge dependent riverflow velocity
    !       Version derived from and analogous to subroutine kasglob
    !
    ! ***** Programmed and developed by Stefan Hagemann
    !
    ! ***** Version 6.0 - Feb. 2016
    ! ***   amod and akdiv need to be calculated in routine and cannot be taken from a_n_kas
    !
    INTEGER,  INTENT(in)           :: steps_per_day  ! number of routine calls per day
    REAL(dp), INTENT(in)           :: a_k(grid_hd%nlon, grid_hd%nlat)    ! array of k-parameter [day]
    REAL(dp), INTENT(in)           :: a_n(grid_hd%nlon, grid_hd%nlat)    ! array of n-parameter
    REAL(dp), INTENT(in)           :: finp(grid_hd%nlon, grid_hd%nlat)   ! input overlandflow/riverflow array
    REAL(dp), INTENT(inout)        :: ymod(grid_hd%nlon, grid_hd%nlat)   ! simulated overlandflow/riverflow
    REAL(dp), INTENT(inout)        :: fmem(:,:,:)    ! intermediate content of reservoir cascade

    TYPE(cart_xidx_2d), INTENT(in) :: a_n_kas(:)     ! kasglob list

    REAL(dp) :: akdiv, fdum, amod, divmm, fmd_sum, arf_k_mod, zfak
    REAL(dp), PARAMETER :: zeps = 1.E-10
    REAL(dp) :: qref_day     ! equivalent of qref per sub time step = qrf_ref / steps_per_day

    INTEGER :: j, jl, jb
    INTEGER :: nn
    INTEGER :: nx, i, extlen, extelem

    divmm = 1._dp/steps_per_day
    qref_day = qrf_ref / steps_per_day  

    nx = SIZE(a_n_kas)
    ! **** Computing modeled value at each grid point

    DO i = 1, nx
      jb = a_n_kas(i)%ilat
      jl = a_n_kas(i)%ilon
      extlen = a_n_kas(i)%extlen
      DO extelem = 1, extlen
        nn = NINT(a_n(jl,jb))
        ! *** Dt=1 day ==> AKDIV = 1./ (AMOD+1.)
!!!        amod = a_n_kas(i)%amod(extelem)
!!!        akdiv = a_n_kas(i)%akdiv(extelem)*divmm

        IF (finp(jl,jb).GT. qref_day) THEN
          zfak = SQRT( SQRT(finp(jl,jb) / qref_day))
          arf_k_mod = a_k(jl,jb) / zfak
        ELSE
          arf_k_mod = a_k(jl,jb)           
        ENDIF
        amod = arf_k_mod * a_n(jl,jb)/AINT(a_n(jl,jb))
        akdiv = 1._dp / (amod + divmm) *divmm

        fdum = finp(jl,jb)

        ! *** Nash-Cascade
        ! *** It is [AMOD] = day ==> AMOD(sec) = AMOD(day) * 1 day
        ! *** Remember: In principle,it is: FDUM = FINP * 1 day
        ! ***           ==> FMEM = x * 1 day
        ! ***           ==> FDUM = x * 1 day * 1 / AMOD(sec)
        ! ***                    = x * 1 day / (AMOD(day) * 1 day)
        ! ***                    = x / AMOD(day)
        ! ***           ==> FMEM = x * 1 day - FDUM * 1 day
        ! ***                    = (x - FDUM) * 1 day
        ! *** Outflow FDUM is computed correctly, Intermediate reservoir unit is
        ! *** a volume flow instead of a volume. This is to avoid
        ! *** back and forth multiplication with factor 1 day = 86400 sec

        DO j = 1, nn
          fmd_sum = fmem(jl,jb,j) + fdum
          fdum = fmd_sum * akdiv
          fmem(jl,jb,j) = fmd_sum - fdum
        END DO
        ymod(jl,jb) = fdum
        jl = jl + 1
      ENDDO
    ENDDO

  END SUBROUTINE kasglob_qrf


  SUBROUTINE hydrology_to_ocean(nlon, nlat, friv, fdir, disch, xresi)

    !
    ! ******* This programs distributes the river discharge from the 0.5
    !         degree inflow points to the considered ocean gridbox
    !
    !  Programmed and developed by Stefan Hagemann, MPI
    !
    ! ***** Version 1.0 -- Oktober 1999
    !
    ! ***** Version 2.0 -- January 2001
    !     ECHAM5- Version incl. Gaussian latitudes
    !
    ! ****** List of Variables
    !
    !  friv = Inflow array on HD model grid
    !  fdir = River direction file that defines river mouthes (destinations) as 0
    !         on HD model grid
    !  xidb = Summation array of inflows, for which no inflowbox into the
    !         ocean was found, e.g. Kaspian Sea and Interior
    !         Drainage Basins
    !  any_ocinflow = Inflow-Point found on ocean grid .TRUE./.FALSE.
    !
    !  nlon = Longitudes of global ocean grid
    !  nlat = Latitudes of global ocean grid
    !  oclorg = Longitudinal origin of global ocean grid
    !  ocscal = Scale/Resolution = Width of Ocean Gridbox in degree
    !  philat = Gaussian latitude of global ocean grid (centre coordinates)
    !
    !   disch = Inflow array on Ocean grid
    !   xresi = Residuum (Runoff+Drainage), which results from different
    !           land sea masks of the atmospheric grid and the 0.5 degree grid.
    !           In the latest HD model version, a start value is passed to the
    !           HD model that may include further residual water terms that should
    !           distributed with the discharge to close the water balance in the
    !           coupled atmosphere ocean system.
    !           XIDB is added to Xresi.
    ! lhd_que = Log-Output switch (.FALSE. = No Log-output to STDOUT)
    !

    INTEGER,  INTENT(in)  :: nlon, nlat
    REAL(dp), INTENT(in)  :: friv(grid_hd%nlon,grid_hd%nlat)
    INTEGER,  INTENT(in)  :: fdir(grid_hd%nlon,grid_hd%nlat)
    REAL(dp), INTENT(in)  :: xresi
    REAL(dp), INTENT(out) :: disch(nlon,nlat)
#if 0
    REAL(dp) :: xjlat
#endif
    REAL(dp) :: xidb
    INTEGER :: jl,jb
    TYPE(cart_idx_2d) :: dest
    LOGICAL :: any_ocinflow

    disch(:,:) = 0.0_dp
    xidb = xresi
    any_ocinflow = .FALSE.

    ! ******* Loop over all inflow points

    DO jb = 1, grid_hd%nlat
       DO jl = 1, grid_hd%nlon

          ! HD ocean cell
          dest = oclook_cache(jl, jb)
          IF (dest%ilon /= -1) THEN
             disch(dest%ilon,dest%ilat) = disch(dest%ilon,dest%ilat) + friv(jl,jb)
             any_ocinflow = .TRUE.
          END IF

          ! internal drainage cells
          IF (fdir(jl,jb) == 5) THEN
             xidb = xidb + friv(jl,jb)
          END IF

      ENDDO
    ENDDO

    ! Distributing the water in XIDB to all Ocean Inflow Points
    ! Applying a weight to treat arid and humid regions differently

    IF (any_ocinflow) THEN
      disch(:,:) = disch(:,:)+disch(:,:)/SUM(disch(:,:)) * xidb
    ELSE
      WRITE(message_text,*) 'error no inflow points on ocean grid found'
      CALL message('hydrology_to_ocean', message_text)

    ENDIF
  END SUBROUTINE hydrology_to_ocean

  SUBROUTINE hydrology_slm_invariants(slm_offline, alake_offline)

    REAL(dp), INTENT(in), OPTIONAL :: slm_offline(:,:)    ! integer land sea mask
    REAL(dp), INTENT(in), OPTIONAL :: alake_offline(:,:)  ! fractional lake mask

    INTEGER  :: il, jb, i
    REAL(dp) :: ra, rb, rh, rd
    REAL(dp), POINTER :: gl(:,:)

    IF (PRESENT(slm_offline) .AND. PRESENT(alake_offline)) THEN
       gl_slm = slm_offline
       gl_alake = alake_offline
    END IF

    WRITE(message_text,*) 'philon= ', philon(1),philon(2),philon(3),' ... ',philon(nlon) 
    CALL message('hydrology_model', message_text)
    WRITE(message_text,*) 'philat= ', philat(1),philat(2),philat(3),' ... ',philat(ngl) 
    CALL message('hydrology_model', message_text)


    IF (p_parallel_io) THEN

       ! setup area in m^2 of the HD model internal grid

       ra = 2.0_dp*api*a*a/REAL(grid_hd%nlon,dp)
       rb = api/REAL(grid_hd%nlat,dp)
       rd = 0.5_dp*api
#ifndef HD_5MIN
       DO i = 1, grid_hd%nlat
          rh = SIN(-rd+(i-1)*rb)-SIN(-rd+i*rb)
          hd_area(i) = ABS(rh)*ra
       END DO
#endif

       ! Gaussian latitudes in degrees

       IF (ldebughd) THEN
          WRITE(message_text,*) 'philon= ', philon(1),philon(2),philon(3),' ... ',philon(nlon) 
          CALL message('hydrology_model', message_text)
          WRITE(message_text,*) 'philat= ', philat(1),philat(2),philat(3),' ... ',philat(ngl) 
          CALL message('hydrology_model', message_text)
       END IF

       ! Definition of ocean land mask, in contrast to slm with lakes represented as land

       ALLOCATE(olm(nlon,ngl))
       olm(:,:) = MERGE(gl_slm(:,:), 1._dp, gl_alake(:,:) < 0.5_dp)

       ! define array with number of reservoirs for the base flow.
       ! There is just 1 baseflow reservoir (linear) at each grid point. The array 
       ! is needed to be able to use routine kasglob.

       WHERE (hd_lsm > 0.5_dp)
          agf_n = 1._dp
       ELSEWHERE
          agf_n = 0._dp
       END WHERE

       ! check consistency of HD input data

       IF (ANY((hd_lsm > 0.5_dp .AND. (fdir == -1 .OR. fdir == 0)) &
            .OR. (hd_lsm < 0.5_dp .AND. (fdir /= -1 .AND. fdir /= 0 .AND. fdir /= 5)))) THEN
          CALL finish('hydrology_slm_invariants', &
               'hd slm and runoff directions do not match')
       END IF

       IF (ANY((hd_lsm > 0.5_dp .AND. arf_n < 0.5_dp .AND. fdir /= 5) &
            .OR. (hd_lsm < 0.5_dp .AND. arf_n > 0.5_dp))) THEN
          CALL finish('hydrology_slm_invariants', &
               'hd slm and riverflow reservoir numbers do not match')
       END IF

       IF (ANY((hd_lsm > 0.5_dp .AND. alf_n < 0.5_dp .AND. fdir /= 5) &
            .OR. (hd_lsm < 0.5_dp .AND. alf_n > 0.5_dp))) THEN
          CALL finish('hydrology_slm_invariants', &
               'hd slm and overlandflow reservoir numbers do not match')
       END IF

       IF (ANY((hd_lsm > 0.5_dp .AND. agf_n <= 0._dp) &
            .OR. (hd_lsm < 0.5_dp .AND. agf_n /= 0._dp))) THEN
          CALL finish('hydrology_slm_invariants', &
               'hd slm and baseflow reservoir numbers do not match')
       END IF

       IF (ANY((fdir == 5._dp) .AND. agf_k > 0._dp )) THEN
          CALL finish('hydrology_slm_invariants', &
               'baseflow over interior darinage basin')
       END IF

       ! prepare HD grid information for locean = .TRUE., i.e. presumably for oasis3.
       ! grid cell centers
       IF (locean) THEN
         DO il = 1, grid_hd%nlon
           lon_hd(il,:) = grid_hd%origin_lon + REAL(il,dp)*grid_hd%resolution - 0.5_dp*grid_hd%resolution
         END DO
         DO jb = 1, grid_hd%nlat
           lat_hd(:,jb) = grid_hd%origin_lat - REAL(jb,dp)*grid_hd%resolution + 0.5_dp*grid_hd%resolution
         END DO
         CALL fill_oclook_caches(olm, fdir)
       ENDIF

       IF (nremap.EQ.1) THEN
         CALL read_remap_matrix
       ENDIF

       CALL create_kasglob_list(alf_n, alf_k, hd_steps_per_day, alf_n_kas)
       CALL create_kasglob_list(agf_n, agf_k, hd_steps_per_day, agf_n_kas)
       CALL create_kasglob_list(arf_n, arf_k, riverflow_steps_per_day, arf_n_kas)
    END IF
  END SUBROUTINE hydrology_slm_invariants

  SUBROUTINE cleanup_hd_slm_invariants
    IF (p_parallel_io) THEN
      CALL cleanup_kasglob_list(alf_n_kas)
      CALL cleanup_kasglob_list(agf_n_kas)
      CALL cleanup_kasglob_list(arf_n_kas)
    END IF
  END SUBROUTINE cleanup_hd_slm_invariants

  SUBROUTINE set_riverflow_timestep
    !*************************************************************************
    !
    ! **** Routine that defines the internal time step used for riverflow
    !      calculations
    !
    ! steps_per_hour = Number of river flow time steps per hour, used for testing 5 Min. vs.
    !
    USE mo_time_control, ONLY: ev_puthd, get_interval_seconds_next

    INTEGER:: steps_per_hour

    ! TODO: if called by JSBACH/ECHAM as subroutine, set
    !       hd_calling_interval = get_interval_seconds_next(ev_puthd)
    !       in the respoective module
    hd_calling_interval = delta_time

#ifdef HD_5MIN
    ! At 5 Min res., the riverflow time step should not be greater than 1h (= 3600 s).
    ! The variable riverflow_timestep defines the number of internal timesteps
    ! per HD model time steps.

    steps_per_hour = 2

    IF (hd_calling_interval <= 1800) THEN
       riverflow_timestep = MAX(1, INT(FLOAT(steps_per_hour)/2. + 0.001))  ! no internal riverflow time steps needed
    ELSE IF (hd_calling_interval <= 3600) THEN
       riverflow_timestep = steps_per_hour       ! internal riverflow time steps may be needed
    ELSE IF (hd_calling_interval <= 7200) THEN
       riverflow_timestep = 2 * steps_per_hour   ! two riverflow time steps per HD time step
    ELSE IF (hd_calling_interval <= 10800) THEN
       riverflow_timestep = 3 * steps_per_hour
    ELSE IF (hd_calling_interval <= 14400) THEN
       riverflow_timestep = 4 * steps_per_hour
    ELSE IF (hd_calling_interval <= 21600) THEN
       riverflow_timestep = 6 * steps_per_hour
    ELSE IF (hd_calling_interval <= 43200) THEN
       riverflow_timestep = 12 * steps_per_hour
    ELSE IF (hd_calling_interval <= 64800) THEN
       riverflow_timestep = 18 * steps_per_hour
    ELSE IF (hd_calling_interval <= 86400) THEN
       riverflow_timestep = 24 * steps_per_hour
    ELSE
       WRITE (message_text,*) 'The hydrology model should be called at least once per day. '&
            //'hd_calling_interval = ', hd_calling_interval, ' seconds.'
       CALL finish ('mo_hydrology: set_riverflow_timestep', message_text)
    END IF
#else
    ! the riverflow time step should not be greater than 6 hours (= 21.600 s).
    ! The variable riverflow_timestep defines the number of internal timesteps
    ! per HD model time steps.

    IF (hd_calling_interval <= 21600) THEN
       riverflow_timestep = 1                  ! no internal riverflow time steps needed
    ELSE IF (hd_calling_interval <= 43200) THEN
       riverflow_timestep = 2                  ! two riverflow time steps per HD time step
    ELSE IF (hd_calling_interval <= 64800) THEN
       riverflow_timestep = 3
    ELSE IF (hd_calling_interval <= 86400) THEN
       riverflow_timestep = 4
    ELSE
       WRITE (message_text,*) 'The hydrology model should be called at least once a day. '&
            //'hd_calling_interval = ', hd_calling_interval, ' seconds.'
       CALL finish ('mo_hydrology: set_riverflow_timestep', message_text)
    END IF
#endif

    hd_steps_per_day = NINT(86400._dp/hd_calling_interval)
    riverflow_steps_per_day = hd_steps_per_day * riverflow_timestep
    div_riverflow_timestep = 1._dp/REAL(riverflow_timestep,dp)

    IF (ldebughd) THEN
       WRITE (message_text,*) 'hd_steps_per_day = ', hd_steps_per_day
       CALL message ('set_riverflow_timestep',message_text)
       WRITE (message_text,*) 'riverflow_steps_per_day = ', riverflow_steps_per_day
       CALL message ('set_riverflow_timestep',message_text)
    END IF

  END SUBROUTINE set_riverflow_timestep

  SUBROUTINE fill_oclook_caches(slm, fdir)

    !*************************************************************************
    !
    ! **** For each HD grid cell the routine finds the closest ocean grid cell
    !      on the ECHAM grid.
    !
    !      First, the Echam grid cell is found, in which the HD grid cell falls.
    !      Then the closest ocean grid cell is seached within a range of ndd
    !      Echam grid cells.

    REAL(dp), INTENT(in) :: slm(nlon, ngl)    ! ocean land mask on ECHAM grid
    INTEGER,  INTENT(in) :: fdir(grid_hd%nlon, grid_hd%nlat)      ! runoff directions on HD grid

    REAL(dp) :: fb, fl                        ! latitudes/longitudes on HD grid
    INTEGER  :: jb, jl                        ! indices of HD grid latitudes/longirudes
    REAL(dp) :: ocscal                        ! ECHAM grid resolution in degree
    INTEGER  :: jlon, jlat                    ! indices of the ECHAM grid cell, in which HD coordinates (jl,jb) fall.
    TYPE(cart_coord_2d) :: coord              ! REAL(dp) indices of HD cells relative to the ECHAM indices. 

    ocscal = fullcirc / REAL(nlon, dp)
    DO jb = 1, grid_hd%nlat
       fb = lat_hd(1,jb)
       coord%latitude = fb
       jlat = dec_monotonic_closest_midpoint(philat, fb, aub=90._dp, alb=-90._dp)
       DO jl = 1, grid_hd%nlon
          fl = MOD(lon_hd(jl,jb) + fullcirc, fullcirc)
          jlon = inc_monotonic_closest_midpoint(philon, fl, alb=0._dp, aub=360._dp)
          coord%longitude = fl
          coord%dir = fdir(jl,jb)
          oclook_cache(jl,jb) = oclook(slm, jlon, jlat, coord)
       END DO
    END DO
  END SUBROUTINE fill_oclook_caches

  SUBROUTINE create_kasglob_list(a_n, a_k, steps_per_day, a_n_kas)
    !
    REAL(dp), INTENT(in) :: a_n(:, :), a_k(:, :)
    INTEGER, INTENT(in) :: steps_per_day
    TYPE(cart_xidx_2d), ALLOCATABLE, INTENT(inout) :: a_n_kas(:)
    REAL(dp) :: divmm
    INTEGER :: size_i, size_j, i, j, num_extents, extent, extlen
    INTEGER :: jl, extelem

    divmm = 1._dp/REAL(steps_per_day,dp)

    size_i = SIZE(a_n, 1)   ! number of longitudes of the HD-grid
    size_j = SIZE(a_n, 2)   ! number of latitudes of the HD-grid
    num_extents = 0
    IF (ALLOCATED(a_n_kas)) THEN
      CALL cleanup_kasglob_list(a_n_kas)
    END IF

    ! Each grid row is devided into pieces with positive reservoir numbers.
    ! The number of these pieces is counted (num_extents).  
    DO j = 1, size_j
      i = 1
      DO WHILE(i <= size_i)
        IF (a_n(i, j) > 0.5_dp) THEN
          num_extents = num_extents + 1
          i = i + 1
          DO WHILE(i <= size_i)
            IF (a_n(i, j) <= 0.5_dp) EXIT
            i = i + 1
          END DO
        ELSE
          i = i + 1
        END IF
      END DO
    END DO

    ! Each stripe with positive reservoir numbers is filled with data:
    !   ilat: y-index of the first cell of the stripe
    !   ilon: x-index of the first cell of the stripe
    !   extlen: number of grid cells in the stripe
    !   amod: amod for each grid cell in the stripe
    !   akdiv: value for each grid cell in the stripe (depending on time step)

    ALLOCATE(a_n_kas(num_extents))
    extent = 0
    DO j = 1, size_j
      i = 1
      DO WHILE(i <= size_i)
        IF (a_n(i, j) > 0.5_dp) THEN
          extent = extent + 1
          extlen = 0
          a_n_kas(extent)%ilat = j
          a_n_kas(extent)%ilon = i
          DO WHILE(i <= size_i)
            IF (a_n(i, j) <= 0.5_dp) EXIT
            i = i + 1
            extlen = extlen + 1
          END DO
          a_n_kas(extent)%extlen = extlen
          ALLOCATE(a_n_kas(extent)%amod(extlen), a_n_kas(extent)%akdiv(extlen))
          DO extelem = 1, extlen
            jl = a_n_kas(extent)%ilon + extelem - 1
!print*, 'vg: ', a_k(jl,j), a_n(jl,j), AINT(a_n(jl,j)) ! Zeile wichtig fuer intel compiler!!
            a_n_kas(extent)%amod(extelem) = a_k(jl,j) * a_n(jl,j) &
                 / AINT(a_n(jl,j))
            a_n_kas(extent)%akdiv(extelem) = 1.0_dp &
                 / (a_n_kas(extent)%amod(extelem) + divmm)
          END DO
        ELSE
          i = i + 1
        END IF
      END DO
    END DO
  END SUBROUTINE create_kasglob_list

  SUBROUTINE cleanup_kasglob_list(a_n_kas)
    TYPE(cart_xidx_2d), ALLOCATABLE, INTENT(inout) :: a_n_kas(:)
    INTEGER :: i, n
    n = SIZE(a_n_kas)
    DO i = 1, n
      DEALLOCATE(a_n_kas(i)%amod, a_n_kas(i)%akdiv)
    END DO
    DEALLOCATE(a_n_kas)
  END SUBROUTINE cleanup_kasglob_list

!******************************************************************************
  SUBROUTINE hd_remap_05degto5min(nlon_src, nlat_src, src_array,         &
                                   nlon_dst, nlat_dst, dst_array)
!******************************************************************************

!     ******* This programs rmaps from the global 0.5 degree grid (org. HD model)
!             to the new 5 Min grid
!             Interpolation is conservative as grid do exactly match.
!
!     *** Programmierung und Entwicklung: Stefan Hagemann
!     *** Version 1.0 -- December 2014

      INTEGER,  INTENT(in)  :: nlon_src, nlat_src           ! source grid dimensions (0.5 deg)
      INTEGER,  INTENT(in)  :: nlon_dst, nlat_dst           ! target grid dimensions (5 Min.)
      REAL(dp), INTENT(in)  :: src_array(nlon_src,nlat_src) ! field on source grid
      REAL(dp), INTENT(out) :: dst_array(nlon_dst,nlat_dst) ! field on target grid
!
      INTEGER :: jl, jb, ilon, ilat
!
      DO jb= 1, nlat_src
      DO jl= 1, nlon_src
        DO ilat=1,6
        DO ilon=1,6
          dst_array( (jl-1)*6+ilon, (jb-1)*6+ilat) = src_array(jl, jb)
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      RETURN
  END SUBROUTINE hd_remap_05degto5min

  SUBROUTINE hd_remap(nlon_src, nlat_src, src_array, nlon_dst, nlat_dst, global_corr, &
                      src_area, dst_area, dst_array)

    ! do SCRIP remapping: first order conservative remapping 

    INTEGER,  INTENT(in)  :: nlon_src, nlat_src           ! source grid dimensions
    INTEGER,  INTENT(in)  :: nlon_dst, nlat_dst           ! target grid dimensions
    REAL(dp), INTENT(in)  :: src_array(nlon_src,nlat_src) ! field on source grid
    REAL(dp), INTENT(in)  :: src_area(nlon_src,nlat_src)  ! grid cell area on source grid
    REAL(dp), INTENT(in)  :: dst_area(nlon_dst,nlat_dst)  ! grid cell area on target grid
    LOGICAL,  INTENT(in)  :: global_corr                  ! switch for global conservation
    REAL(dp), INTENT(out) :: dst_array(nlon_dst,nlat_dst) ! field on target grid

    REAL(dp), ALLOCATABLE :: src(:)         ! source field as one dimensional array
    REAL(dp), ALLOCATABLE :: dst(:)         ! destination field as one dimensional array
    REAL(dp) :: integral_src
    REAL(dp) :: integral_dst

    INTEGER :: n

    ALLOCATE(src(nlon_src*nlat_src))
    ALLOCATE(dst(nlon_dst*nlat_dst))

    src = RESHAPE(src_array,(/nlon_src*nlat_src/))
    dst = 0._dp
    DO n = 1, num_links
       dst(dst_address(n)) = dst(dst_address(n)) + src(src_address(n)) * remap_matrix(1,n)
    END DO

    dst_array = RESHAPE(dst,(/nlon_dst,nlat_dst/))

    ! assure global conservation
    IF (global_corr) THEN
       integral_src = SUM(src_array * src_area)
       integral_dst = SUM(dst_array * dst_area)
       dst_array = dst_array * integral_src/integral_dst
    END IF
    
    DEALLOCATE(src)
    DEALLOCATE(dst)

  END SUBROUTINE hd_remap

  SUBROUTINE read_remap_matrix
    !
    ! read matrix for remapping of data on the ECHAM grid to the HD model grid
    ! The matrix was generated offline with the SCRIP library (cdo gencon).

    TYPE (FILE_INFO)  :: fileinfo
    INTEGER           :: dimid, varid, fileid
    CHARACTER(len=80) :: rmpfile
    LOGICAL           :: lex
    REAL(dp), ALLOCATABLE :: array_dp(:)    !! double precision dummy array

    ! File names

    rmpfile = 'rmp_hd.nc'

    ! Read parameter: Land sea mask, RDF, ...

    INQUIRE (file=rmpfile, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(rmpfile),'>'
      CALL finish ('read_remap_matrix', message_text)
    ENDIF

    fileinfo%opened = .FALSE.
    CALL IO_open (rmpfile, fileinfo, IO_READ)
    WRITE (message_text,*) 'Reading remap matrix from file ', TRIM(rmpfile)
    CALL message('read_remap_matrix', message_text)

    fileID = fileinfo%file_id

    CALL IO_inq_dimid (fileID, 'num_links', dimid)
    CALL IO_inq_dimlen (fileID, dimid, num_links)
    CALL IO_inq_dimid (fileID, 'num_wgts', dimid)
    CALL IO_inq_dimlen (fileID, dimid, num_weights)

    ALLOCATE (src_address(num_links))    !! address indices of the source grid
    ALLOCATE (dst_address(num_links))    !! address indices of the destination grid
    ALLOCATE (remap_matrix(num_weights,num_links))
    ALLOCATE (array_dp(num_links))       !! double precision dummy array

    CALL IO_inq_varid (fileID, 'src_address', varid)
    CALL IO_get_var_double (fileID, varid, array_dp)
    src_address = NINT(array_dp)
    CALL IO_inq_varid (fileID, 'dst_address', varid)
    CALL IO_get_var_double (fileID, varid, array_dp)
    dst_address = NINT(array_dp)
    CALL IO_inq_varid (fileID, 'remap_matrix', varid)
    CALL IO_get_var_double (fileID, varid, remap_matrix)

    CALL IO_close(fileinfo)
    DEALLOCATE (array_dp)
  
  END SUBROUTINE read_remap_matrix

  SUBROUTINE cleanup_remap_matrix

    DEALLOCATE(src_address, dst_address, remap_matrix)

  END SUBROUTINE cleanup_remap_matrix

  FUNCTION oclook(olm, jlon, jlat, coord) RESULT(idx)

    !*************************************************************************
    !
    ! **** Routine that looks for the next ocean gridbox on the ECHAM grid.
    !      The HD grid cell center is represented by the coord structure.

    ! Input

    REAL(dp),            INTENT(in) :: olm(:,:)  ! ocean land mask on ECHAM grid  
    INTEGER,             INTENT(in) :: jlon      ! index of ECHAM grid longitude
    INTEGER,             INTENT(in) :: jlat      ! index of ECHAM grid latitude
    TYPE(cart_coord_2d), INTENT(in) :: coord     ! indices, coordinates and mask
                                                 ! value of the HD grid cell
    ! Result

    TYPE(cart_idx_2d) :: idx                     ! nearest ocean grid cell found

    ! Local parameters

    REAL(dp), PARAMETER :: deg2rad = 1.74532925199432957692e-2  ! Degree to rad: 2pi/360

    REAL(dp) :: dx                           ! distance between ECHAM and HD grid cell centers
    REAL(dp) :: dxmin                        ! distance to the closet ECHAM ocean cell
    INTEGER  :: ndd                          ! maximum search distance (in grid cells)
    INTEGER  :: nlon, nlat                   ! number of ECHAM longitudes, latitudes
    INTEGER  :: ii, i1, i2, j1, j2, i, j     ! longitude, latitude indices
    REAL(dp) :: lon1, lat1, lon2, lat2       ! spherical coordinates   
    REAL(dp) :: x1, y1, z1, x2, y2, z2       ! x,y,z-coordinates   

    !
    nlon = SIZE(olm, 1)   ! number of ECHAM longitudes 
    nlat = SIZE(olm, 2)   ! number of ECHAM latiudes
 
    ! Nothing needs to be done for HD non-ocean cells

    IF (coord%dir /= -1 .AND. coord%dir /= 0) THEN  ! land or internal drainage
       idx%ilon = -1
       idx%ilat = -1
       RETURN
    END IF

    ! Check whether central ECHAM cell is an ocean grid cell

    IF (olm(jlon,jlat) < 0.5_dp) THEN  ! ocean cell
       idx%ilon = jlon
       idx%ilat = jlat
       RETURN
    END IF

    !
    !  Check the surrounding box with a diameter of ndd grid cells
    !
    ndd = INT(nlat/16) ! maximum distance: T31 3, T63 6 grid cells
    
    ! cell indicees for the search
    i1 = jlon - ndd
    i2 = jlon + ndd
    j1 = MAX(1,jlat - ndd)
    j2 = MIN(nlat,jlat + ndd)

    ! initializations
    dxmin = HUGE(dp)
    idx%ilon = -1
    idx%ilat = -1

    lon1 = coord%longitude * deg2rad  ! HD longitude [rad]
    lat1 = coord%latitude  * deg2rad  ! HD latitude [rad]

    ! find the nearest ocean grid cell
    DO ii = i1,i2
       i = MOD(nlon-1 + ii, nlon) + 1
       DO j = j1,j2
          IF (olm(i,j) < 0.5_dp) THEN  ! ocean cell

             lon2 = philon(i) * deg2rad     ! ECHAM longitude [rad] 
             lat2 = philat(j) * deg2rad     ! ECHAM latitude [rad]
             !
             ! Transformation to x,y,z-coordinates
             !
             x1 = cos(lat1)*cos(lon1)
             y1 = cos(lat1)*sin(lon1)
             z1 = sin(lat1)

             x2 = cos(lat2)*cos(lon2)
             y2 = cos(lat2)*sin(lon2)
             z2 = sin(lat2)
             !
             ! Calculation of the distance
             ! 
             ! direct distance:
             dx = SQRT((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

             ! distance along the surface:
             dx = 2*ASIN(dx/2)

             IF (dx < dxmin) THEN
                dxmin = dx
                idx%ilon = i
                idx%ilat = j
             ENDIF
          END IF
       END DO
    END DO

    ! no ocean grid cell found within the search area?
    IF (dxmin == HUGE(dp)) THEN
       WRITE (message_text,*) 'no ocean cell found for ', coord%longitude, coord%latitude
       CALL message('oclook', message_text)
       WRITE (message_text,*) 'search radius ndd =', ndd, ' needs to be increased.'
       CALL finish('oclook', message_text)
    END IF
  END FUNCTION oclook

  !***********************************************************************
  ! Open outflow timeseries file
  !***********************************************************************
  SUBROUTINE hd_open_timeseries (isolog)

    USE mo_filename, ONLY: find_next_free_unit

    INTEGER,      INTENT(in)  :: isolog

    CHARACTER(80)             :: filename

    IF (p_parallel_io) THEN
       isolog_unit = find_next_free_unit (80, 100)
       WRITE (filename,'(a,i2.2,a)') 'hd_outflow_', isolog, '.log'
       OPEN (isolog_unit, file=TRIM(filename))
    ENDIF

  END SUBROUTINE hd_open_timeseries

  !***********************************************************************
  ! Open outflow timeseries file
  !***********************************************************************
  SUBROUTINE hd_close_timeseries

    IF (p_parallel_io) THEN
       CLOSE (isolog_unit)
    ENDIF

  END SUBROUTINE hd_close_timeseries

  !***********************************************************************
  ! Write outflow timeseries of specific grid cells (compare subroutine 
  ! hydrology_diags)
  !***********************************************************************
  SUBROUTINE hd_write_timeseries (f1, f2)

    USE mo_time_control,  ONLY: current_date, get_date_components

    REAL(dp),      INTENT(in) :: f1, f2
    INTEGER                   :: yr, mo, dy, hr, mn, se

    IF (p_parallel_io) THEN
       CALL get_date_components(current_date, yr, mo, dy, hr, mn, se)
       WRITE(isolog_unit,'(i6.4,i2.2,i2.2,a,i2.2,a,i2.2,a,i2.2,f14.4,f14.4)') &
            yr, mo, dy, ' ', hr, ':', mn, ':', se, F1, F2
    ENDIF

  END SUBROUTINE hd_write_timeseries

  !***********************************************************************
  ! Lateral Routing of discharge 
  !    Routine 1: as in original model at every riverflow time step
  !***********************************************************************
  SUBROUTINE routing(overlandplusbaseflow, riverflow, finfl)
  !
  !        *** Conduct routing in a subroutine for easier exchange of methods
  !
  !  overlandplusbaseflow = local input data array for time step istep
  ! riverflow = local output data array for time step istep
  ! finfl = Inflow data array for each gridbox for time step istep+1
  !  fdir = River direction array
  !
  ! **** Indices
  !
  !    jl = Longitudinal index
  !    jb = Latitudinal index
  !    il = relative change in longitude for the routing
  !    ib = relative change in latitude for the routing
  ! jlnew = jl+il
  ! jbnew = jb+ib

    REAL(dp), INTENT(in)  :: riverflow(:,:)
    REAL(dp), INTENT(in)  :: overlandplusbaseflow(:,:)
    REAL(dp), INTENT(out)  :: finfl(:,:)
    INTEGER :: jl, il, jlnew, jb, ib, jbnew, idir
!  
    ! re-initialize finfl
    finfl(:,:) = 0.0_dp

    ! ---------
    !  routing of outflow to finfl ==> new inflow per gridbox
    ! ---------

    DO jb = 1, grid_hd%nlat
       DO jl = 1, grid_hd%nlon
         idir = fdir(jl, jb) ! from HD parameter input file
         IF (idir > 0) THEN  ! internal land

            ! il, ib = relative direction coordinates [-1,0,1]

            ib = 1 - (idir - 1)/3
            il = MOD(idir - 1, 3) - 1

            jlnew = MOD(jl + il - 1 + grid_hd%nlon, grid_hd%nlon) + 1
            jbnew = jb + ib

          ELSE                ! ocean and coast
            jlnew = jl
            jbnew = jb
          END IF

          ! inflow of the new grid cell:  inflow from other cells 
          !                             + outflow from drainage and overlandflow (overlandplusbaseflow)
          !                             + outflow from riverflow calculations (riverflow)

          finfl(jlnew,jbnew) = finfl(jlnew,jbnew) + overlandplusbaseflow(jl,jb)*div_riverflow_timestep + riverflow(jl,jb)

        ENDDO
     ENDDO

  END SUBROUTINE routing

  !***********************************************************************
  ! Lateral Routing of discharge 
  !    Routine 2: routing via index arrays that are read in from HD parameter file
  !***********************************************************************
  SUBROUTINE routing_via_index(overlandplusbaseflow, riverflow, finfl)
  !
  !        *** Conduct routing in a subroutine for easier exchange of methods
  !
  ! overlandplusbaseflow = local input data array for time step istep
  ! riverflow = local input data array for time step istep
  ! finfl = Inflow data array for each gridbox for time step istep+1
  !  fdir = River direction array

    REAL(dp), INTENT(in)  :: riverflow(:,:)
    REAL(dp), INTENT(in)  :: overlandplusbaseflow(:,:)
    REAL(dp), INTENT(out)  :: finfl(:,:)
    INTEGER :: jl, jb, jlnew, jbnew
!  
    ! re-initialize finfl
    finfl(:,:) = 0.0_dp

    ! ---------
    !  routing of outflow to finfl ==> new inflow per gridbox
    ! ---------

    DO jb = 1, grid_hd%nlat
       DO jl = 1, grid_hd%nlon
         jlnew = filnew(jl, jb)
         jbnew = fibnew(jl, jb)

          ! inflow of the new grid cell:  inflow from other cells 
          !                             + outflow from drainage and overlandflow (overlandplusbaseflow)
          !                             + outflow from riverflow calculations (riverflow)

          finfl(jlnew,jbnew) = finfl(jlnew,jbnew) + overlandplusbaseflow(jl,jb)*div_riverflow_timestep + riverflow(jl,jb)
        ENDDO
     ENDDO

  END SUBROUTINE routing_via_index

  !***********************************************************************
  SUBROUTINE redistribute_sinks(icpl_sinks, disch, hd_outflow) 
  !***********************************************************************
  !
  ! The water in sink cells (internal drainage) can be redistributed before transport into the ocean.
  ! This is necessary in global long-term climate applications where no water should get lost
  ! icpl_sinks =  0   No redistribution
  !               1   Distribute relatively equal to all mouth boxes with discharge > Zero, 
  !                   i.e. the added perentage value is the same at all mouth boxes, so that boxes with 
  !                   small discharge only get a small extra amount, and larger inflows get larger addons.
  !               2   Distribute equally to all ocean boxes.

    INTEGER,  INTENT(in)  :: icpl_sinks  
    REAL(dp), INTENT(in)  :: disch(:,:)
    REAL(dp), INTENT(out) :: hd_outflow(:,:)
    REAL(dp)  :: internal_drainage, discharge_at_mouths, fraction_internal
    INTEGER   :: nocean
    REAL(dp), PARAMETER :: zeps = 1.E-10

    ! Distributing the internal drainage water to
    !   1) all Ocean Inflow Points with applying a weight to treat arid and humid regions differently
    !   2) all ocean points

    IF (icpl_sinks.EQ.1) THEN
!
!     *** Discharge sum over mouth grid boxes 
      discharge_at_mouths = SUM(disch(:,:), fdir(:,:) .EQ. 0)
      IF (discharge_at_mouths.LT.zeps) THEN
        WRITE(message_text,*) 'ERROR: no inflow at any HD river mouth found'
        CALL finish('redistribute_sinks', message_text)
      ENDIF
!     *** Discharge sum over sink grid boxes 
      internal_drainage = SUM(disch(:,:), fdir(:,:) .EQ. 5)
      fraction_internal = internal_drainage / discharge_at_mouths 
      WHERE (disch(:,:).GT.zeps)
        WHERE(fdir(:,:) .EQ. 0)
          hd_outflow(:,:) = disch(:,:) * (1._dp + fraction_internal)
        ELSE WHERE(fdir(:,:) .EQ. 5)
          hd_outflow(:,:) = 0._dp
        ELSEWHERE
          hd_outflow(:,:) = disch(:,:)
        END WHERE
      ELSEWHERE
        hd_outflow(:,:) = 0._dp
      END WHERE
    ELSE IF (icpl_sinks.EQ.2) THEN
!     *** Discharge sum over sink grid boxes 
      internal_drainage = SUM(disch(:,:), fdir(:,:) .EQ. 5)
      nocean = COUNT(fdir(:,:) .EQ. 0 .OR. fdir(:,:) .EQ. -1)
      fraction_internal = internal_drainage / float(nocean) 
      WHERE(fdir(:,:) .EQ. 0 .OR. fdir(:,:) .EQ. -1)
        hd_outflow(:,:) = disch(:,:) + fraction_internal
      ELSEWHERE
        hd_outflow(:,:) = 0._dp
      END WHERE      
    ELSE
      hd_outflow(:,:) = disch(:,:)
    ENDIF

  END SUBROUTINE redistribute_sinks

END MODULE mo_hydrology
