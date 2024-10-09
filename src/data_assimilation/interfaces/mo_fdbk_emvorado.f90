! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! 3DVAR/COSMO source module for feedback file interface to COSMO

MODULE mo_fdbk_emvorado

!-------------------------------------------------------------------------------
!
! Description:
!   COSMO interface to write NetCDF feedobs (or feedback) file (FOF, 
!   Common format for FOF in 3DVAR and COSMO), especially for
!   the radar forward OPERATOR emvorado.
!
!   This module contains the following module procedures:
!     - write_report_radar_1
!     - write_report_radar_2
!     - add_data : interface for: add_inte_vala, add_real_vala_1D,
!                                 add_text_vala, add_real_vala_2D
!     - fill_bodybuf_int, fill_bodybuf_real
!
!   This module has been derived from the original "mo_fdbk_cosmo.f90"
!   and contains somewhat optimized clones of the original procedure
!   "write_report", namely "write_report_radar_1" and "write_report_radar_2".
!   These differ from the original in the sense of a perfomance-optimized
!   interface and vectorizable loop structures for radar data.
!
!-------------------------------------------------------------------------------

#ifdef __DACE__
USE mo_kind,          ONLY: sp

USE mo_fdbk,          ONLY: t_fdbk

USE mo_t_netcdf_file, ONLY: nlen

USE mo_netcdf_param   ! provides access to the parameters defined
                      ! in 'netcdf.inc' of the NetCDF package
#endif

IMPLICIT  NONE

#ifdef __DACE__

!================
! public entities
!================
PRIVATE
!-------------------------
! derived type definitions
!-------------------------
PUBLIC :: t_account
PUBLIC :: t_acc_header
PUBLIC :: t_acc_body
PUBLIC :: t_acc_radar
!------------
! subroutines
!------------
!!! PUBLIC :: write_report, write_report_radar_1, write_report_radar_2
!!! The procedure write_report is only for COSMO. The DACE code has its own version,
!!!  so we make it invisible to other modules. Just the radar routines are visible.
PUBLIC :: write_report_radar_1, write_report_radar_2

!Interface block
 INTERFACE add_data
   MODULE PROCEDURE              &
     add_inte_vala,              &
     add_real_vala_1D,           &
     add_real_vala_2D,           &
     add_text_vala
END INTERFACE


TYPE t_acc_header
  INTEGER           ::   i_body            ! body offset in netcdf file
  INTEGER           ::   l_body            ! number of body elements in report
  INTEGER           ::   n_level           ! numbers of level in report
  INTEGER           ::   i_spec            ! spec offset in netcdf file
  INTEGER           ::   l_spec            ! number of spec elements in report
  INTEGER           ::   spec_r_flags      ! bitwise flags for radar sim config
  INTEGER           ::   varno_back        ! 
  INTEGER           ::   data_category     !
  INTEGER           ::   sub_category      !
  INTEGER           ::   center            !
  INTEGER           ::   sub_center        !
  INTEGER           ::   obstype           !
  INTEGER           ::   codetype          !
  INTEGER           ::   ident             !
  CHARACTER(LEN=16) ::   statid            ! in mo_rad the function specifying this is LEN16
  REAL(KIND=sp)     ::   lat               !
  REAL(KIND=sp)     ::   lon               !
  REAL(KIND=sp)     ::   sun_zenit         !
  REAL(KIND=sp)     ::   sat_zenit         !
  REAL(KIND=sp)     ::   vnyquist          !
  INTEGER           ::   time              !
  INTEGER           ::   time_nomi         !
  INTEGER           ::   time_dbase        !
  INTEGER           ::   z_station         !
  INTEGER           ::   z_modsurf         !
  INTEGER           ::   r_state           !
  INTEGER           ::   r_flags           !
  INTEGER           ::   r_check           !
  INTEGER           ::   sta_corr          !
  INTEGER           ::   index_x           !
  INTEGER           ::   index_y           !
  INTEGER           ::   mdlsfc            !
  INTEGER           ::   instype           !
  INTEGER           ::   retrtype          !
  INTEGER           ::   phase             !
  INTEGER           ::   tracking          !
  INTEGER           ::   meas_type         !
  INTEGER           ::   rad_corr          !
  INTEGER           ::   surftype          !
  INTEGER           ::   flg_1dvar         !
  INTEGER           ::   flg_cld           !
  INTEGER           ::   obs_id            !
  INTEGER           ::   source            !
  INTEGER           ::   record            !
  INTEGER           ::   subset            !
  INTEGER           ::   dbkz              !
  INTEGER           ::   index_d           !
  INTEGER           ::   ct_nwc            ! cloud type from NWCSAF
  REAL(KIND=sp)     ::   ch_nwc            ! cloud top height from NWCSAF
END TYPE t_acc_header

TYPE t_acc_body
  INTEGER           ::   varno             !
  REAL(KIND=sp)     ::   obs               !
  REAL(KIND=sp)     ::   bcor              !
  REAL(KIND=sp)     ::   e_o               !
  REAL(KIND=sp)     ::   azimuth           ! GNSS STD: slant azimuth
  REAL(KIND=sp)     ::   level             !
  INTEGER           ::   level_typ         !
  INTEGER           ::   level_sig         !
  INTEGER           ::   state             !
  INTEGER           ::   flags             !
  INTEGER           ::   check             !
  INTEGER           ::   qual              !
  INTEGER           ::   spec_index        !
  REAL(KIND=sp)     ::   plevel            !
  REAL(KIND=sp)     ::   plev_width        ! width of sensitivity function
  REAL(KIND=sp)     ::   accuracy          !
  REAL(KIND=sp)     ::   dlat              !
  REAL(KIND=sp)     ::   dlon              !
  REAL(KIND=sp)     ::   w_qc              !
  REAL(KIND=sp)     ::   obs_par(2)        !
  REAL,    POINTER  ::   veri_data(:) => NULL() !
  !   assume that length of 'veri_data' is <= 1
  !   i.e. simulated obs from at most 1 model run are stored
! REAL(KIND=sp)     ::   veri_data(1) !
  ! for radar output with nveri=1:
  REAL(KIND=sp)     ::   veri_data1     !
END TYPE t_acc_body

! new data structure for accouting of RADAR
TYPE t_acc_radar
  INTEGER           ::   nrange      !
  INTEGER           ::   nbody       !
  REAL(KIND=sp)     ::   azimuth     !
  REAL(KIND=sp)     ::   elevation   !
  REAL(KIND=sp)     ::   range_start !
  REAL(KIND=sp)     ::   drange      !
END TYPE t_acc_radar

TYPE t_account
  INTEGER                         :: len                ! report length
  INTEGER                         :: offset             ! report offset in FOF
  TYPE ( t_acc_header )           :: header             ! report header
  TYPE ( t_acc_body )  , POINTER  :: body(:) => NULL()  ! report body
END TYPE t_account


CONTAINS

!===============================================================================

SUBROUTINE add_inte_vala ( file, ivala, in, start, count, ierror )

!-------------------------------------------------------------------------------
TYPE(t_fdbk),     INTENT(in)    :: file       ! meta file
INTEGER,          INTENT(in)    :: ivala(:)   ! 1D integer buffer to be stored
INTEGER,          INTENT(in)    :: in         ! variable number in meta file
                                              !   to be stored
INTEGER,          INTENT(in)    :: start      ! start position in netcdf file
INTEGER,          INTENT(in)    :: count      ! number of data to be stored
INTEGER,          INTENT(inout) :: ierror     ! netcdf error
!-------------------------------------------------------------------------------

  ierror = nf_put_vara_int ( file% nc% ncid,                &
                             file% nc% vars(in)% varid,     &
                             [start],                       &
                             [count],                       &
                             ivala(1:count)             )

END SUBROUTINE add_inte_vala

!===============================================================================

SUBROUTINE add_real_vala_1D (file, rvala1, in, start, count, ierror )

!-------------------------------------------------------------------------------
TYPE(t_fdbk),     INTENT(in)    :: file       ! meta file
REAL(KIND=sp),    INTENT(in)    :: rvala1(:)  ! 1D real buffer to be stored
INTEGER,          INTENT(in)    :: in         ! variable number in meta file
                                              !   to be stored
INTEGER,          INTENT(in)    :: start      ! start position in netcdf file
INTEGER,          INTENT(in)    :: count      ! number of data to be stored
INTEGER,          INTENT(inout) :: ierror     ! netcdf error
!-------------------------------------------------------------------------------

  ierror = nf_put_vara_real (file% nc% ncid,              &
                             file% nc% vars(in)% varid,   &
                             [start],                     &
                             [count],                     &
                             rvala1(1:count)              )

END SUBROUTINE add_real_vala_1D

!===============================================================================

SUBROUTINE add_real_vala_2D (file,  rvala2, in, nr, start, count, ierror )

!-------------------------------------------------------------------------------
TYPE(t_fdbk),     INTENT(in)    :: file         ! meta file
REAL(KIND=sp),    INTENT(in)    :: rvala2(:,:)  ! 2D real buffer to be stored
INTEGER,          INTENT(in)    :: in        ,& ! variable number in meta file
                                                !   to be stored
                                   nr           ! position data in rvala2
INTEGER,          INTENT(in)    :: start        ! start position in netcdf file
INTEGER,          INTENT(in)    :: count        ! number of data to be stored
INTEGER,          INTENT(inout) :: ierror       ! netcdf error
!-------------------------------------------------------------------------------

  ierror = nf_put_vara_real (file% nc% ncid,              &
                             file% nc% vars(in)% varid,   &
                             (/start, nr            /),   &
                             (/count, size(rvala2,2)/),   &
                             rvala2(:,:)                )

END SUBROUTINE add_real_vala_2D

!===============================================================================

SUBROUTINE add_text_vala (file,  cvala, in,  start, count, ierror )

!-------------------------------------------------------------------------------
TYPE(t_fdbk)    , INTENT(in)    :: file       ! meta file
CHARACTER(LEN=*), INTENT(in)    :: cvala(:)   ! character buffer to be stored
INTEGER         , INTENT(in)    :: in         ! variable number in meta file
                                              !   to be stored
INTEGER,          INTENT(in)    :: start      ! start position in netcdf file
INTEGER,          INTENT(in)    :: count      ! number of data to be stored
INTEGER,          INTENT(inout) :: ierror     ! netcdf error
!-----------------------------------------------------------
INTEGER   ::  cn      ! string length to be stored
!-------------------------------------------------------------------------------
  cn = file% nc% vars(in)% p(1)% dim% len

  ierror = nf_put_vara_text (file% nc% ncid,              &
                             file% nc% vars(in)% varid,   &
                             (/1,    start/),             &
                             (/cn,   count/),             &
                             cvala(1:count)(1:cn)     )

END SUBROUTINE add_text_vala

!===============================================================================



! Version with rep_header(nrep), rep_body(nrep,:) instead of container report(:)
! otherwise similar to write_report above.
!    Optimized on the NEC for the case that number of reports "nrep" is much 
!    smaller than the typical body length of one report.

SUBROUTINE write_report_radar_1 ( file, rep_header, rep_body, spec_radar       &
                        , nrep, dim2_body, dim2_radar, nbody, ihoff, iboff     &
                        , imdi, rmdich, jerr, yerr, yzstatid )

!-------------------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------
  INTEGER,            INTENT (in)    :: nrep       ! number of reports
  INTEGER,            INTENT (in)    :: dim2_body  ! dim2 of rep_body
  INTEGER,            INTENT (in)    :: dim2_radar ! dim2 of spec_radar
  TYPE(t_fdbk),       INTENT (inout) :: file       ! metafile
  TYPE(t_acc_header), INTENT (in)    :: rep_header(:)! report headers to be stored
  TYPE(t_acc_body),   INTENT (in)    :: rep_body(:,:)! report headers to be stored
  TYPE(t_acc_radar),  INTENT (in)    :: spec_radar(:)! report headers to be stored
  INTEGER,            INTENT (in)    :: nbody      ! number of observations for time box
  INTEGER,            INTENT (in)    :: ihoff      ! header offset for time box
  INTEGER,            INTENT (in)    :: iboff      ! body offset for time box
  INTEGER,            INTENT (in)    :: imdi       ! missing data indicator for
                                                   !   integer (2^31-1)
  REAL(KIND=sp)    ,  INTENT (in)    :: rmdich     ! check value for missing
                                                   !   real data (-1.E30)
  INTEGER,            INTENT (inout) :: jerr       ! error status variable
  ! it is assumed here that LEN >= 72 !
  CHARACTER (LEN= *), INTENT (inout) :: yerr       ! error message
  CHARACTER (LEN= *), INTENT (in), OPTIONAL :: yzstatid       
                                                   ! Optional station id valid for all data in rep_header
!---------------------------------------------------------
  INTEGER                       ::  &
    irep    ,& ! loop indices
    jj, iob ,& ! loop indices
    kcase   ,& ! type of variable (int/real, header/body/radar, etc)
!   kd_hdr  ,& ! length of header dimension
!   kd_body ,& ! length of body   dimension
    varid   ,& ! variable id in netcdf file
    start   ,& ! start position in netcdf file
    count   ,& ! number of elements written
    nvdim   ,& ! nr of dimension veri_data field in meta file
    in      ,& ! variable number to be stored in meta file
    nobs    ,& ! number of observations (body length)
    kos     ,& ! index offset obs from current report
    kerr    ,& ! error status variable
    iu         ! helper variable

  CHARACTER (LEN=nlen) :: yvar    ! NetCDF variable name
  CHARACTER (LEN=nlen) :: ydim    ! NetCDF dimension name

  INTEGER ,           ALLOCATABLE :: ivalh(:) ,& ! buffer for int header elem.
                                     ivalb(:)    ! buffer for int body elements
  REAL(KIND=sp)     , ALLOCATABLE :: rvalh(:) ,& ! buffer for real header elem.
                                     rvalb(:)    ! buffer for real body   elem.
  REAL(KIND=sp)     , ALLOCATABLE :: rvala2(:,:) ! buffer for real body   elem.
  CHARACTER(len=100), ALLOCATABLE :: cvala(:)    ! buffer for char header elem.
!-------------------------------------------------------------------------------
 
  jerr = 0
  kerr = 0
  IF (jerr == 0)  ALLOCATE ( ivalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 1
  IF (jerr == 0)  ALLOCATE ( rvalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 2
  IF (jerr == 0)  ALLOCATE ( cvala (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 3
  IF (jerr == 0)  ALLOCATE ( ivalb (MAX(nbody,dim2_radar*nrep)) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 4
  IF (jerr == 0)  ALLOCATE ( rvalb (MAX(nbody,dim2_radar*nrep)) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 5
  IF (jerr /= 0) THEN
    IF (kerr == 1)  WRITE( yerr,'("ERROR in alloc: ivalh (nrep ) ",I8)' ) nrep
    IF (kerr == 2)  WRITE( yerr,'("ERROR in alloc: rvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 3)  WRITE( yerr,'("ERROR in alloc: cvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 4)  WRITE( yerr,'("ERROR in alloc: ivalb (nbody) ",I8)' ) nbody
    IF (kerr == 5)  WRITE( yerr,'("ERROR in alloc: rvalb (nbody) ",I8)' ) nbody
                                                                          RETURN
  ENDIF

  kerr = nf_enddef(file% nc% ncid )
! IF (kerr /= NF_NOERR) THEN
!   jerr = kerr
!   WRITE( yerr,'("ERROR in nf_enddef in write_reports")' )
!                                                                         RETURN
! ENDIF


  write_loop: DO in = 1, file% nc% nvar

    IF (.NOT. file% nc% vars(in)% opt_used)                                CYCLE
    kcase = 0

    yvar(:) = ' '
    yvar  = TRIM(ADJUSTL(file% nc% vars(in)% name))

    ! initialize write-buffers with missing values:
    ivalh(:) = imdi
    rvalh(:) = rmdich *1.1_sp
    cvala(:) = REPEAT(' ',SIZE(cvala,1))
    ivalb(:) = imdi
    rvalb(:) = rmdich *1.1_sp

! 1-dimensional arrays: integer or real, header or body length
! ------------------------------------------------------------

    IF (file% nc% vars(in)% nvdims == 1) THEN

      kcase = 1
      IF (     (file% nc% vars(in)% xtype == NF_FLOAT)                         &
          .OR. (file% nc% vars(in)% xtype == NF_FILL_DOUBLE))  kcase = 2
!     IF (file% nc% vars(in)% xtype == NF_CHAR)  kcase = 5
      ydim  =  file% nc% vars(in)% p(1)% dim% name
      IF (ydim (1:LEN_TRIM( ydim )) == 'd_body') THEN
        kcase = kcase + 2
      ELSEIF (ydim (1:LEN_TRIM( ydim )) == 'd_radar') THEN
        kcase = kcase + 6
      ELSEIF (ydim (1:LEN_TRIM( ydim )) == 'd_veri') THEN
        kcase = 5
      ELSEIF (ydim (1:LEN_TRIM( ydim )) /= 'd_hdr') THEN
        kcase = -6
      ENDIF

      ! integer header elements
      IF (kcase == 1) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('i_body')
          ivalh (1:nrep) = rep_header(1:nrep)% i_body
        CASE ('l_body')
          ivalh (1:nrep) = rep_header(1:nrep)% l_body
        CASE ('n_level')
          ivalh (1:nrep) = rep_header(1:nrep)% n_level
        CASE ('i_spec')
          ivalh (1:nrep) = rep_header(1:nrep)% i_spec
        CASE ('l_spec')
          ivalh (1:nrep) = rep_header(1:nrep)% l_spec
        CASE ('spec_r_flags')
          ivalh (1:nrep) = rep_header(1:nrep)% spec_r_flags
        CASE ('varno_back')
          ivalh (1:nrep) = rep_header(1:nrep)% varno_back
        CASE ('data_category')
          ivalh (1:nrep) = rep_header(1:nrep)% data_category
        CASE ('sub_category')
          ivalh (1:nrep) = rep_header(1:nrep)% sub_category
        CASE ('center')
          ivalh (1:nrep) = rep_header(1:nrep)% center
        CASE ('sub_center')
          ivalh (1:nrep) = rep_header(1:nrep)% sub_center
        CASE ('obstype')
          ivalh (1:nrep) = rep_header(1:nrep)% obstype
        CASE ('codetype')
          ivalh (1:nrep) = rep_header(1:nrep)% codetype
        CASE ('ident')
          ivalh (1:nrep) = rep_header(1:nrep)% ident
        CASE ('time')
          ivalh (1:nrep) = rep_header(1:nrep)% time
        CASE ('time_nomi')
          ivalh (1:nrep) = rep_header(1:nrep)% time_nomi
        CASE ('time_dbase')
          ivalh (1:nrep) = rep_header(1:nrep)% time_dbase
        CASE ('z_station')
          ivalh (1:nrep) = rep_header(1:nrep)% z_station
        CASE ('z_modsurf')
          ivalh (1:nrep) = rep_header(1:nrep)% z_modsurf
        CASE ('r_state')
          ivalh (1:nrep) = rep_header(1:nrep)% r_state
        CASE ('r_flags')
          ivalh (1:nrep) = rep_header(1:nrep)% r_flags
        CASE ('r_check')
          ivalh (1:nrep) = rep_header(1:nrep)% r_check
        CASE ('sta_corr')
          ivalh (1:nrep) = rep_header(1:nrep)% sta_corr
        CASE ('mdlsfc')
          ivalh (1:nrep) = rep_header(1:nrep)% mdlsfc
        CASE ('instype')
          ivalh (1:nrep) = rep_header(1:nrep)% instype
        CASE ('retrtype')
          ivalh (1:nrep) = rep_header(1:nrep)% retrtype
        CASE ('phase')
          ivalh (1:nrep) = rep_header(1:nrep)% phase
        CASE ('tracking')
          ivalh (1:nrep) = rep_header(1:nrep)% tracking
        CASE ('meas_type')
          ivalh (1:nrep) = rep_header(1:nrep)% meas_type
        CASE ('rad_corr')
          ivalh (1:nrep) = rep_header(1:nrep)% rad_corr
        CASE ('surftype')
          ivalh (1:nrep) = rep_header(1:nrep)% surftype
        CASE ('flg_1dvar')
          ivalh (1:nrep) = rep_header(1:nrep)% flg_1dvar
        CASE ('flg_cld')
          ivalh (1:nrep) = rep_header(1:nrep)% flg_cld
        CASE ('index_x')
          ivalh (1:nrep) = rep_header(1:nrep)% index_x
        CASE ('index_y')
          ivalh (1:nrep) = rep_header(1:nrep)% index_y
        CASE ('obs_id')
          ivalh (1:nrep) = rep_header(1:nrep)% obs_id
        CASE ('source')
          ivalh (1:nrep) = rep_header(1:nrep)% source
        CASE ('record')
          ivalh (1:nrep) = rep_header(1:nrep)% record
        CASE ('subset')
          ivalh (1:nrep) = rep_header(1:nrep)% subset
        CASE ('dbkz')
          ivalh (1:nrep) = rep_header(1:nrep)% dbkz
        CASE ('index_d')
          ivalh (1:nrep) = rep_header(1:nrep)% index_d
        CASE DEFAULT
          ivalh (1:nrep) = file% nc% vars(in)% invalid
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (ivalh == imdi)  ivalh = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalh, in, start, count, jerr )
!       =============

      ! real header elements
      ELSEIF (kcase == 2) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('lat')
          rvalh(1:nrep) = rep_header(1:nrep)% lat
        CASE ('lon')
          rvalh(1:nrep) = rep_header(1:nrep)% lon
        CASE ('sun_zenit')
          rvalh(1:nrep) = rep_header(1:nrep)% sun_zenit
        CASE ('vnyquist')
          rvalh(1:nrep) = rep_header(1:nrep)% vnyquist
        CASE DEFAULT
          rvalh(1:nrep) = REAL(file% nc% vars(in)% rinvalid, sp)
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (rvalh < rmdich)  rvalh = REAL(file% nc% vars(in)% rinvalid, sp)
        CALL add_data ( file, rvalh, in, start, count, jerr )
!       =============

      ! integer body elements
      ELSEIF (kcase == 3) THEN
        kos = 0
        DO irep = 1, nrep
          nobs = rep_header(irep)% l_body
          SELECT CASE (yvar (1:LEN_TRIM( yvar )))
          CASE ('varno')
            ivalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% varno
          CASE ('level_typ')
            ivalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% level_typ
          CASE ('level_sig')
            ivalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% level_sig
          CASE ('state')
            ivalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% state
          CASE ('flags')
            ivalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% flags
          CASE ('check')
            ivalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% check
          CASE ('qual')
            ivalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% qual
          CASE ('spec_index')
            ivalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% spec_index
          CASE ('tovs_flag')
            ivalb (kos+1:kos+nobs) = file% nc% vars(in)% invalid
          CASE DEFAULT
            ivalb (kos+1:kos+nobs) = file% nc% vars(in)% invalid
            kcase = -kcase
          END SELECT
          kos = kos + nobs
        ENDDO
        start = iboff
        count = nbody
        WHERE (ivalb == imdi)  ivalb = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalb, in, start, count, jerr )
!       =============

      ! real body elements
      ELSEIF (kcase == 4) THEN
        kos = 0
        DO irep = 1, nrep
          nobs = rep_header(irep)% l_body
!         SELECT CASE (file% nc% vars(in)% name)
          SELECT CASE (yvar (1:LEN_TRIM( yvar )))
          CASE ('obs')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% obs
          CASE ('bcor')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% bcor
          CASE ('e_o')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% e_o
          CASE ('level')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% level
          CASE ('plevel')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% plevel
          CASE ('accuracy')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% accuracy
          CASE ('dlat')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% dlat
          CASE ('dlon')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% dlon
          CASE ('w_qc')
            rvalb (kos+1:kos+nobs) = rep_body(1:nobs,irep)% w_qc
          CASE DEFAULT
            rvalb (kos+1:kos+nobs) = REAL(file% nc% vars(in)% rinvalid, sp)
            kcase = -kcase
          END SELECT
          kos = kos + nobs
        ENDDO
        start  = iboff
        count  = nbody
        WHERE (rvalb < rmdich)  rvalb = REAL(file% nc% vars(in)% rinvalid, sp)
        CALL add_data ( file, rvalb, in, start, count, jerr )
!       =============

      ! integer ray data
      ELSEIF (kcase == 7) THEN
        start  = rep_header(1)% i_spec
        count  = nrep * dim2_radar
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('radar_nrange')
          ivalb (1:count) = spec_radar(1:count)% nrange
        CASE ('radar_nbody')
          ivalb (1:count) = spec_radar(1:count)% nbody
        CASE default
          ivalb (1:count) = file% nc% vars(in)% invalid
          kcase = -kcase
        END SELECT
        CALL add_data ( file, ivalb, in, start, count, jerr )

      ! real ray data
      ELSEIF (kcase == 8) THEN
        start  = rep_header(1)% i_spec
        count  = nrep * dim2_radar
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('radar_azimuth')
          rvalb (1:count) = spec_radar(1:count)% azimuth
        CASE ('radar_elevation')
          rvalb (1:count) = spec_radar(1:count)% elevation
        CASE ('radar_drange')
          rvalb (1:count) = spec_radar(1:count)% drange
        CASE ('radar_range_start')
          rvalb (1:count) = spec_radar(1:count)% range_start
        CASE default
          rvalb (1:count) = REAL(file% nc% vars(in)% rinvalid, sp)
          kcase = -kcase
        END SELECT
        CALL add_data ( file, rvalb, in, start, count, jerr )

      ! integer 1-d verification meta data elements
      ELSEIF (kcase == 5) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('veri_run_type')
        CASE ('veri_run_class')
        CASE ('veri_forecast_time')
        CASE ('veri_ens_member')
        CASE ('veri_exp_id')
        CASE ('veri_operator_flag')
        CASE DEFAULT
          kcase = -kcase
        END SELECT
      ENDIF

! multi-dimensional arrays: 'statid', 'veri_data'
! -----------------------------------------------

    ! 'statid': the only character string (internally: 2-dim array)

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'statid') THEN

      kcase = 11
      IF (PRESENT(yzstatid)) THEN
        DO jj=1,LEN_TRIM(yzstatid)
          DO irep = 1, nrep
            cvala(irep)(jj:jj) =   yzstatid(jj:jj)
          ENDDO
        END DO
      ELSE
        DO irep = 1, nrep
          cvala(irep) =               rep_header(irep)% statid               &
                         (1:LEN_TRIM( rep_header(irep)% statid))
        ENDDO
      END IF
      start = ihoff
      count = nrep
      CALL add_data ( file, cvala, in, start, count, jerr )
!     =============

    ! 'veri_data': the only truly 2-dimensional array

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'veri_data') THEN

      kcase = 12
!     file% nc% vars(in)% p(2)% dim% len = file% n_veri
      ALLOCATE( rvala2 (nbody, file% n_veri), STAT = jerr )
      IF (jerr /= 0) THEN
        WRITE( yerr,'("ERROR in memory alloc: rvala2 (nbody, file% n_veri) "   &
                    &,2I7)' )  nbody, file% n_veri
                                                                          RETURN
      ENDIF
      rvala2(:,:) = rmdich *1.1_sp
      start = iboff
      count = nbody
!!$      DO jj = 1, file% n_veri
        jj = 1
        kos = 0
!CDIR NODEP
        DO irep = 1, nrep
          nobs = rep_header(irep)% l_body
          DO iob = 1, nobs
            rvala2 (kos+iob,jj) = rep_body(iob,irep)% veri_data1
          ENDDO
          kos = kos + nobs
        ENDDO
        WHERE (rvala2(:,jj)< rmdich) rvala2(:,jj) = REAL(file% nc% vars(in)% rinvalid, sp)
        CALL add_data ( file, rvala2, in, jj, start, count, jerr )
!       =============

      DEALLOCATE(rvala2)

    ! other multi-dimensional arrays (or scalars)
    ELSE
      kcase = 10
      SELECT CASE (yvar (1:LEN_TRIM( yvar )))
      CASE ('veri_model')
      CASE ('veri_initial_date')
      CASE ('veri_resolution')
      CASE ('veri_domain_size')
      CASE ('veri_description')
      CASE DEFAULT
        kcase = -kcase
      END SELECT
    ENDIF

! caution message for unknown variables, dimensions etc.
! ------------------------------------------------------

    IF (kcase <= 0) THEN
      PRINT '("CAUTION writing feedobs file:",I3," unknown variable: ",A )'    &
             , kcase, yvar(1:LEN_TRIM( yvar ))
    ENDIF

!   IF ((kcase == 0) .OR. (kcase == 9)) THEN
!     jerr = 27
!     WRITE( yerr,'("ERROR in write_report:",I2,": unknown var: ",A)' )        &
!            kcase, file% nc% vars(in)% name
!                                                                         RETURN
!   ENDIF
  ENDDO write_loop

  DEALLOCATE(ivalh)
  DEALLOCATE(rvalh)
  DEALLOCATE(cvala)
  DEALLOCATE(ivalb)
  DEALLOCATE(rvalb)

END SUBROUTINE write_report_radar_1

!===============================================================================

!!$ UB>> WORKZONE: implement variables for d_radar as above in write_report_radar_1
!!$                (kcase = 7 and 8)

! .. Version similar to write_report_radar_1, but optimized on the NEC for
!    the case that number of reports "nrep" is much larger than the 
!    typical body length of one report.

SUBROUTINE write_report_radar_2 ( file, rep_header, rep_body, rep_offset, rep_len &
                        , nrep, dim2_body, nbody, ihoff, iboff &
                        , imdi, rmdich, jerr, yerr )
 
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------
  INTEGER,            INTENT (in)    :: nrep       ! number of reports
  INTEGER,            INTENT (in)    :: dim2_body  ! dim2 of rep_body
  TYPE(t_fdbk),       INTENT (inout) :: file       ! metafile
  TYPE(t_acc_header), INTENT (in)    :: rep_header(nrep)! report headers to be stored
  TYPE(t_acc_body),   INTENT (in)    :: rep_body(nrep,dim2_body)! report headers to be stored
  INTEGER,            INTENT (in)    :: rep_offset(nrep)! offset of single reports
  INTEGER,            INTENT (in)    :: rep_len(nrep)   ! len of single report bodys
  INTEGER,            INTENT (in)    :: nbody      ! number of observations for time box
  INTEGER,            INTENT (in)    :: ihoff      ! header offset for time box
  INTEGER,            INTENT (in)    :: iboff      ! body offset for time box
  INTEGER,            INTENT (in)    :: imdi       ! missing data indicator for
                                                   !   integer (2^31-1)
  REAL(KIND=sp)    ,  INTENT (in)    :: rmdich     ! check value for missing
                                                   !   real data (-1.E30)
  INTEGER,            INTENT (inout) :: jerr       ! error status variable
  ! it is assumed here that LEN >= 72 !
  CHARACTER (LEN= *), INTENT (inout) :: yerr       ! error message

!---------------------------------------------------------
  INTEGER                       ::  &
    irep    ,& ! loop indices
    jj, iob, iobu ,& ! loop indices
    kcase   ,& ! type of variable (int/real, header/body, etc)
!   kd_hdr  ,& ! length of header dimension
!   kd_body ,& ! length of body   dimension
    varid   ,& ! variable id in netcdf file
    start   ,& ! start position in netcdf file
    count   ,& ! number of elements written
    nvdim   ,& ! nr of dimension veri_data field in meta file
    in      ,& ! variable number to be stored in meta file
    nobs    ,& ! number of observations (body length)
    kos     ,& ! index offset obs from current report
    kerr    ,& ! error status variable
    iu         ! helper variable

  CHARACTER (LEN=nlen) :: yvar    ! NetCDF variable name
  CHARACTER (LEN=nlen) :: ydim    ! NetCDF dimension name

  INTEGER ,           ALLOCATABLE :: ivalh(:) ,& ! buffer for int header elem.
                                     ivalb(:)    ! buffer for int body elements
  REAL(KIND=sp)     , ALLOCATABLE :: rvalh(:) ,& ! buffer for real header elem.
                                     rvalb(:)    ! buffer for real body   elem.
  REAL(KIND=sp)     , ALLOCATABLE :: rvala2(:,:) ! buffer for real body   elem.
  CHARACTER(len=100), ALLOCATABLE :: cvala(:)    ! buffer for char header elem.
  INTEGER           :: ifillbuf(nrep,dim2_body)
  REAL(KIND=sp)     :: rfillbuf(nrep,dim2_body)
!-------------------------------------------------------------------------------
 
  jerr = 0
  kerr = 0
  IF (jerr == 0)  ALLOCATE ( ivalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 1
  IF (jerr == 0)  ALLOCATE ( rvalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 2
  IF (jerr == 0)  ALLOCATE ( cvala (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 3
  IF (jerr == 0)  ALLOCATE ( ivalb (nbody) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 4
  IF (jerr == 0)  ALLOCATE ( rvalb (nbody) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 5
  IF (jerr /= 0) THEN
    IF (kerr == 1)  WRITE( yerr,'("ERROR in alloc: ivalh (nrep ) ",I8)' ) nrep
    IF (kerr == 2)  WRITE( yerr,'("ERROR in alloc: rvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 3)  WRITE( yerr,'("ERROR in alloc: cvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 4)  WRITE( yerr,'("ERROR in alloc: ivalb (nbody) ",I8)' ) nbody
    IF (kerr == 5)  WRITE( yerr,'("ERROR in alloc: rvalb (nbody) ",I8)' ) nbody
                                                                          RETURN
  ENDIF

  kerr = nf_enddef(file% nc% ncid )
! IF (kerr /= NF_NOERR) THEN
!   jerr = kerr
!   WRITE( yerr,'("ERROR in nf_enddef in write_reports")' )
!                                                                         RETURN
! ENDIF

  
  write_loop: DO in = 1, file% nc% nvar

    IF (.NOT. file% nc% vars(in)% opt_used)                                CYCLE
    kcase = 0

    yvar(:) = ' '
    yvar  = TRIM(ADJUSTL(file% nc% vars(in)% name))

    ! initialize write-buffers with missing values:
    ivalh(:) = imdi
    rvalh(:) = rmdich *1.1_sp
    cvala(:) = REPEAT(' ',SIZE(cvala,1))
    ivalb(:) = imdi
    rvalb(:) = rmdich *1.1_sp

! 1-dimensional arrays: integer or real, header or body length
! ------------------------------------------------------------

    IF (file% nc% vars(in)% nvdims == 1) THEN

      kcase = 1
      IF (     (file% nc% vars(in)% xtype == NF_FLOAT)                         &
          .OR. (file% nc% vars(in)% xtype == NF_FILL_DOUBLE))  kcase = 2
!     IF (file% nc% vars(in)% xtype == NF_CHAR)  kcase = 5
      ydim  =  file% nc% vars(in)% p(1)% dim% name
      IF (ydim (1:LEN_TRIM( ydim )) == 'd_body') THEN
        kcase = kcase + 2
      ELSEIF (ydim (1:LEN_TRIM( ydim )) == 'd_veri') THEN
        kcase = 5
      ELSEIF (ydim (1:LEN_TRIM( ydim )) /= 'd_hdr') THEN
        kcase = -6
      ENDIF

      ! integer header elements
      IF (kcase == 1) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('i_body')
          ivalh (1:nrep) = rep_header(1:nrep)% i_body
        CASE ('l_body')
          ivalh (1:nrep) = rep_header(1:nrep)% l_body
        CASE ('n_level')
          ivalh (1:nrep) = rep_header(1:nrep)% n_level
        CASE ('i_spec')
          ivalh (1:nrep) = rep_header(1:nrep)% i_spec
        CASE ('l_spec')
          ivalh (1:nrep) = rep_header(1:nrep)% l_spec
        CASE ('spec_r_flags')
          ivalh (1:nrep) = rep_header(1:nrep)% spec_r_flags
        CASE ('varno_back')
          ivalh (1:nrep) = rep_header(1:nrep)% varno_back
        CASE ('data_category')
          ivalh (1:nrep) = rep_header(1:nrep)% data_category
        CASE ('sub_category')
          ivalh (1:nrep) = rep_header(1:nrep)% sub_category
        CASE ('center')
          ivalh (1:nrep) = rep_header(1:nrep)% center
        CASE ('sub_center')
          ivalh (1:nrep) = rep_header(1:nrep)% sub_center
        CASE ('obstype')
          ivalh (1:nrep) = rep_header(1:nrep)% obstype
        CASE ('codetype')
          ivalh (1:nrep) = rep_header(1:nrep)% codetype
        CASE ('ident')
          ivalh (1:nrep) = rep_header(1:nrep)% ident
        CASE ('time')
          ivalh (1:nrep) = rep_header(1:nrep)% time
        CASE ('time_nomi')
          ivalh (1:nrep) = rep_header(1:nrep)% time_nomi
        CASE ('time_dbase')
          ivalh (1:nrep) = rep_header(1:nrep)% time_dbase
        CASE ('z_station')
          ivalh (1:nrep) = rep_header(1:nrep)% z_station
        CASE ('z_modsurf')
          ivalh (1:nrep) = rep_header(1:nrep)% z_modsurf
        CASE ('r_state')
          ivalh (1:nrep) = rep_header(1:nrep)% r_state
        CASE ('r_flags')
          ivalh (1:nrep) = rep_header(1:nrep)% r_flags
        CASE ('r_check')
          ivalh (1:nrep) = rep_header(1:nrep)% r_check
        CASE ('sta_corr')
          ivalh (1:nrep) = rep_header(1:nrep)% sta_corr
        CASE ('mdlsfc')
          ivalh (1:nrep) = rep_header(1:nrep)% mdlsfc
        CASE ('instype')
          ivalh (1:nrep) = rep_header(1:nrep)% instype
        CASE ('retrtype')
          ivalh (1:nrep) = rep_header(1:nrep)% retrtype
        CASE ('phase')
          ivalh (1:nrep) = rep_header(1:nrep)% phase
        CASE ('tracking')
          ivalh (1:nrep) = rep_header(1:nrep)% tracking
        CASE ('meas_type')
          ivalh (1:nrep) = rep_header(1:nrep)% meas_type
        CASE ('rad_corr')
          ivalh (1:nrep) = rep_header(1:nrep)% rad_corr
        CASE ('surftype')
          ivalh (1:nrep) = rep_header(1:nrep)% surftype
        CASE ('flg_1dvar')
          ivalh (1:nrep) = rep_header(1:nrep)% flg_1dvar
        CASE ('flg_cld')
          ivalh (1:nrep) = rep_header(1:nrep)% flg_cld
        CASE ('index_x')
          ivalh (1:nrep) = rep_header(1:nrep)% index_x
        CASE ('index_y')
          ivalh (1:nrep) = rep_header(1:nrep)% index_y
        CASE ('obs_id')
          ivalh (1:nrep) = rep_header(1:nrep)% obs_id
        CASE ('source')
          ivalh (1:nrep) = rep_header(1:nrep)% source
        CASE ('record')
          ivalh (1:nrep) = rep_header(1:nrep)% record
        CASE ('subset')
          ivalh (1:nrep) = rep_header(1:nrep)% subset
        CASE ('dbkz')
          ivalh (1:nrep) = rep_header(1:nrep)% dbkz
        CASE ('index_d')
          ivalh (1:nrep) = rep_header(1:nrep)% index_d
        CASE DEFAULT
          ivalh (1:nrep) = file% nc% vars(in)% invalid
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (ivalh == imdi)  ivalh = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalh, in, start, count, jerr )
!       =============

      ! real header elements
      ELSEIF (kcase == 2) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('lat')
          rvalh(1:nrep) = rep_header(1:nrep)% lat
        CASE ('lon')
          rvalh(1:nrep) = rep_header(1:nrep)% lon
        CASE ('sun_zenit')
          rvalh(1:nrep) = rep_header(1:nrep)% sun_zenit
        CASE DEFAULT
          rvalh(1:nrep) = REAL(file% nc% vars(in)% rinvalid, sp)
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (rvalh < rmdich)  rvalh = REAL(file% nc% vars(in)% rinvalid, sp)
        CALL add_data ( file, rvalh, in, start, count, jerr )
!       =============

      ! integer body elements
      ELSEIF (kcase == 3) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('varno')
          ifillbuf = rep_body(:,:)% varno
          CALL fill_bodybuf_int (ifillbuf, ivalb)
        CASE ('level_typ')
          ifillbuf = rep_body(:,:)% level_typ
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('level_sig')
          ifillbuf = rep_body(:,:)% level_sig
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('state')
          ifillbuf = rep_body(:,:)% state
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('flags')
          ifillbuf = rep_body(:,:)% flags
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('check')
          ifillbuf = rep_body(:,:)% check
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('qual')
          ifillbuf = rep_body(:,:)% qual
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('spec_index')
          ifillbuf = rep_body(:,:)% spec_index
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('tovs_flag')
          ifillbuf = file% nc% vars(in)% invalid
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE DEFAULT
          ivalb (:) = file% nc% vars(in)% invalid
          kcase = -kcase
        END SELECT
        start = iboff
        count = nbody
        WHERE (ivalb == imdi)  ivalb = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalb, in, start, count, jerr )
!       =============

      ! real body elements
      ELSEIF (kcase == 4) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('obs')
          rfillbuf = rep_body(:,:)%obs
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('bcor')
          rfillbuf = rep_body(:,:)%bcor
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('e_o')
          rfillbuf = rep_body(:,:)%e_o
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('level')
          rfillbuf = rep_body(:,:)%level
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('plevel')
          rfillbuf = rep_body(:,:)%plevel
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('accuracy')
          rfillbuf = rep_body(:,:)%accuracy
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('dlat')
          rfillbuf = rep_body(:,:)%dlat
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('dlon')
          rfillbuf = rep_body(:,:)%dlon
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('w_qc')
          rfillbuf = rep_body(:,:)%w_qc
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE DEFAULT
          rvalb (:) = REAL(file% nc% vars(in)% rinvalid, sp)
          kcase = -kcase
        END SELECT
        start  = iboff
        count  = nbody
        WHERE (rvalb < rmdich)  rvalb = REAL(file% nc% vars(in)% rinvalid, sp)
        CALL add_data ( file, rvalb, in, start, count, jerr )
!       =============

      ! integer 1-d verification meta data elements
      ELSEIF (kcase == 5) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('veri_run_type')
        CASE ('veri_run_class')
        CASE ('veri_forecast_time')
        CASE ('veri_ens_member')
        CASE ('veri_exp_id')
        CASE ('veri_operator_flag')
        CASE DEFAULT
          kcase = -kcase
        END SELECT
      ENDIF

! multi-dimensional arrays: 'statid', 'veri_data'
! -----------------------------------------------

    ! 'statid': the only character string (internally: 2-dim array)

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'statid') THEN

      kcase = 11
      DO irep = 1, nrep
        cvala(irep) =               rep_header(irep)% statid               &
                       (1:LEN_TRIM( rep_header(irep)% statid))
      ENDDO
      start = ihoff
      count = nrep
      CALL add_data ( file, cvala, in, start, count, jerr )
!     =============

    ! 'veri_data': the only truly 2-dimensional array

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'veri_data') THEN

      kcase = 12
!     file% nc% vars(in)% p(2)% dim% len = file% n_veri
      ALLOCATE( rvala2 (nbody, file% n_veri), STAT = jerr )
      IF (jerr /= 0) THEN
        WRITE( yerr,'("ERROR in memory alloc: rvala2 (nbody, file% n_veri) "   &
                    &,2I7)' )  nbody, file% n_veri
                                                                          RETURN
      ENDIF
      rvala2(:,:) = rmdich *1.1_sp
      start = iboff
      count = nbody
      iobu  = rep_offset(1)
!!$      DO jj = 1, file% n_veri
      jj = 1
        DO iob = 1, dim2_body
!CDIR NODEP
          DO irep = 1, nrep
            IF (iob <= rep_len(irep)) THEN
              iu = rep_offset(irep) - iobu
              rvala2 (iu+iob,jj) = rep_body(irep,iob)% veri_data1
            END IF
          END DO
        END DO
        WHERE (rvala2(:,jj)< rmdich) rvala2(:,jj) = REAL(file% nc% vars(in)% rinvalid, sp)
        CALL add_data ( file, rvala2, in, jj, start, count, jerr )
!       =============

      DEALLOCATE(rvala2)

    ! other multi-dimensional arrays (or scalars)
    ELSE
      kcase = 10
      SELECT CASE (yvar (1:LEN_TRIM( yvar )))
      CASE ('veri_model')
      CASE ('veri_initial_date')
      CASE ('veri_resolution')
      CASE ('veri_domain_size')
      CASE ('veri_description')
      CASE DEFAULT
        kcase = -kcase
      END SELECT
    ENDIF

! caution message for unknown variables, dimensions etc.
! ------------------------------------------------------

    IF (kcase <= 0) THEN
      PRINT '("CAUTION writing feedobs file:",I3," unknown variable: ",A )'    &
             , kcase, yvar(1:LEN_TRIM( yvar ))
    ENDIF

!   IF ((kcase == 0) .OR. (kcase == 9)) THEN
!     jerr = 27
!     WRITE( yerr,'("ERROR in write_report:",I2,": unknown var: ",A)' )        &
!            kcase, file% nc% vars(in)% name
!                                                                         RETURN
!   ENDIF
  ENDDO write_loop

  DEALLOCATE(ivalh)
  DEALLOCATE(rvalh)
  DEALLOCATE(cvala)
  DEALLOCATE(ivalb)
  DEALLOCATE(rvalb)

CONTAINS

  SUBROUTINE fill_bodybuf_int (ibodydata, ibuf)

!!$    INTEGER, INTENT(in)  :: ibodydata (nrep,dim2_body)
!!$    INTEGER, INTENT(inout) :: ibuf (nbody)
    INTEGER, INTENT(in)  :: ibodydata (:,:)
    INTEGER, INTENT(inout) :: ibuf (:)

    INTEGER :: iu, iob, iobu, irep

    ibuf = -9999
    iobu = rep_offset(1)
    DO iob = 1, dim2_body
!CDIR NODEP
      DO irep = 1, nrep
        IF (iob <= rep_len(irep)) THEN
          iu = rep_offset(irep) - iobu
          ibuf (iu+iob) = ibodydata(irep,iob)
        END IF
      END DO
    END DO
  
  END SUBROUTINE fill_bodybuf_int

  SUBROUTINE fill_bodybuf_real (rbodydata, rbuf)

!!$    REAL, INTENT(in)  :: rbodydata (nrep,dim2_body)
!!$    REAL, INTENT(inout) :: rbuf (nbody)
    REAL(KIND=sp),  INTENT(in)  :: rbodydata (:,:)
    REAL(KIND=sp),  INTENT(inout) :: rbuf (:)

    INTEGER :: iu, iob, iobu, irep

    rbuf = -9999.99_sp
    iobu = rep_offset(1)
    DO iob = 1, dim2_body
!CDIR NODEP
      DO irep = 1, nrep
        IF (iob <= rep_len(irep)) THEN
          iu = rep_offset(irep) - iobu
          rbuf (iu+iob) = rbodydata(irep,iob)
        END IF
      END DO
    END DO
  
  END SUBROUTINE fill_bodybuf_real

END SUBROUTINE write_report_radar_2

!-------------------------------------------------------------------------------

#endif

END MODULE mo_fdbk_emvorado
