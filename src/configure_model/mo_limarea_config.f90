! Contains the setup of variables related to the lateral boundary
! condition for limited area models
!
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

MODULE mo_limarea_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_io_units,           ONLY: filename_max
  USE mo_util_string,        ONLY: t_keyword_list,                   &
                                   associate_keyword, with_keywords, &
                                   int2string
  USE mo_exception,          ONLY: message, message_text, finish
  USE mtime,                 ONLY: datetime, timedelta, OPERATOR(-)
  USE mo_util_mtime,         ONLY: mtime_utils, FMT_DDDHH, FMT_DDHHMMSS, FMT_HHH
  USE mo_parallel_config,    ONLY: num_prefetch_proc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_latbc_config, latbc_config, configure_latbc, generate_filename
  PUBLIC :: LATBC_TYPE_CONST, LATBC_TYPE_EXT

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_limarea_config'

  ! CONSTANTS (for better readability):
  !
  ! LATBC_TYPE_CONST: constant lateral boundary conditions derived
  ! from the initial conditions
  INTEGER, PARAMETER :: LATBC_TYPE_CONST       = 0

  ! LATBC_TYPE_EXT: time-dependent lateral boundary conditions
  ! provided by an external source (IFS, COSMO-DE or a
  ! coarser-resolution ICON run)
  INTEGER, PARAMETER :: LATBC_TYPE_EXT         = 1



  !>
  !!----------------------------------------------------------------------------
  !! Derived type containing control variables specific to the nonhydrostatic 
  !! atm model
  !------------------------------------------------------------------------
  TYPE t_latbc_config

    ! variables from namelist
    INTEGER                         :: itype_latbc         ! type of limited area boundary nudging
    REAL(wp)                        :: dtime_latbc         ! dt between two consequtive external latbc files
    CHARACTER(LEN=filename_max)     :: latbc_filename      ! prefix of latbc files
    CHARACTER(LEN=MAX_CHAR_LENGTH)  :: latbc_path          ! directory containing external latbc files
    LOGICAL                         :: latbc_contains_qcqi ! latbc data contain qc and qi (=T) or not (=F)
    CHARACTER(LEN=FILENAME_MAX)     :: latbc_boundary_grid ! grid file defining the lateral boundary
    LOGICAL                         :: init_latbc_from_fg  ! take initial lateral boundary conditions from first guess
    LOGICAL                         :: nudge_hydro_pres    ! use hydrostatic pressure for lateral boundary nudging

    ! factor for pressure bias correction of latbc data
    REAL(wp)                        :: fac_latbc_presbiascor 

    ! if LatBC data is unavailable: number of retries
    INTEGER                         :: nretries

    ! if LatBC data is unavailable: idle wait seconds between retries
    INTEGER                         :: retry_wait_sec

    ! dictionary which maps internal variable names onto GRIB2
    ! shortnames or NetCDF var names used in lateral boundary nudging.
    CHARACTER(LEN=filename_max)     :: latbc_varnames_map_file  

    !
    ! settings derived from the namelist parameters above:
    !
    LOGICAL                         :: lsparse_latbc       ! Flag: TRUE if only boundary rows are read.
    TYPE(timedelta), POINTER        :: dtime_latbc_mtime   ! dt between two consequtive external latbc files


  END TYPE t_latbc_config
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  TYPE(t_latbc_config), TARGET :: latbc_config
  !!----------------------------------------------------------------------------


CONTAINS

  SUBROUTINE configure_latbc()
  !--------------------------------------------------------------------------------------
  !  Set up parameters 
  !--------------------------------------------------------------------------------------
    CHARACTER(*), PARAMETER :: routine = &
      "mo_limarea_config::configure_latbc"

    !----------------------------------------------------
    ! Sanity check and Prints
    !----------------------------------------------------

    IF (latbc_config%itype_latbc == LATBC_TYPE_CONST) THEN

       WRITE(message_text,'(a)')'Lateral boundary nudging using constant boundary data from the initial conditions.'
       CALL message(TRIM(routine),message_text)

    ELSE IF (latbc_config%itype_latbc == LATBC_TYPE_EXT) THEN

       WRITE(message_text,'(a)')'Lateral boundary condition using interpolated boundary data.'
       CALL message(TRIM(routine),message_text)

       IF (num_prefetch_proc == 0) THEN
         WRITE(message_text,'(a)') 'Synchronous latBC input has been disabled'
         CALL finish(TRIM(routine),message_text)
       END IF

    ELSE

       WRITE(message_text,'(a,i8)') 'Invalid lateral boundary condition mode:', latbc_config%itype_latbc
       CALL finish(TRIM(routine),message_text)

    END IF

    ! Check whether an mapping file is provided for prefetching boundary data
    ! calls a finish either when the flag is absent
    !
    IF (latbc_config%latbc_varnames_map_file == ' ') THEN
       WRITE(message_text,'(a)') 'no latbc_varnames_map_file provided.'
       CALL message(TRIM(routine),message_text)
    ENDIF

  
  END SUBROUTINE configure_latbc
  !--------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------
  FUNCTION generate_filename(nroot, jlev, latbc_mtime, opt_mtime_begin) RESULT(result_str)
    CHARACTER(MAX_CHAR_LENGTH )                      :: result_str
    INTEGER,          INTENT(IN)                     :: nroot, jlev
    TYPE(datetime),   INTENT(IN)                     :: latbc_mtime
    ! Optional: Start date, which a time span in the filename is related to.
    TYPE(datetime),   INTENT(IN),  POINTER, OPTIONAL :: opt_mtime_begin

    ! Local variables
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER       :: routine = modname//'::generate_filename'
    TYPE (t_keyword_list), POINTER              :: keywords => NULL()
    CHARACTER(MAX_CHAR_LENGTH)                  :: str
    
    WRITE(str,'(i4)')   latbc_mtime%date%year
    CALL associate_keyword("<y>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%month
    CALL associate_keyword("<m>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%day
    CALL associate_keyword("<d>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%hour
    CALL associate_keyword("<h>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%minute
    CALL associate_keyword("<min>",       TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%second !FLOOR(latbc_mtime%time%second)
    CALL associate_keyword("<sec>",       TRIM(str),                        keywords)
      
    CALL associate_keyword("<nroot>",     TRIM(int2string(nroot,'(i1)')),   keywords)
    CALL associate_keyword("<nroot0>",    TRIM(int2string(nroot,'(i2.2)')), keywords)
    CALL associate_keyword("<jlev>",      TRIM(int2string(jlev, '(i2.2)')), keywords)
    CALL associate_keyword("<dom>",       TRIM(int2string(1,'(i2.2)')),     keywords)
    
    IF (PRESENT(opt_mtime_begin)) THEN
      CALL associate_keyword("<ddhhmmss>", &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_DDHHMMSS)), &
        &                    keywords)
      CALL associate_keyword("<dddhh>",    &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_DDDHH)),    &
        &                    keywords)
      CALL associate_keyword("<hhh>",    &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_HHH)),    &
        &                    keywords)
    END IF

    ! replace keywords in latbc_filename
    result_str = TRIM(with_keywords(keywords, TRIM(latbc_config%latbc_filename)))
  END FUNCTION generate_filename
  !--------------------------------------------------------------------------------------


!-----------------------------------------------------------------------
END MODULE mo_limarea_config
