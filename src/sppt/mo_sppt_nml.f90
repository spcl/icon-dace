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

! Namelist for SPPT (Stochastic Pertubation of Physics Tendencies)
!
! these Subroutines are called by control model and construct the
! SPPT scheme

MODULE mo_sppt_nml

  USE mo_kind,                ONLY: wp
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_sppt_config,         ONLY: sppt_config
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_math_constants,      ONLY: deg2rad

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_sppt_namelist

  !-----------------------------------
  ! namelist variables and parameters
  !-----------------------------------

  LOGICAL  :: lsppt     ! forecast with SPPT

  REAL(wp) :: hinc_rn   ! time increment (s) for drawing a new field of random numbers
  REAL(wp) :: dlat_rn   ! random number coarse grid point distance in meridional direction (deg)
  REAL(wp) :: dlon_rn   ! random number coarse grid point distance in zonal direction (deg)
  REAL(wp) :: range_rn  ! max magnitude of random numbers
  REAL(wp) :: stdv_rn   ! standard deviation of the gaussian distribution of random numbers


  NAMELIST/sppt_nml/ lsppt, hinc_rn             , &
       &             dlat_rn, dlon_rn           , &
       &             range_rn                   , &
       &             stdv_rn


  CONTAINS


  ! --------------------------------------------------------------------------
  !
  ! Subroutine: read_sppt_namelist
  !
  ! Purpose: This subroutine
  !          - read the namelist for SPPT (NWP)
  !          - sets default values
  !
  ! --------------------------------------------------------------------------

  SUBROUTINE read_sppt_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: istat, funit
    INTEGER :: iunit
    INTEGER :: jg

    !-----------------------
    ! 1. default settings
    !-----------------------

    lsppt = .FALSE.     ! if TRUE SPPT is switched on

    hinc_rn = 21600_wp    ! time increment for drawing a new field of random numbers        (seconds)

    ! converted from deg to rad when copied to the config state (see below)
    dlat_rn = 1e-1_wp    ! random number coarse grid point distance in meridional direction (deg)
    dlon_rn = 1e-1_wp    ! random number coarse grid point distance in zonal direction      (deg)

    range_rn = 0.8_wp   ! max magnitude of random numbers                                   (-)

    stdv_rn  = 1.0_wp   ! standard deviation of the gaussian distribution of random numbers (-)


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------

    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('sppt_nml')
      READ(funit,NML=sppt_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('sppt_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, sppt_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, sppt_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, sppt_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check (if necessary)
    !----------------------------------------------------

    ! Currently non required keep for consistency with other config namelist routines.

    !-----------------------
    ! 5. Fill the configuration state
    !-----------------------
    DO jg= 1,max_dom
      sppt_config(jg)%lsppt          = lsppt
      sppt_config(jg)%hinc_rn        = hinc_rn
      sppt_config(jg)%dlat_rn        = dlat_rn * deg2rad
      sppt_config(jg)%dlon_rn        = dlon_rn * deg2rad
      sppt_config(jg)%range_rn       = range_rn
      sppt_config(jg)%stdv_rn        = stdv_rn
    ENDDO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------

    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=sppt_nml)
      CALL store_and_close_namelist(funit, 'sppt_nml')
    ENDIF

    !-------------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !-------------------------------------------------------

    IF(my_process_is_stdio()) WRITE(nnml_output,nml=sppt_nml)

  END SUBROUTINE read_sppt_namelist


END MODULE mo_sppt_nml








