! Namelist for Single Column Model
!
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

MODULE mo_scm_nml  
!-------------------------------------------------------------------------  
!  
!-------------------------------------------------------------------------
!
!
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, finish
  USE mo_namelist,             ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH
  USE mo_io_units,             ONLY: nnml, nnml_output
  USE mo_master_control,       ONLY: use_restart_namelists
  USE mo_mpi,                  ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att,  ONLY: open_tmpfile, store_and_close_namelist,     &
    &                                open_and_restore_namelist, close_tmpfile
  
  IMPLICIT NONE  

  PRIVATE 

  PUBLIC ::  read_scm_namelist, scm_sfc_temp, scm_sfc_qv, scm_sfc_mom, i_scm_netcdf,   &
  & lon_scm, lat_scm, lscm_read_tke, lscm_read_z0, lscm_icon_ini, lscm_ls_forcing_ini, &
  & lscm_random_noise

  !----------------------------------!
  ! scm_nml namelist variables       !
  !----------------------------------!

  INTEGER  :: scm_sfc_temp ! surface boundary condition for temperature,
                           ! 0: TERRA
                           ! 1: DIRICHLET (t_g)
                           ! 2: SENSIBLE HEAT FLUX (shfl_s)

  INTEGER  :: scm_sfc_qv   ! surface boundary condition for moisture
                           ! 0: TERRA
                           ! 1: DIRICHLET(qv_s)
                           ! 2: LATENT HEAT FLUX(qhfl_s)

  INTEGER  :: scm_sfc_mom  ! surface boundary condition for momentum
                           ! 0: TERRA,1=DIRICHLET(u_s=v_s=0, gz0)
                           ! 2: friction velocity(ustar->tvm)

  INTEGER  :: i_scm_netcdf ! data read from 
                           ! 0: ASCII
                           ! 1: normal netcdf file
                           ! 2: DEPHY unified format
  LOGICAL  :: lscm_read_tke
  LOGICAL  :: lscm_read_z0
  LOGICAL  :: lscm_icon_ini
  LOGICAL  :: lscm_random_noise
  LOGICAL  :: lscm_ls_forcing_ini=.FALSE. ! apply_ls_forcing already called?

  NAMELIST/scm_nml/ scm_sfc_temp, scm_sfc_qv, scm_sfc_mom,                     &
    &               i_scm_netcdf, lscm_read_tke, lscm_read_z0, lscm_icon_ini,  &
    &               lscm_random_noise

  !non-namelist parameters - read from netcdf file
  REAL(wp) :: lon_scm      ! fix latitude and longitude for computations depending
  REAL(wp) :: lat_scm      ! on geographical coordinates on torus geometry


CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE read_scm_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: i_status, funit
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_scm_nml: read_scm_namelist'

!-----------------------------------------------------------------------

    !-----------------------
    ! 1. default settings
    !    * default external data for Lindenberg grid point
    !    * external data at res. R03B07
    !-----------------------

    scm_sfc_temp     = 0       ! prescribed surface sensible heat flux
    scm_sfc_qv       = 0       ! prescribed surface latent heat flux
    scm_sfc_mom      = 0       ! prescribed friction velocity
    lscm_read_tke    = .FALSE. ! read init. tke from netcdf
    lscm_read_z0     = .FALSE. ! read z0 from netcdf file
    lscm_icon_ini    = .FALSE. ! read initial conditions produced by ICON
    lscm_random_noise= .FALSE. ! initialize with random noise - for LEM runs
   !i_scm_netcdf               ! defaults set in mo_nml_crosscheck

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------

    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('scm_nml')
      READ(funit,NML=scm_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('scm_nml', status=i_status)
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, scm_nml)
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! 4. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=scm_nml)
      CALL store_and_close_namelist(funit, 'scm_nml')
    ENDIF

    !-----------------------------------------------------
    ! 5. Write the contents of the namelist to an ASCII file
    !-----------------------------------------------------

    IF(my_process_is_stdio()) WRITE(nnml_output,nml=scm_nml)

    
  END SUBROUTINE read_scm_namelist

!-------------------------------------------------------------------------
END MODULE mo_scm_nml  
