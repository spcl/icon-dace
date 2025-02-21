! @brief Namelist for synthetic satellite images
!
! These subroutines are called by read_atmo_namelists and set the switches
! for synthetic satellite images
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

MODULE mo_synsat_nml

  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_synsat_config,       ONLY: config_lsynsat    => lsynsat, &
                                    config_nlev_rttov => nlev_rttov
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_synsat_namelist


  !--------------------------------------!
  ! synsat_nml namelist variables !
  !--------------------------------------!


  LOGICAL :: lsynsat(max_dom)     !< main switch

  INTEGER :: nlev_rttov           !< Number of RTTOV levels

  NAMELIST/synsat_nml/ lsynsat, nlev_rttov


CONTAINS


  !-------------------------------------------------------------------------
  !
  !! Read Namelist for NWP ensemble perturbations. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP ensemble perturbations
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  SUBROUTINE read_synsat_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_synsat_nml: read_synsat_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! Main switch
    lsynsat(:) = .FALSE.

    ! Number of RTTOV levels
    nlev_rttov = 54

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('synsat_nml')
      READ(funit,NML=synsat_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('synsat_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, synsat_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, synsat_nml)             ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, synsat_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    config_lsynsat       = lsynsat
    config_nlev_rttov    = nlev_rttov

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=synsat_nml)                    
      CALL store_and_close_namelist(funit, 'synsat_nml')             
    ENDIF

    !--------------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !--------------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=synsat_nml)


    !--------------------------------------------------------
    ! 8. consistency check
    !--------------------------------------------------------


    IF (ANY(lsynsat)) THEN
      CALL finish(routine, 'Switch "lsynsat": Model has not been configured for RTTOV library.')
    END IF


  END SUBROUTINE read_synsat_namelist

END MODULE mo_synsat_nml
