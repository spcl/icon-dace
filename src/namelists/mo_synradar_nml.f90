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

! Namelist reading for synthetic radar data on the model grid

MODULE mo_synradar_nml

#ifdef HAVE_RADARFWO

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_nml_and_att,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
  USE mo_synradar_config,    ONLY: config_synradar_meta           => synradar_meta          , &               
                                 & config_ydir_mielookup_read     => ydir_mielookup_read    , &
                                 & config_ydir_mielookup_write    => ydir_mielookup_write
  
  USE mo_exception,        ONLY: finish
  USE radar_dbzcalc_params_type, ONLY: t_dbzcalc_params, dbz_namlst_d

#endif

  IMPLICIT NONE
  PUBLIC :: read_synradar_namelist

  ! module name
  CHARACTER(*), PARAMETER :: modname = "mo_synradar_nml"
  
CONTAINS
  !! Read Namelist for I/O.
  !!
  !! This subroutine
  !! - reads the Namelist for I/O
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  SUBROUTINE read_synradar_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN)   :: filename

#ifdef HAVE_RADARFWO
    
    CHARACTER(*), PARAMETER :: routine = modname//":read_synradar_namelist"
    INTEGER                        :: istat, funit
    INTEGER                        :: iunit

    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------


    ! Meta data for reflectivity computations (DBZ, DBZ850, DBZ_CMAX, etc.) on the model grid by using advanced methods
    !  from EMVORADO (Mie-scattering, T-matrix):
    TYPE(t_dbzcalc_params)        :: synradar_meta
    CHARACTER(LEN=filename_max) :: ydir_mielookup_read
    CHARACTER(LEN=filename_max) :: ydir_mielookup_write
    
    NAMELIST/synradar_nml/ synradar_meta, ydir_mielookup_read, ydir_mielookup_write

    !-----------------------
    ! 1. default settings
    !-----------------------

    synradar_meta            = dbz_namlst_d
    synradar_meta%itype_refl = 4      ! default: use the established ICON-method (=4) for dbz-calculations
    ydir_mielookup_read(:)   = ' '    ! only relevant for itype_refl /= 4 (EMVORADO-methods)
    ydir_mielookup_write(:)  = ' '    ! only relevant for itype_refl /= 4 (EMVORADO-methods)
    
    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('synradar_nml')
      READ(funit,NML=synradar_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('synradar_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, synradar_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, synradar_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, synradar_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    SELECT CASE (synradar_meta%itype_refl)
    CASE (1, 3, 4, 5, 6)
      CONTINUE
    CASE default
       CALL finish(routine, "Invalid choice of parameter synradar_meta%itype_refl! Allowed are 1, 3, 4 (default), 5, or 6")
    END SELECT

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    config_synradar_meta           = synradar_meta
    config_ydir_mielookup_read     = ydir_mielookup_read
    config_ydir_mielookup_write    = ydir_mielookup_write
    
    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------

    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=synradar_nml)
      CALL store_and_close_namelist(funit, 'synradar_nml')
    ENDIF

    !-----------------------------------------------------
    ! 6. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------

    IF(my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=synradar_nml)
    END IF

#endif
    
  END SUBROUTINE read_synradar_namelist

END MODULE mo_synradar_nml
