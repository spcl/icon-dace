! Contains the setup of the sleve coordinate
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

MODULE mo_sleve_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_sleve_config       , ONLY: config_min_lay_thckn => min_lay_thckn, &
    &                               config_max_lay_thckn => max_lay_thckn, &
    &                               config_htop_thcknlimit => htop_thcknlimit, &
    &                               config_nshift_above_thcklay => nshift_above_thcklay, &
    &                               config_itype_laydistr  => itype_laydistr,  &
    &                               config_top_height    => top_height   , &
    &                               config_decay_scale_1 => decay_scale_1, &
    &                               config_decay_scale_2 => decay_scale_2, &
    &                               config_decay_exp     => decay_exp    , &
    &                               config_flat_height   => flat_height  , &
    &                               config_stretch_fac   => stretch_fac  , &
    &                               config_lread_smt     => lread_smt 
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC:: read_sleve_namelist

  !----------------------------------------------------------------------------
  ! Namelist variables for the SLEVE coordinate
  !---------------------------------------------------------------------------
  !
  ! a) Parameters specifying the distrubution of the coordinate surfaces
  INTEGER :: itype_laydistr  ! Type of analytical function used for computing the coordinate surface distribution
  REAL(wp):: min_lay_thckn   ! Layer thickness of lowermost level
  REAL(wp):: max_lay_thckn   ! Maximum layer thickness below htop_thcknlimit
  REAL(wp):: htop_thcknlimit ! Height below which the layer thickness must not exceed max_lay_thckn
  REAL(wp):: stretch_fac     ! Factor for stretching/squeezing the model layer distribution
  REAL(wp):: top_height      ! Height of model top

  INTEGER :: nshift_above_thcklay ! Shift above constant-thickness layer for further calculation of layer distribution

  ! b) Parameters for SLEVE definition
  REAL(wp):: decay_scale_1    ! Decay scale for large-scale topography component
  REAL(wp):: decay_scale_2    ! Decay scale for small-scale topography component
  REAL(wp):: decay_exp        ! Exponent for decay function
  REAL(wp):: flat_height      ! Height above which the coordinate surfaces are exactly flat
                              ! (not available in the standard SLEVE definition)

  ! c) Parameter for reading in smoothed topography
  LOGICAL :: lread_smt
 
  NAMELIST /sleve_nml/ min_lay_thckn, max_lay_thckn, htop_thcknlimit, top_height,         &
                       decay_scale_1, decay_scale_2, decay_exp, flat_height, stretch_fac, &
                       lread_smt, itype_laydistr, nshift_above_thcklay

CONTAINS
  !-------------------------------------------------------------------------
  !! Read Namelist for SLEVE coordinate. 
  !!
  !! This subroutine 
  !! - reads the Namelist for SLEVE coordinate
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  SUBROUTINE read_sleve_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    CHARACTER(len=*), PARAMETER :: routine = 'mo_sleve_nml:read_sleve_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! a) Parameters determining the distribution of model layers
    !    (if not read in from a table)
    itype_laydistr  = 1           ! stretched cosine function (2 = third-order polynomial)
    min_lay_thckn   = 50._wp      ! Layer thickness of lowermost layer
    max_lay_thckn   = 25000._wp   ! Maximum layer thickness below htop_thcknlimit
    htop_thcknlimit = 15000._wp   ! Height below which the layer thickness must not exceed max_lay_thckn
    nshift_above_thcklay = 0      ! No layer index shift
    top_height      = 23500._wp   ! Height of model top
    stretch_fac     = 1._wp       ! Scaling factor for stretching/squeezing 
                                  ! the model layer distribution

    ! b) Parameters setting up the decay function of the topographic signal
    decay_scale_1   = 4000._wp    ! Decay scale of large-scale topography component
    decay_scale_2   = 2500._wp    ! Decay scale of small-scale topography component
    decay_exp       = 1.2_wp      ! Exponent for decay function
    flat_height     = 16000._wp   ! Height above which the coordinate surfaces are 
                                      ! flat

    ! c) parameter to switch on/off internal topography smoothing
    lread_smt       = .FALSE.     ! read smoothed topography from file (TRUE/FALSE)

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('sleve_nml')
      READ(funit,NML=sleve_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('sleve_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, sleve_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, sleve_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, sleve_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    ! Comment this call if you want to test other values, but be aware that this is not recommended
    IF (nshift_above_thcklay < 0 .OR. nshift_above_thcklay > 1) THEN
      CALL finish(TRIM(routine),'nshift_above_thcklay should be 0 or 1')
    ENDIF

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    config_min_lay_thckn = min_lay_thckn
    config_max_lay_thckn = max_lay_thckn
    config_htop_thcknlimit = htop_thcknlimit
    config_nshift_above_thcklay = nshift_above_thcklay
    config_itype_laydistr  = itype_laydistr
    config_top_height    = top_height
    config_decay_scale_1 = decay_scale_1
    config_decay_scale_2 = decay_scale_2
    config_decay_exp     = decay_exp
    config_flat_height   = flat_height
    config_stretch_fac   = stretch_fac
    config_lread_smt     = lread_smt

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=sleve_nml)                    
      CALL store_and_close_namelist(funit, 'sleve_nml') 
    ENDIF
    ! 6. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=sleve_nml)

  END SUBROUTINE read_sleve_namelist

END MODULE mo_sleve_nml
