!> Contains structures and methods for JSBACH model config
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_config_class
#ifndef __NO_JSBACH__

  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message
  USE mo_jsb_io_netcdf,      ONLY: t_input_file

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_config, t_jsb_config_p, t_jsb_model_config, new_model_config

  ! Type for model configuration
  TYPE t_jsb_model_config
    CHARACTER(len=filename_max)    :: grid_filename
    CHARACTER(len=filename_max)    :: fract_filename
    CHARACTER(len=30)              :: grid_name
    CHARACTER(len=40)              :: usecase
    LOGICAL                        :: use_tmx
    LOGICAL                        :: use_lakes
    LOGICAL                        :: use_glacier
    LOGICAL                        :: use_quincy
    CHARACTER(len=SHORT_NAME_LEN)  :: tpe_scheme       ! For terraplanet setup: open/closed
    CHARACTER(len=SHORT_NAME_LEN)  :: hsm_mode
    LOGICAL                        :: l_compat401      ! Use configuration compatible with jsbach 4.01
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: output_tiles(:)  ! List of tile names for which output should be generated
    LOGICAL                        :: relative_fractions_in_file  ! true: tile fractions in fract_filename are relative to parent tile
    LOGICAL                        :: init_from_ifs    ! Initialize from IFS analysis
    CHARACTER(len=:), ALLOCATABLE  :: ifs_filename
    TYPE(t_input_file)             :: ifs_input_file

  END TYPE t_jsb_model_config

  ! Abstract type for process configurations
  TYPE, ABSTRACT :: t_jsb_config
    LOGICAL :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    CHARACTER(LEN=filename_max) :: namelist_filename
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: output_tiles(:)  ! List of tile names for which output should be generated
  CONTAINS
    PROCEDURE (Init_config), DEFERRED, PASS(config) :: Init
  END TYPE t_jsb_config

  TYPE t_jsb_config_p
    CLASS(t_jsb_config), POINTER :: p
  END TYPE t_jsb_config_p

  ABSTRACT INTERFACE
    SUBROUTINE Init_config(config)
      IMPORT :: t_jsb_config
      CLASS(t_jsb_config), INTENT(inout) :: config
    END SUBROUTINE Init_config
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_config_class'

CONTAINS

  FUNCTION new_model_config(namelist_filename) RESULT(model_config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CHARACTER(len=*),         INTENT(in) :: namelist_filename
    TYPE(t_jsb_model_config), POINTER    :: model_config

    CHARACTER(len=filename_max) :: grid_filename
    CHARACTER(len=30)           :: grid_name

    NAMELIST /jsb_grid_nml/   &
      grid_filename,          &
      grid_name

    LOGICAL                       :: use_tmx
    LOGICAL                       :: use_lakes
    LOGICAL                       :: use_glacier
    LOGICAL                       :: use_quincy
    LOGICAL                       :: l_compat401
    CHARACTER(len=filename_max)   :: fract_filename
    CHARACTER(len=40)             :: usecase
    CHARACTER(len=SHORT_NAME_LEN) :: tpe_scheme
    CHARACTER(len=SHORT_NAME_LEN) :: hsm_mode
    CHARACTER(len=SHORT_NAME_LEN) :: output_tiles(99)
    LOGICAL                       :: relative_fractions_in_file
    LOGICAL                       :: init_from_ifs
    CHARACTER(len=filename_max)   :: ifs_filename

    NAMELIST /jsb_model_nml/  &
      use_tmx, &
      use_lakes, &
      use_glacier, &
      use_quincy, &
      l_compat401, &
      usecase, &
      tpe_scheme, &
      hsm_mode, &
      fract_filename, &
      output_tiles, &
      relative_fractions_in_file, &
      init_from_ifs, &
      ifs_filename

    INTEGER :: nml_handler, nml_unit, istat, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_model_config'

    CALL message(TRIM(routine), 'Starting model configuration from '//TRIM(namelist_filename))

    ALLOCATE(model_config)

    ! Set defaults
    use_tmx                    = .FALSE.
    use_lakes                  = .TRUE.
    use_glacier                = .TRUE.
    use_quincy                 = .FALSE.
    l_compat401                = .FALSE.
    grid_filename              = ''
    grid_name                  = ''
    usecase                    = ''
    tpe_scheme                 = ''
    hsm_mode                   = 'simple'
    fract_filename             = 'bc_land_frac.nc'
    output_tiles(:)            = ''
    relative_fractions_in_file = .TRUE.
    init_from_ifs              = .FALSE.
    ifs_filename               = 'ifs2icon.nc'

    ! Read namelist
    nml_handler = open_nml(TRIM(namelist_filename))

    nml_unit = position_nml('jsb_grid_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_grid_nml)

    ! Write resulting values in config memory
    model_config%grid_filename = grid_filename
    model_config%grid_name     = grid_name

    nml_unit =  position_nml('jsb_model_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_model_nml)
    CALL close_nml(nml_handler)

    model_config%use_tmx                    = use_tmx
    model_config%use_lakes                  = use_lakes
    model_config%use_glacier                = use_glacier
    model_config%use_quincy                 = use_quincy
    model_config%l_compat401                = l_compat401
    model_config%usecase                    = usecase
    model_config%tpe_scheme                 = tpe_scheme
    model_config%hsm_mode                   = hsm_mode
    model_config%fract_filename             = fract_filename
    model_config%relative_fractions_in_file = relative_fractions_in_file
    model_config%init_from_ifs              = init_from_ifs
    model_config%ifs_filename               = TRIM(ifs_filename)

    DO i=1,99
      IF (output_tiles(i) == '') EXIT
    END DO
    i = i - 1  ! Number of tile names in namelist variable output_tiles
    IF (i == 0) THEN             ! Default is 'box' if no tiles specified in namelist
      i = 1
      output_tiles(1) = 'box'
    END IF
    ALLOCATE(model_config%output_tiles(i))
    model_config%output_tiles(1:i) = output_tiles(1:i)

  END FUNCTION new_model_config

#endif
END MODULE mo_jsb_config_class
