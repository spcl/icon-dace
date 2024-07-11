!> Contains structures and methods for soil and snow energy config
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
MODULE mo_sse_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: real2string, int2string
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sse_config

  TYPE, EXTENDS(t_jsb_config) :: t_sse_config
    INTEGER  :: nsnow           !< Maximum number of snow layers
    LOGICAL  :: l_snow          !< T: Use multi-layer snow scheme;
                                !  Note: This affects thermal snow and soil properties:
                                !  l_dynsnow, l_heat_cap_dyn and l_heat_cond_dyn are only effective if l_snow=.true.
    LOGICAL  :: l_dynsnow       !< T: Calculate snow heat cond. and heat cap. dynamically depending on snow density
    LOGICAL  :: l_heat_cap_dyn  !< T: Dynamic calculation of soil heat capacity;
                                !  F: Static map or FAO data, depending on l_heat_cap_map
    LOGICAL  :: l_heat_cond_dyn !< T: Dynamic calculation of soil heat conductivity;
                                !  F: Static map or FAO data, depending on l_heat_cond_map
    !
    LOGICAL  :: l_heat_cap_map  !< T: Use soil heat capacity from input map, F: derive from FAO
    LOGICAL  :: l_heat_cond_map !< T: Use soil heat conductivity from input map, F: derive from heat cap. and FAO thermal diff.
    LOGICAL  :: l_soil_texture  !< T: Deduce mineral soil thermal parameters from soil texture
    LOGICAL  :: l_freeze        !< T: Consider freezing and thawing of soil water
    LOGICAL  :: l_supercool     !< T: Allow for supercooled soil water
    REAL(wp) :: w_soil_critical !< Critical water/ice content in upper soil layer for correction of new surface
                                !! temperature for freezing/melting [m water equivalent]
    !! quincy
    ! INTEGER       :: nsoil_energy         !< number of soil layers for soil energy calculations   ! depreciated feature, but needed for now, see quincy-trac ticket #521
    ! INTEGER       :: nsoil_water          !< number of soil layers for soil moisture calculations ! depreciated feature, but needed for now, see quincy-trac ticket #521
    REAL(wp)      :: soil_depth           !< actual soil and rooting depth
    REAL(wp)      :: wsr_capacity         !< Water holding capacity of ground skin reservoir [m water equivalent]
    REAL(wp)      :: wsn_capacity         !< Water holding capacity of ground snow reservoir (max snow depth) [m water equivalent]
    REAL(wp)      :: soil_awc_prescribe, &
                     soil_theta_prescribe
    REAL(wp)      :: soil_sand            !< soil sand proportion - site specific from forcing info file
    REAL(wp)      :: soil_silt            !< soil silt proportion - site specific from forcing info file
    REAL(wp)      :: soil_clay            !< soil clay proportion - recalculated from: clay = 1.0 -sand -silt
    REAL(wp)      :: bulk_density

  CONTAINS
    PROCEDURE :: Init => Init_sse_config
  END TYPE t_sse_config

  INTEGER, PARAMETER :: max_snow_layers = 20
  INTEGER, PARAMETER :: max_soil_layers = 20

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_config_class'

CONTAINS

  SUBROUTINE Init_sse_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND, ZAXIS_GENERIC
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_sse_constants,      ONLY: snow_depth_min
    USE mo_jsb_math_constants, ONLY: eps8

    CLASS(t_sse_config), INTENT(inout) :: config

    LOGICAL  :: active, l_snow, l_dynsnow, l_heat_cap_dyn, l_heat_cond_dyn, l_heat_cap_map, l_heat_cond_map
    LOGICAL  :: l_freeze, l_supercool, l_soil_texture
    REAL(wp) :: w_soil_critical
    INTEGER  :: nsnow
    REAL(wp) :: dz_snow(max_snow_layers)
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_sse_nml/                      &
      & active,                                 &
      & nsnow, dz_snow, l_snow, l_dynsnow,      &
      & l_heat_cap_dyn, l_heat_cond_dyn,        &
      & l_heat_cap_map, l_heat_cond_map,        &
      & l_soil_texture,                         &
      & l_freeze, l_supercool, w_soil_critical, &
      & ic_filename, bc_filename

    INTEGER :: nml_handler, nml_unit, istat, i
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)
    REAL(wp), ALLOCATABLE :: depths(:), mids(:)
    INTEGER :: nsoil

    TYPE(t_jsb_vgrid), POINTER :: soil_vgrid, snow_vgrid

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_sse_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting soil and snow energy configuration')

    ! Set defaults
    active          = .TRUE.
    nsnow            = 5
    dz_snow(:)       = 0._wp
    dz_snow(1:nsnow) = 0.05_wp
    l_snow           = .TRUE.
    l_dynsnow        = .TRUE.
    l_heat_cap_dyn   = .TRUE.
    l_heat_cond_dyn  = .TRUE.
    l_heat_cap_map   = .FALSE.
    l_heat_cond_map  = .FALSE.
    l_soil_texture   = .FALSE.
    l_freeze         = .TRUE.
    l_supercool      = .TRUE.
    w_soil_critical  = 5.85036E-3_wp
    ic_filename      = 'ic_land_soil.nc'
    bc_filename      = 'bc_land_soil.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_sse_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_sse_nml)

    CALL close_nml(nml_handler)

    config%active        = active
    config%ic_filename   = ic_filename
    config%bc_filename   = bc_filename

    config%nsnow         = nsnow
    IF (l_snow .AND. nsnow > max_snow_layers) THEN
      CALL finish(TRIM(routine), 'Too many snow layers, maximum is '//TRIM(int2string(max_snow_layers)))
    END IF
    IF (l_snow .AND. nsnow < 3) THEN
      CALL finish(TRIM(routine), 'Number of snow layers must be larger than 2.')
    END IF
    config%l_snow        = l_snow
    IF (l_snow) THEN
      CALL message(TRIM(routine), 'Using '//TRIM(int2string(nsnow))//' snow layers.')
    ELSE
      CALL message(TRIM(routine), 'Not using multi-layer snow model')
    END IF

    IF (.NOT. l_snow .AND. l_dynsnow) THEN
      config%l_dynsnow = .FALSE.
      CALL message(TRIM(routine), 'l_dynsnow not in effect for l_snow=.false.')
    ELSE
      config%l_dynsnow       = l_dynsnow
    END IF
    IF (.NOT. l_snow .AND. l_heat_cap_dyn) THEN
      config%l_heat_cap_dyn = .FALSE.
      CALL message(TRIM(routine), 'l_heat_cap_dyn not in effect for l_snow=.false.')
    ELSE
      config%l_heat_cap_dyn  = l_heat_cap_dyn
    END IF
    IF (.NOT. l_snow .AND. l_heat_cond_dyn) THEN
      config%l_heat_cond_dyn = .FALSE.
      CALL message(TRIM(routine), 'l_heat_cond_dyn not in effect for l_snow=.false.')
    ELSE
      config%l_heat_cond_dyn = l_heat_cond_dyn
    END IF
    config%l_heat_cap_map  = l_heat_cap_map
    config%l_heat_cond_map = l_heat_cond_map
    config%l_soil_texture  = l_soil_texture
    IF (l_soil_texture) THEN
      CALL message(TRIM(routine), 'Using soil texture to derive heat capacity and heat conductivity')
    END IF
    config%l_freeze        = l_freeze
    config%l_supercool     = l_supercool

    config%w_soil_critical     = w_soil_critical
    CALL message(TRIM(routine), 'Critical water/ice content in upper soil layer for correction of '// &
      &                         'surface temperature for freezing/melting: '//TRIM(real2string(w_soil_critical)))

    IF (.NOT. active) RETURN

    input_file = jsb_netcdf_open_input(ic_filename)

    ! @todo: At the moment, the soil layers for the energy calculations are the same as for the hydrology
    ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
    nsoil = SIZE(ptr_1D)
    ALLOCATE(depths(nsoil+1))
    ALLOCATE(mids(nsoil))
    depths(1) = 0._wp
    depths(2:nsoil+1) = ptr_1D(1:nsoil)
    mids(1:nsoil) = (depths(1:nsoil) + depths(2:nsoil+1)) / 2._wp
    DEALLOCATE(ptr_1D)

    CALL input_file%Close()

    soil_vgrid  => new_vgrid('soil_depth_energy', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
      & levels  = mids                 (1:nsoil  ),                              &
      & lbounds = depths               (1:nsoil  ),                              &
      & ubounds = depths               (2:nsoil+1),                              &
      & units='m')
    CALL register_vgrid(soil_vgrid)
    WRITE(message_text, *) 'Soil levels in soil energy (upper) [m]: ', soil_vgrid%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in soil energy (mid)   [m]: ', soil_vgrid%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in soil energy (lower) [m]: ', soil_vgrid%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil level depths in soil energy   [m]: ', soil_vgrid%dz
    CALL message(TRIM(routine), message_text)
    DEALLOCATE(depths, mids)

    IF (l_snow) THEN
      IF (dz_snow(1) < snow_depth_min)                                                                 &
        & CALL finish(TRIM(routine), 'Depth of first snow layer should be larger than snow_depth_min=' &
        &                            //real2string(snow_depth_min))
      ALLOCATE(depths(nsnow+1))
      ALLOCATE(mids(nsnow))
      depths(1) = 0._wp
      DO i=1,nsnow
        depths(i+1) = depths(i) + dz_snow(i)
      END DO
      mids(1:nsnow) = (depths(1:nsnow) + depths(2:nsnow+1)) / 2._wp

      snow_vgrid  => new_vgrid('snow_depth_energy', ZAXIS_GENERIC, nsnow,     &
        & levels  = mids                 (1:nsnow  ),                         &
        & lbounds = depths               (1:nsnow  ),                         &
        & ubounds = depths               (2:nsnow+1),                         &
        & units='m')
      CALL register_vgrid(snow_vgrid)
      WRITE(message_text, *) 'Snow levels in soil energy (upper) [m]: ', snow_vgrid%lbounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow levels in soil energy (mid)   [m]: ', snow_vgrid%levels
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow levels in soil energy (lower) [m]: ', snow_vgrid%ubounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow level depths in soil energy   [m]: ', snow_vgrid%dz
      CALL message(TRIM(routine), message_text)
      DEALLOCATE(depths, mids)
    END IF


  END SUBROUTINE Init_sse_config

#endif
END MODULE mo_sse_config_class
