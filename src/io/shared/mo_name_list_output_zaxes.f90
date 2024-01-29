! Module handling the specification of vertical axes for the output
! module.
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
!
!
! --------------------------------
! Details of the implementation
! --------------------------------
!
! Axis creation:
!
! The CDI routines for the creation of the vertical axes are called
! for each output file separately in "mo_name_list_output_init",
! subroutine "setup_output_vlist".
!
! For model level output, all axes (atmospheric, surface, soil etc.)
! are created.  For pressure, height, and isentropic level output
! only the respective axis is created together with a surface axis
! for the grid information.
!
! Note: The vertical axes for the output of the ocean do not support
!       level selection yet.
!
! Derived data type "t_level_selection":
!
! Vertical levels are selected via objects of the derived data type
! "t_level_selection" (module "mo_level_selection").  If such
! a data object has been initialized for the output file, then the
! vertical axis definitions below create CDI axis objects for the
! selected levels only.
! The write routines in the module "mo_name_list_output" skip levels
! which are not part of the selection.
!
! Furthermore, if a "t_level_selection" object is present, then the
! memory windows for the asynchronous one-sided MPI communication are
! adjusted to the reduced level size
! "output_file%level_selection%n_selected".

MODULE mo_name_list_output_zaxes

  USE ISO_C_BINDING,                        ONLY: C_SIGNED_CHAR
  USE mo_kind,                              ONLY: wp, dp
  USE mo_impl_constants,                    ONLY: SUCCESS, inwp
  USE mo_exception,                         ONLY: finish
  USE mo_zaxis_type,                        ONLY: zaxisTypeList,                                                 &
    &                                             ZA_depth_below_sea, ZA_depth_below_sea_half, ZA_GENERIC_ICE,   &
    &                                             ZA_surface, ZA_isentropic, ZA_altitude, ZA_pressure,           &
    &                                             ZA_cloud_base, ZA_cloud_top, ZA_depth_below_land,              &
    &                                             ZA_depth_below_land_p1, ZA_depth_runoff_g, ZA_depth_runoff_s,  &
    &                                             ZA_height_10m, ZA_height_2m, ZA_isotherm_zero, ZA_lake_bottom, &
    &                                             ZA_lake_bottom_half, ZA_meansea, ZA_mix_layer, ZA_pressure_0,  &
    &                                             ZA_pressure_400, ZA_pressure_800, ZA_ATMOSPHERE,               &
    &                                             ZA_PRES_FL_BOT_TOP,                                            &
    &                                             ZA_reference, ZA_reference_half, ZA_reference_half_hhl,        &
    &                                             ZA_sediment_bottom_tw_half, ZA_snow, ZA_snow_half, ZA_toa,     &
    &                                             ZA_OCEAN_SEDIMENT, ZA_height_2m_layer, ZA_ECHOTOP,             &
    &                                             ZA_TROPOPAUSE, ZA_WSHEAR, ZA_PRESSURE_LAPSERATE, ZA_SRH
  USE mo_level_selection_types,             ONLY: t_level_selection
  USE mo_util_vgrid_types,                  ONLY: vgrid_buffer
  USE mo_math_utilities,                    ONLY: set_zlev, t_value_set
  USE mo_run_config,                        ONLY: num_lev, iforcing, nlev
  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxis, t_verticalAxisList
#ifndef __NO_ICON_ATMO__
  USE mo_io_config,                         ONLY: echotop_meta, wshear_uv_heights, srh_heights
  USE mo_nonhydrostatic_config,             ONLY: ivctype
  USE mo_lnd_nwp_config,                    ONLY: nlev_snow, zml_soil
#endif
#ifndef __NO_ICON_OCEAN__
  USE mo_ocean_nml,                         ONLY: dzlev_m,lhamocc, n_zlev
  USE mo_hamocc_nml,                        ONLY: ks, ksp, dzsed
#endif

  IMPLICIT NONE

  PRIVATE

  ! subroutines
  PUBLIC :: setup_ml_axes_atmo
  PUBLIC :: setup_pl_axis_atmo
  PUBLIC :: setup_hl_axis_atmo
  PUBLIC :: setup_il_axis_atmo
  PUBLIC :: setup_zaxes_oce
  PUBLIC :: vertical_axis

  PUBLIC :: single_level_axis

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_zaxes'


CONTAINS

  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axes for output module: Atmosphere component, model levels
  !
  SUBROUTINE setup_ml_axes_atmo(verticalAxisList, level_selection, log_patch_id)
    TYPE(t_verticalAxisList), INTENT(INOUT),TARGET   :: verticalAxisList
    TYPE(t_level_selection),  INTENT(IN), POINTER    :: level_selection
    INTEGER,                  INTENT(IN)             :: log_patch_id
    ! local variables
    REAL(dp), ALLOCATABLE             :: levels(:), lbounds(:), ubounds(:)
    INTEGER                           :: k, nlev, nlevp1, znlev_soil
    INTEGER                           :: n_echotop, n_wshear, n_srh
    TYPE(t_verticalAxisList), POINTER :: it
#ifndef __NO_ICON_ATMO__

    nlev            = num_lev(log_patch_id)
    nlevp1          = num_lev(log_patch_id) + 1

    ! introduce temporary variable znlev_soil, since global variable
    ! nlev_soil is unknown to the I/O-Processor.
    IF ( iforcing == inwp ) THEN
      znlev_soil = SIZE(zml_soil)
    ELSE
      znlev_soil = 0
    ENDIF

    ! --------------------------------------------------------------------------------------
    ! Definitions for single levels --------------------------------------------------------
    ! --------------------------------------------------------------------------------------

    ! surface level
    CALL verticalAxisList%append(single_level_axis(ZA_surface))

    ! CLOUD BASE LEVEL
    CALL verticalAxisList%append(single_level_axis(ZA_cloud_base))

    ! CLOUD TOP LEVEL
    CALL verticalAxisList%append(single_level_axis(ZA_cloud_top))

    ! LEVEL of 0\deg C isotherm
    CALL verticalAxisList%append(single_level_axis(ZA_isotherm_zero))

    ! Define axis for output on mean sea level
    CALL verticalAxisList%append(single_level_axis(ZA_meansea))

    ! Specific soil axis for Runoff_s
    CALL verticalAxisList%append(single_layer_axis(ZA_depth_runoff_s, 0._dp, 0.1_dp, unit="m"))

    ! Specific soil axis for Runoff_g
    CALL verticalAxisList%append(single_layer_axis(ZA_depth_runoff_g, 0.1_dp, 1._dp, unit="m"))

    ! Specified height level above ground: 2m
    CALL verticalAxisList%append(single_level_axis(ZA_height_2m, opt_level_value=2._dp))

    ! Specified height level above ground: 10m
    CALL verticalAxisList%append(single_level_axis(ZA_height_10m, opt_level_value=10._dp))

    ! Top of atmosphere
    CALL verticalAxisList%append(single_level_axis(ZA_toa))

    ! Bottom of sediment layer penetrated by thermal wave (interface, i.e. only typeOfFirstFixedSurface)
    CALL verticalAxisList%append(single_level_axis(ZA_sediment_bottom_tw_half, opt_unit="m"))

    ! Lake bottom half (interface, i.e. only typeOfFirstFixedSurface)
    CALL verticalAxisList%append(single_level_axis(ZA_lake_bottom_half, opt_unit="m"))

    ! for having ice variable in the atmosphere (like AMIP)
    CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_GENERIC_ICE), 1))

    ! for having variable on the tropopause niveau
    CALL verticalAxisList%append(single_level_axis(ZA_TROPOPAUSE, opt_grib2_level_type=7))

    ! --------------------------------------------------------------------------------------
    ! Definitions for single layers --------------------------------------------------------
    ! --------------------------------------------------------------------------------------

    ! Specified height level above ground: 2m
    ! Layered-Version, which allows to re-set typeOfFirst/SeconFixedSurface
    ! without loosing the information stored in scaledValueOfFirst/SecondFixedSurface
    CALL verticalAxisList%append(single_layer_axis(ZA_height_2m_layer, 2._dp,   0._dp, "m"))


    !DR Note that for this particular axis in GME and COSMO,
    !DR   scaleFactorOfSecondFixedSurface = scaledValueOfSecondFixedSurface = missing
    !DR whereas in ICON
    !DR   scaleFactorOfSecondFixedSurface = scaledValueOfSecondFixedSurface = 0
    !DR This is because, CDI does not allow layer-type vertical axis without setting
    !DR scaleFactorOfSecondFixedSurface and scaledValueOfSecondFixedSurface.
    !
    ! Isobaric surface 800 hPa (layer)
    CALL verticalAxisList%append(single_layer_axis(ZA_pressure_800, 800._dp,   0._dp, "hPa"))
    ! Isobaric surface 400 hPa (layer)
    CALL verticalAxisList%append(single_layer_axis(ZA_pressure_400, 400._dp, 800._dp, "hPa"))
    ! Isobaric surface 0 hPa (layer)
    CALL verticalAxisList%append(single_layer_axis(ZA_pressure_0,   0._dp, 400._dp, "hPa"))

    ! Specific vertical axis for Lake-model ------------------------------------------------
    ! Lake bottom (we define it as a layer in order to be able to re-set
    ! either the first- or secondFixedSurfaces if necessary)
    CALL verticalAxisList%append(single_layer_axis(ZA_lake_bottom, 0._dp, 0._dp, "m"))

    ! Mixing layer (we define it as a layer in order to be able to re-set
    ! either the first- or secondFixedSurfaces if necessary)
    CALL verticalAxisList%append(single_layer_axis(ZA_mix_layer, 1._dp, 0._dp, "m"))

    ! ICON-ART products - Maximum total mass or air concentration in flight level range
    !                     defined by pressure layers
    CALL verticalAxisList%append(single_layer_axis(ZA_PRES_FL_BOT_TOP,  50.00_dp, 1013.25_dp, "hPa"))
    ! ICON-ART products - Colummn integrated total mass concentration (entire atmosphere)
    CALL verticalAxisList%append(single_level_axis(ZA_ATMOSPHERE))

    ! --------------------------------------------------------------------------------------
    ! Definitions for echotops of simulated radar reflectivity
    ! --------------------------------------------------------------------------------------
    n_echotop = echotop_meta(log_patch_id)%nechotop
    IF ( n_echotop > 0) THEN
      CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_ECHOTOP), &
                                   n_echotop, &
                                   zaxisLevels=echotop_meta(log_patch_id)%dbzthresh(1:n_echotop), &
                                   zaxisUnits="dBZ") )
    END IF

    ! --------------------------------------------------------------------------------------
    ! Definitions for vertical windshear components
    !   (u/v differences between a height AGL and lowest model level AGL)
    ! --------------------------------------------------------------------------------------
    n_wshear = COUNT(wshear_uv_heights > 0.0_wp)
    IF ( n_wshear > 0) THEN
      CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_WSHEAR), &
           &                                         n_wshear, &
           &                                         zaxisLevels  = wshear_uv_heights(1:n_wshear), &
           &                                         zaxisLbounds = wshear_uv_heights(1:n_wshear), &
           &                                         zaxisUbounds = (/ (0.0_wp, k=1,n_wshear) /),  &
           &                                         zaxisUnits="m"))
    END IF
    
    ! --------------------------------------------------------------------------------------
    ! Definitions for storm relative helicity
    !   (integral between a height AGL and GND)
    ! --------------------------------------------------------------------------------------
    n_srh = COUNT(srh_heights > 0.0_wp)
    IF ( n_srh > 0) THEN
      CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_SRH), &
           &                                         n_srh, &
           &                                         zaxisLevels  = srh_heights(1:n_srh), &
           &                                         zaxisLbounds = srh_heights(1:n_srh), &
           &                                         zaxisUbounds = (/ (0.0_wp, k=1,n_srh) /),&
           &                                         zaxisUnits="m"))
    END IF
    
    ! Isobaric pressure layer 500 hPa - 850 hPa for the temperature lapse_rate:
    CALL verticalAxisList%append(single_layer_axis(ZA_pressure_lapserate, 500._dp, 850._dp, "hPa"))

    ! --------------------------------------------------------------------------------------
    ! Definitions for reference grids (ZAXIS_REFERENCE) ------------------------------------
    ! --------------------------------------------------------------------------------------

    ! REFERENCE
    CALL verticalAxisList%append(vertical_axis(ZA_reference, nlev,                         &
      &                           levels           = (/ ( REAL(k,dp),   k=1,nlevp1 ) /),   &
      &                           level_selection  = level_selection,                      &
      &                           opt_set_bounds   = .TRUE.,                               &
      &                           opt_number       = get_numberOfVgridUsed(ivctype),       &
      &                           opt_uuid         = vgrid_buffer(log_patch_id)%uuid%DATA, &
      &                           opt_nlevref      = nlevp1))

    ! REFERENCE HALF
    CALL verticalAxisList%append(vertical_axis(ZA_reference_half, nlevp1,                  &
      &                           levels          = (/ ( REAL(k,dp),   k=1,nlevp1 ) /),    &
      &                           level_selection = level_selection,                       &
      &                           opt_number      = get_numberOfVgridUsed(ivctype),        &
      &                           opt_uuid        = vgrid_buffer(log_patch_id)%uuid%DATA,  &
      &                           opt_nlevref     = nlevp1 ))

    ! REFERENCE (special version for HHL)
    CALL verticalAxisList%append(vertical_axis(ZA_reference_half_hhl, nlevp1,                   &
      &                           levels                = (/ ( REAL(k,dp),   k=1,nlevp1 ) /),   &
      &                           level_selection       = level_selection,                      &
      &                           opt_set_bounds        = .TRUE.,                               &
      &                           opt_set_ubounds_value = 0._dp,                                &
      &                           opt_number            = get_numberOfVgridUsed(ivctype),       &
      &                           opt_uuid              = vgrid_buffer(log_patch_id)%uuid%DATA, &
      &                           opt_nlevref           = nlevp1 ))


    ! --------------------------------------------------------------------------------------
    ! Axes for soil model (ZAXIS_DEPTH_BELOW_LAND) -----------------------------------------
    ! --------------------------------------------------------------------------------------

    IF ( iforcing == inwp ) THEN
      ALLOCATE(levels(znlev_soil+1))
      levels(1) = 0._dp
      DO k = 1, znlev_soil
        levels(k+1) = REAL(zml_soil(k)*1000._wp,dp)  ! in mm
      END DO
      CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_depth_below_land_p1), &
        &                                         znlev_soil+1, zaxisLevels=levels, zaxisUnits="mm"))
      DEALLOCATE(levels)

      ALLOCATE(lbounds(znlev_soil), ubounds(znlev_soil), levels(znlev_soil))
      lbounds(1) = 0._dp   ! surface
      DO k = 2, znlev_soil
        lbounds(k)   = REAL((zml_soil(k-1) + (zml_soil(k-1) - lbounds(k-1))),dp)
      ENDDO
      DO k = 1, znlev_soil
        ubounds(k) = REAL((zml_soil(k) + (zml_soil(k) - lbounds(k))),dp)
        levels(k)  = REAL(zml_soil(k)*1000._wp,dp)
      ENDDO
      ubounds(:) = ubounds(:) * 1000._dp        ! in mm
      lbounds(:) = lbounds(:) * 1000._dp        ! in mm
      CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_depth_below_land), &
        &                                         znlev_soil, zaxisLevels=levels,              &
        &                                         zaxisLbounds=lbounds, zaxisUbounds=ubounds,  &
        &                                         zaxisUnits="mm"))
      DEALLOCATE(lbounds, ubounds, levels)
    ENDIF

    ! --------------------------------------------------------------------------------------
    ! Axes for multi-layer snow model (ZAXIS_SNOW) -----------------------------------------
    ! --------------------------------------------------------------------------------------

    IF (nlev_snow > 0) THEN
      ! SNOW-layer axis (for multi-layer snow model)
      ALLOCATE(levels(nlev_snow), lbounds(nlev_snow), ubounds(nlev_snow))
      DO k = 1, nlev_snow
        lbounds(k) = REAL(k,dp)
        levels(k)  = REAL(k,dp)
      ENDDO
      DO k = 1, nlev_snow
        ubounds(k) = REAL(k+1,dp)
      ENDDO
      CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_SNOW), nlev_snow, &
        &                                         zaxisLevels=levels, zaxisLbounds=lbounds,     &
        &                                         zaxisUbounds=ubounds))
      DEALLOCATE(levels, lbounds, ubounds)

      ALLOCATE(levels(nlev_snow+1))
      DO k = 1, nlev_snow+1
        levels(k) = REAL(k,dp)
      END DO
      CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_SNOW_HALF),  &
        &                                         nlev_snow+1, zaxisLevels=levels))
      DEALLOCATE(levels)
    END IF

#endif
    ! #ifndef __NO_ICON_ATMO__

    it => verticalAxisList
    DO
      IF (.NOT. ASSOCIATED(it%next))  EXIT
      it => it%next
    END DO
  END SUBROUTINE setup_ml_axes_atmo


  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axis for output module: Pressure levels
  !
  SUBROUTINE setup_pl_axis_atmo(verticalAxisList, levels, level_selection)
    TYPE(t_verticalAxisList), INTENT(INOUT)       :: verticalAxisList
    TYPE (t_value_set),       INTENT(IN)          :: levels
    TYPE(t_level_selection),  INTENT(IN), POINTER :: level_selection

#ifndef __NO_ICON_ATMO__
    ! surface level (required, e.g., for RLON, RLAT)
    CALL verticalAxisList%append(single_level_axis(ZA_surface))

    ! p-axis
    CALL verticalAxisList%append(vertical_axis(ZA_pressure,               &
      &                       SIZE(levels%values),                        &
      &                       levels                = levels%values,      &
      &                       level_selection       = level_selection,    &
      &                       opt_set_vct_as_levels = .TRUE.))
#endif

  END SUBROUTINE setup_pl_axis_atmo


  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axis for output module: Height levels
  !
  SUBROUTINE setup_hl_axis_atmo(verticalAxisList, levels, level_selection)
    TYPE(t_verticalAxisList), INTENT(INOUT)       :: verticalAxisList
    TYPE (t_value_set),       INTENT(IN)          :: levels
    TYPE(t_level_selection),  INTENT(IN), POINTER :: level_selection

#ifndef __NO_ICON_ATMO__
    ! surface level (required, e.g., for RLON, RLAT)
    CALL verticalAxisList%append(single_level_axis(ZA_surface))
    ! Altitude above mean sea level
    CALL verticalAxisList%append(vertical_axis(ZA_altitude,                  &
      &                          SIZE(levels%values),                        &
      &                          levels                = levels%values,      &
      &                          level_selection       = level_selection,    &
      &                          opt_set_vct_as_levels = .TRUE.))
#endif

  END SUBROUTINE setup_hl_axis_atmo


  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axis for output module: Isentropes
  !
  SUBROUTINE setup_il_axis_atmo(verticalAxisList, levels, level_selection)
    TYPE(t_verticalAxisList), INTENT(INOUT)       :: verticalAxisList
    TYPE (t_value_set),       INTENT(IN)          :: levels
    TYPE(t_level_selection),  INTENT(IN), POINTER :: level_selection

#ifndef __NO_ICON_ATMO__
    ! surface level (required, e.g., for RLON, RLAT)
    CALL verticalAxisList%append(single_level_axis(ZA_surface))
    ! i-axis (isentropes)
    CALL verticalAxisList%append(vertical_axis(ZA_isentropic,                &
      &                          SIZE(levels%values),                        &
      &                          levels                = levels%values,      &
      &                          level_selection       = level_selection,    &
      &                          opt_set_vct_as_levels = .TRUE.))

#endif

  END SUBROUTINE setup_il_axis_atmo


  ! --------------------------------------------------------------------------------------
  !> Setup of vertical axes for output module: Ocean component
  !
  SUBROUTINE setup_zaxes_oce(verticalAxisList, level_selection)
  
    TYPE(t_verticalAxisList), INTENT(INOUT) :: verticalAxisList
    TYPE(t_level_selection),  INTENT(IN), POINTER :: level_selection
    ! local variables
    REAL(wp), ALLOCATABLE             :: levels_i(:), levels_m(:)
    REAL(wp), ALLOCATABLE             :: levels_s(:), levels_sp(:)
    
#ifndef __NO_ICON_OCEAN__
    ALLOCATE(levels_i(n_zlev+1), levels_m(n_zlev))
    CALL set_zlev(levels_i, levels_m, n_zlev, dzlev_m)

    CALL verticalAxisList%append(single_level_axis(ZA_surface))
    CALL verticalAxisList%append(vertical_axis(ZA_depth_below_sea, &
      &                          n_zlev,                           &
      &                          levels = REAL(levels_m,wp), &
      &                          level_selection=level_selection))
    CALL verticalAxisList%append(vertical_axis(ZA_depth_below_sea_half, &
      &                          n_zlev+1,                              &
      &                          levels = REAL(levels_i,wp),            &
      &                          level_selection=level_selection))
    CALL verticalAxisList%append(single_level_axis(ZA_GENERIC_ICE))

    DEALLOCATE(levels_i, levels_m)

    CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_GENERIC_ICE), 1))
    if(lhamocc)then
    ! ocean sediment
    ALLOCATE(levels_s(ks), levels_sp(ksp))

    CALL set_zlev(levels_sp, levels_s, ks, dzsed)
    CALL verticalAxisList%append(t_verticalAxis(zaxisTypeList%getEntry(ZA_OCEAN_SEDIMENT), ks, &
      &                                         zaxisLevels=REAL(levels_s,dp)))
    DEALLOCATE(levels_s, levels_sp)
    endif
#endif

  END SUBROUTINE setup_zaxes_oce


  ! --------------------------------------------------------------------------------------
  !> Utility function: defines z-axis with a single level
  !
  FUNCTION single_level_axis(za_type, opt_level_value, opt_unit, opt_grib2_level_type)
    TYPE(t_verticalAxis) :: single_level_axis
    INTEGER,          INTENT(IN)           :: za_type               !< ICON-internal axis ID (see mo_zaxis_type)
    REAL(dp),         INTENT(IN), OPTIONAL :: opt_level_value       !< level value
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_unit              !< axis unit
    INTEGER,          INTENT(in), OPTIONAL :: opt_grib2_level_type  !< level type as defined in GRIB2
    ! local variables
    REAL(dp) :: levels(1)

    levels(1) = 0.0_dp
    IF (PRESENT(opt_level_value)) levels(1) = opt_level_value

    single_level_axis = t_verticalAxis(zaxisTypeList%getEntry(za_type), 1,  &
         &                             zaxisLevels=levels,                  &
         &                             zaxisDefLtype=opt_grib2_level_type)

    IF (PRESENT(opt_unit))  CALL single_level_axis%set(zaxisUnits=TRIM(opt_unit))
  END FUNCTION single_level_axis


  ! --------------------------------------------------------------------------------------
  !> Utility function: defines z-axis with a single *layer*
  !
  FUNCTION single_layer_axis(za_type, lbound, ubound, unit)
    TYPE(t_verticalAxis) :: single_layer_axis
    INTEGER,             INTENT(IN)    :: za_type        !< ICON-internal axis ID (see mo_zaxis_type)
    REAL(dp),            INTENT(IN)    :: lbound, ubound !< lower and upper bound
    CHARACTER(LEN=*),    INTENT(IN)    :: unit           !< axis unit
    ! local variables
    REAL(dp) :: levels(1), lbounds(1), ubounds(1)

    lbounds(1) = lbound
    ubounds(1) = ubound
    levels(1)  = lbound

    single_layer_axis = t_verticalAxis(zaxisTypeList%getEntry(za_type), 1,           &
      &                                zaxisLevels=levels, zaxisLbounds=lbounds,     &
      &                                zaxisUbounds=ubounds, zaxisUnits=TRIM(unit))
  END FUNCTION single_layer_axis


  ! --------------------------------------------------------------------------------------
  !> Utility function: defines z-axis with given levels, lower and upper bounds.
  !
  FUNCTION vertical_axis(za_type, in_nlevs, levels, level_selection,  &
    &                    opt_set_bounds, opt_set_ubounds_value,       &
    &                    opt_name, opt_number, opt_nlevref,           &
    &                    opt_uuid, opt_set_vct_as_levels,             &
    &                    opt_vct)  RESULT(axis)

    TYPE(t_verticalAxis) :: axis

    INTEGER,                 INTENT(IN) :: za_type        !< ICON-internal axis ID (see mo_zaxis_type)
    INTEGER,                 INTENT(IN) :: in_nlevs       !< no. of levels (if no selection)
    REAL(dp),                INTENT(IN) :: levels(:)      !< axis levels
    TYPE(t_level_selection), INTENT(IN), POINTER :: level_selection
    LOGICAL,      INTENT(IN), OPTIONAL :: opt_set_bounds                !< Flag. Set lower/upper bounds if .TRUE.
    REAL(dp),     INTENT(IN), OPTIONAL :: opt_set_ubounds_value         !< Explicit value for ubounds
    CHARACTER(*), INTENT(IN), OPTIONAL :: opt_name                      !< Name of the zaxis
    INTEGER,      INTENT(IN), OPTIONAL :: opt_number                    !< numberOfVGridUsed
    LOGICAL,      INTENT(IN), OPTIONAL :: opt_set_vct_as_levels         !< set VCT to level values
    INTEGER,      INTENT(IN), OPTIONAL :: opt_nlevref                   !< no. of half levels
    REAL(dp),     INTENT(IN), OPTIONAL :: opt_vct(:)                    !< vertical coordinate table
    INTEGER(KIND = C_SIGNED_CHAR), INTENT(IN), OPTIONAL :: opt_uuid(16) !< UUID of vertical grid

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":vertical_axis"

    REAL(dp), ALLOCATABLE :: axis_levels(:), ubounds(:)     !< level list, upper bounds
    INTEGER               :: nlev                           !< no. of axis levels
    LOGICAL               :: set_bounds, select_levs
    INTEGER               :: ierrstat, i, isrc_lev

    ! First, build a level list. If no "level_selection" was
    ! provided, this level lists is simply a copy of the input
    ! parameter "levels(:)":
    select_levs = ASSOCIATED(level_selection)
    IF (.NOT. select_levs) THEN
      nlev = in_nlevs
      IF (nlev <= 0)  CALL finish(routine, "Internal error!")
    ELSE
      ! count the number of feasible levels that have been selected:
      nlev = 0
      DO i=1,level_selection%n_selected
        IF (level_selection%global_idx(i) <= in_nlevs) THEN
          nlev = nlev + 1
        END IF
      END DO
      IF (nlev <= 0)  CALL finish(routine, "Invalid level selection!")
    END IF
    ALLOCATE(axis_levels(nlev),STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    IF (.NOT. select_levs) THEN
      axis_levels(1:nlev) = levels(1:nlev)
    ELSE
      axis_levels(1:nlev) = levels(level_selection%global_idx(1:nlev))
    END IF

    ! create vertical axis object
    axis = t_verticalAxis(zaxisTypeList%getEntry(za_type), nlev)
    CALL axis%set(zaxisLevels=axis_levels) !necessary for GRIB2

    IF (PRESENT(opt_name)) CALL axis%set(zaxisName=TRIM(opt_name))

    set_bounds = .FALSE.
    IF (PRESENT(opt_set_bounds)) set_bounds = opt_set_bounds

    IF (set_bounds) THEN
      ! note: the axes in ICON use the "axis_levels" values as lower bounds
      CALL axis%set(zaxisLbounds=axis_levels) !necessary for GRIB2

      ALLOCATE(ubounds(nlev))
      IF (PRESENT(opt_set_ubounds_value)) THEN
        ubounds(:) = opt_set_ubounds_value
      ELSE
        ubounds(:) = 0._dp
        IF (.NOT. select_levs) THEN
          ubounds(1:nlev) = levels(2:(nlev+1))
        ELSE
          DO i=1,nlev
            isrc_lev = level_selection%global_idx(i)+1
            IF (isrc_lev <= SIZE(levels)) THEN
              ubounds(i) = levels(isrc_lev)
            END IF
          END DO
        END IF
      END IF
      CALL axis%set(zaxisUbounds=ubounds) !necessary for GRIB2
      DEALLOCATE(ubounds)
    END IF

    ! Set numberOfVGridUsed, dependent on the algorithm chosen to
    ! generate the vertical grid (ivctype)
    IF (PRESENT(opt_number))   CALL axis%set(zaxisNumber=opt_number)
    IF (PRESENT(opt_nlevref))  CALL axis%set(zaxisNlevRef=opt_nlevref)
    IF (PRESENT(opt_uuid))     CALL axis%set(zaxisUUID=opt_uuid) ! write vertical grid UUID
    IF (PRESENT(opt_vct))      CALL axis%set(zaxisVct=opt_vct)
    ! Set Vertical Coordinate Table (VCT) to level values
    IF (PRESENT(opt_set_vct_as_levels)) THEN
      IF (opt_set_vct_as_levels)  CALL axis%set(zaxisVct=axis_levels)
    END IF

    DEALLOCATE(axis_levels, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END FUNCTION vertical_axis


  ! --------------------------------------------------------------------------------------
  !> FUNCTION get_numberOfVGridUsed
  !  Depending on the vertical axis chosen for ICON (ivctype), it gives back the value for
  !  the GRIB2-key 'numberOVGridUsed'. Here, we adhere to the COSMO implementation:
  !
  !  |       Description               |  ivctype  |  numberOfVGridUsed  |
  !  =====================================================================
  !  | height based hybrid Gal-Chen    |    1      |       2             |
  !  | height based SLEVE              |    2      |       4             |
  !
  !
  FUNCTION get_numberOfVgridUsed(ivctype)
    INTEGER                 :: get_numberOfVgridUsed
    INTEGER, INTENT(IN)     :: ivctype
    CHARACTER(*), PARAMETER :: routine = modname//":get_numberOfVgridUsed"

    get_numberOfVgridUsed = -1
    SELECT CASE(ivctype)
      CASE(1)
        get_numberOfVgridUsed = 2
      CASE(2)
        get_numberOfVgridUsed = 4
      CASE DEFAULT
        CALL finish(routine, "invalid ivctype! Must be 1 or 2")
    END SELECT
  END FUNCTION get_numberOfVgridUsed

END MODULE mo_name_list_output_zaxes
