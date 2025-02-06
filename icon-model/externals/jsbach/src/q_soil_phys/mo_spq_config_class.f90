!> QUINCY soil-physics process config
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### define soil-physics-quincy config structure, read soil-physics-quincy namelist and init configuration parameters
!>
MODULE mo_spq_config_class
#ifndef __NO_QUINCY__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: real2string
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_spq_config, max_soil_layers

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the spq process, derived from t_jsb_config
  !!
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_spq_config
    LOGICAL                     :: flag_snow            !< on/off snow accumulation
    INTEGER                     :: nsoil_energy         !< number of soil layers for soil energy calculations   ! depreciated feature, but needed for now, see ticket #521
    INTEGER                     :: nsoil_water          !< number of soil layers for soil moisture calculations ! depreciated feature, but needed for now, see ticket #521
    REAL(wp)                    :: soil_depth           !< actual soil and rooting depth
    REAL(wp)                    :: wsr_capacity         !< Water holding capacity of ground skin reservoir [m water equivalent]
    REAL(wp)                    :: wsn_capacity         !< Water holding capacity of ground snow reservoir (max snow depth) [m water equivalent]
    REAL(wp)                    :: soil_awc_prescribe, &
                                   soil_theta_prescribe
    REAL(wp)                    :: soil_sand            !< soil sand proportion - site specific from forcing info file
    REAL(wp)                    :: soil_silt            !< soil silt proportion - site specific from forcing info file
    REAL(wp)                    :: soil_clay            !< soil clay proportion - recalculated from: clay = 1.0 -sand -silt
    REAL(wp)                    :: bulk_density
    CHARACTER(len=filename_max) :: bc_sso_filename          !< elevation and oro_stddev
    CHARACTER(len=filename_max) :: bc_quincy_soil_filename  !< IQ (QUINCY) soil input data
  CONTAINS
    PROCEDURE :: Init => Init_spq_config
  END TYPE t_spq_config

  INTEGER, PARAMETER :: max_snow_layers = 20  ! consistent with JSBACH4
  INTEGER, PARAMETER :: max_soil_layers = 20  ! consistent with JSBACH4

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_spq_config
  !!
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_spq_config(config)
#ifdef __QUINCY_STANDALONE__
    USE mo_namelist,              ONLY: open_nml, POSITIONED, position_nml, close_nml
#else
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
#endif
    USE mo_jsb_model_class,    ONLY: MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC
    USE mo_spq_constants,      ONLY: snow_height_min
    USE mo_jsb_math_constants, ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_spq_config), INTENT(inout) :: config    !< config type for spq
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    LOGICAL                     :: flag_snow
    INTEGER                     :: nsoil_energy    , &
                                   nsoil_water
    REAL(wp)                    :: dz_energy(max_soil_layers), &                        ! vector of soil layer thicknesses
                                   dz_water(max_soil_layers)
    REAL(wp)                    :: soil_layer_profile_ubound_estimate(max_soil_layers)  ! intial estimate of soil profile using the upper bound of each soil layer [m]
    REAL(wp), ALLOCATABLE       :: ubounds_soil_lay_energy(:), &                        ! layer upper bound (note! upper bound > lower bound)
                                   ubounds_soil_lay_water(:)
    REAL(wp)                    :: wsr_capacity    , &
                                   wsn_capacity    , &
                                   soil_depth
    REAL(wp)                    :: soil_awc_prescribe, &
                                   soil_theta_prescribe
    REAL(wp)                    :: soil_sand       , &
                                   soil_silt       , &
                                   soil_clay
    REAL(wp)                    :: bulk_density
    INTEGER                     :: isoil   ! looping
    INTEGER                     :: i       ! looping over snow layers, consistent with jsb4
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    CHARACTER(len=filename_max) :: bc_sso_filename
    CHARACTER(len=filename_max) :: bc_quincy_soil_filename
    REAL(wp)                    :: k_soil_profile       ! for soil layer calculations
    REAL(wp)                    :: min_layer_depth      ! for soil layer calculations
    INTEGER                     :: nsnow
    REAL(wp)                    :: dz_snow(max_snow_layers)   ! width of each snow layer
    REAL(wp), ALLOCATABLE       :: depths(:), mids(:)         ! snow-layer vgrid init
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb              ! soil layers
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_snow_spq             ! snow layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_spq_config'

#ifdef __QUINCY_STANDALONE__
    NAMELIST /spq_ctl/       &
#else
    NAMELIST /jsb_spq_nml/   &
#endif
                          active                  , &
                          ic_filename             , &
                          bc_filename             , &
                          bc_sso_filename         , &
                          bc_quincy_soil_filename , &
                          flag_snow               , &
                          nsoil_energy            , &
                          nsoil_water             , &
                          dz_energy               , &
                          dz_water                , &
                          wsr_capacity            , &
                          wsn_capacity            , &
                          soil_depth              , &
                          soil_awc_prescribe      , &
                          soil_theta_prescribe    , &
                          soil_sand               , &
                          soil_silt               , &
                          soil_clay               , &
                          bulk_density

    ! variables for reading model-options from namelist
    INTEGER :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    IF (config%model_config%model_scheme == MODEL_QUINCY) THEN
      CALL message(TRIM(routine), 'Starting spq configuration')
    ELSE
      CALL message(TRIM(routine), 'NOT starting spq configuration - not running MODEL_QUINCY')
      RETURN
    END IF
    ! ----------------------------------------------------------------------------------------------------- !

    ! variabales for soil layer calculations
    k_soil_profile  = 0.25_wp
    min_layer_depth = 0.065_wp

    ! Set defaults
    active                        = .TRUE.
    ic_filename                   = 'ic_land_spq.nc'
    bc_filename                   = 'bc_land_spq.nc'
    bc_sso_filename               = 'bc_land_sso.nc'
    bc_quincy_soil_filename       = 'bc_quincy_soil.nc'
    flag_snow                     = .TRUE.
    nsoil_energy                  = 15
    nsoil_water                   = 15
    dz_energy(1:5)                = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)
    dz_energy(6:max_soil_layers)  = 0._wp
    dz_water (1:5)                = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)
    dz_water (6:max_soil_layers)  = 0._wp
    wsr_capacity                  = 2.E-4_wp
    wsn_capacity                  = 2.E-4_wp
    soil_depth                    = 9.5_wp
    soil_awc_prescribe            = 300.0_wp
    soil_theta_prescribe          = 1.0_wp
    soil_sand                     = 0.3_wp
    soil_silt                     = 0.4_wp
    soil_clay                     = 0.3_wp
    bulk_density                  = 1500.0_wp

    ! read the namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
#ifdef __QUINCY_STANDALONE__
    nml_unit = position_nml('spq_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, spq_ctl)
#else
    nml_unit = position_nml('jsb_spq_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_spq_nml)
#endif

    CALL close_nml(nml_handler)

    ! check initialisation of number of soil layers
    IF (nsoil_water .LT. 5 .OR. nsoil_water .GT. 20) &
        CALL finish(routine,'Number of soil layers needs to be a number from the range 5 to 20.')
    IF (nsoil_water .NE. nsoil_energy) &
        CALL finish(routine,'Init with different numbers of soil layers for energy and water not yet supported for QUINCY.')

    ! allocate upper-bound soil-layer vectors with dimensions defined in the namelist
    ALLOCATE(ubounds_soil_lay_energy(nsoil_energy))
    ALLOCATE(ubounds_soil_lay_water(nsoil_water))

    ! If not using 5 layers, estimate soil profile (i.e., dz == layer thicknesses) to mimick standard JSB4 profile
    !  assuming dz_water(:) = dz_energy(:)
    IF(nsoil_water > 5) THEN
      ! initial calculation of soil profile (for soil moisture calculations) calculating the upper bound of each layer
      DO isoil=1,nsoil_water
        soil_layer_profile_ubound_estimate(isoil) = min_layer_depth * &
                                                    exp(k_soil_profile * isoil * REAL(max_soil_layers,wp)/REAL(nsoil_water,wp)) - &
                                                    min_layer_depth
      ENDDO
      ! final calculation of soil profile, calculating the layer thickness (dz(:)) and correcting for minimum soil-layer thickness
      dz_water(1) = min_layer_depth
      DO isoil=2,nsoil_water
        dz_water(isoil) = MAX(min_layer_depth, &
                               soil_layer_profile_ubound_estimate(isoil) - soil_layer_profile_ubound_estimate(isoil-1))
      ENDDO
      ! apply the same soil profile to soil profile for soil energy calculations
      dz_energy(:) = dz_water(:)
    ENDIF

    ! run check for soil depth namelist-definition; may not be less/equal zero or larger than covered by soil profile
    IF (soil_depth .LE. 0.0_wp) CALL finish('Init_spq_config','Initialization with non positive soil depth')
    IF (soil_depth .GT. SUM(dz_energy)) &
      & CALL finish('Init_spq_config','Initialization of soil depth larger than soil profile not possible')

    ! calc upper bound of soil layers (water and energy) from layer thickness
    ! note: upper bound is the larger value compared to lower bound of the same layer!
    ubounds_soil_lay_water(1)        = dz_water(1)
    DO isoil=2,nsoil_water
      ubounds_soil_lay_water(isoil)  = ubounds_soil_lay_water(isoil-1) + dz_water(isoil)
    END DO
    ubounds_soil_lay_energy(1)       = dz_energy(1)
    DO isoil=2,nsoil_energy
      ubounds_soil_lay_energy(isoil) = ubounds_soil_lay_energy(isoil-1) + dz_energy(isoil)
    END DO

    ! pass values as read from file
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%bc_sso_filename           = bc_sso_filename
    config%bc_quincy_soil_filename   = bc_quincy_soil_filename
    config%flag_snow                 = flag_snow
    config%nsoil_energy              = nsoil_energy     ! depreciated feature, but needed for now, see ticket #521
    config%nsoil_water               = nsoil_water      ! depreciated feature, but needed for now, see ticket #521
    config%soil_awc_prescribe        = soil_awc_prescribe
    config%soil_theta_prescribe      = soil_theta_prescribe
    config%wsr_capacity              = wsr_capacity
    config%wsn_capacity              = wsn_capacity
    config%soil_depth                = soil_depth
    config%soil_sand                 = soil_sand
    config%soil_silt                 = soil_silt
    !config%soil_clay                = not using the input value, as soil_clay must be calculated from soil_sand & soil_silt to guarantee: 1=sand+silt+clay
    config%soil_clay                 = 1.0_wp - soil_sand - soil_silt
    config%bulk_density              = bulk_density

    ! check for an error in clay content of soil
    IF(config%soil_clay < eps8) THEN
       WRITE(message_text,'(a)') 'invalid proportion of soil_clay (value < eps8)'
       CALL finish(routine, message_text)
    END IF

    !< Create vertical snow-layer axis (vgrid) - consistent with jsb4
    !!
    !! these are "infrastructure" values \n
    !! layer 1 is the one at the ground \n
    !! snow_lay_thickness_snl(:,:) defines how much each of the layers if filled with snow
    !!
    nsnow             = 5       ! by default 5 snow layer (in jsb4 it is a namelist value, which is actually not used)
    dz_snow(:)        = 0._wp   ! jsb4 default
    dz_snow(1:nsnow)  = 0.05_wp ! jsb4 default
    IF (dz_snow(1) < snow_height_min) THEN
      CALL finish(TRIM(routine), 'Depth of first snow layer should be larger than snow_height_min=' &
      &                            //real2string(snow_height_min))
    END IF
    ALLOCATE(depths(nsnow+1))
    ALLOCATE(mids(nsnow))
    depths(1) = 0._wp
    DO i=1,nsnow
      depths(i+1) = depths(i) + dz_snow(i)
    END DO
    mids(1:nsnow) = (depths(1:nsnow) + depths(2:nsnow+1)) / 2._wp

    vgrid_snow_spq  => new_vgrid('snow_layer_spq', ZAXIS_GENERIC, nsnow,     &
      & levels  = mids                 (1:nsnow  ),                         &
      & lbounds = depths               (1:nsnow  ),                         &
      & ubounds = depths               (2:nsnow+1),                         &
      & units='m')
    CALL register_vgrid(vgrid_snow_spq)
    WRITE(message_text, *) 'Snow layer SPQ_ (upper)  [m]: ', vgrid_snow_spq%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow layer SPQ_ (mid)    [m]: ', vgrid_snow_spq%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow layer SPQ_ (lower)  [m]: ', vgrid_snow_spq%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Snow layer thickness SPQ_ [m]: ', vgrid_snow_spq%dz
    CALL message(TRIM(routine), message_text)
    DEALLOCATE(depths, mids)

    !< Create vertical soil-layer axis (vgrid)
    !!
    !! these are "infrastructure" values, ubounds and dz must be corrected/limited by site-specific bedrock depth "soil_depth" !
    !!
    !! site-specific soil-layer thickness values are stored in soil_depth_sl (and calc in spq_init) \n
    !! the function new_vgrid() does: (a) check for 'dz(:) <= zero', and (b) calc levels & lbounds from ubounds and dz \n
    !! levels is defined as the depth at the center of the layer: levels(:) = 0.5_wp * (lbounds(:) + ubounds(:))
    vgrid_soil_sb => new_vgrid('soil_layer_sb', ZAXIS_GENERIC, nsoil_water, &
                      longname='Soil (water) layers from SPQ_ QUINCY', &
                      units='m', &
                      ! levels=  , &      ! is calculated from ubounds and dz
                      ! lbounds=  , &     ! is calculated from ubounds and dz
                      ubounds=ubounds_soil_lay_water(:), &
                      dz=dz_water(:))
    CALL register_vgrid(vgrid_soil_sb)
    WRITE(message_text, *) 'Soil layer SPQ_ (upper)   [m]: ', vgrid_soil_sb%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil layer SPQ_ (mid)     [m]: ', vgrid_soil_sb%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil layer SPQ_ (lower)   [m]: ', vgrid_soil_sb%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil layer thickness SPQ_ [m]: ', vgrid_soil_sb%dz
    CALL message(TRIM(routine), message_text)
  END SUBROUTINE Init_spq_config

#endif
END MODULE mo_spq_config_class
