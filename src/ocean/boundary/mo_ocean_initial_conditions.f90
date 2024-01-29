! Contains the implementation of the initial conditions for the hydrostatic ocean model.
!
! Contains the implementation of the initial conditions for the hydrostatic ocean model.
! This module controls the initial conditions as well as the initialisation of
! test cases, the top and bottom boundary conditions, and the structure of the
! forcing quantities.
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

MODULE mo_ocean_initial_conditions
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_grid_config,        ONLY: grid_sphere_radius, grid_angular_velocity
  USE mo_physical_constants, ONLY: rgrav, tmelt, tf, earth_angular_velocity,inverse_earth_radius! , SItodBar
  USE mo_math_constants,     ONLY: pi, pi_2, rad2deg, deg2rad
  USE mo_parallel_config,    ONLY: nproma
  USE mo_ocean_nml,          ONLY: iswm_oce, n_zlev, no_tracer, i_sea_ice,            &
    & basin_center_lat, basin_center_lon, basin_height_deg,  basin_width_deg,         &
    & initial_temperature_bottom, initial_temperature_top, initial_temperature_shift, &
    & initial_temperature_north, initial_temperature_south,                           &
    & initial_temperature_scale_depth, initial_temperature_VerticalGradient,         &
    & type_3dimrelax_temp, type_3dimrelax_salt,                                      &
    & initial_salinity_top, initial_salinity_bottom, &
    & topography_type, topography_height_reference,  &
    & sea_surface_height_type, initial_temperature_type, initial_salinity_type, &
    & initial_sst_type, initial_velocity_type, initial_velocity_amplitude,      &
    & forcing_temperature_poleLat, InitialState_InputFileName,                  &
    & smooth_initial_height_weights, smooth_initial_salinity_weights,           &
    & smooth_initial_temperature_weights, &
    & initial_perturbation_waveNumber, initial_perturbation_max_ratio,          &
    & smooth_initial_salinity_iterations, &
    & smooth_initial_height_iterations, smooth_initial_temperature_iterations,  &
    & OceanReferenceDensity, LinearThermoExpansionCoefficient,                  &
    & smooth_initial_velocity_iterations, smooth_initial_velocity_weights,      &
    & vert_cor_type
    
  USE mo_sea_ice_nml,        ONLY: use_IceInitialization_fromTemperature

  USE mo_impl_constants,     ONLY: sea_boundary
  USE mo_dynamics_config,    ONLY: nold
  USE mo_math_types,         ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_math_utilities,     ONLY: gvec2cvec
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_model_domain,       ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_types,          ONLY: t_hydro_ocean_state
  USE mo_scalar_product,     ONLY: map_cell2edges_3D
  USE mo_ocean_math_operators,ONLY: smooth_onCells
  USE mo_ape_params,         ONLY: ape_sst
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range

  USE mo_ocean_nudging,       ONLY: ocean_nudge
  
  USE mo_sync,              ONLY: sync_c, sync_e, sync_patch_array
  USE mo_fortran_tools,     ONLY: assign_if_present
  
  USE mo_read_interface,    ONLY: read_2D_1time, read_3D_1time, on_cells, on_edges, t_stream_id, &
    & read_netcdf_broadcast_method, openInputFile, closeFile
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: apply_initial_conditions, init_ocean_bathymetry !,&
  PUBLIC :: tracer_ConstantSurface, varyTracerVerticallyExponentially
  PUBLIC :: init_cell_2D_variable_fromFile, init_cell_3D_variable_fromFile
!   & SST_LinearMeridional, increaseTracerLevelsLinearly
  
  REAL(wp) :: sphere_radius, u0
  REAL(wp), PARAMETER :: aleph = 0.0_wp
  
  CHARACTER(LEN=12), PARAMETER :: module_name = 'oceInitCond'

  TYPE(t_operator_coeff), POINTER :: this_operators_coeff
  TYPE(t_hydro_ocean_state), POINTER :: this_ocean_state

  ! Should be replaced by reading a file
  REAL(wp), PARAMETER :: tprof(20)=&
    & (/ 18.13_wp, 17.80_wp, 17.15_wp, 16.09_wp, 15.04_wp, 13.24_wp, 11.82_wp,  9.902_wp, &
    & 8.484_wp, 7.341_wp, 5.727_wp, 4.589_wp, 3.807_wp, 3.062_wp, 2.481_wp, 2.194_wp, &
    & 1.789_wp, 1.266_wp, 1.070_wp, 0.9211_wp /)

    ! temperature profile for 4-20 layers for testcase 40/45 and similar
    REAL(wp), PARAMETER :: tprof_var(20)= &
      & (/ 25.0_wp, 23.0_wp, 20.0_wp, 15.0_wp, 10.0_wp, 8.0_wp, 6.0_wp, 5.0_wp, 4.0_wp, 3.0_wp,&
      & 2.0_wp,  1.0_wp,  0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)
    !   & (/ 25.0_wp, 18.0_wp, 12.0_wp, 8.0_wp, 6.0_wp, 4.0_wp, 2.0_wp, 1.0_wp, 0.5_wp, 0.0_wp,&
    !   &     0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/)

    REAL(wp) , PARAMETER :: tprof_4layerstommel(4) = (/20.0_wp,10.0_wp,8.0_wp,6.0_wp/)

  
CONTAINS
  

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE apply_initial_conditions(patch_3d, ocean_state, external_data, &
    & operators_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_external_data)                   :: external_data
    TYPE(t_operator_coeff), INTENT(in), TARGET          :: operators_coeff

!   TYPE(t_patch),POINTER                   :: patch_2d
!     REAL(wp), ALLOCATABLE                   :: check_temp(:,:,:), check_salinity(:,:,:)

    this_operators_coeff => operators_coeff
    ! patch_2d => patch_3d%p_patch_2d(1)
    sphere_radius = grid_sphere_radius
    u0 = (2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)
    
    ! the bathymetry initialization is called  after read_external_data and before seting up the sea-land mask
    !CALL init_ocean_bathymetry(patch_3d=patch_3d,  cells_bathymetry=external_data%oce%bathymetry_c(:,:))
    IF(iswm_oce==1)CALL init_ocean_bathymetry(patch_3d=patch_3d,  cells_bathymetry=external_data%oce%bathymetry_c(:,:))
    
    CALL init_ocean_velocity(patch_3d=patch_3d, normal_velocity=ocean_state%p_prog(nold(1))%vn)


    IF (vert_cor_type .EQ. 0) THEN
      CALL init_ocean_surface_height(patch_3d=patch_3d, ocean_height=ocean_state%p_prog(nold(1))%h(:,:))
    ELSE
      CALL init_ocean_surface_height(patch_3d=patch_3d, ocean_height=ocean_state%p_prog(nold(1))%eta_c(:,:))
    END IF

    IF (no_tracer > 0) &
      & CALL init_ocean_temperature(patch_3d=patch_3d, ocean_temperature=ocean_state%p_prog(nold(1))%tracer(:,:,:,1),&
      &ocean_state=ocean_state)
      
    IF (no_tracer > 1) &
      & CALL init_ocean_salinity(patch_3d=patch_3d, ocean_salinity=ocean_state%p_prog(nold(1))%tracer(:,:,:,2))

    !---------------------------------------------------------------------
    ! CALL initialize_diagnostic_fields( patch_2d, patch_3d, ocean_state, operators_coeff)
    ! CALL fill_tracer_x_height(patch_3d, ocean_state)

  END SUBROUTINE apply_initial_conditions
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE init_cell_3D_variable_fromFile(patch_3d, variable, name, has_missValue, missValue)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), INTENT(inout) :: variable(:,:,:)
    CHARACTER(LEN=*) :: name
    LOGICAL  :: has_missValue
    REAL(wp) :: missValue

    TYPE(t_patch),POINTER  :: patch_2d
    TYPE(t_stream_id) :: stream_id
    INTEGER :: block, idx, start_cell_index, end_cell_index, level
    TYPE(t_subset_range), POINTER :: all_cells
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_cell_3D_variable_fromFile '
    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL


    CALL message (method_name, TRIM(name)//"...")
    ! read temperature
    CALL openInputFile(stream_id, initialState_InputFileName, patch_2d)

    CALL read_3D_1Time( stream_id=stream_id, location=on_cells, &
      & variable_name=name, fill_array=variable,               &
      & has_missValue=has_missValue, missValue=missValue)


    CALL closeFile(stream_id)

    ! write(0,*) variable
!     CALL work_mpi_barrier()
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        DO level = patch_3d%p_patch_1d(1)%dolic_c(idx,block) + 1, n_zlev
!           IF ( variable(idx,level,block) /=  0.0_wp) THEN
!             CALL warning(method_name, "non-zero variable on land")
            variable(idx,level,block) = 0.0_wp
!           ENDIF
        ENDDO
      ENDDO
    ENDDO
    
    CALL sync_patch_array(sync_c, patch_2D, variable)

    WRITE(message_text,*) "miss value=", missValue, has_missValue
    CALL message(method_name, message_text)

  END SUBROUTINE init_cell_3D_variable_fromFile
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE init_cell_2D_variable_fromFile(patch_3d, variable, name, has_missValue, missValue)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), INTENT(inout) :: variable(:,:)
    CHARACTER(LEN=*) :: name
    LOGICAL  :: has_missValue
    REAL(wp) :: missValue
    
    INTEGER :: block, idx, start_cell_index, end_cell_index, level
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_stream_id) :: stream_id
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_cell_2D_variable_fromFile '
    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    
    CALL message (method_name, TRIM(name)//"...")
    ! read temperature
    !  - 2011-11-01, >r7005: read one data set, annual mean only
    !  - "T": annual mean temperature
    ! ram: the input has to be POTENTIAL TEMPERATURE!
    CALL openinputfile(stream_id, initialState_InputFileName, patch_2d, &
      &                read_netcdf_broadcast_method)
    
    CALL read_2D_1Time( stream_id=stream_id, location=on_cells, &
      & variable_name=name, fill_array=variable,               &
      & has_missValue=has_missValue, missValue=missValue)

    CALL closeFile(stream_id)

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        ! FIXME: to complicated, better use land sea mask
        DO level = patch_3d%p_patch_1d(1)%dolic_c(idx,block) + 1, 1
!          IF ( variable(idx,block) /=  0.0_wp) THEN
!            CALL warning(method_name, "non-zero variable on land")
            variable(idx,block) = 0.0_wp
!          ENDIF
        ENDDO
      ENDDO
    ENDDO

    CALL sync_patch_array(sync_c, patch_2D, variable)
    
    WRITE(message_text,*) "miss value=", missValue, has_missValue
    CALL message(method_name, message_text)
   
  END SUBROUTINE init_cell_2D_variable_fromFile
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE init_edge_3D_variable_fromFile(patch_3d, variable, name, has_missValue, missValue)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), INTENT(inout) :: variable(:,:,:)
    CHARACTER(LEN=*) :: name
    LOGICAL  :: has_missValue
    REAL(wp) :: missValue

    TYPE(t_patch),POINTER  :: patch_2d
    TYPE(t_stream_id) :: stream_id
    INTEGER :: block, idx, start_edge_index, end_edge_index, level
    TYPE(t_subset_range), POINTER :: all_edges
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_edge_3D_variable_fromFile '
    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL


    CALL message (method_name, TRIM(name)//"...")
    ! read temperature
    CALL openInputFile(stream_id, initialState_InputFileName, patch_2d)

    CALL read_3D_1time( stream_id=stream_id, location=on_edges, &
      & variable_name=name, fill_array=variable,               &
      & has_missValue=has_missValue, missValue=missValue)

    CALL closeFile(stream_id)
    ! write(0,*) variable
!     CALL work_mpi_barrier()
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_edge_index, end_edge_index)
      DO idx = start_edge_index, end_edge_index
        DO level = patch_3d%p_patch_1d(1)%dolic_e(idx,block) + 1, n_zlev
!           IF ( variable(idx,level,block) /=  0.0_wp) THEN
!             CALL warning(method_name, "non-zero variable on land")
            variable(idx,level,block) = 0.0_wp
!           ENDIF
        ENDDO
      ENDDO
    ENDDO
    
    CALL sync_patch_array(sync_e, patch_2D, variable)
  
  END SUBROUTINE init_edge_3D_variable_fromFile
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE init_ocean_bathymetry(patch_3d, cells_bathymetry)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET  :: cells_bathymetry(:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_bathymetry'
    !-------------------------------------------------------------------------

    IF (topography_type < 200) RETURN ! not analytic bathymetry


    SELECT CASE (topography_type)
    CASE (200)
      ! constant depth given by topography_height_reference
      CALL depth_uniform(patch_3d, cells_bathymetry)

    CASE (201)
      CALL depth_mountain_orography_Williamson_test5(patch_3d, cells_bathymetry)

    CASE (202)
      CALL depth_overflow(patch_3d, cells_bathymetry)

    CASE default
      CALL finish(method_name, "unknown topography_type")
    END SELECT

  END SUBROUTINE init_ocean_bathymetry
  !-------------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_salinity(patch_3d, ocean_salinity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    ! replace by file
    REAL(wp), PARAMETER :: sprof(20)=&
      & (/  34.699219_wp, 34.798244_wp, 34.904964_wp, 34.976841_wp, 35.027084_wp, &
      & 35.026825_wp, 34.960835_wp, 34.862324_wp, 34.752468_wp, 34.656761_wp, 34.596603_wp,&
      & 34.594128_wp, 34.628601_wp, 34.678772_wp, 34.717495_wp, 34.738304_wp, 34.741512_wp,&
      & 34.738205_wp, 34.729176_wp, 34.723465_wp /)

    REAL(wp), PARAMETER :: salinity_profile_20levels(20)= &
      & (/ 34.5_wp, 34.6_wp, 34.7_wp, 34.8_wp, 34.9_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp,&
      & 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp, 35.0_wp/)

    REAL(wp) , PARAMETER :: sprof_4layerstommel(4) = &
      & (/34.699219_wp, 34.798244_wp, 34.904964_wp, 34.976841_wp/)

    LOGICAL  :: has_missValue
    REAL(wp) :: missValue
    REAL(wp), ALLOCATABLE :: old_salinity(:,:,:)
    INTEGER :: i
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_salinity'
    !-------------------------------------------------------------------------

    IF (no_tracer < 2) RETURN        ! no salinity

    ocean_salinity(:,:,:) = 0.0_wp
    has_missValue = .false.
    missValue     = -99999999.0_wp
!     IF (initial_salinity_type < 200) RETURN ! not analytic salinity

    SELECT CASE (initial_salinity_type)
    
    CASE (000)
      CALL message(method_name, ' no initialization')

    CASE (001)
      CALL message(method_name, ': init from file')
      CALL init_cell_3D_variable_fromFile(patch_3d, variable=ocean_salinity, name="S", &
        & has_missValue=has_missValue, missValue=missValue)

    ! copy the initial ocean salinity also to the nudging  salinity
      IF (no_tracer>1 .AND. type_3dimrelax_salt >0) &
         ocean_nudge%data_3dimRelax_Salt(:,:,:) = ocean_salinity(:,:,:)
    !------------------------------
    CASE (200)
      ! uniform salinity or vertically linarly increasing
      CALL message(method_name, ': horizontally homogenous, vertically linear')
      CALL tracer_ConstantSurface(patch_3d=patch_3d, ocean_tracer=ocean_salinity, &
        & top_value=initial_salinity_top)
      CALL increaseTracerVerticallyLinearly(patch_3d=patch_3d, ocean_tracer=ocean_salinity,&
        & bottom_value=initial_salinity_bottom)

    !------------------------------
    CASE (201)
      CALL salinity_Uniform_SpecialArea(patch_3d, ocean_salinity)

    !------------------------------
    CASE (202)
      ! simulates the salinity_profile_20levels array, but for any number of levels
      CALL salinity_AnalyticSmoothVerticalProfile(patch_3d, ocean_salinity)

    !------------------------------
    CASE (212)
      CALL tracer_smoothAPE_LinearDepth(patch_3d, ocean_salinity, &
        & top_value=initial_salinity_top, bottom_value=initial_salinity_bottom)

    !------------------------------
    CASE (228) ! salinity dome 
      CALL salinity_GM_idealized(patch_3d,ocean_salinity) 

    CASE (229)  ! linear salinity slope
      CALL salinity_GM_idealized2(patch_3d,ocean_salinity)  

    CASE (230)  ! 2d salinity blubb
      CALL salinity_GM_idealized3(patch_3d,ocean_salinity)  

    CASE(232) ! horizontal constant 
      CALL salinity_Willebrand_test(patch_3d,ocean_salinity)
      ocean_salinity=34.0_wp+0.50_wp*ocean_salinity
 
    CASE (235)  ! const 35
      ocean_salinity(:,:,:) = 35.0_wp 

    CASE (300)
      CALL tracer_bubble(patch_3d, ocean_salinity ,initial_salinity_top, initial_salinity_bottom)

    CASE (301)
      CALL tracer_bubble_disturbed(patch_3d, ocean_salinity ,initial_salinity_top, initial_salinity_bottom)


    CASE (302)
      CALL tracer_double_bubble(patch_3d, ocean_salinity ,initial_salinity_top, initial_salinity_bottom)

    !------------------------------
    CASE (401)
      ! assign from adhoc array values
      IF (n_zlev==4) THEN
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_salinity, VerticalProfileValue=sprof_4layerstommel)
      ELSEIF  (n_zlev <= 20) THEN
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_salinity, VerticalProfileValue=sprof)
      ELSE
        CALL finish(method_name, 'Number of vertical levels to small or to big: >=4 and <=20')
      ENDIF

      
    CASE (402)
      IF  (n_zlev <= 20) THEN
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_salinity, VerticalProfileValue=salinity_profile_20levels)
      ELSE
        CALL finish(method_name, 'Number of vertical levels > 20')
      ENDIF

    !------------------------------
!    CASE (500)
!      ! uniform salinity or vertically linarly increasing on all cells
!      !   it will initialize also land cells
!      !   this is provided only for testing,
!      !   as it should give the same results as
!      !   for initial_salinity_type = 200
!      CALL message(method_name, ': horizontally homogenous, vertically linear INCLUDING LAND')
!      CALL tracer_ConstantSurface_IncludeLand(patch_3d=patch_3d, ocean_tracer=ocean_salinity, &
!        & top_value=initial_salinity_top, bottom_value=initial_salinity_bottom)

    !------------------------------
     CASE default
      CALL finish(method_name, "unknown initial_salinity_type")

    END SELECT
    
    IF (smooth_initial_salinity_iterations > 0) THEN
      ALLOCATE(old_salinity(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks))
      DO i=1,smooth_initial_salinity_iterations
        CALL message(method_name, "Smoothing initial salinity...")
        old_salinity = ocean_salinity
        CALL smooth_onCells(patch_3D=patch_3d, &
          & in_value=old_salinity, out_value=ocean_salinity, &
          & smooth_weights=smooth_initial_salinity_weights, &
          & has_missValue=has_missValue, missValue=missValue)
      ENDDO
      DEALLOCATE(old_salinity)
    ENDIF
    
    CALL fillVerticallyMissingValues(patch_3d=patch_3d, ocean_tracer=ocean_salinity,&
      & has_missValue=has_missValue, missValue=missValue)

    CALL dbg_print('init_ocean_salinity', ocean_salinity(:,:,:), &
      & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)

  END SUBROUTINE init_ocean_salinity
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_temperature(patch_3d, ocean_temperature, ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET                        :: ocean_temperature(:,:,:)
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state

    LOGICAL  :: has_missValue
    REAL(wp) :: missValue
    REAL(wp), ALLOCATABLE :: old_temperature(:,:,:)
    INTEGER ::i
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_temperature'
    !-------------------------------------------------------------------------

    IF (no_tracer < 1) RETURN        ! no temperature

    has_missValue = .false.
    missValue     = -99999999.0_wp
    
    SELECT CASE (initial_temperature_type)
    !------------------------------
    CASE (000)
    
      ocean_temperature(:,:,:) = 0.0_wp
    
      CALL message(method_name, ' zero initialization')

    CASE (001)
      CALL message(method_name, ': init from file')
      !  - "T": annual mean temperature
      ! ram: the input has to be POTENTIAL TEMPERATURE!
      CALL init_cell_3D_variable_fromFile(patch_3d, variable=ocean_temperature, name="T", &
        & has_missValue=has_missValue, missValue=missValue)
      use_IceInitialization_fromTemperature = .true. ! this should be set in the namelist, here only for safety

      IF (no_tracer>=1 .AND. type_3dimrelax_temp >0) &
         ocean_nudge%data_3dimRelax_temp(:,:,:) = ocean_temperature(:,:,:)


    !------------------------------
    CASE (200)
      ! uniform or linearly decreasing temperature
      ! Temperature is homogeneous in each layer.
      CALL message(method_name, ': horizontally homogenous, vertically linear')
      CALL tracer_ConstantSurface(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & top_value=initial_temperature_top)
        
      CALL increaseTracerVerticallyLinearly(patch_3d=patch_3d, ocean_tracer=ocean_temperature,&
        & bottom_value=initial_temperature_bottom)

    !------------------------------
    CASE (201)
      CALL temperature_CollapsingDensityFront_StuhnePeltier(patch_3d, ocean_temperature)

    !------------------------------
    CASE (202)
      CALL temperature_BasinWithVerticalWall(patch_3d, ocean_temperature)

    !------------------------------
    CASE (203)
      CALL temperature_DanilovsMunkGyre(patch_3d, ocean_temperature)

    !------------------------------
    CASE (204)
      CALL temperature_TropicsPolar(patch_3d, ocean_temperature)

    !------------------------------
    CASE (205)
      CALL temperature_CollapsingDensityFront_WeakGrad(patch_3d, ocean_temperature)

    !------------------------------
    CASE (206)
      CALL temperature_Uniform_SpecialArea(patch_3d, ocean_temperature)

    !------------------------------
    CASE (207)
      CALL temperature_APE(patch_3d, ocean_temperature)

    !------------------------------
    CASE (208)
      CALL message(method_name, ': horizontally non-homogenous, local pertubation')

      ! first create linearly vertically decreasing temperature, uniform horizontally
      CALL tracer_ConstantSurface(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & top_value=initial_temperature_top)

      CALL increaseTracerVerticallyLinearly(patch_3d=patch_3d, ocean_tracer=ocean_temperature,&
        & bottom_value=initial_temperature_bottom)
        
      !Add horizontal variation
      CALL temperature_AddHorizontalVariation(patch_3d, ocean_temperature)

      !Add local perturbation
      CALL temperature_AddLocalPerturbation(patch_3d, ocean_temperature)

    !------------------------------
    CASE (209)
      CALL temperature_uniform_SeparationAtLon(patch_3d, ocean_temperature, wallLonDeg=0.0_wp)

    !------------------------------
    CASE (210)
      CALL temperature_uniform_SeparationAtLat(patch_3d, ocean_temperature, wallLatDeg=basin_center_lat)

    !------------------------------
    CASE (211)
      CALL temperature_circularLonLatPerturbation(patch_3d, ocean_temperature)

    !------------------------------
    CASE (212)
      CALL tracer_smoothAPE_LinearDepth(patch_3d, ocean_temperature, &
        & top_value=initial_temperature_top, bottom_value=initial_temperature_bottom)

    !------------------------------
    CASE (213)
      CALL temperature_smoothAPE_LinearLevels(patch_3d, ocean_temperature)

    !------------------------------
    CASE (214)
      ! not used
      CALL finish(method_name, "214 is not used any more")
!       CALL SST_LinearMeridional(patch_3d, ocean_temperature)
!       !  exponential temperature profile following Abernathey et al., 2011
!       ! CALL varyTracerVerticallyExponentially(patch_3d, ocean_temperature, initial_temperature_bottom, &
!       !   &                                    initial_temperature_scale_depth)
    CASE (215)
      ! not used
      CALL finish(method_name, "215 is not used any more")
!      CALL SST_LinearMeridional(patch_3d, ocean_temperature)
!      !  exponential temperature profile following Abernathey et al., 2011
!      CALL increaseTracerLevelsLinearly(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
!        & bottom_value=initial_temperature_bottom)
!      !CALL temperature_AddSinusoidalPerturbation(patch_3d, ocean_temperature)

     CASE (216)

      CALL temperature_front(patch_3d, ocean_temperature)
          
    CASE (217)
      CALL SST_LinearMeridional(patch_3d, ocean_temperature)
      CALL increaseTracerLevelsLinearly(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & increase_gradient=initial_temperature_VerticalGradient)
      CALL perturbeTracer_LonCosinus(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & waveNumber=initial_perturbation_waveNumber, max_ratio=initial_perturbation_max_ratio)
      CALL perturbeTracer_LatCosinus(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & waveNumber=1.0_wp * initial_perturbation_waveNumber, &
        &  max_ratio=0.1_wp * initial_perturbation_max_ratio)

    CASE (218)
      CALL SST_LinearMeridional(patch_3d, ocean_temperature)
      CALL increaseTracerLevelsLinearly(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & bottom_value=initial_temperature_bottom)
      CALL perturbeTracer_LonCosinus(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & waveNumber=initial_perturbation_waveNumber, max_ratio=initial_perturbation_max_ratio)
      CALL perturbeTracer_LatCosinus(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
        & waveNumber=1.0_wp * initial_perturbation_waveNumber, &
        &  max_ratio=0.1_wp * initial_perturbation_max_ratio)

    !CASE (220)
    
    !  CALL tracer_Redi_test(patch_3d=patch_3d, ocean_tracer=ocean_temperature,ocean_state=ocean_state)

    CASE (221)
      ! Abernathey setup 01; initial SST reflects the heat fluxes
      CALL SST_Abernathey_01(patch_3d=patch_3d, ocean_temperature=ocean_temperature, &
!       & BaseTemperature=initial_temperature_top * 0.5_wp, &
        & VariationAmplitude=initial_temperature_south, &
        & VariationLength = basin_height_deg * deg2rad, &
        & VariationWaveNo=2.5_wp, &
        & SouthLat=(basin_center_lat - 0.5_wp * basin_height_deg) * deg2rad)
            
      CALL varyTracerVerticallyExponentially(patch_3d, ocean_temperature, initial_temperature_bottom, &
        &                                    initial_temperature_scale_depth)

    CASE (222)
      ! Abernathey setup 02; initial SST ls linear
      CALL SST_LinearMeridional(patch_3d, ocean_temperature)

      CALL varyTracerVerticallyExponentially(patch_3d, ocean_temperature, initial_temperature_bottom, &
        &                                    initial_temperature_scale_depth)

    CASE(223)
      CALL temperature_dirac_signal(patch_3d, ocean_temperature)

    CASE(224)
      CALL tracer_Redi_test(patch_3d=patch_3d, ocean_tracer=ocean_temperature,ocean_state=ocean_state)      

    CASE(225)
      CALL tracer_GM_test(patch_3d=patch_3d, ocean_tracer=ocean_temperature,ocean_state=ocean_state)      

    CASE(226)
      CALL tracer_Redi_test2(patch_3d=patch_3d, ocean_tracer=ocean_temperature,ocean_state=ocean_state)      
    CASE(227)
      CALL tracer_GMR_slope_test(patch_3d=patch_3d, ocean_tracer=ocean_temperature,ocean_state=ocean_state)
      
    CASE(228) ! temperature dome
      CALL temperature_GM_idealized(patch_3d,ocean_temperature)  

    CASE(229) ! horizontal constant 
      CALL temperature_GM_idealized2(patch_3d,ocean_temperature)  

    CASE(230) ! horizontal constant 
      CALL temperature_GM_idealized3(patch_3d,ocean_temperature)

    CASE(231) ! horizontal constant 
      CALL temperature_GM_idealized4(patch_3d,ocean_temperature)

    CASE(232) ! horizontal constant 
      CALL temperature_Willebrand_test(patch_3d,ocean_temperature)
      ocean_temperature=20.0_wp-10.0_wp*ocean_temperature
    
    CASE (235)
      ocean_temperature(:,:,:) = 10.0_wp

    !------------------------------
    !! Test cases for zstar
    !------------------------------
    CASE(240)
      CALL tracer_advec_test(patch_3d, ocean_temperature)

    CASE(241)
      CALL lock_exchange(patch_3d, ocean_temperature)
    
    CASE(242)
      CALL overflow(patch_3d, ocean_temperature)
    !------------------------------
    !------------------------------


    CASE(300)
     CALL message(method_name, 'Temperature Kelvin-Helmholtz Test ')
     CALL temperature_KelvinHelmholtzTest(patch_3d, ocean_temperature,&
     &top_value=initial_temperature_top,bottom_value=initial_temperature_bottom )

    !------------------------------
    CASE (307)
      CALL tracer_bubble(patch_3d, ocean_temperature ,initial_temperature_top, initial_temperature_bottom)
    !------------------------------
    CASE (301)
      CALL tracer_bubble_disturbed(patch_3d, ocean_temperature ,initial_temperature_top, initial_temperature_bottom)

    !------------------------------
    CASE (302)
      CALL tracer_double_bubble(patch_3d, ocean_temperature ,initial_temperature_top, initial_temperature_bottom)

    CASE (303)
      CALL tracer_layer(patch_3d, ocean_temperature ,initial_temperature_top, initial_temperature_bottom)

    CASE (304) 
      CALL tracer_bubbles_side_by_side(patch_3d, ocean_temperature ,initial_temperature_top, initial_temperature_bottom)

    CASE (305) 
      CALL Roberts_tracer_bubble(patch_3d, ocean_temperature ,initial_temperature_top, initial_temperature_bottom)

    CASE (306) 
      CALL inclined_layer(patch_3d, ocean_temperature ,initial_temperature_top, initial_temperature_bottom)

    !------------------------------
    CASE (401)
      ! assign from adhoc array values
      IF(n_zlev==4)THEN
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_temperature, VerticalProfileValue=tprof_4layerstommel)
      ELSE
        CALL fill_FromVerticalArrayProfile(patch_3d, ocean_temperature, VerticalProfileValue=tprof)
      ENDIF

    !------------------------------
    CASE (402)
       CALL fill_FromVerticalArrayProfile(patch_3d, ocean_temperature, VerticalProfileValue=tprof_var)

    !------------------------------
!    CASE (500)
!      ! uniform salinity or vertically linarly increasing on all cells
!      !   it will initialize also land cells
!      !   this is provided only for testing,
!      !   as it should give the same results as
!      !   for initial_salinity_type = 200
!      CALL message(method_name, ': horizontally homogenous, vertically linear INCLUDING LAND')
!      CALL tracer_ConstantSurface_IncludeLand(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
!        & top_value=initial_temperature_top, bottom_value=initial_temperature_bottom)

     !------------------------------
     CASE default
      CALL finish(method_name, "unknown initial_temperature_type")

    END SELECT
    
    IF (smooth_initial_temperature_iterations > 0) THEN
      ALLOCATE(old_temperature(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks))
      DO i=1,smooth_initial_temperature_iterations
        CALL message(method_name, "Smoothing temperature...")
        old_temperature = ocean_temperature
        CALL smooth_onCells(patch_3D=patch_3d, &
          & in_value=old_temperature, out_value=ocean_temperature,  &
          & smooth_weights=smooth_initial_temperature_weights,      &
          & has_missValue=has_missValue, missValue=missValue)
      ENDDO
      DEALLOCATE(old_temperature)
    ENDIF

    CALL fillVerticallyMissingValues(patch_3d=patch_3d, ocean_tracer=ocean_temperature, &
      & has_missValue=has_missValue, missValue=missValue)
      
    CALL dbg_print('init_ocean_temperature', ocean_temperature(:,:,:), &
      & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)

  END SUBROUTINE init_ocean_temperature
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_velocity_uv_fromFile(patch_3d, normal_velocity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: normal_velocity(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_cartesian_coordinates), ALLOCATABLE ::cellVelocity_cc(:,:,:)
    REAL(wp), ALLOCATABLE :: u(:,:,:), v(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET :: suv(:,:,:)
    INTEGER :: i, level
    INTEGER :: cell_block, cell_index
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: x1, x2, x3, lon, lat
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    LOGICAL  :: has_missValue
    REAL(wp) :: missValue

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_velocity_uv_fromFile'
    !-------------------------------------------------------------------------

    ! CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    ! fisrt calculate th velocyt at cell centers
    ALLOCATE(u(nproma,n_zlev, patch_2d%alloc_cell_blocks), &
             v(nproma,n_zlev, patch_2d%alloc_cell_blocks))
!     ALLOCATE(u(nproma,n_zlev, patch_2d%alloc_cell_blocks), &
!              v(nproma,n_zlev, patch_2d%alloc_cell_blocks))
    !-------------------------------------------------------------------------
    CALL message(method_name, ': init from file')
    CALL init_cell_3D_variable_fromFile(patch_3d, variable=u, name="u", &
      & has_missValue=has_missValue, missValue=missValue)
    CALL init_cell_3D_variable_fromFile(patch_3d, variable=v, name="v", &
      & has_missValue=has_missValue, missValue=missValue)
      
    IF (smooth_initial_velocity_iterations > 0) THEN
      ALLOCATE(suv(nproma,n_zlev, patch_2d%alloc_cell_blocks))
      DO i=1,smooth_initial_velocity_iterations
        WRITE(message_text,*) "Smoothing velocities..., miss value=", missValue
        CALL message(method_name, message_text)

        suv = u
        CALL smooth_onCells(patch_3D=patch_3d, &
          & in_value=suv, out_value=u,  &
          & smooth_weights=smooth_initial_velocity_weights,      &
          & has_missValue=has_missValue, missValue=missValue)
        suv = v
        CALL smooth_onCells(patch_3D=patch_3d, &
          & in_value=suv, out_value=v,  &
          & smooth_weights=smooth_initial_velocity_weights,      &
          & has_missValue=has_missValue, missValue=missValue)
      ENDDO
      DEALLOCATE(suv)
    
    ENDIF
    
    CALL message(method_name, "fillVerticallyMissingValues...")
    CALL fillVerticallyMissingValues(patch_3d=patch_3d, ocean_tracer=u, &
      & has_missValue=has_missValue, missValue=missValue)
    CALL fillVerticallyMissingValues(patch_3d=patch_3d, ocean_tracer=v, &
      & has_missValue=has_missValue, missValue=missValue)
 
    CALL message(method_name, "gvec2cvec ...")
    ALLOCATE(cellVelocity_cc(nproma,n_zlev, patch_2d%alloc_cell_blocks))

    DO cell_block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, cell_block, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index
        lon = patch_2d%cells%center(cell_index,cell_block)%lon
        lat = patch_2d%cells%center(cell_index,cell_block)%lat
 
        DO level = 1, n_zlev
          CALL gvec2cvec (u(cell_index, level, cell_block), &
                  &  v(cell_index, level, cell_block), &
                  &  lon, lat,                         &
                  & cellVelocity_cc(cell_index, level, cell_block)%x(1), &
                  & cellVelocity_cc(cell_index, level, cell_block)%x(2), &
                  & cellVelocity_cc(cell_index, level, cell_block)%x(3), &
                  & patch_2D%geometry_info)
    
        ENDDO
      ENDDO
    ENDDO
        
    CALL message(method_name, "map_cell2edges_3D ...")
    CALL map_cell2edges_3D( patch_3D, cellVelocity_cc, normal_velocity, this_operators_coeff)
    CALL sync_patch_array(sync_e, patch_2D, normal_velocity)

    CALL message(method_name, "DEALLOCATE ...")
    DEALLOCATE(cellVelocity_cc)
    DEALLOCATE(u, v)
 
  END SUBROUTINE init_ocean_velocity_uv_fromFile
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_velocity_vn_fromFile(patch_3d, normal_velocity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: normal_velocity(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    LOGICAL  :: has_missValue
    REAL(wp) :: missValue

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_velocity_vn_fromFile'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    CALL message(method_name, ': init from file')
    CALL init_edge_3D_variable_fromFile(patch_3d, variable=normal_velocity, name="vn", &
      & has_missValue=has_missValue, missValue=missValue)
    
!     CALL sync_patch_array(sync_e, patch_2D, normal_velocity) done in the file read

  END SUBROUTINE init_ocean_velocity_vn_fromFile
  !-------------------------------------------------------------------------------

        
  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_velocity(patch_3d, normal_velocity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: normal_velocity(:,:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_velocity'
    !-------------------------------------------------------------------------

    normal_velocity(:,:,:) = 0.0_wp

    ! IF (initial_velocity_type < 200) RETURN ! not analytic velocity

    SELECT CASE (initial_velocity_type)

    CASE (000)
      RETURN
 
    CASE (001)
      CALL init_ocean_velocity_uv_fromFile(patch_3d, normal_velocity)

    CASE (002)
      CALL init_ocean_velocity_vn_fromFile(patch_3d, normal_velocity)
        
        
    !------------------------------
    CASE (200)
      ! uniform velocity
      normal_velocity(:,:,:) = 0.0_wp
      CALL message(method_name, ': uniform zero velocity')

    !------------------------------
    CASE (201)
      CALL velocity_LauterRotation(patch_3d, normal_velocity)

    !------------------------------
    CASE (202)
      CALL message(method_name, 'Williamson Test 2 ')
      CALL velocity_WilliamsonTest_2_5(patch_3d, normal_velocity, velocity_amplitude=u0)

    !------------------------------
    CASE (203)
      CALL message(method_name, 'Williamson Test 5 ')
      CALL velocity_WilliamsonTest_2_5(patch_3d, normal_velocity, velocity_amplitude=initial_velocity_amplitude)

    CASE (206)
      CALL message(method_name, 'Williamson Test 6 ')
      CALL velocity_WilliamsonTest_2_6(patch_3d, normal_velocity, velocity_amplitude=initial_velocity_amplitude)
  
    CASE (207)
      CALL message(method_name, 'Galewsky Test ')
      CALL velocity_GalewskyTest(patch_3d, normal_velocity)
 
    CASE (220)
      CALL message(TRIM(method_name), 'Uniform horz velocity')
      CALL velocity_uniHorz(patch_3d, normal_velocity)
 

     CASE (300)
      CALL message(method_name, 'Velocity Kelvin-Helmholtz Test ')
      CALL velocity_KelvinHelmholtzTest(patch_3d, normal_velocity, velocity_amplitude=initial_velocity_amplitude)

 
    !------------------------------
     CASE default
      CALL finish(method_name, "unknown initial_velocity_type")

    END SELECT

    CALL dbg_print('init_ocean_velocity', normal_velocity(:,:,:), &
      & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%edges%owned)

  END SUBROUTINE init_ocean_velocity
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE init_ocean_surface_height(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    REAL(wp), ALLOCATABLE :: old_height(:,:)
    LOGICAL  :: has_missValue
    REAL(wp) :: missValue

    INTEGER :: i
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_ocean_surface_height'
    !-------------------------------------------------------------------------

    ocean_height(:,:) = 0.0_wp
    has_missValue = .false.
    missValue     = -99999999.0_wp

    patch_2d => patch_3d%p_patch_2d(1)

    ! needs to be written with calls !
    SELECT CASE (sea_surface_height_type)
    !------------------------------
    CASE (000)
      CALL message(method_name, ' no initialization')

    CASE (001)
      CALL message(method_name, ': init from file')
      CALL init_cell_2D_variable_fromFile(patch_3d, variable=ocean_height, name="h", &
         & has_missValue=has_missValue, missValue=missValue)
     
    CASE (200)
      ! 0 height, this is the initialization value,
      ! so no need to explicilty define this case
      ! the whole grid is considered sea
      ocean_height(:,:) = 0.0_wp

    CASE (201)
      CALL height_sinLon_cosLat(patch_3d, ocean_height)

    CASE (202)
      CALL height_exponentialDistance(patch_3d, ocean_height)

    CASE (203)
      CALL height_LauterRotation(patch_3d, ocean_height)

    CASE (204)
      CALL height_WilliamsonTest2(patch_3d, ocean_height)

    CASE (205)
      CALL height_WilliamsonTest5(patch_3d, ocean_height)
      
    CASE (206)
      CALL height_WilliamsonTest6(patch_3d, ocean_height)

    CASE (207)
      CALL height_GalewskyTest(patch_3d, ocean_height)

    CASE (221)
      CALL height_quads_checkerboard(patch_3d=patch_3d, ocean_height=ocean_height, base_value=1.0_wp, variation=1.0_wp)

    CASE (235)
      ocean_height(:,:) = 10.0_wp

    CASE default
      CALL finish(method_name, "unknown sea_surface_height_type")

    END SELECT

    IF (smooth_initial_height_iterations > 0) THEN
      ALLOCATE(old_height(nproma,patch_2d%alloc_cell_blocks))
      DO i=1,smooth_initial_height_iterations
        CALL message(method_name, "Smoothing initial height...")
        old_height = ocean_height
        CALL smooth_onCells(patch_3D=patch_3d, &
          & in_value=old_height, out_value=ocean_height,  &
          & smooth_weights=smooth_initial_height_weights, &
          & has_missValue=has_missValue, missValue=missValue)
      ENDDO
      DEALLOCATE(old_height)
    ENDIF
    
    CALL dbg_print('init_ocean_surface_height', ocean_height, module_name,  1, &
        & in_subset=patch_2d%cells%owned)
  END SUBROUTINE init_ocean_surface_height
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE height_sinLon_cosLat(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':height_sinLon_cosLat'
    !-------------------------------------------------------------------------
    ! CASE (201)
    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! #slo#: simple elevation between 30W and 30E (pi/3.)
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        IF ( patch_3d%lsm_c(idx,1,block) <= sea_boundary ) THEN

          ocean_height(idx,block) = 10.0_wp * &
            & SIN(cell_center(idx, block)%lon * 6.0_wp) * COS(cell_center(idx, block)%lat * 3.0_wp)

        ENDIF
      END DO
    END DO

  END SUBROUTINE height_sinLon_cosLat
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE height_exponentialDistance(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: distan
    REAL(wp):: perturbation_lat, perturbation_lon

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':height_exponentialDistance'
    !-------------------------------------------------------------------------
    ! CASE (202)
    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Add elevation perturbation at new values - 35N; 10W
    ! not clear yet
    perturbation_lat = basin_center_lat + 0.1_wp * basin_height_deg
    perturbation_lon = basin_center_lon + 0.1_wp * basin_width_deg
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        IF (patch_3d%p_patch_1d(1)%dolic_c(idx,block) > 0) THEN

          distan=SQRT((cell_center(idx, block)%lat - perturbation_lat * deg2rad)**2 + &
            & (cell_center(idx, block)%lon - perturbation_lon * deg2rad)**2)
          !IF(distan<=15.5_wp*deg2rad) cycle
          IF(distan < 10.0_wp * deg2rad) THEN
            ocean_height(idx,block) = 0.5_wp & !ocean_state%p_prog(nold(1))%h(idx,block)&
              & + 0.3_wp * EXP(-(distan/(2.2_wp*deg2rad))**2)
          ENDIF

        ENDIF

      END DO
    END DO

  END SUBROUTINE height_exponentialDistance
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! Initial datum for height h, test case unsteady solid body
  ! rotation of L\"auter et al.(2007).
  SUBROUTINE height_LauterRotation(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':height_LauterRotation'
    !-------------------------------------------------------------------------
    ! CASE (203)

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! test_usbr_h
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        ocean_height(idx,block)  = test_usbr_h( cell_center(idx, block)%lon, cell_center(idx, block)%lat, 0.0_wp)
      END DO
    END DO

  END SUBROUTINE height_LauterRotation
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE height_WilliamsonTest2(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index

!   INTEGER :: v1_idx, v1_blk, v2_idx, v2_blk, v3_idx, v3_blk
!   TYPE(t_cartesian_coordinates), POINTER :: vertex_cartesian(:,:)
!   TYPE(t_geographical_coordinates), POINTER :: vertex_lonlat(:,:)
!   TYPE(t_cartesian_coordinates) :: barycenter_cartesian
!   TYPE(t_geographical_coordinates) :: barycenter_lonlat

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':height_WilliamsonTest2'
    !-------------------------------------------------------------------------
    ! CASE (204)

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center
!   vertex_cartesian => patch_2d%verts%cartesian
!   vertex_lonlat => patch_2d%verts%vertex

    ! test2_h
    CALL message(method_name, ' h for Williamson Test 2')
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
  !         v1_idx = patch_2d%cells%vertex_idx(idx, block, 1)
  !         v1_blk = patch_2d%cells%vertex_blk(idx, block, 1)
  !         v2_idx = patch_2d%cells%vertex_idx(idx, block, 2)
  !         v2_blk = patch_2d%cells%vertex_blk(idx, block, 2)
  !         v3_idx = patch_2d%cells%vertex_idx(idx, block, 3)
  !         v3_blk = patch_2d%cells%vertex_blk(idx, block, 3)
  !         min_lat = MIN( vertex_lonlat(v1_idx, v1_blk)%lat, vertex_lonlat(v2_idx, v2_blk)%lat, &
  !           & vertex_lonlat(v3_idx, v3_blk)%lat)
  !         max_lat = MAX( vertex_lonlat(v1_idx, v1_blk)%lat, vertex_lonlat(v2_idx, v2_blk)%lat, &
  !           & vertex_lonlat(v3_idx, v3_blk)%lat)
  !         lat = (min_lat+max_lat) * 0.5_wp
  !         ocean_height(idx,block) = test2_h( cell_center(idx, block)%lon, lat, 0.0_wp)
  !         barycenter_cartesian%x = (vertex_cartesian(v1_idx, v1_blk)%x + &
  !                                   vertex_cartesian(v2_idx, v2_blk)%x + &
  !                                   vertex_cartesian(v3_idx, v3_blk)%x) / 3.0_wp
  !         barycenter_lonlat = cc2gc(barycenter_cartesian)
  !         ocean_height(idx,block) = test2_h( barycenter_lonlat%lon, barycenter_lonlat%lat, 0.0_wp)

          ! this is the correct one
        ocean_height(idx,block) = test2_h( cell_center(idx, block)%lon, cell_center(idx, block)%lat)
      END DO
    END DO

  END SUBROUTINE height_WilliamsonTest2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE height_WilliamsonTest5(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':height_WilliamsonTest5'
    !-------------------------------------------------------------------------
    ! CASE (205)
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! test5_h
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        ocean_height(idx,block) = test5_h( cell_center(idx, block)%lat)
      END DO
    END DO

  END SUBROUTINE height_WilliamsonTest5
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE height_WilliamsonTest6(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':height_WilliamsonTest6'
    !-------------------------------------------------------------------------
    ! CASE (205)
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! test6_h
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        ocean_height(idx,block) = Williamson_test6_h( cell_center(idx, block)%lon, cell_center(idx, block)%lat, 0.0_wp)
      END DO
    END DO
write(0,*)'Williamson-Test6:h', maxval(ocean_height),minval(ocean_height)
  END SUBROUTINE height_WilliamsonTest6
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE height_GalewskyTest(patch_3d, ocean_height)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: h_perturb, phi_2, alpha,beta

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':height_GalewskyTest'
    !-------------------------------------------------------------------------
    ! CASE (205)
    patch_2d    => patch_3d%p_patch_2d(1)
    all_cells   => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Basic height
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        ocean_height(idx,block) = galewsky_h(cell_center(idx, block)%lat)
      END DO
    END DO
 write(0,*)'Galewsky-Test:h', maxval(ocean_height),minval(ocean_height)   
    !Add perturbation
    phi_2=0.25_wp*pi
    beta=1.0_wp/15.0_wp
    alpha=1.0_wp/3.0_wp
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        IF(cell_center(idx,block)%lon<= pi .and. cell_center(idx,block)%lon> -pi)THEN
        h_perturb=120.0_wp*cos(cell_center(idx, block)%lat)&
        &*exp(-(cell_center(idx,block)%lon/alpha)**2)*exp(-((phi_2-cell_center(idx, block)%lat)/beta)**2)
        ocean_height(idx,block) = ocean_height(idx,block)+h_perturb
! write(123,*)'perturb',ocean_height(idx,block)-h_perturb,ocean_height(idx,block), h_perturb,&
! &cos(cell_center(idx, block)%lat),&
! &exp(-(cell_center(idx,block)%lon/alpha)**2),exp(-((phi_2-cell_center(idx, block)%lat)/beta)**2)
        ENDIF
      
      END DO
    END DO
    
    !
write(0,*)'Galewsky-Test:h', maxval(ocean_height),minval(ocean_height)
  END SUBROUTINE height_GalewskyTest
  !-------------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------------
  !> Initial datum for zonal velocity u, test case unsteady solid body
  ! rotation of L\"auter et al.(2007).
  !-------------------------------------------------------------------------
  SUBROUTINE velocity_LauterRotation(patch_3d, vn)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: vn(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges

    INTEGER :: edge_block, edge_index, level
    INTEGER :: start_edges_index, end_edges_index
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: t       ! point of time
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: angle1, angle2, edge_vn, COS_angle1, SIN_angle1

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':velocity_usbr_u'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL

    t = 0.0_wp
    angle1 = .25_wp * pi
    COS_angle1 = COS(angle1)
    SIN_angle1 = SIN(angle1)

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat

        angle2 = point_lon + grid_angular_velocity * t
        uu = COS(point_lat) * COS_angle1
        uu = uu + COS(angle2) * SIN(point_lat) * SIN_angle1
        uu = u0 * uu

        vv = SIN(angle2) * SIN_angle1
        vv = -1._wp * u0 * vv

        edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
              & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          vn(edge_index, level, edge_block) = edge_vn
        ENDDO

      ENDDO
    ENDDO

  END SUBROUTINE  velocity_LauterRotation
  !-----------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------
  SUBROUTINE velocity_WilliamsonTest_2_5(patch_3d, vn, velocity_amplitude)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: vn(:,:,:)
    REAL(wp), INTENT(in) :: velocity_amplitude

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_cartesian_coordinates), POINTER ::cellVelocity_cc(:,:,:)

    INTEGER :: level
    INTEGER :: cell_block, cell_index
    INTEGER :: start_cells_index, end_cells_index
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: x1, x2, x3
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':velocity_WilliamsonTest_2_5'
    !-------------------------------------------------------------------------

    ! CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    ! fisrt calculate th velocyt at cell centers
    ALLOCATE(cellVelocity_cc(nproma,n_zlev, patch_2d%alloc_cell_blocks))
    DO cell_block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, cell_block, start_cells_index, end_cells_index)
      DO cell_index = start_cells_index, end_cells_index
        point_lon = patch_2d%cells%center(cell_index,cell_block)%lon
        point_lat = patch_2d%cells%center(cell_index,cell_block)%lat

        uu = COS(point_lat) * COS(aleph)
        uu = uu + COS(point_lon) * SIN(point_lat) * SIN(aleph)
        uu = velocity_amplitude * uu

        vv = SIN(point_lon) * SIN(aleph)
        vv = -1._wp * velocity_amplitude * vv

        CALL gvec2cvec(  uu, vv, point_lon, point_lat, x1, x2, x3)
        DO level = 1, n_zlev
          cellVelocity_cc(cell_index, level, cell_block)%x(1) = x1 
          cellVelocity_cc(cell_index, level, cell_block)%x(2) = x2
          cellVelocity_cc(cell_index, level, cell_block)%x(3) = x3
        ENDDO

      ENDDO
    ENDDO

    ! map velocity to edge centers
    CALL map_cell2edges_3D( patch_3D, cellVelocity_cc, vn, this_operators_coeff)
    CALL sync_patch_array(sync_e, patch_2D, vn)
    
    DEALLOCATE(cellVelocity_cc)


    ! calculate the velocity directly on edges
!     DO edge_block = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
!       DO edge_index = start_edges_index, end_edges_index
!         point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
!         point_lat = patch_2d%edges%center(edge_index,edge_block)%lat
! 
!         uu = COS(point_lat) * COS(aleph)
!         uu = uu + COS(point_lon) * SIN(point_lat) * SIN(aleph)
!         uu = velocity_amplitude * uu
! 
!         vv = SIN(point_lon) * SIN(aleph)
!         vv = -1._wp * velocity_amplitude * vv
! 
!         edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
!               & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2
! 
!         DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
!           vn(edge_index, level, edge_block) = edge_vn
!         ENDDO
! 
!       ENDDO
!     ENDDO

  END SUBROUTINE  velocity_WilliamsonTest_2_5
  !-----------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------
  SUBROUTINE velocity_WilliamsonTest_2_6(patch_3d, vn, velocity_amplitude)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: vn(:,:,:)
    REAL(wp), INTENT(in) :: velocity_amplitude

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges

    INTEGER :: edge_block, edge_index, level
    INTEGER :: start_edges_index, end_edges_index
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: edge_vn

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':velocity_WilliamsonTest_2_6'
    !-------------------------------------------------------------------------

    ! CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat

        uu = Williamson_test6_u(point_lon,point_lat,velocity_amplitude)

        vv = Williamson_test6_v(point_lon,point_lat,velocity_amplitude)

        edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
              & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          vn(edge_index, level, edge_block) = edge_vn
        ENDDO

      ENDDO
    ENDDO
write(0,*)'Williamson-Test6:vn', maxval(vn),minval(vn)
  END SUBROUTINE  velocity_WilliamsonTest_2_6
  !-----------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------
  SUBROUTINE velocity_GalewskyTest(patch_3d, vn)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: vn(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges

    INTEGER :: edge_block, edge_index, level
    INTEGER :: start_edges_index, end_edges_index
    REAL(wp) :: point_lat     ! latitude of point
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: edge_vn

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':velocity_WilliamsonTest_2_5'
    !-------------------------------------------------------------------------

    ! CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat

        uu = galewsky_u(point_lat)

        vv = galewsky_v()

        edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
              & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          vn(edge_index, level, edge_block) = edge_vn
        ENDDO

      ENDDO
    ENDDO
 !write(0,*)'Galewsky-Test6:vn', maxval(vn),minval(vn)
  END SUBROUTINE  velocity_GalewskyTest
  !-----------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !> Uniform velocity
  !-------------------------------------------------------------------------
  SUBROUTINE velocity_uniHorz(patch_3d, vn)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: vn(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges

    INTEGER :: edge_block, edge_index, level
    INTEGER :: start_edges_index, end_edges_index
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: t       ! point of time
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: angle1, angle2, edge_vn, COS_angle1, SIN_angle1

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':velocity_usbr_u'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
        uu = 10.0_wp

        vv =  0._wp

        edge_vn = uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
              & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          vn(edge_index, level, edge_block) = edge_vn
        ENDDO

      ENDDO
    ENDDO

  END SUBROUTINE  velocity_uniHorz
  !-----------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------
  SUBROUTINE velocity_KelvinHelmholtzTest(patch_3d, vn, velocity_amplitude)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: vn(:,:,:)
    REAL(wp), INTENT(in) :: velocity_amplitude

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges

    INTEGER :: edge_block, edge_index
    INTEGER :: start_edges_index, end_edges_index
    REAL(wp) :: point_lon, point_lat     ! latitude of point
    REAL(wp) :: uu, vv      ! zonal,  meridional velocity
    REAL(wp) :: edge_vn!, COS_angle1, SIN_angle1
    REAL(wp) :: shear_depth
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':velocity_KelvinHelmholtz'
    !-------------------------------------------------------------------------

    ! CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_edges => patch_2d%edges%ALL
    
    shear_depth  = 0.05_wp
    
    edge_vn = 0.1_wp

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index
      
         point_lon = patch_2d%edges%center(edge_index,edge_block)%lon* rad2deg
         point_lat = patch_2d%edges%center(edge_index,edge_block)%lat* rad2deg
!          IF(point_lat>=basin_center_lat)THEN   
!            !uu=tanh((point_lat-0.025)*300.0_wp) 
!            uu=tanh((point_lat-shear_depth)*300.0_wp) 
!          ELSEIF(point_lat<basin_center_lat)THEN
!            !uu=tanh((0.75_wp-point_lat)*300.0_wp) 
!            uu=tanh((shear_depth-point_lat)*300.0_wp) 
!          ENDIF
          IF(point_lat>=basin_center_lat)THEN
            !uu=tanh((point_lat-0.025)*300.0_wp) 
            uu=tanh((point_lat+shear_depth)*300.0_wp) 
          ELSEIF(point_lat<basin_center_lat)THEN
            !uu=tanh((0.75_wp-point_lat)*300.0_wp) 
            uu=tanh((shear_depth-point_lat)*300.0_wp) 
          ENDIF
        vv=0.1_wp*sin(2.0_wp*pi*point_lon)
   
        edge_vn =(uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1 &
           & + vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2) 
           vn(edge_index, 1:n_zlev, edge_block) = velocity_amplitude*edge_vn

      ENDDO
    ENDDO

  END SUBROUTINE  velocity_KelvinHelmholtzTest
  !-----------------------------------------------------------------------------------
 
  
  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_uniform_SeparationAtLon(patch_3d, ocean_temperature, wallLonDeg)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)
    REAL(wp), INTENT(in) :: wallLonDeg

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    
    REAL(wp) :: wallLonRads, basin_center_lon_rads

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_uniform_SeparationAtLon'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    wallLonRads = wallLonDeg * deg2rad
    basin_center_lon_rads = basin_center_lon * deg2rad

    !Add horizontal variation
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        IF((patch_2d%cells%center(idx,block)%lon - basin_center_lon_rads) >= wallLonRads) THEN
          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
            ocean_temperature(idx,level,block) = 10.0_wp
          ENDDO
        ELSE
          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
            ocean_temperature(idx,level,block) = 5.0_wp
          ENDDO
        ENDIF

      END DO
    END DO

  END SUBROUTINE temperature_uniform_SeparationAtLon
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_uniform_SeparationAtLat(patch_3d, ocean_temperature, wallLatDeg)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)
    REAL(wp), INTENT(in) :: wallLatDeg

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: wallLatRads, basin_center_lat_rads

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_uniform_SeparationAtLon'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    wallLatRads = wallLatDeg * deg2rad
    basin_center_lat_rads = basin_center_lat * deg2rad

    !Add horizontal variation
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        IF((patch_2d%cells%center(idx,block)%lat - basin_center_lat_rads) >= wallLatRads)THEN
          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
            ocean_temperature(idx,level,block) = 10.0_wp
          ENDDO
        ELSE
          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
            ocean_temperature(idx,level,block) = 5.0_wp
          ENDDO
        ENDIF

      END DO
    END DO

  END SUBROUTINE temperature_uniform_SeparationAtLat
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_AddLocalPerturbation(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_AddLocalPerturbation'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !Add horizontal variation
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        lat_deg = patch_2d%cells%center(idx,block)%lat*rad2deg
        lon_deg = patch_2d%cells%center(idx,block)%lon*rad2deg

        IF(ABS(lon_deg) < 2.5_wp .AND. &
           ABS(lat_deg-basin_center_lat) < 0.25_wp*basin_height_deg) THEN

          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
            ocean_temperature(idx,level,block) = ocean_temperature(idx,level,block) * 0.75_wp
          END DO

        ENDIF

      END DO
    END DO
  END SUBROUTINE temperature_AddLocalPerturbation
  !-------------------------------------------------------------------------------


!   !-------------------------------------------------------------------------------
!   SUBROUTINE temperature_AddSinusoidalPerturbation(patch_3d, ocean_temperature)
!     TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
!     REAL(wp), TARGET :: ocean_temperature(:,:,:)
! 
!     TYPE(t_patch),POINTER   :: patch_2d
!     TYPE(t_subset_range), POINTER :: all_cells
! 
!     INTEGER :: block, idx, level
!     INTEGER :: start_cell_index, end_cell_index
!     REAL(wp):: lat_deg, lon_deg
! 
!     CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_AddLocalPerturbation'
!     !-------------------------------------------------------------------------
! 
!     CALL message(method_name, ' ')
! 
!     patch_2d => patch_3d%p_patch_2d(1)
!     all_cells => patch_2d%cells%ALL
! 
!     !Add horizontal variation
!     DO block = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!       DO idx = start_cell_index, end_cell_index
!         lat_deg = patch_2d%cells%center(idx,block)%lat*rad2deg
!         lon_deg = patch_2d%cells%center(idx,block)%lon*rad2deg
! 
! !        IF(ABS(lon_deg) < 2.5_wp .AND. &
! !          ABS(lat_deg-basin_center_lat) < 0.25_wp*basin_height_deg) THEN
! !        write(123,*)'t-perturb',ocean_temperature(idx,1,block),ocean_temperature(idx,2,block)*&
! !            &0.01_wp*sin(50.0_wp*patch_2d%cells%center(idx,block)%lon)
!           IF(  lat_deg<= basin_center_lat+ 0.5_wp*basin_height_deg)THEN    
!           DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
!             ocean_temperature(idx,level,block) = ocean_temperature(idx,level,block) +ocean_temperature(idx,level,block)*&
!             &0.01_wp*sin(50.0_wp*patch_2d%cells%center(idx,block)%lon)
!           END DO
! 
!         ENDIF
! 
!       END DO
!     END DO
!       
!   END SUBROUTINE temperature_AddSinusoidalPerturbation
!   !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_AddSinusoidalPerturbation(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_AddLocalPerturbation'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !Add horizontal variation
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        lat_deg = patch_2d%cells%center(idx,block)%lat*rad2deg

!        IF(ABS(lon_deg) < 2.5_wp .AND. &
!          ABS(lat_deg-basin_center_lat) < 0.25_wp*basin_height_deg) THEN
!        write(123,*)'t-perturb',ocean_temperature(idx,1,block),ocean_temperature(idx,2,block)*&
!            &0.01_wp*sin(50.0_wp*patch_2d%cells%center(idx,block)%lon)
          IF(  lat_deg<= basin_center_lat+ 0.5_wp*basin_height_deg)THEN    
          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
            ocean_temperature(idx,level,block) = ocean_temperature(idx,level,block) +ocean_temperature(idx,level,block)*&
            &0.01_wp*sin(50.0_wp*patch_2d%cells%center(idx,block)%lon)
          END DO

        ENDIF

      END DO
    END DO
  END SUBROUTINE temperature_AddSinusoidalPerturbation
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_AddHorizontalVariation(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, z_temp_max

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_AddHorizontalVariation'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !Add horizontal variation
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        lat_deg = patch_2d%cells%center(idx,block)%lat
        z_temp_max=0.01_wp* (lat_deg-basin_center_lat) * (lat_deg-basin_center_lat)

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

          ocean_temperature(idx,level,block) = ocean_temperature(idx,level,block) * &
            & EXP(-z_temp_max/basin_height_deg)!(1.0_wp-exp(-z_temp_max/basin_height_deg))

        END DO
      END DO
    END DO

  END SUBROUTINE temperature_AddHorizontalVariation
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_Uniform_SpecialArea(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg
    REAL(wp):: z_lat1, z_lat2, z_lon1, z_lon2

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_Uniform_SpecialArea'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! 2012-10-31: Indonesian Archipelago - connected to Atlantic? (4 cpu)
    z_lat1  =  -5.0_wp
    z_lat2  =  10.0_wp
    z_lon1  = 115.0_wp
    z_lon2  = 135.0_wp
    ! 2012-10-31: corresponding NAtl: 1N 20W, 3 levels
    z_lon1  = 145.0_wp
    z_lon2  = 160.0_wp
    ! 2012-11-08: Test homogen with one cell differ at 10N; 70E
    z_lat1  =   0.0_wp
    z_lat2  =  15.0_wp
    z_lon1  =  60.0_wp
    z_lon2  =  80.0_wp

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = cell_center(idx, block)%lat * rad2deg
        lon_deg = cell_center(idx, block)%lon * rad2deg

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

          ocean_temperature(idx,level,block) = 5.0_wp

          IF ( (lat_deg >= z_lat1 .AND. lat_deg <= z_lat2) .AND. &
            & (lon_deg >= z_lon1 .AND. lon_deg <= z_lon2) .AND. &
            & ( level <= 1 ) ) THEN

            ocean_temperature(idx,level,block) =  6.0_wp

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE temperature_Uniform_SpecialArea
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_front(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, width

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_front'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    width = 0.0_wp
    ocean_temperature(:,:,:)=0.0_wp

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
  
        lat_deg = patch_2d%cells%center(idx,block)%lat*rad2deg

         DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
           IF( lat_deg>= basin_center_lat+ width)THEN    

            ocean_temperature(idx,level,block) =  initial_temperature_north
          !lower channel boudary
          ELSEIF( lat_deg<= basin_center_lat- width )THEN 
            ocean_temperature(idx,level,block) =  initial_temperature_south  
          ELSE!channel interior  
          ENDIF
        END DO
      END DO
    END DO

  END SUBROUTINE temperature_front
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE salinity_Uniform_SpecialArea(patch_3d, ocean_salinity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg
    REAL(wp):: z_lat1, z_lat2, z_lon1, z_lon2

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':salinity_Uniform_SpecialArea'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! 2012-10-31: Indonesian Archipelago - connected to Atlantic? (4 cpu)
    z_lat1  =  -5.0_wp
    z_lat2  =  10.0_wp
    z_lon1  = 115.0_wp
    z_lon2  = 135.0_wp
    ! 2012-10-31: corresponding NAtl: 1N 20W, 3 levels
    z_lon1  = 145.0_wp
    z_lon2  = 160.0_wp
    ! 2012-11-08: Test homogen with one cell differ at 10N; 70E
    z_lat1  =   0.0_wp
    z_lat2  =  15.0_wp
    z_lon1  =  60.0_wp
    z_lon2  =  80.0_wp

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = cell_center(idx, block)%lat * rad2deg
        lon_deg = cell_center(idx, block)%lon * rad2deg

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

          ocean_salinity(idx,level,block) = 35.0_wp

          IF ( (lat_deg >= z_lat1 .AND. lat_deg <= z_lat2) .AND. &
             & (lon_deg >= z_lon1 .AND. lon_deg <= z_lon2) .AND. &
             & ( level <= 1 ) ) THEN

             ocean_salinity(idx,level,block) = 34.8_wp

          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE salinity_Uniform_SpecialArea
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_APE(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_APE'
    !-------------------------------------------------------------------------
    ! Important:
    !   use initial_temperature_top=27.0 initial_temperature_bottom=0.0
    !   to be consistent with the old setup
    CALL message(method_name, ' using sst1')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
         DO level=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(idx,block))
          ocean_temperature(idx,level,block) = MIN(MAX( &
            & ape_sst(initial_sst_type, patch_2d%cells%center(idx,block)%lat) - tmelt,  & ! SST in Celsius
            & initial_temperature_bottom), initial_temperature_top)
        END DO
      END DO
    END DO

    CALL increaseTracerVerticallyLinearly(patch_3d, ocean_tracer=ocean_temperature, &
      & bottom_value=initial_temperature_bottom)

  END SUBROUTINE temperature_APE
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_smoothAPE_LinearLevels(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: temperature_difference, poleLat, waveNo

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_APE'
    !-------------------------------------------------------------------------
    CALL message(method_name, ' using smoothAPE')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    temperature_difference = initial_temperature_top - initial_temperature_bottom
    poleLat = ABS(forcing_temperature_poleLat * deg2rad)
    waveNo = pi_2 / poleLat
    
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        DO level=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(idx,block))
          ocean_temperature(idx,level,block) = MAX(initial_temperature_bottom + &
            & (COS(waveNo * MIN(ABS(patch_2d%cells%center(idx,block)%lat), poleLat))**2) * temperature_difference, &
            & initial_temperature_bottom)
        END DO
      END DO
    END DO

    ! add meridional temperature slope over all latitudes
    
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        DO level=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(idx,block))
          ocean_temperature(idx,level,block) = MAX(ocean_temperature(idx,level,block) + &
            & (patch_2d%cells%center(idx,block)%lat + poleLat) / pi_2 * initial_temperature_shift, &
            & initial_temperature_bottom)
        END DO
      END DO
    END DO

    ! decrease of temperature from top to bottom
    CALL increaseTracerLevelsLinearly(patch_3d, ocean_tracer=ocean_temperature, &
      & bottom_value=initial_temperature_bottom)

  END SUBROUTINE temperature_smoothAPE_LinearLevels
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE tracer_smoothAPE_LinearDepth(patch_3d, ocean_tracer, top_value, bottom_value)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in) :: top_value, bottom_value

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: tracer_difference, poleLat, waveNo

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_APE'
    !-------------------------------------------------------------------------
    CALL message(method_name, ' using smoothAPE')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    tracer_difference = top_value - bottom_value
    ! #slo#: Caution, there is a mixture of forcing and initialization!
    poleLat = ABS(forcing_temperature_poleLat * deg2rad)
    waveNo = pi_2 / poleLat

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        DO level=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(idx,block))
          ocean_tracer(idx,level,block) = bottom_value + &
            & (COS(waveNo * MIN(ABS(patch_2d%cells%center(idx,block)%lat), poleLat))**2) * tracer_difference
        END DO
      END DO
    END DO

    ! add meridional tracer slope over all latitudes
    IF (initial_temperature_shift /= 0.0_wp) THEN
      DO block = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
        DO idx = start_cell_index, end_cell_index
          DO level=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(idx,block))
            ocean_tracer(idx,level,block) = ocean_tracer(idx,level,block) + &
              & (patch_2d%cells%center(idx,block)%lat + poleLat) / pi_2 * initial_temperature_shift
          END DO
        END DO
      END DO
    ENDIF

    ! decrease of tracer from top to bottom
    CALL increaseTracerVerticallyLinearly(patch_3d, ocean_tracer=ocean_tracer, &
      & bottom_value=bottom_value)

  END SUBROUTINE tracer_smoothAPE_LinearDepth
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE SST_LinearMeridional(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: temperature_difference, basin_northBoundary, basin_southBoundary, lat_diff
    REAL(wp) :: lat(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':SST_LinearMeridional'
    !-------------------------------------------------------------------------
    CALL message(method_name, ' using meridional gradient over basin height')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

 !  CALL assign_if_present(length,length_opt)

    lat(:,:) = patch_2d%cells%center(:,:)%lat

    ocean_temperature(:,:,:) = initial_temperature_south

    temperature_difference = initial_temperature_north - initial_temperature_south
    basin_northBoundary    = (basin_center_lat + 0.5_wp*basin_height_deg) * deg2rad
    basin_southBoundary    = (basin_center_lat - 0.5_wp*basin_height_deg) * deg2rad
    lat_diff               = basin_northBoundary - basin_southBoundary  !  basin_height_deg*deg2rad
    
!     write(0,*) "*****************************"
!     write(0,*) "init SST_LinearMeridional"
!     write(0,*) "T:", initial_temperature_north, initial_temperature_south, temperature_difference
    
    
    level=1
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        ocean_temperature(idx,level,block) = &
          & initial_temperature_north - temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
        ocean_temperature(idx,level,block) = &
          & MERGE(ocean_temperature(idx,level,block), initial_temperature_north, lat(idx,block)<basin_northBoundary)
        ocean_temperature(idx,level,block) = &
          & MERGE(ocean_temperature(idx,level,block), initial_temperature_south, lat(idx,block)>basin_southBoundary)
          
!         write(0,*) "T=", ocean_temperature(idx,level,block), " at:", lat(idx,block)*rad2deg
      END DO
    END DO

  END SUBROUTINE SST_LinearMeridional
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE SST_Abernathey_01(patch_3d, ocean_temperature, VariationAmplitude, VariationLength, VariationWaveNo, SouthLat)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)
    REAL(wp), INTENT(in) :: VariationAmplitude, VariationLength, VariationWaveNo
    REAL(wp), INTENT(in) :: SouthLat
    
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: y_lat
              
    CALL SST_LinearMeridional(patch_3d, ocean_temperature)

    ! add a perturbation analogous to the forcing
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%all

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        IF (patch_3D%lsm_c(idx,1,block) <= sea_boundary) THEN

          y_lat = patch_2d%cells%center(idx,block)%lat - SouthLat

          ocean_temperature(idx,1,block) = &
            & ocean_temperature(idx,1,block) - VariationAmplitude * COS(VariationWaveNo * pi * y_lat/VariationLength)

        ENDIF
      END DO
    END DO
    
!     CALL SST_constant(patch_3d=patch_3d, ocean_temperature=ocean_temperature, &
!       & constant_temperature=NorthTemperature, &
!       & LowerLat=NorthLat)
    
          
  END SUBROUTINE SST_Abernathey_01
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE SST_constant(patch_3d, ocean_temperature, constant_temperature, LowerLat)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)
    REAL(wp), INTENT(in) :: constant_temperature
    REAL(wp), OPTIONAL :: LowerLat
    
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: lower_lat

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':SST_constant'
    !-------------------------------------------------------------------------
    CALL message(method_name, ' using meridional gradient over basin height')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    IF (PRESENT(LowerLat)) THEN
      lower_lat = LowerLat
    ELSE
      lower_lat = -100.0
    ENDIF
    
    level=1
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        IF (patch_2d%cells%center(idx,block)%lat >= lower_lat) THEN
          ocean_temperature(idx,level,block) = constant_temperature
        ENDIF
      END DO
    END DO

  END SUBROUTINE SST_constant
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_dirac_signal(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg
    REAL(wp):: z_lat1, z_lat2, z_lon1, z_lon2
    LOGICAL :: set_single_triangle

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_Uniform_SpecialArea'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! 2012-10-31: Indonesian Archipelago - connected to Atlantic? (4 cpu)
    z_lat1  =  -5.0_wp
    z_lat2  =  10.0_wp
    z_lon1  = 115.0_wp
    z_lon2  = 135.0_wp
    ! 2012-10-31: corresponding NAtl: 1N 20W, 3 levels
    z_lon1  = 145.0_wp
    z_lon2  = 160.0_wp
    ! 2012-11-08: Test homogen with one cell differ at 10N; 70E
    z_lat1  =   0.0_wp
    z_lat2  =  15.0_wp
    z_lon1  =  60.0_wp
    z_lon2  =  80.0_wp
    
    z_lat1=basin_center_lat
    z_lat2=basin_center_lat +1.0
    
    z_lon1=basin_center_lon
    z_lon2=basin_center_lon+1.0

    set_single_triangle=.false.
    
    ocean_temperature = 5.0_wp !-5.0_wp
    
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = cell_center(idx, block)%lat * rad2deg
        lon_deg = cell_center(idx, block)%lon * rad2deg

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
       
          IF ( (lat_deg >= z_lat1 .AND. lat_deg <= z_lat2) .AND. &
            & (lon_deg >= z_lon1 .AND. lon_deg <= z_lon2) .AND. &
            & ( level <= 1 )) THEN

            IF(.NOT.set_single_triangle)THEN
              ocean_temperature(idx,1:5,block) =  6.0_wp!-4.0_wp
              !set_single_triangle=.true.
!write(1020,*)'indices',idx,block              
            ELSE
            
            ENDIF


          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE temperature_dirac_signal
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !-Square blob for testing advection
  SUBROUTINE tracer_advec_test(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg
    REAL(wp):: z_lat1, z_lat2, z_lon1, z_lon2
    LOGICAL :: set_single_triangle

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_Uniform_SpecialArea'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    z_lat1  =  - 5._wp
    z_lat2  =    5._wp
    z_lon1  =  - 5._wp
    z_lon2  =    5._wp
    
    set_single_triangle=.false.
    
    ocean_temperature = 5.0_wp !-5.0_wp
    
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = cell_center(idx, block)%lat * rad2deg
        lon_deg = cell_center(idx, block)%lon * rad2deg

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

!            write(*,*) level, lat_deg, lon_deg            
       
          IF ( (lat_deg >= z_lat1 .AND. lat_deg <= z_lat2) .AND. &
            & (lon_deg >= z_lon1 .AND. lon_deg <= z_lon2) ) THEN 

            IF(.NOT.set_single_triangle)THEN
              ocean_temperature(idx, level, block) =  7.5_wp!-4.0_wp
              !set_single_triangle=.true.
!write(1020,*)'indices',idx,block              
            ELSE
            
            ENDIF


          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE tracer_advec_test 
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !-Lock Exchange test case
  !! Approximating the following reference
  !! https://doi.org/10.1016/j.ocemod.2014.12.004
  SUBROUTINE lock_exchange(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg
    REAL(wp):: z_lat1, z_lat2, z_lon1, z_lon2
    LOGICAL :: set_single_triangle

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_Uniform_SpecialArea'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    z_lon1  =    0._wp
    z_lat1  =    0._wp
    
    set_single_triangle=.false.
    
    ocean_temperature = 5.0_wp 
    
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = cell_center(idx, block)%lat * rad2deg
        lon_deg = cell_center(idx, block)%lon * rad2deg

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

!          IF (lon_deg >= z_lon1 )  THEN 
         IF (lat_deg >= z_lat1 )  THEN 

            IF(.NOT.set_single_triangle)THEN
              ocean_temperature(idx, level, block) =  30._wp
            ENDIF


          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE lock_exchange 
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !-overflow test case
  !! Approximating the following reference
  !! https://doi.org/10.1016/j.ocemod.2014.12.004
  SUBROUTINE overflow(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: y_loc, y_loc0 
    LOGICAL :: set_single_triangle

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_Uniform_SpecialArea'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' ')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    set_single_triangle=.false.

    y_loc0 = 20000._wp
    
    ocean_temperature = 10.0_wp 
    
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        y_loc = patch_2D%cells%cartesian_center(idx, block)%x(2)

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

         IF (y_loc >= y_loc0 )  THEN 

            IF(.NOT.set_single_triangle)THEN
              ocean_temperature(idx, level, block) =  20._wp
            ENDIF


          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE overflow 
  !-------------------------------------------------------------------------------




  !-------------------------------------------------------------------------------
  SUBROUTINE tracer_Redi_test(patch_3d, ocean_tracer,ocean_state)
  !
  !This test is for testsuite use: it reuqires the density field to be stationary!
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
   TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg
    REAL(wp),POINTER :: density(:,:,:)
    REAL(wp):: slope_parameter =0.5_wp
    REAL(wp) :: x_coord, z_coord
    !REAL(wp) :: upper_level, middle_level, lower_level
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_Redi_test'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' tracer_Redi_test')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    
    density=> ocean_state%p_diag%rho(:,:,:)
    ocean_tracer=0.0_wp
    density=0.0_wp
    slope_parameter =0.00001_wp
!    slope_parameter =0.15_wp !0.5
    
!write(*,*)'surface gradient',temperature_difference
!    lat(:,:) = patch_2d%cells%center(:,:)%lat    
!    level=1
!    DO block = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!      DO idx = start_cell_index, end_cell_index
!        tracer(idx,level,block) = &
!          & initial_temperature_south - temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
!          
!      END DO
!    END DO
!    bottom_value=initial_temperature_bottom 
!
!    DO block = all_cells%start_block, all_cells%end_block
!     CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!     DO idx = start_cell_index, end_cell_index
!        
!       IF (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1) /= 0.0_wp) THEN
!         linear_increase = (ocean_tracer(idx,1,block)- bottom_value  ) / & 
!           & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))
!       ELSE
!         linear_increase = 0.0_wp
!       ENDIF
!       !linear_increase=slope_parameter*linear_increase
!        
!       !First third of levels
!       DO level = 2,7!INT(n_zlev/3.0_wp)+1!patch_3d%p_patch_1d(1)%dolic_c(idx,block) !2,7
!         ocean_tracer(idx,level,block) &
!           & = ocean_tracer(idx,level-1,block) - linear_increase * &
!           &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!       END DO
!
!       !Second third of levels
!       level=8!INT(n_zlev/3.0_wp)+2
!       ocean_tracer(:,level,:)  = 16.5_WP!max(ocean_tracer(:,level-1,:),0.0_wp)! maxval(ocean_tracer(:,level-1,:))
!       DO level = 9,12!INT(n_zlev/3.0_wp)+3,2*INT(n_zlev/3.0_wp)!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
!         ocean_tracer(idx,level,block) &
!           & = ocean_tracer(idx,level-1,block) -0.1_wp! linear_increase * &
!           !            &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!       END DO
!       
!       level=20!n_zlev
!        tracer(idx,level,block) = &
!          & initial_temperature_bottom + temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
!          
!        DO level =19,13,-1
!!n_zlev-1,2*INT(n_zlev/3.0_wp)+1,-1!INT(2.0_wp*n_zlev/3.0_wp)+1, n_zlev!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
!          ocean_tracer(idx,level,block) &
!            & = ocean_tracer(idx,level+1,block) + linear_increase * &
!            &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!        END DO
!
!
!      END DO
!    END DO
!
!    level=1
!    DO block = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!      DO idx = start_cell_index, end_cell_index
!        tracer(idx,level,block) = &
!          & initial_temperature_south - temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
!          
!        density(idx,level,block)=OceanReferenceDensity - LinearThermoExpansionCoefficient*tracer(idx,level,block)  
!
!      END DO
!    END DO
!
!    bottom_value=initial_temperature_bottom 
!
!    DO block = all_cells%start_block, all_cells%end_block
!     CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!     DO idx = start_cell_index, end_cell_index
!        
!       IF (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1) /= 0.0_wp) THEN
!         linear_increase = (ocean_tracer(idx,1,block)- bottom_value  ) / & 
!           & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))
!       ELSE
!         linear_increase = 0.0_wp
!       ENDIF
!       !linear_increase=slope_parameter*linear_increase
!        
!       DO level = 2,n_zlev!INT(n_zlev/3.0_wp)+1!patch_3d%p_patch_1d(1)%dolic_c(idx,block) !2,7
!         !ocean_tracer(idx,level,block) &
!         !  & = ocean_tracer(idx,level-1,block) - linear_increase * &
!         !  &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!
!         density(idx,level,block) &
!           & = density(idx,level-1,block) + 0.5_wp*linear_increase * &
!           &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!           
!       END DO
!
!      END DO
!    END DO
 
 
    density(:,:,:)=1023.0_wp
    ocean_tracer(:,:,:)=2.0_wp
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(idx,block)%lat * rad2deg
        lon_deg = patch_2d%cells%center(idx,block)%lon * rad2deg
        
        x_coord = (lon_deg - (basin_center_lon -0.5_wp*basin_width_deg))/(basin_width_deg)
        IF(lat_deg>=basin_center_lat-0.5.AND.lat_deg<basin_center_lat+0.5)THEN

        DO level = 1,n_zlev

          z_coord =&
          & (patch_3d%p_patch_1d(1)%zlev_m(1)- patch_3d%p_patch_1d(1)%zlev_m(level))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  

           !density(idx,level,block)=density(idx,level,block)&
           !&-tanh(slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !

          !density(idx,level,block)=density(idx,level,block)+tanh(z_coord)*density(idx,level,block)         
          !density(idx,level,block)=1023+density(idx,level,block)
!  write(2040,*)'dens',z_coord,level,density(idx,level,block),tanh(slope_parameter*(z_coord-tanh(x_coord))),&
!&tanh(z_coord),tanh(z_coord)*density(idx,level,block)
!           ENDIF
!          density(idx,level,block)=tanh(5.0_wp*(z_coord-0.0_wp-slope_parameter*8.0_wp*(pi**3)*(x_coord**3)&
!          &*(sin(pi*x_coord)-0.5_wp*sin(2_wp*pi*x_coord))**2)) 
!          density(idx,level,block)=tanh(5.0_wp*(x_coord-0.0_wp-slope_parameter*8.0_wp*(pi**3)*(z_coord**3)&
!          &*(sin(pi*z_coord)-0.5_wp*sin(2_wp*pi*z_coord))**2)) 
          !density(idx,level,block)=-tanh(x_coord)!(5.0_wp*(z_coord-0.25_wp-slope_parameter*8.0_wp*(pi**3)*(x_coord**3)))!&
!          &*(sin(pi*x_coord)-0.5_wp*sin(2_wp*pi*x_coord)))) 

          IF(x_coord<=0.5_wp.and.x_coord>=0.3_wp.and.-z_coord<=0.5_wp.and.-z_coord>=0.3_wp)THEN
          !ocean_tracer(idx,level,block)=0.25_wp*(cos((20_wp*z_coord-5)*pi/3.0_wp)+1.0_wp)&
          !&*(cos((20_wp*x_coord-5_wp)*pi/3.0_wp)+1.0_wp)
          
          ocean_tracer(idx,level,block)=ocean_tracer(idx,level,block)&
          !&exp(-0.01*((z_coord+0.45)**2+(x_coord-0.55)**2))
          &+0.5_wp*(cos((20_wp*(z_coord+0.4_wp))*pi/3.0_wp))&
          &*(cos((20_wp*(x_coord-0.4_wp))*pi/3.0_wp))

  !write(2040,*)'density',x_coord,z_coord,level!,density(idx,level,block),ocean_tracer(idx,level,block)         
          ENDIF
!IF(density(idx,level,block)/=1.0_wp.and. density(idx,level,block)/=-1.0_wp)THEN         
! write(2040,*)'density',x_coord,z_coord,level,density(idx,level,block),ocean_tracer(idx,level,block)  
!ENDIF                                  
        END DO
        ENDIF          
      END DO
    END DO
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(idx,block)%lat * rad2deg
        lon_deg = patch_2d%cells%center(idx,block)%lon * rad2deg
        
        x_coord = (lon_deg - (basin_center_lon -0.5_wp*basin_width_deg))/(basin_width_deg)
        IF(lat_deg>=basin_center_lat-0.5.AND.lat_deg<basin_center_lat+0.5)THEN

        DO level = 1,n_zlev

          z_coord =&
          & (patch_3d%p_patch_1d(1)%zlev_m(1)- patch_3d%p_patch_1d(1)%zlev_m(level))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  

           !density(idx,level,block)=density(idx,level,block)&
           !&-tanh(slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !
!To be preserved: this is used for Redi test
          IF(x_coord+z_coord>=0.25_wp.OR.x_coord+z_coord<=-0.25_wp)THEN
           density(idx,level,block)=density(idx,level,block)+0.000115_wp*density(idx,level,block)!&
          ! &-tanh(0.5*slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !          
          !ENDIF
          !ocean_tracer(idx,level,block)=ocean_tracer(idx,level,block)-0.5_wp*ocean_tracer(idx,level,block)
          !ELSEIF(x_coord+z_coord<=-0.25_wp)THEN
          ! density(idx,level,block)=density(idx,level,block)-0.000115_wp*density(idx,level,block)!&
          ! &-tanh(0.5*slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !          
          ELSE
          density(idx,level,block)=1023_wp&
           !&-tanh(slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) 
           &-tanh(slope_parameter*(z_coord+tanh(x_coord+z_coord)))!this introduces another variation in density
          
          ENDIF

        END DO
        ENDIF          
      END DO
    END DO

DO level = 1, n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
!write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
    CALL dbg_print('trac_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
!     CALL dbg_print('rho init', density(:,level,:), method_name, 3, in_subset=all_cells)
END DO
DO level = 1, n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
!write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
!    CALL dbg_print('trac_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
     CALL dbg_print('rho init', density(:,level,:), method_name, 3, in_subset=all_cells)
END DO
!write(*,*)'leave init'
!    CALL dbg_print('aft. AdvIndivTrac: trac_old', ocean_tracer(:,2,:), method_name, 3, in_subset=all_cells)

!stop
  END SUBROUTINE tracer_Redi_test
  !-------------------------------------------------------------------------------

 !-------------------------------------------------------------------------------
  SUBROUTINE tracer_Redi_test2(patch_3d, ocean_tracer,ocean_state)
  !
  !This test is for testsuite use: it reuqires the density field to be stationary!
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
   TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg
    REAL(wp),POINTER :: density(:,:,:)
    REAL(wp):: slope_parameter =0_wp
    REAL(wp) :: x_coord, z_coord,x_coord_prime
    !REAL(wp) :: upper_level, middle_level, lower_level
    REAL(wp) :: inv_cell_characteristic_length,cell_aspect_ratio
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_Redi_test'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' tracer_Redi_test')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    
    density=> ocean_state%p_diag%rho(:,:,:)
    ocean_tracer=0.0_wp
    density=0.0_wp
!    slope_parameter =0.00001_wp
    slope_parameter =5.0_wp !1.5_wp !0.75_wp !0.15
    
    density(:,:,:)=0.0_wp
    ocean_tracer(:,:,:)=0.0_wp
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
      
        inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(idx,block))
      

        lat_deg = patch_2d%cells%center(idx,block)%lat * rad2deg
        
        x_coord = (lat_deg +0.5*basin_height_deg)/(basin_height_deg)
        x_coord_prime=x_coord !1.0_wp-x_coord
!write(1234,*)'x-coord',x_coord, lat_deg        
        !IF(lat_deg>=basin_center_lat-0.5.AND.lat_deg<basin_center_lat+0.5)THEN

        DO level = 1,n_zlev
          cell_aspect_ratio=patch_3d%p_patch_1d(1)%prism_thick_c(idx,level,block) &
          & * inv_cell_characteristic_length
                
          z_coord =&
          &1+((patch_3d%p_patch_1d(1)%zlev_m(1)- patch_3d%p_patch_1d(1)%zlev_m(level))&
          &/patch_3d%p_patch_1d(1)%zlev_m(n_zlev))
          !z_coord =z_coord-0.5   

          density(idx,level,block)=-TANH(5.0_wp*(z_coord-0.25_wp&
          &+cell_aspect_ratio*slope_parameter*(pi**3)*(x_coord_prime**3)&
          &*(sin(pi*x_coord_prime)-0.5_wp*sin(2.0_wp*pi*x_coord_prime))**2))
          
          !IF((x_coord>=0.1_wp.AND.x_coord<=0.3_wp).AND.(z_coord>=0.1_wp.AND.z_coord<=0.3_wp))THEN
          IF((x_coord>=0.6_wp.AND.x_coord<=0.9_wp).AND.(z_coord>=0.1_wp.AND.z_coord<=0.3_wp))THEN          
            ocean_tracer(idx,level,block)=max(0.25*cos((20_wp*z_coord-5_wp)*pi/3_wp+1)*cos((20_wp*x_coord-5_wp)*pi/3_wp+1),0.0_wp)        
          ELSE
            ocean_tracer(idx,level,block)=0.0_wp
          ENDIF
          ocean_state%p_diag%rho_GM(idx,level,block)=density(idx,level,block)          
          
        END DO
        !ENDIF          
      END DO
    END DO

DO level = 1, n_zlev
    CALL dbg_print('trac_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
END DO
DO level = 1, n_zlev
     CALL dbg_print('rho init', density(:,level,:), method_name, 3, in_subset=all_cells)
END DO

!stop
  END SUBROUTINE tracer_Redi_test2
  !-------------------------------------------------------------------------------



 !-------------------------------------------------------------------------------
  SUBROUTINE tracer_GMR_slope_test(patch_3d, ocean_tracer,ocean_state)
  !
  !This test is for testsuite use: it reuqires the density field to be stationary!
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
   TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg
    REAL(wp),POINTER :: density(:,:,:)
    REAL(wp):: slope_parameter =0.5_wp
    REAL(wp) :: lat_diff
    !REAL(wp) :: upper_level, middle_level, lower_level
    REAL(wp) :: temperature_difference,basin_southBoundary
    REAL(wp), POINTER :: tracer(:,:,:) 
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_Redi_test'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' tracer_Redi_test')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    
    tracer =>ocean_tracer(:,:,:)
    density=> ocean_state%p_diag%rho(:,:,:)
    ocean_tracer=0.0_wp
    density=0.0_wp

!    slope_parameter =0.00001_wp
    slope_parameter =0.5_wp !0.15

    temperature_difference = (initial_temperature_south - initial_temperature_north)
    

    basin_southBoundary    = (basin_center_lat - 0.5_wp*basin_height_deg) * deg2rad
    lat_diff               = basin_height_deg*deg2rad!basin_northBoundary - basin_southBoundary  !  basin_height_deg*deg2rad
 
 
    density(:,:,:)=1023.0_wp
    ocean_tracer(:,:,:)=0.0_wp
    
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(idx,block)%lat !* rad2deg
        slope_parameter=(lat_deg-basin_southBoundary)/lat_diff

        tracer(idx,1,block)=initial_temperature_south- slope_parameter*temperature_difference       

!write(123,*)'data',tracer(idx,1,block),initial_temperature_south, (lat_deg-basin_southBoundary)*temperature_difference/lat_diff,&
!&(lat_deg-basin_southBoundary),slope_parameter,temperature_difference/lat_diff 
        
        DO level = 2,n_zlev
          tracer(idx,level,block)=tracer(idx,level-1,block)&
          &-1.0_wp/(patch_3d%p_patch_1d(1)%zlev_m(level)-patch_3d%p_patch_1d(1)%zlev_m(level-1))
          !z_coord =&
          !& (patch_3d%p_patch_1d(1)%zlev_m(1)- patch_3d%p_patch_1d(1)%zlev_m(level))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  

         END DO
!write(1234,*)'data',tracer(idx,:,block)
         
      END DO
    END DO
 
DO level = 1, 1!n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
!write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
    CALL dbg_print('trac_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
!     CALL dbg_print('rho init', density(:,level,:), method_name, 3, in_subset=all_cells)
END DO
!DO level = 1, n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
!write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
!    CALL dbg_print('trac_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
!     CALL dbg_print('rho init', density(:,level,:), method_name, 3, in_subset=all_cells)
!END DO
!write(*,*)'leave init'
!    CALL dbg_print('aft. AdvIndivTrac: trac_old', ocean_tracer(:,2,:), method_name, 3, in_subset=all_cells)


  END SUBROUTINE tracer_GMR_slope_test
  !-------------------------------------------------------------------------------



  !-------------------------------------------------------------------------------
  SUBROUTINE tracer_GM_test(patch_3d, ocean_tracer,ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
   TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp),POINTER :: density(:,:,:)
    REAL(wp):: slope_parameter =0.5_wp
    REAL(wp) :: x_coord, z_coord
    !REAL(wp) :: upper_level, middle_level, lower_level
!   REAL(wp) :: temperature_difference,basin_northBoundary,basin_southBoundary,lat_diff
    REAL(wp) :: lat_deg, lon_deg
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_Redi_test'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' tracer_Redi_test')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    
    density=> ocean_state%p_diag%rho(:,:,:)
    ocean_tracer=0.0_wp
    density=0.0_wp
!    slope_parameter =0.00001_wp
    slope_parameter =0.05
!   temperature_difference = slope_parameter*(initial_temperature_south - initial_temperature_north)
    
!   basin_northBoundary    = (basin_center_lat + 0.5_wp*basin_height_deg) * deg2rad
!   basin_southBoundary    = (basin_center_lat - 0.5_wp*basin_height_deg) * deg2rad
!   lat_diff               = basin_northBoundary - basin_southBoundary  !  basin_height_deg*deg2rad
!write(*,*)'surface gradient',temperature_difference
!    lat(:,:) = patch_2d%cells%center(:,:)%lat    
!    level=1
!    DO block = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!      DO idx = start_cell_index, end_cell_index
!        tracer(idx,level,block) = &
!          & initial_temperature_south - temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
!          
!      END DO
!    END DO
!    bottom_value=initial_temperature_bottom 
!
!    DO block = all_cells%start_block, all_cells%end_block
!     CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!     DO idx = start_cell_index, end_cell_index
!        
!       IF (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1) /= 0.0_wp) THEN
!         linear_increase = (ocean_tracer(idx,1,block)- bottom_value  ) / & 
!           & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))
!       ELSE
!         linear_increase = 0.0_wp
!       ENDIF
!       !linear_increase=slope_parameter*linear_increase
!        
!       !First third of levels
!       DO level = 2,7!INT(n_zlev/3.0_wp)+1!patch_3d%p_patch_1d(1)%dolic_c(idx,block) !2,7
!         ocean_tracer(idx,level,block) &
!           & = ocean_tracer(idx,level-1,block) - linear_increase * &
!           &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!       END DO
!
!       !Second third of levels
!       level=8!INT(n_zlev/3.0_wp)+2
!       ocean_tracer(:,level,:)  = 16.5_WP!max(ocean_tracer(:,level-1,:),0.0_wp)! maxval(ocean_tracer(:,level-1,:))
!       DO level = 9,12!INT(n_zlev/3.0_wp)+3,2*INT(n_zlev/3.0_wp)!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
!         ocean_tracer(idx,level,block) &
!           & = ocean_tracer(idx,level-1,block) -0.1_wp! linear_increase * &
!           !            &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!       END DO
!       level=20!n_zlev
!        tracer(idx,level,block) = &
!          & initial_temperature_bottom + temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
!          
!        DO level =19,13,-1
!!n_zlev-1,2*INT(n_zlev/3.0_wp)+1,-1!INT(2.0_wp*n_zlev/3.0_wp)+1, n_zlev!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
!          ocean_tracer(idx,level,block) &
!            & = ocean_tracer(idx,level+1,block) + linear_increase * &
!            &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!        END DO
!      END DO
!    END DO
!
!    level=1
!    DO block = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!      DO idx = start_cell_index, end_cell_index
!        tracer(idx,level,block) = &
!          & initial_temperature_south - temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
!          
!        density(idx,level,block)=OceanReferenceDensity - LinearThermoExpansionCoefficient*tracer(idx,level,block)  
!
!      END DO
!    END DO
!
!    bottom_value=initial_temperature_bottom 
!
!    DO block = all_cells%start_block, all_cells%end_block
!     CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!     DO idx = start_cell_index, end_cell_index
!        
!       IF (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1) /= 0.0_wp) THEN
!         linear_increase = (ocean_tracer(idx,1,block)- bottom_value  ) / & 
!           & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))
!       ELSE
!         linear_increase = 0.0_wp
!       ENDIF
!       !linear_increase=slope_parameter*linear_increase
!        
!       DO level = 2,n_zlev!INT(n_zlev/3.0_wp)+1!patch_3d%p_patch_1d(1)%dolic_c(idx,block) !2,7
!         !ocean_tracer(idx,level,block) &
!         !  & = ocean_tracer(idx,level-1,block) - linear_increase * &
!         !  &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!
!         density(idx,level,block) &
!           & = density(idx,level-1,block) + 0.5_wp*linear_increase * &
!           &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!           
!       END DO
!      END DO
!    END DO
 
 
    density(:,:,:)=OceanReferenceDensity
    ocean_tracer(:,:,:)=2.0_wp
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(idx,block)%lat * rad2deg
        lon_deg = patch_2d%cells%center(idx,block)%lon * rad2deg
        
        x_coord = (lon_deg - (basin_center_lon -0.5_wp*basin_width_deg))/(basin_width_deg)
        IF(lat_deg>=basin_center_lat-0.5.AND.lat_deg<basin_center_lat+0.5)THEN

        DO level = 1,n_zlev

          z_coord =&
          & (patch_3d%p_patch_1d(1)%zlev_m(1)- patch_3d%p_patch_1d(1)%zlev_m(level))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  

           !density(idx,level,block)=density(idx,level,block)&
           !&-tanh(slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !

          !density(idx,level,block)=density(idx,level,block)+tanh(z_coord)*density(idx,level,block)         
          !density(idx,level,block)=1023+density(idx,level,block)
!  write(2040,*)'dens',z_coord,level,density(idx,level,block),tanh(slope_parameter*(z_coord-tanh(x_coord))),&
!&tanh(z_coord),tanh(z_coord)*density(idx,level,block)
!           ENDIF
!          density(idx,level,block)=tanh(5.0_wp*(z_coord-0.0_wp-slope_parameter*8.0_wp*(pi**3)*(x_coord**3)&
!          &*(sin(pi*x_coord)-0.5_wp*sin(2_wp*pi*x_coord))**2)) 
!          density(idx,level,block)=tanh(5.0_wp*(x_coord-0.0_wp-slope_parameter*8.0_wp*(pi**3)*(z_coord**3)&
!          &*(sin(pi*z_coord)-0.5_wp*sin(2_wp*pi*z_coord))**2)) 
          !density(idx,level,block)=-tanh(x_coord)!(5.0_wp*(z_coord-0.25_wp-slope_parameter*8.0_wp*(pi**3)*(x_coord**3)))!&
!          &*(sin(pi*x_coord)-0.5_wp*sin(2_wp*pi*x_coord)))) 

           ocean_tracer(idx,level,block)=20.0_wp&
           &-tanh(slope_parameter*(-z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) 
           !&-tanh(slope_parameter*(z_coord+tanh(x_coord+z_coord)))!this introduces another variation in density
 

!          IF(x_coord<=0.5_wp.and.x_coord>=0.3_wp.and.-z_coord<=0.5_wp.and.-z_coord>=0.3_wp)THEN
!          !ocean_tracer(idx,level,block)=0.25_wp*(cos((20_wp*z_coord-5)*pi/3.0_wp)+1.0_wp)&
!          !&*(cos((20_wp*x_coord-5_wp)*pi/3.0_wp)+1.0_wp)
!          
!          ocean_tracer(idx,level,block)=ocean_tracer(idx,level,block)&
!          !&exp(-0.01*((z_coord+0.45)**2+(x_coord-0.55)**2))
!          &+0.5_wp*(cos((20_wp*(z_coord+0.4_wp))*pi/3.0_wp))&
!          &*(cos((20_wp*(x_coord-0.4_wp))*pi/3.0_wp))
!          !ocean_tracer(idx,level,block)=max(ocean_tracer(idx,level,block),0.0_wp)
!  !write(2040,*)'density',x_coord,z_coord,level!,density(idx,level,block),ocean_tracer(idx,level,block)         
!          ENDIF

        END DO
        ENDIF          
      END DO
    END DO
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(idx,block)%lat * rad2deg
        lon_deg = patch_2d%cells%center(idx,block)%lon * rad2deg
        
        x_coord = (lon_deg - (basin_center_lon -0.5_wp*basin_width_deg))/(basin_width_deg)
        IF(lat_deg>=basin_center_lat-0.5.AND.lat_deg<basin_center_lat+0.5)THEN

        DO level = 1,n_zlev

          z_coord =&
          & (patch_3d%p_patch_1d(1)%zlev_m(1)- patch_3d%p_patch_1d(1)%zlev_m(level))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  

           !density(idx,level,block)=density(idx,level,block)&
           !&-tanh(slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !
         IF(x_coord+1.5_wp+4*z_coord>=0.5)THEN
           ocean_tracer(idx,level,block)=ocean_tracer(idx,level,block)+0.115_wp*ocean_tracer(idx,level,block)!&
!          ! &-tanh(0.5*slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !          
!          !ENDIF
          ELSEIF(x_coord+1.5_wp+4_wp*z_coord<=-0.5_wp)THEN
           ocean_tracer(idx,level,block)=ocean_tracer(idx,level,block)-0.115_wp*ocean_tracer(idx,level,block)!&
!          ! &-tanh(0.5*slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !          
!!          !ELSE
!!          !density(idx,level,block)=1023_wp&
!!          ! !&-tanh(slope_parameter*(z_coord+tanh(x_coord))) !-tanh(x_coord)*(1.0_wp+z_coord) !
!!          ! &-tanh(slope_parameter*(z_coord+tanh(x_coord+z_coord)))
          ENDIF


        END DO
        ENDIF          
      END DO
    END DO

density=OceanReferenceDensity - LinearThermoExpansionCoefficient * ocean_tracer
DO level = 1, n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
!write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
    CALL dbg_print('trac_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
!     CALL dbg_print('rho init', density(:,level,:), method_name, 3, in_subset=all_cells)
END DO
DO level = 1, n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
!write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
!    CALL dbg_print('trac_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
     CALL dbg_print('rho init', density(:,level,:), method_name, 3, in_subset=all_cells)
END DO
!write(*,*)'leave init'
!    CALL dbg_print('aft. AdvIndivTrac: trac_old', ocean_tracer(:,2,:), method_name, 3, in_subset=all_cells)

!stop
  END SUBROUTINE tracer_GM_test
  !-------------------------------------------------------------------------------


 !-------------------------------------------------------------------------------
  SUBROUTINE tracer_Redi_test_withdensity(patch_3d, ocean_tracer,ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp),POINTER :: density(:,:,:)
    REAL(wp):: slope_parameter =0.5_wp
    REAL(wp) :: temperature_difference,basin_northBoundary,basin_southBoundary,lat_diff,bottom_value
    REAL(wp) :: lat(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks),linear_increase  
    REAL(wp), POINTER :: tracer(:,:,:) 
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_Redi_test'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' tracer_Redi_test')

    patch_2d    => patch_3d%p_patch_2d(1)
    all_cells   => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center
    
    tracer  => ocean_tracer(:,:,:)
    density => ocean_state%p_diag%rho(:,:,:)
    ocean_tracer=5.0_wp
    density=OceanReferenceDensity
    slope_parameter =0.5_wp
    temperature_difference = slope_parameter*(initial_temperature_south - initial_temperature_north)
    
    basin_northBoundary    = (basin_center_lat + 0.5_wp*basin_height_deg) * deg2rad
    basin_southBoundary    = (basin_center_lat - 0.5_wp*basin_height_deg) * deg2rad
    lat_diff               = basin_northBoundary - basin_southBoundary  !  basin_height_deg*deg2rad
!write(*,*)'surface gradient',temperature_difference
    lat(:,:) = patch_2d%cells%center(:,:)%lat    
    level=1
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        tracer(idx,level,block) = &
          & initial_temperature_south - temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
          
        density(idx,level,block)=OceanReferenceDensity - LinearThermoExpansionCoefficient*tracer(idx,level,block)  
!write(123,*)'density:temp', density(idx,level,block), tracer(idx,level,block)      
      END DO
    END DO
! stop       
    bottom_value=initial_temperature_bottom 

    DO block = all_cells%start_block, all_cells%end_block
     CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
     DO idx = start_cell_index, end_cell_index
        
       IF (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1) /= 0.0_wp) THEN
         linear_increase = (ocean_tracer(idx,1,block)- bottom_value  ) / & 
           & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))
       ELSE
         linear_increase = 0.0_wp
       ENDIF
       !linear_increase=slope_parameter*linear_increase
        
       DO level = 2,n_zlev!INT(n_zlev/3.0_wp)+1!patch_3d%p_patch_1d(1)%dolic_c(idx,block) !2,7
         !ocean_tracer(idx,level,block) &
         !  & = ocean_tracer(idx,level-1,block) - linear_increase * &
         !  &     patch_3d%p_patch_1d(1)%del_zlev_i(level)

         density(idx,level,block) &
           & = density(idx,level-1,block) + 0.5_wp*linear_increase * &
           &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
           
       END DO

      END DO
    END DO
 
   ocean_tracer(:,:,:)=5.0_wp
   DO block = all_cells%start_block, all_cells%end_block
     CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
     DO idx = start_cell_index, end_cell_index
        
        
         level = 3!INT(n_zlev/2.0_wp)
         !This criterion singles out one cell: a dirac signal: 0.03 only one cell
         IF(     abs(cell_center(idx,block)%lat-basin_center_lat*deg2rad)<0.003_wp&
          & .AND.abs(cell_center(idx,block)%lon-basin_center_lon*deg2rad)<0.003_wp )THEN

         ocean_tracer(idx,level,block)= ocean_tracer(idx,level,block) + 0.5_wp
         ENDIF

!         DO level = 4,n_zlev!INT(n_zlev/3.0_wp)+1!patch_3d%p_patch_1d(1)%dolic_c(idx,block) !2,7
!         ocean_tracer(idx,level,block) &
!           & = ocean_tracer(idx,level-1,block) - linear_increase * &
!           &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
!         END DO
           

      END DO
    END DO
!

 !rho(1:levels) = OceanReferenceDensity - LinearThermoExpansionCoefficient * t(1:levels)
DO level = 1, n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
!write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
    CALL dbg_print('Init trac_old', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
END DO
!write(*,*)'leave init'
!    CALL dbg_print('aft. AdvIndivTrac: trac_old', ocean_tracer(:,2,:), method_name, 3, in_subset=all_cells)
DO level = 1, n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
!write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
    CALL dbg_print('Init: density', density(:,level,:), method_name, 3, in_subset=all_cells)
END DO


  END SUBROUTINE tracer_Redi_test_withdensity
  !-------------------------------------------------------------------------------



  !-------------------------------------------------------------------------------
  SUBROUTINE tracer_Redi_test_withdensity0(patch_3d, ocean_tracer,ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
   TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp),POINTER :: density(:,:,:)
    !REAL(wp) :: x_coord, z_coord,linear_increase,linear_decrease
    !REAL(wp) :: upper_level, middle_level, lower_level
    REAL(wp) :: temperature_difference,basin_eastBoundary,lat_diff,bottom_value
    REAL(wp) :: basin_northBoundary,linear_increase
    REAL(wp) :: lat(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: lon(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)  
    REAL(wp), POINTER :: tracer(:,:,:) 
    REAL(wp) :: x_1, x_3,xi
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_Redi_test'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' tracer_Redi_test')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    
    tracer =>ocean_tracer(:,:,:)
    density => ocean_state%p_diag%rho

    ocean_tracer=0.0_wp
    density     =0.0_wp
    xi=0.01_wp

    
    basin_eastBoundary    = (basin_center_lat + 0.5_wp*basin_width_deg) * deg2rad

!write(*,*)'surface gradient',temperature_difference
    lat(:,:) = patch_2d%cells%center(:,:)%lat    
    lon(:,:) = patch_2d%cells%center(:,:)%lon  

    basin_eastBoundary = maxval(lon)
    basin_northBoundary    = (basin_center_lat + 0.5_wp*basin_height_deg) * deg2rad
        
    
    
    level=1
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        tracer(idx,level,block) = &
          & initial_temperature_south - temperature_difference*((basin_northBoundary-lat(idx,block))/lat_diff)
          
        density(idx,level,block)=OceanReferenceDensity - LinearThermoExpansionCoefficient*tracer(idx,level,block)  
!write(123,*)'density:temp', density(idx,level,block), tracer(idx,level,block)      
      END DO
    END DO
! stop       
    bottom_value=initial_temperature_bottom 

    DO block = all_cells%start_block, all_cells%end_block
     CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
     DO idx = start_cell_index, end_cell_index
        
       IF (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1) /= 0.0_wp) THEN
         linear_increase = (ocean_tracer(idx,1,block)- bottom_value  ) / & 
           & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))
       ELSE
         linear_increase = 0.0_wp
       ENDIF
       !linear_increase=slope_parameter*linear_increase
        
       DO level = 2,n_zlev!INT(n_zlev/3.0_wp)+1!patch_3d%p_patch_1d(1)%dolic_c(idx,block) !2,7
         !ocean_tracer(idx,level,block) &
         !  & = ocean_tracer(idx,level-1,block) - linear_increase * &
         !  &     patch_3d%p_patch_1d(1)%del_zlev_i(level)

         density(idx,level,block) &
           & = density(idx,level-1,block) + 1.0_wp*linear_increase * &
           &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
           
       END DO

      END DO
    END DO
    
    
    
    
          
!write(0,*)'max',maxval(lon),minval(lon),basin_westBoundary,basin_eastBoundary
    DO block = all_cells%start_block, all_cells%end_block
     CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
     DO idx = start_cell_index, end_cell_index
        
        
       DO level = 1,n_zlev
       
         x_3 = 1.0_wp-(patch_3d%p_patch_1d(1)%zlev_m(level)-patch_3d%p_patch_1d(1)%zlev_m(1))/&
         &(patch_3d%p_patch_1d(1)%zlev_m(n_zlev)-patch_3d%p_patch_1d(1)%zlev_m(1))
         
!        x_1 = (lon(idx,block)-basin_westBoundary)/(basin_eastBoundary-basin_westBoundary)
         x_1 = lon(idx,block)/basin_eastBoundary   
                
         density(idx,level,block) =-tanh( 5.0_wp*(x_3-0.25_wp&
         &-xi*8.0_wp*(pi*x_1**3)*(sin(pi*x_1)-0.5_wp*sin(2.0_wp*pi*x_1))**2))
         IF(( 0.1<x_1) .and. (x_1<0.4) .and. (0.1<x_3) .and. (x_3<0.4))THEN
         tracer(idx,level,block)=0.25_wp*(cos(20.0_wp*x_3-5)*(pi/3.0_wp)+1)&
         &*(cos(20.0_wp*x_1-5)*(pi/3.0_wp)+1)
         ELSE
         tracer(idx,level,block)=0.0_wp
         ENDIF
write(12345,*)'init', level,x_1**3,x_3,&
!&patch_3d%p_patch_1d(1)%zlev_m(level),&
&(sin(pi*x_1)-0.5_wp*sin(2.0_wp*pi*x_1))**2,&
&(x_3-0.25_wp-xi*8.0_wp*(pi**3)*(x_1**3)*(sin(pi*x_1)-0.5_wp*sin(2.0_wp*pi*x_1))**2),&        
&ocean_tracer(idx,level,block),density(idx,level,block)
       END DO

       

      END DO
    END DO
    
DO level = 1, n_zlev
!z_coord = (patch_3d%p_patch_1d(1)%zlev_m(level)- patch_3d%p_patch_1d(1)%zlev_m(1))/patch_3d%p_patch_1d(1)%zlev_m(n_zlev)  
write(*,*)'temp',level,maxval(ocean_tracer(:,level,:)),minval(ocean_tracer(:,level,:))                
    CALL dbg_print('Initial: density', density(:,level,:), method_name, 3, in_subset=all_cells)
    CALL dbg_print('Initial: tracer', tracer(:,level,:), method_name, 3, in_subset=all_cells)
END DO
!write(*,*)'leave init'
!    CALL dbg_print('aft. AdvIndivTrac: trac_old', ocean_tracer(:,2,:), method_name, 3, in_subset=all_cells)

!stop
  END SUBROUTINE tracer_Redi_test_withdensity0
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_Willebrand_test(patch_3d, ocean_tracer)
  !
  !This test is for testsuite use: it reuqires the density field to be stationary!
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg, z_tmp
    REAL(wp),POINTER :: density(:,:,:)
    REAL(wp):: slope_parameter =0_wp
    real(wp):: centerline
    REAL(wp) :: x_coord, z_coord,z_coord_prime,x_coord_prime,width
    REAL(wp) :: left_basin_boundary_lon, right_basin_boundary_lon
    !REAL(wp) :: upper_level, middle_level, lower_level
    REAL(wp) :: temperature_difference,basin_northBoundary,basin_southBoundary,lat_diff,bottom_value
    REAL(wp) :: lat(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: inv_cell_characteristic_length,cell_characteristic_length, cell_aspect_ratio
    REAL(wp), POINTER :: tracer(:,:,:) 
    REAL(wp) :: center_point
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_Redi_test'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' temperature_Willebrand_test')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center
    lat(:,:) = patch_2d%cells%center(:,:)%lat    
    
    tracer =>ocean_tracer(:,:,:)
    ocean_tracer(:,:,:)=0.0_wp
    
    basin_northBoundary    = (basin_center_lat + 0.5_wp*basin_height_deg) * deg2rad
    basin_southBoundary    = (basin_center_lat - 0.5_wp*basin_height_deg) * deg2rad
    lat_diff               = basin_northBoundary - basin_southBoundary  !  basin_height_deg*deg2rad
 
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
      
        inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(idx,block))
        cell_characteristic_length     = SQRT(patch_2D%cells%area(idx,block))
      
        lat = patch_2d%cells%center(idx,block)%lat 
        
        x_coord=(lat(idx,block)-basin_southBoundary)/lat_diff
        
        x_coord_prime=1.0_wp-x_coord!1.0_wp-x_coord
        
        DO level = 1,n_zlev!-1

           z_coord =1.0_wp+&
           &((patch_3d%p_patch_1d(1)%zlev_m(1)- patch_3d%p_patch_1d(1)%zlev_m(level))&
           &/patch_3d%p_patch_1d(1)%zlev_m(n_zlev))
           
           z_coord_prime=1.0_wp-z_coord
           !width=0.099_wp !0.15_wp
           width=0.15_wp
           centerline=0.5_wp+  0.25_wp*tanh((x_coord-0.5_wp)/1.5_wp) !For Redi test =0.6, for GMR shallower=1.5  0.4_wp 0.13 0.2: steep slope; 0.5 small slope !!x_o=0.5
           center_point=centerline!0.5
           ocean_tracer(idx,level,block)=exp(-( (z_coord-center_point)**2 )*40_wp) !For Redi test=70, for GMR=40

        END DO
      END DO
    END DO

DO level = 1, n_zlev
    CALL dbg_print('temperature_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
END DO
!stop


  END SUBROUTINE temperature_Willebrand_test
  !-------------------------------------------------------------------------------



  !-------------------------------------------------------------------------------
!  SUBROUTINE tracer_ConstantSurface_IncludeLand(patch_3d, ocean_tracer, top_value, bottom_value)
!    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!    REAL(wp), TARGET :: ocean_tracer(:,:,:)
!    REAL(wp), INTENT(in) :: top_value, bottom_value
!
!    TYPE(t_patch),POINTER   :: patch_2d
!    TYPE(t_subset_range), POINTER :: all_cells
!
!    INTEGER :: block, idx, level
!    INTEGER :: start_cell_index, end_cell_index
!    REAL(wp) :: linear_increase
!
!    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_ConstantSurface'
!    !-------------------------------------------------------------------------
!
!    patch_2d => patch_3d%p_patch_2d(1)
!    all_cells => patch_2d%cells%ALL
!
!    DO block = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!      DO idx = start_cell_index, end_cell_index
!          ocean_tracer(idx,1,block) = top_value
!      END DO
!    END DO
!
!    linear_increase = (bottom_value - top_value) / (REAL(n_zlev,wp)-1.0_wp)
!
!    ! write(0,*) n_zlev
!    ! write(0,*) bottom_value, top_value, linear_increase
!
!    DO block = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
!      DO idx = start_cell_index, end_cell_index
!
!        DO level = 2, n_zlev
!          ocean_tracer(idx,level,block) &
!            & = ocean_tracer(idx,level-1,block) + linear_increase
!        !  write(0,*) ocean_tracer(idx,level,block), ocean_tracer(idx,level-1,block)
!        END DO
!      END DO
!    END DO
!
!    ! ocean_tracer(:,:,:) = top_value
!
!  END SUBROUTINE tracer_ConstantSurface_IncludeLand
  !-------------------------------------------------------------------------------
  SUBROUTINE tracer_quads_checkerboard(patch_3d, ocean_tracer, base_value, variation)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in) :: base_value, variation

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER  :: block, idx, level, start_cell_index, end_cell_index
    INTEGER  :: checkerboard_top_mod
    REAL(wp) :: checkerboard_mod

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_quads_checkerboard'

    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        checkerboard_top_mod = MODULO((block - all_cells%start_block) * nproma + idx, 2)
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
          checkerboard_mod = (REAL(MODULO(checkerboard_top_mod + level, 2),wp) - 0.5_wp) * 2.0_wp ! this is -1,+1
          ocean_tracer(idx,level,block) = base_value + checkerboard_mod * variation
        ENDDO
      END DO
    END DO

  END SUBROUTINE tracer_quads_checkerboard
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE height_quads_checkerboard(patch_3d, ocean_height, base_value, variation)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_height(:,:)
    REAL(wp), INTENT(in) :: base_value, variation

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER  :: block, idx, level, start_cell_index, end_cell_index
    INTEGER  :: checkerboard_top_mod
    REAL(wp) :: checkerboard_mod, column_sign

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':height_quads_checkerboard'

    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    column_sign = -1
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        IF (MODULO((block - all_cells%start_block) * nproma + idx - 1, 4) == 0) THEN
          column_sign = -column_sign
        ENDIF
        checkerboard_top_mod = MODULO((block - all_cells%start_block) * nproma + idx, 2)
        checkerboard_mod = column_sign * (REAL(checkerboard_top_mod,wp) - 0.5_wp) * 2.0_wp ! this is -1,+1
        
        write(0,*) block, idx, " checkerboard_mod=", checkerboard_mod
        DO level = 1, MIN(patch_3d%p_patch_1d(1)%dolic_c(idx,block),1)
          ocean_height(idx,block) = base_value + checkerboard_mod * variation
        ENDDO
      END DO
    END DO

  END SUBROUTINE height_quads_checkerboard
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE tracer_ConstantSurface(patch_3d, ocean_tracer, top_value)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in) :: top_value

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: block, idx, level, start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_ConstantSurface'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        DO level = 1, MIN(1,patch_3d%p_patch_1d(1)%dolic_c(idx,block))
          ocean_tracer(idx,level,block) = top_value
        ENDDO
      END DO
    END DO

  END SUBROUTINE tracer_ConstantSurface
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  SUBROUTINE fillVerticallyMissingValues(patch_3d, ocean_tracer, has_missValue, missValue)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    LOGICAL :: has_missValue
    REAL(wp) :: missValue
      
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    !-------------------------------------------------------------------------
    IF (.NOT. has_missValue) RETURN
    
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index        
        DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
        
          IF (ocean_tracer(idx,level,block) == missValue) &
            & ocean_tracer(idx,level,block) = ocean_tracer(idx,level-1,block)
        
        END DO
      END DO
    END DO

  END SUBROUTINE fillVerticallyMissingValues
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! decrease tvertically linerarly the given tracer based on the top level value
  ! of the tracer and using a decres of (top_value - bottom_value) / (n_zlev - 1)
  SUBROUTINE increaseTracerVerticallyLinearly(patch_3d, ocean_tracer, bottom_value, start_level,end_level)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in), OPTIONAL :: bottom_value
    INTEGER, OPTIONAL    :: start_level
    INTEGER, OPTIONAL    :: end_level

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: linear_increase
    INTEGER    :: top_level
    INTEGER    :: bottom_level
    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

   IF(PRESENT(start_level))THEN
      top_level=start_level
    ELSE
      top_level=2
    ENDIF
    IF(PRESENT(end_level))THEN
      bottom_level=end_level
    ELSE
      bottom_level=n_zlev
    ENDIF

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        
        IF (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1) /= 0.0_wp) THEN
          linear_increase = (bottom_value - ocean_tracer(idx,1,block) ) / & 
            & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))
        ELSE
          linear_increase = 0.0_wp
        ENDIF

        DO level = top_level, bottom_level!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
          ocean_tracer(idx,level,block) &
            & = ocean_tracer(idx,level-1,block) + linear_increase * &
            &     patch_3d%p_patch_1d(1)%del_zlev_i(level)
        END DO

      END DO
    END DO

  END SUBROUTINE increaseTracerVerticallyLinearly
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  ! decrease tvertically linerarly the given tracer based on the top level value
  ! of the tracer and using a decres of (top_value - bottom_value) / (n_zlev - 1)
  SUBROUTINE de_increaseTracerVertically(patch_3d, ocean_tracer,&
  & decrease_start_level,decrease_end_level, increase_start_level,increase_end_level)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    INTEGER  :: decrease_start_level
    INTEGER :: decrease_end_level
    INTEGER :: increase_start_level
    INTEGER :: increase_end_level

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: linear_increase

    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL


    !decrease
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        linear_increase = 0.001_wp
        !(ocean_tracer(idx,decrease_start_level,block) - 0.1_wp*ocean_tracer(idx,decrease_end_level,block) ) / & 
        !& (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))

        DO level = decrease_start_level, decrease_end_level !n_zlev!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
          ocean_tracer(idx,level,block) &
            & = ocean_tracer(idx,level,block) - linear_increase *ocean_tracer(idx,level,block)
        END DO
        
        DO level = decrease_end_level,increase_start_level !n_zlev!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
          ocean_tracer(idx,level,block) &
            & = ocean_tracer(idx,level,block) 
        END DO
        
        !linear_increase = 0.5_wp*(ocean_tracer(idx,decrease_start_level,block) - &
        !  0.1_wp*ocean_tracer(idx,decrease_end_level,block) ) / & 
        !  & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))
       
        DO level = increase_end_level, increase_start_level, -1 !n_zlev!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
          ocean_tracer(idx,level,block) &
            & = ocean_tracer(idx,level,block) + linear_increase*ocean_tracer(idx,level,block) 
        END DO
        
        
      END DO
    END DO

    DO level = 1, n_zlev !n_zlev!patch_3d%p_patch_1d(1)%dolic_c(idx,block)
     write(0,*)'in',level,maxval( ocean_tracer(:,level,:)),minval( ocean_tracer(:,level,:))
      END DO

!stop

  END SUBROUTINE de_increaseTracerVertically
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! pertrube the tracer according to the cos(lon)
  SUBROUTINE perturbeTracer_LonCosinus(patch_3d, ocean_tracer, waveNumber, max_ratio)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in), OPTIONAL :: waveNumber, max_ratio

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index

    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index      
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
        
          ocean_tracer(idx,level,block) = ocean_tracer(idx,level,block) * &
            & (1.0_wp + max_ratio * COS(waveNumber * patch_2d%cells%center(idx,block)%lon))
          
        END DO
      END DO
    END DO

  END SUBROUTINE perturbeTracer_LonCosinus
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! pertrube the tracer according to the cos(lon)
  SUBROUTINE perturbeTracer_LatCosinus(patch_3d, ocean_tracer, waveNumber, max_ratio)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in), OPTIONAL :: waveNumber, max_ratio

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index

    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL


    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

          ocean_tracer(idx,level,block) = ocean_tracer(idx,level,block) * &
            & (1.0_wp + max_ratio * COS(waveNumber * patch_2d%cells%center(idx,block)%lat))

        END DO
      END DO
    END DO

  END SUBROUTINE perturbeTracer_LatCosinus
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! decrease tvertically linerarly the given tracer based on the top level value
  ! of the tracer and using a decres of (top_value - bottom_value) / (n_zlev - 1)
  SUBROUTINE increaseTracerLevelsLinearly(patch_3d, ocean_tracer, bottom_value, increase_gradient)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in), OPTIONAL :: bottom_value, increase_gradient

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: linear_increase

    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    linear_increase = 0.0_wp
    

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        IF (PRESENT(bottom_value)) &
          & linear_increase = (bottom_value - ocean_tracer(idx,1,block) ) / (REAL(n_zlev,wp)-1.0_wp)
      
        DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
          IF (PRESENT(increase_gradient)) &
            & linear_increase = increase_gradient * patch_3d%p_patch_1D(1)%prism_center_dist_c(idx,level,block)
         
          ocean_tracer(idx,level,block) = ocean_tracer(idx,level-1,block) + linear_increase
        END DO

      END DO
    END DO

  END SUBROUTINE increaseTracerLevelsLinearly
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !!  exponential temperature profile following Abernathey et al., 2011
  SUBROUTINE varyTracerVerticallyExponentially(patch_3d, ocean_tracer, bottom_value, scale_depth)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), INTENT(in) :: bottom_value
    REAL(wp), INTENT(in) :: scale_depth

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: bottom_depth, temperature_difference, exp_neghoverH, exp_negzoverH

    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    bottom_depth = patch_3d%p_patch_1d(1)%zlev_m(n_zlev)      !  below H is T=T_bot
    exp_neghoverH = exp((-1.0_wp)*bottom_depth/scale_depth)   !  constant: 1-e**(-H/h)

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        temperature_difference = ocean_tracer(idx,1,block) - bottom_value

        DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

          exp_negzoverH = exp((-1.0_wp)*patch_3d%p_patch_1d(1)%zlev_m(level)/scale_depth)
          ocean_tracer(idx,level,block) = bottom_value + &
            &  temperature_difference*(exp_negzoverH - exp_neghoverH) / (1.0_wp - exp_neghoverH)
     !      &   EXP((-1.0_wp)*bottom_depth/scale_depth) /                          &
     !      &   (1.0_wp - exp((-1.0_wp)*bottom_depth/scale_depth)))
     !      &  (EXP((-1.0_wp)*patch_3d%p_patch_1d(1)%del_zlev_i(level)/scale_depth) - &
     !      &   EXP((-1.0_wp)*bottom_depth/scale_depth) /                          &
     !      &   (1.0_wp - exp((-1.0_wp)*bottom_depth/scale_depth)))
        END DO

      END DO
    END DO

  END SUBROUTINE varyTracerVerticallyExponentially
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_CollapsingDensityFront_WeakGrad(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, z_tmp
    REAL(wp):: z_ldiff, z_ltrop, z_lpol
    REAL(wp):: z_ttrop, z_tpol, z_tdiff

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_CollapsingDensityFront_WeakGrad'
    !-------------------------------------------------------------------------

    CALL message(method_name, ' Collapsing density front with weaker gradient')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    ! Temperature profile in first layer depends on latitude only
    ! Construct temperature profile
    !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
    z_ttrop = 10.0_wp      ! tropical temperature
    z_tpol  =  5.0_wp      ! polar temperature
    z_ttrop =  5.0_wp      ! 2011-09-02: instable stratification
    z_tpol  = 10.0_wp      ! 2011-09-02: instable stratification
    z_lpol  = 70.0_wp      ! polar boundary latitude of transition zone

    z_ttrop = 25.0_wp      ! 2011-09-05: stable stratification
    z_tpol  = 10.0_wp      ! 2011-09-05: stable stratification
    z_ltrop = 15.0_wp      ! tropical boundary latitude of transition zone
    z_lpol  = 60.0_wp      ! polar boundary latitude of transition zone
    z_tdiff = z_ttrop - z_tpol
    z_ldiff = z_lpol  - z_ltrop

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(idx,block)%lat * rad2deg

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

              IF(ABS(lat_deg)>=z_lpol)THEN

                ocean_temperature(idx,level,block) =  z_tpol

                ! ocean_state%p_diag%temp_insitu(idx,level,block) = z_tpol
                !            ocean_state%p_prog(nold(1))%tracer(idx,level,block,1)&
                !             &= 30.0_wp!convert_insitu2pot_temp_func(ocean_state%p_diag%temp_insitu(idx,level,block),&
                !                     &ocean_state%p_prog(nold(1))%tracer(idx,level,block,2),&
                !                     &sfc_press_bar)
                !SItodBar*OceanReferenceDensity*v_base%zlev_m(level))!1013.0_wp)SItodBar*101300.0_wp)!

              ELSEIF (ABS(lat_deg)<=z_ltrop) THEN

                ! ocean_state%p_diag%temp_insitu(idx,level,block) = z_ttrop
                ocean_temperature(idx,level,block) = z_ttrop


              ELSE ! IF(ABS(lat_deg)<z_lpol .AND. ABS(lat_deg)>z_ltrop)THEN
                !   z_tmp = pi*((abs(lat_deg) - z_lpol)/z_ldiff)
                !   ocean_state%p_diag%temp_insitu(idx,level,block) = z_tpol + 0.5_wp*z_tdiff*(1.0_wp+cos(z_tmp))
                z_tmp = 0.5_wp*pi*((ABS(lat_deg) - z_ltrop)/z_ldiff)
                ! ocean_state%p_diag%temp_insitu(idx,level,block) = z_ttrop - z_tdiff*SIN(z_tmp)
                ocean_temperature(idx,level,block) = z_ttrop - z_tdiff*SIN(z_tmp)
                !      if (level==1) write(*,*) 'zlat,ztmp(deg),temp', &
                !   &  block,idx,lat_deg,(abs(lat_deg)-z_lpol)/z_ldiff,ocean_state%p_diag%temp_insitu(idx,level,block)
              ENDIF
!            ELSE
!              ! ocean_state%p_diag%temp_insitu(idx,level,block) = z_tdeep
!              ocean_temperature(idx,level,block) = z_tdeep
!            ENDIF  ! level=1
        END DO
      END DO
    END DO

  END SUBROUTINE temperature_CollapsingDensityFront_WeakGrad
  !-------------------------------------------------------------------------------



  !-------------------------------------------------------------------------------
  ! Temperature profile depends on latitude and depth
  ! Construct temperature profile
  !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
  !   for maximum tropical temperature see tprof
  SUBROUTINE temperature_KelvinHelmholtzTest(patch_3d, ocean_temperature, top_value, bottom_value)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)
    REAL(wp) :: top_value, bottom_value

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp) :: shear_depth, shear_center, shear_top,shear_bottom 
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_KelvinHelmholtzTest'
    !-------------------------------------------------------------------------

    patch_2d    => patch_3d%p_patch_2d(1)
    all_cells   => patch_2d%cells%ALL

    
    shear_depth  = 4.0_wp
    shear_center = INT(0.5_wp*n_zlev)
    shear_top    = shear_center-INT(0.5_wp*shear_depth)
    shear_bottom = shear_center+INT(0.5_wp*shear_depth)
    
    top_value    = 10.0_wp
    bottom_value = 5.0_wp
    
!         linear_increase = (bottom_value - top_value ) / & 
!           & (patch_3d%p_patch_1d(1)%zlev_m(n_zlev) - patch_3d%p_patch_1d(1)%zlev_m(1))


    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
      
      
        !Above shear layer
        DO level = 1, INT(shear_top)
          ocean_temperature(idx,level,block) = top_value
        ENDDO
        
        !Shear layer
        DO level = INT(shear_top+1),INT(shear_bottom-1)
          ocean_temperature(idx,level,block) = ocean_temperature(idx,level-1,block)-(top_value-bottom_value)/shear_depth

        ENDDO
        
        !Below shear layer
        DO level = INT(shear_bottom),n_zlev
          ocean_temperature(idx,level,block) = bottom_value
        ENDDO
      END DO
    END DO
!      DO level = 1, n_zlev
!       write(*,*)'ocean_temperature', level, maxval(ocean_temperature(:, level,:))   
!      ENDDO

  END SUBROUTINE temperature_KelvinHelmholtzTest
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! Temperature profile depends on latitude and depth
  ! Construct temperature profile
  !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
  !   for maximum tropical temperature see tprof
  SUBROUTINE temperature_TropicsPolar(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, z_tmp
    REAL(wp):: z_ldiff, z_ltrop, z_lpol
    REAL(wp):: z_ttrop, z_tpol, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_TropicsPolar'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    ! Temperature profile depends on latitude and depth
    ! Construct temperature profile
    !   ttrop for lat<ltrop; tpol for lat>lpol; cos for transition zone
    !   for maximum tropical temperature see values above
    CALL message(method_name, ': simple tropics-pol/vertical temperature profile')

    IF (i_sea_ice == 0) THEN
      z_tpol  =  5.0_wp      ! polar temperature
    ELSE
      z_tpol = tf
    ENDIF
    z_ltrop = 15.0_wp      ! tropical latitude for temperature gradient
    z_lpol  = 60.0_wp      ! polar latitude for temperature gradient
    z_ldiff = z_lpol  - z_ltrop

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        lat_deg = patch_2d%cells%center(idx,block)%lat * rad2deg

        ! bugfix: z_tpol was 0 for 10 levels since level was inner loop
        !         does not effect 4 levels
        z_tpols = z_tpol
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

          z_ttrop = tprof(level)
          z_tpols = MIN(z_tpols,tprof(level))

          IF(ABS(lat_deg)>=z_lpol)THEN

             ocean_temperature(idx,level,block) = z_tpols

          ELSEIF(ABS(lat_deg)<=z_ltrop)THEN

            ocean_temperature(idx,level,block) = z_ttrop

          ELSE ! IF(ABS(lat_deg)<z_lpol .AND. ABS(lat_deg)>z_ltrop)THEN

            z_tdiff = z_ttrop - z_tpols
            z_tmp = 0.5_wp*pi*((ABS(lat_deg) - z_ltrop)/z_ldiff)
            ocean_temperature(idx,level,block) = z_ttrop - z_tdiff * SIN(z_tmp)

          ENDIF

        END DO
      END DO
    END DO

  END SUBROUTINE temperature_TropicsPolar
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_circularLonLatPerturbation(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: distan
    REAL(wp):: perturbation_lat, perturbation_lon, perturbation_width
    REAL(wp):: temperature

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_circularLonLatPerturbation'
    !-------------------------------------------------------------------------

    CALL message(method_name, ':')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center


    perturbation_lat= pi / 6._wp*rad2deg
    perturbation_lon=-pi_2*rad2deg
    !perturbation_lat = basin_center_lat! + 0.1_wp*basin_height_deg!             !45.5_wp
    !perturbation_lon =  0.0_wp!0.1_wp*basin_width_deg                           !4.5_wp
    !max_perturbation  = 20.0_wp            !20.1_wp
    perturbation_width  =  7.0_wp*pi/64.0_wp !10.0_wp!5.0_wp!1.5_wp

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        distan=SQRT((cell_center(idx,block)%lat - perturbation_lat * deg2rad)**2 + &
          & (cell_center(idx,block)%lon - perturbation_lon*deg2rad)**2)

        !Local hot perturbation
        IF(distan<=perturbation_width)THEN
          temperature = (1.0_wp+COS(pi*distan/perturbation_width)) * 0.5_wp +2.0_wp
        ELSE
          temperature = 0.0_wp
        ENDIF

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
          ocean_temperature(idx, level, block) = temperature
        END DO

      END DO
    END DO
    
  END SUBROUTINE temperature_circularLonLatPerturbation
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_DanilovsMunkGyre(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: levels
    REAL(wp):: distan
    REAL(wp):: perturbation_lat, perturbation_lon,  max_perturbation, perturbation_width
    ! REAL(wp):: z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_DanilovsMunkGyre'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    !------------------------------
    CALL message(method_name, ': Danilovs Munk gyre flow')

    perturbation_lat = basin_center_lat !- 0.1_wp * basin_height_deg
    perturbation_lon = basin_center_lon !- 0.1_wp * basin_width_deg
    max_perturbation  = 1.0!0.1_wp!20.1_wp
    perturbation_width  = 2.0_wp!1.5_wp

    ! Next update 2011-05-24: due to Danilov the perturbation should be -1 Kelvin, width 3.0
    ! 05-25: max and width larger: -2.0 and 5.0
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        levels = patch_3d%p_patch_1d(1)%dolic_c(idx,block)

        IF (levels > 0) THEN
          ! level=1:  250m  T= 20 - 0.9375 = 19.0625
          ! level=2:  750m  T= 20 - 2.8125 = 17.1875
          ! level=3: 1250m  T= 20 - 4.6875 = 15.3125
          ! level=4: 1750m  T= 20 - 6.5625 = 13.4375
          ocean_temperature(idx,1:levels,block) = 20.0_wp
                  !stratification
          DO level = 1, levels
            ocean_temperature(idx,level,block) =          &
              & ocean_temperature(idx,level,block)-0.5_wp*level          
          END DO
                  
                  
          distan = SQRT((cell_center(idx,block)%lat - perturbation_lat * deg2rad)**2 + &
            & (cell_center(idx,block)%lon - perturbation_lon * deg2rad)**2)

  !Commented out PK 4/2015        !Local cold perturbation
          IF(distan<=5.0_wp)THEN
            DO level = 1, 1!levels
              ocean_temperature(idx,level,block) =          &
                & ocean_temperature(idx,level,block)        &
                & - max_perturbation*EXP(-(distan/(perturbation_width*deg2rad))**2) !&
             !                &   * sin(pi*v_base%zlev_m(level)/4000.0_wp)!&
             !   & * SIN(pi*patch_3d%p_patch_1d(1)%zlev_m(level) / patch_3d%p_patch_1d(1)%zlev_i(levels+1))
!write(123,*)'perturb',max_perturbation*EXP(-(distan/(perturbation_width*deg2rad))**2)                   
            END DO
          ENDIF !Local hot perturbation

        END IF !(levels > 0)
      END DO
    END DO
  END SUBROUTINE temperature_DanilovsMunkGyre
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_BasinWithVerticalWall(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_BasinWithVerticalWall'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Adjusting density front in a basin: vertical wall at basin_center_lon
    CALL message(method_name,' Adjusting density front in a basin with vertical wall')

    !Impose temperature profile. Profile
    !depends on latitude only and is uniform across
    !all vertical layers
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        IF (cell_center(idx,block)%lon >= basin_center_lon) THEN
          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
            ocean_temperature(idx,1:n_zlev,block) = 30.0_wp
          ENDDO
        ELSE
          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
            ocean_temperature(idx,1:n_zlev,block) = 25.0_wp
          ENDDO
        ENDIF
      END DO
    END DO

  END SUBROUTINE temperature_BasinWithVerticalWall
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_CollapsingDensityFront_StuhnePeltier(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, z_tmp
    ! REAL(wp):: perturbation_lat, perturbation_lon,  z_ltrop, z_lpol
    ! REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_CollapsingDensityFront_StuhnePeltier'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    CALL message(method_name, ': Collapsing density front, Stuhne-Peltier')

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(idx,block)%lat * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across
        !all vertical layers
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

          IF (ABS(lat_deg) >= 40.0_wp) THEN

            ocean_temperature(idx,level,block) = 5.0_wp

          ELSEIF (ABS(lat_deg) <= 20.0_wp) THEN

            ocean_temperature(idx,level,block) =  30.0_wp

          ELSE ! IF (ABS(lat_deg) < 40.0_wp .AND. ABS(lat_deg) > 20.0_wp)THEN

            z_tmp = pi*((ABS(lat_deg) -20.0_wp)/20.0_wp)

            ocean_temperature(idx,level,block) = &
              & 5.0_wp + 0.5_wp * 25.0_wp * (1.0_wp + COS(z_tmp))

          ENDIF

        END DO
      END DO
    END DO

   END SUBROUTINE temperature_CollapsingDensityFront_StuhnePeltier
  !-------------------------------------------------------------------------------



  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_GM_idealized(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: BLOCK, idx, level, ll
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg!, lon_deg, z_tmp
    ! REAL(wp):: perturbation_lat, perturbation_lon,  z_ltrop, z_lpol
    ! REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_GM_idealized'
    !-------------------------------------------------------------------------
    REAL(wp) :: scal, delta_t_back, tano
    
    ! initialisation with stable background stratification and a latitude dependend t and/or s  anomaly
    ! amplitude of the anomaly is decreasing with depth

    delta_t_back=1.0_wp ! increase per level
    tano=0.0_wp


    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    CALL message(method_name, ': Collapsing density front, Stuhne-Peltier')

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(idx,block)%lat * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across
        !all vertical layers
        DO level = 1, n_zlev

          ll=n_zlev+1-level

          ocean_temperature(idx,level,BLOCK)=0.0_wp
          IF (ABS(lat_deg) <= 45.0_wp) THEN

           scal=(COS(lat_deg/45.0_wp * pi) +1.0_wp) *0.5_wp
           ocean_temperature(idx,level,BLOCK) =0.0  + delta_t_back*ll + scal*ll*tano

          ELSE

           ocean_temperature(idx,level,BLOCK) =0.0  + delta_t_back*ll
 
          ENDIF

        END DO
      END DO
    END DO

   END SUBROUTINE temperature_GM_idealized
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE salinity_GM_idealized(patch_3d, ocean_salinity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: BLOCK, idx, level, ll
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg!, lon_deg, z_tmp
    ! REAL(wp):: perturbation_lat, perturbation_lon,  z_ltrop, z_lpol
    ! REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':salinity_GM_idealized'
    !-------------------------------------------------------------------------
     REAL(wp) :: scal, delta_s_back, sano
    
    ! initialisation with stable background stratification and a latitude dependend t and/or s  anomaly
    ! amplitude of the anomaly is decreasing with depth

    delta_s_back=0.1_wp ! increase per level
    sano=0.01_wp


    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    CALL message(method_name, ': Collapsing density front, Stuhne-Peltier')

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(idx,block)%lat * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across
        !all vertical layers
        DO level = 1, n_zlev

          ll=n_zlev+1-level

          ocean_salinity(idx,level,BLOCK)=0.0_wp
          IF (ABS(lat_deg) <= 45.0_wp) THEN

           scal=(COS(lat_deg/45.0_wp *pi) +1.0_wp) *0.5_wp
           !scal=(COS((lat_deg*deg2rad)/5.0_wp ) +1.0_wp) *0.5_wp
           ocean_salinity(idx,level,BLOCK) =35.0  + delta_s_back*ll + scal*ll*sano

          ELSE

           ocean_salinity(idx,level,BLOCK) =35.0  + delta_s_back*ll
 
          ENDIF

        END DO
      END DO
    END DO

   END SUBROUTINE salinity_GM_idealized
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE temperature_GM_idealized2(patch_3d, ocean_temperature)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: BLOCK, idx, level, ll
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg!, lon_deg, z_tmp
    ! REAL(wp):: perturbation_lat, perturbation_lon,  z_ltrop, z_lpol
    ! REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_GM_idealized'
    !-------------------------------------------------------------------------
    REAL(wp) :: delta_t_back, north, south

    ! initialisation with stable background stratification and a latitude dependend t and/or s  anomaly

    delta_t_back=1.0_wp ! increase per level

    north=45.0_wp
    south=-45.0_wp

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    CALL message(method_name, ': Collapsing density front, Stuhne-Peltier')

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(idx,block)%lat * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across
        !all vertical layers
        DO level = 1, n_zlev

          ll=n_zlev+1-level

          ocean_temperature(idx,level,BLOCK)=0.0_wp
          IF (lat_deg < south) THEN

            ocean_temperature(idx,level,BLOCK) =0.0  + delta_t_back*ll

          ELSEIF(lat_deg > north) THEN

            ocean_temperature(idx,level,BLOCK) =0.0  + delta_t_back*ll

          ELSEIF(lat_deg <= north .AND. lat_deg >= south) THEN

            ocean_temperature(idx,level,BLOCK) =0.0  + delta_t_back*ll

          ENDIF

        END DO
      END DO
    END DO

   END SUBROUTINE temperature_GM_idealized2
  !-------------------------------------------------------------------------------

  SUBROUTINE salinity_GM_idealized2(patch_3d, ocean_salinity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: BLOCK, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg!, lon_deg, z_tmp
    ! REAL(wp):: perturbation_lat, perturbation_lon,  z_ltrop, z_lpol
    ! REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':salinity_GM_idealized'
    !-------------------------------------------------------------------------
    REAL(wp) :: north, south, sano, scal

    ! initialisation with stable background stratification and a latitude dependend t and/or s  anomaly

    sano=0.1
    north=45.0_wp
    south=-45.0_wp

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    CALL message(method_name, ': Collapsing density front, Stuhne-Peltier')

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(idx,block)%lat * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across
        !all vertical layers
        DO level = 1, n_zlev

          ocean_salinity(idx,level,BLOCK)=0.0_wp
          IF (lat_deg < south) THEN

            ocean_salinity(idx,level,BLOCK) =35.0_wp

          ELSEIF(lat_deg > north) THEN

            ocean_salinity(idx,level,BLOCK) =35.0_wp + sano

          ELSEIF(lat_deg <= north .AND. lat_deg >= south) THEN

            scal=(north - lat_deg) / (north - south)
            ocean_salinity(idx,level,BLOCK) =  35.0_wp + sano - scal*sano

          ENDIF

        END DO
      END DO
    END DO

   END SUBROUTINE salinity_GM_idealized2
  !-------------------------------------------------------------------------------

  SUBROUTINE salinity_GM_idealized3(patch_3d, ocean_salinity)



    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: BLOCK, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg!, distan !, z_tmp
    ! REAL(wp):: perturbation_lat, perturbation_lon,  z_ltrop, z_lpol
    ! REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':salinity_GM_idealized'
    !-------------------------------------------------------------------------
    REAL(wp) :: a,b,c,xlon,alon_0,alat_0,height,sssu

    ! initialisation with stable background stratification and a latitude dependend t and/or s  anomaly

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    CALL message(method_name, ': Collapsing density front, Stuhne-Peltier')

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(idx,block)%lat * rad2deg
        lon_deg = cell_center(idx,block)%lon * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across

        alat_0=-30.0_wp
        !alat_0=45.0_wp        
        alon_0=250_wp
        !alon_0=60_wp

        !all vertical layers
        DO level = 1, n_zlev


          ocean_salinity(idx,level,BLOCK)=0.0_wp

          a=(lat_deg-alat_0)**2
          xlon=MERGE(lon_deg,lon_deg+360,lon_deg.GE.0)

          b= (xlon-alon_0)**2
          !b= (lon_deg-alon_0)**2          
          c= 10_wp**2

          height=0.0_wp
          sssu=0.8_wp*EXP(- ( a + b ) / c)

          IF (patch_3d%p_patch_1d(1)%zlev_m(level) .LE. height) THEN

            ocean_salinity(idx,level,BLOCK) =34.0_wp

          ELSEIF(patch_3d%p_patch_1d(1)%zlev_m(level) .ge. 1400.0_wp) THEN

            ocean_salinity(idx,level,BLOCK) =35.0_wp

          ELSE


!         distan=SQRT((cell_center(idx, block)%lat*rad2deg + 5.0_wp)**2 + &
!            & (xlon - 180_wp)**2)
 !write(1020,*)'dist', distan,lat_deg,lon_deg,10.0_wp * deg2rad          
          !IF(distan < 10_wp) THEN
 

            ocean_salinity(idx,level,BLOCK) =  34.1_wp + sssu*patch_3d%p_patch_1d(1)%zlev_m(level)/1400.0_wp
          ! ENDIF 
          ENDIF

        END DO
      END DO
    END DO
!stop
   END SUBROUTINE salinity_GM_idealized3
  !-------------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------------
  SUBROUTINE salinity_Willebrand_test(patch_3d, ocean_tracer)
  !
  !This test is for testsuite use: it reuqires the density field to be stationary!
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg, z_tmp
    REAL(wp),POINTER :: density(:,:,:)
    REAL(wp):: slope_parameter =0_wp
    real(wp):: centerline
    REAL(wp) :: x_coord, z_coord,z_coord_prime,x_coord_prime,width
    REAL(wp) :: left_basin_boundary_lon, right_basin_boundary_lon
    !REAL(wp) :: upper_level, middle_level, lower_level
    REAL(wp) :: temperature_difference,basin_northBoundary,basin_southBoundary,lat_diff,bottom_value
    REAL(wp) :: lat(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: inv_cell_characteristic_length,cell_characteristic_length, cell_aspect_ratio
    REAL(wp), POINTER :: tracer(:,:,:) 
    REAL(wp) :: center_point
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_Redi_test'
    !-------------------------------------------------------------------------

    CALL message(TRIM(method_name), ' temperature_Willebrand_test')

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center
    lat(:,:) = patch_2d%cells%center(:,:)%lat    
    
    tracer =>ocean_tracer(:,:,:)
    ocean_tracer(:,:,:)=0.0_wp
    
    basin_northBoundary    = (basin_center_lat + 0.5_wp*basin_height_deg) * deg2rad
    basin_southBoundary    = (basin_center_lat - 0.5_wp*basin_height_deg) * deg2rad
    lat_diff               = basin_northBoundary - basin_southBoundary  !  basin_height_deg*deg2rad
 
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
      
        inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(idx,block))
        cell_characteristic_length     = SQRT(patch_2D%cells%area(idx,block))
      
        lat = patch_2d%cells%center(idx,block)%lat 
        
        x_coord=(lat(idx,block)-basin_southBoundary)/lat_diff
        
        x_coord_prime=1.0_wp-x_coord!1.0_wp-x_coord
        
        DO level = 1,n_zlev!-1

           z_coord =1.0_wp+&
           &((patch_3d%p_patch_1d(1)%zlev_m(1)- patch_3d%p_patch_1d(1)%zlev_m(level))&
           &/patch_3d%p_patch_1d(1)%zlev_m(n_zlev))
           
           z_coord_prime=1.0_wp-z_coord
           !width=0.099_wp !0.15_wp
           width=0.15_wp
           centerline=0.5_wp+  0.25_wp*tanh((x_coord-0.5_wp)/1.5_wp) !For Redi test =0.6, for GMR shallower=1.5  0.4_wp 0.13 0.2: steep slope; 0.5 small slope !!x_o=0.5
           center_point=centerline!0.5
           ocean_tracer(idx,level,block)=exp(-( (z_coord-center_point)**2 )*40_wp) !For Redi test=70, for GMR=40

        END DO
      END DO
    END DO

DO level = 1, n_zlev
    CALL dbg_print('temperature_init', ocean_tracer(:,level,:), method_name, 3, in_subset=all_cells)
END DO
!stop


  END SUBROUTINE salinity_Willebrand_test
  !-------------------------------------------------------------------------------


  SUBROUTINE temperature_GM_idealized3(patch_3d, ocean_temperature)



    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: BLOCK, idx, level, ll
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg !, z_tmp
    ! REAL(wp):: perturbation_lat, perturbation_lon,  z_ltrop, z_lpol
    ! REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_GM_idealized'
    !-------------------------------------------------------------------------
    REAL(wp) :: a,b,c,xlon,alon_0,alat_0,sssu,delta_t_back,tano

    ! initialisation with stable background stratification and a latitude dependend t and/or s  anomaly

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    delta_t_back=0.2_wp ! increase per level
    tano=0.0_wp


    CALL message(method_name, ': Collapsing density front, Stuhne-Peltier')

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(idx,block)%lat * rad2deg
        lon_deg = cell_center(idx,block)%lon * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across



        alat_0=-30.0_wp
        alon_0=250.0_wp


        !all vertical layers
        DO level = 1, n_zlev

          ll=n_zlev+1-level

          ocean_temperature(idx,level,BLOCK)=0.0_wp

          a=(lat_deg-alat_0)**2
          xlon=MERGE(lon_deg,lon_deg+360,lon_deg.GE.0)
          b= (xlon-alon_0)**2
          c= 10_wp**2

          sssu=0.8_wp*EXP(- ( a + b ) / c)

          ocean_temperature(idx,level,BLOCK) =8.0_wp  + delta_t_back*ll+tano*sssu

        END DO
      END DO
    END DO

   END SUBROUTINE temperature_GM_idealized3
  !-------------------------------------------------------------------------------
  
 !-------------------------------------------------------------------------------
  SUBROUTINE temperature_GM_idealized4(patch_3d, ocean_temperature)


    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_temperature(:,:,:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: BLOCK, idx, level, ll
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lat_deg, lon_deg !, z_tmp
    ! REAL(wp):: perturbation_lat, perturbation_lon,  z_ltrop, z_lpol
    ! REAL(wp):: z_ttrop, z_tpol, z_tdeep, z_tdiff, z_tpols

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':temperature_GM_idealized'
    !-------------------------------------------------------------------------
    REAL(wp) :: a,b,c,xlon,alon_0,alat_0,sssu,delta_t_back,tano

    ! initialisation with stable background stratification and a latitude dependend t and/or s  anomaly

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    delta_t_back=2.0_wp * 10.0/real(n_zlev) ! increase per level
    tano=0.0_wp


    CALL message(method_name, ': Collapsing density front, Stuhne-Peltier')

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transer to latitude in degrees
        lat_deg = cell_center(idx,block)%lat * rad2deg
        lon_deg = cell_center(idx,block)%lon * rad2deg
        !Impose emperature profile. Profile
        !depends on latitude only and is uniform across



        alat_0=-30.0_wp
        !alat_0=45.0_wp
        alon_0=250.0_wp
        !alon_0=60_wp

        !all vertical layers
        DO level = 1, n_zlev

          ll=n_zlev+1-level

          ocean_temperature(idx,level,BLOCK)=0.0_wp

          a=(lat_deg-alat_0)**2
          xlon=MERGE(lon_deg,lon_deg+360,lon_deg.GE.0)
          b= (xlon-alon_0)**2
          c=  10_wp**2

          sssu=0.8_wp*EXP(- ( a + b ) / c)

          ocean_temperature(idx,level,BLOCK) =0.0_wp  + delta_t_back*ll+tano*sssu


        END DO
      END DO
    END DO

  END SUBROUTINE temperature_GM_idealized4
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE salinity_AnalyticSmoothVerticalProfile(patch_3d, ocean_salinity)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_salinity(:,:,:)

    REAL(wp) :: salinity_profile(n_zlev)
    INTEGER :: level

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':salinity_AnalyticSmoothVerticalProfile'
    !-------------------------------------------------------------------------

    CALL message(method_name,' Creating analytic profile:')

    DO level=1,n_zlev

      salinity_profile(level) = &
        & MIN(34.1_wp + LOG(1.3_wp + SQRT(patch_3d%p_patch_1d(1)%zlev_m(level)) * 0.05), 35.0_wp)

      ! write(0,*) level, patch_3D%p_patch_1D(1)%zlev_m(level), " salinity:", salinity_profile(level)
      WRITE(message_text,*) level, patch_3d%p_patch_1d(1)%zlev_m(level), " salinity:", salinity_profile(level)
      CALL message(method_name, message_text)
    ENDDO

    CALL fill_FromVerticalArrayProfile(patch_3d=patch_3d, ocean_tracer=ocean_salinity, VerticalProfileValue=salinity_profile)

  END SUBROUTINE salinity_AnalyticSmoothVerticalProfile
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE fill_FromVerticalArrayProfile(patch_3d, ocean_tracer, VerticalProfileValue)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET :: ocean_tracer(:,:,:)
    REAL(wp), TARGET :: VerticalProfileValue(:)

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index

    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    !------------------------------
    ! assign from adhoc array values
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
          ocean_tracer(idx,level,block) = VerticalProfileValue(level)
        END DO
      END DO
    END DO

  END SUBROUTINE fill_FromVerticalArrayProfile
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  ! Initial datum for height, test case 2 of Williamson et al.(1992).
  !
  !
  FUNCTION test2_h( point_lon, point_lat) result(p_hh)
    REAL(wp), PARAMETER :: h0 = 2.94e4_wp * rgrav  ! maximum height
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lon     ! longitude of point
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh      ! height
    
    ! !LOCAL VARIABLES:
    REAL(wp)             :: z_fact1   ! 1st factor
    REAL(wp)             :: z_fact2   ! 2nd factor
    
    ! 1st factor
    
    z_fact1 = sphere_radius * grid_angular_velocity
    z_fact1 = z_fact1 + 0.5_wp * u0
    z_fact1 = z_fact1 * u0 * rgrav
    
    ! 2nd factor
    
    z_fact2 = SIN(point_lat) * COS(aleph)
    z_fact2 = z_fact2 - COS(point_lon) * COS(point_lat) * SIN(aleph)
    z_fact2 = z_fact2 * z_fact2
    
    ! height
    
    p_hh = h0 - z_fact1 * z_fact2
    
  END FUNCTION test2_h
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  ! Initial datum for height, test case 5 of Williamson et al.(1992).
  !
  FUNCTION test5_h( point_lat) result(p_hh)
    REAL(wp), PARAMETER :: h0    = 5960._wp  ! maximum height
    REAL(wp), PARAMETER :: uzero = 20._wp    ! maximum velocity
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh      ! height
    
    ! !LOCAL VARIABLES:
    REAL(wp)             :: z_fact1   ! 1st factor
    REAL(wp)             :: z_fact2   ! 2nd factor
    

    ! 1st factor
    
    z_fact1 = sphere_radius * grid_angular_velocity
    z_fact1 = z_fact1 + 0.5_wp * uzero
    z_fact1 = z_fact1 * uzero * rgrav
    
    ! 2nd factor
    
    z_fact2 = SIN(point_lat)
    z_fact2 = z_fact2 * z_fact2
    
    ! height
    
    p_hh = h0 - z_fact1 * z_fact2
    
  END FUNCTION test5_h
  !-------------------------------------------------------------------------
  

  !-----------------------------------------------------------------------------------
  !
  ! Initial datum for orography, test case 5 of Williamson et al.(1992).
  !
  !
  SUBROUTINE depth_mountain_orography_Williamson_test5(patch_3d, cells_bathymetry)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET, INTENT(inout)  :: cells_bathymetry(:,:)

    REAL(wp), PARAMETER :: h_s0  = 2000._wp  ! maximum height of mountain

    REAL(wp) :: point_lon     ! longitude of point
    REAL(wp) :: point_lat     ! latitude of point

    REAL(wp)             :: point_height      ! orography

    REAL(wp)             :: z_lon_mc  ! Mountain center, longitude ...
    REAL(wp)             :: z_lat_mc  !          ... and latitude
    REAL(wp)             :: z_rad_mt  ! radius of mountain
    REAL(wp)             :: z_dist_mc ! distance from mountain center
    REAL(wp)             :: z_diff    ! difference of coordinates
    REAL(wp)             :: z_min_dist_sq ! min of square of distances

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    ! Local Variables
    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':depth_mountain_orography_Williamson_test5'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    ! center and radius of mountain
    z_lon_mc = -pi_2
    z_lat_mc = pi / 6._wp
    z_rad_mt = pi / 9._wp

 !   patch_3d%lsm_c(:,:,:) = sea
 !   patch_3d%lsm_e(:,:,:) = sea

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        point_lon = patch_2d%cells%center(idx, block)%lon
        point_lat = patch_2d%cells%center(idx, block)%lat

        ! square of distance (in geographical coordinate sense) of point
        ! from mountain center

        z_diff = point_lon - z_lon_mc
        z_diff = z_diff * z_diff
        z_dist_mc = z_diff

        z_diff = point_lat - z_lat_mc
        z_diff = z_diff * z_diff
        z_dist_mc = z_dist_mc + z_diff

        ! if point inside mountain range take its distance, else take mountain radius

        z_diff = z_rad_mt * z_rad_mt
        z_min_dist_sq = MIN ( z_diff, z_dist_mc)
        z_dist_mc = SQRT( z_min_dist_sq)

        ! conical shape of mountain, depending on distance from mountain center

        point_height = z_dist_mc / z_rad_mt
        point_height = 1._wp - point_height

        cells_bathymetry(idx, block) = h_s0 * point_height

      ENDDO
    ENDDO

  END SUBROUTINE depth_mountain_orography_Williamson_test5
  !-----------------------------------------------------------------------------------

  
  !-----------------------------------------------------------------------------------
  !
  ! Initial datum for overflow test case 
  !
  SUBROUTINE depth_overflow(patch_3d, cells_bathymetry)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET, INTENT(inout)  :: cells_bathymetry(:,:)

    REAL(wp)             :: point_height      ! orography

    REAL(wp)             :: d1, d2, d, sigma, x0, x 

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    ! Local Variables
    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':depth_mountain_orography_Williamson_test5'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    d1 =  500._wp 
    d2 = 2000._wp 
    
    sigma =  7000._wp
    x0    = 40000._wp

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        x = patch_2D%cells%cartesian_center(idx, block)%x(2)

        d = d1 + 0.5_wp*(d2 - d1)*(1.0_wp + tanh( (x - x0)/sigma) )
        
        IF ( cells_bathymetry(idx, block) < 0.0_wp ) THEN ! Only touch it if its not land
          cells_bathymetry(idx, block) = -d
        ENDIF


      ENDDO
    ENDDO

  END SUBROUTINE depth_overflow
  !-----------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------
  SUBROUTINE depth_uniform(patch_3d, cells_bathymetry)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    REAL(wp), TARGET, INTENT(inout)  :: cells_bathymetry(:,:)


    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx
    INTEGER :: start_cell_index, end_cell_index

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':depth_uniform'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        cells_bathymetry(idx, block) = topography_height_reference

      ENDDO
    ENDDO
    
  END SUBROUTINE depth_uniform
  !-----------------------------------------------------------------------------------
  

  !-----------------------------------------------------------------------------------
  FUNCTION test_usbr_h( point_lon, point_lat, p_t) result(p_hh)
    !
    ! !DESCRIPTION:
    ! Initial datum for height h, test case unsteady solid body
    ! rotation of L\"auter et al.(2007).
    
    !
    ! !DEFINED PARAMETERS:
    REAL(wp), PARAMETER :: d0    = 133681.0_wp  ! additive constant
    
    ! !INPUT PARAMETERS:
    REAL(wp), INTENT(in) :: point_lon     ! longitude of point
    REAL(wp), INTENT(in) :: point_lat     ! latitude of point
    REAL(wp), INTENT(in) :: p_t       ! point of time
    
    ! !RETURN VALUE:
    REAL(wp)             :: p_hh      ! height
    
    ! !LOCAL VARIABLES:
    REAL(wp)             :: z_phi_t_k ! 1st summand
    REAL(wp)             :: z_summand ! 2nd summand
    REAL(wp)             :: z_fact    ! factor
    REAL(wp)             :: z_angle1  ! 1st angle
    REAL(wp)             :: z_angle2  ! 2nd angle
    !-----------------------------------------------------------------------
    
    z_angle1 = .25_wp * pi
    z_angle2 = point_lon + grid_angular_velocity * p_t
    
    ! 1st summand: \phi_t(\vec c) \cdot \vec k
    
    z_phi_t_k = SIN(point_lat) * COS(z_angle1)
    z_phi_t_k = z_phi_t_k - COS(z_angle2) * COS(point_lat) * SIN(z_angle1)
    z_phi_t_k = u0 * z_phi_t_k
    
    ! 2nd summand: r_e \grid_angular_velocity \sin \varphi
    
    z_summand = sphere_radius * grid_angular_velocity * SIN(point_lat)
    
    ! one factor
    
    z_fact    = .5_wp *  z_phi_t_k + z_summand
    
    ! height
    
    p_hh      = d0 - z_phi_t_k *  z_fact
    p_hh      = p_hh * rgrav
    !write(*,*)'param:', u0, pi, rgrav,sphere_radius, grid_angular_velocity
    !stop
  END FUNCTION test_usbr_h
  !-------------------------------------------------------------------------
    

!   !-------------------------------------------------------------------------
!   !
!   ! Performs  numerical integration between -$\frac{\pi}{2}$ and $\frac{\pi}{2}$
!   ! to compute geostrophically balanced initial state used
!   ! in test 3.
!   !
!   FUNCTION geostrophic_balance( point_lat, func)  result(p_hh)
!     
!     INTERFACE                        ! selected function
!       
!       FUNCTION func(p_t) result(p_vv)
!         
!         USE mo_kind, ONLY: wp
!         
!         REAL(wp), INTENT(in) :: p_t
!         REAL(wp)             :: p_vv
!         
!       END FUNCTION func
!       
!     END INTERFACE
!     
!     ! !INPUT PARAMETERS:
!     REAL(wp), INTENT(in) :: point_lat           ! rotated latitude
!     
!     ! !RETURN VALUE:
!     REAL(wp)             :: p_hh            ! balanced height
!     
!     ! !LOCAL VARIABLES:
!     INTEGER :: j               ! loop index
!     
!     REAL(wp)             :: z_a             ! left bound
!     REAL(wp)             :: z_b             ! right bound
!     REAL(wp)             :: cell_lat           ! latitude in loop
!     REAL(wp)             :: z_step          ! step
!     REAL(wp)             :: z_val, z_val2   ! intermediate values
!     !-----------------------------------------------------------------------
!     
!     z_a = -1._wp * pi_2
!     z_b = point_lat
!     
!     z_step = 0.02_wp * ( z_b - z_a)
!     
!     p_hh = 0._wp
!     
!     cell_lat = z_a - 0.5_wp * z_step
!     
!     DO j = 1, 50
!       cell_lat = cell_lat + z_step
!       
!       z_val = func(cell_lat)
!       
!       z_val2 = 2._wp * grid_angular_velocity * SIN(cell_lat)
!       z_val2 = z_val2 + z_val * TAN(cell_lat)* sphere_radius
!       z_val2 = z_val * z_val2
!       
!       p_hh = p_hh + z_val2 * z_step
!       
!     ENDDO
!     
!   END FUNCTION geostrophic_balance
!   !-------------------------------------------------------------------------
!     
!   
!   !-----------------------------------------------------------------------------------
!   ! Performs  numerical integration between -$\frac{\pi}{2}$ and $\frac{\pi}{2}$
!   ! to compute geostrophically balanced initial state used
!   ! in test 3.
!   !
!   !
!   FUNCTION geostrophic_balance_11( phi, func)  result(p_hh)
!     
!     INTERFACE                        ! selected function
!       
!       FUNCTION func(p_t) result(p_vv)
!         
!         USE mo_kind, ONLY: wp
!         
!         REAL(wp), INTENT(in) :: p_t
!         REAL(wp)             :: p_vv
!         
!       END FUNCTION func
!       
!     END INTERFACE
!     
!     ! !INPUT PARAMETERS:
!     REAL(wp), INTENT(in) :: phi           ! rotated latitude
!     ! !RETURN VALUE:
!     REAL(wp)             :: p_hh            ! balanced height
!     ! !LOCAL VARIABLES:
!     INTEGER :: j               ! loop index
!     REAL(wp)             :: phi_a             ! left bound
!     REAL(wp)             :: phi_b             ! right bound
!     REAL(wp)             :: phidash           ! latitude in loop
!     REAL(wp)             :: dphi          ! step
!     REAL(wp)             :: u, temp   ! intermediate values
! 
!     !-----------------------------------------------------------------------
! 
!     phi_a = -0.5_wp * pi
!     phi_b = phi
!     
!     dphi = 0.01_wp * ( phi_b - phi_a)
!     
!     p_hh = 0._wp
!     
!     phidash = phi_a - 0.5_wp * dphi
!     
!     DO j = 1, 100
!       phidash = phidash + dphi
!       
!       u = func(phidash)
!       
!       temp = 2._wp * grid_angular_velocity * SIN(phidash)
!       temp = temp + ( u * TAN(phidash)* sphere_radius)
!       temp = sphere_radius *rgrav * u * temp
!       
!       p_hh = p_hh + temp * dphi
!       
!     ENDDO
!     
!     p_hh = 10000._wp - p_hh
!     !     print*, "phh", INT(360*phi/pi), INT(p_hh)
!     
!   END FUNCTION geostrophic_balance_11
!   !-----------------------------------------------------------------------------------
!   
!   !-----------------------------------------------------------------------------------
! !   SUBROUTINE fill_tracer_x_height(patch_3d, ocean_state)
! !     TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
! !     TYPE(t_hydro_ocean_state), TARGET :: ocean_state
! !     
! !     TYPE(t_subset_range), POINTER :: all_cells
! !     INTEGER :: tracer_idx, level, block, idx, start_cell_idx, end_cell_idx
! !     
! !     IF (.NOT. use_tracer_x_height) RETURN
! !     
! !     all_cells => patch_3d%p_patch_2d(1)%cells%ALL
! !     
! !     DO tracer_idx = 1, no_tracer
! !       ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration_x_height(:, :, :) = 0.0_wp
! !       DO block = all_cells%start_block, all_cells%end_block
! !         CALL get_index_range(all_cells, block, start_cell_idx, end_cell_idx)
! !         DO idx = start_cell_idx, end_cell_idx
! !           DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
! !             
! !             ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration_x_height(idx, level, block) = &
! !               & ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration(idx,level,block)   *      &
! !               & patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(idx, level, block)
! !             
! !           ENDDO
! !         ENDDO
! !       ENDDO
! !     ENDDO
! !     
! !   END SUBROUTINE fill_tracer_x_height
!   !-----------------------------------------------------------------------------------
!   
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  Williamson_test6_h
!  
! !FUNCTION INTERFACE: 
  FUNCTION Williamson_test6_h(p_lon, p_lat, p_t) RESULT(p_hh)
!
! !DESCRIPTION:
!
! Initial datum for geopotential h, test case 6 of Williamson et al.(1992).
!

! !DEFINED PARAMETERS:  
!    REAL (wp), PARAMETER  :: h0 = 8000._wp, re_omg_kk = 50._wp
    REAL (wp), PARAMETER  :: h0 = 8000._wp, omg_kk = 7.848e-6_wp !(re * omg_kk is not 50.)
                                                                 ! pripodas 
    INTEGER,   PARAMETER  :: r = 4

! !INPUT PARAMETERS:  
    REAL(wp) , INTENT(in) :: p_lon, p_lat, p_t
 
! !RETURN VALUE:  
    REAL(wp)              :: p_hh

! !LOCAL VARIABLES:  
   ! REAL(wp)              :: z_omg, z_phia, z_phib, z_phic, z_r_earth_angular_velocity 
    REAL(wp)              :: z_phia, z_phib, z_phic, z_r_omega , z_re_omg_kk 
    REAL(wp)              :: z_cosfi, z_cosfi2, z_cosfir, z_cosfir2, z_cosfir2m2
    REAL(wp)              :: z_cosdl, z_cosd2l, z_dlon, z_rr1r2
    REAL(wp)              :: z_val

    !INTEGER               :: i_r1, i_r1r1, i_r1r2, i_r2, j
    INTEGER               :: j
    REAL(wp)              :: z_r, z_r1, z_r1r1, z_r1r2, z_r2   !pripodas, better transform to real values

!EOP  
!-----------------------------------------------------------------------  
!BOC

    z_r_omega  = sphere_radius * earth_angular_velocity
    !z_omg     = re_omg_kk / re
    z_re_omg_kk= sphere_radius * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
   
    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r1    = z_r1 * z_r1
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2

    z_dlon    = omg_kk * z_r * (3._wp+z_r) - 2.0_wp * earth_angular_velocity
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (p_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)  
    z_cosd2l  = COS(2._wp * z_dlon)  

    z_cosfi   = COS(p_lat)
    z_cosfi2  = z_cosfi  * z_cosfi    ! cos^2(lat)

    z_cosfir  = z_cosfi
    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    z_cosfir2m2 = z_cosfir 
    z_cosfir2m2 = z_cosfir2m2 * z_cosfir2m2   ! cos^{2*r1-2}(lat)

    z_cosfir  = z_cosfir * z_cosfi    ! cos^{r1}(lat)
    z_cosfir2 = z_cosfir * z_cosfir   ! cos^{2*r1}(lat)

    z_val  = -.25_wp + z_r
    z_val  = 2._wp * z_val * z_val - 2.125_wp   ! 2r^2 - r -2

    z_phia = z_val * z_cosfi2
    
    z_val  = 2._wp * r * z_r

    z_phia = z_phia - z_val 
    z_val  = z_cosfi2 * z_cosfi2 * z_r1
    z_phia = z_phia + z_val 

    z_phia = .25_wp * z_re_omg_kk * z_re_omg_kk * z_cosfir2m2 * z_phia 
    z_val  = .5_wp * z_re_omg_kk * (2._wp * z_r_omega + z_re_omg_kk) * z_cosfi2
    z_phia = z_val + z_phia 
    
    z_phib = -1._wp * z_cosfi2 * z_r1r1 + z_r1r1 + 1._wp 
    z_phib = z_re_omg_kk * (z_r_omega + z_re_omg_kk) * z_cosfir * z_phib
    z_phib = 2._wp * z_rr1r2 * z_phib
 
    z_phic = z_r1 * z_cosfi2 - 1._wp * z_r2
    z_phic = .25_wp * z_re_omg_kk * z_re_omg_kk * z_cosfir2 * z_phic

    p_hh   = (z_phia + z_phib * z_cosdl + z_phic * z_cosd2l) * rgrav
    p_hh   = h0 + p_hh

  END FUNCTION Williamson_test6_h

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  Williamson_test6_u
!  
! !FUNCTION INTERFACE: 
  FUNCTION Williamson_test6_u( p_lon, p_lat, p_t) RESULT(p_uu)
!
! !DESCRIPTION:
!
! Initial datum for zonal velocity u, test case 6 of Williamson et al.(1992) .
!

! !DEFINED PARAMETERS:  
   ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER  ::  omg_kk = 7.848e-6_wp !(re * omg_kk is not 50.)
                                                                 ! pripodas 
    INTEGER,   PARAMETER  :: r = 4

! !INPUT PARAMETERS:  
    REAL(wp) , INTENT(in) :: p_lon, p_lat, p_t
 
! !RETURN VALUE:  
    REAL(wp)              :: p_uu

! !LOCAL VARIABLES:  
    !REAL(wp)              :: z_omg, z_r_omega 
    REAL(wp)              :: z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_cosfi2, z_sinfi, z_sinfi2
    REAL(wp)              :: z_cosfir, z_cosfirm1
    REAL(wp)              :: z_cosdl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val

!    INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER               :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2    !pripodas, better transform to real values

!EOP  
!-----------------------------------------------------------------------  
!BOC

    !z_omg     = re_omg_kk / re
    z_re_omg_kk= sphere_radius * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
    z_r       = REAL(r, wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2

    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * earth_angular_velocity
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (p_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)  

    z_sinfi   = SIN(p_lat)
    z_sinfi2  = z_sinfi * z_sinfi

    z_cosfi   = COS(p_lat)
    z_cosfi2  = z_cosfi * z_cosfi

    z_cosfir  = z_cosfi
    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    z_cosfirm1 = z_cosfir 

    z_val      = z_r * z_sinfi2 - z_cosfi2
    z_val      = z_cosfirm1 * z_val * z_cosdl
    z_val      = z_cosfi + z_val
    p_uu       = z_re_omg_kk * z_val 

  END FUNCTION Williamson_test6_u

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  Williamson_test6_v
!  
! !FUNCTION INTERFACE: 
  FUNCTION Williamson_test6_v( p_lon, p_lat, p_t) RESULT(p_vv)
!
! !DESCRIPTION:
!
! Initial datum for meridional velocity v, test case 6 of Williamson 
! et al.(1992).
!

! !DEFINED PARAMETERS:  
   ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER  ::  omg_kk = 7.848e-6_wp !(re * omg_kk is not 50.)
                                                                 ! pripodas 
    INTEGER,   PARAMETER  :: r = 4

! !INPUT PARAMETERS:  
    REAL(wp) , INTENT(in) :: p_lon, p_lat, p_t
 
! !RETURN VALUE:  
    REAL(wp)              :: p_vv

! !LOCAL VARIABLES:  
    !REAL(wp)              :: z_omg, z_r_omega   !pripodas, we use omg_kk and not re_omg_kk
    REAL(wp)              :: z_re_omg_kk
    REAL(wp)              :: z_cosfi, z_sinfi
    REAL(wp)              :: z_cosfir, z_cosfirm1
    REAL(wp)              :: z_sindl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val

!    INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER               :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2    !pripodas, better transform to real values

!EOP  
!-----------------------------------------------------------------------  
!BOC

    !z_omg     = re_omg_kk / re
    z_re_omg_kk= sphere_radius * omg_kk  !pripodas, the initial parameter is omg_kk and not re_omg_kk
    
    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2 

    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * earth_angular_velocity
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (p_lon - z_dlon) * z_r
    z_sindl   = SIN(z_dlon)  

    z_sinfi   = SIN(p_lat)

    z_cosfi   = COS(p_lat)

    z_cosfir  = z_cosfi
    DO j= 2, r-1
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO
    z_cosfirm1 = z_cosfir 

    z_val      = z_cosfirm1 * z_sinfi * z_sindl
    p_vv       = -1._wp * z_re_omg_kk * z_r * z_val 

  END FUNCTION Williamson_test6_v

!EOC  
!-------------------------------------------------------------------------  
!BOP
!
! !IROUTINE:  Williamson_test6_vort
!  
! !FUNCTION INTERFACE: 
   FUNCTION Williamson_test6_vort( p_lon, p_lat, p_t) RESULT(p_vt)
!
! !DESCRIPTION:
!
! Initial datum for relative vorticity, test case 6 of Williamson et al.(1992).
!

! !DEFINED PARAMETERS:  
   ! REAL (wp), PARAMETER  :: re_omg_kk = 50._wp
    REAL (wp), PARAMETER  ::  omg_kk = 7.848e-6_wp !(re * omg_kk is not 50.)
                                                                 ! pripodas 
    INTEGER,   PARAMETER  :: r = 4

! !INPUT PARAMETERS:  
    REAL(wp) , INTENT(in) :: p_lon, p_lat, p_t
 
! !RETURN VALUE:  
    REAL(wp)              :: p_vt

! !LOCAL VARIABLES:  
    !REAL(wp)              :: z_omg, z_r_omega   !pripodas, we use omg_kk and not re_omg_kk
    REAL(wp)              :: z_cosfi, z_sinfi
    REAL(wp)              :: z_cosfir
    REAL(wp)              :: z_cosdl, z_dlon, z_rr1r2
    REAL(wp)              :: z_val

    !INTEGER               :: i_r1, i_r1r2, i_r2, j
    INTEGER               :: j
    REAL(wp)              :: z_r, z_r1, z_r1r2, z_r2   !pripodas, better transform to real values

!EOP  
!-----------------------------------------------------------------------  
!BOC

    !z_omg     = re_omg_kk / re
    
    z_r       = REAL(r,wp)
    z_r1      = z_r + 1._wp
    z_r2      = z_r + 2._wp
    z_r1r2    = z_r1 * z_r2
    z_rr1r2   = 1._wp / z_r1r2

    z_dlon    = z_r * (3._wp+z_r) * omg_kk - 2.0_wp * earth_angular_velocity
    z_dlon    = z_dlon * z_rr1r2 * p_t
    z_dlon    = (p_lon - z_dlon) * z_r
    z_cosdl   = COS(z_dlon)  

    z_sinfi   = SIN(p_lat)

    z_cosfi   = COS(p_lat)

    z_cosfir  = z_cosfi
    DO j= 2, r
      z_cosfir = z_cosfir * z_cosfi   ! cos^{j}(lat)
    ENDDO

    z_val     = z_cosfir * z_r1 * z_r2 * z_cosdl
    z_val     = 2._wp - z_val
    p_vt      = omg_kk * z_sinfi * z_val

  END FUNCTION Williamson_test6_vort

!EOC  
  FUNCTION sphere_h( p_lon, p_lat) RESULT(h_site)

    !This test case situates a cone on the sphere, which declines with time.
    
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
    REAL(wp), INTENT(in) :: p_lat     ! latitude of point

    ! !RETURN VALUE:  
    REAL(wp) :: h_site, radius  

    REAL(wp) :: x2, y2, r, d


    d = 20000.0           !Height of the cone
    radius = 0.5_wp        !Radius of the cone
!   radius = 1.5*L_x        !Radius of the cone on the plane

    h_site=2000.0_wp !2000.0       !Depth of the ocean
    
    !  The cone can be situated at (0,0):
    x2=(p_lon)*(p_lon)
    y2=(p_lat)*(p_lat)

    r=sqrt(x2+y2)          !Distance of the specific point to the center of the cone.

    if(r.lt.(radius))then
       h_site=h_site-d/(radius*radius*radius)*(r*r*r)+3.0*d/(radius*radius)*(r*r) &
            & -3.0*d/(radius)*(r)+d
    end if

   
  END FUNCTION sphere_h


  FUNCTION sphere_u() RESULT(u_site)
! !RETURN VALUE:  
    REAL(wp)             :: u_site      ! meridional velocity

    u_site = 0.0_wp

  END FUNCTION sphere_u



  FUNCTION sphere_v() RESULT(v_site)
! !RETURN VALUE:  
    REAL(wp)             :: v_site      ! meridional velocity


    v_site = 0.0_wp

  END FUNCTION sphere_v


  FUNCTION sphere_oro() RESULT(p_or)

! !RETURN VALUE:  
    REAL(wp)             :: p_or      ! orography

!EOP  
!-----------------------------------------------------------------------  
!BOC

    p_or = 0.0


  END FUNCTION sphere_oro


  FUNCTION sphere_wind(p_lat, direction) RESULT(wd)
    REAL(wp), INTENT(IN) :: p_lat
    INTEGER, INTENT(in) :: direction

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':sphere_wind'

! !RETURN VALUE:  
    REAL(wp)             :: wd        ! wind

!EOP  
!-----------------------------------------------------------------------  
!BOC

    IF(direction==1)THEN
       wd = -0.1_wp*cos(9.0_wp*p_lat/4.0_wp)
    ELSE IF(direction==2)THEN
       wd = 0.0_wp
    ELSE
      CALL finish(method_name, "Wrong wind direction")
    END IF

  END FUNCTION  sphere_wind




  FUNCTION vortex_h() RESULT(h_site)

    ! !RETURN VALUE:  
    REAL(wp) :: h_site 

!    REAL(wp) :: x2, y2, r, d, radius 



!!$    d = 2000.0           !Height of the cone
!!$    radius = 0.25*L_x        !Radius of the cone
!!$!   radius = 1.5*L_x        !Radius of the cone
!!$
!!$    h_site=10200.0_wp+1.0/4000.0*L_x       !Height of the normal water columns
!!$    
!!$    x2=(p_lon-0.5*L_x)*(p_lon-0.5*L_x)
!!$    y2=(p_lat-0.5*L_y)*(p_lat-0.5*L_y)
!!$    r=sqrt(x2+y2)          !Distance of the specific point to the center of the cone.
!!$    
!!$    if(r.lt.(radius))then
!!$       h_site=h_site-d/(radius*radius*radius)*(r*r*r)+3.0*d/(radius*radius)*(r*r) &
!!$            & -3.0*d/(radius)*(r)+d
!!$    end if


    h_site=5.0_wp 
   
  END FUNCTION vortex_h


  FUNCTION vortex_u() RESULT(u_site)

! !RETURN VALUE:  
    REAL(wp)             :: u_site      ! meridional velocity

!   REAL(wp) :: d,radius 
!   REAL(wp) :: L_x = 5000000.0
!   REAL(wp) :: L_y = 4330000.0


!   d = 150.0           !Height of the cone
!   radius = 0.25*L_x        !Radius of the cone
!   radius = 1.5*L_x        !Radius of the cone

!!$    u_site=0.0_wp
!!$    
!!$    x2=(p_lon-0.5*L_x)*(p_lon-0.5*L_x)
!!$    y2=(p_lat-0.5*L_y)*(p_lat-0.5*L_y)
!!$    r=sqrt(x2+y2)          !Distance of the specific point to the center of the cone.
!!$    
!!$    if(r.lt.(radius))then
!!$       u_site=-d/(radius*radius*radius)*(r*r*r)+3.0*d/(radius*radius)*(r*r) &
!!$            & -3.0*d/(radius)*(r)+d
!!$    end if



!!$    abs = sqrt((p_lon-0.2_wp*L_x)**2.0_wp+&
!!$           & (p_lat-0.5_wp*L_y)**2.0_wp)
!!$      
!!$    IF(abs.lt.0.2_wp*L_y)THEN
!!$       u_site = 1.0_wp
!!$    ELSE
    u_site = 0.0_wp!1.57_wp
!    END IF
    

  END FUNCTION vortex_u



  FUNCTION vortex_v() RESULT(v_site)
! !RETURN VALUE:  
    REAL(wp)             :: v_site      ! meridional velocity


    v_site = 0.0_wp
  END FUNCTION vortex_v



  FUNCTION vortex_vort() RESULT(vort_site)
! !RETURN VALUE:  
    REAL(wp)             :: vort_site      ! meridional velocity
    

    vort_site = 1.0_wp
  END FUNCTION vortex_vort


  FUNCTION vortex_wind(direction) RESULT(wd)
 
! !INPUT PARAMETERS:  
    INTEGER, INTENT(in) :: direction

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':vortex_wind'


! !RETURN VALUE:  
    REAL(wp)             :: wd        ! wind

!EOP  
!-----------------------------------------------------------------------  
!BOC

    IF(direction==1)THEN
       wd = 1.1_wp
    ELSE IF(direction==2)THEN
       wd = 0.0_wp
    ELSE
       CALL finish(method_name, 'Wrong wind direction')
    END IF

  END FUNCTION vortex_wind




  FUNCTION vortex_oro(p_lon) RESULT(p_or)

! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: p_lon     ! longitude of point
!   REAL(wp), INTENT(in) :: p_lat     ! latitude of point

! !RETURN VALUE:  
    REAL(wp)             :: p_or      ! orography

    REAL(wp) :: x2, d
    REAL(wp) :: L_x = 5000000.0
!   REAL(wp) :: L_y = 4330000.0

!EOP  
!-----------------------------------------------------------------------  
!BOC



!!$    d = 0.3           !Height of the cone
!!$  
!!$    x2=(p_lon-0.5*L_x)
!!$    y2=(p_lat-0.5*L_y)
!!$    
!!$    p_or=d/(1.0_wp+x2*x2+y2*y2)**(3.0_wp/2.0_wp)
!!$
!!$
    d = 0.3           !Height of the cone
  
    x2=(p_lon-0.5*L_x)
!    y2=(p_lat-0.5*L_y)
    
    p_or=d/(1.0_wp+x2*x2)**(3.0_wp/2.0_wp)

!!$    IF(p_lon.ge.0.3_wp*L_x.and.p_lon.lt.0.4_wp*L_x)THEN
!!$       p_or=0.5_wp*(p_lon-0.3_wp*L_x)/(0.1_wp*L_x)
!!$    ELSE IF(p_lon.ge.0.4_wp*L_x.and.p_lon.lt.0.6_wp*L_x)THEN
!!$       p_or=0.5_wp
!!$    ELSE IF(p_lon.ge.0.6_wp*L_x.and.p_lon.lt.0.7_wp*L_x)THEN
!!$       p_or=0.5_wp*(1.0_wp-(p_lon-0.6_wp*L_x)/(0.1_wp*L_x))
!!$    END IF
    
    p_or=0.0_wp

  END FUNCTION vortex_oro



  FUNCTION vortex_sphere_h() RESULT(h_site)

    ! !RETURN VALUE:  
    REAL(wp) :: h_site 

    h_site=10000.0_wp 
   
  END FUNCTION vortex_sphere_h


  FUNCTION vortex_sphere_u() RESULT(u_site)

! !RETURN VALUE:  
    REAL(wp)             :: u_site      ! meridional velocity
 

    u_site = 0.0_wp    

  END FUNCTION vortex_sphere_u



  FUNCTION vortex_sphere_v() RESULT(v_site)

! !RETURN VALUE:  
    REAL(wp)             :: v_site      ! meridional velocity

    v_site = 0.0_wp
  END FUNCTION vortex_sphere_v

  FUNCTION vortex_wind_sphere() RESULT(wd)

! !RETURN VALUE:  
    REAL(wp)             :: wd        ! wind

    wd = 0.0_wp

  END FUNCTION vortex_wind_sphere



  FUNCTION galewsky_h( p_lat) RESULT(h_site)

   IMPLICIT NONE
   REAL(wp) , INTENT(in):: p_lat
   REAL(wp) :: h_site

   h_site = 10166.0_wp-geostr_balance11(p_lat, galewsky_uu)
   
  END FUNCTION galewsky_h


  FUNCTION galewsky_u(p_lat) RESULT(u_site)

   IMPLICIT NONE
   REAL(wp) , INTENT(in):: p_lat
   REAL(wp) ::  u_site
   REAL(wp) ::  phi0, phi1, umax, en

   phi0=pi/7._wp
   phi1=pi/2._wp - phi0
   en=exp(-4._wp/(phi0-phi1)**2)
   umax=80._wp

   if ((p_lat.gt.phi0).and.(p_lat.lt.phi1))then
        u_site=umax/en*exp(1._wp/(p_lat-phi0)/(p_lat-phi1))
        if (u_site.lt. 0.001) then
              u_site  = 0.0d0
        end if
   else 
        u_site=0._wp
   endif

  END FUNCTION galewsky_u

  FUNCTION galewsky_uu(p_lat) RESULT(u_site)

   IMPLICIT NONE
   REAL(wp) , INTENT(in):: p_lat
   REAL(wp) ::  u_site
   REAL(wp) ::  phi0, phi1, umax, en

   phi0=pi/7._wp
   phi1=pi/2._wp - phi0
   en=exp(-4._wp/(phi0-phi1)**2)
   umax=80._wp

   if ((p_lat.gt.phi0).and.(p_lat.lt.phi1))then
        u_site=umax/en*exp(1._wp/(p_lat-phi0)/(p_lat-phi1))
        if (u_site.lt. 0.001) then
              u_site  = 0.0d0
        end if
   else 
        u_site=0._wp
   endif

  END FUNCTION galewsky_uu


  FUNCTION galewsky_v() RESULT(v_site)

    REAL(wp)             :: v_site      ! meridional velocity


    v_site = 0.0_wp
    

  END FUNCTION galewsky_v

  FUNCTION galewsky_oro() RESULT(p_or)

    REAL(wp)             :: p_or      ! orography

!EOP  
!-----------------------------------------------------------------------  
!BOC

!    p_or = 10166.0_wp
  p_or = 0.0_wp

  END FUNCTION galewsky_oro

  FUNCTION geostr_balance11( phi, func)  RESULT(p_hh)
!
    ! !DESCRIPTION:
! Performs  numerical integration between -$\frac{\pi}{2}$ and $\frac{\pi}{2}$
! to compute geostrophically balanced initial state used
! in test 3.

! !INTERFACE:
    INTERFACE                        ! selected function

      FUNCTION func(p_lat) RESULT(p_vv)  

        USE mo_kind, ONLY: wp
               
        REAL(wp), INTENT(in) :: p_lat
        REAL(wp)             :: p_vv

      END FUNCTION func
       
    END INTERFACE
   
! !INPUT PARAMETERS:  
    REAL(wp), INTENT(in) :: phi           ! rotated latitude

! !RETURN VALUE:  
    REAL(wp)             :: p_hh            ! balanced height

! !LOCAL VARIABLES:  
    INTEGER              :: j,nint               ! loop index

    REAL(wp)             :: phi_a             ! left bound
    REAL(wp)             :: phi_b             ! right bound
    REAL(wp)             :: phidash           ! latitude in loop
    REAL(wp)             :: dphi          ! step
    REAL(wp)             :: temp   ! intermediate values


!EOP  
!-----------------------------------------------------------------------  
!BOC

!!$    !Midpoint rule:
!!$
!!$    phi_a = -0.5_wp * pi
!!$    phi_b = phi
!!$
!!$    nint = 1000
!!$
!!$    dphi = 1.0_wp/real(nint) * ( phi_b - phi_a)
!!$
!!$    p_hh = 0._wp
!!$
!!$    phidash = phi_a - 0.5_wp * dphi
!!$
!!$    DO j = 1, nint
!!$       phidash = phidash + dphi
!!$       
!!$       u = func(phidash)
!!$      
!!$       temp = 2._wp * omega * SIN(phidash)
!!$       temp = temp + ( u * TAN(phidash)* rre)
!!$       temp = re *rgrav * u * temp
!!$
!!$       p_hh = p_hh + temp * dphi
!!$
!!$    ENDDO



    !Composite Simpson rule (following wikipedia)
    phi_a = -0.5_wp * pi
    phi_b = phi

    nint = 10000    !2*nint is the number of intervals

    dphi = 1.0_wp/(2.0_wp*real(nint)) * ( phi_b - phi_a)

    p_hh = 0._wp

    !First value:
    phidash = phi_a

    temp = h_value(phidash)
    p_hh = p_hh + temp * dphi/3.0_wp

    !and last value:

    phidash = phi_b

    temp = h_value(phidash)
    p_hh = p_hh + temp * dphi/3.0_wp

    !First sum:

    phidash = phi_a
    


    DO j = 1, (nint-1)
       phidash = phidash + 2.0_wp*dphi

       temp = h_value(phidash)
       p_hh = p_hh + 2.0_wp*temp * dphi/3.0_wp

    ENDDO

    !Last sum:
     phidash = phi_a - dphi
    
    DO j = 1, nint
       phidash = phidash + 2.0_wp*dphi
       
       temp = h_value(phidash)
       p_hh = p_hh + 4.0_wp*temp * dphi/3.0_wp


    ENDDO


    CONTAINS
    
      FUNCTION h_value(phi)  RESULT(temp)  
        REAL(wp) :: phi
        REAL(wp) :: temp
        REAL(wp) :: u

        u = func(phi)
       
        temp = 2._wp * earth_angular_velocity * SIN(phi)
        temp = temp + ( u * TAN(phi)* inverse_earth_radius)
        temp = sphere_radius *rgrav * u * temp
        
      END FUNCTION h_value

  END FUNCTION geostr_balance11

   SUBROUTINE tracer_bubble(patch_3d, tracer,bubble_inside, bubble_outside, lat_bubble_opt, lon_bubble_opt,&
  & radius_bubble_opt, layers_above_bubble_opt)
! This subroutine places an ellipsoid of defined salinity or temperature, which has its maximum/minimum 
! value at the midpoint and approaches the value of its environment linearly with radius
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d 
    REAL(wp), TARGET :: tracer(:,:,:)
    REAL(wp),intent(in):: bubble_inside, bubble_outside
    REAL(wp),intent(in),optional:: lat_bubble_opt, lon_bubble_opt, radius_bubble_opt
    INTEGER,intent(in),optional::layers_above_bubble_opt

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    INTEGER ::layers_above_bubble, layers_bubble
    REAL(wp):: lat_deg, lon_deg, z_tmp, test
    REAL(wp):: dist, dist_layer
    REAL(wp):: lat_bubble, lon_bubble, radius_bubble

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_bubble'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Coordinates midpoint of bubble in South Pacific:
    lat_bubble= 20.0_wp
    lon_bubble  = -40.0_wp
    radius_bubble = 2.5_wp  !5.0_wp
    layers_above_bubble = 78 !31
    layers_bubble = 2 !46 
    dist_layer=layers_bubble/2.0_wp !radius in z direction
    CALL assign_if_present(lat_bubble,lat_bubble_opt)
    CALL assign_if_present(lon_bubble,lon_bubble_opt)
    CALL assign_if_present(radius_bubble,radius_bubble_opt)
    CALL assign_if_present(layers_above_bubble,layers_above_bubble_opt)
    CALL assign_if_present(layers_bubble,layers_above_bubble_opt)


    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transfer to latitude in degrees
        lat_deg = cell_center(idx,block)% lat * rad2deg
        lon_deg = cell_center(idx,block)% lon * rad2deg

        !to determine the closest point on grid to given midpoint of the bubble
        dist = SQRT( (lat_bubble-lat_deg)**2.0_wp + (lon_bubble-lon_deg)**2.0_wp )

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

            IF (      (dist <= radius_bubble) &
              & .AND. (level <= layers_above_bubble + layers_bubble) &
              & .AND. (level >= layers_above_bubble)) THEN
              
              test = (ABS(level - (layers_above_bubble + dist_layer)) / dist_layer)
              z_tmp = (dist/radius_bubble) + test

              IF (z_tmp <= 1) THEN

                tracer(idx,level,block) = bubble_inside !&
                   !& bubble_outside + (bubble_inside-bubble_outside) * (1.0_wp - z_tmp) * (1.0_wp - test)

              ELSE

                tracer(idx,level,block) = bubble_outside

              END IF

            ELSE

                tracer(idx,level,block) = bubble_outside

            END IF

        END DO 
      END DO
    END DO


   END SUBROUTINE tracer_bubble


 SUBROUTINE tracer_bubble_disturbed(patch_3d, tracer,bubble_inside, bubble_outside, lat_bubble_opt,&
  & lon_bubble_opt, radius_bubble_opt, layers_above_bubble_opt)
! This subroutine places an ellipsoid of defined salinity or temperature, which has its maximum/minimum 
! value at the midpoint and approaches the value of its environment linearly with radius. There is a
! disturbance at the upper side of the bubble if it's lighter than the environment and vise versa.
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d 
    REAL(wp), TARGET :: tracer(:,:,:)
    REAL(wp),intent(in):: bubble_inside, bubble_outside
    REAL(wp),intent(in),optional:: lat_bubble_opt, lon_bubble_opt, radius_bubble_opt
    INTEGER,intent(in),optional::layers_above_bubble_opt

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    INTEGER ::layers_above_bubble, layers_bubble, layers_perturbation
    REAL(wp):: lat_deg, lon_deg, z_tmp, test, amplitude_perturbation
    REAL(wp):: dist, dist_layer
    REAL(wp):: lat_bubble, lon_bubble, radius_bubble

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_bubble_disturbed'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Coordinates midpoint of bubble (in South Pacific):
    lat_bubble= 20.0_wp
    lon_bubble  = -40.0_wp
    radius_bubble = 15.0_wp
    layers_above_bubble = 35 !15 !35
    layers_bubble = 40 !20 !40 
    dist_layer=layers_bubble/2.0_wp  !radius in z direction
    layers_perturbation = 1
    amplitude_perturbation = 0.10_wp
    CALL assign_if_present(lat_bubble,lat_bubble_opt)
    CALL assign_if_present(lon_bubble,lon_bubble_opt)
    CALL assign_if_present(radius_bubble,radius_bubble_opt)
    CALL assign_if_present(layers_above_bubble,layers_above_bubble_opt)
    CALL assign_if_present(layers_bubble,layers_above_bubble_opt)

    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transfer to latitude in degrees
        lat_deg = cell_center(idx,block)% lat * rad2deg
        lon_deg = cell_center(idx,block)% lon * rad2deg

        !to determine the closest point on grid to given midpoint of the bubble
        dist = SQRT( (lat_bubble-lat_deg)**2.0_wp + (lon_bubble-lon_deg)**2.0_wp )

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

            IF (      (dist <= radius_bubble) &
              & .AND. (level <= layers_above_bubble + layers_bubble) &
              & .AND. (level >= layers_above_bubble)) THEN
              
              test = (ABS(level - (layers_above_bubble + dist_layer)) / dist_layer)
              z_tmp = (dist/radius_bubble) + test

               IF (z_tmp <= 1) THEN
 
                 tracer(idx,level,block) = &
                    & bubble_outside + (bubble_inside-bubble_outside) * (1.0_wp - z_tmp) * (1.0_wp - test)
 
                      IF( bubble_inside < bubble_outside .AND. & 
                          & level >= layers_above_bubble + layers_bubble - layers_perturbation ) THEN
    
                         tracer(idx,level,block) = bubble_outside + (bubble_inside-bubble_outside) *&
                                            & ( 1.0_wp - z_tmp ) * ( 1.0_wp - test )- &
                                            & amplitude_perturbation + &
                                            & SIN( 2.0_wp * pi * dist / radius_bubble ) * & 
                                            !& COS( pi * dist / ( 2.0_wp * radius_bubble ) ) * & 
                                            & amplitude_perturbation
      
                      END IF
 
                     IF( bubble_inside > bubble_outside .AND. &
                         & level <= layers_above_bubble + layers_perturbation ) THEN
    
                        tracer(idx,level,block) = bubble_outside + (bubble_inside-bubble_outside) *&
                                           & ( 1.0_wp - z_tmp ) * ( 1.0_wp - test ) + &
                                           & amplitude_perturbation - &
                                           & SIN( 2.0_wp * pi * dist / radius_bubble ) * & 
                                           & amplitude_perturbation
     
                     END IF

              ELSE

                tracer(idx,level,block) = bubble_outside

              END IF

            ELSE

                tracer(idx,level,block) = bubble_outside

            END IF

        END DO 
      END DO
    END DO


   END SUBROUTINE tracer_bubble_disturbed

 SUBROUTINE tracer_double_bubble(patch_3d, tracer,bubble_inside, bubble_outside, lat_bubble_opt,&
  & lon_bubble_opt, radius_bubble_opt, layers_above_bubble_opt)
! This subroutine places an ellipsoid of defined salinity or temperature, which has its minimum/maximum 
! value at the midpoint and approaches the value of its environment linearly with radius. Below the
! bubble is another smaller bubble with opposite density.
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d 
    REAL(wp), TARGET :: tracer(:,:,:)
    REAL(wp),intent(in):: bubble_inside, bubble_outside
    REAL(wp),intent(in),optional:: lat_bubble_opt, lon_bubble_opt, radius_bubble_opt
    INTEGER,intent(in),optional::layers_above_bubble_opt

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    INTEGER ::layers_above_bubble, layers_bubble
    REAL(wp):: lat_deg, lon_deg, z_tmp, test
    REAL(wp):: dist, dist_layer, lat_small_bubble, dist_small
    REAL(wp):: lat_bubble, lon_bubble, radius_bubble, small_bubble_inside

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_double_bubble'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Coordinates midpoint of bubble in South Pacific:
    lat_bubble= 20.0_wp
    lon_bubble  = -40.0_wp
    radius_bubble = 15.0_wp
    layers_above_bubble = 5
    layers_bubble = 35 
    lat_small_bubble= 20.0_wp

    dist_layer=layers_bubble/2.0_wp  !radius in z direction
    small_bubble_inside = 2.0_wp * ( bubble_outside - bubble_inside ) + bubble_inside
    CALL assign_if_present(lat_bubble,lat_bubble_opt)
    CALL assign_if_present(lon_bubble,lon_bubble_opt)
    CALL assign_if_present(radius_bubble,radius_bubble_opt)
    CALL assign_if_present(layers_above_bubble,layers_above_bubble_opt)
    CALL assign_if_present(layers_bubble,layers_above_bubble_opt)


    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transfer to latitude in degrees
        lat_deg = cell_center(idx,block)% lat * rad2deg
        lon_deg = cell_center(idx,block)% lon * rad2deg

        !to determine the closest point on grid to given midpoint of the bubble
        dist = SQRT( (lat_bubble-lat_deg)**2.0_wp + (lon_bubble-lon_deg)**2.0_wp )
        dist_small = SQRT( (lat_small_bubble-lat_deg)**2.0_wp + (lon_bubble-lon_deg)**2.0_wp )

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

            IF (      (dist <= radius_bubble) &
              & .AND. (level <= layers_above_bubble + layers_bubble) &
              & .AND. (level >= layers_above_bubble)) THEN
              
              test = (ABS(level - (layers_above_bubble + dist_layer)) / dist_layer)
              z_tmp = (dist/radius_bubble) + test

              IF (z_tmp <= 1) THEN

                tracer(idx,level,block) = &
                   & bubble_outside + (bubble_inside-bubble_outside) * (1.0_wp - z_tmp) * (1.0_wp - test)

              ELSE

                tracer(idx,level,block) = bubble_outside

              END IF

            ! small bubble beneath the bigger one 
            ELSEIF (  ( dist_small <= radius_bubble / 3.0_wp ) &
              & .AND. ( level <= layers_above_bubble + layers_bubble + 5 +  layers_bubble - 10 ) &
              & .AND. ( level >= layers_above_bubble + layers_bubble + 5 ) ) THEN

              test = (ABS(level - (layers_above_bubble  + layers_bubble + 5 + dist_layer - 5)) / dist_layer)
              z_tmp = ( 4.0_wp * dist_small / radius_bubble ) + test

              IF (z_tmp <= 1) THEN

                tracer(idx,level,block) = &
                   & bubble_outside + ( small_bubble_inside - bubble_outside ) * &
                   & ( 1.0_wp - z_tmp ) * ( 1.0_wp - test )

              ELSE

                tracer(idx,level,block) = bubble_outside

              END IF

            ELSE

                tracer(idx,level,block) = bubble_outside

            END IF

        END DO 
      END DO
    END DO


   END SUBROUTINE tracer_double_bubble


  SUBROUTINE tracer_layer(patch_3d, tracer,tracer_top_opt, tracer_bottom_opt)
! This subroutine places an layer of defined salinity or temperature above an 
! layer of another salinity or temperature. At a given latitude and longitude
! there is a disturbance.
     TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d 
     REAL(wp), TARGET :: tracer(:,:,:)
     REAL(wp),intent(in):: tracer_top_opt, tracer_bottom_opt
 
     TYPE(t_patch),POINTER   :: patch_2d
     TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
     TYPE(t_subset_range), POINTER :: all_cells
 
     INTEGER :: block, idx, level
     INTEGER :: start_cell_index, end_cell_index
     REAL(wp) :: lat_deg, lon_deg, amplitude_perturbation, dist
     REAL(wp) :: lat_disturbance, lon_disturbance, z_tmp, layers_top
     REAL(wp) :: tracer_top, tracer_bottom, test, radius_bubble
 
     CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_layer'
     !-------------------------------------------------------------------------
 
 
     patch_2d => patch_3d%p_patch_2d(1)
     all_cells => patch_2d%cells%ALL
     cell_center => patch_2d%cells%center

     amplitude_perturbation = 3.0
     lat_disturbance = 20.0_wp
     lon_disturbance = -40.0_wp
     layers_top = 10
 
     CALL assign_if_present(tracer_top,tracer_top_opt)
     CALL assign_if_present(tracer_bottom,tracer_bottom_opt)
 
     DO block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
       DO idx = start_cell_index, end_cell_index
 
         !transfer to latitude in degrees
         lat_deg = cell_center(idx,block)% lat * rad2deg
         lon_deg = cell_center(idx,block)% lon * rad2deg
         dist = SQRT( (lat_disturbance-lat_deg)**2.0_wp + (lon_disturbance-lon_deg)**2.0_wp )
         layers_top = 10.0 + 2 * ( SIN( ( 2.0 * pi * lon_deg / 180) ) + SIN( (2.0 * pi * lat_deg) / 360 ) )

         DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
 
            IF ( level <= layers_top ) THEN

                tracer(idx,level,block) = tracer_top
 
            ELSE
 
                tracer(idx,level,block) = tracer_bottom
 
            END IF


            IF (      (dist <= radius_bubble) &
              & .AND. (level <= layers_top + amplitude_perturbation) &
              & .AND. (level >= layers_top)) THEN
              
              test = (ABS(level - (layers_top + amplitude_perturbation)) / amplitude_perturbation)
              z_tmp = (dist/radius_bubble) + test

              IF (z_tmp <= 1) THEN

                tracer(idx,level,block) = tracer_top

              END IF
! 
! !            IF ( ( level == layers_top + 1 ) .AND. dist <= 10 ) THEN
!              IF ( dist <= 10 .AND. level <= layers_top + COS( pi * (level -layers_top) / 2.0_wp ) &
!                 & * amplitude_perturbation) THEN
!                 
!                    tracer(idx,level,block) = tracer_top !+ amplitude_perturbation - &
! !                                    & COS( pi * dist / 20.0_wp ) * &
! !                                    & amplitude_perturbation 
 
            END IF 

         END DO 
       END DO
     END DO
 
 
    END SUBROUTINE tracer_layer



 SUBROUTINE tracer_bubbles_side_by_side(patch_3d, tracer,bubble_inside, bubble_outside, lat_bubble_opt, lon_bubble_opt,&
  & radius_bubble_opt, layers_above_bubble_opt)
! This subroutine places two ellipsoid of defined salinity or temperature next to each other, which has its maximum/minimum 
! value at the midpoint and approaches the value of its environment linearly with radius
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d 
    REAL(wp), TARGET :: tracer(:,:,:)
    REAL(wp),intent(in):: bubble_inside, bubble_outside
    REAL(wp),intent(in),optional:: lat_bubble_opt, lon_bubble_opt, radius_bubble_opt
    INTEGER,intent(in),optional::layers_above_bubble_opt

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    INTEGER ::layers_above_bubble, layers_bubble
    REAL(wp):: lat_deg, lon_deg, z_tmp, test
    REAL(wp):: dist, dist_layer,dist2, bubble_inside2
    REAL(wp):: lat_bubble, lon_bubble, radius_bubble

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tracer_bubbles_side_by_side'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Coordinates midpoint of bubble in South Pacific:
    lat_bubble= 20.0_wp
    lon_bubble  = -40.0_wp
    radius_bubble = 15.0_wp
    layers_above_bubble = 35
    layers_bubble = 40 
    dist_layer = layers_bubble/2.0_wp !radius in z direction
    bubble_inside2 = bubble_outside + bubble_inside
    CALL assign_if_present(lat_bubble,lat_bubble_opt)
    CALL assign_if_present(lon_bubble,lon_bubble_opt)
    CALL assign_if_present(radius_bubble,radius_bubble_opt)
    CALL assign_if_present(layers_above_bubble,layers_above_bubble_opt)
    CALL assign_if_present(layers_bubble,layers_above_bubble_opt)


    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transfer to latitude in degrees
        lat_deg = cell_center(idx,block)% lat * rad2deg
        lon_deg = cell_center(idx,block)% lon * rad2deg

        !to determine the closest point on grid to given midpoint of the bubble
        dist = SQRT( (lat_bubble-lat_deg)**2.0_wp + (lon_bubble-lon_deg)**2.0_wp )
        dist2 = SQRT( (lat_bubble+radius_bubble-lat_deg)**2.0_wp + (lon_bubble-lon_deg)**2.0_wp )

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

            IF (      (dist <= radius_bubble) &
              & .AND. (level <= layers_above_bubble + layers_bubble) &
              & .AND. (level >= layers_above_bubble)) THEN
              
              test = (ABS(level - (layers_above_bubble + dist_layer)) / dist_layer)
              z_tmp = (dist/radius_bubble) + test

              IF (z_tmp <= 1) THEN

                tracer(idx,level,block) = &
                   & bubble_outside + (bubble_inside-bubble_outside) * (1.0_wp - z_tmp) * (1.0_wp - test)

              ELSE

                tracer(idx,level,block) = bubble_outside

              END IF

            ELSEIF (  (dist2 <= radius_bubble) &
              & .AND. (level <= layers_above_bubble + layers_bubble) &
              & .AND. (level >= layers_above_bubble)) THEN
              
              test = (ABS(level - (layers_above_bubble + dist_layer)) / dist_layer)
              z_tmp = (dist2/radius_bubble) + test

              IF (z_tmp <= 1) THEN

                tracer(idx,level,block) = &
                   & bubble_outside + (bubble_inside2-bubble_outside) * (1.0_wp - z_tmp) * (1.0_wp - test)

              ELSE

                tracer(idx,level,block) = bubble_outside

              END IF

            ELSE

                tracer(idx,level,block) = bubble_outside

            END IF

        END DO 
      END DO
    END DO

  END SUBROUTINE tracer_bubbles_side_by_side


  SUBROUTINE Roberts_tracer_bubble(patch_3d, tracer,bubble_inside, bubble_outside, lat_bubble_opt, lon_bubble_opt,&
    & radius_bubble_opt, layers_above_bubble_opt)
! The difference to the upper bubble cases is that the radius is given in meter and  the temperature decreases exponential
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d 
    REAL(wp), TARGET :: tracer(:,:,:)
    REAL(wp),intent(in):: bubble_inside, bubble_outside
    REAL(wp),intent(in),optional:: lat_bubble_opt, lon_bubble_opt, radius_bubble_opt
    INTEGER,intent(in),optional::layers_above_bubble_opt

    TYPE(t_patch),POINTER   :: patch_2d
    TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: block, idx, level
    INTEGER :: start_cell_index, end_cell_index
    INTEGER ::layers_above_bubble, layers_bubble
    REAL(wp):: lat_deg, lon_deg
    REAL(wp):: dist_xy, dist_z, dist_layer
    REAL(wp):: lat_bubble, lon_bubble, radius_bubble
    REAL(wp):: layerthickness, a, s, r, b,c, dLat, dLon

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':Roberts_tracer_bubble'
    !-------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    cell_center => patch_2d%cells%center

    ! Coordinates midpoint of bubble in South Pacific:
    layerthickness = 30000.0 !m
    lat_bubble= 20.0_wp
    lon_bubble  = -40.0_wp
    radius_bubble = 15.0_wp
    layers_above_bubble =47  !45
    layers_bubble = 20 
    dist_layer=layers_bubble/2.0_wp !radius in z directions
    a = 250000.0_wp  ! inner radius in meter
    s = 250000.0_wp  ! outer radius in meter
  

    CALL assign_if_present(lat_bubble,lat_bubble_opt)
    CALL assign_if_present(lon_bubble,lon_bubble_opt)
    CALL assign_if_present(radius_bubble,radius_bubble_opt)
    CALL assign_if_present(layers_above_bubble,layers_above_bubble_opt)
    CALL assign_if_present(layers_bubble,layers_above_bubble_opt)


    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
      DO idx = start_cell_index, end_cell_index

        !transfer to latitude in degrees
        lat_deg = cell_center(idx,block)% lat * rad2deg
        lon_deg = cell_center(idx,block)% lon * rad2deg


        ! transform latitude to distance to midpoint in meters
        dLat = (lat_bubble-lat_deg) * deg2rad
        dLon = (lon_bubble-lon_deg) * deg2rad
        b = sin( dLat / 2.0_wp ) * sin( dLat / 2.0_wp ) &
          & + cos( lat_deg * deg2rad ) *  cos( lat_bubble * deg2rad ) &
          & * sin( dLon / 2.0_wp ) * sin( dLon / 2.0_wp )
        c = 2.0_wp * atan2( sqrt( b ) , sqrt( 1.0_wp - b ) )
        dist_xy = 6378.137_wp * c * 1000.0_wp !1609.00_wp

        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

             ! transform number of layers to distance to midpoint in meters
             dist_z = abs(level - (layers_above_bubble + dist_layer)) &
                    & * layerthickness 

              r = sqrt( dist_xy**2.0_wp + dist_z**2.0_wp )


              IF ( r <= a ) THEN

                tracer(idx,level,block) = bubble_inside

              ELSEIF ( r > a ) THEN

                tracer(idx,level,block) = &
                   & bubble_outside + (bubble_inside-bubble_outside) &
                   & * exp( - ( r - a )**2.0_wp / s**2.0_wp )

              END IF

        END DO 
      END DO
    END DO

  END SUBROUTINE Roberts_tracer_bubble


  SUBROUTINE inclined_layer(patch_3d, tracer,tracer_top_opt, tracer_bottom_opt)

     TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d 
     REAL(wp), TARGET :: tracer(:,:,:)
     REAL(wp),intent(in):: tracer_top_opt, tracer_bottom_opt
 
     TYPE(t_patch),POINTER   :: patch_2d
     TYPE(t_geographical_coordinates), POINTER :: cell_center(:,:)
     TYPE(t_subset_range), POINTER :: all_cells
 
     INTEGER :: block, idx, level
     INTEGER :: start_cell_index, end_cell_index
     REAL(wp) :: lat_deg
     REAL(wp) :: h, depth, lat_neu
     REAL(wp) :: tracer_top, tracer_bottom
 
     CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':inclined_layer'
     !-------------------------------------------------------------------------
 
     depth = 80.0_wp
 
     patch_2d => patch_3d%p_patch_2d(1)
     all_cells => patch_2d%cells%ALL
     cell_center => patch_2d%cells%center

     CALL assign_if_present(tracer_top,tracer_top_opt)
     CALL assign_if_present(tracer_bottom,tracer_bottom_opt)
 
     DO block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, block, start_cell_index, end_cell_index)
       DO idx = start_cell_index, end_cell_index
 
         !transfer to latitude in degrees
         lat_deg = cell_center(idx,block)% lat * rad2deg

         IF ( abs( lat_deg ) < 65.0_wp ) THEN
  
          lat_neu = lat_deg + 65.0_wp
          h = lat_neu / 130.0_wp * depth

          DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)

             IF ( level > h  ) THEN

                 tracer(idx,level,block) = tracer_bottom
 
             ELSE
 
                 tracer(idx,level,block) = tracer_top
 
             END IF

         END DO 

         ELSEIF ( lat_deg > 65.0_wp ) THEN
          
          tracer(idx,level,block) = tracer_bottom
        
         ELSE

          tracer(idx,level,block) = tracer_top
         
         END IF

       END DO
     END DO
 
 
    END SUBROUTINE inclined_layer




  !-----------------------------------------------------------------------------------
!   SUBROUTINE fill_tracer_x_height(patch_3d, ocean_state)
!     TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
!     TYPE(t_hydro_ocean_state), TARGET :: ocean_state
!     
!     TYPE(t_subset_range), POINTER :: all_cells
!     INTEGER :: tracer_idx, level, block, idx, start_cell_idx, end_cell_idx
!     
!     IF (.NOT. use_tracer_x_height) RETURN
!     
!     all_cells => patch_3d%p_patch_2d(1)%cells%ALL
!     
!     DO tracer_idx = 1, no_tracer
!       ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration_x_height(:, :, :) = 0.0_wp
!       DO block = all_cells%start_block, all_cells%end_block
!         CALL get_index_range(all_cells, block, start_cell_idx, end_cell_idx)
!         DO idx = start_cell_idx, end_cell_idx
!           DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(idx,block)
!             
!             ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration_x_height(idx, level, block) = &
!               & ocean_state%p_prog(nold(1))%ocean_tracers(tracer_idx)%concentration(idx,level,block)   *      &
!               & patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(idx, level, block)
!             
!           ENDDO
!         ENDDO
!       ENDDO
!     ENDDO
!     
!   END SUBROUTINE fill_tracer_x_height
  !-----------------------------------------------------------------------------------
  
END MODULE mo_ocean_initial_conditions
