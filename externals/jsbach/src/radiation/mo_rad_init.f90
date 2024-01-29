!> Initialization of the the radiation memory.
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
MODULE mo_rad_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid
  USE mo_jsb_grid,          ONLY: Get_grid
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,            ONLY: missval

  dsl4jsb_Use_processes RAD_
  dsl4jsb_Use_config(RAD_)
  dsl4jsb_Use_memory(RAD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: rad_init

  TYPE t_rad_init_vars
    REAL(wp), POINTER ::  &
      & alb_background(:,:), &
      & alb_vis_soil  (:,:), &
      & alb_nir_soil  (:,:), &
      & alb_vis_can   (:,:), &
      & alb_nir_can   (:,:)
  END TYPE t_rad_init_vars

  TYPE(t_rad_init_vars) :: rad_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_init'

CONTAINS

  !
  !> Intialize radiation process (after memory has been set up)
  !
  SUBROUTINE rad_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':rad_init'

    model => get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL rad_read_init_vars(tile)
    END IF

    CALL rad_init_bc(tile)
    CALL rad_init_ic(tile)
    !! quincy
    IF (model%config%use_quincy .AND. tile%contains_vegetation) THEN
      CALL rad_init_q_ic(tile)
    ENDIF

    IF (tile%Is_last_process_tile(RAD_)) THEN
      CALL rad_finalize_init_vars()
    END IF

  END SUBROUTINE rad_init

  SUBROUTINE rad_finalize_init_vars

    DEALLOCATE( &
      & rad_init_vars%alb_background, &
      & rad_init_vars%alb_vis_soil,   &
      & rad_init_vars%alb_nir_soil,   &
      & rad_init_vars%alb_vis_can,    &
      & rad_init_vars%alb_nir_can     &
      & )
      
  END SUBROUTINE rad_finalize_init_vars

  SUBROUTINE rad_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    INTEGER :: nproma, nblks

    dsl4jsb_Def_config(RAD_)

    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER  :: ptr_2D(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':rad_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()


    ALLOCATE( &
      & rad_init_vars%alb_background(nproma, nblks), &
      & rad_init_vars%alb_vis_soil  (nproma, nblks), &
      & rad_init_vars%alb_nir_soil  (nproma, nblks), &
      & rad_init_vars%alb_vis_can   (nproma, nblks), &
      & rad_init_vars%alb_nir_can   (nproma, nblks)  &
      & )
      
    rad_init_vars%alb_background(:,:) = missval
    rad_init_vars%alb_vis_soil  (:,:) = missval
    rad_init_vars%alb_nir_soil  (:,:) = missval
    rad_init_vars%alb_vis_can   (:,:) = missval
    rad_init_vars%alb_nir_can   (:,:) = missval

    dsl4jsb_Get_config(RAD_)

    IF (.NOT. dsl4jsb_Config(RAD_)%active) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting radiation init vars')

    ! TODO: adapt names of variables in input files to names in memory

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(RAD_)%bc_filename), model%grid_id)

    ! background albedo (without snow or canopy)
    ptr_2D => input_file%Read_2d( &
      & variable_name='albedo',           &
      & fill_array = rad_init_vars%alb_background)
    ptr_2D = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    ! soil albedo
    IF (tile%contains_soil) THEN
      ptr_2D => input_file%Read_2d( &
        & variable_name='albedo_soil_vis',  &
        & fill_array = rad_init_vars%alb_vis_soil)
      ptr_2D = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

      ptr_2D => input_file%Read_2d( &
        & variable_name='albedo_soil_nir',  &
        & fill_array = rad_init_vars%alb_nir_soil)
      ptr_2D = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
    END IF

    IF (tile%contains_vegetation) THEN
      ! Get canopy albedo
      IF (dsl4jsb_Config(RAD_)%use_alb_canopy .OR. tile%lcts(1)%lib_id == 0) THEN
        ptr_2D => input_file%Read_2d( &
          & variable_name='albedo_veg_vis',   &
          & fill_array = rad_init_vars%alb_vis_can)
        ptr_2D = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
      ELSE
        rad_init_vars%alb_vis_can = dsl4jsb_Lctlib_param(AlbedoCanopyVIS)
      END IF
      IF (dsl4jsb_Config(RAD_)%use_alb_canopy .OR. tile%lcts(1)%lib_id == 0) THEN
        ptr_2D => input_file%Read_2d( &
          & variable_name='albedo_veg_nir',   &
          & fill_array = rad_init_vars%alb_nir_can)
        ptr_2D = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
      ELSE
        rad_init_vars%alb_nir_can = dsl4jsb_Lctlib_param(AlbedoCanopyNIR)
      END IF

    END IF

    CALL input_file%Close()

    NULLIFY(ptr_2D)

  END SUBROUTINE rad_read_init_vars

  SUBROUTINE rad_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! TYPE(t_jsb_model), POINTER :: model

    ! dsl4jsb_Def_memory(RAD_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':rad_init_bc'

    IF (.NOT. tile%Is_process_active(RAD_)) RETURN

    ! model => get_model(tile%owner_model_id)

    ! dsl4jsb_Get_memory(RAD_)

    ! IF (debug_on()) CALL message(routine, 'Setting  radiation boundary conditions for tile '//TRIM(tile%name))

    ! nothing to do in the moment ...

  END SUBROUTINE rad_init_bc

  SUBROUTINE rad_init_ic(tile)

    USE mo_rad_constants,  ONLY: AlbedoCanopySnow_age, AlbedoCanopySnow_temp, AlbedoVisInitial, AlbedoNirInitial, &
      &                          AlbedoGlacierVisMax, AlbedoGlacierNirMax, AlbedoLakeWater

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_memory(RAD_)

    REAL(wp) :: albedo_age_weight

    CHARACTER(len=*), PARAMETER :: routine = modname//':rad_init_ic'

    IF (.NOT. tile%Is_process_active(RAD_)) RETURN

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_config(RAD_)

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing radiation memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !
    ! TBD: adapt names of variables in input files to names in memory

    ! background albedo (without snow or canopy)
    dsl4jsb_var2D_onDomain(RAD_,alb_background) = rad_init_vars%alb_background

    IF (tile%is_lake) THEN
      dsl4jsb_var2D_onDomain(RAD_,alb_vis) = AlbedoLakeWater
      dsl4jsb_var2D_onDomain(RAD_,alb_nir) = AlbedoLakeWater
    ELSE IF (tile%is_glacier) THEN
      dsl4jsb_var2D_onDomain(RAD_,alb_vis)     = AlbedoGlacierVisMax
      dsl4jsb_var2D_onDomain(RAD_,alb_nir)     = AlbedoGlacierNirMax
      dsl4jsb_var2D_onDomain(RAD_,alb_vis_lnd) = AlbedoGlacierVisMax
      dsl4jsb_var2D_onDomain(RAD_,alb_nir_lnd) = AlbedoGlacierNirMax
    ELSE
      dsl4jsb_var2D_onDomain(RAD_,alb_vis)     = AlbedoVisInitial
      dsl4jsb_var2D_onDomain(RAD_,alb_nir)     = AlbedoNirInitial
      dsl4jsb_var2D_onDomain(RAD_,alb_vis_lnd) = AlbedoVisInitial
      dsl4jsb_var2D_onDomain(RAD_,alb_nir_lnd) = AlbedoNirInitial
    END IF

    !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(RAD_,alb_vis)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(RAD_,alb_nir)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(RAD_,alb_vis_lnd)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(RAD_,alb_nir_lnd)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(RAD_,alb_background))

    ! soil albedo
    IF (tile%contains_soil) THEN
      dsl4jsb_var2D_onDomain(RAD_, alb_vis_soil) = rad_init_vars%alb_vis_soil
      dsl4jsb_var2D_onDomain(RAD_, alb_nir_soil) = rad_init_vars%alb_nir_soil
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(RAD_,alb_nir_soil), dsl4jsb_var2D_onDomain(RAD_, alb_vis_soil))
    END IF

    ! Note: If use_alb_canopy=FALSE and PFTs are used, alb_vis_can and alb_nir_can will be different on the tiles
    !       where the rad process runs (PFTs => values from lctlib) and the parent tiles (values from ini file). These
    !       two variables are not aggregated at the moment and should therefore only be used on the tiles where the
    !       rad process runs! (TODO)
    IF (tile%contains_vegetation) THEN
      ! Get canopy albedo
      IF (dsl4jsb_Config(RAD_)%use_alb_canopy .OR. tile%lcts(1)%lib_id == 0) THEN
        dsl4jsb_var2D_onDomain(RAD_, alb_vis_can) = rad_init_vars%alb_vis_can
      ELSE
        dsl4jsb_var2D_onDomain(RAD_, alb_vis_can) = dsl4jsb_Lctlib_param(AlbedoCanopyVIS)
      END IF
      IF (dsl4jsb_Config(RAD_)%use_alb_canopy .OR. tile%lcts(1)%lib_id == 0) THEN
        dsl4jsb_var2D_onDomain(RAD_, alb_nir_can) = rad_init_vars%alb_nir_can
      ELSE
        dsl4jsb_var2D_onDomain(RAD_, alb_nir_can) = dsl4jsb_Lctlib_param(AlbedoCanopyNIR)
      END IF
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(RAD_, alb_nir_can), dsl4jsb_var2D_onDomain(RAD_, alb_vis_can))

      ! The albedo of snow on canopy:
      ! R: In the model the albedo of snow is weighted with albedo_age_weight between two albedo schemes (age and temperature).
      !    For snow on soil this is done with calculations of snow age and temperature and with differentiation between vis and nir
      !    (see mo_rad_process in subroutine calc_snow_albedo).
      !    However, for snow ON CANOPY the albedo from the age scheme and the albedo from the temperature scheme are simply assumed
      !    to be constant (see the following lines).
      !    Therefore as all inputs to this equation are constants, it was not necessary to calculate it for each time step and block
      !    again (as it was in JSBACH3 in mo_land_surface: update_albedo_snowage_temp). Therefore it put it here and not in
      !    merge_albedos_of_vegtile.
      albedo_age_weight = dsl4jsb_Config(RAD_)%albedo_age_weight

      ! Weight between age and temperature impact on snow albedo for snow on canopy
      dsl4jsb_Config(RAD_)%AlbedoCanopySnow = &
        & albedo_age_weight * AlbedoCanopySnow_age + (1._wp - albedo_age_weight) * AlbedoCanopySnow_temp

    END IF

  END SUBROUTINE rad_init_ic

  !-----------------------------------------------------------------------------------------------------
  !> Intialize radiation process (after memory has been set up)
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE rad_init_q_ic(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    !USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_rad_constants,         ONLY: rfr_ratio_toc

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< one tile with data structure for one lct

    !TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':rad_init_q_ic'

    dsl4jsb_Def_memory(RAD_)

    dsl4jsb_Real2D_onDomain      :: rfr_ratio_boc
    dsl4jsb_Real2D_onDomain      :: rfr_ratio_boc_tvegdyn_mavg
    dsl4jsb_Real3D_onDomain      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_shaded_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_shaded_tcnl_mavg_cl

    IF (.NOT. tile%Is_process_active(RAD_)) RETURN

    IF (debug_on()) CALL message(routine, 'Setting initial conditions of radiation memory (quincy) for tile '// &
      &                                   TRIM(tile%name))

    !model  => Get_model(tile%owner_model_id)

    !dsl4jsb_Get_config(RAD_)

    dsl4jsb_Get_memory(RAD_)

    dsl4jsb_Get_var2D_onDomain(RAD_, rfr_ratio_boc)
    dsl4jsb_Get_var2D_onDomain(RAD_, rfr_ratio_boc_tvegdyn_mavg)
    dsl4jsb_Get_var3D_onDomain(RAD_, ppfd_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(RAD_, ppfd_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(RAD_, ppfd_shaded_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(RAD_, ppfd_shaded_tcnl_mavg_cl)
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    rfr_ratio_boc(:,:)                  = rfr_ratio_toc
    rfr_ratio_boc_tvegdyn_mavg(:,:)     = rfr_ratio_boc(:,:) 

    ppfd_sunlit_tfrac_mavg_cl(:,:,:)    = 150.0_wp
    ppfd_sunlit_tcnl_mavg_cl(:,:,:)     = 150.0_wp
    
    ppfd_shaded_tfrac_mavg_cl(:,:,:)    = 10.0_wp
    ppfd_shaded_tcnl_mavg_cl(:,:,:)     = 10.0_wp

  END SUBROUTINE rad_init_q_ic

#endif
END MODULE mo_rad_init
