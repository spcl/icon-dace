!> Initialization of the the soil and snow memory.
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
MODULE mo_sse_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message
  USE mo_jsb_control,        ONLY: debug_on

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string, soil_init_from_texture

  dsl4jsb_Use_processes SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sse_init, init_soil_properties

  REAL(wp), PARAMETER :: &
    & fao_thermal_diff(0:6) = (/7.4e-7_wp ,8.7e-7_wp ,8.0e-7_wp ,7.4e-7_wp ,7.1e-7_wp ,6.7e-7_wp ,-1._wp/), &
    & fao_vol_hcap    (0:6) = (/2.25e+6_wp,1.93e+6_wp,2.10e+6_wp,2.25e+6_wp,2.36e+6_wp,2.48e+6_wp,-1._wp/)

  ! heat capacity [J/Km^3] (s. table 2.2 in Reick et al. (2020) Reports on Earth System Science 240)
  REAL(wp), PARAMETER ::  heat_cap_sand = 1.93e+06_wp
  REAL(wp), PARAMETER ::  heat_cap_silt = 2.25e+06_wp
  REAL(wp), PARAMETER ::  heat_cap_clay = 2.48e+06_wp
  REAL(wp), PARAMETER ::  heat_cap_oc = 2.5e+06_wp ! see vol_hcap_org in mo_sse_constants.f90

  ! heat conductivity (s. table 2.2 in Reick et al. (2020) Reports on Earth System Science 240)
  REAL(wp), PARAMETER ::  heat_cond_sand = 1.6791_wp
  REAL(wp), PARAMETER ::  heat_cond_silt = 1.665_wp
  REAL(wp), PARAMETER ::  heat_cond_clay = 1.6616_wp
  REAL(wp), PARAMETER ::  heat_cond_oc = 0.25_wp ! see hcond_org in mo_sse_constants.f90

  TYPE t_sse_init_vars
    REAL(wp), POINTER ::     &
      & fao         (:,:  ) => NULL(), &
      & vol_heat_cap(:,:  ) => NULL(), &
      & heat_cond   (:,:  ) => NULL(), &
      & fr_sand     (:,:  ) => Null(), &
      & fr_silt     (:,:  ) => Null(), &
      & fr_clay     (:,:  ) => Null(), &
      & fr_sand_deep(:,:  ) => Null(), &
      & fr_silt_deep(:,:  ) => Null(), &
      & fr_clay_deep(:,:  ) => Null(), &
      & t_srf_mon   (:,:,:) => NULL(), &
      & ifs_tsoil   (:,:,:) => NULL(), &
      & ifs_tskin   (:,:  ) => NULL()
  END TYPE t_sse_init_vars

  TYPE(t_sse_init_vars) :: sse_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_init'

CONTAINS

  !
  !> Intialize soil and snow energy process (after memory has been set up)
  !
  SUBROUTINE sse_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init'

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL sse_read_init_vars(tile)
    END IF

    CALL sse_init_bc(tile)
    CALL sse_init_ic(tile)

    IF (tile%Is_last_process_tile(SSE_)) THEN
      CALL sse_finalize_init_vars()
    END IF

  END SUBROUTINE sse_init

  SUBROUTINE sse_finalize_init_vars

    DEALLOCATE(                     &
      & sse_init_vars%vol_heat_cap, &
      & sse_init_vars%heat_cond     &
      & )

    IF (ASSOCIATED(sse_init_vars%fao))          DEALLOCATE(sse_init_vars%fao)
    IF (ASSOCIATED(sse_init_vars%fr_sand))      DEALLOCATE(sse_init_vars%fr_sand)
    IF (ASSOCIATED(sse_init_vars%fr_silt))      DEALLOCATE(sse_init_vars%fr_silt)
    IF (ASSOCIATED(sse_init_vars%fr_clay))      DEALLOCATE(sse_init_vars%fr_clay)
    IF (ASSOCIATED(sse_init_vars%fr_sand_deep)) DEALLOCATE(sse_init_vars%fr_sand_deep)
    IF (ASSOCIATED(sse_init_vars%fr_silt_deep)) DEALLOCATE(sse_init_vars%fr_silt_deep)
    IF (ASSOCIATED(sse_init_vars%fr_clay_deep)) DEALLOCATE(sse_init_vars%fr_clay_deep)

    IF (ASSOCIATED(sse_init_vars%t_srf_mon))    DEALLOCATE(sse_init_vars%t_srf_mon)

    IF (ASSOCIATED(sse_init_vars%ifs_tsoil))    DEALLOCATE(sse_init_vars%ifs_tsoil)
    IF (ASSOCIATED(sse_init_vars%ifs_tskin))    DEALLOCATE(sse_init_vars%ifs_tskin)

  END SUBROUTINE sse_finalize_init_vars

  SUBROUTINE sse_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(HYDRO_)

    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< temporary pointers

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks, isoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_e => Get_vgrid('soil_depth_energy')

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Reading soil energy boundary conditions for tile '//TRIM(tile%name))

    ! Get config parameters
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    ALLOCATE(                                          &
      & sse_init_vars%vol_heat_cap(nproma, nblks    ), &
      & sse_init_vars%heat_cond   (nproma, nblks    )  &
      & )

    sse_init_vars%vol_heat_cap(:,:)   = missval
    sse_init_vars%heat_cond   (:,:)   = missval

    IF (      .NOT. (dsl4jsb_Config(SSE_)%l_heat_cap_map .AND. dsl4jsb_Config(SSE_)%l_heat_cond_map) &
      & .AND. .NOT. dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      ALLOCATE(sse_init_vars%fao(nproma, nblks))
      sse_init_vars%fao(:,:)   = missval
    END IF

    IF (dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      ALLOCATE(                                          &
        & sse_init_vars%fr_sand     (nproma, nblks    ), &
        & sse_init_vars%fr_silt     (nproma, nblks    ), &
        & sse_init_vars%fr_clay     (nproma, nblks    ), &
        & sse_init_vars%fr_sand_deep(nproma, nblks    ), &
        & sse_init_vars%fr_silt_deep(nproma, nblks    ), &
        & sse_init_vars%fr_clay_deep(nproma, nblks    )  &
        & )

      sse_init_vars%fr_sand     (:,:)   = missval
      sse_init_vars%fr_silt     (:,:)   = missval
      sse_init_vars%fr_clay     (:,:)   = missval
      sse_init_vars%fr_sand_deep(:,:)   = missval
      sse_init_vars%fr_silt_deep(:,:)   = missval
      sse_init_vars%fr_clay_deep(:,:)   = missval
    END IF ! l_soil_texture

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SSE_)%bc_filename), model%grid_id)

    ! soil texture information
    IF (dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      ptr_2D => input_file%Read_2d(             &
        & variable_name='FR_SAND',              &
        & fill_array = sse_init_vars%fr_sand)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) > 0._wp)

      ptr_2D => input_file%Read_2d(             &
        & variable_name='FR_SILT',              &
        & fill_array = sse_init_vars%fr_silt)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) > 0._wp)

      ptr_2D => input_file%Read_2d(             &
        & variable_name='FR_CLAY',              &
        & fill_array = sse_init_vars%fr_clay)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) > 0._wp)

      ptr_2D => input_file%Read_2d(             &
        & variable_name='SUB_FR_SAND',              &
        & fill_array = sse_init_vars%fr_sand_deep)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) > 0._wp)

      ptr_2D => input_file%Read_2d(             &
        & variable_name='SUB_FR_SILT',              &
        & fill_array = sse_init_vars%fr_silt_deep)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) > 0._wp)

      ptr_2D => input_file%Read_2d(             &
        & variable_name='SUB_FR_CLAY',              &
        & fill_array = sse_init_vars%fr_clay_deep)
      ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) > 0._wp)

      ! calculate soil heat capacity [J/(m^3*K)] incl. organic fraction
      CALL soil_init_from_texture( &
         heat_cap_sand, heat_cap_silt, heat_cap_clay, heat_cap_oc,                                          &
         sse_init_vars%fr_sand(:,:), sse_init_vars%fr_silt(:,:), sse_init_vars%fr_clay(:,:),                &
         sse_init_vars%fr_sand_deep(:,:), sse_init_vars%fr_silt_deep(:,:), sse_init_vars%fr_clay_deep(:,:), &
         sse_init_vars%vol_heat_cap(:,:))

      ! calculate soil heat conductivity incl. organic fraction
      CALL soil_init_from_texture( &
         heat_cond_sand, heat_cond_silt, heat_cond_clay, heat_cond_oc,                                      &
         sse_init_vars%fr_sand(:,:), sse_init_vars%fr_silt(:,:), sse_init_vars%fr_clay(:,:),                &
         sse_init_vars%fr_sand_deep(:,:), sse_init_vars%fr_silt_deep(:,:), sse_init_vars%fr_clay_deep(:,:), &
         sse_init_vars%heat_cond(:,:))
    ELSE
      IF (.NOT. dsl4jsb_Config(SSE_)%l_heat_cap_map .OR. .NOT. dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
        ! FAO
        ptr_2D => input_file%Read_2d(variable_name='fao')
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)
        sse_init_vars%fao(:,:) = MAX(REAL(LBOUND(fao_thermal_diff,1),wp), &
                                     MIN(REAL(UBOUND(fao_thermal_diff,1),wp), ptr_2D(:,:)))
      END IF

      IF (dsl4jsb_Config(SSE_)%l_heat_cap_map) THEN
        ! Volumetric heat capacity of dry soil [J/(m^3*K)]

        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                   &
            & variable_name='heat_capacity_mineral',      &
            & fill_array = sse_init_vars%vol_heat_cap)
        ELSE
          ptr_2D => input_file%Read_2d(                   &
            & variable_name='heat_capacity',              &
            & fill_array = sse_init_vars%vol_heat_cap)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 2.e6_wp, ptr_2D(:,:) > 0._wp)
      END IF

      IF (dsl4jsb_Config(SSE_)%l_heat_cond_map) THEN
        ! Heat conductivity of mineral soil
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                   &
            & variable_name='heat_conductivity_mineral',  &
            & fill_array = sse_init_vars%heat_cond)
        ELSE
          ptr_2D => input_file%Read_2d(                   &
            & variable_name='heat_conductivity',          &
            & fill_array = sse_init_vars%heat_cond)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) > 0._wp)
      END IF

    END IF ! l_soil_texture

    CALL input_file%Close()

    IF (model%config%init_from_ifs) THEN

      ALLOCATE(sse_init_vars%ifs_tsoil(nproma,ifs_nsoil,nblks))
      sse_init_vars%ifs_tsoil(:,:,:) = missval
      ALLOCATE(sse_init_vars%ifs_tskin(nproma,nblks))
      sse_init_vars%ifs_tskin(:,:) = missval
      DO isoil=1,ifs_nsoil
        ptr_2D => model%config%ifs_input_file%Read_2d_1lev_1time( &
          & variable_name='STL'//int2string(isoil))
        sse_init_vars%ifs_tsoil(:,isoil,:) = ptr_2D(:,:)
      END DO
      ptr_3D => model%config%ifs_input_file%Read_2d_time( &
        & variable_name='SKT', start_time_step=1, end_time_step=1)
      sse_init_vars%ifs_tskin(:,:) = ptr_3D(:,:,1)

    ELSE

      ALLOCATE(sse_init_vars%t_srf_mon(nproma, nblks, 12))
      sse_init_vars%t_srf_mon(:,:,:) = missval

      input_file = jsb_netcdf_open_input(dsl4jsb_Config(SSE_)%ic_filename, model%grid_id)

      ! Read climatological surface temperature
      ptr_3D => input_file%Read_2d_time( &
        & variable_name='surf_temp')
      sse_init_vars%t_srf_mon(:,:,:) = MERGE(ptr_3D, 280._wp, ptr_3D > 0._wp)

      CALL input_file%Close()

    END IF

  END SUBROUTINE sse_read_init_vars

  SUBROUTINE sse_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)

    INTEGER :: iblk

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting soil energy boundary conditions for tile '//TRIM(tile%name))

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)

    ! Volumetric heat capacity of dry soil [J/(m^3*K)]
    IF (dsl4jsb_Config(SSE_)%l_heat_cap_map .OR. dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      dsl4jsb_var_ptr(SSE_, vol_heat_cap) = sse_init_vars%vol_heat_cap
    ELSE
      DO iblk=1,SIZE(sse_init_vars%fao,2)
        dsl4jsb_var_ptr(SSE_,vol_heat_cap)(:,iblk) = fao_vol_hcap(NINT(sse_init_vars%fao(:,iblk)))
      END DO
    END IF

    ! Heat conductivity of mineral soil
    IF (dsl4jsb_Config(SSE_)%l_heat_cond_map .OR. dsl4jsb_Config(SSE_)%l_soil_texture) THEN
      dsl4jsb_var_ptr(SSE_, heat_cond) = sse_init_vars%heat_cond
    ELSE
      ! Thermal diffusivity of soil from FAO
      DO iblk=1,SIZE(sse_init_vars%fao,2)
        dsl4jsb_var_ptr(SSE_,thermal_diffusivity)(:,iblk) = fao_thermal_diff(NINT(sse_init_vars%fao(:,iblk)))
      END DO
      dsl4jsb_var2D_onDomain(SSE_, heat_cond) = &
        & dsl4jsb_var2D_onDomain(SSE_, vol_heat_cap) * dsl4jsb_var2D_onDomain(SSE_, thermal_diffusivity)
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(SSE_,thermal_diffusivity)) ASYNC(1)
    END IF

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SSE_,vol_heat_cap)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SSE_,heat_cond))

  END SUBROUTINE sse_init_bc

  SUBROUTINE sse_init_ic(tile)

    USE mo_sse_process, ONLY: init_soil_temperature
    USE mo_util,        ONLY: ifs2soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)

    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':sse_init_ic'

    model => get_model(tile%owner_model_id)
    soil_e => Get_vgrid('soil_depth_energy')

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)

    ! Get soil config
    dsl4jsb_Get_config(SSE_)

    IF (.NOT. tile%Is_process_active(SSE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing soil energy memory for tile '//TRIM(tile%name))

    ! Soil and surface temperatures
    !
    IF (model%config%init_from_ifs) THEN

      ! CALL ifs2soil(sse_init_vars%ifs_tsoil(:,:,:), ifs_soil_depth(:), dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,:), soil_e%ubounds(:), &
      !   & ifs_sfc=sse_init_vars%ifs_tskin(:,:))
      CALL ifs2soil(sse_init_vars%ifs_tsoil(:,:,:), ifs_soil_depth(:), dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,:), soil_e%ubounds(:))

    ELSE
      DO iblk=1,SIZE(sse_init_vars%t_srf_mon, DIM=2)
        CALL init_soil_temperature(                   &
          & soil_e%levels(:),                         &
          & sse_init_vars%t_srf_mon(:,iblk,1:12),     &
          & dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,iblk))
      END DO
    END IF

    !$ACC UPDATE DEVICE(dsl4jsb_var3D_onDomain(SSE_,t_soil_sl)) ASYNC(1)

    ! IF (TRIM(tile%name) /= 'box') THEN
    IF ( model%config%use_tmx .OR. TRIM(tile%name) == 'land') THEN
      dsl4jsb_var2D_onDomain(SEB_,t)        = dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,1,:)
      dsl4jsb_var2D_onDomain(SEB_,t_unfilt) = dsl4jsb_var2D_onDomain(SEB_,t)

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SEB_,t_unfilt)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(SEB_,t))
    END IF

    ! Note: called from model_init now
    ! CALL init_soil_properties(tile)

  END SUBROUTINE sse_init_ic

  SUBROUTINE init_soil_properties(tile)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min, tmelt, alf, rhoh2o
    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, &   ! vol_hcap_snow_min,
      &                                  vol_hcap_org_top, vol_hcap_org_below, &
      &                                  hcond_ice, hcond_snow, &         !, hcond_snow_min, hcond_org
      &                                  hcond_org_top, hcond_dry_org_top, hcond_org_below, hcond_dry_org_below
    USE mo_hydro_constants,        ONLY: vol_porosity_org_top, vol_porosity_org_below
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity, Get_liquid_max
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & vol_heat_cap, &
      & heat_cond, &
      & vol_heat_cap_snow, &
      & heat_cond_snow
    dsl4jsb_Real3D_onDomain :: &
      & soil_depth_sl, &
      & wtr_soil_sl, &
      & ice_soil_sl, &
      & t_soil_sl, &
      & matric_pot_sl, &
      & bclapp_sl, &
      & fract_org_sl, &
      & vol_heat_cap_sl, &
      & heat_cond_sl, &
      & vol_porosity_sl

    ! Locally allocated vectors
    !
    dsl4jsb_Real3D_onDomain :: &
      & liquid_max,            &
      & fract_water,           &
      & fract_ice,             &
      & vol_heat_cap_tmp,      &
      & vol_heat_cap_org_tmp,  &
      & heat_cond_tmp,         &
      & vol_porosity_org_tmp,  &
      & hcond_org_tmp,         &
      & hcond_dry_org_tmp

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    INTEGER  :: nsoil, nproma, nblks, ic, is, ib
    REAL(wp), POINTER :: soil_e_dz(:)
    REAL(wp) ::heat_cap_sl, ws_max_sl

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    model => Get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels
    soil_e_dz => soil_e%dz

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for domain
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var3D_onDomain(SSE_, vol_heat_cap_sl)
    dsl4jsb_Get_var3D_onDomain(SSE_, heat_cond_sl)
    dsl4jsb_Get_var2D_onDomain(SSE_, vol_heat_cap_snow)
    dsl4jsb_Get_var2D_onDomain(SSE_, heat_cond_snow)

    nproma = SIZE(heat_cond_snow,1)
    nblks  = SIZE(heat_cond_snow,2)

    ! On glacier tile, set variables and return
    IF (tile%is_glacier) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(3)
      DO is=1,nsoil
        DO ib=1,nblks
          DO ic=1,nproma
            vol_heat_cap_sl(ic,is,ib) = vol_hcap_ice
            heat_cond_sl   (ic,is,ib) = hcond_ice
          END DO
        END DO
      END DO
      !$ACC END LOOP
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO ib=1,nblks
        DO ic=1,nproma
          vol_heat_cap_snow(ic  ,ib) = vol_hcap_snow
          heat_cond_snow   (ic  ,ib) = hcond_snow
        END DO
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL
      RETURN
    END IF

    dsl4jsb_Get_var2D_onDomain(SSE_,   vol_heat_cap)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)

    IF (model%config%init_from_ifs .AND. dsl4jsb_Config(SSE_)%l_freeze) THEN

      dsl4jsb_Get_var3D_onDomain(SSE_,   t_soil_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_, matric_pot_sl)
      dsl4jsb_Get_var3D_onDomain(HYDRO_, bclapp_sl)

      ALLOCATE(liquid_max (nproma, nsoil, nblks))
      !$ACC ENTER DATA CREATE(liquid_max)

      IF (dsl4jsb_Config(SSE_)%l_supercool) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1) PRIVATE(ws_max_sl)
        DO is = 1, nsoil
          DO ib = 1, nblks
            DO ic = 1, nproma
              ws_max_sl = soil_depth_sl(ic,is,ib) * vol_porosity_sl(ic,is,ib)
              liquid_max(ic,is,ib) = Get_liquid_max( &
                & t_soil_sl    (ic,is,ib), &
                & ws_max_sl              , &
                & matric_pot_sl(ic,is,ib), &
                & bclapp_sl    (ic,is,ib)  &
                & )
            END DO
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
        DO is = 1, nsoil
          DO ib = 1, nblks
            DO ic = 1, nproma
              liquid_max(ic,is,ib) = 0._wp
            END DO
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            heat_cap_sl = soil_e_dz(is) * vol_heat_cap(ic,ib)
            ice_soil_sl(ic,is,ib) = 0._wp
            IF (t_soil_sl(ic,is,ib) < tmelt .AND. wtr_soil_sl(ic,is,ib) > liquid_max(ic,is,ib)) THEN
              ice_soil_sl(ic,is,ib) = MIN(wtr_soil_sl(ic,is,ib) - liquid_max(ic,is,ib),                &
                &                         heat_cap_sl * (tmelt - t_soil_sl(ic,is,ib)) / (alf * rhoh2o) &
                &                        )
            END IF
            wtr_soil_sl(ic,is,ib) = wtr_soil_sl(ic,is,ib) - ice_soil_sl(ic,is,ib)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      ! t_soil_sl(:,:,:) = t_soil_sl(:,:,:) + ice_soil_sl(:,:,:) * alf * rhoh2o / heat_cap_sl(:,:,:)
      ! IF (TRIM(tile%name) == 'land') THEN
      !   dsl4jsb_var_ptr(SEB_,t)(:,:)        = t_soil_sl(:,1,:)
      !   dsl4jsb_var_ptr(SEB_,t_unfilt)(:,:) = t_soil_sl(:,1,:)
      ! END IF

      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(liquid_max)
      DEALLOCATE(liquid_max)

      ! dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl) = &
      !   & MIN(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl), dsl4jsb_var3D_onDomain(HYDRO_, w_soil_fc_sl))
      ! !$acc update device(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl))

    END IF

    dsl4jsb_Get_var2D_onDomain(SSE_,   heat_cond)

    IF (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onDomain(HYDRO_, fract_org_sl)
    ELSE
      ALLOCATE(fract_org_sl(nproma, nsoil, nblks))
      !$ACC ENTER DATA CREATE(fract_org_sl)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            fract_org_sl(ic,is,ib) = 0._wp
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ALLOCATE(vol_heat_cap_tmp(nproma,nsoil,nblks), heat_cond_tmp(nproma,nsoil,nblks))
    !$ACC ENTER DATA CREATE(vol_heat_cap_tmp, heat_cond_tmp)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
    DO is = 1, nsoil
      DO ib = 1, nblks
        DO ic = 1, nproma
          vol_heat_cap_tmp(ic,is,ib) = vol_heat_cap(ic,ib)
          heat_cond_tmp   (ic,is,ib) = heat_cond(ic,ib)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_dyn) THEN
      ALLOCATE(fract_water(nproma, nsoil, nblks), fract_ice(nproma, nsoil, nblks), &
        &      vol_heat_cap_org_tmp(nproma, nsoil, nblks))
      !$ACC ENTER DATA CREATE(fract_water, fract_ice, vol_heat_cap_org_tmp)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO ic = 1, nproma
        DO is = 1, nsoil
          DO ib = 1, nblks
            IF (soil_depth_sl(ic,is,ib) > 0._wp) THEN
              fract_water(ic,is,ib) = wtr_soil_sl(ic,is,ib) / soil_depth_sl(ic,is,ib)
              fract_ice  (ic,is,ib) = ice_soil_sl(ic,is,ib) / soil_depth_sl(ic,is,ib)
            ELSE
              fract_water(ic,is,ib) = 0._wp
              fract_ice  (ic,is,ib) = 0._wp
            END IF
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC UPDATE HOST(fract_water, fract_ice)
      !$ACC WAIT(1)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO ib = 1, nblks
        DO ic = 1, nproma
           vol_heat_cap_org_tmp(ic,1,ib) = vol_hcap_org_top
        END DO
      END DO
      !$ACC END PARALLEL LOOP
     !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 2, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            vol_heat_cap_org_tmp(ic,is,ib) = vol_hcap_org_below
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            vol_heat_cap_sl (ic,is,ib) = calc_vol_heat_capacity( &
              &                                 vol_heat_cap_tmp (ic,is,ib)      , &
              &                                 vol_heat_cap_org_tmp(ic,is,ib)   , &
              &                                 vol_porosity_sl  (ic,is,ib)      , &
              &                                 fract_org_sl     (ic,is,ib)      , &
              &                                 fract_water      (ic,is,ib)      , &
              &                                 fract_ice        (ic,is,ib)      , &
              &                                 soil_depth_sl    (ic,is,ib)      , &
              &                                 soil_e_dz        (   is   )        &
              & )
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(fract_water, fract_ice, vol_heat_cap_org_tmp)
      DEALLOCATE(fract_water, fract_ice, vol_heat_cap_org_tmp)
    ELSE

      !$ACC UPDATE DEVICE(vol_heat_cap) ASYNC(1)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib=1,nblks
          DO ic = 1, nproma
            vol_heat_cap_sl(ic,is,ib) = vol_heat_cap(ic,ib)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_dyn) THEN
      ALLOCATE( &
        & vol_porosity_org_tmp(nproma,nsoil,nblks), &
        & hcond_org_tmp(nproma,nsoil,nblks),        &
        & hcond_dry_org_tmp(nproma,nsoil,nblks))
      !$ACC ENTER DATA CREATE(vol_porosity_org_tmp, hcond_org_tmp, hcond_dry_org_tmp)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO ib = 1, nblks
        DO ic = 1, nproma
          vol_porosity_org_tmp(ic,1,ib) = vol_porosity_org_top
          hcond_org_tmp(ic,1,ib)        = hcond_org_top
          hcond_dry_org_tmp(ic,1,ib)    = hcond_dry_org_top
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 2, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            vol_porosity_org_tmp(ic,is,ib) = vol_porosity_org_below
            hcond_org_tmp(ic,is,ib)        = hcond_org_below
            hcond_dry_org_tmp(ic,is,ib)    = hcond_dry_org_below
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 2, nsoil
        DO ib = 1, nblks
          DO ic = 1, nproma
            heat_cond_sl(ic,is,ib) = calc_thermal_conductivity( &
              &                                  heat_cond_tmp       (ic,is,ib)   , &
              &                                  wtr_soil_sl         (ic,is,ib)   , &
              &                                  ice_soil_sl         (ic,is,ib)   , &
              &                                  vol_porosity_sl     (ic,is,ib)   , &
              &                                  fract_org_sl        (ic,is,ib)   , &
              &                                  vol_porosity_org_tmp(ic,is,ib)   , &
              &                                  hcond_org_tmp       (ic,is,ib)   , &
              &                                  hcond_dry_org_tmp   (ic,is,ib)   , &
              &                                  soil_depth_sl       (ic,is,ib)   , &
              &                                  soil_e_dz           (   is   )     &
              & )
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(vol_porosity_org_tmp, hcond_org_tmp, hcond_dry_org_tmp)
      DEALLOCATE(vol_porosity_org_tmp,hcond_org_tmp,hcond_dry_org_tmp)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1)
      DO is = 1, nsoil
        DO ib=1,nblks
          DO ic = 1, nproma
            heat_cond_sl(ic,is,ib) = heat_cond(ic,ib)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (.NOT. (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil)) THEN
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      IF (model%config%l_compat401) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO ib=1,nblks
          DO ic = 1, nproma
            vol_heat_cap_snow(ic,ib) = 634500._wp
            heat_cond_snow   (ic,ib) = 0.31_wp
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
          ! Use snow density values for fresh snow for initialization
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO ib=1,nblks
            DO ic = 1, nproma
              vol_heat_cap_snow(ic,ib) = ci * dens_snow_min              ! Use specific heat of ice (Verseghy 1991)
              heat_cond_snow   (ic,ib) =   2.5e-6_wp * dens_snow_min**2._wp &
                &                        - 0.123E-3_wp *dens_snow_min + 0.024_wp ! Calonne et al. (2011)
              ! vol_heat_cap_snow(ic,ib) = vol_hcap_snow_min
              ! heat_cond_snow   (ic,ib) = hcond_snow_min
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO ib=1,nblks
            DO ic = 1, nproma
              vol_heat_cap_snow(ic,ib) = vol_hcap_snow
              heat_cond_snow   (ic,ib) = hcond_snow
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      END IF
    END IF

    !$ACC WAIT(1)

  END SUBROUTINE init_soil_properties

#endif
END MODULE mo_sse_init
