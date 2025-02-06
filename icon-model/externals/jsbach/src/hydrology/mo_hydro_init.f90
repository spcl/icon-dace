!> Initialization of the the hydrology memory.
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
MODULE mo_hydro_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, warning, finish
  USE mo_jsb_control,        ONLY: debug_on, jsbach_runs_standalone

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string, soil_init_from_texture

  dsl4jsb_Use_processes HYDRO_, PHENO_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)

  IMPLICIT NONE
  PRIVATE

  ! parameter taken from TERRA for sand, loam (here called silt), clay, peat (here called oc), respectively
  ! (s. sfc_terra_data.f90)

  ! pore volume (fraction of volume, called cporv in TERRA)
  REAL(wp), PARAMETER ::  porosity_sand = 0.364_wp
  REAL(wp), PARAMETER ::  porosity_silt = 0.455_wp
  REAL(wp), PARAMETER ::  porosity_clay = 0.507_wp
  REAL(wp), PARAMETER ::  porosity_oc = 0.863_wp

  ! field capacity (fraction of volume, called cfcap in TERRA)
  REAL(wp), PARAMETER ::  field_cap_sand = 0.196_wp
  REAL(wp), PARAMETER ::  field_cap_silt = 0.340_wp
  REAL(wp), PARAMETER ::  field_cap_clay = 0.463_wp
  REAL(wp), PARAMETER ::  field_cap_oc = 0.763_wp

  ! hydr. conductivity at saturation ([m/s], called ckw0 in TERRA)
  REAL(wp), PARAMETER ::  hyd_cond_sand = 4.79e-5_wp
  REAL(wp), PARAMETER ::  hyd_cond_silt = 5.31e-6_wp
  REAL(wp), PARAMETER ::  hyd_cond_clay = 8.50e-8_wp
  REAL(wp), PARAMETER ::  hyd_cond_oc = 5.80e-8_wp

  ! wilting point (fraction of volume, called cpwp in TERRA)
  REAL(wp), PARAMETER ::  wilt_sand = 0.042_wp
  REAL(wp), PARAMETER ::  wilt_silt = 0.11_wp
  REAL(wp), PARAMETER ::  wilt_clay = 0.257_wp
  REAL(wp), PARAMETER ::  wilt_oc = 0.265_wp

  ! pore size index
  REAL(wp), PARAMETER :: pore_size_index_sand = 0.35_wp
  REAL(wp), PARAMETER :: pore_size_index_silt = 0.2_wp
  REAL(wp), PARAMETER :: pore_size_index_clay = 0.13_wp
  REAL(wp), PARAMETER :: pore_size_index_oc = 0.65_wp

  ! Clapp & Hornberger (s. Beringer et al 2001)
  REAL(wp), PARAMETER ::  bclapp_sand = 3.39_wp   ! (type 1 in Beringer)
  REAL(wp), PARAMETER ::  bclapp_silt = 4.98_wp   ! (type 5 in Beringer)
  REAL(wp), PARAMETER ::  bclapp_clay = 10.38_wp  ! (type 10 in Beringer)
  REAL(wp), PARAMETER ::  bclapp_oc = 4._wp       ! (type 12 in Beringer)

  ! soil-matric potential (s. Beringer et al 2001)
  REAL(wp), PARAMETER ::  matric_pot_sand = -0.04729_wp ! (type 1 in Beringer)
  REAL(wp), PARAMETER ::  matric_pot_silt = -0.45425_wp ! (type 5 in Beringer)
  REAL(wp), PARAMETER ::  matric_pot_clay = -0.633_wp   ! (type 10 in Beringer)
  REAL(wp), PARAMETER ::  matric_pot_oc = -0.12_wp      ! (type 12 in Beringer)

  ! factor to compensate profile of hyd_cond_sat with depth in TERRA
  REAL(wp), PARAMETER ::  hyd_cond_sat_profile = 0.432332_wp

  ! residual soil moisture (fraction of volume; Maidment, Handbook of Hydrology, 1993)
  REAL(wp), PARAMETER ::  wres_sand = 0.020_wp  ! (type sand)
  REAL(wp), PARAMETER ::  wres_silt = 0.015_wp  ! (type silt loam)
  REAL(wp), PARAMETER ::  wres_clay = 0.090_wp  ! (type clay)
  REAL(wp), PARAMETER ::  wres_oc   = 0.150_wp  ! (Letts et al., 2000)

  PUBLIC :: hydro_init, hydro_sanitize_state

  TYPE t_hydro_init_vars
    REAL(wp), POINTER ::                  &
      & elevation        (:,:) => NULL(), &
      & oro_stddev       (:,:) => NULL(), &
      & soil_depth       (:,:) => NULL(), &
      & vol_porosity     (:,:) => NULL(), &
      & hyd_cond_sat     (:,:) => NULL(), &
      & matric_pot       (:,:) => NULL(), &
      & bclapp           (:,:) => NULL(), &
      & pore_size_index  (:,:) => NULL(), &
      & vol_field_cap    (:,:) => NULL(), &
      & vol_p_wilt       (:,:) => NULL(), &
      & vol_wres         (:,:) => NULL(), &
      & root_depth       (:,:) => NULL(), &
      & weq_snow_soil    (:,:) => NULL(), &
      & fract_pond_max   (:,:) => NULL(), &
      & depth_pond_max   (:,:) => NULL(), &
      & wtr_soil_sl    (:,:,:) => NULL(), &
      & fract_org_sl   (:,:,:) => NULL(), &
      & fr_sand          (:,:) => NULL(), &
      & fr_silt          (:,:) => NULL(), &
      & fr_clay          (:,:) => NULL(), &
      & fr_sand_deep     (:,:) => NULL(), &
      & fr_silt_deep     (:,:) => NULL(), &
      & fr_clay_deep     (:,:) => NULL(), &
      & ifs_smi_sl     (:,:,:) => NULL(), &
      & ifs_weq_snow     (:,:) => NULL()
  END TYPE t_hydro_init_vars

  TYPE(t_hydro_init_vars) :: hydro_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_init'

CONTAINS

  !
  !> Intialize hydrology process (after memory has been set up)
  !
  SUBROUTINE hydro_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init'

    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL hydro_read_init_vars(tile)
    END IF

    CALL hydro_init_bc(tile)
    CALL hydro_init_ic(tile)

    IF (tile%Is_last_process_tile(HYDRO_)) THEN
      CALL hydro_finalize_init_vars()
    END IF

  END SUBROUTINE hydro_init

  SUBROUTINE hydro_finalize_init_vars

    DEALLOCATE( &
      & hydro_init_vars%oro_stddev      ,       &
      & hydro_init_vars%soil_depth      ,       &
      & hydro_init_vars%vol_porosity    ,       &
      & hydro_init_vars%hyd_cond_sat    ,       &
      & hydro_init_vars%matric_pot      ,       &
      & hydro_init_vars%bclapp          ,       &
      & hydro_init_vars%pore_size_index ,       &
      & hydro_init_vars%vol_field_cap   ,       &
      & hydro_init_vars%vol_p_wilt      ,       &
      & hydro_init_vars%vol_wres        ,       &
      & hydro_init_vars%root_depth      ,       &
      & hydro_init_vars%fract_org_sl    ,       &
      & hydro_init_vars%fract_pond_max  ,       &
      & hydro_init_vars%depth_pond_max  ,       &
      & hydro_init_vars%weq_snow_soil)

    IF (ASSOCIATED(hydro_init_vars%elevation))    DEALLOCATE(hydro_init_vars%elevation)
    IF (ASSOCIATED(hydro_init_vars%fr_sand))      DEALLOCATE(hydro_init_vars%fr_sand)
    IF (ASSOCIATED(hydro_init_vars%fr_silt))      DEALLOCATE(hydro_init_vars%fr_silt)
    IF (ASSOCIATED(hydro_init_vars%fr_clay))      DEALLOCATE(hydro_init_vars%fr_clay)
    IF (ASSOCIATED(hydro_init_vars%fr_sand_deep)) DEALLOCATE(hydro_init_vars%fr_sand_deep)
    IF (ASSOCIATED(hydro_init_vars%fr_silt_deep)) DEALLOCATE(hydro_init_vars%fr_silt_deep)
    IF (ASSOCIATED(hydro_init_vars%fr_clay_deep)) DEALLOCATE(hydro_init_vars%fr_clay_deep)
    IF (ASSOCIATED(hydro_init_vars%wtr_soil_sl))  DEALLOCATE(hydro_init_vars%wtr_soil_sl)
    IF (ASSOCIATED(hydro_init_vars%ifs_smi_sl))   DEALLOCATE(hydro_init_vars%ifs_smi_sl)
    IF (ASSOCIATED(hydro_init_vars%ifs_weq_snow)) DEALLOCATE(hydro_init_vars%ifs_weq_snow)

  END SUBROUTINE hydro_finalize_init_vars

  SUBROUTINE hydro_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_)

    REAL(wp), POINTER :: &
      & ptr_2D(:,  :),   &
      & ptr_3D(:,:,:)
#ifdef __ICON__
    LOGICAL,  POINTER   :: is_in_domain(:,:) ! T: cell in domain (not halo)
#endif

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w
    TYPE(t_input_file) :: input_file
    INTEGER :: isoil, nproma, nblks, nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_read_init_vars'


    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Reading/setting hydrology init vars')

    model => Get_model(tile%owner_model_id)

    ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels
#ifdef __ICON__
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
#endif

    ALLOCATE( &
      & hydro_init_vars%oro_stddev      (nproma, nblks),       &
      & hydro_init_vars%soil_depth      (nproma, nblks),       &
      & hydro_init_vars%vol_porosity    (nproma, nblks),       &
      & hydro_init_vars%hyd_cond_sat    (nproma, nblks),       &
      & hydro_init_vars%matric_pot      (nproma, nblks),       &
      & hydro_init_vars%bclapp          (nproma, nblks),       &
      & hydro_init_vars%pore_size_index (nproma, nblks),       &
      & hydro_init_vars%vol_field_cap   (nproma, nblks),       &
      & hydro_init_vars%vol_p_wilt      (nproma, nblks),       &
      & hydro_init_vars%vol_wres        (nproma, nblks),       &
      & hydro_init_vars%root_depth      (nproma, nblks),       &
      & hydro_init_vars%fract_pond_max  (nproma, nblks),       &
      & hydro_init_vars%depth_pond_max  (nproma, nblks),       &
      & hydro_init_vars%weq_snow_soil   (nproma, nblks)        &
      & )

      hydro_init_vars%oro_stddev     (:,:)   = missval
      hydro_init_vars%soil_depth     (:,:)   = missval
      hydro_init_vars%vol_porosity   (:,:)   = missval
      hydro_init_vars%hyd_cond_sat   (:,:)   = missval
      hydro_init_vars%matric_pot     (:,:)   = missval
      hydro_init_vars%bclapp         (:,:)   = missval
      hydro_init_vars%pore_size_index(:,:)   = missval
      hydro_init_vars%vol_field_cap  (:,:)   = missval
      hydro_init_vars%vol_p_wilt     (:,:)   = missval
      hydro_init_vars%vol_wres       (:,:)   = missval
      hydro_init_vars%root_depth     (:,:)   = missval
      hydro_init_vars%fract_pond_max (:,:)   = missval
      hydro_init_vars%depth_pond_max (:,:)   = missval
      hydro_init_vars%weq_snow_soil  (:,:)   = missval

    ALLOCATE(hydro_init_vars%fract_org_sl(nproma, nsoil, nblks))
    hydro_init_vars%fract_org_sl(:,:,:) = missval

    IF (jsbach_runs_standalone()) THEN
       ALLOCATE(hydro_init_vars%elevation(nproma,nblks))
       hydro_init_vars%elevation(:,:) = missval
    END IF
    IF (model%config%init_from_ifs) THEN
      ALLOCATE(hydro_init_vars%ifs_weq_snow(nproma,nblks))
      hydro_init_vars%ifs_weq_snow(:,:) = missval
    END IF
    IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
      ALLOCATE(hydro_init_vars%ifs_smi_sl(nproma,ifs_nsoil,nblks))
      hydro_init_vars%ifs_smi_sl(:,:,:) = missval
    ELSE
      ALLOCATE(hydro_init_vars%wtr_soil_sl(nproma, nsoil, nblks))
      hydro_init_vars%wtr_soil_sl(:,:,:) = missval
    END IF
    IF (dsl4jsb_Config(HYDRO_)%l_soil_texture) THEN
      ALLOCATE( &
        & hydro_init_vars%fr_sand         (nproma, nblks),       &
        & hydro_init_vars%fr_silt         (nproma, nblks),       &
        & hydro_init_vars%fr_clay         (nproma, nblks),       &
        & hydro_init_vars%fr_sand_deep    (nproma, nblks),       &
        & hydro_init_vars%fr_silt_deep    (nproma, nblks),       &
        & hydro_init_vars%fr_clay_deep    (nproma, nblks)        &
        & )
      hydro_init_vars%fr_sand        (:,:)   = missval
      hydro_init_vars%fr_silt        (:,:)   = missval
      hydro_init_vars%fr_clay        (:,:)   = missval
      hydro_init_vars%fr_sand_deep   (:,:)   = missval
      hydro_init_vars%fr_silt_deep   (:,:)   = missval
      hydro_init_vars%fr_clay_deep   (:,:)   = missval
    END IF ! l_soil_texture

    !----------

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_sso_filename), model%grid_id)

    ! Orography ...
    IF (jsbach_runs_standalone()) THEN
       ptr_2D => input_file%Read_2d(  &
         & variable_name='elevation', &
         & fill_array = hydro_init_vars%elevation)
       hydro_init_vars%elevation = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
    END IF

    ptr_2D => input_file%Read_2d(  &
      & variable_name='orostd',    &
      & fill_array = hydro_init_vars%oro_stddev)
    hydro_init_vars%oro_stddev = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    ! Read in maximum fraction and depth of local surface depressions
    IF (dsl4jsb_Config(HYDRO_)%l_ponds .AND. .NOT. tile%is_lake .AND. .NOT. tile%is_glacier) THEN
      IF (input_file%Has_var('surf_depr_fract') .AND. input_file%Has_var('surf_depr_depth')) THEN
        ! Use depression area fraction as maximum pond fraction
        ptr_2D => input_file%Read_2d(                                 &
          & variable_name='surf_depr_fract',                           &
          & fill_array=hydro_init_vars%fract_pond_max)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)
        ! Use depression depth as maximum pond depth (Attention: Values are negative!)
        ptr_2D => input_file%Read_2d(                                 &
          & variable_name='surf_depr_depth',                           &
          & fill_array=hydro_init_vars%depth_pond_max)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) <= 0._wp)
      ELSE
        CALL finish(TRIM(routine), '*** Error: surf_depr_fract and/or surf_depr_depth not found in bc_sso file. '// &
          & 'Pond scheme cannot run without these fields')
      END IF
    END IF

    CALL input_file%Close()

    IF (tile%contains_soil) THEN

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_filename), model%grid_id)

      ! Total soil depth until bedrock
      ptr_2D => input_file%Read_2d(   &
        & variable_name='soil_depth', &
        & fill_array = hydro_init_vars%soil_depth)
      hydro_init_vars%soil_depth = MAX(0.1_wp, MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp))

      !> reads maps if not using the soil texture for calculating the soil properties
      !>
      IF (.NOT. dsl4jsb_Config(HYDRO_)%l_soil_texture) THEN
        ! Volumetric soil porosity
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & fill_array = hydro_init_vars%vol_porosity, &
            & variable_name='soil_porosity_mineral')
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & fill_array = hydro_init_vars%vol_porosity, &
            & variable_name='soil_porosity')
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp)

        ! Saturated hydraulic conductivity
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='hyd_cond_sat_mineral',      &
            & fill_array = hydro_init_vars%hyd_cond_sat)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='hyd_cond_sat',              &
            & fill_array = hydro_init_vars%hyd_cond_sat)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 4.e-6_wp, ptr_2D(:,:) >= 0._wp)

        ! Volumetric field capacity [m/m]
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='soil_field_cap_mineral',    &
            & fill_array = hydro_init_vars%vol_field_cap)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='soil_field_cap',            &
            & fill_array = hydro_init_vars%vol_field_cap)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Volumetric permanent wilting point
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='wilting_point_mineral',     &
            & fill_array=hydro_init_vars%vol_p_wilt)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='wilting_point',             &
            & fill_array=hydro_init_vars%vol_p_wilt)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Volumetric residual water content (if missing deduce from porosity)
        IF (input_file%Has_var('residual_water')) THEN
          IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
            ptr_2D => input_file%Read_2d(                &
              & variable_name='residual_water_mineral',  &
              & fill_array=hydro_init_vars%vol_wres)
          ELSE
            ptr_2D => input_file%Read_2d(                &
              & variable_name='residual_water',          &
              & fill_array=hydro_init_vars%vol_wres)
          END IF
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)
        ELSE
          WRITE (message_text,*) 'BC File does not contain data on residual water content. Will use approximation instead!'
          CALL message (TRIM(routine), message_text)
          hydro_init_vars%vol_wres = hydro_init_vars%vol_porosity * 0.2_wp
        END IF

        ! Soil matric potential
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='moisture_pot_mineral',      &
            & fill_array = hydro_init_vars%matric_pot)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='moisture_pot',              &
            & fill_array = hydro_init_vars%matric_pot)
        END IF
        ! @todo There is still old bc data around where the matric potential is given with a positiv sign while
        !       it should be negativ from a soil physics perspective and is expected as such from the soil hydrology
        !       routines. Until this bc data is not used anymore, the following check is required to ensure the
        !       sign of the soil matric potential is (converted) correct.
#ifdef __ICON__
        ! Exclude values outside of the domain (not available in the echam infrastructure)
        WHERE(.NOT. is_in_domain(:,:))
          ptr_2D(:,:) = 0._wp
        ENDWHERE
#endif
        ! Exclude (standard) missing values
        WHERE(ABS(ptr_2D(:,:)) > 1.0e+20_wp)
          ptr_2D(:,:) = 0._wp
        ENDWHERE
        IF (ALL(ptr_2D(:,:) < 1.0e-10_wp)) THEN
          ! Correct: all values are negative
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), -0.2_wp, ptr_2D(:,:) <= 0._wp)
        ELSE IF (ALL(ptr_2D(:,:) > -1.0e-10_wp)) THEN
          ! False: all values are positive --> convert to negative values
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp) * (-1._wp)
        ELSE
          ! Totally wrong: dataset contains positive and negativ values
          WRITE (message_text,*) 'Found positive and negative values in soil moisture matric potential: ', &
            & MINVAL(ptr_2D(:,:)),' <-> ',MAXVAL(ptr_2D(:,:))
          CALL finish (TRIM(routine), message_text)
        END IF

        ! Exponent B in Clapp and Hornberger
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='bclapp_mineral',            &
            & fill_array = hydro_init_vars%bclapp)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='bclapp',                    &
            & fill_array = hydro_init_vars%bclapp)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 6._wp, ptr_2D(:,:) >= 0._wp)

        ! Pore size distribution index
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='pore_size_index_mineral',   &
            & fill_array = hydro_init_vars%pore_size_index)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='pore_size_index',           &
            & fill_array = hydro_init_vars%pore_size_index)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.25_wp, ptr_2D(:,:) >= 0.05_wp)

      ELSE   ! l_soil_texture

        ! Soil texture information

        ! fraction of sand upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_SAND',                   &
          & fill_array=hydro_init_vars%fr_sand)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! fraction of silt upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_SILT',                   &
          & fill_array=hydro_init_vars%fr_silt)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! fraction of clay upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_CLAY',                   &
          & fill_array=hydro_init_vars%fr_clay)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! fraction of sand deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_SAND',               &
          & fill_array=hydro_init_vars%fr_sand_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! fraction of silt deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_SILT',               &
          & fill_array=hydro_init_vars%fr_silt_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! fraction of clay deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_CLAY',               &
          & fill_array=hydro_init_vars%fr_clay_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! calculate volumetric soil porosity of mineral soil
        CALL soil_init_from_texture( &
           porosity_sand, porosity_silt, porosity_clay, porosity_oc,                                                &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%vol_porosity(:,:))

        ! calculate saturated hydraulic conductivity of mineral soil
        CALL soil_init_from_texture( &
           hyd_cond_sand, hyd_cond_silt, hyd_cond_clay, hyd_cond_oc,                                                &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%hyd_cond_sat(:,:))
        hydro_init_vars%hyd_cond_sat(:,:) = hydro_init_vars%hyd_cond_sat(:,:) * hyd_cond_sat_profile

        ! calculate volumetric field capacity [m/m] of mineral soil
        CALL soil_init_from_texture( &
           field_cap_sand, field_cap_silt, field_cap_clay, field_cap_oc,                                            &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%vol_field_cap(:,:))

        ! calculate plant wilting point of mineral soil
        CALL soil_init_from_texture( &
           wilt_sand, wilt_silt, wilt_clay, wilt_oc,                                                                &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%vol_p_wilt(:,:))

        ! calculate residual soil water content of mineral soil
        CALL soil_init_from_texture( &
           wres_sand, wres_silt, wres_clay, wres_oc,                                                                &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%vol_wres(:,:))

        ! calculate matric potential of mineral soil
        CALL soil_init_from_texture( &
           matric_pot_sand, matric_pot_silt, matric_pot_clay, matric_pot_oc,                                        &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%matric_pot(:,:))

        ! calculate pore_size_index of mineral soil
        CALL soil_init_from_texture( &
           pore_size_index_sand, pore_size_index_silt, pore_size_index_clay, pore_size_index_oc,                    &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%pore_size_index(:,:))

        ! calculate exponent b in Clapp & Hornberger of mineral soil
        CALL soil_init_from_texture( &
           bclapp_sand, bclapp_silt, bclapp_clay, bclapp_oc,                                                        &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%bclapp(:,:))

      END IF ! l_soil_texture

      ! Read soil organic carbon fractions for each soil layer
      IF (dsl4jsb_Config(HYDRO_)%l_socmap) THEN
        DO isoil=1,nsoil
          ptr_3D => input_file%Read_2d_extdim( &
            & variable_name='fract_org_sl',     &
            & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
          hydro_init_vars%fract_org_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
        END DO
      END IF

      ! Root depth [m]
      ! TODO root depth should become a function of PFTs eventually, instead of being a soil property
      ptr_2D => input_file%Read_2d(   &
        & variable_name='root_depth', &
        & fill_array = hydro_init_vars%root_depth)
      hydro_init_vars%root_depth = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
      ! Make sure there is root depth equivalent to 0.2 m weq in all cells
      !   that contain field capacity data but no root depth
      WHERE (hydro_init_vars%vol_field_cap(:,:) > 1.0e-10_wp)
        hydro_init_vars%root_depth(:,:) = MAX(hydro_init_vars%root_depth(:,:), &
          & 0.2_wp / hydro_init_vars%vol_field_cap(:,:))
      ELSEWHERE
        hydro_init_vars%root_depth(:,:) = 0._wp
      ENDWHERE
      ! Constrain root depth because in the bc file it may be larger then soil depth
      hydro_init_vars%root_depth = MIN(hydro_init_vars%root_depth, hydro_init_vars%soil_depth)

      CALL input_file%Close()

    END IF

    IF (tile%contains_soil) THEN

      IF (.NOT. model%config%init_from_ifs .OR. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
        input_file = jsb_netcdf_open_input(dsl4jsb_Config(HYDRO_)%ic_filename, model%grid_id)
      END IF

      ! Initialize hydrology variables from ic file

      IF (model%config%init_from_ifs) THEN

        ptr_3D => model%config%ifs_input_file%Read_2d_time( &
          & variable_name='W_SNOW', start_time_step=1, end_time_step=1)
        hydro_init_vars%ifs_weq_snow(:,:) = ptr_3D(:,:,1)

      ELSE

        ! Initial snow depth
        IF (input_file%Has_var('snow')) THEN
          ptr_2D => input_file%Read_2d( &
            & variable_name='snow',     &
            & fill_array = hydro_init_vars%weq_snow_soil)
          hydro_init_vars%weq_snow_soil = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
        ELSE
          hydro_init_vars%weq_snow_soil = 0._wp
          CALL message(TRIM(routine), 'Initializing snow to zero')
        END IF

      END IF

      IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN

        DO isoil=1,ifs_nsoil
          ptr_2D => model%config%ifs_input_file%Read_2d_1lev_1time( &
            & variable_name='SMIL'//int2string(isoil))
          hydro_init_vars%ifs_smi_sl(:,isoil,:) = ptr_2D(:,:)
        END DO

      ELSE

        DO isoil=1,nsoil
          ptr_3D => input_file%Read_2d_extdim( &
            & variable_name='layer_moist',     &
            & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
          hydro_init_vars%wtr_soil_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
        END DO

      END IF

      IF (.NOT. model%config%init_from_ifs .OR. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
        CALL input_file%Close()
      END IF

    END IF

  END SUBROUTINE hydro_read_init_vars

  SUBROUTINE hydro_init_bc(tile)

    USE mo_util,                  ONLY: soil_depth_to_layers_2d
    USE mo_hydro_process,         ONLY: calc_orographic_features

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    TYPE(t_jsb_model), POINTER        :: model
    TYPE(t_jsb_grid),  POINTER        :: grid
    TYPE(t_jsb_vgrid), POINTER        :: soil_w
    INTEGER                           :: i, nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_bc'

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Setting  hydrology boundary conditions for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)

    grid   => Get_grid(model%grid_id)
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    ! Get hydrology memory of the tile
    dsl4jsb_Get_memory(HYDRO_)

    !----------

    ! Orography ...
    IF (jsbach_runs_standalone()) dsl4jsb_var2D_onDomain(HYDRO_,elevation) = hydro_init_vars%elevation

    dsl4jsb_var2D_onDomain(HYDRO_,oro_stddev) = hydro_init_vars%oro_stddev

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_,elevation)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_,oro_stddev))

    ! Compute steepness parameter for surface runoff
    DO i = 1, grid%Get_nblks()
      CALL calc_orographic_features(                    &
        ! in
        & grid%Get_nproma(),                            &
        & grid%nlat_g,                                  & !< (Effective) number of latitudes for steepness parameter
        & dsl4jsb_var_ptr(HYDRO_, oro_stddev ) (:,i),   &
        ! out
        & dsl4jsb_var_ptr(HYDRO_, steepness  ) (:,i)    &
        & )
    END DO

    IF (.NOT. tile%is_lake .AND. .NOT. tile%is_glacier) THEN
      ! Maximum depth and extent of surface water ponds
      IF (dsl4jsb_Config(HYDRO_)%l_ponds) THEN
        ! Note, surface depression depth is given as negative value in the data, but represent
        !   a storage on top of the soil --> needs to be converted to positive values!
        dsl4jsb_var_ptr(HYDRO_, fract_pond_max) = hydro_init_vars%fract_pond_max(:,:)
        dsl4jsb_var_ptr(HYDRO_, weq_pond_max)   = hydro_init_vars%depth_pond_max(:,:) * (-1._wp)
      ELSE
        ! Note, that ponds are implicitely disabled by setting the maximum allowed pond fraction
        !   and depth to zero. Thus, we can avoid a large number of l_ponds conditions for surface
        !   and soil hydrology processes.
        dsl4jsb_var_ptr(HYDRO_, fract_pond_max) = 0._wp
        dsl4jsb_var_ptr(HYDRO_, weq_pond_max)   = 0._wp
      END IF
    END IF

    IF (tile%contains_soil) THEN

      ! Initialize hydrology variables from bc file

      ! Total soil depth until bedrock
      dsl4jsb_var2D_onDomain(HYDRO_, soil_depth) = hydro_init_vars%soil_depth

      ! Compute actual soil layer depths based on the fixed layer thicknesses from input file and
      ! the bedrock depth just read in.
      dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) = soil_depth_to_layers_2d( &
        & dsl4jsb_var2D_onDomain(HYDRO_, soil_depth),                          & ! Total soil depth until bedrock (from textures)
        & soil_w%dz(:))                                                          ! Soil layer thicknesses

      ! Volumetric soil porosity
      dsl4jsb_var_ptr(HYDRO_, vol_porosity) = hydro_init_vars%vol_porosity

      ! Saturated hydraulic conductivity
      dsl4jsb_var_ptr(HYDRO_, hyd_cond_sat) = hydro_init_vars%hyd_cond_sat

      ! Matric potential
      dsl4jsb_var_ptr(HYDRO_, matric_pot) = hydro_init_vars%matric_pot

      ! Exponent B in Clapp and Hornberger
      dsl4jsb_var_ptr(HYDRO_, bclapp) = hydro_init_vars%bclapp

      ! Pore size distribution index
      dsl4jsb_var_ptr(HYDRO_, pore_size_index) = hydro_init_vars%pore_size_index

      ! Volumetric field capacity [m/m]
      dsl4jsb_var_ptr(HYDRO_, vol_field_cap) = hydro_init_vars%vol_field_cap

      ! Get water content at field capacity by multiplying volumetric field capacity with soil layer thicknesses
      ! Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,i,:) = hydro_init_vars%vol_field_cap(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Volumetric permanent wilting point
      ! @todo: Do we need this for bare soil tiles?
      dsl4jsb_var_ptr(HYDRO_, vol_p_wilt) = hydro_init_vars%vol_p_wilt

      ! Get permanent wilting point by multiplying with soil layer thicknesses
      ! Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl) (:,i,:) = hydro_init_vars%vol_p_wilt(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Volumetric residual water contant
      dsl4jsb_var_ptr(HYDRO_, vol_wres) = hydro_init_vars%vol_wres

      ! Get absolute residual water content by multiplying with soil layer thicknesses
      ! Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_res_sl) (:,i,:) = hydro_init_vars%vol_wres(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Get water holding capacity (water content at saturation) by multiplying volumetric porosity with soil layer thicknesses
      ! Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_sat_sl) (:,i,:) = dsl4jsb_var2D_onDomain(HYDRO_,vol_porosity) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Soil organic carbon fraction from map
      IF (dsl4jsb_Config(HYDRO_)%l_socmap) THEN
        dsl4jsb_var_ptr(HYDRO_, fract_org_sl) = hydro_init_vars%fract_org_sl
      END IF

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_sat_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_pwp_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_fc_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_res_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_field_cap)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_p_wilt)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_wres)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, pore_size_index)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, bclapp)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, matric_pot)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, hyd_cond_sat)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_porosity)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, fract_org_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, fract_pond_max)), &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, weq_pond_max)), &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, soil_depth))

    END IF

    IF (tile%contains_vegetation) THEN

      dsl4jsb_var_ptr(HYDRO_, root_depth) = hydro_init_vars%root_depth
      dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl) = soil_depth_to_layers_2d( &
        & dsl4jsb_var2D_onDomain(HYDRO_, root_depth),                          & ! Total rooting depth
        & soil_w%dz(:))    ! Soil layer thicknesses

      ! Maximum root zone soil moisture
      ! Note: weq_rootzone_max corresponds to water content at field capacity.
      dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) = &
        & hydro_init_vars%vol_field_cap * hydro_init_vars%root_depth
      ! Constrain rootzone capacity if upper limit is given to soil water capacity
      IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
        dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) = MIN( &
          & dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max), dsl4jsb_Config(HYDRO_)%w_soil_limit)
      END IF

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, root_depth)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max))
    END IF

  END SUBROUTINE hydro_init_bc


  SUBROUTINE hydro_init_ic(tile)

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_hydro_util,          ONLY: get_amount_in_rootzone
    USE mo_util,                ONLY: ifs2soil

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)

    TYPE(t_jsb_model), POINTER              :: model
    TYPE(t_jsb_vgrid), POINTER              :: soil_w
    INTEGER                                 :: i
    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_ic'

    !-----------------------------------------------------------------------
    ! CONTENT

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing hydrology memory for tile '//TRIM(tile%name))

    model  => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')

    ! Get hydrology config
    dsl4jsb_Get_config(HYDRO_)

    ! Get memory of the tile
    dsl4jsb_Get_memory(HYDRO_)
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    !----------

    IF (tile%contains_soil) THEN

      ! Initial snow depth
      IF (model%config%init_from_ifs) THEN
        dsl4jsb_var_ptr(HYDRO_, weq_snow_soil)(:,:) = hydro_init_vars%ifs_weq_snow(:,:)
        dsl4jsb_var2D_onDomain(HYDRO_, weq_snow) = dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)
      ELSE
        dsl4jsb_var_ptr(HYDRO_, weq_snow_soil) = hydro_init_vars%weq_snow_soil
        dsl4jsb_var2D_onDomain(HYDRO_, weq_snow) = dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)
      END IF

      ! Initial soil moisture
      IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN

        CALL ifs2soil(hydro_init_vars%ifs_smi_sl(:,:,:), ifs_soil_depth(:),        &
          &           dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:), soil_w%ubounds(:))

        ! Convert soil moisture index (SMI) in IFS file to water content
        WHERE (dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) < -999._wp)
          ! Fill missing values
          dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) = &
            &  0.5_wp * (dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:) + dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,:,:))
        ELSE WHERE
          dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) = &
            & MIN(dsl4jsb_var_ptr(HYDRO_,wtr_soil_sat_sl)(:,:,:), &
            &     MAX(0._wp, &
            &         dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) &
            &         * (dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,:,:) - dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:)) &
            &         + dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:) &
            &        ) &
            &    )
        END WHERE

      ELSE

        dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl) (:,:,:) = hydro_init_vars%wtr_soil_sl(:,:,:)

      END IF

      ! Initial soil water overflow (used in terra planet experiments)
      dsl4jsb_var2D_onDomain(HYDRO_, tpe_overflow) = 0._wp

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, weq_snow))      &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl))   &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, tpe_overflow))

    END IF

    IF (tile%contains_soil) THEN
      ! Initialize organic soil layer fractions and modify saturated, field capacity and perm. wilting point capacities accordingly
      ! Note: at the moment, fract_org_sl doesn't change over time, but it could
      IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
        IF (.NOT. dsl4jsb_Config(HYDRO_)%l_socmap) THEN
          ! Compute soil organic carbon fraction from forest fraction instead of using map from bc file
          dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl) = 0._wp
          IF (tile%is_vegetation) THEN
            ! For general VEG tile
            dsl4jsb_var_ptr(HYDRO_, fract_org_sl)(:,1,:) = dsl4jsb_var2D_onDomain(PHENO_, fract_forest)
            IF (tile%lcts(1)%lib_id /= 0)  THEN ! only if the present tile is a PFT (otherwise lctlib not defined)
              ! For PFT tiles
              ! Note: in current setup hydrology does not run on PFTs
              IF (dsl4jsb_Lctlib_param(ForestFlag)) THEN
                dsl4jsb_var_ptr(HYDRO_, fract_org_sl)(:,1,:) = dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max)
              END IF
            END IF
          END IF
        END IF

        WHERE (dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl) <= 1.0e-10_wp)
          dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl) = 0._wp
        ENDWHERE
        !$ACC UPDATE DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)) ASYNC(1)
      END IF

      ! update soil properties with organic layer fractions
      CALL init_soil_properties(tile)

      ! Make sure that soil water content is not larger than the maximum soil water capacity.
      ! Note: the volumetric porosity represents the maximum water capacity of the specific soil layer. Reduce
      ! soil water and ice if this limit is exceeded.
      IF (ANY(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)     + dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl) > &
        &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl))) THEN
        ! First limit soil water to maximum soil water capacity
        dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl) = MAX(0._wp, &
          & MIN(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl),     &
          &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)))
        ! Then limit soil ice to the remaining capacity (without soil water)
        dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl) = MAX(0._wp, &
          & MIN(dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl),     &
          &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) &
          &     - dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)))
      END IF

      ! The following can change later depending on the organic layers and soil ice content.
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sat_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl)  * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl) =  MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_field_cap_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_pwp_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_p_wilt_sl)    * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_res_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_wres_sl)      * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl))      &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))      &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sat_sl))  &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl))   &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_pwp_sl))  &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_res_sl))
    END IF

    IF (tile%contains_vegetation) THEN   ! Assuming there is also soil and roots

      ! Re-compute maximum root zone moisture based on soil parameters modified by organic fractions
      ! Note: All plant related computations use weq_rootzone_max reduced by ice and supercooled water content,
      !       which represents the unfrozen part of the soil.
      !       This assumes that plants retract/extents their roots immediately if the soil freezes/thaws.
      DO i=1,SIZE(dsl4jsb_var_ptr(HYDRO_, weq_rootzone_max), 2)
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, wtr_soil_fc_sl)   (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, weq_rootzone_max) (:,  i)  )
      END DO
      !$ACC UPDATE HOST(dsl4jsb_var_ptr(HYDRO_, weq_rootzone_max)) ASYNC(1)
      !$ACC WAIT(1)
      ! Set maximum root zone moisture to 0.2 for glacier and ocean areas (where weq_rootzone_max == 0)
      dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) = MERGE(       &
          & dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max), 0.2_wp, &
          & dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) > 0._wp )
      IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
        dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) = MIN( &
        & dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max), dsl4jsb_Config(HYDRO_)%w_soil_limit)
      END IF
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, weq_rootzone_max)) ASYNC(1)

      ! Initialization of relative root zone soil moisture (only water, not ice).
      !   This is needed for the LoGro-P phenology, as phenology is updated before hydrology.

      ! wtr_rootzone calculated here is only used to calculate relative soil moisture for phenology.
      ! It is re-calculated in update_soil_hydrology.
      DO i=1,SIZE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone), 2)
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, wtr_soil_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)   (:,  i)  )
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, ice_soil_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, ice_rootzone)   (:,  i)  )
      END DO
      !$ACC UPDATE ASYNC(1) &
      !$ACC HOST(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)) &
      !$ACC HOST(dsl4jsb_var_ptr(HYDRO_, ice_rootzone))
      !$ACC WAIT(1)

      ! Make sure wtr_rootzone is smaller than weq_rootzone_max - ice_rootzone.
      ! Note: this corresponds to the maximum possible amount of liquid water in the root zone.
      WHERE (dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) > &
        & dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) = MAX(0._wp, &
        & dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
      END WHERE
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)) ASYNC(1)

      ! Initialization of wtr_rootzone_rel
      !   This is re-calculated every time step in update_water_stress.
      WHERE (dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone) > 0._wp)
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone_rel) = &
          &  dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) / &
          & (dsl4jsb_var2D_onDomain(HYDRO_, weq_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
      ELSE WHERE
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone_rel) = 0._wp
      END WHERE

      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone_rel)) ASYNC(1)

    END IF ! Tile with vegetation

  END SUBROUTINE hydro_init_ic


  SUBROUTINE init_soil_properties(tile)
    ! This routine is needed because soil snow and energy (SSE) processes are calculated before
    ! the soil properties are calculated the first time within the hydrology in update_soil_properties.

    USE mo_hydro_constants, ONLY: vol_porosity_org_top, hyd_cond_sat_org_top, bclapp_org_top, matric_pot_org_top, &
      & pore_size_index_org_top, vol_field_cap_org_top, vol_p_wilt_org_top, vol_wres_org_top, &
      & vol_porosity_org_below, hyd_cond_sat_org_below, bclapp_org_below, matric_pot_org_below, &
      & pore_size_index_org_below, vol_field_cap_org_below, vol_p_wilt_org_below, vol_wres_org_below

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & hyd_cond_sat,          &
      & vol_porosity,          &
      & bclapp,                &
      & matric_pot,            &
      & pore_size_index,       &
      & vol_field_cap,         &
      & vol_p_wilt,            &
      & vol_wres
    dsl4jsb_Real3D_onDomain :: &
      & fract_org_sl,          &
      & hyd_cond_sat_sl,       &
      & vol_porosity_sl,       &
      & bclapp_sl,             &
      & matric_pot_sl,         &
      & pore_size_index_sl,    &
      & vol_field_cap_sl,      &
      & vol_p_wilt_sl,         &
      & vol_wres_sl

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    REAL(wp) ::  hyd_cond_sat_org, vol_porosity_org, bclapp_org, matric_pot_org, pore_size_index_org, &
      & vol_field_cap_org, vol_p_wilt_org, vol_wres_org

    INTEGER  :: nsoil, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (tile%is_lake) RETURN

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, hyd_cond_sat)       ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_porosity)       ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, bclapp)             ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, matric_pot)         ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, pore_size_index)    ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_field_cap)      ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_p_wilt)         ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_wres)           ! in
    IF (dsl4jsb_Config(HYDRO_)%l_organic)   &
     &   dsl4jsb_Get_var3D_onDomain(HYDRO_, fract_org_sl)  ! in
    dsl4jsb_Get_var3D_onDomain(HYDRO_, hyd_cond_sat_sl)    ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)    ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, bclapp_sl)          ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, matric_pot_sl)      ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, pore_size_index_sl) ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_field_cap_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_p_wilt_sl)      ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_wres_sl)        ! out

    hyd_cond_sat_sl   (:,:,:) = SPREAD(hyd_cond_sat   (:,:), DIM=2, ncopies=nsoil)
    vol_porosity_sl   (:,:,:) = SPREAD(vol_porosity   (:,:), DIM=2, ncopies=nsoil)
    bclapp_sl         (:,:,:) = SPREAD(bclapp         (:,:), DIM=2, ncopies=nsoil)
    matric_pot_sl     (:,:,:) = SPREAD(matric_pot     (:,:), DIM=2, ncopies=nsoil)
    pore_size_index_sl(:,:,:) = SPREAD(pore_size_index(:,:), DIM=2, ncopies=nsoil)
    vol_field_cap_sl  (:,:,:) = SPREAD(vol_field_cap  (:,:), DIM=2, ncopies=nsoil)
    vol_p_wilt_sl     (:,:,:) = SPREAD(vol_p_wilt     (:,:), DIM=2, ncopies=nsoil)
    vol_wres_sl       (:,:,:) = SPREAD(vol_wres       (:,:), DIM=2, ncopies=nsoil)
    IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
      ! Update soil parameters with organic fractions for deep and top soil layers
      ! Attention: hyd_cond_sat_sl is calculated differently (correct) in update_soil_properties
      DO i=1,nsoil
        IF (i == 1) THEN
          hyd_cond_sat_org    = hyd_cond_sat_org_top
          vol_porosity_org    = vol_porosity_org_top
          bclapp_org          = bclapp_org_top
          matric_pot_org      = matric_pot_org_top
          pore_size_index_org = pore_size_index_org_top
          vol_field_cap_org   = vol_field_cap_org_top
          vol_p_wilt_org      = vol_p_wilt_org_top
          vol_wres_org        = vol_wres_org_top
        ELSE
          hyd_cond_sat_org    = hyd_cond_sat_org_below
          vol_porosity_org    = vol_porosity_org_below
          bclapp_org          = bclapp_org_below
          matric_pot_org      = matric_pot_org_below
          pore_size_index_org = pore_size_index_org_below
          vol_field_cap_org   = vol_field_cap_org_below
          vol_p_wilt_org      = vol_p_wilt_org_below
          vol_wres_org        = vol_wres_org_below
        END IF

        hyd_cond_sat_sl   (:,i,:) = (1 - fract_org_sl(:,i,:)) * hyd_cond_sat_sl   (:,i,:) &
          & + fract_org_sl(:,i,:) * hyd_cond_sat_org
        vol_porosity_sl   (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_porosity_sl   (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_porosity_org
        bclapp_sl         (:,i,:) = (1 - fract_org_sl(:,i,:)) * bclapp_sl         (:,i,:) &
          & + fract_org_sl(:,i,:) * bclapp_org
        matric_pot_sl     (:,i,:) = (1 - fract_org_sl(:,i,:)) * matric_pot_sl     (:,i,:) &
          & + fract_org_sl(:,i,:) * matric_pot_org
        pore_size_index_sl(:,i,:) = (1 - fract_org_sl(:,i,:)) * pore_size_index_sl(:,i,:) &
          & + fract_org_sl(:,i,:) * pore_size_index_org
        vol_field_cap_sl  (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_field_cap_sl  (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_field_cap_org
        vol_p_wilt_sl     (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_p_wilt_sl     (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_p_wilt_org
        vol_wres_sl       (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_wres_sl       (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_wres_org
      END DO
    END IF

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl)         &
    !$ACC   DEVICE(matric_pot_sl, pore_size_index_sl, vol_field_cap_sl) &
    !$ACC   DEVICE(vol_p_wilt_sl, vol_wres_sl)

  END SUBROUTINE init_soil_properties


  !> Sanitize hydro variables.
  !!
  !! Loading a model state from a file with lower than double precision can lead to model aborts because the state
  !! violates model constraints due to rounding errors in soil water content and root depth. This routine sanitizes
  !! the current state by limiting soil water and ice to the field capacity, and by recomputing the root depth.
  SUBROUTINE hydro_sanitize_state(tile)

    USE mo_util, ONLY: soil_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_sanitize_state'

    !> If water content exceeds field capacity by more than this factor prior to sanitization, a warning is issued.
    REAL(wp), PARAMETER :: warn_level_weq_soil_excess = 1.05_wp
    REAL(wp), PARAMETER :: zfcmin = 1.e-10_wp              ! minimum field capacity (or rather epsilon?)

    INTEGER :: maxloc_weq_soil_excess(3)
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Real2D_onDomain :: &
      & root_depth,            &
      & soil_depth

    dsl4jsb_Real3D_onDomain :: &
      & root_depth_sl,         &
      & soil_depth_sl,         &
      & vol_porosity_sl,       &
      & wtr_soil_sat_sl,       &
      & wtr_soil_sl,           &
      & ice_soil_sl

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, root_depth)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, soil_depth)

    dsl4jsb_Get_var3D_onDomain(HYDRO_, root_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sat_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, ice_soil_sl)

    IF (tile%contains_soil) THEN

      ! Issue a warning if current soil water anywherer exceeds saturation capacity by more than a factor of
      ! `warn_level_weq_soil_excess`.
      IF (ANY(wtr_soil_sl(:,:,:) > warn_level_weq_soil_excess * wtr_soil_sat_sl(:,:,:))) THEN
        maxloc_weq_soil_excess(:) = MAXLOC(wtr_soil_sl(:,:,:) - wtr_soil_sat_sl(:,:,:))
        WRITE (message_text,*) 'liquid water content above saturation capacity at ', maxloc_weq_soil_excess(:), &
          & ': wtr_soil_sl = ',                                                                                 &
          & wtr_soil_sl    (maxloc_weq_soil_excess(1), maxloc_weq_soil_excess(2), maxloc_weq_soil_excess(3)),   &
          & ', wtr_soil_sat_sl = ',                                                                             &
          & wtr_soil_sat_sl(maxloc_weq_soil_excess(1), maxloc_weq_soil_excess(2), maxloc_weq_soil_excess(3))
        CALL warning(routine, message_text)
      END IF
      ! Limit soil water and soil ice to saturation amounts.
      IF (ANY((wtr_soil_sl(:,:,:) + ice_soil_sl(:,:,:))                &
        &  > MAX(vol_porosity_sl (:,:,:) * soil_depth_sl(:,:,:), 0._wp))) THEN
        ! Limit liquid soil water to saturation amount
        wtr_soil_sl(:,:,:) = MIN(wtr_soil_sl(:,:,:), wtr_soil_sat_sl(:,:,:))

        ! Limit ice content such that soil water + soil ice <= saturation.
        ice_soil_sl(:,:,:) = MAX(0._wp, MIN(ice_soil_sl(:,:,:), &
          &                             MAX(vol_porosity_sl (:,:,:) * soil_depth_sl(:,:,:), 0._wp) - wtr_soil_sl(:,:,:)))

        !$ACC UPDATE DEVICE(wtr_soil_sl, ice_soil_sl) ASYNC(1)
      END IF
    END IF

    IF (tile%contains_vegetation) THEN
      ! Due to rounding, the root depth in a layer might differ from the soil depth even though the layer is contained
      ! in the root zone. The might not add up to the total root depth for the same reason. Fix by recomputing.
      soil_w => Get_vgrid('soil_depth_water')
      root_depth(:,:) = MIN(root_depth(:,:), soil_depth(:,:))
      root_depth_sl(:,:,:) = soil_depth_to_layers_2d(root_depth(:,:), soil_w%dz(:))

      !$ACC UPDATE DEVICE(root_depth, root_depth_sl) ASYNC(1)
    END IF

  END SUBROUTINE hydro_sanitize_state

#endif
END MODULE mo_hydro_init
