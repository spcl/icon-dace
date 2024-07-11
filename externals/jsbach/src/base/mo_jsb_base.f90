!> Contains main JSBACH structure and methods.
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
MODULE mo_jsb_base
#ifndef __NO_JSBACH__

  USE mo_jsb_class,   ONLY: jsbach
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message, finish
  USE mo_util,        ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_setup_models, jsbach_setup_tiles

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Setup the ICON-Land models
  !>
  !>   Routine is called by both:
  !>      ICON-Land standalone: jsbach: src/drivers/mo_jsbach_model.f90: jsbach_model
  !>      host model (ICON):    icon: src/drivers/mo_atmo_model.f90: construct_atmo_model
  !>
  SUBROUTINE jsbach_setup_models(master_namelist_filename)

    USE mo_jsb_control,         ONLY: jsb_models_nml, get_no_of_models, init_jsb_master_control, &
      &                               jsbach_runs_standalone, jsbach_is_restarted
    USE mo_jsb_parallel,        ONLY: init_parallel_iface, Get_omp_no_of_threads
    USE mo_jsb_model_class,     ONLY: t_jsb_model, new_model, MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_subset,          ONLY: jsbach_subsets
    USE mo_jsb_config_class,    ONLY: new_model_config
    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process
    USE mo_jsb_time,            ONLY: init_time
    USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,            ONLY: Register_vgrid, Get_vgrid
    USE mo_jsb_io,              ONLY: ZAXIS_SURFACE
    USE mo_pheno_parameters,    ONLY: init_pheno_parameters
#ifndef __NO_QUINCY__
    USE mo_q_rad_parameters,    ONLY: init_q_rad_parameters
    USE mo_atmland_constants,   ONLY: init_atmland_constants
    USE mo_q_assimi_parameters, ONLY: init_q_assimi_parameters
    USE mo_veg_constants,       ONLY: init_veg_constants
    USE mo_sb_constants,        ONLY: init_sb_constants
    USE mo_q_pheno_parameters,  ONLY: init_q_pheno_parameters
    USE mo_spq_constants,       ONLY: init_spq_constants
    USE mo_lnd_bgcm_idx,        ONLY: set_bgcm_element_idx
#endif

    dsl4jsb_Use_processes SEB_, ASSIMI_, PHENO_, TURB_, SSE_, HYDRO_, RAD_, &
      &                   CARBON_, DISTURB_, FUEL_, PPLCC_, WLCC_, FLCC_, FAGE_

    dsl4jsb_Use_config(SSE_)
    dsl4jsb_Use_config(HYDRO_)
    dsl4jsb_Use_config(RAD_)
    dsl4jsb_Use_config(CARBON_)
    dsl4jsb_Use_config(DISTURB_)
    dsl4jsb_Use_config(FUEL_)
    dsl4jsb_Use_config(FAGE_)
    dsl4jsb_Use_config(FLCC_)
    dsl4jsb_Use_config(WLCC_)
    dsl4jsb_Use_config(PPLCC_)

    CHARACTER(len=*), INTENT(in) :: master_namelist_filename

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_config(DISTURB_)
    dsl4jsb_Def_config(FUEL_)
    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_config(FLCC_)
    dsl4jsb_Def_config(WLCC_)
    dsl4jsb_Def_config(PPLCC_)

    TYPE(t_jsb_model), POINTER :: model => NULL()
    TYPE(t_jsb_vgrid), POINTER :: soil_e, soil_w, snow_e

    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid

    INTEGER :: no_of_models, model_idx, iproc
    INTEGER :: nsoil_e, nsoil_w, nsnow_e
    INTEGER :: model_scheme

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:jsbach_setup_models'

    CALL message(TRIM(routine), 'starting basic setup')

    CALL init_parallel_iface

    CALL init_jsb_master_control(master_namelist_filename)

    !IF (jsbach_runs_standalone()) CALL init_mpi_communicators !TODO-check-with-RMS: remove: done in the driver!

    no_of_models = get_no_of_models()
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models = no_of_models

    ALLOCATE(jsbach_subsets(no_of_models))

    ! Create models according to master namelist
    DO model_idx=1,no_of_models
      NULLIFY(model)
      model => new_model(jsb_models_nml(model_idx)%model_id,                &
        &                  jsb_models_nml(model_idx)%model_name,              &
        &                  jsb_models_nml(model_idx)%model_shortname,         &
        &                  jsb_models_nml(model_idx)%model_description,       &
        &                  jsb_models_nml(model_idx)%model_namelist_filename)
      jsbach%models(jsb_models_nml(model_idx)%model_id)%m => model
      ALLOCATE(jsbach_subsets(model_idx)%sub(Get_omp_no_of_threads()))
    END DO

    ! Configure models according to model namelist(s)
    DO model_idx=1,no_of_models
       NULLIFY(model)
       model => jsbach%models(model_idx)%m
       model%config => new_model_config(model%namelist_filename)

       ! Require that only one model_scheme is used!
       IF (model_idx == 1) THEN
         model_scheme = model%config%model_scheme
       ELSE IF (model%config%model_scheme /= model_scheme) THEN
          CALL finish(routine, 'Only one model_scheme allowed!')
       END IF

#ifndef __NO_QUINCY__
       !> set the bgcm index var used in the scientific routines
       !>   it depends on the model config and is identical across models
       !>   assuming the namelist options are identical across models
       !>
       IF (model_idx == 1) THEN
         CALL set_bgcm_element_idx(model%config%elements_index_map)
       END IF
#endif

       !> parameter init
       !>
       IF (model_idx == 1) THEN
         ! assuming the shared/constants/ modules are available here already (impl, math, physical)
         SELECT CASE (model_scheme)
         CASE (MODEL_QUINCY)
           CALL message(TRIM(routine), 'starting init of QUINCY process parameters')
#ifndef __NO_QUINCY__
           CALL init_q_rad_parameters
           CALL init_q_assimi_parameters
           CALL init_veg_constants
           CALL init_spq_constants
           CALL init_sb_constants
           CALL init_q_pheno_parameters
           CALL init_atmland_constants
#else
           CALL finish(routine, 'Model has been compiled without support for QUINCY, use --enable-quincy')
#endif
         CASE (MODEL_JSBACH)
           CALL message(TRIM(routine), 'starting init of JSBACH process parameters')
           CALL init_pheno_parameters
         END SELECT
       END IF

       CALL message(TRIM(routine), 'New jsbach model: '//TRIM(model%shortname))
       CALL message('     ', '... ID: '//TRIM(int2string(model%id)))
       CALL message('     ', '... Usecase: '//TRIM(model%config%usecase))
       CALL message('     ', '... Using configuration from namelist: '//TRIM(model%namelist_filename))
       CALL message('     ', '... Using tile fractions from: '//TRIM(model%config%fract_filename))

       IF (jsbach_runs_standalone()) THEN
         CALL init_time(TRIM(model%namelist_filename), TRIM(model%shortname), jsbach_is_restarted())
       ELSE
         CALL init_time()
       END IF

      CALL model%Set_mode(model%config%hsm_mode)

      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(surface)

      ALLOCATE(model%processes(max_no_of_processes))
      DO iproc=1,max_no_of_processes
        model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename))
      END DO

      ! Point from each associated process%config to the model_config
      DO iproc = 1,max_no_of_processes
        IF (ASSOCIATED(model%processes(iproc)%p)) THEN
          IF (ASSOCIATED(model%processes(iproc)%p%config)) THEN
            model%processes(iproc)%p%config%model_config => model%config
          ENDIF
        ENDIF
      ENDDO

      ! Configure processes: run 'Init_PROC_config'
      CALL model%Configure_processes()

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        dsl4jsb_Get_config(SSE_)
        dsl4jsb_Get_config(HYDRO_)
        dsl4jsb_Get_config(RAD_)
        dsl4jsb_Get_config(CARBON_)
        dsl4jsb_Get_config(FUEL_)
        dsl4jsb_Get_config(DISTURB_)
        dsl4jsb_Get_config(FAGE_)
        dsl4jsb_Get_config(WLCC_)
        dsl4jsb_Get_config(FLCC_)
        dsl4jsb_Get_config(PPLCC_)

        ! Make number of vertical layers available from the options type
        IF (model%Is_process_enabled(HYDRO_)) THEN
          soil_w => Get_vgrid('soil_depth_water')
          nsoil_w = soil_w%n_levels
        END IF
        IF (model%Is_process_enabled(SSE_)) THEN
          soil_e => Get_vgrid('soil_depth_energy')
          nsoil_e = soil_e%n_levels
          IF (dsl4jsb_Config(SSE_)%l_snow) THEN
            snow_e => Get_vgrid('snow_depth_energy')
            nsnow_e = snow_e%n_levels
          END IF
        END IF
        CALL model%Set_options(nsoil_e=nsoil_e, nsoil_w=nsoil_w, nsnow_e=nsnow_e)

        ! Crosscheck process configs
        IF (model%Is_process_enabled(HYDRO_) .AND. model%Is_process_enabled(SSE_) &
          & .AND. (nsoil_e /= nsoil_w)) THEN
          CALL finish(routine, &
            & 'The number of soil layers used for water and energy currently need to be identical:'  &
            & //' nsoil_w = '//TRIM(int2string(nsoil_w))//', nsoil_e = '//TRIM(int2string(nsoil_e)))
        END IF

      END SELECT

      IF (TRIM(model%config%tpe_scheme) /= '') THEN
        CALL message(TRIM(routine), 'Running terraplanet experiment with scheme '//TRIM(model%config%tpe_scheme))
        IF (model%config%use_lakes) THEN
          CALL message('    ', '... using lakes')
        ELSE
          CALL message('    ', '... not using lakes')
        END IF
        IF (model%config%use_glacier) THEN
          CALL message('    ', '... using glaciers')
        ELSE
          CALL message('    ', '... not using glaciers')
        END IF
      END IF

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
        ! There are several jsbach processes for which active is set to TRUE as default, however, these are not required for IQ:
        model%processes(SSE_)%p%config%active = .FALSE.
        model%processes(SEB_)%p%config%active = .FALSE.
        model%processes(RAD_)%p%config%active = .FALSE.
        model%processes(TURB_)%p%config%active = .FALSE.
        model%processes(ASSIMI_)%p%config%active = .FALSE.
        model%processes(PHENO_)%p%config%active = .FALSE.
      CASE (MODEL_JSBACH)
        IF (model%config%l_compat401) THEN
          CALL message(TRIM(routine), 'Using settings compatible to jsbach4.01:')
          CALL message('    ', '... using soil heat capacity from mineral soil map')
          CALL message('    ', '... using soil heat conductivity from mineral soil map')
          CALL message('    ', '... using fixed snow heat capacity and snow heat conductivity')
          CALL message('    ', '... not using multi-layer snow model')
          CALL message('    ', '... not using dynamical snow density and thickness')
          CALL message('    ', '... not using frozen ground')
          CALL message('    ', '... not using organic soil fractions')
          CALL message('    ', '... using canopy albedo from ini file')
          dsl4jsb_Config(SSE_)%l_heat_cap_dyn = .FALSE.
          dsl4jsb_Config(SSE_)%l_heat_cond_dyn = .FALSE.
          dsl4jsb_Config(SSE_)%l_heat_cap_map = .FALSE.
          dsl4jsb_Config(SSE_)%l_heat_cond_map = .FALSE.
          dsl4jsb_Config(SSE_)%l_snow = .FALSE.
          dsl4jsb_Config(SSE_)%l_dynsnow = .FALSE.
          dsl4jsb_Config(SSE_)%l_freeze = .FALSE.
          dsl4jsb_Config(HYDRO_)%l_organic = .FALSE.
          dsl4jsb_Config(RAD_)%use_alb_canopy = .TRUE.
        END IF

        IF (model%Do_fractions_change()) THEN
          IF (.NOT. model%Is_process_enabled(PPLCC_)) THEN
            CALL message(TRIM(routine), 'Activating PPLCC_ process since some process changes fractions')
          END IF
          dsl4jsb_Config(PPLCC_)%active = .TRUE.
        END IF

        IF (model%config%init_from_ifs) THEN
          CALL message(TRIM(routine), 'Initialize JSBACH from IFS analysis file '//TRIM(model%config%ifs_filename))
        END IF

        IF (.NOT. dsl4jsb_Config(SSE_)%l_snow) THEN
          ! Don't allow infiltration at sub-zero soil temperature for old snow scheme
          dsl4jsb_Config(HYDRO_)%l_infil_subzero = .FALSE.
          CALL message(TRIM(routine), 'Not allowing infiltration at sub-zero soil temperature (l_snow=false)')
        ELSE IF (.NOT. dsl4jsb_Config(HYDRO_)%l_infil_subzero) THEN
          ! Print message if infiltration at sub-zero temperature is switched off in namelist
          CALL message(TRIM(routine), 'Not allowing infiltration at sub-zero soil temperature')
        END IF

        IF (dsl4jsb_Config(RAD_)%use_alb_canopy .OR. TRIM(model%config%usecase) == 'jsbach_lite') THEN
          CALL message(TRIM(routine), 'Using canopy albedo from ini file '//TRIM(dsl4jsb_Config(RAD_)%bc_filename)//&
            &                         ' ... therefore no influence of vegetation on albedo')
        ELSE
          CALL message(TRIM(routine), 'Using canopy albedo from lctlib file')
        END IF

        IF (model%Is_process_enabled(FLCC_) .AND. .NOT. model%Is_process_enabled(DISTURB_)) THEN
          dsl4jsb_Config(DISTURB_)%active = .TRUE.
          CALL message(TRIM(routine), 'Switching on disturbance process since fire lcc process is active')
        END IF
        IF (model%Is_process_enabled(WLCC_) .AND. .NOT. model%Is_process_enabled(DISTURB_)) THEN
          dsl4jsb_Config(DISTURB_)%active = .TRUE.
          CALL message(TRIM(routine), 'Switching on disturbance process since wind lcc process is active')
        END IF
        IF (model%Is_process_enabled(DISTURB_) .AND. .NOT. model%Is_process_enabled(FUEL_)) THEN
          dsl4jsb_Config(FUEL_)%active = .TRUE.
          CALL message(TRIM(routine), 'Switching on fuel process since disturbance process is active')
        END IF
        IF (model%Is_process_enabled(DISTURB_) .AND. .NOT. model%Is_process_enabled(CARBON_)) THEN
          dsl4jsb_Config(CARBON_)%active = .TRUE.
          CALL message(TRIM(routine), 'Switching on carbon process since disturbance process is active')
        END IF

        IF(model%Is_process_enabled(FAGE_) &
            & .AND. .NOT. (TRIM(model%config%usecase) == 'jsbach_forest_age_classes')) THEN
          CALL finish(routine, &
            & 'The FAGE process should only be used with the jsbach_forest_age_classes usecase! Please check!')
        ENDIF

        IF (model%Is_process_enabled(FAGE_) .AND. model%Is_process_enabled(DISTURB_) &
            & .AND. .NOT. model%Is_process_enabled(FLCC_)) THEN
          CALL finish(routine, &
            & 'When FAGE is enabled and disturbances are active, flcc should be active, too! Please check!')
        END IF
        IF (model%Is_process_enabled(FAGE_) .AND. model%Is_process_enabled(DISTURB_) &
            & .AND. .NOT. model%Is_process_enabled(WLCC_)) THEN
          CALL finish(routine, &
            & 'When FAGE is enabled and disturbances are active, wlcc should be active, too! Please check!')
        END IF

      END SELECT

    END DO

    CALL message(TRIM(routine), 'done basic setup (models, process configs, QUINCY process constants init)')

  END SUBROUTINE jsbach_setup_models

  SUBROUTINE jsbach_setup_tiles(model_id)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_model_usecases, ONLY: init_usecase
    USE mo_jsb_lctlib_class,   ONLY: Read_lctlib

    INTEGER, INTENT(in) :: model_id

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(LEN=30)           :: lctlib_filename
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:jsbach_setup_tiles'

    model => Get_model(model_id)

    ! select lctlib file depending on the usecase
    ! @TODO should always link to the same lctlib filename and do the rest by scripting environment
    SELECT CASE (TRIM(model%config%usecase))
    ! quincy model with 11 PFT tiles incl. bare soil
    CASE ('quincy_11_pfts')
      lctlib_filename = 'lctlib_quincy_nlct14.def'
    ! jsbach_lite & jsbach_pfts
    CASE DEFAULT
      lctlib_filename = 'lctlib_nlct21.def'
    END SELECT

    ! import lctlib
    ! use model_scheme_char, instead of the ENUM model_scheme, to differentiate between model schemes
    !  because mo_jsb_model_class is not used in this function
    model%lctlib => Read_lctlib(TRIM(lctlib_filename), &
      &                         TRIM(model%config%model_scheme_char), &
      &                         TRIM(model%config%usecase))
    CALL message(TRIM(routine), 'Imported lctlib from '//TRIM(lctlib_filename)// &
      &                         ' using the model scheme '//model%config%model_scheme_char)

    ! initialize usecase if defined via namelist
    IF (model%config%usecase /= '') THEN
      ! Shortcut model configuration through usecase
      CALL init_usecase(model)
    ELSE
      ! Explicit model configuration through namelist
      ! TODO
    END IF

  END SUBROUTINE jsbach_setup_tiles

#endif
END MODULE mo_jsb_base
