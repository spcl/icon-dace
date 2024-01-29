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
  USE mo_exception,   ONLY: message
  USE mo_util,        ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach
  PUBLIC :: jsbach_setup_models, jsbach_setup_tiles

CONTAINS

  SUBROUTINE jsbach_setup_models(master_namelist_filename)

    USE mo_jsb_control,         ONLY: jsb_models_nml, get_no_of_models, init_jsb_master_control, &
      &                               jsbach_runs_standalone, jsbach_is_restarted
    USE mo_jsb_parallel,        ONLY: init_parallel_iface, Get_omp_no_of_threads
    USE mo_jsb_model_class,     ONLY: t_jsb_model, new_model
    USE mo_jsb_subset,          ONLY: jsbach_subsets
    USE mo_jsb_config_class,    ONLY: new_model_config
    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process
    USE mo_jsb_time,            ONLY: init_time
    USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,            ONLY: Register_vgrid
    USE mo_jsb_io,              ONLY: ZAXIS_SURFACE
    USE mo_rad_constants,       ONLY: init_rad_constants
    USE mo_atmland_constants,   ONLY: init_atmland_constants
    USE mo_assimi_constants,    ONLY: init_assimi_constants
    USE mo_veg_constants,       ONLY: init_veg_constants
    ! USE mo_sse_constants,       ONLY: init_sse_constants
    ! USE mo_sb_constants,        ONLY: init_sb_constants
    USE mo_pheno_constants,     ONLY: init_pheno_constants

    dsl4jsb_Use_processes SSE_, HYDRO_, RAD_, CARBON_, DISTURB_, FUEL_, PPLCC_
    dsl4jsb_Use_config(SSE_)
    dsl4jsb_Use_config(HYDRO_)
    dsl4jsb_Use_config(RAD_)
    dsl4jsb_Use_config(CARBON_)
    dsl4jsb_Use_config(DISTURB_)
    dsl4jsb_Use_config(FUEL_)
    dsl4jsb_Use_config(PPLCC_)

    CHARACTER(len=*), INTENT(in) :: master_namelist_filename

    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_config(DISTURB_)
    dsl4jsb_Def_config(FUEL_)
    dsl4jsb_Def_config(PPLCC_)

    TYPE(t_jsb_model), POINTER :: model => NULL()

    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid

    INTEGER :: no_of_models, model_idx, iproc

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

      CALL model%Configure_processes()

      dsl4jsb_Get_config(SSE_)
      dsl4jsb_Get_config(HYDRO_)
      dsl4jsb_Get_config(RAD_)
      dsl4jsb_Get_config(CARBON_)
      dsl4jsb_Get_config(FUEL_)
      dsl4jsb_Get_config(DISTURB_)
      dsl4jsb_Get_config(PPLCC_)

      ! Crosscheck process configs
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
        IF (.NOT. dsl4jsb_Config(PPLCC_)%active) THEN
          CALL message(TRIM(routine), 'Activating PPLCC_ process since some process changes fractions')
        END IF
        dsl4jsb_Config(PPLCC_)%active = .TRUE.
      END IF

      IF (model%config%init_from_ifs) THEN
        CALL message(TRIM(routine), 'Initialize JSBACH from IFS analysis file '//TRIM(model%config%ifs_filename))
      END IF

      IF (dsl4jsb_Config(RAD_)%use_alb_canopy .OR. TRIM(model%config%usecase) == 'jsbach_lite') THEN
        CALL message(TRIM(routine), 'Using canopy albedo from ini file '//TRIM(dsl4jsb_Config(RAD_)%bc_filename)//&
          &                         ' ... therefore no influence of vegetation on albedo')
      ELSE
        CALL message(TRIM(routine), 'Using canopy albedo from lctlib file')
      END IF

      IF (dsl4jsb_Config(DISTURB_)%active .AND. .NOT. dsl4jsb_Config(FUEL_)%active) THEN
        dsl4jsb_Config(FUEL_)%active = .TRUE.
        CALL message(TRIM(routine), 'Switching on fuel process since disturbance process is active')
      END IF
      IF (dsl4jsb_Config(DISTURB_)%active .AND. .NOT. dsl4jsb_Config(CARBON_)%active) THEN
        dsl4jsb_Config(CARBON_)%active = .TRUE.
        CALL message(TRIM(routine), 'Switching on carbon process since disturbance process is active')
      END IF

    END DO

    !> quincy constants init
    !! assuming the shared/constants/ modules are available here already (impl, math, physical)
    IF (model%config%use_quincy) THEN
      CALL message(TRIM(routine), 'Starting init of QUINCY process constants')
      CALL init_rad_constants      
      CALL init_atmland_constants  
      CALL init_assimi_constants   
      CALL init_veg_constants      
      ! CALL init_sse_constants
      ! CALL init_sb_constants
      CALL init_pheno_constants    
      IF (debug_on()) CALL message(TRIM(routine), 'initialized quincy constants)')
    ENDIF

    IF (debug_on()) CALL message(TRIM(routine), 'done basic setup (models, process configs)')

  END SUBROUTINE jsbach_setup_models

  SUBROUTINE jsbach_setup_tiles(model_id)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_model_usecases, ONLY: init_usecase
    USE mo_jsb_lctlib_class,   ONLY: Read_lctlib

    INTEGER, INTENT(in) :: model_id

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:jsbach_setup_tiles'

    model => Get_model(model_id)

    model%lctlib => Read_lctlib('lctlib_nlct21.def', model%config%use_quincy)
    CALL message(TRIM(routine), 'Imported lctlib from '//'lctlib_nlct21.def')

    IF (model%config%usecase /= '') THEN
      ! Shortcut model configuration through usecase
      CALL init_usecase(model)
    ELSE
      ! Explicit model configuration through namelist
      ! TBD
    END IF

  END SUBROUTINE jsbach_setup_tiles

#endif
END MODULE mo_jsb_base
