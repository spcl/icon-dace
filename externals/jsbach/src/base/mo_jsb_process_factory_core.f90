!> Contains process factory.
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
MODULE mo_jsb_process_factory_core
#ifndef __NO_JSBACH__

  USE mo_jsb_process_class,    ONLY: t_jsb_process, Get_process_name
  USE mo_jsb_memory_class,     ONLY: t_jsb_memory
  !
  USE mo_a2l_memory_class,     ONLY: t_a2l_memory
  USE mo_l2a_memory_class,     ONLY: t_l2a_memory
  USE mo_seb_memory_class,     ONLY: t_seb_memory
  USE mo_turb_memory_class,    ONLY: t_turb_memory
  USE mo_sse_memory_class,     ONLY: t_sse_memory
  USE mo_hydro_memory_class,   ONLY: t_hydro_memory
  USE mo_assimi_memory_class,  ONLY: t_assimi_memory
  USE mo_pheno_memory_class,   ONLY: t_pheno_memory
  USE mo_rad_memory_class,     ONLY: t_rad_memory
  USE mo_carbon_memory_class,  ONLY: t_carbon_memory
  USE mo_disturb_memory_class, ONLY: t_disturb_memory
  USE mo_fuel_memory_class,    ONLY: t_fuel_memory
  USE mo_pplcc_memory_class,   ONLY: t_pplcc_memory
  USE mo_fage_memory_class,    ONLY: t_fage_memory
  USE mo_alcc_memory_class,    ONLY: t_alcc_memory
  USE mo_nlcc_memory_class,    ONLY: t_nlcc_memory
#ifndef __NO_JSBACH_HD__
  USE mo_hd_memory_class,      ONLY: t_hd_memory
#endif
#ifndef __NO_QUINCY__
  USE mo_q_assimi_memory_class,ONLY: t_q_assimi_memory
  USE mo_q_pheno_memory_class, ONLY: t_q_pheno_memory
  USE mo_q_rad_memory_class,   ONLY: t_q_rad_memory
  USE mo_spq_memory_class,     ONLY: t_spq_memory
  USE mo_veg_memory_class,     ONLY: t_veg_memory
  USE mo_sb_memory_class,      ONLY: t_sb_memory
#endif

  USE mo_rad_interface,        ONLY: Register_rad_tasks
  USE mo_rad_config_class,     ONLY: t_rad_config
  USE mo_seb_interface,        ONLY: Register_seb_tasks
  USE mo_seb_config_class,     ONLY: t_seb_config
  USE mo_turb_interface,       ONLY: Register_turb_tasks
  USE mo_turb_config_class,    ONLY: t_turb_config
  USE mo_sse_interface,        ONLY: Register_sse_tasks
  USE mo_sse_config_class,     ONLY: t_sse_config
  USE mo_hydro_interface,      ONLY: Register_hydro_tasks
  USE mo_hydro_config_class,   ONLY: t_hydro_config
  USE mo_assimi_interface,     ONLY: Register_assimi_tasks
  USE mo_assimi_config_class,  ONLY: t_assimi_config
  USE mo_pheno_interface,      ONLY: Register_pheno_tasks
  USE mo_pheno_config_class,   ONLY: t_pheno_config
#ifndef __NO_JSBACH_HD__
  USE mo_hd_interface,         ONLY: Register_hd_tasks
  USE mo_hd_config_class,      ONLY: t_hd_config
#endif
  USE mo_carbon_interface,      ONLY: Register_carbon_tasks
  USE mo_carbon_config_class,   ONLY: t_carbon_config
  USE mo_disturb_interface,     ONLY: Register_disturb_tasks
  USE mo_disturb_config_class,  ONLY: t_disturb_config
  USE mo_fuel_interface,        ONLY: Register_fuel_tasks
  USE mo_fuel_config_class,     ONLY: t_fuel_config
  USE mo_pplcc_interface,       ONLY: Register_pplcc_tasks
  USE mo_pplcc_config_class,    ONLY: t_pplcc_config
  USE mo_fage_interface,        ONLY: Register_fage_tasks
  USE mo_fage_config_class,     ONLY: t_fage_config
  USE mo_alcc_interface,        ONLY: Register_alcc_tasks
  USE mo_alcc_config_class,     ONLY: t_alcc_config
  USE mo_flcc_interface,        ONLY: Register_flcc_tasks
  USE mo_flcc_config_class,     ONLY: t_flcc_config
  USE mo_wlcc_interface,        ONLY: Register_wlcc_tasks
  USE mo_wlcc_config_class,     ONLY: t_wlcc_config
  USE mo_nlcc_interface,        ONLY: Register_nlcc_tasks
  USE mo_nlcc_config_class,     ONLY: t_nlcc_config
#ifndef __NO_QUINCY__
  USE mo_q_rad_interface,       ONLY: Register_q_rad_tasks
  USE mo_q_rad_config_class,    ONLY: t_q_rad_config
  USE mo_q_assimi_interface,    ONLY: Register_q_assimi_tasks
  USE mo_q_assimi_config_class, ONLY: t_q_assimi_config
  USE mo_q_pheno_interface,     ONLY: Register_q_pheno_tasks
  USE mo_q_pheno_config_class,  ONLY: t_q_pheno_config
  USE mo_q_veg_interface,       ONLY: Register_veg_tasks_quincy
  USE mo_veg_config_class,      ONLY: t_veg_config
  USE mo_q_sb_interface,        ONLY: Register_sb_tasks_quincy
  USE mo_sb_config_class,       ONLY: t_sb_config
  USE mo_q_spq_interface,       ONLY: Register_spq_tasks_quincy
  USE mo_spq_config_class,      ONLY: t_spq_config
#endif

  USE mo_jsb_process_class, ONLY: &
    & A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, HD_, RAD_, ASSIMI_, PHENO_, CARBON_ , DISTURB_, FUEL_, &
    & PPLCC_, FAGE_, ALCC_, FLCC_, WLCC_, NLCC_, VEG_, SB_, SPQ_, &
    & Q_PHENO_, Q_ASSIMI_, Q_RAD_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Create_process_memory
  PUBLIC :: Create_process, max_no_of_processes

  INTEGER, PARAMETER :: max_no_of_processes = 42

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_factory_core'

CONTAINS

  FUNCTION Create_process(iproc, model_id, namelist_filename) RESULT(return_ptr)

    INTEGER, INTENT(in) :: iproc
    INTEGER, INTENT(in) :: model_id
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CLASS(t_jsb_process), POINTER :: return_ptr

    ALLOCATE(t_jsb_process::return_ptr)

    return_ptr%id             = iproc
    return_ptr%name           = Get_process_name(iproc)
    return_ptr%owner_model_id = model_id

    ! Register tasks and set some process properties
    ! e.g. the land cover change (lcc) processes should set l_changes_fractions to .TRUE. here
    ! and processes that do not have a memory module also need to specify this here.
    SELECT CASE (iproc)
    CASE (A2L_, L2A_)
      ! Do nothing for A2L and L2A: these processes currently have no tasks and no config
    CASE (RAD_)
      CALL Register_rad_tasks(return_ptr, model_id)
      ALLOCATE(t_rad_config::return_ptr%config)
    CASE (PHENO_)
      CALL Register_pheno_tasks(return_ptr, model_id)
      ALLOCATE(t_pheno_config::return_ptr%config)
    CASE (SEB_)
      CALL Register_seb_tasks(return_ptr, model_id)
      ALLOCATE(t_seb_config::return_ptr%config)
    CASE (TURB_)
      CALL Register_turb_tasks(return_ptr, model_id)
      ALLOCATE(t_turb_config::return_ptr%config)
    CASE (SSE_)
      CALL Register_sse_tasks(return_ptr, model_id)
      ALLOCATE(t_sse_config::return_ptr%config)
    CASE (HYDRO_)
      CALL Register_hydro_tasks(return_ptr, model_id)
      ALLOCATE(t_hydro_config::return_ptr%config)
    CASE (ASSIMI_)
      CALL Register_assimi_tasks(return_ptr, model_id)
      ALLOCATE(t_assimi_config::return_ptr%config)
#ifndef __NO_JSBACH_HD__
    CASE (HD_)
      CALL Register_hd_tasks(return_ptr, model_id)
      ALLOCATE(t_hd_config::return_ptr%config)
#endif
    CASE (CARBON_)
      CALL Register_carbon_tasks(return_ptr, model_id)
      ALLOCATE(t_carbon_config::return_ptr%config)
    CASE (DISTURB_)
      CALL Register_disturb_tasks(return_ptr, model_id)
      ALLOCATE(t_disturb_config::return_ptr%config)
    CASE (FUEL_)
      CALL Register_fuel_tasks(return_ptr, model_id)
      ALLOCATE(t_fuel_config::return_ptr%config)
    CASE (PPLCC_)
      CALL Register_pplcc_tasks(return_ptr, model_id)
      ALLOCATE(t_pplcc_config::return_ptr%config)
    CASE (FAGE_)
      return_ptr%l_changes_fractions = .TRUE.
      CALL Register_fage_tasks(return_ptr, model_id)
      ALLOCATE(t_fage_config::return_ptr%config)
    CASE (ALCC_)
      return_ptr%l_changes_fractions = .TRUE.
      CALL Register_alcc_tasks(return_ptr, model_id)
      ALLOCATE(t_alcc_config::return_ptr%config)
    CASE (FLCC_)
      return_ptr%l_changes_fractions = .TRUE.
      return_ptr%has_memory = .FALSE.
      CALL Register_flcc_tasks(return_ptr, model_id)
      ALLOCATE(t_flcc_config::return_ptr%config)
    CASE (WLCC_)
      return_ptr%l_changes_fractions = .TRUE.
      return_ptr%has_memory = .FALSE.
      CALL Register_wlcc_tasks(return_ptr, model_id)
      ALLOCATE(t_wlcc_config::return_ptr%config)
    CASE (NLCC_)
      return_ptr%l_changes_fractions = .TRUE.
      CALL Register_nlcc_tasks(return_ptr, model_id)
      ALLOCATE(t_nlcc_config::return_ptr%config)
#ifndef __NO_QUINCY__
    CASE (Q_RAD_)
      CALL Register_q_rad_tasks(return_ptr, model_id)
      ALLOCATE(t_q_rad_config::return_ptr%config)
    CASE (Q_ASSIMI_)
      CALL Register_q_assimi_tasks(return_ptr, model_id)
      ALLOCATE(t_q_assimi_config::return_ptr%config)
    CASE (Q_PHENO_)
      CALL Register_q_pheno_tasks(return_ptr, model_id)
      ALLOCATE(t_q_pheno_config::return_ptr%config)
    CASE (VEG_)
      CALL Register_veg_tasks_quincy(return_ptr, model_id)
      ALLOCATE(t_veg_config::return_ptr%config)
    CASE (SB_)
      CALL Register_sb_tasks_quincy(return_ptr, model_id)
      ALLOCATE(t_sb_config::return_ptr%config)
    CASE (SPQ_)
      CALL Register_spq_tasks_quincy(return_ptr, model_id)
      ALLOCATE(t_spq_config::return_ptr%config)
#endif
    CASE DEFAULT
      return_ptr => NULL()
    END SELECT

    IF (ASSOCIATED(return_ptr)) THEN
      IF (ASSOCIATED(return_ptr%config)) THEN
        return_ptr%config%namelist_filename = TRIM(namelist_filename)
      END IF
    END IF

  END FUNCTION Create_process

  FUNCTION Create_process_memory(iproc) RESULT(return_ptr)

    USE mo_a2l_memory_class,     ONLY: max_no_of_a2l_vars => max_no_of_vars
    USE mo_l2a_memory_class,     ONLY: max_no_of_l2a_vars => max_no_of_vars
    USE mo_seb_memory_class,     ONLY: max_no_of_seb_vars => max_no_of_vars
    USE mo_turb_memory_class,    ONLY: max_no_of_turb_vars => max_no_of_vars
    USE mo_sse_memory_class,     ONLY: max_no_of_sse_vars => max_no_of_vars
    USE mo_hydro_memory_class,   ONLY: max_no_of_hydro_vars => max_no_of_vars
    USE mo_assimi_memory_class,  ONLY: max_no_of_assimi_vars => max_no_of_vars
    USE mo_pheno_memory_class,   ONLY: max_no_of_pheno_vars => max_no_of_vars
    USE mo_rad_memory_class,     ONLY: max_no_of_rad_vars => max_no_of_vars
    USE mo_carbon_memory_class,  ONLY: max_no_of_carbon_vars => max_no_of_vars
    USE mo_disturb_memory_class, ONLY: max_no_of_disturb_vars => max_no_of_vars
    USE mo_fuel_memory_class,    ONLY: max_no_of_fuel_vars => max_no_of_vars
    USE mo_pplcc_memory_class,   ONLY: max_no_of_pplcc_vars => max_no_of_vars
    USE mo_fage_memory_class,    ONLY: max_no_of_fage_vars => max_no_of_vars
    USE mo_alcc_memory_class,    ONLY: max_no_of_alcc_vars => max_no_of_vars
    USE mo_nlcc_memory_class,    ONLY: max_no_of_nlcc_vars => max_no_of_vars
#ifndef __NO_JSBACH_HD__
    USE mo_hd_memory_class,      ONLY: max_no_of_hd_vars => max_no_of_vars
#endif
#ifndef __NO_QUINCY__
    USE mo_q_assimi_memory_class,ONLY: max_no_of_q_assimi_vars => max_no_of_vars
    USE mo_q_pheno_memory_class, ONLY: max_no_of_q_pheno_vars => max_no_of_vars
    USE mo_q_rad_memory_class,   ONLY: max_no_of_q_rad_vars => max_no_of_vars
    USE mo_veg_memory_class,     ONLY: max_no_of_veg_vars => max_no_of_vars
    USE mo_sb_memory_class,      ONLY: max_no_of_sb_vars => max_no_of_vars
    USE mo_spq_memory_class,     ONLY: max_no_of_spq_vars => max_no_of_vars
#endif
    INTEGER, INTENT(in) :: iproc
    CLASS(t_jsb_memory), POINTER   :: return_ptr

    SELECT CASE(iproc)
    CASE (A2L_)
      ALLOCATE(t_a2l_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_a2l_vars
    CASE (L2A_)
      ALLOCATE(t_l2a_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_l2a_vars
    CASE (SEB_)
      ALLOCATE(t_seb_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_seb_vars
    CASE (TURB_)
      ALLOCATE(t_turb_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_turb_vars
    CASE (SSE_)
      ALLOCATE(t_sse_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_sse_vars
    CASE (HYDRO_)
      ALLOCATE(t_hydro_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_hydro_vars
    CASE (ASSIMI_)
      ALLOCATE(t_assimi_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_assimi_vars
    CASE (PHENO_)
      ALLOCATE(t_pheno_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_pheno_vars
    CASE (RAD_)
      ALLOCATE(t_rad_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_rad_vars
#ifndef __NO_JSBACH_HD__
    CASE (HD_)
      ALLOCATE(t_hd_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_hd_vars
#endif
    CASE (CARBON_)
      ALLOCATE(t_carbon_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_carbon_vars
    CASE (DISTURB_)
      ALLOCATE(t_disturb_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_disturb_vars
    CASE (FUEL_)
      ALLOCATE(t_fuel_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_fuel_vars
    CASE (PPLCC_)
      ALLOCATE(t_pplcc_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_pplcc_vars
    CASE (FAGE_)
      ALLOCATE(t_fage_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_fage_vars
    CASE (ALCC_)
      ALLOCATE(t_alcc_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_alcc_vars
    CASE (NLCC_)
      ALLOCATE(t_nlcc_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_nlcc_vars
#ifndef __NO_QUINCY__
    CASE (Q_ASSIMI_)
      ALLOCATE(t_q_assimi_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_q_assimi_vars
    CASE (Q_PHENO_)
      ALLOCATE(t_q_pheno_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_q_pheno_vars
    CASE (Q_RAD_)
      ALLOCATE(t_q_rad_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_q_rad_vars
    CASE (VEG_)
      ALLOCATE(t_veg_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_veg_vars
    CASE (SB_)
      ALLOCATE(t_sb_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_sb_vars
    CASE (SPQ_)
      ALLOCATE(t_spq_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_spq_vars
#endif
    END SELECT

    ALLOCATE(return_ptr%vars(return_ptr%max_no_of_vars))

  END FUNCTION Create_process_memory

#endif
END MODULE mo_jsb_process_factory_core
