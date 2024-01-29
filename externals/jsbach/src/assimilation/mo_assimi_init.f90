!> Initialization of the the assimilation memory
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
MODULE mo_assimi_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: get_model

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: assimi_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Run assimilation init
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE assimi_init(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_process_class,     ONLY: ASSIMI_
    USE mo_assimi_constants,      ONLY: t_jmax_opt_min

    ! Use of process configurations (t_PROC_config)
    !dsl4jsb_Use_config(ASSIMI_)    
    
    ! USE PROC_memory
    dsl4jsb_Use_memory(ASSIMI_)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':assimi_init'

    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    !dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onDomain      :: beta_air
    dsl4jsb_Real2D_onDomain      :: beta_soa
    dsl4jsb_Real2D_onDomain      :: soa_tsoa_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soa_tphen_mavg
    dsl4jsb_Real2D_onDomain      :: beta_air_daytime
    dsl4jsb_Real2D_onDomain      :: beta_air_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: beta_air_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps_daytime
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_ps_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_daytime
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_daytime_dacc
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_tphen_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_tfrac_mavg
    dsl4jsb_Real2D_onDomain      :: beta_soil_gs_tcnl_mavg
    dsl4jsb_Real2D_onDomain      :: t_jmax_opt
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    model => Get_model(tile%owner_model_id)
    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN
    IF (.NOT. model%config%use_quincy) RETURN
    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of assimi memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
    !dsl4jsb_Get_config(ASSIMI_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soa)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, soa_tsoa_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soa_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air_daytime)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps_daytime)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_daytime)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, beta_soil_gs_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(ASSIMI_, t_jmax_opt)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    beta_air(:,:)                     = 1.0_wp
    beta_soa(:,:)                     = 1.0_wp
    soa_tsoa_mavg(:,:)                = -10.0_wp
    beta_soa_tphen_mavg(:,:)          = 1.0_wp
    beta_soil_ps(:,:)                 = 1.0_wp
    beta_soil_gs(:,:)                 = 1.0_wp
    
    beta_air_daytime(:,:)             = 1.0_wp
    beta_air_daytime_dacc(:,:)        = 1.0_wp
    beta_soil_ps_daytime(:,:)         = 1.0_wp
    beta_soil_ps_daytime_dacc(:,:)    = 1.0_wp
    beta_soil_gs_daytime(:,:)         = 1.0_wp
    beta_soil_gs_daytime_dacc(:,:)    = 1.0_wp
    
    beta_air_tfrac_mavg(:,:)          = 1.0_wp
    beta_air_tcnl_mavg(:,:)           = 1.0_wp
    beta_soil_ps_tfrac_mavg(:,:)      = 1.0_wp
    beta_soil_ps_tcnl_mavg(:,:)       = 1.0_wp
    beta_soil_gs_tphen_mavg(:,:)      = 1.0_wp
    beta_soil_gs_tfrac_mavg(:,:)      = 1.0_wp
    beta_soil_gs_tcnl_mavg(:,:)       = 1.0_wp

    t_jmax_opt(:,:)                   = t_jmax_opt_min 

  END SUBROUTINE assimi_init

#endif
END MODULE mo_assimi_init
