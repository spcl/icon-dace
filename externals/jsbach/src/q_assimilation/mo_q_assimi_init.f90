!> QUINCY assimilation variables init
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
!>#### initialization of assimilation memory variables using, e.g., ic & bc input files
!>
MODULE mo_q_assimi_init
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: q_assimi_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_assimi_init'

CONTAINS

  ! ======================================================================================================= !
  !> Run assimilation init
  !>
  SUBROUTINE q_assimi_init(tile)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_process_class,     ONLY: Q_ASSIMI_
    USE mo_q_assimi_parameters,   ONLY: t_jmax_opt_min
    dsl4jsb_Use_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':q_assimi_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
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
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(Q_ASSIMI_)) RETURN
    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of assimi memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_air)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soa)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, soa_tsoa_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soa_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_air_daytime)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_air_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_ps)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_ps_daytime)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_gs_daytime)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_gs_daytime_dacc)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_gs_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_gs_tfrac_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, beta_soil_gs_tcnl_mavg)
    dsl4jsb_Get_var2D_onDomain(Q_ASSIMI_, t_jmax_opt)
    ! ----------------------------------------------------------------------------------------------------- !

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
  END SUBROUTINE q_assimi_init

#endif
END MODULE mo_q_assimi_init
