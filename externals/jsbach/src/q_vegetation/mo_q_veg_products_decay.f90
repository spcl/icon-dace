!> QUINCY calculate decay of product pools
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
!>#### calculate decay of product pools
!>
MODULE mo_q_veg_products_decay
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_PP_FUEL_ID, VEG_BGCM_PP_PAPER_ID, VEG_BGCM_PP_FIBERBOARD_ID,      &
    &       VEG_BGCM_PP_OIRW_ID, VEG_BGCM_PP_PV_ID, VEG_BGCM_PP_SAWNWOOD_ID, VEG_BGCM_FPROD_DECAY_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_products_decay

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_products_decay'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Calculate the decay of the product pools
  !>
  SUBROUTINE update_products_decay(tile, options)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: VEG_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_constants,         ONLY: tau_pp_fuelwood, tau_pp_paper, tau_pp_fiberboard, &
      &                                 tau_pp_oirw, tau_pp_pv, tau_pp_sawnwood
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model                  !< the model
    TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_lctlib_element),   POINTER :: lctlib                 !< land-cover-type library - parameter across pft's
    REAL(wp)                          :: dtime                  !< timestep length
    INTEGER                           :: iblk, ics, ice, nc     !< dimensions
    CHARACTER(len=*), PARAMETER       :: routine = TRIM(modname)//':update_products_decay'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt1L2D :: veg_pp_fuel_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_paper_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_fiberboard_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_oirw_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_pv_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_sawnwood_mt
    dsl4jsb_Def_mt1L2D :: fprod_decay_mt
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_FUEL_ID, veg_pp_fuel_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_PAPER_ID, veg_pp_paper_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_FIBERBOARD_ID, veg_pp_fiberboard_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_OIRW_ID, veg_pp_oirw_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_PV_ID, veg_pp_pv_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_SAWNWOOD_ID, veg_pp_sawnwood_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_FPROD_DECAY_ID, fprod_decay_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    ! TODO: remove later, when having a harvest process
    !       for now: fill one matrix for testing purposes
    IF(ANY(veg_pp_sawnwood_mt == 0._wp)) THEN
      veg_pp_sawnwood_mt(:,:) = 6._wp
    ENDIF

    fprod_decay_mt(:,:) = 0._wp
    CALL calc_products_decay(nc, dtime, model%config%is_element_used(:), tau_pp_fuelwood, veg_pp_fuel_mt, fprod_decay_mt)
    CALL calc_products_decay(nc, dtime, model%config%is_element_used(:), tau_pp_paper, veg_pp_paper_mt, fprod_decay_mt)
    CALL calc_products_decay(nc, dtime, model%config%is_element_used(:), tau_pp_fiberboard, veg_pp_fiberboard_mt, fprod_decay_mt)
    CALL calc_products_decay(nc, dtime, model%config%is_element_used(:), tau_pp_oirw, veg_pp_oirw_mt, fprod_decay_mt)
    CALL calc_products_decay(nc, dtime, model%config%is_element_used(:), tau_pp_pv, veg_pp_pv_mt, fprod_decay_mt)
    CALL calc_products_decay(nc, dtime, model%config%is_element_used(:), tau_pp_sawnwood, veg_pp_sawnwood_mt, fprod_decay_mt)

  END SUBROUTINE update_products_decay

  ! ======================================================================================================= !
  !>
  !> calculate product pool decay -- for now using a 'simple' exponential decay
  !> -- in the future probably to be updated following suggestions in Nuetzel (2021, LMU master thesis)
  !>
  !> C from decaying wood products enters the atmospheric CO2 balance
  !> N and P from decaying wood products are simply recorded as flux but don't go anywhere (15N idem)
  !> 13C and 14C are following 12C
  !> 14C: radioactive decay (in veg_update_pools)
  !>
  SUBROUTINE calc_products_decay(nc, dtime, l_element_used , tau, veg_pp_mt, fprod_decay_mt)
    USE mo_jsb_math_constants,     ONLY: one_day, one_year
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,      INTENT(in) :: nc                  !< block dimension
    REAL(wp),     INTENT(in) :: dtime               !< timestep length
    LOGICAL,      INTENT(in) :: l_element_used(:)   !< indicates which elements are in use
    REAL(wp),     INTENT(in) :: tau                 !< decay time of decaying product pool
    REAL(wp),  INTENT(inout) :: veg_pp_mt(:,:)      !< decaying bgcm product pool
    REAL(wp),  INTENT(inout) :: fprod_decay_mt(:,:) !< bgcm flux of product decay
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                          :: fract_turnover    !< turnover fraction
    INTEGER                           :: ic, id_elem       !< loop counter
    CHARACTER(len=*), PARAMETER       :: routine = TRIM(modname)//':calc_products_decay'
    ! ----------------------------------------------------------------------------------------------------- !
    ! turnover fraction
    fract_turnover = 1.0_wp / tau / one_day / one_year * dtime

    DO ic = 1,nc
      DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
        IF (l_element_used(id_elem)) THEN
          fprod_decay_mt(id_elem, ic) = fprod_decay_mt(id_elem, ic) + (veg_pp_mt(id_elem, ic) * fract_turnover)
          veg_pp_mt(id_elem, ic) = veg_pp_mt(id_elem, ic) * (1._wp - fract_turnover)
        END IF
      END DO
    END DO

  END SUBROUTINE calc_products_decay

#endif
END MODULE mo_q_veg_products_decay
