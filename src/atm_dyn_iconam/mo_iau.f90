!
! Routines and functions for the Incremental Analysis Update (IAU)
!
! This module contains routines and functions for the Incremental Analysis
! Update (IAU), including the iterative IAU. The characteristics of the
! iterative IAU are described in the following using the example of a global
! NWP forecast run:
!
! IAU iteration
!
!                     input
!                       /
!                      /
!                     /
!          ........../
!         /
!        /
!       /
!      /
!     /
!  -90min               0min              90min
! ---|------------------|------------------|------------->
!    |//////////////////| - - - - - - - - - - - - - - - ->
!                               free forecast (iteration = false)
!    \________IAU_______/
!                       |
!                       /
!                      /
!                     /
!          ........../
!         /   reset
!        /
!       /
!      /
!     /
!  -90min               0min              90min
! ---|------------------|------------------|------------->
!    |//////////////////|//////////////////| free forecast
!
!    \_________________IAU________________/
!
! @Literature:
! Bloom, S. C., Takacs, L. L., da Silva, A. M., & Ledvina, D. (1996).
! Data Assimilation Using Incremental Analysis Updates, Monthly Weather Review,
! 124(6), 1256-1271
! Polavarapu, S., Ren, S., Clayton, A. M., Sankey, D., & Rochon, Y. (2004).
! On the Relationship between Incremental Analysis Updating and
! Incremental Digital Filtering, Monthly Weather Review, 132(10), 2495-2502
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

MODULE mo_iau

  USE mo_kind,                    ONLY: wp
  USE mo_impl_constants,          ONLY: LSS_JSBACH, SUCCESS
  USE mo_exception,               ONLY: finish
  USE mo_math_constants,          ONLY: pi
  USE mo_save_restore,            ONLY: save_var_group_state, restore_var_group_state, &
    &                                   reinit_var_group_state
  USE mo_model_domain,            ONLY: t_patch
  USE mo_nonhydro_types,          ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_parallel_config,         ONLY: nproma
  USE mo_grid_config,             ONLY: n_dom
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config
  USE mo_assimilation_config,     ONLY: assimilation_config
  USE mo_initicon_config,         ONLY: type_iau_wgt, is_iau_active, iau_wgt_dyn, iau_wgt_adv, &
    &                                   qcana_mode, qiana_mode, qrsgana_mode
  USE mo_run_config,              ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqh, &
    &                                   iqm_max, iqni, iqnc, iqnr, iqns, iqng, iqnh, iqbin, iqb_i, iqb_e, &
    &                                   ldass_lhn
  USE mo_dynamics_config,         ONLY: nnow, nnew, nnow_rcf, nnew_rcf
  USE mo_advection_config,        ONLY: advection_config
  USE mo_nonhydrostatic_config,   ONLY: kstart_moist
  USE mo_fortran_tools,           ONLY: assert_acc_device_only, assert_acc_host_only
  USE mo_hash_table,              ONLY: t_HashTable, hashTable_make
  USE mo_util_texthash,           ONLY: text_hash, text_isEqual
  USE mo_satad,                   ONLY: qsat_rho
  USE mo_nh_diagnose_pres_temp,   ONLY: diag_pres, diag_temp
  USE mo_timer,                   ONLY: ltimer, timer_iau_save_restore, timer_start, timer_stop

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_iau'


  ! variables
  PUBLIC :: t_saveinit_state     ! state for saving initial state for double IAU runs

  ! subroutines
  !
  PUBLIC :: save_initial_state
  PUBLIC :: reset_to_initial_state
  PUBLIC :: compute_iau_wgt
  PUBLIC :: iau_update_tracer


  ! state for saving initial state
  TYPE :: t_saveinit_state

    TYPE(t_HashTable), POINTER :: fields     => NULL()
    TYPE(t_HashTable), POINTER :: fields_jsb => NULL()

  CONTAINS
    PROCEDURE :: init     => saveinit_state_init     !< constructor
    PROCEDURE :: finalize => saveinit_state_finalize !< destructor
  END TYPE t_saveinit_state


  TYPE(t_saveinit_state), ALLOCATABLE  :: saveinit(:)

CONTAINS


  SUBROUTINE saveinit_state_init(saveinit_data)
    CLASS(t_saveinit_state), INTENT(INOUT) :: saveinit_data

    saveinit_data%fields     => hashTable_make(text_hash, text_isEqual)
    saveinit_data%fields_jsb => hashTable_make(text_hash, text_isEqual)

  END SUBROUTINE saveinit_state_init


  SUBROUTINE saveinit_state_finalize(saveinit_data)
    CLASS(t_saveinit_state), INTENT(INOUT) :: saveinit_data

    IF (ASSOCIATED(saveinit_data%fields)) THEN
      CALL saveinit_data%fields%destruct
      !
      DEALLOCATE(saveinit_data%fields)
    END IF

    IF (ASSOCIATED(saveinit_data%fields_jsb)) THEN
      CALL saveinit_data%fields_jsb%destruct
      !
      DEALLOCATE(saveinit_data%fields_jsb)
    END IF

  END SUBROUTINE saveinit_state_finalize


  !----------------------------------------------------------------------------
  !>
  !! Saves the initial state of NWP applications for the IAU iteration mode.
  !!
  !! Running ICON in iterative IAU mode requires the model to be
  !! reinitialized after the first IAU iteration step.
  !! Technically, the reinitialization is based on a save/restore mechanism.
  !! A well chosen subset of ICON's prognostic and diagnostic fields is saved
  !! prior to the first iteration and restored thereafter, just before the
  !! second iteration.
  !! For efficiency reasons the mechanism distinguishes between two types of fields:
  !! I ) fields that must be saved and restored (such as the first guess fields)
  !! II) fields for which a reinitialization by its initial value field%info%initval
  !!     is sufficient.
  !!
  !! These fields are collected in one of the following variable groups:
  !!
  !! IAU_RESTORE_VARS ! Dynamical core + NWP physics variables that require save/restore
  !! IAU_INIT_VARS    ! Dynamical core + NWP physics variables that require reinitialization
  !!                    by a constant.
  !! JSB_INIT_VARS    ! JSBACH specific variables
  !!
  !! Additional fields may be added or deleted if necessary. To do so, please use the
  !! 'in_group' optional argument of the add_var/add_ref calls and assign the desired
  !! field(s) to one of the above mentioned groups.
  !!
  SUBROUTINE save_initial_state(p_patch, lacc)

    TYPE(t_patch),    INTENT(IN) :: p_patch(:)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jg
    INTEGER :: ierrstat
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':save_initial_state'

    ! make sure that accelerator is not used
    CALL assert_acc_host_only(routine, lacc)

    IF (ltimer) CALL timer_start(timer_iau_save_restore)


    ALLOCATE(saveinit(n_dom), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, 'Allocation failed for saveinit')

    DO jg = 1, n_dom

      IF(.NOT. p_patch(jg)%ldom_active) CYCLE

      CALL saveinit(jg)%init

      ! save dynamics + NWP physics fields
      CALL save_var_group_state('iau_restore_vars', p_patch(jg), saveinit(jg)%fields)

      ! save jsbach specific fields
      IF (atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH) THEN
        CALL save_var_group_state('jsb_init_vars', p_patch(jg), saveinit(jg)%fields_jsb)
      END IF

    ENDDO  ! jg

    IF (ltimer) CALL timer_stop(timer_iau_save_restore)

  END SUBROUTINE save_initial_state


  !----------------------------------------------------------------------------
  !>
  !! Restores the initial state of NWP applications for the IAU iteration mode.
  !!
  SUBROUTINE restore_initial_state(p_patch, p_nh, lacc)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),          INTENT(INOUT) :: p_nh(:)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jg, ic, je, jb
    INTEGER :: ierrstat

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':restore_initial_state'

    ! make sure that accelerator is not used
    CALL assert_acc_host_only(routine, lacc)

    IF (ltimer) CALL timer_start(timer_iau_save_restore)

    DO jg = 1, n_dom

      IF(.NOT. p_patch(jg)%ldom_active) CYCLE

      ! restore dynamics  + NWP physics fields
      CALL restore_var_group_state('iau_restore_vars', p_patch(jg), saveinit(jg)%fields)

      ! reinitialize fields which are included in the group IAU_INIT_VARS
      CALL reinit_var_group_state('iau_init_vars', p_patch(jg))

      ! restore jsbach specific fields
      IF (atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH) THEN
        CALL restore_var_group_state('jsb_init_vars', p_patch(jg), saveinit(jg)%fields_jsb)
      END IF

      ! For the limited-area mode and one-way nesting, we also need to reset grf_tend_vn on the nudging points
      !
      DO ic = 1, p_nh(jg)%metrics%nudge_e_dim
        je = p_nh(jg)%metrics%nudge_e_idx(ic)
        jb = p_nh(jg)%metrics%nudge_e_blk(ic)
        p_nh(jg)%diag%grf_tend_vn(je,:,jb) = 0._wp
      ENDDO

      ! deallocate
      CALL saveinit(jg)%finalize

    ENDDO

    DEALLOCATE(saveinit, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, 'Deallocation failed for saveinit')

    IF (ltimer) CALL timer_stop(timer_iau_save_restore)

  END SUBROUTINE restore_initial_state


  !----------------------------------------------------------------------------
  !>
  !! Wrapper routine for restore_initial_state, which in addition restores
  !! several control fields such as
  !! * time level arrays
  !! * mtime events for NWP physics, LHN, and DACE
  !!
  SUBROUTINE reset_to_initial_state(p_patch, p_nh)

    TYPE(t_patch),        INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),     INTENT(INOUT) :: p_nh(:)

    INTEGER :: jg

    ! reinitialize mtime events
    !
    ! NOTE: reinitialization of DACE event group 'mec_Events ' missing.
    !       I am not sure whether a reinitialization is strictly needed.
    !
    DO jg = 1, n_dom
      ! for NWP physics
      CALL atm_phy_nwp_config(jg)%phyProcs%reinitEvents()
      IF (ldass_lhn) THEN
        ! for latent heat nudging
        CALL assimilation_config(jg)%dass_g%reinitEvents()
      ENDIF
    ENDDO
    !
    ! reinitialize time level arrays
    nnow(:)     = 1
    nnow_rcf(:) = 1
    nnew(:)     = 2
    nnew_rcf(:) = 2
    !
    ! reset model fields to its initial state
    ! Note that this must happen after the reinitialization
    ! of the time level arrays.
     CALL restore_initial_state(p_patch(1:), p_nh)

  END SUBROUTINE reset_to_initial_state


  !>
  !! Compute weights for incremental analysis update
  !!
  !! Compute weights for incremental analysis update.
  !! 2 weights are provided:
  !! - iau_wgt_dyn can be used for all fields that need to be updated
  !!   every (fast) dynamics time step
  !! - iau_wgt_adv can be used for all fields that need to be updated
  !!   every (slow) advection time step.
  !!
  SUBROUTINE compute_iau_wgt(sim_time, dt, dt_iau, lreset_wgt_adv)

    REAL(wp)        , INTENT(IN)  :: sim_time          !< Simulation time since model
                                                       !< start
    REAL(wp)        , INTENT(IN)  :: dt                !< time step
    REAL(wp)        , INTENT(IN)  :: dt_iau            !< width of IAU window
    LOGICAL         , INTENT(IN)  :: lreset_wgt_adv    !< If true, reset the accumulated weight for the advective time step

    ! local variables
    REAL(wp)  :: time_iau_elapsed                      !< elapsed time since IAU start [s]
    REAL(wp)  :: fct_eval                              !< result of top-hat or sin2 function

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_init_utils:compute_iau_wgt'
    !-------------------------------------------------------------------------

    ! initialize
    IF (lreset_wgt_adv) iau_wgt_adv   = 0._wp

    ! compute elapsed time (in s) since IAU start
    !
    ! trivial so far, however will be changed to mtime when the functionality of
    ! computing the timedelta between two dates becomes available.
    time_iau_elapsed = sim_time


    IF (time_iau_elapsed <= dt_iau) THEN
      is_iau_active = .TRUE.

      SELECT CASE (type_iau_wgt)
        CASE(1)  ! top-hat function
          fct_eval = iau_top_hat(dt_iau,time_iau_elapsed)

        CASE(2)  ! sin2 function
          fct_eval = iau_sin2   (dt_iau,time_iau_elapsed)

        CASE(3)  ! sin function
          fct_eval = iau_sin    (dt_iau,time_iau_elapsed)

        CASE default
          CALL finish(routine,&
                      'Invalid IAU weighting function. Must be 1, 2 or 3.')
      END SELECT

      ! compute weights by multiplying with the time step
      iau_wgt_dyn = fct_eval * dt
      iau_wgt_adv = iau_wgt_adv + iau_wgt_dyn

    ELSE
      is_iau_active = .FALSE.
      iau_wgt_dyn   = 0._wp
      iau_wgt_adv   = 0._wp
    ENDIF

!!$write(0,*) "sim_time, is_iau_active, iau_wgt_dyn, iau_wgt_adv: ", &
!!$  & sim_time, is_iau_active, iau_wgt_dyn, iau_wgt_adv

  END SUBROUTINE compute_iau_wgt


  !>
  !! Evaluates top-hat function at a particular point in time
  !!
  !! Evaluates top-hat function at a particular point in time
  !! Top-hat function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  FUNCTION iau_top_hat (dt,cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = 1._wp/dt
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_top_hat


  !>
  !! Evaluates SIN2 function at a particular point in time
  !!
  !! Evaluates SIN2 function at a particular point in time
  !! SIN2 function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  FUNCTION iau_sin2 (dt,cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = (2._wp/dt) * SIN(pi*cur_time/dt)**2
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_sin2

  !>
  !! Evaluates SIN function at a particular point in time
  !!
  !! Evaluates SIN function at a particular point in time
  !! SIN function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  FUNCTION iau_sin (dt, cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = ((PI/2._wp)/dt) * SIN(PI*cur_time/dt)
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_sin


  !
  ! Add IAU increment to qv during IAU phase
  !
  ! Add analysis increments from data assimilation to qv
  !
  ! Initial revision by Daniel Reinert, DWD (2018-05-18)
  ! Previously, this code snippet was part of nh_update_tracer_phy
  ! 
  SUBROUTINE iau_update_tracer( pt_prog, p_metrics, pt_diag, pt_prog_rcf, &
    &                     jg, jb, i_startidx, i_endidx, kend, lacc )

    TYPE(t_nh_prog)    ,INTENT(IN)   :: pt_prog      !< NH prog state at dynamic time step
    TYPE(t_nh_metrics) ,INTENT(IN)   :: p_metrics    !< NH metrics variables
    TYPE(t_nh_diag)    ,INTENT(INOUT):: pt_diag      !< the diagnostic variables
    TYPE(t_nh_prog)    ,INTENT(INOUT):: pt_prog_rcf  !< the tracer field at
                                                      !< reduced calling frequency
    INTEGER            ,INTENT(IN)   :: jg           !< domain ID
    INTEGER            ,INTENT(IN)   :: jb           !< block index
    INTEGER            ,INTENT(IN)   :: i_startidx   !< hor. start idx
    INTEGER            ,INTENT(IN)   :: i_endidx     !< hor. end idx
    INTEGER            ,INTENT(IN)   :: kend         !< vert. end idx
    LOGICAL, OPTIONAL  ,INTENT(IN)   :: lacc         ! If true, use openacc

    ! Local variables
    INTEGER  :: jk,jc
    INTEGER  :: iqb
    REAL(wp) :: zqin
    REAL(wp) :: zrhw(nproma, kend) ! relative humidity w.r.t. water


    CALL assert_acc_device_only("iau_update_tracer", lacc)

    ! add analysis increments from data assimilation to qv
    !
    ! Diagnose pressure and temperature for subsequent calculations
    CALL diag_temp (pt_prog, pt_prog_rcf, advection_config(jg)%trHydroMass%list, pt_diag, &
                    jb, i_startidx, i_endidx, 1, kstart_moist(jg), kend)
    CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, 1, kend)

    !$ACC DATA CREATE(zrhw) PRESENT(pt_prog, p_metrics, pt_diag, pt_prog_rcf, atm_phy_nwp_config)

    ! Compute relative humidity w.r.t. water
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jk = 1, kend
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        zrhw(jc,jk) = pt_prog_rcf%tracer(jc,jk,jb,iqv)/qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
      ENDDO
    ENDDO

    ! GZ: This loop needs to be split for correct vectorization because rhoc_incr is allocated for qcana_mode >= 1 only;
    !     otherwise, the NEC runs into a segfault. Likewise, the remaining case selections need to be done outside the
    !     vectorized loops in order to avoid invalid memory accesses.
    !$ACC LOOP SEQ
    DO jk = 1, kend
      IF (qcana_mode >= 1) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqin)
          DO jc = i_startidx, i_endidx
          IF (qcana_mode == 2 .AND. pt_prog_rcf%tracer(jc,jk,jb,iqc) > 0._wp) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + &
              iau_wgt_adv*pt_diag%rhov_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb)
            pt_prog_rcf%tracer(jc,jk,jb,iqc) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqc) + &
              iau_wgt_adv*pt_diag%rhoc_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          ELSE 
            zqin = (pt_diag%rhov_incr(jc,jk,jb)+pt_diag%rhoc_incr(jc,jk,jb))/pt_prog%rho(jc,jk,jb)
            ! DA increments of humidity are limited to positive values if p > 150 hPa and RH < 2% or QV < 5.e-7
            IF (pt_diag%pres(jc,jk,jb) > 15000._wp .AND. zrhw(jc,jk) < 0.02_wp .OR. &
              pt_prog_rcf%tracer(jc,jk,jb,iqv) < 5.e-7_wp) zqin = MAX(0._wp, zqin)
            pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqv) + iau_wgt_adv*zqin)
          ENDIF
        ENDDO
      ELSE
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqin)
          DO jc = i_startidx, i_endidx
          zqin = pt_diag%rhov_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb)
          ! DA increments of humidity are limited to positive values if p > 150 hPa and RH < 2% or QV < 5.e-7
          IF (pt_diag%pres(jc,jk,jb) > 15000._wp .AND. zrhw(jc,jk) < 0.02_wp .OR. &
            pt_prog_rcf%tracer(jc,jk,jb,iqv) < 5.e-7_wp) zqin = MAX(0._wp, zqin)
          pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqv) + iau_wgt_adv*zqin)
        ENDDO
      ENDIF

      IF (qiana_mode > 0) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqi) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqi) + &
            iau_wgt_adv*pt_diag%rhoi_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (qrsgana_mode > 0) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqr) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqr) + &
            iau_wgt_adv * pt_diag%rhor_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          pt_prog_rcf%tracer(jc,jk,jb,iqs) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqs) + &
            iau_wgt_adv * pt_diag%rhos_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (qrsgana_mode > 0 .AND. iqg <= iqm_max) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqg) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqg) + &
            iau_wgt_adv * pt_diag%rhog_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%l2moment) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          IF (qcana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqnc) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnc) + &
                 iau_wgt_adv * pt_diag%rhonc_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
          IF (qiana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqni) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqni) + &
                 iau_wgt_adv * pt_diag%rhoni_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
          IF (qrsgana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqh) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqh) + &
                 iau_wgt_adv * pt_diag%rhoh_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqnr) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnr) + &
                 iau_wgt_adv * pt_diag%rhonr_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqns) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqns) + &
                 iau_wgt_adv * pt_diag%rhons_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqng) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqng) + &
                 iau_wgt_adv * pt_diag%rhong_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqnh) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnh) + &
                 iau_wgt_adv * pt_diag%rhonh_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%lsbm) THEN
        IF (qrsgana_mode > 0) THEN
          !$ACC LOOP SEQ
          DO iqb = iqb_i, iqb_e
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO jc = i_startidx, i_endidx
              pt_prog_rcf%tracer(jc,jk,jb,iqbin(iqb)) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqbin(iqb)) + &
                   iau_wgt_adv * pt_diag%rhonh_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            END DO
          ENDDO
        ENDIF
      ENDIF

    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE iau_update_tracer

END MODULE mo_iau
