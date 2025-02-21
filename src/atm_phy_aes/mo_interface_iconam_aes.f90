!
!
!------------------------------------------------------------------------------
! This module provides the interface between dynamics&transport and physics.
!
!------------------------------------------------------------------------------
!
! Interface between ICONAM dynamics&transport and AES physics
!
! The provisional atmospheric state existing after dynamics&transport is used
! as input for the phyiscs parameterizations. The resulting physical tendencies
! are used to update and thus obtain the final atmospheric state of the time
! step. Thus the entire physics is treated as "fast" physics.
!
! Procedure:
! - Prepare the input fields from the prognostic state,
! - Call a subroutine to prepare the physics boundary conditions,
! - Call a subroutine to calculate the physics tendencies,
! - Call a subroutine to couple atmosphere and ocean, if needed,
! - Update the prognostic state with the physics tendencies,
! - Synchronize the prognostic state
!
! Memory:
! - Global atmospheric state fields needed in phyiscs are bundled in the
!   'field' (f) data type, using references to the memory of the original
!   fields organized in a number of data types used in the dynamics and
!   transport components.
! - Tendencies dx/dt resulting from physics are stored in the 'tend' (t) data
!   type, which has its own memory. These physics tendencies are used for the
!   final update of the prognostic variables before the end of the interface.
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

!----------------------------
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

!----------------------------

MODULE mo_interface_iconam_aes

  USE mo_kind                  ,ONLY: wp, vp
  USE mo_exception             ,ONLY: print_value

  USE mo_timer                 ,ONLY: ltimer, timer_start, timer_stop, &
    &                                 timer_dyn2phy, timer_d2p_sync,   &
    &                                 timer_aes_bcs, timer_aes_phy,    &
    &                                 timer_phy2dyn, timer_p2d_sync

  USE mo_run_config            ,ONLY: ntracer, iqv
  USE mo_advection_config      ,ONLY: advection_config
  USE mo_aes_phy_config        ,ONLY: aes_phy_config, aes_phy_tc, dt_zero
  USE mo_aes_vdf_config        ,ONLY: aes_vdf_config

  USE mtime                    ,ONLY: t_datetime  => datetime , newDatetime , deallocateDatetime     ,&
    &                                 t_timedelta => timedelta, newTimedelta, deallocateTimedelta    ,&
    &                                 max_timedelta_str_len   , getPTStringFromSeconds               ,&
    &                                 OPERATOR(+), OPERATOR(>)

  USE mo_model_domain          ,ONLY: t_patch
  USE mo_intp_data_strc        ,ONLY: t_int_state
  USE mo_nonhydro_types        ,ONLY: t_nh_prog, t_nh_diag
  USE mo_aes_phy_memory        ,ONLY: t_aes_phy_field, prm_field, &
    &                                 t_aes_phy_tend, prm_tend

  USE mo_sync                  ,ONLY: sync_c, sync_e, sync_patch_array_mult
  USE mo_intp_rbf              ,ONLY: rbf_vec_interpol_cell
  USE mo_loopindices           ,ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants        ,ONLY: min_rlcell_int, grf_bdywidth_c, &
    &                                 min_rledge_int, grf_bdywidth_e

  USE mo_math_constants        ,ONLY: rad2deg
  USE mo_physical_constants    ,ONLY: rd, p0ref, rd_o_cpd, vtmpc1

  USE mo_aes_phy_bcs           ,ONLY: aes_phy_bcs
  USE mo_aes_phy_main          ,ONLY: aes_phy_main

  USE mo_jsb_interface         ,ONLY: jsbach_start_timestep, jsbach_finish_timestep
  
  

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interface_iconam_aes

CONTAINS
  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE interface_iconam_aes(dt,            & !in
    &                             datetime_new,  & !in
    &                             patch,         & !in
    &                             int_state,     & !in
    &                             dyn_now,       & !inout
    &                             dyn_new,       & !inout
    &                             adv_now,       & !inout
    &                             adv_new,       & !inout
    &                             diag)            !inout

    !
    !> Arguments
    !
    REAL(wp)              , INTENT(in)            :: dt           !< advective time step
    TYPE(t_datetime)      , POINTER               :: datetime_new !< date and time at the end of this time step

    TYPE(t_patch)         , INTENT(inout), TARGET :: patch        !< grid/patch info
    TYPE(t_int_state)     , INTENT(in)   , TARGET :: int_state    !< interpolation state

    TYPE(t_nh_diag)       , INTENT(inout), TARGET :: diag         !< diagnostic variables
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: dyn_now      !< progn. vars in    dynamics  for wind, temp. rho, ...
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: dyn_new      !< progn. vars after dynamics  for wind, temp. rho, ...
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: adv_now      !< progn. vars in    advection for tracers
    TYPE(t_nh_prog)       , INTENT(inout), TARGET :: adv_new      !< progn. vars after advection for tracers

    !> Local variables

    CHARACTER(len=max_timedelta_str_len) :: neg_dt_string !< negative time delta as string
    TYPE(t_timedelta)     , POINTER      :: neg_dt_mtime  !< negative time delta as mtime variable
    TYPE(t_datetime)      , POINTER      :: datetime      !< date and time at the beginning of this time step

    TYPE(t_aes_phy_field) , POINTER :: f
    TYPE(t_aes_phy_tend)  , POINTER :: t

    INTEGER               , POINTER :: trHydroMass_list(:)

    INTEGER  :: ncd          !< number of child patches
    INTEGER  :: nlev         !< number of full levels

    INTEGER  :: rls_c, rle_c
    INTEGER  :: jbs_c, jbe_c !< start and end indices for rows of cells
    INTEGER  :: jcs,jce      !< cell start and end indices

    INTEGER  :: rls_e, rle_e
    INTEGER  :: jbs_e, jbe_e !< start and end indices for rows of edges
    INTEGER  :: jes,jee      !< edge start and end indices

    INTEGER  :: jcn,jbn      !< jc and jb of neighbor cells sharing an edge je

    INTEGER  :: jg           !< grid   index
    INTEGER  :: jt           !< tracer index
    INTEGER  :: jb           !< block  index
    INTEGER  :: jk           !< level  index
    INTEGER  :: jc           !< cell   index
    INTEGER  :: je           !< edge   index

    REAL(wp) :: vn1, vn2, tend_vn_phy !< for computing dvn/dt|phy

    REAL(wp) :: inv_dt       !< 1/dt

    !=====================================================================================
    !

    ! Compute the datetime, from which this time step started.
    ! -> This assures correct events for the physics.
    ! -> Physics boundary conditions for the datetime
    !    from which the time step started.
    !
    CALL getPTStringFromSeconds(-dt, neg_dt_string)
    neg_dt_mtime => newTimedelta(neg_dt_string)
    datetime     => newDatetime(datetime_new)
    datetime     =  datetime_new + neg_dt_mtime

    ! grid index
    jg    = patch%id

    ! number of full levels on this patch
    nlev  = patch%nlev
    !
    ncd   = MAX(1,patch%n_childdom)
    !
    ! cell index ranges
    rls_c = grf_bdywidth_c+1
    rle_c = min_rlcell_int
    jbs_c = patch%cells%start_blk(rls_c,  1)
    jbe_c = patch%cells%  end_blk(rle_c,ncd)

    ! edge index ranges
    rls_e = grf_bdywidth_e+1
    rle_e = min_rledge_int
    jbs_e = patch%edges%start_blk(rls_e,  1)
    jbe_e = patch%edges%  end_blk(rle_e,ncd)
    !
    ! inverse time step
    inv_dt = 1.0_wp/dt
    !
    ! associate pointers
    f => prm_field(jg)
    t => prm_tend (jg)
    !
    trHydroMass_list => advection_config(jg)%trHydroMass%list
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! Memory references for the atmospheric state
    !
    ! dynamic pointers to the "new" time slice of the prognostic state
    f%rho      => dyn_new%rho
    f%qtrc_dyn => adv_new%tracer
    !
    ! dynamic pointers to the "now" time slice of the prognostic state,
    ! used here to provide work space for internal updating in physics
    f%qtrc_phy => adv_now%tracer
    f%ta       => dyn_now%theta_v ! in physics, otherwise => diag%temp
    f%wa       => dyn_now%w
    !
    ! static pointers to diagnostics
    ! f%pfull  => diag%pres
    ! f%phalf  => diag%pres_ifc
    ! f%ua     => diag%u
    ! f%va     => diag%v
    ! f%ta     => diag%temp       ! before and after physics
    ! f%tv     => diag%tempv
    ! f%mair   => diag%airmass_new
    ! f%mtrcvi => diag%tracer_vi
    !
    ! static pointers to metrics fields
    ! f%zhalf  => metrics%z_ifc
    ! f%zfull  => metrics%z_mc
    ! f%dzhalf => metrics%ddqz_z_full
    ! f%geom
    ! f%geoi   => metrics%
    !
    ! copies from the patch information
    ! f%clon      = p_patch%cells%center%lon
    ! f%clat      = p_patch%cells%center%lat
    ! f%areacella = p_patch%cells%area
    ! f%coriol    = p_patch%cells%f_c
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! Manage GPU memory
    !
    !$ACC DATA &
    !$ACC   PRESENT(dyn_new%vn, dyn_new%w, dyn_new%rho) &
    !$ACC   PRESENT(dyn_new%exner, dyn_new%theta_v, adv_new%tracer) &
    !$ACC   PRESENT(dyn_now%exner, dyn_now%theta_v, adv_now%tracer) &
    !$ACC   PRESENT(diag%u, diag%v, diag%temp, diag%tempv) &
    !$ACC   PRESENT(diag%exner_pr, diag%exner_dyn_incr) &
    !$ACC   PRESENT(diag%ddt_exner_phy, diag%ddt_tracer_adv) &
    !$ACC   PRESENT(int_state%c_lin_e, trHydroMass_list) &
    !$ACC   PRESENT(patch%edges%cell_idx, patch%edges%cell_blk, patch%edges%primal_normal_cell) &
    !$ACC   PRESENT(f%clon, f%clat, f%pfull) &
    !$ACC   PRESENT(f%rho, f%mair, f%ta, f%wa) &
    !$ACC   PRESENT(f%qtrc_dyn, f%qtrc_phy, f%mtrcvi) &
    !$ACC   PRESENT(t%ua_phy, t%va_phy, t%wa_phy, t%ta_phy) &
    !$ACC   PRESENT(t%qtrc_phy, t%mtrcvi_phy) &
    !$ACC   COPYIN(aes_phy_config(jg:jg))
    !
    !=====================================================================================

    IF (ltimer) CALL timer_start(timer_dyn2phy)

    !=====================================================================================
    !
    ! Prepare the input fields from the prognostic state and initialize output fields
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jt,jb,jk,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb, jbs_c, jbe_c, jcs, jce, rls_c, rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev
        DO jc = jcs, jce
          !
          ! copy temperature for in/out work in physics
          f%ta    (jc,jk,jb) = diag%temp(jc,jk,jb)
          !
          ! reset physics tendencies
          t%ua_phy(jc,jk,jb) = 0.0_wp
          t%va_phy(jc,jk,jb) = 0.0_wp
          t%ta_phy(jc,jk,jb) = 0.0_wp
          !
        END DO !jc
      END DO !jk
      !$ACC END PARALLEL
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev+1
        DO jc = jcs, jce
          !
          ! copy vertical velocity for in/out work
          f%wa    (jc,jk,jb) = dyn_new%w(jc,jk,jb)
          !
          ! reset physics tendencies
          t%wa_phy(jc,jk,jb) = 0.0_wp
          !
        END DO !jc
      END DO !jk
      !$ACC END PARALLEL
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jt = 1,ntracer
        DO jk = 1,nlev
          DO jc = jcs, jce
            !
            ! handle negative tracer mass fractions resulting from dynamics
            !
            IF (aes_phy_config(jg)%iqneg_d2p /= 0) THEN
                IF (adv_new%tracer(jc,jk,jb,jt) < 0.0_wp) THEN

                  IF (aes_phy_config(jg)%iqneg_d2p == 1 .OR. aes_phy_config(jg)%iqneg_d2p == 3) THEN
                     CALL print_value('d2p:index of grid   jg',jg)
                     CALL print_value('d2p:index of block  jb',jb)
                     CALL print_value('d2p:index of tracer jt',jt)
                     CALL print_value('d2p:index of level  jk',jk)
                     CALL print_value('d2p:index of cell   jc',jc)
                     CALL print_value('d2p:pressure      [Pa]',f%pfull(jc,jk,jb))
                     CALL print_value('d2p:longitude    [deg]',f%clon(jc,jb)*rad2deg)
                     CALL print_value('d2p:latitude     [deg]',f%clat(jc,jb)*rad2deg)
                     CALL print_value('d2p:tracer(jt) [kg/kg]',adv_new%tracer(jc,jk,jb,jt))
                  END IF

                  IF (aes_phy_config(jg)%iqneg_d2p == 2 .OR. aes_phy_config(jg)%iqneg_d2p == 3) THEN
                     diag%ddt_tracer_adv(jc,jk,jb,jt) = diag%ddt_tracer_adv(jc,jk,jb,jt) - adv_new%tracer(jc,jk,jb,jt)*inv_dt
                     adv_new%tracer(jc,jk,jb,jt) = 0.0_wp
                  END IF
               END IF
            END IF
            !
            ! copy tracer mass fraction for in/out work in physics
            f%qtrc_phy(jc,jk,jb,jt)  = adv_new%tracer(jc,jk,jb,jt)
            !
            ! reset physics tendencies
            t%qtrc_phy(jc,jk,jb,jt)  = 0.0_wp
            !
          END DO ! jc
        END DO ! jk
      END DO ! jt
      !$ACC END PARALLEL
      !
    END DO ! jb
!$OMP END DO
!$OMP END PARALLEL

    !$ACC WAIT(1)

    ! only if turbulent diffusion parameterizations are used
    IF (aes_phy_tc(jg)%dt_vdf > dt_zero) THEN
      !
      IF (ltimer) CALL timer_start(timer_d2p_sync)
      CALL sync_patch_array_mult(SYNC_E, patch, 1, dyn_new%vn)
      IF (ltimer) CALL timer_stop(timer_d2p_sync)
      !
      ! interpolate vn -> (u,v)
      CALL rbf_vec_interpol_cell(dyn_new%vn, patch, int_state,       &! in
        &                        f%ua, f%va,                         &! out
        &                        opt_rlstart=rls_c, opt_rlend=rle_c, &! in
        &                        opt_acc_async=.TRUE.)                ! in
      !
      ! only if the Smagorinsky turbulent diffusion parameterizations is used
      IF (aes_vdf_config(jg)%turb == 2) THEN
        !
        ! synchronize input fields to allow horizontal operations
        IF (ltimer) CALL timer_start(timer_d2p_sync)
        CALL sync_patch_array_mult(SYNC_C, patch, 4, f%rho, f%ua, f%va, f%wa)
        IF (ltimer) CALL timer_stop(timer_d2p_sync)
        !
      END IF
      !
    END IF
    !
    !=====================================================================================

    IF (ltimer)  CALL timer_stop (timer_dyn2phy)

    !=====================================================================================
    !
    ! Prepare the physics boundary conditions
    !
    IF (ltimer) CALL timer_start(timer_aes_bcs)
    CALL aes_phy_bcs(patch, datetime, dt)
    IF (ltimer) CALL timer_stop(timer_aes_bcs)
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! Calculate the physics tendencies
    !

    IF (aes_phy_config(jg)%ljsb) CALL jsbach_start_timestep(jg, datetime, dt)


    IF (ltimer) CALL timer_start(timer_aes_phy)
    CALL aes_phy_main(patch, datetime, dt)
    IF (ltimer) CALL timer_stop(timer_aes_phy)


    IF (aes_phy_config(jg)%ljsb) CALL jsbach_finish_timestep(jg, datetime, dt)

    !
    !=====================================================================================

    !=====================================================================================
    !
    ! Couple atmosphere and ocean, if needed
    !







    !
    !=====================================================================================

    IF (ltimer) CALL timer_start(timer_phy2dyn)

    !=====================================================================================
    !
    ! Update the prognostic state with the physics tendencies
    !
    ! 
    ! Only if turbulent diffusion parameterizations are used
    IF (aes_phy_tc(jg)%dt_vdf > dt_zero) THEN
      !
      IF (ltimer) CALL timer_start(timer_p2d_sync)
      CALL sync_patch_array_mult(SYNC_C, patch, 2, t%ua_phy, t%va_phy)
      IF (ltimer) CALL timer_stop(timer_p2d_sync)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, je, jes, jee, jcn, jbn, vn1, vn2, tend_vn_phy) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs_e,jbe_e
        !
        CALL get_indices_e(patch, jb, jbs_e, jbe_e, jes, jee, rls_e, rle_e)
        IF (jes>jee) CYCLE
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(jcn, jbn, vn1, vn2, tend_vn_phy) COLLAPSE(2)
        DO jk = 1,nlev
          DO je = jes,jee
            !
            jcn =   patch%edges%cell_idx(je,jb,1)
            jbn =   patch%edges%cell_blk(je,jb,1)
            vn1 =   t%ua_phy(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v1 &
              &   + t%va_phy(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v2
            !
            jcn =   patch%edges%cell_idx(je,jb,2)
            jbn =   patch%edges%cell_blk(je,jb,2)
            vn2 =   t%ua_phy(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v1 &
              &   + t%va_phy(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v2
            !
            tend_vn_phy = int_state%c_lin_e(je,1,jb)*vn1 + int_state%c_lin_e(je,2,jb)*vn2
            !
            ! new normal wind
            dyn_new%vn(je,jk,jb) = dyn_new%vn(je,jk,jb) + dt*tend_vn_phy
            !
            ! Set physics forcing to zero so that it is not re-applied in the dynamical core
            diag%ddt_vn_phy(je,jk,jb) = 0._vp
            !
          END DO ! je
        END DO ! jk
        !$ACC END PARALLEL
        !
      END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      !
    END IF

    ! Only if the Smagorinsky turbulent diffusion parameterizations is used
    IF (aes_phy_tc(jg)%dt_vdf > dt_zero .AND. aes_vdf_config(jg)%turb == 2) THEN
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, jcs, jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs_c,jbe_c
        !
        CALL get_indices_c(patch, jb, jbs_c, jbe_c, jcs, jce, rls_c, rle_c)
        IF (jcs>jce) CYCLE
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1,nlev
          DO jc = jcs, jce
            !
            ! new vertical wind
            dyn_new%w(jc,jk,jb) = dyn_new%w(jc,jk,jb) + dt*t%wa_phy(jc,jk,jb)
            !
          END DO ! jc
        END DO ! jk
        !$ACC END PARALLEL
        !
      END DO ! jb
!$OMP END DO
!$OMP END PARALLEL
      !
    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jt, jk, jc, jcs, jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs_c,jbe_c
      !
      CALL get_indices_c(patch, jb, jbs_c, jbe_c, jcs, jce, rls_c, rle_c)
      IF (jcs>jce) CYCLE
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jt = 1,ntracer
        DO jk = 1,nlev
          DO jc = jcs, jce
            !
            ! new tracer mass fraction
            adv_new%tracer(jc,jk,jb,jt) = adv_new%tracer(jc,jk,jb,jt) + dt*t%qtrc_phy(jc,jk,jb,jt)
            !
            ! handle negative tracer mass fractions resulting from physics
            IF (aes_phy_config(jg)%iqneg_p2d /= 0) THEN
              IF (adv_new%tracer(jc,jk,jb,jt) < 0.0_wp) THEN

                IF (aes_phy_config(jg)%iqneg_p2d == 1 .OR. aes_phy_config(jg)%iqneg_p2d == 3) THEN
                  CALL print_value('p2d:index of grid   jg',jg)
                  CALL print_value('p2d:index of block  jb',jb)
                  CALL print_value('p2d:index of tracer jt',jt)
                  CALL print_value('p2d:index of level  jk',jk)
                  CALL print_value('p2d:index of cell   jc',jc)
                  CALL print_value('p2d:pressure      [Pa]',f%pfull(jc,jk,jb))
                  CALL print_value('p2d:longitude    [deg]',f%clon(jc,jb)*rad2deg)
                  CALL print_value('p2d:latitude     [deg]',f%clat(jc,jb)*rad2deg)
                  CALL print_value('p2d:tracer(jt) [kg/kg]',adv_new%tracer(jc,jk,jb,jt))
                END IF

                IF (aes_phy_config(jg)%iqneg_p2d == 2 .OR. aes_phy_config(jg)%iqneg_p2d == 3) THEN
                  t%qtrc_phy(jc,jk,jb,jt) = t%qtrc_phy(jc,jk,jb,jt) - adv_new%tracer(jc,jk,jb,jt)*inv_dt
                  adv_new%tracer(jc,jk,jb,jt) = 0.0_wp
                END IF
              END IF
            END IF
            !
          END DO ! jc
        END DO ! jk
      END DO ! jt
      !$ACC END PARALLEL
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1,nlev
        DO jc = jcs, jce
          !
          ! new temperature
          diag%temp (jc,jk,jb) = diag%temp(jc,jk,jb) + dt*t%ta_phy(jc,jk,jb)
          !
          ! new virtual temperature
          diag%tempv(jc,jk,jb) = diag%temp(jc,jk,jb)*(1._wp + vtmpc1*adv_new%tracer(jc,jk,jb,iqv)      &
            &                                         - SUM(adv_new%tracer(jc,jk,jb,trHydroMass_list)))
          !
          ! save provisional "new" exner from the slow-physics-forced dynamics
          dyn_now%exner(jc,jk,jb) = dyn_new%exner(jc,jk,jb)
          !
          ! new exner
          dyn_new%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd/p0ref*dyn_new%rho(jc,jk,jb)*diag%tempv(jc,jk,jb)))
          !
          ! add exner change from fast physics to exner_pr in order to avoid unphysical sound wave generation
          diag%exner_pr(jc,jk,jb) = diag%exner_pr(jc,jk,jb) + dyn_new%exner(jc,jk,jb) - dyn_now%exner(jc,jk,jb)
          !
          dyn_new%theta_v(jc,jk,jb) = diag%tempv(jc,jk,jb)/dyn_new%exner(jc,jk,jb)
          !
          ! set physics forcing to zero so that it is not re-applied in the dynamical core
          diag%ddt_exner_phy(jc,jk,jb) = 0._vp
          !
          ! set the dynamical exner increment to zero, accumulated in solve_nh, must be zero for next step
          diag%exner_dyn_incr(jc,jk,jb) = 0._vp
          !
        END DO ! jc
      END DO ! jk
      !$ACC END PARALLEL
      !
    END DO !jb
!$OMP END DO
!$OMP END PARALLEL
    !
    !=====================================================================================

    !=====================================================================================
    !
    ! Synchronize the prognostic state
    !
    IF (ltimer) CALL timer_start(timer_p2d_sync)
    CALL sync_patch_array_mult(SYNC_E, patch, 1,         &
      &                        dyn_new%vn)
    CALL sync_patch_array_mult(SYNC_C, patch, 5+ntracer, &
      &                        dyn_new%w,                &
      &                        dyn_new%rho,              &
      &                        dyn_new%exner,            &
      &                        dyn_new%theta_v,          &
      &                        diag%exner_pr,            &
      &                        f4din=adv_new%tracer)
    IF (ltimer) CALL timer_stop(timer_p2d_sync)
    !
    !$ACC WAIT(1)
    !$ACC END DATA
    !
    !=====================================================================================

    ! physics is finished: point 'ta' back to 'temp'
    f%ta => diag%temp

    NULLIFY(f)
    NULLIFY(t)

    CALL deallocateDatetime(datetime)
    CALL deallocateTimedelta(neg_dt_mtime)

    IF (ltimer) CALL timer_stop(timer_phy2dyn)

    !=====================================================================================
    !
    ! Now the final new state is ready.
    !
    !=====================================================================================

  END SUBROUTINE interface_iconam_aes
  !----------------------------------------------------------------------------

END MODULE mo_interface_iconam_aes
