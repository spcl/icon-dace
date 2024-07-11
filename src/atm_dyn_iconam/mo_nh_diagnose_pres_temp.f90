! Contains utilities for diagnose pressure, temperature and air mass in nh model
!
!
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
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_diagnose_pres_temp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog
  USE mo_run_config,          ONLY: iqv, iforcing, lforcing
  USE mo_impl_constants,      ONLY: min_rlcell, iheldsuarez
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_physical_constants,  ONLY: rd, grav, vtmpc1, p0ref, rd_o_cpd
  USE mo_timer,               ONLY: timers_level, timer_start, timer_stop, timer_diagnose_pres_temp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_advection_config,    ONLY: advection_config
  USE mo_dynamics_config,     ONLY: ldeepatmo
  USE mo_grid_config,         ONLY: grid_sphere_radius
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE

  REAL(wp), PARAMETER :: cpd_o_rd  = 1._wp / rd_o_cpd
  REAL(wp), PARAMETER :: grav_o_rd = grav / rd

  PUBLIC :: diagnose_pres_temp
  PUBLIC :: diag_temp
  PUBLIC :: diag_pres
  PUBLIC :: calc_qsum
  PUBLIC :: compute_airmass

  CONTAINS

  ! moved here from mo_nh_stepping to avoid circular dependencies
  !>
  !! diagnose_pres_temp
  !!
  !! Diagnoses pressure and temperature from NH prognostic fields
  !!
  SUBROUTINE diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, pt_diag, pt_patch, &
    &                            opt_calc_temp, opt_calc_pres, opt_calc_temp_ifc,    &
    &                            lnd_prog, opt_slev, opt_rlend                       )


!!$    CHARACTER(len=*), PARAMETER ::  &
!!$      &  routine = 'mo_nh_diagnose_pres_temp:diagnose_pres_temp' 

    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: pt_prog      !!the prognostic variables
    TYPE(t_nh_prog),    INTENT(IN)    :: pt_prog_rcf  !!the prognostic variables which are
                                                      !! treated with reduced calling frequency

    TYPE(t_lnd_prog),   INTENT(IN), OPTIONAL :: lnd_prog 
    TYPE(t_nh_diag),    INTENT(INOUT) :: pt_diag      !!the diagnostic variables


    TYPE(t_patch),      INTENT(IN)    :: pt_patch    ! Patch

    LOGICAL, INTENT(IN), OPTIONAL   :: opt_calc_temp, opt_calc_pres, opt_calc_temp_ifc

    INTEGER, INTENT(IN), OPTIONAL :: opt_slev, opt_rlend 

    INTEGER  :: jb,jk,jc,jg
    INTEGER  :: nlev, nlevp1              !< number of full levels
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: slev, slev_moist

    LOGICAL  :: l_opt_calc_temp, l_opt_calc_pres, l_opt_calc_temp_ifc


    IF (timers_level > 8) CALL timer_start(timer_diagnose_pres_temp)

    
    ! Check for optional arguments

    IF ( PRESENT(opt_calc_temp_ifc ) ) THEN
      l_opt_calc_temp_ifc = opt_calc_temp_ifc
    ELSE
      l_opt_calc_temp_ifc = .FALSE.
    ENDIF

    IF ( PRESENT(opt_calc_temp ) ) THEN
      l_opt_calc_temp = opt_calc_temp
    ELSE
      l_opt_calc_temp = .TRUE.
    ENDIF

    IF ( PRESENT(opt_calc_pres ) ) THEN
      l_opt_calc_pres = opt_calc_pres
    ELSE
      l_opt_calc_pres = .TRUE.
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell
    ENDIF

    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    ENDIF

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    ! Nest boundaries are always included
    i_rlstart = 1

    jg = pt_patch%id
    ! start index for moisture variables other than QV
    slev_moist = MAX(kstart_moist(jg),slev)

    i_startblk = pt_patch%cells%start_block(i_rlstart)
    i_endblk   = pt_patch%cells%end_block(i_rlend)


    !$ACC DATA PRESENT(advection_config(jg)%trHydroMass%list) IF(i_am_accel_node)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( pt_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      IF ( l_opt_calc_temp) THEN
        IF ( lforcing .AND. iforcing /= iheldsuarez  ) THEN

          CALL diag_temp (pt_prog, pt_prog_rcf, advection_config(jg)%trHydroMass%list, &
            &            pt_diag, jb, i_startidx, i_endidx, slev, slev_moist, nlev)


        ELSE ! .NOT. lforcing or Held-Suarez test forcing

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = slev, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
               pt_diag%tempv(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) * pt_prog%exner(jc,jk,jb)
               pt_diag%temp(jc,jk,jb)  = pt_diag%tempv  (jc,jk,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL

        ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !> diagnose temperature on interface levels
      !-------------------------------------------------------------------------
      
      IF ( l_opt_calc_temp_ifc ) THEN
        
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = MAX(slev+1,2), nlev
!DIR$ IVDEP
          DO jc =  i_startidx, i_endidx
            pt_diag%temp_ifc(jc,jk,jb) = &
              p_metrics%wgtfac_c(jc,jk,jb)*pt_diag%temp(jc,jk,jb) +      &
              (1._wp-p_metrics%wgtfac_c(jc,jk,jb))*pt_diag%temp(jc,jk-1,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        IF ( PRESENT(lnd_prog) ) THEN
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG VECTOR
          DO jc =  i_startidx, i_endidx
            pt_diag%temp_ifc(jc,     1,jb) = pt_diag%temp (jc,1,jb)
            pt_diag%temp_ifc(jc,nlevp1,jb) = lnd_prog%t_g (jc,jb)
          ENDDO
          !$ACC END PARALLEL
        ELSE
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG VECTOR
          DO jc =  i_startidx, i_endidx
            pt_diag%temp_ifc(jc,     1,jb) = pt_diag%temp (jc,1,jb)
            pt_diag%temp_ifc(jc,nlevp1,jb) = pt_diag%temp (jc,nlev,jb)
          ENDDO
          !$ACC END PARALLEL
        ENDIF

      ENDIF !l_opt_calc_temp_ifc

      !-------------------------------------------------------------------------
      !> diagnose pressure on main and interface levels
      !!    and   pressure thickness
      !!
      !-------------------------------------------------------------------------

      IF ( l_opt_calc_pres ) CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, slev, nlev)

    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA
    
    IF (timers_level > 8) CALL timer_stop(timer_diagnose_pres_temp)

  END SUBROUTINE diagnose_pres_temp


  !!
  !! Reduced version for pressure diagnosis to be called from within a block loop
  !! Diagnoses
  !! - pressure on full levels
  !! - pressure on half levels
  !! - pressure thickness
  !! - surface pressure
  !! Note that the pressure is diagnosed by vertical integration of the 
  !! hydrostatic equation!
  !!
  SUBROUTINE diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, slev, nlev)


    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: pt_prog      !!the prognostic variables
 
    TYPE(t_nh_diag),    INTENT(INOUT) :: pt_diag      !!the diagnostic variables


    INTEGER, INTENT(IN) :: jb, i_startidx, i_endidx, slev, nlev

    INTEGER  :: jk,jc

    REAL(wp) :: dz1, dz2, dz3

    !$ACC DATA PRESENT(pt_diag%tempv, pt_diag%pres_sfc, pt_diag%pres_ifc, pt_diag%pres) &
    !$ACC   PRESENT(pt_diag%dpres_mc, pt_prog%theta_v, pt_prog%exner, p_metrics%ddqz_z_full) &
    !$ACC   IF(i_am_accel_node)

    IF (.NOT. ldeepatmo) THEN

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR PRIVATE(dz1, dz2, dz3)
      DO jc = i_startidx, i_endidx
        ! Height differences between surface and third-lowest main level
        dz1 = p_metrics%ddqz_z_full(jc,nlev,jb)
        dz2 = p_metrics%ddqz_z_full(jc,nlev-1,jb)
        dz3 = 0.5_wp*p_metrics%ddqz_z_full(jc,nlev-2,jb)
        
        ! Compute surface pressure starting from third-lowest level; this is done
        ! in order to avoid contamination by sound-wave activity in the presence of strong latent heating
        pt_diag%pres_sfc(jc,jb) = p0ref * EXP( cpd_o_rd*LOG(pt_prog%exner(jc,nlev-2,jb)) + &
          grav_o_rd*(dz1/pt_diag%tempv(jc,nlev,jb) + dz2/pt_diag%tempv(jc,nlev-1,jb) +     &
          dz3/pt_diag%tempv(jc,nlev-2,jb)) )
        
        pt_diag%pres_ifc(jc,nlev+1,jb) = pt_diag%pres_sfc(jc,jb)
      ENDDO
      !$ACC END PARALLEL

      
      !-------------------------------------------------------------------------
      !> diagnose pressure for physics parameterizations
      !! this is accomplished by vertical integration of the hydrostatic equation
      !! because the physics schemes actually need the air mass represented 
      !! by a given model layer
      !-------------------------------------------------------------------------
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP SEQ
      DO jk = nlev, slev,-1
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          
          ! pressure at interface levels
          pt_diag%pres_ifc(jc,jk,jb) = pt_diag%pres_ifc(jc,jk+1,jb)                  &
            & *EXP(-grav_o_rd*p_metrics%ddqz_z_full(jc,jk,jb)/pt_diag%tempv(jc,jk,jb))
          
          ! pressure at main levels
          pt_diag%pres(jc,jk,jb) = SQRT(pt_diag%pres_ifc(jc,jk,jb) * &
                                        pt_diag%pres_ifc(jc,jk+1,jb) )
          
          ! layer thickness with respect to pressure
          pt_diag%dpres_mc(jc,jk,jb) = pt_diag%pres_ifc(jc,jk+1,jb) &
                                     - pt_diag%pres_ifc(jc,jk  ,jb)
          
        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ELSE

      CALL diag_pres_deepatmo(z_mc=p_metrics%z_mc(:,:,jb), z_ifc=p_metrics%z_ifc(:,:,jb),                  &
        & exner=pt_prog%exner(:,:,jb), tempv=pt_diag%tempv(:,:,jb), pres_sfc=pt_diag%pres_sfc(:,jb),       & 
        & pres_ifc=pt_diag%pres_ifc(:,:,jb), pres=pt_diag%pres(:,:,jb), dpres_mc=pt_diag%dpres_mc(:,:,jb), &
        & start_indices=[i_startidx,slev], end_indices=[i_endidx,nlev],                                    &
#ifdef _OPENACC
        & lacc=i_am_accel_node)
#else
        & lacc=.FALSE.)
#endif

    ENDIF ! IF (.NOT. ldeepatmo)

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE diag_pres



  !!
  !! Reduced version for temperature diagnosis to be called from within a block loop
  !! Diagnoses 
  !! - virtual temperature
  !! - temperature
  !!
  !!
  SUBROUTINE diag_temp (pt_prog, pt_prog_rcf, condensate_list, pt_diag, &
    &                   jb, i_startidx, i_endidx, slev, slev_moist, nlev)


    TYPE(t_nh_prog),    INTENT(IN)    :: pt_prog       !!the prognostic variables
    TYPE(t_nh_prog),    INTENT(IN)    :: pt_prog_rcf   !!the prognostic variables which are
                                                       !! treated with reduced calling frequency
    INTEGER        ,    INTENT(IN)    :: &             !! IDs of all tracers containing 
      &  condensate_list(:)                            !! prognostic condensate. Required for
                                                       !! computing the water loading term.  
    TYPE(t_nh_diag),    INTENT(INOUT) :: pt_diag       !!the diagnostic variables


    INTEGER, INTENT(IN) :: jb, i_startidx, i_endidx, slev, slev_moist, nlev 

    INTEGER  :: jk,jc
    REAL(wp) :: z_qsum(nproma,nlev)

    !$ACC DATA PRESENT(pt_prog_rcf, pt_diag, pt_prog) &
    !$ACC   CREATE(z_qsum) &
    !$ACC   IF(i_am_accel_node)

    CALL calc_qsum (pt_prog_rcf%tracer, z_qsum, condensate_list, jb, i_startidx, i_endidx, slev, slev_moist, nlev)

    !$ACC PARALLEL DEFAULT(PRESENT) ATTACH(pt_prog_rcf%tracer) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev, nlev
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
        pt_diag%tempv(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) * pt_prog%exner(jc,jk,jb)
        pt_diag%temp(jc,jk,jb)  = pt_diag%tempv(jc,jk,jb) /            &
          ( 1._wp+vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)-z_qsum(jc,jk) )
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE diag_temp


  !!
  !! Calculates the sum of selected tracer mass fractions
  !! Extracted from diag_temp (see above) in order to encapsulate the code duplication needed for vectorization
  !!
  SUBROUTINE calc_qsum (tracer, qsum, tracer_list, jb, i_startidx, i_endidx, slev, slev_moist, nlev)

    REAL(wp), INTENT(IN)    :: tracer(:,:,:,:)       !! tracer array
    REAL(wp), INTENT(INOUT) :: qsum(:,:)             !! output: sum of condensates [kg kg-1]
    INTEGER,  INTENT(IN)    :: tracer_list(:)        !! IDs of tracers for which the sum is taken
    INTEGER,  INTENT(IN)    :: jb, i_startidx, i_endidx, slev, slev_moist, nlev

    INTEGER  :: jk,jc
#ifdef __SX__
    INTEGER  :: jl, jt
#endif


#ifdef __SX__
    qsum(:,MIN(slev,slev_moist):nlev) = 0._wp
    DO jl = 1, SIZE(tracer_list)
      jt = tracer_list(jl)
      DO jk = slev_moist, nlev
        DO jc = i_startidx, i_endidx
          qsum(jc,jk) = qsum(jc,jk) + tracer(jc,jk,jb,jt)
        ENDDO
      ENDDO
    ENDDO
#else
    !$ACC DATA NO_CREATE(qsum, tracer, tracer_list)

    IF (slev < slev_moist) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, slev_moist-1
        DO jc = i_startidx, i_endidx
          qsum(jc,jk) = 0._wp
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = slev_moist, nlev
      DO jc = i_startidx, i_endidx
        qsum(jc,jk) = SUM(tracer(jc,jk,jb,tracer_list))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !DA: wait here is due to mo_nh_interface_nwp
    !$ACC WAIT
    !$ACC END DATA
#endif


  END SUBROUTINE calc_qsum


  !>
  !! Compute air mass within grid cell
  !!
  !! Compute air mass within grid cell. Note that here, the air mass is defined
  !! as \rho*\Delta z [kg m-2]. Computing the true grid cell air mass 
  !! requires an additional multiplication with the grid cell area.
  !!
  SUBROUTINE compute_airmass (p_patch, p_metrics, rho, airmass)

    TYPE(t_patch),      INTENT(IN   ) :: p_patch
    TYPE(t_nh_metrics), INTENT(IN   ) :: p_metrics
    REAL(wp),           INTENT(IN   ) :: rho(:,:,:)      ! air density [kg m-3]
    REAL(wp),           INTENT(INOUT) :: airmass(:,:,:)  ! air mass    [kg m-2]

    INTEGER :: nlev                  ! number of vertical levels
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jk,jb
  !---------------------------------------------------------!

    ! number of vertical levels
    nlev = p_patch%nlev

    ! halo points must be included !
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)


    !$ACC DATA PRESENT(rho, airmass, p_metrics%ddqz_z_full, p_metrics%deepatmo_vol_mc) IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          airmass(jc,jk,jb) = rho(jc,jk,jb)*p_metrics%ddqz_z_full(jc,jk,jb)*p_metrics%deepatmo_vol_mc(jk)
        ENDDO  ! jc
      ENDDO  ! jk
      !$ACC END PARALLEL

    ENDDO ! jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE compute_airmass


  !>
  !! Deep-atmosphere variant of diag_pres (for use within block loop)
  !!
  SUBROUTINE diag_pres_deepatmo( z_mc,          & !in
    &                            z_ifc,         & !in
    &                            exner,         & !in
    &                            tempv,         & !in
    &                            pres_sfc,      & !inout
    &                            pres_ifc,      & !inout
    &                            pres,          & !inout
    &                            dpres_mc,      & !inout
    &                            start_indices, & !in
    &                            end_indices,   & !in
    &                            lacc)            !optin

    ! In/out variables
    REAL(wp), INTENT(IN)              :: z_mc(:,:)        ! Height of cell centers
    REAL(wp), INTENT(IN)              :: z_ifc(:,:)       ! Height of cell interfaces
    REAL(wp), INTENT(IN)              :: exner(:,:)       ! Exner pressure
    REAL(wp), INTENT(IN)              :: tempv(:,:)       ! Virtual temperature
    REAL(wp), INTENT(INOUT)           :: pres_sfc(:)      ! Surface pressure
    REAL(wp), INTENT(INOUT)           :: pres_ifc(:,:)    ! Pressure at cell interfaces
    REAL(wp), INTENT(INOUT)           :: pres(:,:)        ! Pressure at cell centers
    REAL(wp), INTENT(INOUT)           :: dpres_mc(:,:)    ! Pressure difference between lower and upper cell interfaces
    INTEGER,  INTENT(IN)              :: start_indices(2) ! Start indices: [i_startidx,slev]
    INTEGER,  INTENT(IN)              :: end_indices(2)   ! End indices:   [i_endidx,  nlev]
    LOGICAL,  INTENT(IN),    OPTIONAL :: lacc             ! Optional flag for use of OpenACC

    ! Local variables
    REAL(wp) :: inv_grid_sphere_radius
    REAL(wp) :: zgpot_ifc_nlevp1, zgpot_ifc_nlev, zgpot_ifc_nlevm1, zgpot_ifc_nlevm2
    REAL(wp) :: dzgpot_mc_nlev, dzgpot_mc_nlevm1, dzgpot_mc_nlevm2_halved
    REAL(wp) :: zgpot_mc, dzgpot_mc, zgpot_lifc, zgpot_uifc
    INTEGER  :: jk, jc
    INTEGER  :: nlev
    LOGICAL  :: lzacc

    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! Rudimentary consistency checks
    IF (ANY(start_indices(:) < 1) .OR. ANY(end_indices(:) < 1)) THEN
      RETURN
    ELSEIF (ANY(start_indices(:) > end_indices(:))) THEN
      RETURN
    ELSEIF (ANY(end_indices(:) > SHAPE(z_ifc))) THEN
      RETURN
    ENDIF

    nlev = end_indices(2)

    ! Inverse Earth radius
    inv_grid_sphere_radius = 1._wp / grid_sphere_radius

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR &
    !$ACC   PRIVATE(zgpot_ifc_nlevp1, zgpot_ifc_nlev, zgpot_ifc_nlevm1, zgpot_ifc_nlevm2) &
    !$ACC   PRIVATE(dzgpot_mc_nlev, dzgpot_mc_nlevm1, dzgpot_mc_nlevm2_halved)
    DO jc = start_indices(1), end_indices(1)

      ! Auxiliary geopotential heights
      zgpot_ifc_nlevp1 = z_ifc(jc,nlev+1) / (1._wp + inv_grid_sphere_radius * z_ifc(jc,nlev+1))
      zgpot_ifc_nlev   = z_ifc(jc,nlev)   / (1._wp + inv_grid_sphere_radius * z_ifc(jc,nlev))
      zgpot_ifc_nlevm1 = z_ifc(jc,nlev-1) / (1._wp + inv_grid_sphere_radius * z_ifc(jc,nlev-1))
      zgpot_ifc_nlevm2 = z_ifc(jc,nlev-2) / (1._wp + inv_grid_sphere_radius * z_ifc(jc,nlev-2))

      ! Auxiliary geopotential height differences
      dzgpot_mc_nlev          = zgpot_ifc_nlev   - zgpot_ifc_nlevp1
      dzgpot_mc_nlevm1        = zgpot_ifc_nlevm1 - zgpot_ifc_nlev
      dzgpot_mc_nlevm2_halved = 0.5_wp * (zgpot_ifc_nlevm2 - zgpot_ifc_nlevm1)

      ! Surface pressure
      pres_sfc(jc) = p0ref * EXP(cpd_o_rd * LOG(exner(jc,nlev-2)) + grav_o_rd * (dzgpot_mc_nlev / tempv(jc,nlev) &
        &          + dzgpot_mc_nlevm1 / tempv(jc,nlev-1) + dzgpot_mc_nlevm2_halved / tempv(jc,nlev-2)))

      ! Start value for vertical integration of hydrostatic balance
      pres_ifc(jc,nlev+1) = pres_sfc(jc)

    ENDDO  !jc
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO jk = end_indices(2), start_indices(2), -1 
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR &
      !$ACC   PRIVATE(zgpot_mc, zgpot_lifc, zgpot_uifc, dzgpot_mc)
      DO jc = start_indices(1), end_indices(1)

        ! Auxiliary geopotential heights and height differences
        zgpot_mc   = z_mc(jc,jk)    / (1._wp + inv_grid_sphere_radius * z_mc(jc,jk))
        zgpot_lifc = z_ifc(jc,jk+1) / (1._wp + inv_grid_sphere_radius * z_ifc(jc,jk+1))
        zgpot_uifc = z_ifc(jc,jk)   / (1._wp + inv_grid_sphere_radius * z_ifc(jc,jk))
        dzgpot_mc  = zgpot_uifc - zgpot_lifc

        ! Integrate hydrostatic balance: cell interfaces
        pres_ifc(jc,jk) = pres_ifc(jc,jk+1) * EXP(-grav_o_rd * dzgpot_mc / tempv(jc,jk))

        ! Integrate hydrostatic balance: cell centers
        ! (since gpot(0.5*(z1+z2)) /= 0.5*(gpot(z1)+gpot(z2)), where gpot(z)=zgpot, 
        ! the pressure at main levels cannot be computed from the geometric mean, 
        ! so we have to integrate again)
        pres(jc,jk) = pres_ifc(jc,jk+1) * EXP(-grav_o_rd * (zgpot_mc - zgpot_lifc) / tempv(jc,jk))

        dpres_mc(jc,jk) = pres_ifc(jc,jk+1) - pres_ifc(jc,jk)

      ENDDO  !jc
    ENDDO  !jk
    !$ACC END PARALLEL

  END SUBROUTINE diag_pres_deepatmo

END MODULE  mo_nh_diagnose_pres_temp

