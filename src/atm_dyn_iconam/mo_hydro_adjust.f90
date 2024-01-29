!
! Routines for computing hydrostatically balanced initial conditions
!
! A hydrostatically balanced state is achieved by integrating the discretized (!)
! third equation of motion, assuming dw/dt=0.
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
#include "omp_definitions.inc"
!----------------------------

MODULE mo_hydro_adjust

  USE mo_kind,                  ONLY: wp
  USE mo_model_domain,          ONLY: t_patch
  USE mo_exception,             ONLY: finish
  USE mo_impl_constants,        ONLY: SUCCESS, min_rlcell
  USE mo_physical_constants,    ONLY: rd, rdv, cpd, cvd_o_rd, p0ref, vtmpc1, o_m_rdv
  USE mo_parallel_config,       ONLY: nproma
  USE mo_run_config,            ONLY: msg_level
  USE mo_nonhydro_types,        ONLY: t_nh_metrics
  USE mo_loopindices,           ONLY: get_indices_c
  USE mo_satad,                 ONLY: sat_pres_water
  USE mo_util_string,           ONLY: int2string, real2string
  USE mo_util_table,            ONLY: t_table, initialize_table, add_table_column, &
    &                                 set_table_entry, print_table, finalize_table
  USE mo_sync,                  ONLY: global_max
  USE mo_mpi,                   ONLY: my_process_is_stdio, p_pe_work_only

  IMPLICIT NONE

  PRIVATE

  ! subroutines
  !
  PUBLIC :: hydro_adjust
  PUBLIC :: hydro_adjust_const_thetav
  PUBLIC :: hydro_adjust_iterative

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nh_init_utils'

CONTAINS

  !------------------------------------------------------------------------------
  !>
  !! SUBROUTINE hydro_adjust
  !! Computes hydrostatically balanced initial condition by bottom-up integration
  !! Virtual temperature is kept constant during the adjustment process
  !!
  !! Input/Output: density, Exner pressure, virtual potential temperature
  !!
  SUBROUTINE hydro_adjust(p_patch, p_nh_metrics, rho, exner, theta_v )


    TYPE(t_patch),      INTENT(IN)       :: p_patch
    TYPE(t_nh_metrics), INTENT(IN)       :: p_nh_metrics

    ! Thermodynamic fields - all defined at full model levels
    REAL(wp), INTENT(INOUT) :: rho(:,:,:)        ! density (kg/m**3)
    REAL(wp), INTENT(INOUT) :: exner(:,:,:)      ! Exner pressure
    REAL(wp), INTENT(INOUT) :: theta_v(:,:,:)    ! virtual potential temperature (K)


    ! LOCAL VARIABLES
    REAL(wp) :: temp_v(nproma,p_patch%nlev) ! virtual temperature
    REAL(wp), DIMENSION(nproma) :: z_fac1, z_fac2, z_fac3, za, zb, zc

    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,temp_v,z_fac1,z_fac2,z_fac3,za,zb,zc) ICON_OMP_DEFAULT_SCHEDULE

    ! The full model grid including the lateral boundary interpolation zone of
    ! nested domains and MPI-halo points is processed; depending on the setup
    ! of the parallel-read routine, the input fields may need to be synchronized
    ! before entering this routine.

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      ! Compute virtual temperature
      DO jk = 1, nlev
        DO jc = 1, nlen
          temp_v(jc,jk) = theta_v(jc,jk,jb)*exner(jc,jk,jb)
        ENDDO
      ENDDO

      ! Now compute hydrostatically balanced prognostic fields:
      ! The following expressions are derived from the discretized (!) third
      ! equation of motion, assuming dw/dt = 0, and solved for the exner pressure.
      ! Because the vertical discretization differs between the triangular and
      ! hexagonal NH cores, a case discrimination is needed here
      DO jk = nlev-1, 1, -1
        DO jc = 1, nlen
          z_fac1(jc) = p_nh_metrics%wgtfac_c(jc,jk+1,jb)*(temp_v(jc,jk+1) &
            - p_nh_metrics%theta_ref_mc(jc,jk+1,jb)*exner(jc,jk+1,jb))    &
            - (1._wp-p_nh_metrics%wgtfac_c(jc,jk+1,jb))                   &
            * p_nh_metrics%theta_ref_mc(jc,jk,jb)*exner(jc,jk+1,jb)

          z_fac2(jc) = (1._wp-p_nh_metrics%wgtfac_c(jc,jk+1,jb))*temp_v(jc,jk) &
            *exner(jc,jk+1,jb)

          z_fac3(jc) = p_nh_metrics%exner_ref_mc(jc,jk+1,jb)     &
            -p_nh_metrics%exner_ref_mc(jc,jk,jb)-exner(jc,jk+1,jb)

          za(jc) = (p_nh_metrics%theta_ref_ic(jc,jk+1,jb)                     &
            *exner(jc,jk+1,jb)+z_fac1(jc))/p_nh_metrics%ddqz_z_half(jc,jk+1,jb)

          zb(jc) = -(za(jc)*z_fac3(jc)+z_fac2(jc)/p_nh_metrics%ddqz_z_half(jc,jk+1,jb) &
            + z_fac1(jc)*p_nh_metrics%d_exner_dz_ref_ic(jc,jk+1,jb))

          zc(jc) = -(z_fac2(jc)*z_fac3(jc)/p_nh_metrics%ddqz_z_half(jc,jk+1,jb) &
            + z_fac2(jc)*p_nh_metrics%d_exner_dz_ref_ic(jc,jk+1,jb))
        ENDDO !jc

        DO jc = 1, nlen
          exner(jc,jk,jb)      = (zb(jc)+SQRT(zb(jc)**2+4._wp*za(jc)*zc(jc)))/(2._wp*za(jc))
          theta_v(jc,jk,jb)    = temp_v(jc,jk)/exner(jc,jk,jb)
          rho(jc,jk,jb)        = exner(jc,jk,jb)**cvd_o_rd*p0ref/(rd*theta_v(jc,jk,jb))
        ENDDO

      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE hydro_adjust


  !------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE hydro_adjust_const_thetav
  !! Computes hydrostatically balanced initial condition by either top-down or bottom-up 
  !! integration.
  !! In contrast to the routine hydro_adjust, virtual potential temperature is kept constant
  !! during the adjustment, leading to a simpler formula
  !!
  !! Input/Output: density, Exner pressure, virtual potential temperature
  !!
  SUBROUTINE hydro_adjust_const_thetav(p_patch, p_nh_metrics, lintegrate_topdown, rho, exner, theta_v)


    TYPE(t_patch),      INTENT(IN)       :: p_patch
    TYPE(t_nh_metrics), INTENT(IN)       :: p_nh_metrics
    LOGICAL,            INTENT(IN)       :: lintegrate_topdown   ! TRUE : top-down integration for exner
                                                                 ! FALSE: bottom-up integration for exner

    ! Thermodynamic fields - all defined at full model levels
    REAL(wp), INTENT(INOUT) :: rho(:,:,:)        ! density (kg/m**3)
    REAL(wp), INTENT(INOUT) :: exner(:,:,:)      ! Exner pressure
    REAL(wp), INTENT(INOUT) :: theta_v(:,:,:)    ! virtual potential temperature (K)


    ! LOCAL VARIABLES
    REAL(wp), DIMENSION(nproma) :: theta_v_pr_ic

    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,theta_v_pr_ic) ICON_OMP_DEFAULT_SCHEDULE

    ! The full model grid including the lateral boundary interpolation zone of
    ! nested domains and MPI-halo points is processed; depending on the setup
    ! of the parallel-read routine, the input fields may need to be synchronized
    ! before entering this routine.

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      ! Now compute hydrostatically balanced prognostic fields:
      ! The following expressions are derived from the discretized (!) third
      ! equation of motion, assuming dw/dt = 0, and solved for the exner pressure.
      !
      IF (lintegrate_topdown) THEN
        ! top-down integration
        !
        DO jk = 2, nlev
          DO jc = 1, nlen
            theta_v_pr_ic(jc) = p_nh_metrics%wgtfac_c(jc,jk,jb) *        &
             (theta_v(jc,jk,jb) - p_nh_metrics%theta_ref_mc(jc,jk,jb)) + &
             (1._wp-p_nh_metrics%wgtfac_c(jc,jk,jb)) *                   &
             (theta_v(jc,jk-1,jb)-p_nh_metrics%theta_ref_mc(jc,jk-1,jb)  )

            exner(jc,jk,jb) = exner(jc,jk-1,jb) + p_nh_metrics%exner_ref_mc(jc,jk,jb) -      &
              p_nh_metrics%exner_ref_mc(jc,jk-1,jb) + p_nh_metrics%ddqz_z_half(jc,jk,jb)*    &
              theta_v_pr_ic(jc)*p_nh_metrics%d_exner_dz_ref_ic(jc,jk,jb)/(theta_v_pr_ic(jc)+ &
              p_nh_metrics%theta_ref_ic(jc,jk,jb))

          ENDDO
        ENDDO

      ELSE
        ! bottom-up integration
        !
        DO jk = nlev-1,1,-1
          DO jc = 1, nlen

            theta_v_pr_ic(jc) = p_nh_metrics%wgtfac_c(jc,jk+1,jb) *        &
             (theta_v(jc,jk+1,jb) - p_nh_metrics%theta_ref_mc(jc,jk+1,jb)) + &
             (1._wp-p_nh_metrics%wgtfac_c(jc,jk+1,jb)) *                   &
             (theta_v(jc,jk,jb)-p_nh_metrics%theta_ref_mc(jc,jk,jb)  )

            exner(jc,jk,jb) = exner(jc,jk+1,jb)                                                &
              &  + (p_nh_metrics%exner_ref_mc(jc,jk,jb) - p_nh_metrics%exner_ref_mc(jc,jk+1,jb)) &
              &  - p_nh_metrics%ddqz_z_half(jc,jk+1,jb)/(theta_v_pr_ic(jc)+p_nh_metrics%theta_ref_ic(jc,jk+1,jb)) &
              &  * theta_v_pr_ic(jc) * p_nh_metrics%d_exner_dz_ref_ic(jc,jk+1,jb)

          ENDDO
        ENDDO

      ENDIF  ! lintegrate_topdown

      DO jk = 1, nlev
        DO jc = 1, nlen
          rho(jc,jk,jb) = exner(jc,jk,jb)**cvd_o_rd*p0ref/(rd*theta_v(jc,jk,jb))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE hydro_adjust_const_thetav



  !--------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE hydro_adjust_iterative
  !! Compute hydrostatically balanced initial state in terms of
  !!   theta_v, exner, rho, qv
  !! by integrating the discretized (!) third equation of motion, assuming dw/dt = 0
  !!
  !! Input:  
  !!   - prescribed profile of t and rel_hum
  !!   - reference state exner_ref, theta_v_ref
  !!   - boundary condition for exner
  !! If the integration is prescribed bottom-up (top-down) a lower (upper) 
  !! boundary condition for exner must be provided. 
  !!
  !! Output: exner, theta_v, rho, qv
  !!
  SUBROUTINE hydro_adjust_iterative(p_patch, p_nh_metrics, temp_ini, rh_ini,   &
    &                               exner, theta_v, rho, qv, luse_exner_fg,    &
    &                               opt_exner_lbc, opt_exner_ubc)

    CHARACTER(len=*), PARAMETER :: &
       &  routine = modname//':hydro_adjust_iterative'

    TYPE(t_patch),      INTENT(IN)   :: p_patch
    TYPE(t_nh_metrics), INTENT(IN)   :: p_nh_metrics

    REAL(wp),           INTENT(IN)   :: temp_ini(:,:,:) ! prescribed temperature profile

    REAL(wp),           INTENT(IN)   :: rh_ini(:,:,:)   ! prescribed relative humidity profile

    REAL(wp), TARGET,   INTENT(INOUT):: exner(:,:,:)    ! Exner pressure
                                                        ! exner(:,nlev,:) is taken as lower boundary condition
    REAL(wp),           INTENT(INOUT):: theta_v(:,:,:)  ! virtual potential temperature (K)
    REAL(wp),           INTENT(INOUT):: rho(:,:,:)      ! density (kg/m**3)
    REAL(wp),           INTENT(INOUT):: qv(:,:,:)       ! water vapour mass fraction (kg/kg)

    LOGICAL,            INTENT(IN)   :: luse_exner_fg   ! Determines which exner field will be used as first guess
                                                        ! TRUE: use exner field provided via argument list
                                                        ! FALSE: use exner reference profile

    REAL(wp), OPTIONAL, INTENT(IN)   :: opt_exner_ubc(:,:)! exner pressure upper boundary condiction (uppermost full level)
                                                          ! dim: (nproma, nblks_c)
    REAL(wp), OPTIONAL, INTENT(IN)   :: opt_exner_lbc(:,:)! exner pressure lower boundary condition (lowermost full level)
                                                          ! dim: (nproma, nblks_c)

    ! local
    REAL(wp) :: theta_v_ic(nproma,p_patch%nlev+1)    ! interpolated face value
    REAL(wp) :: theta_v_pr_ic(nproma,p_patch%nlev+1) ! interpolated face value (perturbation)

    !
    REAL(wp) :: exner_new  (nproma,p_patch%nlev)     ! iterative solution for exner
    REAL(wp) :: theta_v_new(nproma,p_patch%nlev)     ! iterative solution for theta_v
    REAL(wp) :: rho_new    (nproma,p_patch%nlev)     ! iterative solution for rho
    REAL(wp) :: qv_new     (nproma,p_patch%nlev)     ! iterative solution for qv
    REAL(wp) :: z_theta_v_pr(nproma,p_patch%nlev)    ! perturbation virtual potential temperature

    REAL(wp) :: tempv                                ! virtual temperature
    REAL(wp) :: pres                                 ! pressure
    REAL(wp) :: es                                   ! saturation vapour pressure

    REAL(wp), POINTER :: ptr_exner_fg(:,:,:)         ! exner first guess
                                                     ! will be either exner, or p_nh_metrics%exner_ref_mc 
                                                     ! depending on luse_exner_fg

    REAL(wp), ALLOCATABLE, TARGET :: &               ! wp precision for exner_ref_mc
      &  exner_ref_mc_wp(:,:,:)                                     

    LOGICAL  :: lintegrate_topdown                   ! TRUE : use opt_exner_ubc
                                                     ! FALSE: opt_exner_lbc

    INTEGER :: jb, jk, jc, iter                      ! loop indices
    INTEGER :: nlev
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: ist

    INTEGER :: irow                                   ! row to fill

    ! convergence monitoring
    REAL(wp), ALLOCATABLE :: diff_max_theta_v(:,:)
    REAL(wp), ALLOCATABLE :: diff_max_exner(:,:)
    REAL(wp), ALLOCATABLE :: diff_max_theta_v_tot(:)
    REAL(wp), ALLOCATABLE :: diff_max_exner_tot(:)
    !
    ! table for monitoring convergence
    TYPE(t_table)   :: table
    CHARACTER(LEN = *), PARAMETER :: iterationCntCol  = "Iteration cnt",    &
      &                              maxDiffThetavCol = "max diff THETA_V [K]", &
      &                              maxDiffExnerCol  = "max diff EXNER"

    ! table for depicting final vertical profiles 
    TYPE(t_table)   :: table_profiles
    CHARACTER(LEN = *), PARAMETER :: levelCol  = "full level height [m]", &
      &                              thetavCol = "theta_v [K]",           &
      &                              qvCol     = "qv [kg/kg]",            &
      &                              presCol   = "pres [hPa]",            &
      &                              refPresCol= "reference pres [hPa]"

    INTEGER, PARAMETER :: max_iteration = 25  ! maximum number of iterations
    !-----------------------------------------------------------------------

    ! check whether the upper or lower boundary condition for exner is provided
    IF (PRESENT(opt_exner_ubc)) THEN
      lintegrate_topdown = .TRUE.
    ELSE IF (PRESENT(opt_exner_lbc)) THEN
      lintegrate_topdown=.FALSE.
    ELSE
      CALL finish(routine, "either opt_exner_ubc or opt_exner_lbc must be given")
    ENDIF

    ! choose exner first guess
    IF (luse_exner_fg) THEN
      ptr_exner_fg => exner(:,:,:)
    ELSE
      ! this detour is necessary, as exner_ref_mc can be either single precision, 
      ! or double precision
      ALLOCATE(exner_ref_mc_wp(SIZE(p_nh_metrics%exner_ref_mc,1), &
        &                      SIZE(p_nh_metrics%exner_ref_mc,2), &
        &                      SIZE(p_nh_metrics%exner_ref_mc,3) ), STAT=ist)
      IF (ist /= SUCCESS) THEN
        CALL finish (routine, &
          &  'allocation for exner_ref_mc_wp failed '     )
      ENDIF
      exner_ref_mc_wp(:,:,:) = REAL(p_nh_metrics%exner_ref_mc(:,:,:),wp)

      ptr_exner_fg => exner_ref_mc_wp(:,:,:)
    ENDIF


    nlev = p_patch%nlev

    ALLOCATE(diff_max_theta_v(max_iteration,p_patch%nblks_c),  &
      &      diff_max_exner(max_iteration,p_patch%nblks_c),    &
      &      diff_max_theta_v_tot(max_iteration),              &
      &      diff_max_exner_tot(max_iteration), STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine, &
        &  'allocation for diff_max_theta_v, diff_max_exner, ' // &
        &  'diff_max_theta_v_tot, diff_max_exner_tot failed '     )
    ENDIF



    ! Strategy:
    ! Given the vertical profile of temperture t and relative humidity rh, 
    ! compute hydrostatically balanced profiles for exner, rho, theta_v and qv.
    !
    ! 0) Set first guess for exner and qv
    !    qv   : qv=0
    !    exner: either reference, or user-specific profile
    ! 
    ! 1) Compute t_v from t and qv
    ! 2) Compute theta_v=t_v/exner.
    ! 3) Compute exner from the discretized third equation of motion
    ! 4) Compute specific humidity qv(rh,p,es)
    ! 5) Compute rho = f(exner, theta_v) from the equation of state
    ! 6) Iterate on steps 1-5 until a desired convergence criterion is met.

    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,iter,exner_new,qv_new,theta_v_new,rho_new, &
!$OMP            tempv,z_theta_v_pr,theta_v_ic,theta_v_pr_ic,es,pres)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, rl_start, rl_end)


      ! initialize exner_new and qv_new
      !
      !
      ! set first guess for qv
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          qv_new(jc,jk)  = 0._wp
        ENDDO
      ENDDO
      !
      IF (lintegrate_topdown) THEN
        ! set exner upper boundary condition
        DO jc = i_startidx, i_endidx
          exner_new(jc,1)  = opt_exner_ubc(jc,jb)
        ENDDO
        !
        ! set first guess
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
            exner_new(jc,jk) = ptr_exner_fg(jc,jk,jb)
          ENDDO
        ENDDO

      ELSE
        ! set exner lower boundary condition
        DO jc = i_startidx, i_endidx
          exner_new(jc,nlev)  = opt_exner_lbc(jc,jb)
        ENDDO
        !
        ! set first guess
        DO jk = 1, nlev-1
          DO jc = i_startidx, i_endidx
            exner_new(jc,jk) = ptr_exner_fg(jc,jk,jb)
          ENDDO
        ENDDO
      ENDIF


      iterate: DO iter = 1,max_iteration

        !**********************************************************
        ! (1) compute virtual temperature t_v(t,qv)
        ! (2) compute theta_v(tempv,exner)
        !**********************************************************
        !
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ! virtual temperature
            tempv = temp_ini(jc,jk,jb) * (1._wp + vtmpc1*qv_new(jc,jk))
            !
            ! virtual potential temperature
            theta_v_new(jc,jk)= tempv/exner_new(jc,jk)
            !
            ! perturbation value
            z_theta_v_pr(jc,jk) = theta_v_new(jc,jk) - p_nh_metrics%theta_ref_mc(jc,jk,jb)
          ENDDO           
        ENDDO



        !***********************************************************
        ! (3) solve the discretized third equation of motion for exner, 
        !     under the assumption dw/dt=0 
        !***********************************************************
        !
        ! linear cell to face interpolation of theta_v_new and z_theta_v_pr 
        !  
        DO jk=2,nlev
          DO jc = i_startidx, i_endidx
            theta_v_ic(jc,jk) = p_nh_metrics%wgtfac_c(jc,jk,jb) * theta_v_new(jc,jk) &
              &               + (1._wp-p_nh_metrics%wgtfac_c(jc,jk,jb))*theta_v_new(jc,jk-1)

            theta_v_pr_ic(jc,jk) = p_nh_metrics%wgtfac_c(jc,jk,jb) * z_theta_v_pr(jc,jk)  &
              &               + (1._wp-p_nh_metrics%wgtfac_c(jc,jk,jb))* z_theta_v_pr(jc,jk-1)   
          ENDDO
        ENDDO


        ! Now compute hydrostatically balanced prognostic fields:
        ! The following expression is derived from the discretized third
        ! equation of motion, assuming dw/dt = 0, and solved for the exner pressure.
        IF (lintegrate_topdown) THEN
          DO jk = 2,nlev
            DO jc = i_startidx, i_endidx

              exner_new(jc,jk) = exner_new(jc,jk-1)                                                &
                &  + (p_nh_metrics%exner_ref_mc(jc,jk,jb) - p_nh_metrics%exner_ref_mc(jc,jk-1,jb)) &
                &  + p_nh_metrics%ddqz_z_half(jc,jk,jb)/theta_v_ic(jc,jk)                          &
                &  * theta_v_pr_ic(jc,jk) * p_nh_metrics%d_exner_dz_ref_ic(jc,jk,jb)

            ENDDO
          ENDDO
        ELSE
          DO jk = nlev-1,1,-1
            DO jc = i_startidx, i_endidx

              exner_new(jc,jk) = exner_new(jc,jk+1)                                                &
                &  + (p_nh_metrics%exner_ref_mc(jc,jk,jb) - p_nh_metrics%exner_ref_mc(jc,jk+1,jb)) &
                &  - p_nh_metrics%ddqz_z_half(jc,jk+1,jb)/theta_v_ic(jc,jk+1)                      &
                &  * theta_v_pr_ic(jc,jk+1) * p_nh_metrics%d_exner_dz_ref_ic(jc,jk+1,jb)

            ENDDO
          ENDDO
        ENDIF

        !**********************************************************
        ! (4) compute saturation vapour pressure es(temp)
        !     and specific humidity qv
        !     * makes use of updated exner pressure
        !**********************************************************
        !
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ! saturation vapour pressure
            es   = sat_pres_water(temp_ini(jc,jk,jb))
            pres = p0ref * EXP((cpd/rd)*LOG(exner_new(jc,jk)))
            ! specific humidity
            qv_new(jc,jk) =  rdv* rh_ini(jc,jk,jb) * es / (pres - o_m_rdv*rh_ini(jc,jk,jb)*es)
          ENDDO           
        ENDDO

        !**********************************************************
        ! (5) update rho
        !**********************************************************
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            rho_new(jc,jk) = exner_new(jc,jk)**cvd_o_rd*p0ref/(rd*theta_v_new(jc,jk))
          ENDDO
        ENDDO


        !***********************************************************
        ! check for convergence
        !***********************************************************
        diff_max_theta_v(iter,jb) = MAXVAL(ABS(theta_v_new(i_startidx:i_endidx,1:nlev)   &
          &                          - theta_v(i_startidx:i_endidx,1:nlev,jb)) )
        diff_max_exner(iter,jb)   = MAXVAL(ABS(exner_new(i_startidx:i_endidx,1:nlev)   &
          &                          - exner(i_startidx:i_endidx,1:nlev,jb)) )


        !***********************************************************
        ! copy results back to state variables
        !***********************************************************
        DO jk = 1,nlev
          DO jc = i_startidx, i_endidx
            exner(jc,jk,jb)   = exner_new(jc,jk)
            theta_v(jc,jk,jb) = theta_v_new(jc,jk)
            rho(jc,jk,jb)     = rho_new(jc,jk)
            qv(jc,jk,jb)      = qv_new(jc,jk)
          ENDDO
        ENDDO

      END DO iterate

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


    ! Take maximum over all PEs
    !
    DO iter = 1,max_iteration
      diff_max_theta_v_tot(iter) = MAXVAL(diff_max_theta_v(iter,i_startblk:i_endblk))
      diff_max_exner_tot(iter)   = MAXVAL(diff_max_exner(iter,i_startblk:i_endblk))

      diff_max_theta_v_tot(iter) = global_max(diff_max_theta_v_tot(iter))
      diff_max_exner_tot(iter)   = global_max(diff_max_exner_tot(iter))
    ENDDO

    ! convergence control output 
    !
    IF ((p_pe_work_only == 0) .AND. .NOT. my_process_is_stdio() .AND. msg_level >= 10) THEN

      ! print table for monitoring convergence behaviour
      !
      write(0,*) ""
      write(0,*) "hydrostatic adjustment: convergence control"
      CALL initialize_table(table)
      !
      CALL add_table_column(table, iterationCntCol)
      CALL add_table_column(table, maxDiffThetavCol)
      CALL add_table_column(table, maxDiffExnerCol)

      irow = 0
      DO iter=1, max_iteration
        irow = irow + 1
        CALL set_table_entry(table,irow,iterationCntCol, ADJUSTL(TRIM(int2string(iter))))
        CALL set_table_entry(table,irow,maxDiffThetavCol, ADJUSTL(TRIM(real2string(diff_max_theta_v_tot(iter)))))
        CALL set_table_entry(table,irow,maxDiffExnerCol, ADJUSTL(TRIM(real2string(diff_max_exner_tot(iter)))))
      ENDDO
      !
      CALL print_table(table, opt_delimiter=' | ')
      CALL finalize_table(table)


      ! print table with final hydrostatically balanced profiles
      !
      write(0,*) ""
      write(0,*) "Hydrostatically balanced profiles"
      CALL initialize_table(table_profiles)
      !
      CALL add_table_column(table, levelCol)
      CALL add_table_column(table, thetavCol)
      CALL add_table_column(table, qvCol)
      CALL add_table_column(table, presCol)
      CALL add_table_column(table, refPresCol)
      !
      DO jk=1, nlev
        CALL set_table_entry(table_profiles,jk,levelCol  , ADJUSTL(TRIM(real2string(p_nh_metrics%z_mc(1,jk,i_startblk)))))
        CALL set_table_entry(table_profiles,jk,thetavCol , ADJUSTL(TRIM(real2string(theta_v(1,jk,i_startblk)))))
        CALL set_table_entry(table_profiles,jk,qvCol     , ADJUSTL(TRIM(real2string(qv(1,jk,i_startblk)))))
        pres = 0.01_wp*p0ref * EXP((cpd/rd)*LOG(exner(1,jk,i_startblk)))
        CALL set_table_entry(table_profiles,jk,presCol   , ADJUSTL(TRIM(real2string(pres))))
        pres = 0.01_wp*p0ref * EXP((cpd/rd)*LOG(p_nh_metrics%exner_ref_mc(1,jk,i_startblk)))
        CALL set_table_entry(table_profiles,jk,refPresCol, ADJUSTL(TRIM(real2string(pres))))
      ENDDO
      !
      CALL print_table(table_profiles, opt_delimiter=' | ')
      CALL finalize_table(table_profiles)

    ENDIF  ! (p_pe_work_only == 0) .AND. NOT. my_process_is_stdio()

    ! cleanup
    DEALLOCATE(diff_max_theta_v, diff_max_exner, diff_max_theta_v_tot, diff_max_exner_tot)
    IF (ALLOCATED(exner_ref_mc_wp)) THEN
      NULLIFY(ptr_exner_fg)
      DEALLOCATE(exner_ref_mc_wp)
    ENDIF

  END SUBROUTINE hydro_adjust_iterative

END MODULE mo_hydro_adjust
