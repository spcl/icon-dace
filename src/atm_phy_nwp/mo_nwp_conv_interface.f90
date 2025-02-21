!OPTION! -cont
! this command should fix the problem of copying arrays in a subroutine call
!
! This module is the interface between nwp_nh_interface to the
! convection parameterisation(s):
! inwp_conv == 1 == Tiedtke-Bechtold convection
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

MODULE mo_nwp_conv_interface

  USE mo_kind,                 ONLY: wp
  USE mo_parallel_config,      ONLY: nproma
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag,&
    &                                t_nh_metrics
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend, &
       &                             t_nwp_phy_stochconv, t_ptr_cloud_ensemble
  USE mo_nwp_phy_state,        ONLY: phy_params
  USE mo_run_config,           ONLY: iqv, iqc, iqi, iqr, iqs, nqtendphy, lart
  USE mo_physical_constants,   ONLY: grav, alf, cvd, cpd
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_cumaster,             ONLY: cumastrn
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_comin_config,         ONLY: comin_config, t_comin_tracer_info
  USE mo_art_config,           ONLY: art_config
  USE mo_util_phys,            ONLY: nwp_con_gust
  USE mo_exception,            ONLY: finish, message_text

  ! for stochastic convection
  USE mo_sync,                 ONLY: sync_patch_array,sync_patch_array_mult,SYNC_C
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_math_constants,       ONLY: rad2deg
  USE mtime,                   ONLY: datetime
  USE mo_gribout_config,       ONLY: gribout_config
  USE mo_fortran_tools,        ONLY: init, t_ptr_tracer, set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE



  PUBLIC  ::  nwp_convection

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_conv_interface'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_convection ( tcall_conv_jg,             & !>input
    &                         linit,                     & !>input
    &                         p_patch,p_metrics,         & !>input
    &                         ext_data,                  & !>input
    &                         p_prog,                    & !>input
    &                         p_prog_rcf,                & !>inout
    &                         mtime_datetime,            & !>input
    &                         p_diag ,                   & !>inout
    &                         prm_diag,                  & !>inout
    &                         prm_nwp_tend,              & !>inout
    &                         prm_nwp_stochconv,         & !>inout
    &                         pt_int_state,              & !
    &                         lacc                       ) !>in

    TYPE(t_patch)               ,INTENT(in)   :: p_patch          !<grid/patch info.
    TYPE(t_external_data)       ,INTENT(in)   :: ext_data         !< external data
    TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
    TYPE(t_nh_prog)             ,INTENT(in)   :: p_prog           !<the dyn prog vars
    TYPE(datetime), POINTER     ,INTENT(in)   :: mtime_datetime
    TYPE(t_nh_prog), TARGET     ,INTENT(inout):: p_prog_rcf       !<call freq
    TYPE(t_nh_diag), TARGET     ,INTENT(inout):: p_diag           !<the dyn diag vars
    TYPE(t_nwp_phy_diag)        ,INTENT(inout):: prm_diag         !<the atm phys vars
    TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend     !< atm tend vars
    TYPE(t_nwp_phy_stochconv)   ,INTENT(inout):: prm_nwp_stochconv!< stoch conv vars
    TYPE(t_int_state)           ,INTENT(IN)   :: pt_int_state

    REAL(wp)                    ,INTENT(in)   :: tcall_conv_jg    !< time interval for
                                                                  !< convection
    LOGICAL                     ,INTENT(in)   :: linit            !< .TRUE. if initialization call
    LOGICAL, OPTIONAL           ,INTENT(in)   :: lacc             !< to run on GPUs

    CHARACTER(*), PARAMETER    :: routine = modname//"::nwp_convection"
    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full and half levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !< blocks
    INTEGER :: i_startidx, i_endidx    !< slices

    REAL(wp)         :: z_omega_p(nproma,p_patch%nlev) !< vertical velocity in p-system
    REAL(wp)         :: z_plitot (nproma,p_patch%nlev) !< cloud water + cloud ice
    REAL(wp), TARGET :: z_qhfl (nproma,p_patch%nlevp1) !< 3D moisture flux( convection)
    REAL(wp), TARGET :: z_shfl (nproma,p_patch%nlevp1) !< 3D sensible heat flux "-"
    REAL(wp)         :: z_dtdqv  (nproma,p_patch%nlev) !< 3D moisture convergence
                                                       !< on output, the convection scheme adds the convective
                                                       !< qv-tendency. It is, however no longer used.
                                                       !< Instead we make use of the \rho*qv-tendency
                                                       !< ptenrhoq
    REAL(wp) :: z_dtdt   (nproma,p_patch%nlev)         !< temporal temperature tendency
    REAL(wp) :: z_dtdt_sv(nproma,p_patch%nlev)         !< save array for temperature tendency
    REAL(wp) :: z_ddspeed(nproma)                      !< maximum downdraft speed at the surface

    ! stochastic cloud ensemble variables:
    ! variables to contain profiles/sfc fluxes averaged over halo, only used for stoch conv
    ! could be made allocatable, but then would have to allocate within block loop
    REAL(wp), TARGET, DIMENSION(nproma,p_patch%nlev)    :: qvmean, tempmean, presmean, umean, vmean
    REAL(wp), TARGET, DIMENSION(nproma,p_patch%nlevp1)  :: shfl_avg,qhfl_avg
    ! pointers to point either at averaged profiles/sfc fluxes (stoch) or "normal" profiles (default)
    REAL(wp), POINTER   :: p_qv(:,:), p_temp(:,:), p_pres(:,:), p_u(:,:), p_v(:,:)
    REAL(wp), POINTER   :: p_shfl_avg(:,:),p_qhfl_avg(:,:)

    INTEGER       :: iseed(nproma) ! seed for random number generator in stoch schemes
    ! logicals switching various stochastic convection options
    LOGICAL       :: lstoch_expl,lvvcouple,lvv_shallow_deep,lstoch_sde,lstoch_deep,lspinup
    TYPE(t_ptr_cloud_ensemble) :: p_cloud_ensemble

    ! Local scalars:
    INTEGER  :: jk,jc,jb,jg,jt,l,jc2,jb2   !< block indices
    INTEGER  :: zk850, zk950               !< level indices
    REAL(wp) :: u850, u950, v850, v950     !< zonal and meridional velocity at specific heights
    REAL(wp) :: ticeini, lfocvd, wfac, cpdocvd, area_norm
    INTEGER  :: iqrd, iqsd
    LOGICAL  :: lcompute_lpi               !< compute lpi_con, mlpi_con, koi, lpi_con_max and mlpi_con_max
    LOGICAL  :: lcompute_lfd               !< compute lfd_con, lfd_con_max
    LOGICAL  :: lzacc                      !< to check, if lacc is present

    ! Tracer specific variables:
    INTEGER :: nconv_tracer_tot, nconv
    ! General pointers to tracers/tendencies for convection
    TYPE(t_ptr_tracer), POINTER :: ptr_conv_tracer_tend(:)
    TYPE(t_ptr_tracer), POINTER :: ptr_conv_tracer(:)
    ! Pointers to structures in memory allocated within ART
    TYPE(t_ptr_tracer), POINTER :: ptr_conv_tracer_tend_art(:)
    TYPE(t_ptr_tracer), POINTER :: ptr_conv_tracer_art(:)
    ! Extended list of tracers by ComIn (incl. ART if present)
    TYPE(t_ptr_tracer), ALLOCATABLE, TARGET :: ptr_conv_tracer_tend_comin(:)
    TYPE(t_ptr_tracer), ALLOCATABLE, TARGET :: ptr_conv_tracer_comin(:)
    TYPE(t_comin_tracer_info), POINTER :: this_info => NULL()

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA &
    !$ACC   PRESENT(kstart_moist(jg)) &
    !$ACC   CREATE(z_ddspeed, z_dtdqv, z_dtdt, z_dtdt_sv) &
    !$ACC   CREATE(z_omega_p, z_plitot, z_qhfl, z_shfl) &
    !$ACC   IF(lzacc)

    ! local variables related to the blocking
    jg        = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    lfocvd  = alf/cvd
    cpdocvd = cpd/cvd
    ticeini = 256.15_wp

    ! IDs for optional arguments for detrainment of rain and snow
    IF (atm_phy_nwp_config(jg)%ldetrain_conv_prec) THEN
      iqrd = iqr
      iqsd = iqs
    ELSE
      iqrd = nqtendphy
      iqsd = nqtendphy
    ENDIF

    nconv_tracer_tot = comin_config%comin_icon_domain_config(jg)%nconv_tracer
    IF (lart) nconv_tracer_tot = nconv_tracer_tot + art_config(jg)%nconv_tracer
    IF (comin_config%comin_icon_domain_config(jg)%nconv_tracer > 0) THEN
      ALLOCATE( ptr_conv_tracer_tend_comin(nconv_tracer_tot), ptr_conv_tracer_comin(nconv_tracer_tot) )
    ENDIF

    ! switch on stochastic convection scheme with namelist parameter, default is false/off
    lstoch_expl=atm_phy_nwp_config(jg)%lstoch_expl
    ! switch on stochastic differential equations. default is off. If lstoch_expl=T lstoch_sde must be false.
    lstoch_sde=atm_phy_nwp_config(jg)%lstoch_sde
    ! switch on stochastic scheme for deep convection. default is off.
    lstoch_deep=atm_phy_nwp_config(jg)%lstoch_deep
    ! decide whether to use vertical velocity at 650hPa to couple shallow Cu to resolved deep convection
    !   default: false/off
    lvvcouple=atm_phy_nwp_config(jg)%lvvcouple
    ! decide whether to use vertical velocity at 650hP to distinguish between shallow/deep convective
    !   grid points within the convection code
    lvv_shallow_deep=atm_phy_nwp_config(jg)%lvv_shallow_deep
    ! spinup cloud ensemble only during inital slow physics call, and when spinup
    ! is requested via namelist parameter

    IF (atm_phy_nwp_config(jg)%lstoch_spinup .AND. linit) THEN
       lspinup=.TRUE.
    ELSE
       lspinup=.FALSE.
    ENDIF

    IF (lstoch_expl .or. lstoch_sde .or. lstoch_deep) THEN
       ! initialize stochstic diagnostic variables:
       CALL init(prm_diag%mf_b)
       CALL init(prm_diag%mf_p)
       CALL init(prm_diag%mf_num)
    ENDIF

    IF (lstoch_expl .or. lstoch_sde) THEN
      CALL sync_patch_array_mult(SYNC_C,p_patch,5,p_diag%temp,p_prog_rcf%tracer(:,:,:,iqv),p_diag%pres,p_diag%u,p_diag%v)
      CALL sync_patch_array(SYNC_C,p_patch,prm_diag%shfl_s)
      CALL sync_patch_array(SYNC_C,p_patch,prm_diag%qhfl_s)
    ENDIF

    ! compute lpi_con(_max) only if all relevant fields are allocated (non-dummy).
    ! This is only the case, if any of the lpi-fields is requested in the output_nml.
    ! GZ: taking the size product of all 4 fields causes an integer overflow for nproma > 215 (2**7.75)
    lcompute_lpi = SIZE(prm_diag%lpi_con_max,1) * SIZE(prm_diag%lpi_con,1) > 0 .AND.  &
      &            SIZE(prm_diag%mlpi_con_max,1)* SIZE(prm_diag%mlpi_con,1) > 0
    !
    ! compute lfd_con(_max) only if all relevant fields are allocated (non-dummy).
    lcompute_lfd = SIZE(prm_diag%lfd_con_max,1) * SIZE(prm_diag%lfd_con,1) > 0

!FIXME: PGI + OpenMP produce deadlock in this loop. Compiler bug suspected
!$OMP PARALLEL DO PRIVATE(jb,jc,jk,jt,i_startidx,i_endidx,z_omega_p,z_plitot,z_qhfl,z_shfl,z_dtdqv,&
!$OMP            z_dtdt,z_dtdt_sv,zk850,zk950,u850,u950,v850,v950,wfac,z_ddspeed,&
!$OMP            iseed,presmean,umean,vmean,qvmean,tempmean,qhfl_avg,shfl_avg,l,jc2,jb2,area_norm, &
!$OMP            p_pres,p_u,p_v,p_qv,p_temp,p_qhfl_avg,p_shfl_avg,p_cloud_ensemble), ICON_OMP_GUIDED_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)


      IF( atm_phy_nwp_config(jg)%inwp_convection == 1 ) THEN

        !>
        !! define convective-related fields
        !! NOTE: Heat fluxes are defined negative upwards in the convection scheme.
        !!       In the turbulences scheme (1,2,3) they defined either positive or
        !!       negative upward.
        !! Thus pass fluxes to cumastrn that are negative when upwards!!!

        SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)
        CASE (0)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx,i_endidx
            z_qhfl(jc,nlevp1) = - 4.79846_wp*1.e-5_wp !> moisture flux kg/m2/s
            z_shfl(jc,nlevp1) = - 17._wp              !! sens. heat fl W/m**2
          ENDDO
          !$ACC END PARALLEL

        CASE DEFAULT

          ! In turb1,turb2 and turb3, the flux is positive downwards / negative upwards

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx,i_endidx
            z_qhfl(jc,nlevp1) = prm_diag%qhfl_s(jc,jb)
            z_shfl(jc,nlevp1) = prm_diag%shfl_s(jc,jb)
          ENDDO
          !$ACC END PARALLEL

        END SELECT

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG
        DO jk = 1, kstart_moist(jg)-1
          !$ACC LOOP VECTOR
          DO jc = 1, nproma
            z_omega_p (jc,jk) = 0._wp
            z_dtdqv   (jc,jk) = 0._wp
            z_dtdt    (jc,jk) = 0._wp
            z_plitot  (jc,jk) = 0._wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL


        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG
        DO jk = kstart_moist(jg),nlev
          !$ACC LOOP VECTOR
          DO jc = i_startidx,i_endidx
            ! vertical velocity in p-system
            z_omega_p(jc,jk)= -0.5_wp*(p_prog%w(jc,jk,jb)+p_prog%w(jc,jk+1,jb)) &
                            * p_prog%rho(jc,jk,jb)*grav
            ! moisture convergence
            z_dtdqv(jc,jk) =                                                                 &
              p_diag%ddt_tracer_adv(jc,jk,jb,iqv) + prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqv)

            ! temperature tendencies from other physical processes (used for mass flux closure)
            z_dtdt(jc,jk) =                              &
              &    prm_nwp_tend%ddt_temp_radsw(jc,jk,jb) &
              &  + prm_nwp_tend%ddt_temp_radlw(jc,jk,jb) &
              &  + prm_nwp_tend%ddt_temp_turb (jc,jk,jb) &
              &  + p_diag%ddt_temp_dyn(jc,jk,jb)
            z_dtdt_sv(jc,jk) = z_dtdt(jc,jk)

            ! cloud water + cloud ice for entrainment computation
            z_plitot(jc,jk) = p_prog_rcf%tracer(jc,jk,jb,iqc) &
                            + p_prog_rcf%tracer(jc,jk,jb,iqi)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ! The following input fields must be reset to zero because the convective
        ! tendencies are added to them
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG
        DO jk = 1, nlev
          !$ACC LOOP VECTOR
          DO jc = 1, nproma
            prm_nwp_tend%ddt_u_pconv (jc,jk,jb) = 0._wp
            prm_nwp_tend%ddt_v_pconv (jc,jk,jb) = 0._wp

            prm_diag%rain_con_rate_3d(jc,jk,jb) = 0._wp
            prm_diag%snow_con_rate_3d(jc,jk,jb) = 0._wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        IF ( lart .AND. art_config(jg)%nconv_tracer > 0 ) THEN
          DO jt=1,art_config(jg)%nconv_tracer
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = 1, nproma
                prm_nwp_tend%conv_tracer_tend(jb,jt)%ptr(jc,jk) = 0._wp
              END DO
            END DO
            !$ACC END PARALLEL
          ENDDO
          ptr_conv_tracer_tend_art => prm_nwp_tend%conv_tracer_tend(jb,1:art_config(jg)%nconv_tracer)
          ptr_conv_tracer_art      => p_prog_rcf%conv_tracer(jb,1:art_config(jg)%nconv_tracer) 
        ELSE
          ptr_conv_tracer_art      => NULL()
          ptr_conv_tracer_tend_art => NULL()
        ENDIF

        !-------------------------------------------------------------------------
        !> Convection
        !-------------------------------------------------------------------------

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP VECTOR
        DO jc = i_startidx,i_endidx                         ! ldshcv is set in mo_nwp_turb_sfc_interface.f90
          prm_diag%ldshcv(jc,jb) = .TRUE.                   ! here: option to overwrite DUALM choice
        ENDDO
        !$ACC END PARALLEL

        ! Preparing fields for stochastic convection routines
        IF (lstoch_expl .or. lstoch_sde .or. lstoch_deep) THEN

           DO jc = i_startidx, i_endidx
              !  Generate seed for random number generator
              iseed(jc) = (INT(mtime_datetime%time%hour) * 60) + INT(mtime_datetime%time%minute)+ &
        &                  INT(mtime_datetime%date%day) +  NINT(p_patch%cells%center(jc,jb)%lon*108000000._wp + &
        &                  (p_patch%cells%center(jc,jb)%lat * rad2deg * 100._wp) )+ &
        &                  gribout_config(1)%perturbationNumber*10000
              ! Initalise  fields to contain averaged surface fluxes
              shfl_avg(jc,nlevp1)=0.0_wp
              qhfl_avg(jc,nlevp1)=0.0_wp
           END DO
        ENDIF

        ! Average input profiles and surface fluxes over halo for stochastic convection schemes.
        ! NOTE! Currently only applied to shallow convection schemes (explicit or SDE), not to
        ! deep stochastic scheme (averaging impacts overall global convective activity too much).
        IF (lstoch_expl .or. lstoch_sde) THEN !don't average for deep

           ! average surface fluxes
           DO l=1, pt_int_state%cell_environ%max_nmbr_nghbr_cells
              DO jc = i_startidx, i_endidx
                 jc2       = pt_int_state%cell_environ%idx      ( jc, jb, l)
                 jb2       = pt_int_state%cell_environ%blk      ( jc, jb, l)
                 area_norm = pt_int_state%cell_environ%area_norm( jc, jb, l)

                 shfl_avg(jc,nlevp1) = shfl_avg(jc,nlevp1) + area_norm * prm_diag%shfl_s(jc2,jb2)
                 qhfl_avg(jc,nlevp1) = qhfl_avg(jc,nlevp1) + area_norm * prm_diag%qhfl_s(jc2,jb2)
              END DO! i_startidx
           END DO! l cell neighbours
           ! Exclude cell row adjacent to the boundary interpolation zone from averaging
           DO jc = i_startidx, i_endidx
             IF (p_patch%cells%refin_ctrl(jc,jb) == grf_bdywidth_c+1) THEN
               shfl_avg(jc,nlevp1) = prm_diag%shfl_s(jc,jb)
               qhfl_avg(jc,nlevp1) = prm_diag%qhfl_s(jc,jb)
             ENDIF
           ENDDO

           ! average profiles for T, qv, u, v, pressure
           DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                 presmean( jc, jk) = 0.0_wp
                  tempmean( jc, jk) = 0.0_wp
                  umean( jc, jk) = 0.0_wp
                  vmean( jc, jk) = 0.0_wp
                  qvmean( jc, jk) = 0.0_wp
               enddo
               DO l=1, pt_int_state%cell_environ%max_nmbr_nghbr_cells
                  DO jc = i_startidx, i_endidx
                     jc2       = pt_int_state%cell_environ%idx      ( jc, jb, l)
                     jb2       = pt_int_state%cell_environ%blk      ( jc, jb, l)
                     area_norm = pt_int_state%cell_environ%area_norm( jc, jb, l)

                     presmean(jc,jk) = presmean(jc,jk) + area_norm * p_diag%pres(jc2,jk,jb2)
                     tempmean(jc,jk) = tempmean(jc,jk) + area_norm * p_diag%temp(jc2,jk,jb2)
                     qvmean(jc,jk)   = qvmean(jc,jk)   + area_norm * p_prog_rcf%tracer(jc2,jk,jb2,iqv)
                     umean(jc,jk)    = umean(jc,jk)    + area_norm * p_diag%u(jc2,jk,jb2)
                     vmean(jc,jk)    = vmean(jc,jk)    + area_norm * p_diag%v(jc2,jk,jb2)
                  END DO! i_startidx
              END DO! l cell neighbours
           END DO! nlev

           ! Pointer assignment
           p_pres     => presmean
           p_temp     => tempmean
           p_qv       => qvmean
           p_u        => umean
           p_v        => vmean
           p_shfl_avg => shfl_avg
           p_qhfl_avg => qhfl_avg
        ELSE !lstoch_
          ! without stochastic schemes, point at "normal" profiles, surface fluxes
          p_pres     => p_diag%pres(:,:,jb)
          p_temp     => p_diag%temp(:,:,jb)
          p_u        => p_diag%u(:,:,jb)
          p_v        => p_diag%v(:,:,jb)
          p_qv       => p_prog_rcf%tracer(:,:,jb,iqv)
          p_shfl_avg => z_shfl(:,:)
          p_qhfl_avg => z_qhfl(:,:)
        ENDIF

        IF (lstoch_expl) THEN
           p_cloud_ensemble%mf_i       => prm_nwp_stochconv%mf_i(:,:,jb)
           p_cloud_ensemble%time_i     => prm_nwp_stochconv%time_i(:,:,jb)
           p_cloud_ensemble%life_i     => prm_nwp_stochconv%life_i(:,:,jb)
           p_cloud_ensemble%area_i     => prm_nwp_stochconv%area_i(:,:,jb)
           p_cloud_ensemble%type_i     => prm_nwp_stochconv%type_i(:,:,jb)
           p_cloud_ensemble%ktype_i    => prm_nwp_stochconv%ktype_i(:,:,jb)
           p_cloud_ensemble%depth_i    => prm_nwp_stochconv%depth_i(:,:,jb)
           p_cloud_ensemble%base_i     => prm_nwp_stochconv%base_i(:,:,jb)
           p_cloud_ensemble%used_cell  => prm_nwp_stochconv%used_cell(:,:,jb)
        ELSE
           NULLIFY(p_cloud_ensemble%mf_i)
           NULLIFY(p_cloud_ensemble%time_i)
           NULLIFY(p_cloud_ensemble%life_i)
           NULLIFY(p_cloud_ensemble%area_i)
           NULLIFY(p_cloud_ensemble%type_i)
           NULLIFY(p_cloud_ensemble%ktype_i)
           NULLIFY(p_cloud_ensemble%depth_i)
           NULLIFY(p_cloud_ensemble%base_i)
           NULLIFY(p_cloud_ensemble%used_cell)          
        ENDIF

        IF ( comin_config%comin_icon_domain_config(jg)%nconv_tracer > 0 ) THEN
          ! Check if ART tracers are present, fill the first part of ptr_conv_tracer_comin / ptr_conv_tracer_tend_comin with ART lists
          IF ( lart .AND. art_config(jg)%nconv_tracer > 0 ) THEN
            ptr_conv_tracer_comin     (1:art_config(jg)%nconv_tracer) = p_prog_rcf%conv_tracer       (jb,1:art_config(jg)%nconv_tracer)
            ptr_conv_tracer_tend_comin(1:art_config(jg)%nconv_tracer) = prm_nwp_tend%conv_tracer_tend(jb,1:art_config(jg)%nconv_tracer)
            nconv = art_config(jg)%nconv_tracer
          ELSE
            nconv = 0
          ENDIF

          ! Append ComIn tracers
          this_info => comin_config%comin_icon_domain_config(jg)%tracer_info_head
          DO WHILE (ASSOCIATED(this_info))
            IF (this_info%idx_conv > 0) THEN
              nconv = nconv + 1
              ptr_conv_tracer_comin     (nconv)%ptr => p_prog_rcf%tracer(:,:,jb,this_info%idx_tracer)
              ptr_conv_tracer_tend_comin(nconv)%ptr => prm_nwp_tend%ddt_tracer_pconv(:,:,jb,this_info%idx_conv)
            ENDIF !this_info%idx_conv > 0
            this_info => this_info%next
          ENDDO
          IF (nconv /= nconv_tracer_tot) THEN
            WRITE (message_text,'(A,A,I3,A,I3)') "Number of convective tracer is inconsistent: ",  &
              &                               "nconv = ",nconv,", nconv_tracer_tot = ",nconv_tracer_tot
            CALL finish(routine, message_text)
          ENDIF !nconv /= nconv_tracer_tot

          ! Overwrite pointers
          ptr_conv_tracer      => ptr_conv_tracer_comin(1:nconv_tracer_tot)
          ptr_conv_tracer_tend => ptr_conv_tracer_tend_comin(1:nconv_tracer_tot)

        ELSE ! points to NULL() if no ART tracers are present
          ptr_conv_tracer      => ptr_conv_tracer_art
          ptr_conv_tracer_tend => ptr_conv_tracer_tend_art
        ENDIF !comin_config%comin_icon_domain_config(jg)%nconv_tracer > 0

        CALL cumastrn &
&         (kidia  = i_startidx            , kfdia  = i_endidx               ,& !> IN
&          klon   = nproma ,     ktdia  = kstart_moist(jg)  , klev = nlev   ,& !! IN
&          ldland = ext_data%atm%llsm_atm_c(:,jb), ptsphy = tcall_conv_jg   ,& !! IN
&          ldlake = ext_data%atm%llake_c(:,jb), k950 = prm_diag%k950(:,jb)  ,& !! IN
&          phy_params = phy_params(jg),trop_mask=prm_diag%tropics_mask(:,jb),& !! IN
&          mtnmask=p_metrics%mask_mtnpoints(:,jb), pten=p_temp(:,:)         ,& !! IN
&          pqen   = p_qv(:,:)                                               ,& !! IN
&          puen   = p_u(:,:)              , pven   = p_v(:,:)               ,& !! IN
&          plitot = z_plitot              , pvervel= z_omega_p              ,& !! IN
&          plen  =p_prog_rcf%tracer(:,:,jb,iqc), pien=p_prog_rcf%tracer(:,:,jb,iqi),& !!IN
&          shfl_s = prm_diag%shfl_s(:,jb) , qhfl_s=prm_diag%qhfl_s(:,jb)    ,& !! IN
&          pqhfl  = p_qhfl_avg            , pahfs  = p_shfl_avg             ,& !! IN
&          pap    = p_pres(:,:)       , paph   = p_diag%pres_ifc(:,:,jb)    ,& !! IN
&          pgeo   = p_metrics%geopot_agl (:,:,jb)                           ,& !! IN
&          pgeoh  = p_metrics%geopot_agl_ifc(:,:,jb)                        ,& !! IN
&          zdph   = p_diag%dpres_mc     (:,:,jb)                            ,& !! IN
&          zdgeoh = p_metrics%dgeopot_mc(:,:,jb)                            ,& !! IN
&          pcloudnum = prm_diag%cloud_num(:,jb)                             ,& !! IN
&          ptenta = p_diag%ddt_temp_dyn(:,:,jb)                             ,& !! IN
&          ptenqa = p_diag%ddt_tracer_adv(:,:,jb,iqv)                       ,& !! IN
&          ptent  = z_dtdt                                                  ,& !! INOUT
&          ptenu  = prm_nwp_tend%ddt_u_pconv     (:,:,jb)                   ,& !! OUT
&          ptenv  = prm_nwp_tend%ddt_v_pconv     (:,:,jb)                   ,& !! OUT
&          ptenq  = z_dtdqv                                                 ,& !! INOUT
&          ptenrhoq  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqv)            ,& !! OUT
&          ptenrhol  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqc)            ,& !! OUT
&          ptenrhoi  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqi)            ,& !! OUT
&          ptenrhor  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqrd)           ,& !! OUT
&          ptenrhos  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqsd)           ,& !! OUT
&          ldcum  = prm_diag%locum   (:,jb)                                 ,& !! OUT
&          ktype  = prm_diag%ktype   (:,jb)                                 ,& !! OUT
&          kcbot  = prm_diag%mbas_con(:,jb)                                 ,& !! OUT
&          kctop  = prm_diag%mtop_con(:,jb)                                 ,& !! OUT
&          LDSHCV = prm_diag%ldshcv  (:,jb)                                 ,& !! IN
&          pmfu   =      prm_diag%con_udd(:,:,jb,1)                         ,& !! OUT
&          pmfd   =      prm_diag%con_udd(:,:,jb,2)                         ,& !! OUT
&          pmfude_rate = prm_diag%con_udd(:,:,jb,3)                         ,& !! OUT
&          pmfdde_rate = prm_diag%con_udd(:,:,jb,4)                         ,& !! OUT
&          ptu    =      prm_diag%con_udd(:,:,jb,5)                         ,& !! OUT
&          pqu    =      prm_diag%con_udd(:,:,jb,6)                         ,& !! OUT
&          plu    =      prm_diag%con_udd(:,:,jb,7)                         ,& !! OUT
&          pcore  =      prm_diag%con_udd(:,:,jb,8)                         ,& !! OUT
&          pmflxr =      prm_diag%rain_con_rate_3d(:,:,jb)                  ,& !! OUT
&          pmflxs =      prm_diag%snow_con_rate_3d(:,:,jb)                  ,& !! OUT
&          prain  =      prm_diag%rain_upd (:,jb)                           ,& !! OUT
&          pdtke_con =   prm_nwp_tend%ddt_tke_pconv(:,:,jb)                 ,& !! OUT
&          pcape  =      prm_diag%cape     (:,jb)                           ,& !! OUT
&          pvddraf =     z_ddspeed(:)                                       ,& !! OUT
&          pcen     =    ptr_conv_tracer                                    ,& !! IN
&          ptenrhoc =    ptr_conv_tracer_tend                               ,& !! OUT
&          l_lpi  =      lcompute_lpi                                       ,& !! IN
&          l_lfd  =      lcompute_lfd                                       ,& !! IN
&          lpi    =      prm_diag%lpi_con(:,jb)                             ,& !! OUT
&          mlpi   =      prm_diag%mlpi_con(:,jb)                            ,& !! OUT
&          koi    =      prm_diag%koi(:,jb)                                 ,& !! OUT
&          lfd    =      prm_diag%lfd_con(:,jb)                             ,& !! OUT
&          peis   =      prm_diag%conv_eis(:,jb)                            ,& !! OUT
&          lspinup      = lspinup                                           ,& !! IN
&          k650=         prm_diag%k650(:,jb)                                ,& !! IN
&          k700=         prm_diag%k700(:,jb)                                ,& !! IN
&          temp_s =      p_diag%temp(:,nlev,jb)                             ,& !! IN
&          cell_area    = p_patch%cells%area(:,jb)                          ,& !! IN
&          iseed        = iseed                                             ,& !! IN
&          mf_bulk=      prm_diag%mf_b(:,jb)                                ,& !! OUT  
&          mf_perturb =  prm_diag%mf_p(:,jb)                                ,& !! OUT
&          mf_num =      prm_diag%mf_num(:,jb)                              ,& !! OUT
&          p_cloud_ensemble = p_cloud_ensemble                              ,& !! INOUT
&          pclnum_a     = prm_nwp_stochconv%clnum_a(:,jb)                   ,& !! INOUT
&          pclmf_a      = prm_nwp_stochconv%clmf_a(:,jb)                    ,& !! INOUT
&          pclnum_p     = prm_nwp_stochconv%clnum_p(:,jb)                   ,& !! INOUT
&          pclmf_p      = prm_nwp_stochconv%clmf_p(:,jb)                    ,& !! INOUT
&          pclnum_d     = prm_nwp_stochconv%clnum_d(:,jb)                   ,& !! INOUT
&          pclmf_d      = prm_nwp_stochconv%clmf_d(:,jb)                    ,& !! INOUT
&          lacc   = lzacc                                                   )  !! IN


        ! Postprocessing on some fields

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)

        ! Conversion from temperature tendencies at constant pressure to constant volume is now done here
        !$ACC LOOP SEQ
        DO jk = kstart_moist(jg),nlev
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx,i_endidx
            prm_nwp_tend%ddt_temp_pconv  (jc,jk,jb) =  &
              &  ( z_dtdt   (jc,jk) - z_dtdt_sv(jc,jk) ) * cpdocvd
          ENDDO
        ENDDO


        ! Convert detrained cloud ice into cloud water if the temperature is only slightly below freezing
        ! and convective cloud top is not cold enough for substantial ice initiation
        !$ACC LOOP SEQ
        DO jk = kstart_moist(jg),nlev
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wfac)
          DO jc = i_startidx,i_endidx
            IF (prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi) > 0._wp .AND. p_diag%temp(jc,jk,jb) > ticeini) THEN
              wfac = MAX(0._wp, MIN(1._wp,0.25_wp*(p_diag%temp(jc,jk,jb)-ticeini)) + &
                     0.25_wp*MIN(0._wp,p_diag%temp(jc,prm_diag%mtop_con(jc,jb),jb)-ticeini) )
              prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc) = prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc) + &
                wfac*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)

              prm_nwp_tend%ddt_temp_pconv(jc,jk,jb) = prm_nwp_tend%ddt_temp_pconv(jc,jk,jb) - &
                lfocvd*wfac*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)/p_prog%rho(jc,jk,jb)

              prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi) = (1._wp-wfac)*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)
            ENDIF
          ENDDO
        ENDDO


        ! convective contribution to wind gust
        ! (based on simple parameterization by Peter Bechthold)
        !
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zk850, zk950, u850, u950, v850, v950)
        DO jc=i_startidx,i_endidx
          IF ( prm_diag%ktype(jc,jb) == 1 )  THEN   ! penetrative convection
            zk850 = prm_diag%k850(jc,jb)
            zk950 = prm_diag%k950(jc,jb)
            ! We take the arithmetic mean of u(jc,zk850,jb) and u(jc,zk850-1,jb)
            ! as well as v(jc,zk850,jb) and v(jc,zk850-1,jb), since the levels
            ! zk850 and zk950 are located just below the respective threshold heights.
            u850 = 0.5_wp * (p_diag%u(jc,zk850,jb) + p_diag%u(jc,zk850-1,jb))
            u950 = 0.5_wp * (p_diag%u(jc,zk950,jb) + p_diag%u(jc,zk950-1,jb))
            v850 = 0.5_wp * (p_diag%v(jc,zk850,jb) + p_diag%v(jc,zk850-1,jb))
            v950 = 0.5_wp * (p_diag%v(jc,zk950,jb) + p_diag%v(jc,zk950-1,jb))

            prm_diag%con_gust(jc,jb) = nwp_con_gust( u850, u950, v850, v950 ) * MIN(1._wp,0.5_wp*z_ddspeed(jc))
          ELSE
            prm_diag%con_gust(jc,jb) = 0._wp
          ENDIF
        ENDDO  ! jc


        IF (lcompute_lpi) THEN
          ! Store the maximum of lpi_con and mlpi_con
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx,i_endidx
            prm_diag%lpi_con_max(jc,jb)=MAX(prm_diag%lpi_con_max(jc,jb),      &
              &                             prm_diag%lpi_con    (jc,jb))
            prm_diag%mlpi_con_max(jc,jb)=MAX(prm_diag%mlpi_con_max(jc,jb),    &
              &                              prm_diag%mlpi_con    (jc,jb))
          ENDDO
        ENDIF

        IF (lcompute_lfd) THEN
          ! Store the maximum of lfd_con
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx,i_endidx
            prm_diag%lfd_con_max(jc,jb)=MAX(prm_diag%lfd_con_max(jc,jb),      &
              &                             prm_diag%lfd_con    (jc,jb))
          ENDDO
        ENDIF

        !$ACC END PARALLEL

      ENDIF !inwp_conv

    ENDDO  ! jb
    !$ACC WAIT(1)

!$OMP END PARALLEL DO


    IF ( ALLOCATED( ptr_conv_tracer_tend_comin) ) DEALLOCATE(ptr_conv_tracer_tend_comin)
    IF ( ALLOCATED( ptr_conv_tracer_comin) )      DEALLOCATE(ptr_conv_tracer_comin)

    !$ACC END DATA

  END SUBROUTINE nwp_convection

END MODULE mo_nwp_conv_interface

