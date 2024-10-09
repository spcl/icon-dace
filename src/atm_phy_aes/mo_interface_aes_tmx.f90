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

! Subroutine interface_aes_tmx calls the turbulent mixing scheme
! and the surface schemes using tmx.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_aes_tmx

  USE mo_kind                ,ONLY: wp, sp
  USE mtime                  ,ONLY: datetime, OPERATOR(>)

  USE mo_exception           ,ONLY: finish, warning

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: ntracer
  USE mo_aes_phy_config      ,ONLY: aes_phy_config, aes_phy_tc, dt_zero
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field, &
    &                               t_aes_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_vdf

  USE mo_ccycle_config       ,ONLY: ccycle_config
  USE mo_physical_constants  ,ONLY: amco2, amd
  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2mmr

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqnc, iqni, iqt, ico2
  USE mo_vdf                 ,ONLY: t_vdf, new_vdf
  
  USE mo_aes_sfc_indices     ,ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_surface_diag        ,ONLY: nsurf_diag
  USE mo_run_config          ,ONLY: lart
  USE mo_aes_vdf_config      ,ONLY: aes_vdf_config
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_impl_constants_grf  ,ONLY: grf_bdywidth_c
  USE mo_impl_constants      ,ONLY: min_rlcell_int, min_rlcell,max_dom
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_nh_testcases_nml    ,ONLY: nh_test_name

  USE mo_aes_thermo,          ONLY: potential_temperature
  
#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_aes_tmx, init_tmx, vdf_dom

  TYPE :: t_vdf_p
    TYPE(t_vdf), POINTER :: p
  END TYPE t_vdf_p

  TYPE(t_vdf_p) :: vdf_dom(max_dom)

  LOGICAL, SAVE :: l_init_or_restart = .TRUE.

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_aes_tmx'

CONTAINS

  SUBROUTINE interface_aes_tmx  (patch                ,&
       &                         is_in_sd_ed_interval ,&
       &                         is_active            ,&
       &                         datetime_old         ,&
       &                         dtime               )

    USE mo_variable, ONLY: unbind_variable, bind_variable
    USE mo_vdf_atmo, ONLY: t_vdf_atmo_inputs
    USE mo_vdf_sfc, ONLY: t_vdf_sfc_inputs

    ! Arguments
    !
    TYPE(t_patch)   ,TARGET ,INTENT(inout) :: patch
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: dtime

    ! Shortcuts
    !
    LOGICAL :: lparamcpl, l2moment, ldtrad_gt0
    INTEGER :: fc_vdf
    TYPE(t_aes_phy_field)   ,POINTER    :: field
    TYPE(t_aes_phy_tend)    ,POINTER    :: tend

    TYPE(t_vdf), POINTER :: vdf

    ! Local variables
    !
    INTEGER  :: jg, jc, jk, jsfc, jt, jice, copy_nblks_c
    INTEGER  :: jb,jbs,jbe,jcs,jce,ncd,rls,rle
    INTEGER  :: nlev, nlevm1, nlevp1
    INTEGER  :: ntrac
    INTEGER  :: nice        ! for simplicity (ice classes)
    !
    ! Pointers to results (nproma,patch%nlev,patch%nblks_c)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & tend_ta_vdf, tend_ua_vdf, tend_va_vdf, &
      & tend_qv_vdf, tend_qc_vdf, tend_qi_vdf
    ! Pointers to results (nproma,patch%nlev,patch%nblks_c,ntracer)
    ! REAL(wp), POINTER, DIMENSION(:,:,:,:) :: &
    !   & tend_qtrc_vdf
    ! Pointers to results (nproma,patch%nlev+1,patch%nblks_c)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & tend_wa_vdf
    ! Pointers to results (nproma,patch%nblks_c,nsfc_type)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & tend_ts

    ! Local varaibles
    REAL(wp) :: q_rlw_impl(nproma,patch%nblks_c), &
      &         tend_ta_rlw_impl(nproma,patch%nblks_c)
    REAL(wp), POINTER :: ptr_r2d(:,:), ptr_r3d(:,:,:)

    !
    INTEGER, POINTER :: turb
    LOGICAL :: l_use_rad

    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_aes_tmx'

    IF (ltimer) CALL timer_start(timer_vdf)

    jg           = patch%id
    nlev         = patch%nlev

    ! associate pointers
    ! lparamcpl =  aes_phy_config(jg)%lparamcpl
    ! l2moment  =  aes_phy_config(jg)%l2moment
    ! ldtrad_gt0=  aes_phy_tc(jg)%dt_rad > dt_zero
    ! fc_vdf    =  aes_phy_config(jg)%fc_vdf
    field     => prm_field(jg)
    tend      => prm_tend (jg)
    vdf       => vdf_dom  (jg)%p

    nlevm1 = nlev-1
    nlevp1 = nlev+1
    ntrac  = ntracer-iqt+1  ! number of tracers excluding water vapour and hydrometeors
 
    nice   = prm_field(jg)%kice
    turb => aes_vdf_config(jg)%turb

    IF ( is_in_sd_ed_interval ) THEN
      !
      IF ( is_active ) THEN

        l_use_rad = aes_phy_tc(jg)%dt_rad > dt_zero

        rls = grf_bdywidth_c + 1
        rle = min_rlcell_int

        jbs = patch%cells%start_block(rls)
        jbe = patch%cells%end_block  (rle)

        IF (.NOT. aes_vdf_config(jg)%use_tmx) THEN
          CALL finish(routine, 'ERROR: namelist parameter aes_vdf_config%use_tmx=.FALSE.!')
        END IF

        ! After init or restart, hand over some diagnostics
        IF (l_init_or_restart) THEN
          ! vdf%sfc%Get_diagnostic_r3d('roughness length momentum, tile') = field%z0m_tile(:,:,:)
          ! vdf%sfc%Get_diagnostic_r3d('roughness length heat, tile')     = field%z0h_tile(:,:,:)
          l_init_or_restart = .FALSE.
        END IF

        !
        ! Some variable pointers are swapped in mo_interface_iconam_aes.f90:interface_iconam_aes
        ! before and after the physics is called so that they are potentially pointing to targets
        ! that are different from the targets set during initialization of tmx.
        ! Here we need to make sure that the tmx variables point to the correct targets at each time step.
        !
        SELECT TYPE (ins => vdf%atmo%inputs)
        TYPE IS (t_vdf_atmo_inputs)
    
          ptr_r3d => field% qtrc_phy(:,:,:,iqv)
          CALL unbind_variable(vdf%atmo%states%search('water vapor'))
          CALL bind_variable(vdf%atmo%states%search('water vapor'), ptr_r3d)
          CALL unbind_variable(vdf%atmo%inputs%list%Search('water vapor'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('water vapor'), ptr_r3d)
          ins%pqm1 => ins%list%Get_ptr_r3d('water vapor')
          __acc_attach(ins%pqm1)
      
          ptr_r3d => field% qtrc_phy(:,:,:,iqc)
          CALL unbind_variable(vdf%atmo%states%search('cloud water'))
          CALL bind_variable(vdf%atmo%states%search('cloud water'), ptr_r3d)
          CALL unbind_variable(vdf%atmo%inputs%list%Search('cloud water'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('cloud water'), ptr_r3d)
          ins%pxlm1 => ins%list%Get_ptr_r3d('cloud water')
          __acc_attach(ins%pxlm1)

          ptr_r3d => field% qtrc_phy(:,:,:,iqi)
          CALL unbind_variable(vdf%atmo%states%search('cloud ice'))
          CALL bind_variable(vdf%atmo%states%search('cloud ice'), ptr_r3d)
          CALL unbind_variable(vdf%atmo%inputs%list%Search('cloud ice'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('cloud ice'), ptr_r3d)
          ins%pxim1 => ins%list%Get_ptr_r3d('cloud ice')
          __acc_attach(ins%pxim1)
      
          ptr_r3d => field% qtrc_phy(:,:,:,iqr)
          CALL unbind_variable(vdf%atmo%inputs%list%Search('rain'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('rain'), ptr_r3d)
          ins%pxrm1 => ins%list%Get_ptr_r3d('rain')
          __acc_attach(ins%pxrm1)
      
          ptr_r3d => field% qtrc_phy(:,:,:,iqs)
          CALL unbind_variable(vdf%atmo%inputs%list%Search('snow'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('snow'), ptr_r3d)
          ins%pxsm1 => ins%list%Get_ptr_r3d('snow')
          __acc_attach(ins%pxsm1)
      
          ptr_r3d => field% qtrc_phy(:,:,:,iqg)
          CALL unbind_variable(vdf%atmo%inputs%list%Search('graupel'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('graupel'), ptr_r3d)
          ins%pxgm1 => ins%list%Get_ptr_r3d('graupel')
          __acc_attach(ins%pxgm1)
      
          CALL unbind_variable(vdf%atmo%states%search('temperature'))
          CALL bind_variable(vdf%atmo%states%search('temperature'), field%ta)
          CALL unbind_variable(vdf%atmo%inputs%list%Search('temperature'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('temperature'), field%ta)
          ins%ptm1 => ins%list%Get_ptr_r3d('temperature')
          __acc_attach(ins%ptm1)
      
          CALL unbind_variable(vdf%atmo%states%search('vertical velocity'))
          CALL bind_variable(vdf%atmo%states%search('vertical velocity'), field%wa)
          CALL unbind_variable(vdf%atmo%inputs%list%Search('vertical wind'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('vertical wind'), field%wa)
          ins%pwp1 => ins%list%Get_ptr_r3d('vertical wind')
          __acc_attach(ins%pwp1)
            ! CALL vdf%atmo%inputs%Set_pointers()

          CALL unbind_variable(vdf%atmo%inputs%list%Search('air density'))
          CALL bind_variable(vdf%atmo%inputs%list%Search('air density'), field%rho)
          ins%rho => ins%list%Get_ptr_r3d('air density')
          __acc_attach(ins%rho)
  
        END SELECT

        SELECT TYPE (ins => vdf%sfc%inputs)
        TYPE IS (t_vdf_sfc_inputs)

          ptr_r2d => field% qtrc_phy(:,nlev,:,iqv)
          CALL unbind_variable(vdf%sfc%inputs%list%Search('atm total water'))
          CALL bind_variable(vdf%sfc%inputs%list%Search('atm total water'), ptr_r2d)
          ins%qa => ins%list%Get_ptr_r2d('atm total water')
          __acc_attach(ins%qa)

          ptr_r2d => field%ta(:,nlev,:)
          CALL unbind_variable(vdf%sfc%inputs%list%Search('atm temperature'))
          CALL bind_variable(vdf%sfc%inputs%list%Search('atm temperature'), ptr_r2d)
          ins%ta => ins%list%Get_ptr_r2d('atm temperature')
          __acc_attach(ins%ta)

        END SELECT

        !
        ! Move forward one time step
        !
        CALL vdf%Compute(datetime_old)

        ! Retrieve computed tendency for temperature
        tend_ta_vdf => vdf%atmo%Get_tendency_r3d('temperature')
        ! Retrieve computed tendency for water vapor
        tend_qv_vdf => vdf%atmo%Get_tendency_r3d('water vapor')
        ! Retrieve computed tendency for cloud water
        tend_qc_vdf => vdf%atmo%Get_tendency_r3d('cloud water')
        ! Retrieve computed tendency for cloud ice
        tend_qi_vdf => vdf%atmo%Get_tendency_r3d('cloud ice')
        ! Retrieve computed tendencies for horizontal velocity
        tend_ua_vdf => vdf%atmo%Get_tendency_r3d('eastward wind')
        tend_va_vdf => vdf%atmo%Get_tendency_r3d('northward wind')
        ! Retrieve computed tendency for vertical velocity
        tend_wa_vdf => vdf%atmo%Get_tendency_r3d('vertical velocity')
        ! Retrieve computed tendency for surface temperature on tiles (for output only)
        tend_ts => vdf%sfc%Get_tendency_r3d('surface temperature')

        !$ACC DATA CREATE(q_rlw_impl, tend_ta_rlw_impl)

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,jk,jsfc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs, jbe

          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
      
          !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(field%qtrc_phy) ASYNC(1)
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = jcs, jce

              ! Add tendency from turbulent transport to physics tendency
              ! (is used in iconam_aes interface)
              tend%ta_phy(jc,jk,jb) = tend%ta_phy(jc,jk,jb) + tend_ta_vdf(jc,jk,jb)
              ! Update physics state
              field%ta(jc,jk,jb) = field%ta(jc,jk,jb) + tend_ta_vdf(jc,jk,jb) * dtime

              ! Add tendency from turbulent transport to physics tendency
              ! (is used in iconam_aes interface)
              tend%qtrc_phy (jc,jk,jb,iqv) = tend%qtrc_phy (jc,jk,jb,iqv) + tend_qv_vdf(jc,jk,jb)
              ! Update physics state
              field%qtrc_phy(jc,jk,jb,iqv) = field%qtrc_phy(jc,jk,jb,iqv) + tend_qv_vdf(jc,jk,jb) * dtime

              ! Add tendency from turbulent transport to physics tendency
              ! (is used in iconam_aes interface)
              tend%qtrc_phy(jc,jk,jb,iqc) = tend%qtrc_phy(jc,jk,jb,iqc) + tend_qc_vdf(jc,jk,jb)
              ! Update physics state
              field%qtrc_phy(jc,jk,jb,iqc) = field%qtrc_phy(jc,jk,jb,iqc) + tend_qc_vdf(jc,jk,jb) * dtime

              ! Add tendency from turbulent transport to physics tendency
              ! (is used in iconam_aes interface)
              tend%qtrc_phy(jc,jk,jb,iqi) = tend%qtrc_phy(jc,jk,jb,iqi) + tend_qi_vdf(jc,jk,jb)
              ! Update physics state
              field%qtrc_phy(jc,jk,jb,iqi) = field%qtrc_phy(jc,jk,jb,iqi) + tend_qi_vdf(jc,jk,jb) * dtime
              ! 
              tend%ua_phy(jc,jk,jb) = tend%ua_phy(jc,jk,jb) + tend_ua_vdf(jc,jk,jb)
              tend%va_phy(jc,jk,jb) = tend%va_phy(jc,jk,jb) + tend_va_vdf(jc,jk,jb)
              ! Update physics state
              field%ua(jc,jk,jb) = field%ua(jc,jk,jb) + tend_ua_vdf(jc,jk,jb) * dtime
              field%va(jc,jk,jb) = field%va(jc,jk,jb) + tend_va_vdf(jc,jk,jb) * dtime
              !
              tend%wa_phy(jc,jk,jb) = tend%wa_phy(jc,jk,jb) + tend_wa_vdf(jc,jk,jb)
              ! Update physics state
              field%wa(jc,jk,jb) = field%wa(jc,jk,jb) + tend_wa_vdf(jc,jk,jb) * dtime

            END DO
          END DO
          !$ACC END LOOP

          IF (ASSOCIATED(tend%ta_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                tend%ta_vdf(jc,jk,jb) = tend_ta_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF

          IF (ASSOCIATED(tend%ua_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                tend%ua_vdf(jc,jk,jb) = tend_ua_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF

          IF (ASSOCIATED(tend%va_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                tend%va_vdf(jc,jk,jb) = tend_va_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF

          IF (ASSOCIATED(tend%qtrc_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = jcs, jce
                tend%qtrc_vdf(jc,jk,jb,iqv) = tend_qv_vdf(jc,jk,jb)
                tend%qtrc_vdf(jc,jk,jb,iqc) = tend_qc_vdf(jc,jk,jb)
                tend%qtrc_vdf(jc,jk,jb,iqi) = tend_qi_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF
          !$ACC END PARALLEL

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          IF (ASSOCIATED(tend%wa_vdf)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO jk = 1, nlevp1
              DO jc = jcs, jce
                tend%wa_vdf(jc,jk,jb) = tend_wa_vdf(jc,jk,jb)
              END DO
            END DO
            !$ACC END LOOP
          END IF

          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jsfc = 1, nsfc_type
            DO jc = jcs, jce
              field%ts_tile(jc,jb,jsfc) = field%ts_tile(jc,jb,jsfc) + tend_ts(jc,jb,jsfc) * dtime
            END DO
          END DO
          !$ACC END LOOP

          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = jcs, jce

            q_rlw_impl(jc,jb) =                                               &
              &  ( (field%rld_rt(jc,nlev,jb)-field%rlu_rt(jc,nlev,jb))  & ! ( rln  from "radiation", at top of layer nlev
              &   -(field%rlds  (jc,jb)     -field%rlus  (jc,jb)     )) & !  -rlns from "radheating" and "update_surface")
              & -field%q_rlw_nlev(jc,jb)                                       ! -old heating in layer nlev from "radheating"

            ! Correction related to implicitness, due to the fact that surface model only used
            ! part of longwave radiation to compute new surface temperature
            !  
            IF (ASSOCIATED(field%q_rlw_impl)) THEN
              field%q_rlw_impl(jc,jb) = q_rlw_impl(jc,jb)
            END IF

            ! convert heating
            tend_ta_rlw_impl(jc,jb) = q_rlw_impl(jc,jb) / field%cvair(jc,nlev,jb)
            !
            IF (ASSOCIATED(tend%ta_rlw_impl)) THEN
              tend%ta_rlw_impl(jc,jb) = tend_ta_rlw_impl(jc,jb)
            END IF

            ! for output: accumulate heating
            IF (ASSOCIATED(field% q_phy)) THEN
              field%q_phy(jc,nlev,jb) = field%q_phy(jc,nlev,jb) + q_rlw_impl(jc,jb)
            END IF
            IF (ASSOCIATED(field% q_phy_vi)) THEN
              field%q_phy_vi(jc,jb) = field%q_phy_vi(jc,jb) + q_rlw_impl(jc,jb)
            END IF
            !
            ! accumulate surface lw increment for later updating the model state
            ! but only if radiative processes are active.
            IF (l_use_rad) THEN
              ! use tendency to update the model state
              tend%ta_phy(jc,nlev,jb) = tend%ta_phy(jc,nlev,jb) + tend_ta_rlw_impl(jc,jb)
              !
              ! update physics state for input to the next physics process
              ! use tendency to update the physics state
              field%ta(jc,nlev,jb) = field%ta(jc,nlev,jb) + tend_ta_rlw_impl(jc,jb) * dtime
            END IF

          ENDDO
          !$ACC END LOOP

          !$ACC END PARALLEL

        END DO

!$OMP END PARALLEL DO

        !$no EXIT DATA COPYOUT(tend_qtrc_vdf(iqv)%p)
        !$no EXIT DATA COPYOUT(tend_qtrc_vdf(iqc)%p)
        !$no EXIT DATA COPYOUT(tend_qtrc_vdf(iqi)%p)
        !$no EXIT DATA COPYOUT(tend_qtrc_vdf(1:3))

        !$ACC WAIT(1)

        !$ACC END DATA

        ! Retrieve diagnostics (mostly for output or restart only)
        ! IF (ASSOCIATED(field%q_vdf)) field%q_vdf(:,:,:) = vdf%atmo%Get_diagnostic_r3d('layer heating')

        ! field%cptgz(:,:,:) = vdf%atmo%Get_diagnostic_r3d('static energy')

        ! field%ts (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc temperature')
        ! !
        ! field%lhflx (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc latent heat flux')
        ! field%shflx (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc sensible heat flux')
        ! field%u_stress (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc zonal wind stress')
        ! field%v_stress (:,:) = vdf%sfc%Get_diagnostic_r2d('sfc mer. wind stress')
        ! !
        ! field%lhflx_tile (:,:,:) = vdf%sfc%Get_diagnostic_r3d('sfc latent heat flux, tile')
        ! field%shflx_tile (:,:,:) = vdf%sfc%Get_diagnostic_r3d('sfc sensible heat flux, tile')
        ! field%u_stress_tile (:,:,:) = vdf%sfc%Get_diagnostic_r3d('sfc zonal wind stress, tile')
        ! field%v_stress_tile (:,:,:) = vdf%sfc%Get_diagnostic_r3d('sfc mer. wind stress, tile')
        ! !
        ! field%cfm (:,:,:) = vdf%atmo%Get_diagnostic_r3d('exchange coefficient momentum')
        ! field%cfh (:,:,:) = vdf%atmo%Get_diagnostic_r3d('drag coefficient scalar')
        ! !
        ! field%z0m(:,:) = vdf%sfc%Get_diagnostic_r2d('roughness length momentum')
        ! field%z0h(:,:) = vdf%sfc%Get_diagnostic_r2d('roughness length heat')
      
      END IF
    END IF

    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_vdf)

  END SUBROUTINE interface_aes_tmx

  !New vdf
  SUBROUTINE init_tmx(p_patch, dtime)

    USE mo_variable, ONLY: bind_variable, t_variable
    USE mo_vdf,      ONLY: heat_type, momentum_type
    USE mo_vdf_atmo, ONLY: t_vdf_atmo_inputs
    ! USE mo_vdf_sfc,  ONLY: t_vdf_sfc_diagnostics
    USE mo_tmx_field_class, ONLY: isfc_oce, isfc_ice, isfc_lnd
    ! USE mo_vdf_diag_smag

    USE mo_nonhydro_state,     ONLY: p_nh_state
    USE mo_nonhydro_types,     ONLY: t_nh_metrics, t_nh_diag, t_nh_prog
    USE mo_dynamics_config,    ONLY: nnow
    USE mo_physical_constants, ONLY: cpd, cpv, cvd, cvv, Tf, tmelt, vmr_to_mmr_co2

    USE mo_master_config, ONLY: isRestart

    TYPE(t_patch), INTENT(inout), TARGET :: p_patch
    REAL(wp),      INTENT(in)    :: dtime

    TYPE(t_aes_phy_field),POINTER :: field
    TYPE(t_aes_phy_tend), POINTER :: tend
    TYPE(t_nh_metrics),   POINTER :: p_nh_metrics
    TYPE(t_nh_diag),      POINTER :: p_nh_diag
    TYPE(t_nh_prog),      POINTER :: p_nh_prog

    TYPE(t_vdf), POINTER :: vdf

    TYPE(t_patch), POINTER :: patch
    INTEGER :: nlev, nlevp1, jg
    INTEGER :: rls, rle, jbs, jbe, jcs, jce, jb, jc
    REAL(wp), POINTER :: ptr_r2d(:,:), ptr_r3d(:,:,:)
    REAL(sp), POINTER :: ptr_s2d(:,:), ptr_s3d(:,:,:)
    INTEGER, ALLOCATABLE :: sfc_types(:)

    REAL(wp), POINTER :: dz_srf(:,:)
    REAL(wp), POINTER :: zco2(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_tmx'

    patch => p_patch

    jg = patch%id

    nlev = patch%nlev
    nlevp1 = nlev + 1

    rls = grf_bdywidth_c + 1
    rle = min_rlcell_int

    jbs = patch%cells%start_block(rls)
    jbe = patch%cells%end_block  (rle)

    ! ALLOCATE(dz_srf(nproma,patch%nblks_c))
    
    ! TODO simple hack here
    ALLOCATE(zco2(nproma,patch%nblks_c))
    !$ACC ENTER DATA CREATE(zco2)
!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs, jbe

      CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
      DO jc = jcs, jce
        zco2(jc,jb) = 348.0e-06_wp * vmr_to_mmr_co2
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    field     => prm_field(jg)
    tend      => prm_tend (jg)
    p_nh_metrics => p_nh_state(jg)%metrics
    p_nh_diag => p_nh_state(jg)%diag
    p_nh_prog => p_nh_state(jg)%prog(nnow(jg))

    ! Question: use fields from AES field or from e.g. p_nh_state_lists(jg)%metrics p_nh_state_lists(jg)%diag? !!!!!!!!!!

    ! The order of sfc_types must be consistent with how variables are defined in aes memory!
    ALLOCATE(sfc_types(0))
    IF (iwtr <= nsfc_type) sfc_types = [sfc_types, isfc_oce]
    IF (iice <= nsfc_type) sfc_types = [sfc_types, isfc_ice]
    IF (ilnd <= nsfc_type) sfc_types = [sfc_types, isfc_lnd]

    vdf => new_vdf(patch, nproma, nlev=nlev, nsfc_tiles=nsfc_type, sfc_types=sfc_types, dt=dtime)
    __acc_attach(vdf)

    ! CALL vdf%atmo%Add_state('dry static energy', type=heat_type, field=field%cptgz)
    ! CALL vdf%atmo%Add_state('temperature', heat_type, dims=[nproma,nlev,patch%nblks_c])
    CALL vdf%atmo%Add_state('temperature',          type=heat_type,     field=field%ta)
    ptr_r3d => field% qtrc_phy(:,:,:,iqv)
    CALL vdf%atmo%Add_state('water vapor',          type=heat_type,     field=ptr_r3d)
    ptr_r3d => field% qtrc_phy(:,:,:,iqc)
    CALL vdf%atmo%Add_state('cloud water',          type=heat_type,     field=ptr_r3d)
    ptr_r3d => field% qtrc_phy(:,:,:,iqi)
    CALL vdf%atmo%Add_state('cloud ice',            type=heat_type,     field=ptr_r3d)
    CALL vdf%atmo%Add_state('eastward wind',        type=momentum_type, field=field%ua)
    CALL vdf%atmo%Add_state('northward wind',       type=momentum_type, field=field%va)
    CALL vdf%atmo%Add_state('vertical velocity',    type=momentum_type, field=field%wa)

    ! Bind variables to atmo config list
    CALL bind_variable(vdf%atmo%config%list%Search('cpd'), cpd)
    CALL bind_variable(vdf%atmo%config%list%Search('cvd'), cvd)
    CALL bind_variable(vdf%atmo%config%list%Search('Smagorinsky constant'), aes_vdf_config(jg)%smag_constant)
    CALL bind_variable(vdf%atmo%config%list%Search('maximum turbulence length scale'), aes_vdf_config(jg)%max_turb_scale)
    CALL bind_variable(vdf%atmo%config%list%Search('minimum Km'),aes_vdf_config(jg)%km_min)
    CALL bind_variable(vdf%atmo%config%list%Search('reverse prandtl number'),aes_vdf_config(jg)%rturb_prandtl)
    CALL bind_variable(vdf%atmo%config%list%Search('prandtl number'),aes_vdf_config(jg)%turb_prandtl)
    CALL bind_variable(vdf%atmo%config%list%Search('switch to activate Louis formula'),aes_vdf_config(jg)%use_louis)
    CALL bind_variable(vdf%atmo%config%list%Search('Louis constant b'),aes_vdf_config(jg)%louis_constant_b)
    CALL bind_variable(vdf%atmo%config%list%Search('time step'),dtime)
    CALL bind_variable(vdf%atmo%config%list%Search('solver type'), aes_vdf_config(jg)%solver_type)
    CALL bind_variable(vdf%atmo%config%list%Search('energy type'), aes_vdf_config(jg)%energy_type)
    CALL bind_variable(vdf%atmo%config%list%Search('dissipation factor'), aes_vdf_config(jg)%dissipation_factor)

    ! Bind variables to atmo input list
    ! 3d
    CALL bind_variable(vdf%atmo%inputs%list%Search('temperature'),         field%ta)
    CALL bind_variable(vdf%atmo%inputs%list%Search('zonal wind'),          field%ua)
    CALL bind_variable(vdf%atmo%inputs%list%Search('meridional wind'),     field%va)
    CALL bind_variable(vdf%atmo%inputs%list%Search('vertical wind'),       field%wa)
    CALL bind_variable(vdf%atmo%inputs%list%Search('virtual temperature'), field%tv)
    CALL bind_variable(vdf%atmo%inputs%list%Search('air density'),         field%rho)
    CALL bind_variable(vdf%atmo%inputs%list%Search('moist air mass'),      field%mair)
    CALL bind_variable(vdf%atmo%inputs%list%Search('specific heat of air at constant volume'),  field%cvair)
    CALL bind_variable(vdf%atmo%inputs%list%Search('geometric height full'), field%zf)
    CALL bind_variable(vdf%atmo%inputs%list%Search('geometric height half'), field%zh)
    !
    ptr_r3d => field% qtrc_phy(:,:,:,iqv)
    CALL bind_variable(vdf%atmo%inputs%list%Search('water vapor'), ptr_r3d)
    ptr_r3d => field% qtrc_phy(:,:,:,iqc)
    CALL bind_variable(vdf%atmo%inputs%list%Search('cloud water'), ptr_r3d)
    ! TODO: check for tracer existance
    ptr_r3d => field% qtrc_phy(:,:,:,iqi)
    CALL bind_variable(vdf%atmo%inputs%list%Search('cloud ice'), ptr_r3d)
    ptr_r3d => field% qtrc_phy(:,:,:,iqr)
    CALL bind_variable(vdf%atmo%inputs%list%Search('rain'), ptr_r3d)
    ptr_r3d => field% qtrc_phy(:,:,:,iqs)
    CALL bind_variable(vdf%atmo%inputs%list%Search('snow'), ptr_r3d)
    ptr_r3d => field% qtrc_phy(:,:,:,iqg)
    CALL bind_variable(vdf%atmo%inputs%list%Search('graupel'), ptr_r3d)
    !
    CALL bind_variable(vdf%atmo%inputs%list%Search('full level pressure'), field%pfull)
    CALL bind_variable(vdf%atmo%inputs%list%Search('half level pressure'), field%phalf)
    CALL bind_variable(vdf%atmo%inputs%list%Search('layer thickness'), field%dz)
    CALL bind_variable(vdf%atmo%inputs%list%Search('layer thickness full'), p_nh_metrics%ddqz_z_full)
    CALL bind_variable(vdf%atmo%inputs%list%Search('inverse layer thickness full'), p_nh_metrics%inv_ddqz_z_full)
    CALL bind_variable(vdf%atmo%inputs%list%Search('inverse layer thickness half'), p_nh_metrics%inv_ddqz_z_half)
#ifdef __MIXED_PRECISION
    ptr_s3d => p_nh_metrics%ddqz_z_half(:,:,:)
    CALL bind_variable(vdf%atmo%inputs%list%Search('layer thickness half'), ptr_s3d)
#else
    ptr_r3d => p_nh_metrics%ddqz_z_half(:,:,:)
    CALL bind_variable(vdf%atmo%inputs%list%Search('layer thickness half'), ptr_r3d)
#endif
    CALL bind_variable(vdf%atmo%inputs%list%Search('geopotential above groundlevel at interface and cell center'), p_nh_metrics%geopot_agl_ifc)
    ! 2d
    ! CALL bind_variable(vdf%atmo%inputs%list%Search('area fraction with wet land surface'), field%csat)
    ! CALL bind_variable(vdf%atmo%inputs%list%Search('area fraction with wet land surface (air)'), field%cair)

    !
    ! Surface
    !
    CALL vdf%sfc%Add_state('saturation specific humidity', type=heat_type, dims=[nproma,patch%nblks_c,SIZE(sfc_types)])
    CALL vdf%sfc%Add_state('surface temperature', type=heat_type, field=field%ts_tile)
    
    ! Bind variables to sfc config list
    CALL bind_variable(vdf%sfc%config%list%Search('cpd'), cpd)
    CALL bind_variable(vdf%sfc%config%list%Search('cvd'), cvd)
    CALL bind_variable(vdf%sfc%config%list%Search('cvv'), cvv)
    CALL bind_variable(vdf%sfc%config%list%Search('time step'),dtime)
    CALL bind_variable(vdf%sfc%config%list%Search('minimum surface wind speed'), aes_vdf_config(jg)%min_sfc_wind)
    CALL bind_variable(vdf%sfc%config%list%Search('ocean roughness length'),     aes_vdf_config(jg)%z0m_oce)
    CALL bind_variable(vdf%sfc%config%list%Search('ice roughness length'),       aes_vdf_config(jg)%z0m_ice)
    CALL bind_variable(vdf%sfc%config%list%Search('minimal roughness length'),   aes_vdf_config(jg)%z0m_min)
    CALL bind_variable(vdf%sfc%config%list%Search('weight for interpolation to surface_layer mid level'), aes_vdf_config(jg)%fsl)
    CALL bind_variable(vdf%sfc%config%list%Search('number of sea ice thickness classes'), field%kice)

    ! Bind variables to sfc input list
    ptr_r2d => field%ta(:,nlev,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm temperature'), ptr_r2d)
    ptr_r2d => field%tv(:,nlev,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm virtual temperature'), ptr_r2d)
    ptr_r2d => field%ua(:,nlev,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm zonal wind'),     ptr_r2d)
    ptr_r2d => field%va(:,nlev,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm meridional wind'), ptr_r2d)
    ptr_r2d => field% qtrc_phy(:,nlev,:,iqv)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm total water'), ptr_r2d)
    ptr_r2d => field%rho(:,nlev,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm density'), ptr_r2d)
    ptr_r2d => field%pfull(:,nlev,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm full level pressure'), ptr_r2d)
    ptr_r2d => field%phalf(:,nlevp1,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('surface pressure'), ptr_r2d)
    ! CALL bind_variable(vdf%sfc%inputs%list%Search('surface pressure'), p_nh_diag%pres_sfc)
    ptr_r2d => field%zf(:,nlev,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm geometric height full'), ptr_r2d)
    ptr_r2d => field%zh(:,nlevp1,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('sfc geometric height half'), ptr_r2d)
    !
    CALL bind_variable(vdf%sfc%inputs%list%Search('surface rain flux, large-scale'), field%rsfl)
    CALL bind_variable(vdf%sfc%inputs%list%Search('surface snow flux, large-scale'), field%ssfl)
    !
    CALL bind_variable(vdf%sfc%inputs%list%Search('surface downward longwave radiation'),                field%rlds)
    CALL bind_variable(vdf%sfc%inputs%list%Search('surface downward shortwave radiation'),               field%rsds)
    CALL bind_variable(vdf%sfc%inputs%list%Search('all-sky surface downward direct visible radiation'),  field%rvds_dir)
    CALL bind_variable(vdf%sfc%inputs%list%Search('all-sky surface downward direct near-IR radiation'),  field%rnds_dir)
    CALL bind_variable(vdf%sfc%inputs%list%Search('all-sky surface downward direct PAR radiation'),      field%rpds_dir)
    CALL bind_variable(vdf%sfc%inputs%list%Search('all-sky surface downward diffuse visible radiation'), field%rvds_dif)
    CALL bind_variable(vdf%sfc%inputs%list%Search('all-sky surface downward diffuse near-IR radiation'), field%rnds_dif)
    CALL bind_variable(vdf%sfc%inputs%list%Search('all-sky surface downward diffuse PAR radiation'),     field%rpds_dif)
    !
    CALL bind_variable(vdf%sfc%inputs%list%Search('surface temperature, tile'), field%ts_tile)
    CALL bind_variable(vdf%sfc%inputs%list%Search('grid box fraction of tiles'), field%frac_tile)
    !
    ! TODO: lw surface emissivity should be tile-specific and, for land, should be returned from land model
    CALL bind_variable(vdf%sfc%inputs%list%Search('longwave surface emissivity'), field%emissivity)
    CALL bind_variable(vdf%sfc%inputs%list%Search('u-component of ocean current'), field%ocu)
    CALL bind_variable(vdf%sfc%inputs%list%Search('v-component of ocean current'), field%ocv)
    !
    ptr_r2d => field%hi(:,1,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('thickness of sea ice'), ptr_r2d)
    !
    CALL bind_variable(vdf%sfc%inputs%list%Search('cosine of zenith angle'), field%cosmu0)
    CALL bind_variable(vdf%sfc%inputs%list%Search('atm CO2 concentration'), zco2) ! TODO carbon cycle
    !
#ifdef __MIXED_PRECISION
    ptr_s2d => p_nh_metrics%ddqz_z_half(:,nlevp1,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('reference height in surface layer times 2'), ptr_s2d)
#else
    ptr_r2d => p_nh_metrics%ddqz_z_half(:,nlevp1,:)
    CALL bind_variable(vdf%sfc%inputs%list%Search('reference height in surface layer times 2'), ptr_r2d)
#endif
    ! dz_srf(:,:) = 2._wp * (field%zh(:,nlev,:) - field%zh(:,nlevp1,:))
    ! dz_srf(:,:) = 2._wp * (p_nh_metrics%z_mc(:,nlev,:) - p_nh_metrics%z_ifc(:,nlevp1,:))
    ! CALL bind_variable(vdf%sfc%inputs%list%Search('reference height in surface layer times 2'), dz_srf)
    !
    CALL bind_variable(vdf%atmo%diagnostics%list%Search('static energy'),                               field%cptgz)
    IF (ASSOCIATED(field%cptgzvi)) THEN
      CALL bind_variable(vdf%atmo%diagnostics%list%Search('static energy, vert. int.'),                 field%cptgzvi)
    END IF
    IF (ASSOCIATED(field%q_vdf)) THEN
      CALL bind_variable(vdf%atmo%diagnostics%list%Search('layer heating'),                             field%q_vdf)
    END IF
    IF (ASSOCIATED(field%kedisp)) THEN
      CALL bind_variable(vdf%atmo%diagnostics%list%Search('dissipation of kinetic energy, vert. int.'), field%kedisp)
    END IF
    IF (ASSOCIATED(field%utmxvi)) THEN
      CALL bind_variable(vdf%atmo%diagnostics%list%Search('moist internal energy after tmx, vert. int.'), field%utmxvi)
    END IF
    IF (ASSOCIATED(tend%utmxvi)) THEN
      CALL bind_variable(vdf%atmo%diagnostics%list%Search('tendency of vert. int. moist internal energy'), tend%utmxvi)
    END IF
    CALL bind_variable(vdf%atmo%diagnostics%list%Search('exchange coefficient momentum full'),          field%cfm)
    CALL bind_variable(vdf%atmo%diagnostics%list%Search('exchange coefficient scalar full'),            field%cfh)
    !
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc temperature'),                              field%ts)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc radiative temperature'),                    field%ts_rad)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc longwave net flux, tile'),                  field%lwflxsfc_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc shortwave net flux, tile'),                 field%swflxsfc_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc longwave upward flux'),                     field%rlus)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc shortwave upward flux'),                    field%rsus)
    !
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc evapotranspiration'),                       field%evap)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc latent heat flux'),                         field%lhflx)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc sensible heat flux'),                       field%shflx)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('energy flux at surface from thermal exchange'), field%ufts)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('energy flux at surface from vapor exchange'),   field%ufvs)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc zonal wind stress'),                        field%u_stress)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc mer. wind stress'),                         field%v_stress)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc evapotranspiration, tile'),                 field%evap_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc latent heat flux, tile'),                   field%lhflx_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc sensible heat flux, tile'),                 field%shflx_tile) 
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc zonal wind stress, tile'),                  field%u_stress_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('sfc mer. wind stress, tile'),                   field%v_stress_tile)
    IF (ASSOCIATED(field%z0m)) &
      & CALL bind_variable(vdf%sfc%diagnostics%list%Search('roughness length momentum'),                    field%z0m)
    IF (ASSOCIATED(field%z0h)) &
      & CALL bind_variable(vdf%sfc%diagnostics%list%Search('roughness length heat'),                        field%z0h)
    IF (ASSOCIATED(field%z0m_tile)) &
      & CALL bind_variable(vdf%sfc%diagnostics%list%Search('roughness length momentum, tile'),              field%z0m_tile)
    IF (ASSOCIATED(field%z0h_tile)) &
      & CALL bind_variable(vdf%sfc%diagnostics%list%Search('roughness length heat, tile'),                  field%z0h_tile)
    IF (ASSOCIATED(field%cfm_tile)) &
      & CALL bind_variable(vdf%sfc%diagnostics%list%Search('exchange coefficient for momentum, tile'),  field%cfm_tile)
    IF (ASSOCIATED(field%cfh_tile)) &
      & CALL bind_variable(vdf%sfc%diagnostics%list%Search('exchange coefficient for scalar, tile'),    field%cfh_tile)
    !
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo VIS direct, tile'),                      field%albvisdir_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo VIS diffuse, tile'),                     field%albvisdif_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo NIR direct, tile'),                      field%albnirdir_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo NIR diffuse, tile'),                     field%albnirdif_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo, tile'),                                 field%albedo_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo VIS direct'),                            field%albvisdir)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo VIS diffuse'),                           field%albvisdif)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo NIR direct'),                            field%albnirdir)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo NIR diffuse'),                           field%albnirdif)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('albedo'),                                       field%albedo)

    ptr_r2d => field%Qtop(:,1,:)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('energy flux available for surface melting of sea ice'), ptr_r2d)
    ptr_r2d => field%Qbot(:,1,:)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('energy flux at ice-ocean interface'),                   ptr_r2d)
    ptr_r2d => field%hs(:,1,:)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('thickness of snow on sea ice'),                         ptr_r2d)

    CALL bind_variable(vdf%sfc%diagnostics%list%Search('2m temperature'),                               field%tas)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('2m temperature, tile'),                         field%tas_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('2m specific humidity'),                         field%qv2m)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('2m specific humidity, tile'),                   field%qv2m_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('2m dewpoint temperature'),                      field%dew2)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('2m dewpoint temperature, tile'),                field%dew2_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('10m wind speed'),                               field%sfcwind)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('10m zonal wind'),                               field%uas)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('10m meridional wind'),                          field%vas)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('10m wind speed, tile'),                         field%sfcwind_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('10m zonal wind, tile'),                         field%uas_tile)
    CALL bind_variable(vdf%sfc%diagnostics%list%Search('10m meridional wind, tile'),                    field%vas_tile)

    CALL vdf%Lock_variable_sets()

    vdf_dom(jg)%p => vdf
    
    CALL vdf%atmo%Init()

  END SUBROUTINE init_tmx

END MODULE mo_interface_aes_tmx
