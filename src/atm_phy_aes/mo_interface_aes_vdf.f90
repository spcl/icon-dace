!
! Subroutine interface_aes_vdf calls the vertical diffusion and the surface schemes.
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

MODULE mo_interface_aes_vdf

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: t_datetime => datetime, OPERATOR(>)

  USE mo_exception           ,ONLY: finish

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: ntracer

  USE mo_aes_phy_config      ,ONLY: aes_phy_config, aes_phy_tc, dt_zero
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field, &
    &                               t_aes_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_vdf, &
    &                               timer_vdf_dn, timer_vdf_sf, timer_vdf_up

  USE mo_ccycle_config       ,ONLY: ccycle_config
  USE mo_physical_constants  ,ONLY: vmr_to_mmr_co2
  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2mmr

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqnc, iqni, iqt, ico2
  USE mo_turb_vdiff          ,ONLY: vdiff_down, vdiff_up, nvar_vdiff, nmatrix, imh, &
    &                               imqv, ih_vdiff=>ih, iqv_vdiff=>iqv
  USE mo_turb_vdiff_params   ,ONLY: VDIFF_TURB_3DSMAGORINSKY

  USE mo_aes_sfc_indices     ,ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_surface             ,ONLY: update_surface
  USE mo_surface_diag        ,ONLY: nsurf_diag
  USE mo_aes_vdf_config      ,ONLY: aes_vdf_config
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_impl_constants_grf  ,ONLY: grf_bdywidth_c
  USE mo_impl_constants      ,ONLY: min_rlcell_int
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_nh_testcases_nml    ,ONLY: is_dry_cbl, isrfc_type


  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_aes_vdf

CONTAINS

  SUBROUTINE interface_aes_vdf(patch)

    ! Arguments
    !
    TYPE(t_patch)   ,TARGET ,INTENT(inout) :: patch

    ! Shortcuts
    !
    LOGICAL :: l2moment, ldtrad_gt0
    INTEGER :: fc_vdf
    TYPE(t_aes_phy_field)   ,POINTER    :: field
    TYPE(t_aes_phy_tend)    ,POINTER    :: tend

    ! Local variables
    !
    TYPE(t_datetime), POINTER :: datetime
    REAL(wp) :: pdtime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER  :: jg, jl, jk, jsfc, jt, jice
    INTEGER  :: jb,jbs,jbe,jcs,jce,ncd,rls,rle
    INTEGER  :: nlev, nlevm1, nlevp1
    INTEGER  :: ntrac
    INTEGER  :: nice        ! for simplicity (ice classes)
    !
    REAL(wp) :: zxt_emis(nproma,ntracer-iqt+1,patch%nblks_c)   !< tracer tendency due to surface emission

    REAL(wp) :: zcpt_sfc_tile(nproma,patch%nblks_c,nsfc_type)  !< dry static energy at surface

    ! Coefficient matrices and right-hand-side vectors for the turbulence solver
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)
    REAL(wp) :: ri_tile(nproma,patch%nblks_c,nsfc_type)           !< Richardson number
    REAL(wp) :: zaa    (nproma,patch%nlev,3,nmatrix,patch%nblks_c)      !< coeff. matrices, all variables
    REAL(wp) :: zaa_btm(nproma,3,nsfc_type,imh:imqv,patch%nblks_c)!< last row of coeff. matrix of heat and moisture
    REAL(wp) :: zbb    (nproma,patch%nlev,nvar_vdiff,patch%nblks_c)             !< r.h.s., all variables
    REAL(wp) :: zbb_btm(nproma,nsfc_type,ih_vdiff:iqv_vdiff,patch%nblks_c)!< last row of r.h.s. of heat and moisture

    ! Temporary arrays used by VDIFF, JSBACH
    REAL(wp) :: zfactor_sfc(nproma,patch%nblks_c)
    REAL(wp) :: zco2       (nproma,patch%nblks_c)          !< co2 value passed on to jsbach

    REAL(wp) :: zqx      (nproma,patch%nlev,patch%nblks_c)       !< total cloud condensate
    REAL(wp) :: zthvvar  (nproma,patch%nlev,patch%nblks_c)       !< intermediate value of thvvar
    REAL(wp) :: dummy    (nproma,patch%nlev,patch%nblks_c)       !< to replace thvvar
    REAL(wp) :: dummyx   (nproma,patch%nlev,patch%nblks_c)       !< to replace xvar
    REAL(wp) :: ztottevn (nproma,patch%nlev,patch%nblks_c)       !< intermediate value of TTE
    REAL(wp) :: zch_tile (nproma,patch%nblks_c,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zchn_tile(nproma,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zcdn_tile(nproma,nsfc_type)  !< for "nsurf_diag"
!!$    REAL(wp) :: zcfnc_tile(nproma,nsfc_type) !< for "nsurf_diag"
    REAL(wp) :: zbn_tile (nproma,patch%nblks_c,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbhn_tile(nproma,patch%nblks_c,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbm_tile (nproma,patch%nblks_c,nsfc_type)  !< for "nsurf_diag"
    REAL(wp) :: zbh_tile (nproma,patch%nblks_c,nsfc_type)  !< for "nsurf_diag"

    ! inout variables of vdiff
    !
    REAL(wp) :: ts_tile(nproma,patch%nblks_c,nsfc_type)
    REAL(wp) :: z0m_tile(nproma,patch%nblks_c,nsfc_type)

    ! inout variables of vdiff_down
    !
    REAL(wp) :: ustar     (nproma,patch%nblks_c)
    REAL(wp) :: wstar_tile(nproma,patch%nblks_c,nsfc_type)

    ! output variables of vdiff_down
    !
    REAL(wp) :: wstar      (nproma,patch%nblks_c)
    REAL(wp) :: qs_sfc_tile(nproma,patch%nblks_c,nsfc_type)
    REAL(wp) :: hdtcbl     (nproma,patch%nblks_c)
    REAL(wp) :: ri_atm     (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: mixlen     (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: cfm        (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: cfm_tile   (nproma,patch%nblks_c,nsfc_type)
    REAL(wp) :: cfh        (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: cfh_tile   (nproma,patch%nblks_c,nsfc_type)
    REAL(wp) :: cfv        (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: cftotte    (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: cfthv      (nproma,patch%nlev,patch%nblks_c)

    ! inout variables of update_surface
    !
    REAL(wp) :: rlus            (nproma,patch%nblks_c)
    REAL(wp) :: albvisdir_ice   (nproma,prm_field(patch%id)%kice,patch%nblks_c)
    REAL(wp) :: albnirdir_ice   (nproma,prm_field(patch%id)%kice,patch%nblks_c)
    REAL(wp) :: albvisdif_ice   (nproma,prm_field(patch%id)%kice,patch%nblks_c)
    REAL(wp) :: albnirdif_ice   (nproma,prm_field(patch%id)%kice,patch%nblks_c)

    ! output variables of update_surface
    !
    REAL(wp) :: q_snocpymlt     (nproma,patch%nblks_c)
    REAL(wp) :: tend_ta_sfc     (nproma,patch%nblks_c)
    REAL(wp) :: q_rlw_impl      (nproma,patch%nblks_c)
    REAL(wp) :: tend_ta_rlw_impl(nproma,patch%nblks_c)
    REAL(wp) :: ts              (nproma,patch%nblks_c)
    REAL(wp) :: ts_rad          (nproma,patch%nblks_c)
    REAL(wp) :: evap            (nproma,patch%nblks_c)
    REAL(wp) :: cair            (nproma,patch%nblks_c)
    REAL(wp) :: csat            (nproma,patch%nblks_c)
    REAL(wp) :: z0h_lnd         (nproma,patch%nblks_c)

    REAL(wp) :: albvisdir       (nproma,patch%nblks_c)
    REAL(wp) :: albnirdir       (nproma,patch%nblks_c)
    REAL(wp) :: albvisdif       (nproma,patch%nblks_c)
    REAL(wp) :: albnirdif       (nproma,patch%nblks_c)
    REAL(wp) :: albvisdir_tile  (nproma,patch%nblks_c,nsfc_type)
    REAL(wp) :: albnirdir_tile  (nproma,patch%nblks_c,nsfc_type)
    REAL(wp) :: albvisdif_tile  (nproma,patch%nblks_c,nsfc_type)
    REAL(wp) :: albnirdif_tile  (nproma,patch%nblks_c,nsfc_type)
    REAL(wp) :: albedo          (nproma,patch%nblks_c)
    REAL(wp) :: albedo_tile     (nproma,patch%nblks_c,nsfc_type)

    ! output variables of vdiff_up
    !
    REAL(wp) :: kedisp       (nproma,patch%nblks_c)
    REAL(wp) :: q_vdf        (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: tend_ta_vdf  (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: tend_ua_vdf  (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: tend_va_vdf  (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: tend_wa_vdf  (nproma,patch%nlev+1,patch%nblks_c)
    REAL(wp), TARGET :: tend_qtrc_vdf(nproma,patch%nlev,patch%nblks_c,ntracer)
    REAL(wp), POINTER, CONTIGUOUS :: tend_qtrc_vdf_iqt(:,:,:,:)
    REAL(wp), TARGET :: tend_qtrc_vdf_dummy(nproma,patch%nlev,patch%nblks_c,0)
    REAL(wp) :: mmr_co2

    REAL(wp) :: ddt_u(nproma,patch%nlev,patch%nblks_c), ddt_v(nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: ta_hori_tend(nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: qv_hori_tend(nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: ql_hori_tend(nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: qi_hori_tend(nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: qnc_hori_tend(nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: qni_hori_tend(nproma,patch%nlev,patch%nblks_c)

    INTEGER, POINTER :: turb

    IF (ltimer) CALL timer_start(timer_vdf)

    jg   = patch%id
    nlev = patch%nlev

    datetime             => aes_phy_tc(jg)%datetime
    pdtime               =  aes_phy_tc(jg)%dt_phy_sec
    is_in_sd_ed_interval =  aes_phy_tc(jg)%is_in_sd_ed_interval_vdf
    is_active            =  aes_phy_tc(jg)%is_active_vdf

    ! associate pointers
    l2moment  =  aes_phy_config(jg)%l2moment
    ldtrad_gt0=  aes_phy_tc(jg)%dt_rad > dt_zero
    fc_vdf    =  aes_phy_config(jg)%fc_vdf
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    nlevm1 = nlev-1
    nlevp1 = nlev+1
    ntrac  = ntracer-iqt+1  ! number of tracers excluding water vapour and hydrometeors

    nice   = prm_field(jg)%kice
    turb => aes_vdf_config(jg)%turb

    !$ACC DATA CREATE(zxt_emis) IF(ntrac > 0)
    !$ACC DATA CREATE(ta_hori_tend, qv_hori_tend, ql_hori_tend, qi_hori_tend) IF(turb == 2)
    !$ACC DATA PRESENT(field, tend, ccycle_config) &
    !$ACC   PRESENT(field%qtrc_phy) & ! ACCWA (nvhpc on levante): to prevent illegal address during kernel execution
    !$ACC   CREATE(zcpt_sfc_tile, ri_tile, zqx, zbn_tile) &
    !$ACC   CREATE(zbhn_tile, zbm_tile, zbh_tile, dummy, dummyx) &
    !$ACC   CREATE(wstar, qs_sfc_tile, hdtcbl, ri_atm, mixlen, cfm) &
    !$ACC   CREATE(cfm_tile, cfh, cfh_tile, cfv, cftotte, cfthv, zaa) &
    !$ACC   CREATE(zaa_btm, zbb, zbb_btm, zfactor_sfc, ddt_u, ddt_v) &
    !$ACC   CREATE(zthvvar, ztottevn, zch_tile, kedisp, tend_ua_vdf) &
    !$ACC   CREATE(tend_va_vdf, tend_wa_vdf, q_vdf, tend_qtrc_vdf, q_snocpymlt, zco2) &
#if !(defined(_CRAYFTN) && _RELEASE_MAJOR <= 16)
    !ACCWA zero sized arrays are not properly supported (CAST-33010) and cause here 
    ! weird "present" error
    !$ACC   CREATE(tend_qtrc_vdf_dummy) &
#endif
    !$ACC   CREATE(tend_ta_sfc, q_rlw_impl, tend_ta_rlw_impl, tend_ta_vdf) &
    !$ACC   CREATE(ts_tile, z0m_tile, ustar, wstar_tile, rlus) &
    !$ACC   CREATE(albvisdir_ice, albnirdir_ice, albvisdif_ice) &
    !$ACC   CREATE(albnirdif_ice, ts, ts_rad, evap, cair, csat, z0h_lnd) &
    !$ACC   CREATE(albvisdir, albnirdir, albvisdif, albnirdif) &
    !$ACC   CREATE(albvisdir_tile, albnirdir_tile, albvisdif_tile) &
    !$ACC   CREATE(albnirdif_tile, albedo, albedo_tile) &
    !$ACC   CREATE(qnc_hori_tend, qni_hori_tend)

    IF ( is_dry_cbl ) THEN
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
      field% qtrc_phy(:,:,:,iqv) = 0._wp
      field% qtrc_phy(:,:,:,iqi) = 0._wp
      field% qtrc_phy(:,:,:,iqc) = 0._wp
      !$ACC END KERNELS
    END IF

    !$NOser verbatim zaa = 0._wp
    !$NOser verbatim zbb = 0._wp
    !$NOser verbatim !$ACC UPDATE DEVICE( zaa, zbb ) ASYNC(1)

!!$    ! Emission of aerosols or other tracers (not implemented yet)
!!$    IF (ntrac>0) THEN
!!$       CALL tracer_emission()
!!$    ENDIF
!!$    !
!!$    ! Dry deposition of aerosols or other tracers (not implemented yet)
!!$    CALL dry_deposition()


    IF ( is_in_sd_ed_interval ) THEN
      !
      IF ( is_active ) THEN

        rls = grf_bdywidth_c+1
        rle = min_rlcell_int

        ncd = MAX(1,patch%n_childdom)
        jbs = patch%cells%start_blk(rls,  1)
        jbe = patch%cells%end_blk  (rle,ncd)

!$OMP PARALLEL DO PRIVATE(jb,jl,jk,jt,jsfc,jice,jcs,jce,mmr_co2)
        !## jb loop1, is_in_sd_ed_interval == T .and. is_active == T
        DO jb = jbs, jbe
          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
          IF (jcs>jce) CYCLE
          !
          ! Set dummy values to zero to prevent invalid floating point operations:
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              dummy (jl,jk,jb) = 0._wp
              dummyx(jl,jk,jb) = 0._wp
            END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
          ! Surface emissions of extra tracers if ntrac > 0
          !
          ! - default is no emission
          IF (ntrac > 0) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(ntrac > 0)
            DO jt = 1,ntrac
              DO jl = jcs,jce
                zxt_emis(jl,jt,jb) = 0._wp
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          !
          ! Carbon cycle
          ! ------------
          ! (see also after vdiff_up)
          ! If the c-cycle is used then the surface co2 flux and the
          ! atmospheric co2 concentration must be prepared for the
          ! vertical diffusion.
          !
          IF (ccycle_config(jg)%iccycle /= 0 .AND. (ico2 < iqt .OR. ntracer < ico2)) THEN
             CALL finish('interface_aes_vdf','The C-cycle cannot be used without CO2 tracer (ico2<iqt or ntracer<ico2)')
          END IF
          !
          ! DA: fuse all the 1D copies in a single ACC kernel
          !$ACC PARALLEL NO_CREATE(zxt_emis) DEFAULT(PRESENT) ASYNC(1)
          SELECT CASE (ccycle_config(jg)%iccycle)
             !
          CASE (0) ! no c-cycle
             !
             ! This should have no meaning, but strangely changing zco2(:) changes the result
             ! of interface_aes_vdf/update_surface/jsbach.
             ! The original parameter aes_rad_config(jg)%vmr_co2 is replaced by its value:
             !
             !$ACC LOOP GANG VECTOR
             DO jl = jcs,jce
               zco2(jl,jb) = 348.0e-06_wp * vmr_to_mmr_co2
             END DO
             !
          CASE (1) ! c-cycle with interactive atm. co2 concentration
             !
             ! co2 flux at the surface
             !$ACC LOOP GANG VECTOR
             DO jl = jcs,jce
               zxt_emis(jl,ico2-(iqt-1),jb) = field%fco2nat(jl,jb)
             END DO
             !
             ! co2 concentration in the lowermost layer
             !$ACC LOOP GANG VECTOR
             DO jl = jcs,jce
               zco2(jl,jb) = field% qtrc_phy(jl,nlev,jb,ico2)
             END DO
             !
          CASE (2) ! c-cycle with prescribed atm. co2 concentration
             !
             ! co2 flux at the surface
             !$ACC LOOP GANG VECTOR
             DO jl = jcs,jce
               zxt_emis(jl,ico2-(iqt-1),jb) = field%fco2nat(jl,jb)
             END DO
             !
             ! co2 concentration in the lowermost layer
             SELECT CASE (ccycle_config(jg)%ico2conc)
             CASE (2)
                mmr_co2 = ccycle_config(jg)%vmr_co2 * vmr_to_mmr_co2
                !
                !$ACC LOOP GANG VECTOR
                DO jl = jcs,jce
                  zco2(jl,jb) = mmr_co2
                END DO
             CASE (4)
                !$ACC LOOP GANG VECTOR
                DO jl = jcs,jce
                  zco2(jl,jb) = ghg_co2mmr
                END DO
             END SELECT
             !
          END SELECT
          !$ACC END PARALLEL
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              zqx(jl,jk,jb) =  field%qtrc_phy(jl,jk,jb,iqc) + field%qtrc_phy(jl,jk,jb,iqi)
            END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO jl = jcs,jce
            ustar   (jl,jb) = field% ustar   (jl,jb)
            cair    (jl,jb) = field% cair    (jl,jb)
            csat    (jl,jb) = field% csat    (jl,jb)
            z0h_lnd (jl,jb) = field% z0h_lnd (jl,jb)
            rlus    (jl,jb) = field% rlus    (jl,jb)
          END DO
          !$ACC END PARALLEL LOOP
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jsfc=1,nsfc_type
            DO jl=jcs,jce
              ts_tile   (jl,jb,jsfc) = field% ts_tile   (jl,jb,jsfc)
              wstar_tile(jl,jb,jsfc) = field% wstar_tile(jl,jb,jsfc)
              z0m_tile  (jl,jb,jsfc) = field% z0m_tile  (jl,jb,jsfc)
            END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jice=1,nice
            DO jl=jcs,jce
              albvisdir_ice (jl,jice,jb) = field% albvisdir_ice(jl,jice,jb)
              albnirdir_ice (jl,jice,jb) = field% albnirdir_ice(jl,jice,jb)
              albvisdif_ice (jl,jice,jb) = field% albvisdif_ice(jl,jice,jb)
              albnirdif_ice (jl,jice,jb) = field% albnirdif_ice(jl,jice,jb)
            END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
        END DO !## jb loop 1
        !$ACC WAIT
!$OMP END PARALLEL DO

        IF (ltimer) CALL timer_start(timer_vdf_dn)
        !
        ! Turbulent mixing, part I:
        ! - computation of exchange coefficients in the atmosphere and at the surface;
        ! - build up the tridiagonal linear algebraic system;
        ! - downward sweep (Gaussian elimination from top till level nlev-1)
        !
        !----------------------------------------------------------------------------------------

        CALL vdiff_down( nproma, patch%nblks_c,             &! in
              &          patch%nblks_v, patch%nblks_e,      &! in
              &          nlev, nlevm1, nlevp1,              &! in
              &          ntrac, nsfc_type,                  &! in
              &          iwtr, iice, ilnd,                  &! in, indices of different surface types
              &          pdtime,                            &! in, time step
              &          field%coriol(:,:),                 &! in, Coriolis parameter
              ! Smagorinsky
              &          patch,                             &! in
              &          l2moment,                          &! in, l2moment switch for horizontal tendencies
              !
              &          field%   zf(:,:,:),                &! in, geopot. height above sea level, full level
              &          field%   zh(:,:,:),                &! in, geopot. height above sea level, half level
              &          field%frac_tile(:,:,:),            &! in, area fraction of each sfc type
              &                 ts_tile(:,:,:),             &! in, surface temperature
              &          field% ocu (:,:),                  &! in, ocean sfc velocity, u-component
              &          field% ocv (:,:),                  &! in, ocean sfc velocity, v-component
              &          field% phalf(:,nlevp1,:),          &! in, sfc pressure
              &          field%   ua(:,:,:),                &! in, um1
              &          field%   va(:,:,:),                &! in, vm1
              &          field%   wa(:,:,:),                &! in, wp1
              &          field%   ta(:,:,:),                &! in, tm1
              &          field% qtrc_phy(:,:,:,iqv),        &! in, qm1
              &          field% qtrc_phy(:,:,:,iqc),        &! in, xlm1
              &          field% qtrc_phy(:,:,:,iqi),        &! in, xim1
              &                  zqx(:,:,:),                &! in, xlm1 + xim1
              &          field% qtrc_phy(:,:,:,iqt:),       &! in, xtm1
              &          field% mair(:,:,:),                &! in, moist     air mass
              &          field% rho(:,:,:),                 &! in, air density
              &          field% phalf(:,:,:),               &! in, aphm1
              &          field% pfull(:,:,:),               &! in, apm1
              &          field%   tv(:,:,:),                &! in, virtual temperaturea
              &          field% aclc(:,:,:),                &! in, cloud fraction
              &          zxt_emis(:,:,:),                   &! in, zxtems
              &          dummy(:,:,:),                      &! in, variance of theta_v at step t-dt
              &          dummyx(:,:,:),                     &! in
              &                 z0m_tile(:,:,:),            &! in
              &          field% tottem1(:,:,:),             &! in, TTE at step t-dt
              &          aes_vdf_config(jg),                &! in
              &                 ustar  (:,:),               &! inout
              &                 wstar  (:,:),               &! out, convective velocity scale
              &                 wstar_tile(:,:,:),          &! inout, convective velocity scale (each sfc type)
              &                 qs_sfc_tile(:,:,:),         &! out, sfc specific humidity at saturation
              &                 hdtcbl  (:,:),              &! out, for output
              &                 ri_atm  (:,:,:),            &! out, for output
              &                 ri_tile (:,:,:),            &! out, for nsurf_diag
              &                 mixlen  (:,:,:),            &! out, for output
              &                 cfm     (:,:,:),            &! out, for output
              &                 cfm_tile(:,:,:),            &! out, for output and "vdiff_up"
              &                 cfh     (:,:,:),            &! out, for output
              &                 cfh_tile(:,:,:),            &! out, for output and "vdiff_up"
              &                 cfv     (:,:,:),            &! out, for output
              &                 cftotte (:,:,:),            &! out, for output
              &                 cfthv   (:,:,:),            &! out, for output
              &          zaa(:,:,:,:,:), zaa_btm(:,:,:,:,:),&
              &          zbb(:,:,:,:), zbb_btm(:,:,:,:),    &! out, for "vdiff_up"
              ! for Smagorinsky
              &          ddt_u(:,:,:), ddt_v(:,:,:),        &! out
              &          tend_wa_vdf(:,:,:),                &! out
              &          ta_hori_tend(:,:,:),               &! out
              &          qv_hori_tend(:,:,:),               &! out
              &          ql_hori_tend(:,:,:),               &! out
              &          qi_hori_tend(:,:,:),               &! out
              ! for 2 moment scheme and Smagorinsky
              &          qnc_hori_tend(:,:,:),              &! out
              &          qni_hori_tend(:,:,:),              &! out
              !
              &          zfactor_sfc(:,:),                  &! out, for "vdiff_up"
              &          zcpt_sfc_tile(:,:,:),              &! out, for "vdiff_up"
              &          field%cptgz(:,:,:),                &! out, for "vdiff_up"
              &          zthvvar(:,:,:),                    &! out, for "vdiff_up"
              &          ztottevn (:,:,:),                  &! out, for "vdiff_up"
              &          zch_tile(:,:,:),                   &! out, for "nsurf_diag"
!!$               &          zchn_tile(:,:),                  &! out, for "nsurf_diag"
!!$               &          zcdn_tile(:,:),                  &! out, for "nsurf_diag"
!!$               &          zcfnc_tile(:,:),                 &! out, for "nsurf_diag"
              &          zbn_tile(:,:,:),                   &! out, for "nsurf_diag"
              &          zbhn_tile(:,:,:),                  &! out, for "nsurf_diag"
              &          zbm_tile(:,:,:),                   &! out, for "nsurf_diag"
              &          zbh_tile(:,:,:),                   &! out, for "nsurf_diag"
              &          pcsat = csat(:,:),                 &! in, optional, area fraction with wet land surface
              &          pcair = cair(:,:),                 &! in, optional, area fraction with wet land surface (air)
              &          paz0lh = z0h_lnd(:,:)          )    ! in, optional, roughness length for heat over land
        !
        !
        ! DA: vdiff_down has its own ACC WAIT due to automatic arrays in ACC data sections
        !----------------------------------------------------------------------------------------

        rls = grf_bdywidth_c+1
        rle = min_rlcell_int

        ncd = MAX(1,patch%n_childdom)
        jbs = patch%cells%start_blk(rls,  1)
        jbe = patch%cells%end_blk  (rle,ncd)

!$OMP PARALLEL DO PRIVATE(jb,jl,jk,jsfc,jcs,jce)
        !## jb loop2, is_in_sd_ed_interval == T .and. is_active == T
        DO jb = jbs, jbe
          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
          IF (jcs>jce) CYCLE
          !
          ! store in memory for output, if requested
          !
          ! DA: fuse all the 1D copies in a single ACC kernel
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          IF (ASSOCIATED(field% wstar)) THEN
            !$ACC LOOP GANG VECTOR
            DO jl=jcs,jce
              field% wstar(jl,jb) = wstar(jl,jb)
            END DO
          END IF
          IF (ASSOCIATED(field% hdtcbl)) THEN
            !$ACC LOOP GANG VECTOR
            DO jl=jcs,jce
              field% hdtcbl(jl,jb) = hdtcbl(jl,jb)
            END DO
          END IF
          !$ACC END PARALLEL
          !
          IF (ASSOCIATED(field% qs_sfc_tile)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jsfc=1,nsfc_type
              DO jl=jcs,jce
                field% qs_sfc_tile(jl,jb,jsfc) = qs_sfc_tile(jl,jb,jsfc)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% ri_atm)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk=1,nlev
              DO jl=jcs,jce
                field% ri_atm(jl,jk,jb) = ri_atm(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% mixlen)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk=1,nlev
              DO jl=jcs,jce
                field% mixlen(jl,jk,jb) = mixlen(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% cfm)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk=1,nlev
              DO jl=jcs,jce
                field% cfm(jl,jk,jb) = cfm(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% cfm_tile)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jsfc=1,nsfc_type
              DO jl=jcs,jce
                field% cfm_tile(jl,jb,jsfc) = cfm_tile(jl,jb,jsfc)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% cfh)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk=1,nlev
              DO jl=jcs,jce
                field% cfh(jl,jk,jb) = cfh(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% cfh_tile)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jsfc=1,nsfc_type
              DO jl=jcs,jce
                field% cfh_tile(jl,jb,jsfc) = cfh_tile(jl,jb,jsfc)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% cfv)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk=1,nlev
              DO jl=jcs,jce
                field% cfv(jl,jk,jb) = cfv(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% cftotte)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk=1,nlev
              DO jl=jcs,jce
                field% cftotte(jl,jk,jb) = cftotte(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(field% cfthv)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk=1,nlev
              DO jl=jcs,jce
                field% cfthv(jl,jk,jb) = cfthv(jl,jk,jb)
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END IF
          IF (ASSOCIATED(tend% wa_vdf)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk=1,nlev+1
              DO jl=jcs,jce
                tend% wa_vdf(jl,jk,jb) = tend_wa_vdf(jl,jk,jb)
              END DO
            END DO
          END IF
          !
        END DO !## jb loop 2
        !$ACC WAIT
!$OMP END PARALLEL DO
        !
        !
        IF (ltimer) CALL timer_stop(timer_vdf_dn)
        !
        rls = grf_bdywidth_c+1
        rle = min_rlcell_int

        ncd = MAX(1,patch%n_childdom)
        jbs = patch%cells%start_blk(rls,  1)
        jbe = patch%cells%end_blk  (rle,ncd)

!$OMP PARALLEL DO PRIVATE(jb,jl,jk,jt,jcs,jce)
        !## jb loop2-1, is_in_sd_ed_interval == T .and. is_active == T
        DO jb = jbs, jbe
          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
          IF (jcs>jce) CYCLE
          !
          ! Surface processes that provide time-dependent lower boundary
          ! condition for wind, temperature, tracer concentration, etc.
          !
          !$ACC WAIT
          !
          IF (ltimer) CALL timer_start(timer_vdf_sf)

          CALL update_surface(jg, jcs, jce, nproma, field%kice,               &! in
               &              nlev, nsfc_type,                                &! in
               &              iwtr, iice, ilnd,                               &! in, indices of surface types
               &              datetime,                                       &! in, date and time at beginning of this time step
               &              pdtime,                                         &! in, time step
               &              field%frac_tile(:,jb,:),                        &! in, area fraction
               &                     cfh_tile(:,jb,:),                        &! in, from "vdiff_down"
               &                     cfm_tile(:,jb,:),                        &! in, from "vdiff_down"
               &              zfactor_sfc(:,jb),                              &! in, from "vdiff_down"
               &              field% ocu (:,jb),                              &! in, ocean sfc velocity, u-component
               &              field% ocv (:,jb),                              &! in, ocean sfc velocity, v-component
               &              zaa(:,:,:,:,jb), zaa_btm(:,:,:,:,jb),           &! inout
               &              zbb(:,:,:,jb),   zbb_btm(:,:,:,jb),             &! inout
               &               zcpt_sfc_tile(:,jb,:),                         &! inout, from "vdiff_down", for "vdiff_up"
               &                 qs_sfc_tile(:,jb,:),                         &! inout, from "vdiff_down", for "vdiff_up"
               &                     ts_tile(:,jb,:),                         &! inout
               &              field%u_stress    (:,  jb),                     &! out
               &              field%v_stress    (:,  jb),                     &! out
               &              field% lhflx      (:,  jb),                     &! out
               &              field% shflx      (:,  jb),                     &! out
               &                     evap       (:,jb),                       &! out, for "cucall"
               &              field%u_stress_tile  (:,jb,:),                  &! out
               &              field%v_stress_tile  (:,jb,:),                  &! out
               &              field% lhflx_tile    (:,jb,:),                  &! out
               &              field% shflx_tile    (:,jb,:),                  &! out
               &              field%  evap_tile    (:,jb,:),                  &! out
               &              field%  fco2nat      (:,  jb),                  &! out
               &              nblock = jb,                                    &! in
               &              lsm = field%lsmask(:,jb),                       &! in, land-sea mask
               &              alake = field%alake(:,jb),                      &! in, lake fraction
               &              pu    = field% ua(:,nlev,jb),                   &! in, um1
               &              pv    = field% va(:,nlev,jb),                   &! in, vm1
               &              ptemp = field% ta(:,nlev,jb),                   &! in, tm1
               &              pq = field% qtrc_phy(:,nlev,jb,iqv),            &! in, qm1
               &              pco2 = zco2(:,jb),                              &! in, co2, lowest level or fixed value
               &              prsfl = field% rsfl(:,jb),                      &! in, rain surface large scale (from cloud)
               &              pssfl = field% ssfl(:,jb),                      &! in, snow surface large scale (from cloud)
               &              rlds        = field% rlds (:,jb),               &! in,  downward surface  longwave flux [W/m2]
               &              rlus        =        rlus (:,jb),               &! inout, upward surface  longwave flux [W/m2]
               &              rsds        = field% rsds (:,jb),               &! in,  downward surface shortwave flux [W/m2]
               &              rsus        = field% rsus (:,jb),               &! in,  upward surface shortwave flux [W/m2]
               !
               &              rvds_dir   = field%rvds_dir   (:,jb),           &! in, all-sky downward direct visible radiation at surface
               &              rpds_dir   = field%rpds_dir   (:,jb),           &! in, all-sky downward direct PAR     radiation at surface
               &              rnds_dir   = field%rnds_dir   (:,jb),           &! in, all-sky downward direct near-IR radiation at surface
               &              rvds_dif   = field%rvds_dif   (:,jb),           &! in, all-sky downward diffuse visible radiation at surface
               &              rpds_dif   = field%rpds_dif   (:,jb),           &! in, all-sky downward diffuse PAR     radiation at surface
               &              rnds_dif   = field%rnds_dif   (:,jb),           &! in, all-sky downward diffuse near-IR radiation at surface
               !
               &              pmair = field%mair(:,:,jb),                     &! in, moist air mass
               &              ps = field% phalf(:,nlevp1,jb),                 &! in, paphm1, half level pressure
               &              pcosmu0 = field% cosmu0(:,jb),                  &! in, amu0_x, cos of zenith angle
               &              pch_tile = zch_tile(:,jb,:),                    &! in, from "vdiff_down" for JSBACH
               &              pcsat = csat(:,jb),                             &! inout, area fraction with wet land surface
               &              pcair = cair(:,jb),                             &! inout, area fraction with wet land surface (air)
               &              q_snocpymlt    = q_snocpymlt(:,jb),             &! out, heating  by melting snow on the canopy [W/m2]
               &              z0m_tile = z0m_tile(:,jb,:),                    &! inout, roughness length for momentum over tiles
               &              z0h_lnd  = z0h_lnd (:,jb),                      &! out, roughness length for heat over land
               &              albvisdir      = albvisdir     (:,jb)  ,        &! out
               &              albnirdir      = albnirdir     (:,jb)  ,        &! out
               &              albvisdif      = albvisdif     (:,jb)  ,        &! out
               &              albnirdif      = albnirdif     (:,jb)  ,        &! out
               &              albvisdir_tile = albvisdir_tile(:,jb,:),        &! out
               &              albnirdir_tile = albnirdir_tile(:,jb,:),        &! out
               &              albvisdif_tile = albvisdif_tile(:,jb,:),        &! out
               &              albnirdif_tile = albnirdif_tile(:,jb,:),        &! out
               &              albedo         = albedo        (:,jb)  ,        &! out
               &              albedo_tile    = albedo_tile   (:,jb,:),        &! out
               &              emissivity     = field% emissivity    (:,jb)  , &! in
               &              pco2_flux_tile = field% co2_flux_tile(:,jb,:),  &! inout
               &              ptsfc     = ts    (:,jb),                       &! out
               &              ptsfc_rad = ts_rad(:,jb),                       &! out
               &              rlns_tile = field%lwflxsfc_tile(:,jb,:),        &! out (for coupling)
               &              rsns_tile = field%swflxsfc_tile(:,jb,:),        &! out (for coupling)
               &              lake_ice_frc = field%lake_ice_frc(:,jb),        &! out
               &              Tsurf = field% Tsurf(:,:,jb),                   &! inout, for sea ice
               &              T1    = field% T1   (:,:,jb),                   &! inout, for sea ice
               &              T2    = field% T2   (:,:,jb),                   &! inout, for sea ice
               &              hi    = field% hi   (:,:,jb),                   &! in   , for sea ice
               &              hs    = field% hs   (:,:,jb),                   &! inout, for sea ice
               &              conc  = field% conc (:,:,jb),                   &! in   , for sea ice
               &              Qtop  = field% Qtop (:,:,jb),                   &! out  , for sea ice
               &              Qbot  = field% Qbot (:,:,jb),                   &! out  , for sea ice
               &              albvisdir_ice = albvisdir_ice(:,:,jb),          &! inout ice albedos
               &              albnirdir_ice = albnirdir_ice(:,:,jb),          &! inout
               &              albvisdif_ice = albvisdif_ice(:,:,jb),          &! inout
               &              albnirdif_ice = albnirdif_ice(:,:,jb))           ! inout
          !
          !
          ! DA: update_surface has its own ACC WAIT due to automatic arrays in ACC data sections
          !----------------------------------------------------------------------------------------
          IF (ltimer) CALL timer_stop(timer_vdf_sf)
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(field% q_snocpymlt)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
            DO jl = jcs,jce
              field% q_snocpymlt(jl, jb)   = q_snocpymlt(jl,jb)
            END DO
          END IF

          !
          !
          !
          ! Turbulent mixing, part II:
          ! - Elimination for the lowest model level using boundary conditions
          !   provided by the surface model(s);
          ! - Back substitution to get solution of the tridiagonal system;
          ! - Compute tendencies and additional diagnostics.
          !
          IF (ltimer) CALL timer_start(timer_vdf_up)

          IF(ntracer >= iqt) THEN
            tend_qtrc_vdf_iqt => tend_qtrc_vdf(:,:,:,iqt:)
          ELSE
            tend_qtrc_vdf_iqt => tend_qtrc_vdf_dummy
          ENDIF
          !$ACC WAIT
          !
          CALL vdiff_up(jcs, jce, nproma, nlev, nlevm1,  &! in
               &        ntrac, nsfc_type,                &! in
               &        iwtr,                            &! in, indices of different sfc types
               &        pdtime,                          &! in, time steps
               &        field%frac_tile(:,jb,:),         &! in, area fraction of each sfc type
               &               cfm_tile(:,jb,:),         &! in
               &        zaa(:,:,:,:,jb),                 &! in, from "vdiff_down"
               &        field%cptgz(:,:,jb),             &! in, from "vdiff_down"
               &        field%   ua(:,:,jb),             &! in, um1
               &        field%   va(:,:,jb),             &! in, vm1
               &        field%   ta(:,:,jb),             &! in, tm1
               &        field% mair(:,:,jb),             &! in, moist     air mass [kg/m2]
               &        field% qtrc_phy(:,:,jb,iqv),     &! in, qm1
               &        field% qtrc_phy(:,:,jb,iqc),     &! in, xlm1
               &        field% qtrc_phy(:,:,jb,iqi),     &! in, xim1
               &        field% qtrc_phy(:,:,jb,iqt:),    &! in, xtm1
               &        field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
               &             ztottevn(:,:,jb),           &! in, TTE at intermediate time step
               &        aes_vdf_config(jg),              &! in
               &        zbb(:,:,:,jb),                   &! inout
               &        zthvvar(:,:,jb),                 &! in
               &        dummyx(:,:,jb),                  &! inout
               &        z0m_tile(:,jb,:),                &! inout
               &                 kedisp(:,jb),           &! out, vert. integr. diss. kin. energy [W/m2]
               &          tend_ua_vdf(:,:,jb),           &! out
               &          tend_va_vdf(:,:,jb),           &! out
               &                q_vdf(:,:,jb),           &! out   heating W/m2
               &        tend_qtrc_vdf(:,:,jb,iqv),       &! out
               &        tend_qtrc_vdf(:,:,jb,iqc),       &! out
               &        tend_qtrc_vdf(:,:,jb,iqi),       &! out
!               &        tend_qtrc_vdf(:,:,jb,iqt:),     &! out
               &        tend_qtrc_vdf_iqt(:,:,jb,:),     & ! out
               &        field%   z0m   (:,  jb),         &! out, for the next step
               &        dummy(:,:,jb),                   &!
               &        field%      totte(:,:,jb)        )! out
!!$               &        field%      totte(:,:,jb),       &! out
!!$               &        field%   sh_vdiff(:,  jb),       &! out, for energy diagnostic
!!$               &        field%   qv_vdiff(:,  jb)        )! out, for energy diagnostic
          !
          !
          ! DA: vdiff_up has its own ACC WAIT due to automatic arrays in ACC data sections in vdiff_tendencies
          !----------------------------------------------------------------------------------------

          IF ( turb == VDIFF_TURB_3DSMAGORINSKY ) THEN ! Smagorinksy
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs, jce
                tend_ua_vdf(jl,jk,jb) = tend_ua_vdf(jl,jk,jb) + ddt_u(jl,jk,jb)
                tend_va_vdf(jl,jk,jb) = tend_va_vdf(jl,jk,jb) + ddt_v(jl,jk,jb)
              END DO
            END DO
          END IF

          !
          ! Carbon -cycle
          ! -------------
          ! (see also before vdiff_down)
          ! For now do not change the co2 tracer if a prescribed co2 concentration
          ! is used as input to land (and ocean).
          IF (ccycle_config(jg)%iccycle == 2) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs,jce
                tend_qtrc_vdf(jl,jk,jb,ico2) = 0.0_wp
              END DO
            END DO
          END IF
          !
          !----------------------------------------------------------------------------------------
          !
          IF (ltimer) CALL timer_stop(timer_vdf_up)
          !
          ! store in memory for output or recycling
          !
          ! DA: fuse all the 1D copies in a single ACC kernel
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          IF (ASSOCIATED(field% kedisp  )) THEN
            !$ACC LOOP GANG VECTOR
            DO jl = jcs,jce
              field% kedisp  (jl, jb)   = kedisp (jl,jb)
            END DO
          END IF
          !
          IF (ASSOCIATED(field% q_vdf_vi)) THEN
            !$ACC LOOP GANG VECTOR
            DO jl = jcs,jce
              field% q_vdf_vi(jl,jb) = SUM(q_vdf(jl,:,jb))
            END DO
          END IF
          !$ACC END PARALLEL
          !
          IF (ASSOCIATED(field% q_vdf)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs,jce
                field% q_vdf(jl,jk,jb) = q_vdf(jl,jk,jb)
              END DO
            END DO
          END IF
          !
          IF (ASSOCIATED(tend% ua_vdf)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs,jce
                tend% ua_vdf(jl,jk,jb) = tend_ua_vdf(jl,jk,jb)
              END DO
            END DO
          END IF
          IF (ASSOCIATED(tend% va_vdf)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs,jce
                tend% va_vdf(jl,jk,jb) = tend_va_vdf(jl,jk,jb)
              END DO
            END DO
          END IF
          !
          IF (ASSOCIATED(tend% qtrc_vdf )) THEN
            IF (l2moment) THEN
              !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
              DO jk = 1,nlev
                DO jl = jcs,jce
                  tend% qtrc_vdf(jl,jk,jb,iqv) = tend_qtrc_vdf(jl,jk,jb,iqv)
                  tend% qtrc_vdf(jl,jk,jb,iqc) = tend_qtrc_vdf(jl,jk,jb,iqc)
                  tend% qtrc_vdf(jl,jk,jb,iqi) = tend_qtrc_vdf(jl,jk,jb,iqi)
                  tend% qtrc_vdf(jl,jk,jb,iqnc) = tend_qtrc_vdf(jl,jk,jb,iqnc)
                  tend% qtrc_vdf(jl,jk,jb,iqni) = tend_qtrc_vdf(jl,jk,jb,iqni)
                END DO
              END DO
            ELSE
              !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
              DO jk = 1,nlev
                DO jl = jcs,jce
                  tend% qtrc_vdf(jl,jk,jb,iqv) = tend_qtrc_vdf(jl,jk,jb,iqv)
                  tend% qtrc_vdf(jl,jk,jb,iqc) = tend_qtrc_vdf(jl,jk,jb,iqc)
                  tend% qtrc_vdf(jl,jk,jb,iqi) = tend_qtrc_vdf(jl,jk,jb,iqi)
                  !$ACC LOOP SEQ
                  DO jt = iqt,ntracer
                    !iqt is the index of the first non-water tracer
                    tend% qtrc_vdf(jl,jk,jb,jt) = tend_qtrc_vdf(jl,jk,jb,jt)
                  END DO
                END DO
              END DO
            END IF
          END IF
        !
        ENDDO !## jb loop2-1 END
        !$ACC WAIT
!$OMP END PARALLEL DO

      ELSE  ! if ( is_active )

        rls = grf_bdywidth_c+1
        rle = min_rlcell_int

        ncd = MAX(1,patch%n_childdom)
        jbs = patch%cells%start_blk(rls,  1)
        jbe = patch%cells%end_blk  (rle,ncd)

!$OMP PARALLEL DO PRIVATE(jb,jl,jk,jt,jcs,jce)
        !## jb loop3, is_in_sd_ed_interval == T .and. is_active == F
        DO jb = jbs, jbe
          CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
          IF (jcs>jce) CYCLE
          !
          ! retrieve from memory for recycling
          !
          IF (ASSOCIATED(field% q_snocpymlt)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
            DO jl = jcs,jce
              q_snocpymlt(jl,jb) = field% q_snocpymlt(jl,  jb)
            END DO
          END IF
          !
          IF (ASSOCIATED(field% q_vdf)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs,jce
                q_vdf(jl,jk,jb) = field% q_vdf(jl,jk,jb)
              END DO
            END DO
          END IF
          !
          IF (ASSOCIATED(tend% ua_vdf)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs,jce
                tend_ua_vdf(jl,jk,jb) = tend% ua_vdf(jl,jk,jb)
              END DO
            END DO
          END IF
          IF (ASSOCIATED(tend% va_vdf)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs,jce
                tend_va_vdf(jl,jk,jb) = tend% va_vdf(jl,jk,jb)
              END DO
            END DO
          END IF
          IF (ASSOCIATED(tend% wa_vdf)) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev+1
              DO jl = jcs,jce
                tend_wa_vdf(jl,jk,jb) = tend% wa_vdf(jl,jk,jb)
              END DO
            END DO
          END IF
          !
          IF (ASSOCIATED(tend% qtrc_vdf )) THEN
            IF (l2moment) THEN
              !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
              DO jk = 1,nlev
                DO jl = jcs,jce
                  tend_qtrc_vdf(jl,jk,jb,iqv) = tend% qtrc_vdf(jl,jk,jb,iqv)
                  tend_qtrc_vdf(jl,jk,jb,iqc) = tend% qtrc_vdf(jl,jk,jb,iqc)
                  tend_qtrc_vdf(jl,jk,jb,iqi) = tend% qtrc_vdf(jl,jk,jb,iqi)
                  tend_qtrc_vdf(jl,jk,jb,iqnc) = tend% qtrc_vdf(jl,jk,jb,iqnc)
                  tend_qtrc_vdf(jl,jk,jb,iqni) = tend% qtrc_vdf(jl,jk,jb,iqni)
                END DO
              END DO
            ELSE
              !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
              DO jk = 1,nlev
                DO jl = jcs,jce
                  tend_qtrc_vdf(jl,jk,jb,iqv) = tend% qtrc_vdf(jl,jk,jb,iqv)
                  tend_qtrc_vdf(jl,jk,jb,iqc) = tend% qtrc_vdf(jl,jk,jb,iqc)
                  tend_qtrc_vdf(jl,jk,jb,iqi) = tend% qtrc_vdf(jl,jk,jb,iqi)
                  !$ACC LOOP SEQ
                  DO jt = iqt,ntracer
                    tend_qtrc_vdf(jl,jk,jb,jt) = tend% qtrc_vdf(jl,jk,jb,jt)
                  END DO
                END DO
              END DO
            END IF
          END IF
        !
        ENDDO !## jb loop3 END
        !$ACC WAIT
!$OMP END PARALLEL DO

      END IF ! if ( is_active )
      !
      !
      rls = grf_bdywidth_c+1
      rle = min_rlcell_int

      ncd = MAX(1,patch%n_childdom)
      jbs = patch%cells%start_blk(rls,  1)
      jbe = patch%cells%end_blk  (rle,ncd)

!## jb loop4, is_in_sd_ed_interval == T
!$OMP PARALLEL DO PRIVATE(jb,jl,jk,jt,jsfc,jice,jcs,jce)
      DO jb = jbs, jbe
        CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
        IF (jcs>jce) CYCLE
        !
        ! Surface effect on the lowermost layer
        !
        IF (aes_phy_config(jg)%ljsb) THEN
          !
          ! convert    heating
          ! q_snocpymlt = heating for melting of snow on canopy
          !             = cooling of atmosphere --> negative sign
          !
          ! DA: fuse all the 1D copies in a single ACC kernel
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            tend_ta_sfc(jl,jb) = -q_snocpymlt(jl,jb) * field% qconv(jl,nlev,jb)
          END DO
          !
          IF (ASSOCIATED(tend% ta_sfc)) THEN
            !$ACC LOOP GANG VECTOR
            DO jl = jcs, jce
              tend% ta_sfc(jl,jb) = tend_ta_sfc(jl,jb)
            END DO
          END IF

          ! for output: accumulate heating
          IF (ASSOCIATED(field% q_phy)) THEN
            !$ACC LOOP GANG VECTOR
            DO jl = jcs, jce
              field% q_phy(jl,nlev,jb) = field% q_phy(jl,nlev,jb) - q_snocpymlt(jl,jb)
            END DO
          END IF
          IF (ASSOCIATED(field% q_phy_vi)) THEN
            !$ACC LOOP GANG VECTOR
            DO jl = jcs, jce
              field% q_phy_vi(jl, jb) = field% q_phy_vi(jl, jb) - q_snocpymlt(jl,jb)
            END DO
          END IF
          !
          ! accumulate tendencies for later updating the model state
          SELECT CASE(fc_vdf)
          CASE(0)
             ! diagnostic, do not use tendency
          CASE(1)
             ! use tendency to update the model state
             !$ACC LOOP GANG VECTOR
             DO jl = jcs, jce
               tend% ta_phy(jl,nlev,jb) = tend% ta_phy(jl,nlev,jb) + tend_ta_sfc(jl,jb)
             END DO
          END SELECT
          !
          ! update physics state for input to the next physics process
          SELECT CASE(fc_vdf)
          CASE(0)
             ! diagnostic, do not use tendency
          CASE(1)
             ! use tendency to update the physics state
             !$ACC LOOP GANG VECTOR
             DO jl = jcs, jce
               field% ta(jl,nlev,jb) = field% ta(jl,nlev,jb) + tend_ta_sfc(jl,jb)*pdtime
            END DO
          END SELECT
          !$ACC END PARALLEL

        END IF
        !
        !
        ! Vertical diffusion effect on the atmospheric column
        !
        ! convert    heating
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO jk = 1,nlev
          DO jl = jcs, jce
            tend_ta_vdf(jl,jk,jb) = q_vdf(jl,jk,jb) * field% qconv(jl,jk,jb)
          END DO
        END DO
        !
        ! for Smagorinsky
        IF ( turb == VDIFF_TURB_3DSMAGORINSKY ) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs, jce
              tend_ta_vdf(jl,jk,jb) = tend_ta_vdf(jl,jk,jb) + ta_hori_tend(jl,jk,jb)
              tend_qtrc_vdf(jl,jk,jb,iqv) = tend_qtrc_vdf(jl,jk,jb,iqv) + qv_hori_tend(jl,jk,jb)
              tend_qtrc_vdf(jl,jk,jb,iqc) = tend_qtrc_vdf(jl,jk,jb,iqc) + ql_hori_tend(jl,jk,jb)
              tend_qtrc_vdf(jl,jk,jb,iqi) = tend_qtrc_vdf(jl,jk,jb,iqi) + qi_hori_tend(jl,jk,jb)
            END DO
          END DO
          IF (l2moment) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs, jce
                tend_qtrc_vdf(jl,jk,jb,iqnc) = tend_qtrc_vdf(jl,jk,jb,iqnc) + qnc_hori_tend(jl,jk,jb)
                tend_qtrc_vdf(jl,jk,jb,iqni) = tend_qtrc_vdf(jl,jk,jb,iqni) + qni_hori_tend(jl,jk,jb)
              END DO
            END DO
          END IF ! l2moment
        END IF

        !
        IF (ASSOCIATED(tend% ta_vdf)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs, jce
              tend% ta_vdf(jl,jk,jb) = tend_ta_vdf(jl,jk,jb)
            END DO
          END DO
        END IF

        ! for output: accumulate heating
        IF (ASSOCIATED(field% q_phy)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs, jce
              field% q_phy(jl,jk,jb) = field% q_phy(jl,jk,jb) + q_vdf(jl,jk,jb)
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% q_phy_vi)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO jl = jcs, jce
            field% q_phy_vi(jl, jb) = field% q_phy_vi(jl, jb) + SUM(q_vdf(jl,:,jb))
          END DO
        END IF
        !
        ! accumulate tendencies for later updating the model state
        SELECT CASE(fc_vdf)
        CASE(0)
          ! diagnostic, do not use tendency
        CASE(1)
          ! use tendency to update the model state
          IF (l2moment) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs, jce
                tend%   ua_phy(jl,jk,jb)      = tend%   ua_phy(jl,jk,jb)      + tend_ua_vdf  (jl,jk,jb)
                tend%   va_phy(jl,jk,jb)      = tend%   va_phy(jl,jk,jb)      + tend_va_vdf  (jl,jk,jb)
                tend%   wa_phy(jl,jk,jb)      = tend%   wa_phy(jl,jk,jb)      + tend_wa_vdf  (jl,jk,jb)
                tend%   ta_phy(jl,jk,jb)      = tend%   ta_phy(jl,jk,jb)      + tend_ta_vdf  (jl,jk,jb)
                tend% qtrc_phy(jl,jk,jb,iqv)  = tend% qtrc_phy(jl,jk,jb,iqv)  + tend_qtrc_vdf(jl,jk,jb,iqv)
                tend% qtrc_phy(jl,jk,jb,iqc)  = tend% qtrc_phy(jl,jk,jb,iqc)  + tend_qtrc_vdf(jl,jk,jb,iqc)
                tend% qtrc_phy(jl,jk,jb,iqi)  = tend% qtrc_phy(jl,jk,jb,iqi)  + tend_qtrc_vdf(jl,jk,jb,iqi)
                tend% qtrc_phy(jl,jk,jb,iqnc) = tend% qtrc_phy(jl,jk,jb,iqnc) + tend_qtrc_vdf(jl,jk,jb,iqnc)
                tend% qtrc_phy(jl,jk,jb,iqni) = tend% qtrc_phy(jl,jk,jb,iqni) + tend_qtrc_vdf(jl,jk,jb,iqni)
              END DO
            END DO
          ELSE
            ! use tendency to update the model state
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs, jce
                tend%   ua_phy(jl,jk,jb)      = tend%   ua_phy(jl,jk,jb)      + tend_ua_vdf  (jl,jk,jb)
                tend%   va_phy(jl,jk,jb)      = tend%   va_phy(jl,jk,jb)      + tend_va_vdf  (jl,jk,jb)
                tend%   wa_phy(jl,jk,jb)      = tend%   wa_phy(jl,jk,jb)      + tend_wa_vdf  (jl,jk,jb)
                tend%   ta_phy(jl,jk,jb)      = tend%   ta_phy(jl,jk,jb)      + tend_ta_vdf  (jl,jk,jb)
                tend% qtrc_phy(jl,jk,jb,iqv)  = tend% qtrc_phy(jl,jk,jb,iqv)  + tend_qtrc_vdf(jl,jk,jb,iqv)
                tend% qtrc_phy(jl,jk,jb,iqc)  = tend% qtrc_phy(jl,jk,jb,iqc)  + tend_qtrc_vdf(jl,jk,jb,iqc)
                tend% qtrc_phy(jl,jk,jb,iqi)  = tend% qtrc_phy(jl,jk,jb,iqi)  + tend_qtrc_vdf(jl,jk,jb,iqi)
                !$ACC LOOP SEQ
                DO jt = iqt,ntracer
                  tend% qtrc_phy(jl,jk,jb,jt) = tend% qtrc_phy(jl,jk,jb,jt) + tend_qtrc_vdf(jl,jk,jb,jt)
                END DO
              END DO
            END DO
          END IF
          IF ( is_dry_cbl ) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs, jce
                tend% qtrc_phy(jl,jk,jb,iqv)  = 0._wp
                tend% qtrc_phy(jl,jk,jb,iqc)  = 0._wp
                tend% qtrc_phy(jl,jk,jb,iqi)  = 0._wp
              END DO
            END DO
          END IF

        END SELECT
        !
        ! update physics state for input to the next physics process
        SELECT CASE(fc_vdf)
        CASE(0)
          ! diagnostic, do not use tendency
        CASE(1)
          ! use tendency to update the physics state
          ! prognostic
          IF (l2moment) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs, jce
                field%       ua(jl,jk,jb)      = field%       ua(jl,jk,jb)      + tend_ua_vdf  (jl,jk,jb)     *pdtime
                field%       va(jl,jk,jb)      = field%       va(jl,jk,jb)      + tend_va_vdf  (jl,jk,jb)     *pdtime
                field%       wa(jl,jk,jb)      = field%       wa(jl,jk,jb)      + tend_wa_vdf  (jl,jk,jb)     *pdtime
                field%       ta(jl,jk,jb)      = field%       ta(jl,jk,jb)      + tend_ta_vdf  (jl,jk,jb)     *pdtime
                field% qtrc_phy(jl,jk,jb,iqv)  = field% qtrc_phy(jl,jk,jb,iqv)  + tend_qtrc_vdf(jl,jk,jb,iqv) *pdtime
                field% qtrc_phy(jl,jk,jb,iqc)  = field% qtrc_phy(jl,jk,jb,iqc)  + tend_qtrc_vdf(jl,jk,jb,iqc) *pdtime
                field% qtrc_phy(jl,jk,jb,iqi)  = field% qtrc_phy(jl,jk,jb,iqi)  + tend_qtrc_vdf(jl,jk,jb,iqi) *pdtime
                field% qtrc_phy(jl,jk,jb,iqnc) = field% qtrc_phy(jl,jk,jb,iqnc) + tend_qtrc_vdf(jl,jk,jb,iqnc)*pdtime
                field% qtrc_phy(jl,jk,jb,iqni) = field% qtrc_phy(jl,jk,jb,iqni) + tend_qtrc_vdf(jl,jk,jb,iqni)*pdtime
              END DO
            END DO
          ELSE
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs, jce
                field%       ua(jl,jk,jb)      = field%       ua(jl,jk,jb)      + tend_ua_vdf  (jl,jk,jb)     *pdtime
                field%       va(jl,jk,jb)      = field%       va(jl,jk,jb)      + tend_va_vdf  (jl,jk,jb)     *pdtime
                field%       wa(jl,jk,jb)      = field%       wa(jl,jk,jb)      + tend_wa_vdf  (jl,jk,jb)     *pdtime
                field%       ta(jl,jk,jb)      = field%       ta(jl,jk,jb)      + tend_ta_vdf  (jl,jk,jb)     *pdtime
                field% qtrc_phy(jl,jk,jb,iqv)  = field% qtrc_phy(jl,jk,jb,iqv)  + tend_qtrc_vdf(jl,jk,jb,iqv) *pdtime
                field% qtrc_phy(jl,jk,jb,iqc)  = field% qtrc_phy(jl,jk,jb,iqc)  + tend_qtrc_vdf(jl,jk,jb,iqc) *pdtime
                field% qtrc_phy(jl,jk,jb,iqi)  = field% qtrc_phy(jl,jk,jb,iqi)  + tend_qtrc_vdf(jl,jk,jb,iqi) *pdtime
                !$ACC LOOP SEQ
                DO jt = iqt,ntracer
                  field% qtrc_phy(jl,jk,jb,jt) = field% qtrc_phy(jl,jk,jb,jt)   + tend_qtrc_vdf(jl,jk,jb,jt)  *pdtime
                END DO
              END DO
            END DO
          END IF
          IF ( is_dry_cbl ) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
            DO jk = 1,nlev
              DO jl = jcs, jce
                tend% qtrc_phy(jl,jk,jb,iqv)  = 0._wp
                tend% qtrc_phy(jl,jk,jb,iqc)  = 0._wp
                tend% qtrc_phy(jl,jk,jb,iqi)  = 0._wp
              END DO
            END DO
          END IF 
          !
          ! diagnostic
          ! 2-tl-scheme
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs, jce
              field% tottem1   (jl,jk,jb) = field% totte (jl,jk,jb)
            END DO
          END DO
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO jl = jcs, jce
            field% ts        (jl,jb)   = ts        (jl,jb)
            field% ts_rad    (jl,jb)   = ts_rad    (jl,jb)
            field% ustar     (jl,jb)   = ustar     (jl,jb)
            field% evap      (jl,jb)   = evap      (jl,jb)
            field% cair      (jl,jb)   = cair      (jl,jb)
            field% csat      (jl,jb)   = csat      (jl,jb)
            field% z0h_lnd   (jl,jb)   = z0h_lnd   (jl,jb)
            field% rlus      (jl,jb)   = rlus      (jl,jb)
            !
            field% albvisdir (jl,jb)   = albvisdir (jl,jb)
            field% albnirdir (jl,jb)   = albnirdir (jl,jb)
            field% albvisdif (jl,jb)   = albvisdif (jl,jb) 
            field% albnirdif (jl,jb)   = albnirdif (jl,jb)
            field% albedo    (jl,jb)   = albedo    (jl,jb)
          END DO
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jsfc = 1,nsfc_type
            DO jl = jcs,jce
              field% ts_tile   (jl,jb,jsfc) = ts_tile   (jl,jb,jsfc)
              field% wstar_tile(jl,jb,jsfc) = wstar_tile(jl,jb,jsfc)
              field% z0m_tile  (jl,jb,jsfc) = z0m_tile  (jl,jb,jsfc)
              !
              field% albvisdir_tile (jl,jb,jsfc) = albvisdir_tile(jl,jb,jsfc)
              field% albnirdir_tile (jl,jb,jsfc) = albnirdir_tile(jl,jb,jsfc)
              field% albvisdif_tile (jl,jb,jsfc) = albvisdif_tile(jl,jb,jsfc)
              field% albnirdif_tile (jl,jb,jsfc) = albnirdif_tile(jl,jb,jsfc)
              field% albedo_tile    (jl,jb,jsfc) = albedo_tile   (jl,jb,jsfc)
              !
            END DO
          END DO
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jice=1,nice
            DO jl=jcs,jce
              field% albvisdir_ice(jl,jice,jb) = albvisdir_ice(jl,jice,jb)
              field% albnirdir_ice(jl,jice,jb) = albnirdir_ice(jl,jice,jb)
              field% albvisdif_ice(jl,jice,jb) = albvisdif_ice(jl,jice,jb)
              field% albnirdif_ice(jl,jice,jb) = albnirdif_ice(jl,jice,jb)
            END DO
          END DO
          !
        END SELECT
        !
        !
        ! Correction related to implicitness, due to the fact that surface model only used
        ! part of longwave radiation to compute new surface temperature
        !
        ! DA: fuse all the 1D copies in a single ACC kernel
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        IF ( isrfc_type >= 5 ) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            q_rlw_impl(jl,jb) = 0._wp
          END DO
        ELSE
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            q_rlw_impl(jl,jb) =                                               &
                &  ( (field%rld_rt(jl,nlev,jb)-field%rlu_rt(jl,nlev,jb))  & ! ( rln  from "radiation", at top of layer nlev
                &   -(field%rlds  (jl,jb)     -field%rlus  (jl,jb)     )) & !  -rlns from "radheating" and "update_surface")
                & -field%q_rlw_nlev(jl,jb)                                       ! -old heating in layer nlev from "radheating"
          END DO
        END IF
        !
        IF (ASSOCIATED(field%q_rlw_impl)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field%q_rlw_impl(jl,jb) = q_rlw_impl(jl,jb)
          END DO
        END IF

        ! convert    heating
        !$ACC LOOP GANG VECTOR
        DO jl = jcs,jce
          tend_ta_rlw_impl(jl,jb) = q_rlw_impl(jl,jb) * field% qconv(jl,nlev,jb)
        END DO
        !
        IF (ASSOCIATED(tend%ta_rlw_impl)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            tend%ta_rlw_impl(jl,jb) = tend_ta_rlw_impl(jl,jb)
          END DO
        END IF

        ! for output: accumulate heating
        IF (ASSOCIATED(field% q_phy)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% q_phy(jl,nlev,jb) = field% q_phy(jl,nlev,jb) + q_rlw_impl(jl,jb)
          END DO
        END IF
        IF (ASSOCIATED(field% q_phy_vi)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% q_phy_vi(jl,jb) = field% q_phy_vi(jl,jb) + q_rlw_impl(jl,jb)
          END DO
        END IF
        !
        ! accumulate surface lw increment for later updating the model state
        ! but only if radiative processes are active.
        IF ( ldtrad_gt0 ) THEN
          SELECT CASE(fc_vdf)
          CASE(1)
            ! use tendency to update the model state
            !$ACC LOOP GANG VECTOR
            DO jl = jcs,jce
              tend%ta_phy(jl,nlev,jb) = tend%ta_phy(jl,nlev,jb) + tend_ta_rlw_impl(jl,jb)
            END DO
          END SELECT
          !
          ! update physics state for input to the next physics process
          SELECT CASE(fc_vdf)
          CASE(0)
            ! diagnostic, do not use tendency
          CASE(1)
            ! use tendency to update the physics state
            !$ACC LOOP GANG VECTOR
            DO jl = jcs,jce
              field% ta(jl,nlev,jb) = field% ta(jl,nlev,jb) + tend_ta_rlw_impl(jl,jb)*pdtime
            END DO
          END SELECT
        END IF
        !$ACC END PARALLEL
        !$ACC WAIT

        ! Turbulent mixing, part III:
        ! - Further diagnostics.
        !----------------------------------------------------------------------------------------
        ! TODO: async
        CALL nsurf_diag(jcs, jce, nproma, nsfc_type,     &! in
            &          ilnd,                            &! in
            &          field%frac_tile(:,jb,:),         &! in
            &          field% qtrc_phy(:,nlev,jb,iqv),  &! in humidity qm1
            &          field%    ta(:,nlev,jb),         &! in tm1
            &          field% pfull(:,nlev,jb),         &! in, apm1
            &          field% phalf(:,nlevp1,jb),       &! in, aphm1
            &                  zqx(:,nlev,jb),          &! in, xlm1 + xim1
            &          field%   ua(:,nlev,jb),          &! in, um1
            &          field%   va(:,nlev,jb),          &! in, vm1
            &          field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
            &          field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
            &          field% zf  (:,nlev  ,jb),        &! in, height of lowermost full level (m)
            &          field% zh  (:,nlev+1,jb),        &! in, surface height    (m)
            &          field%cptgz(:,nlev,jb),          &! in dry static energy
            &          zcpt_sfc_tile(:,jb,:),           &! in dry static energy
            &          zbn_tile(:,jb,:),                &! in for diagnostic
            &          zbhn_tile(:,jb,:),               &! in for diagnostic
            &          zbh_tile(:,jb,:),                &! in for diagnostic
            &          zbm_tile(:,jb,:),                &! in for diagnostic
            &          ri_tile(:,jb,:),                 &! in
            &          field%sfcWind(:,  jb),           &! out 10m windspeed
            &          field%    tas(:,  jb),           &! out temperature in 2m
            &          field%   dew2(:,  jb),           &! out dew point temperature in 2m
            &          field%    uas(:,  jb),           &! out zonal wind in 10m
            &          field%    vas(:,  jb),           &! out meridional wind in 10m
            &          field%tasmax (:,  jb),           &! out max 2m temperature
            &          field%tasmin (:,  jb),           &! out min 2m temperature
            &          field%sfcWind_tile(:,jb,:),      &! out 10m windspeed on tiles
            &          field%    tas_tile(:,jb,:),      &! out temperature in 2m on tiles
            &          field%   dew2_tile(:,jb,:),      &! out dew point temperature in 2m on tiles
            &          field%    uas_tile(:,jb,:),      &! out zonal wind in 10m on tiles
            &          field%    vas_tile(:,jb,:)       )! out meridional wind in 10m on tiles

      ENDDO !## jb loop4 END
!$OMP END PARALLEL DO
      !
      !
    ELSE ! if ( is_in_sd_ed_interval )
      !
      rls = grf_bdywidth_c+1
      rle = min_rlcell_int

      ncd = MAX(1,patch%n_childdom)
      jbs = patch%cells%start_blk(rls,  1)
      jbe = patch%cells%end_blk  (rle,ncd)

!$OMP PARALLEL DO PRIVATE(jb,jl,jk,jt,jsfc,jice,jcs,jce)
      !## jb loop5, is_in_sd_ed_interval == F
      DO jb = jbs, jbe
        CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
        IF (jcs>jce) CYCLE
        !
        !
        ! vdiff_down
        ! DA: fuse all the 1D copies in a single ACC kernel
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jl = jcs,jce
          field% ustar(jl,jb) = 0.0_wp
        END DO
        IF (ASSOCIATED(field% wstar)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% wstar(jl,jb) = 0.0_wp
          END DO
        END IF
        IF (ASSOCIATED(field% hdtcbl  )) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% hdtcbl(jl, jb) = 0.0_wp
          END DO
        END IF
        !$ACC END PARALLEL
        !
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO jsfc = 1,nsfc_type
          DO jl = jcs,jce
            field% wstar_tile(jl,jb,jsfc) = 0.0_wp
          END DO
        END DO
        IF (ASSOCIATED(field% qs_sfc_tile)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jsfc = 1,nsfc_type
            DO jl = jcs,jce
              field% qs_sfc_tile (jl,jb,jsfc) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% ri_atm  )) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              field% ri_atm(jk,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% mixlen)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              field% mixlen(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% cfm)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              field% cfm(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% cfm_tile)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jsfc = 1,nsfc_type
            DO jl = jcs,jce
              field% cfm_tile(jl,jb,jsfc) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% cfh)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              field% cfh(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% cfh_tile)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jsfc = 1,nsfc_type
            DO jl = jcs,jce
              field% cfh_tile(jl,jb,jsfc) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% cfv)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              field% cfv(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% cftotte )) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              field% cftotte(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(field% cfthv)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              field% cfthv(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        !
        ! update_surface
        ! DA: fuse all the 1D copies in a single ACC kernel
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jl = jcs,jce
          field% u_stress  (jl,jb) = 0.0_wp
          field% v_stress  (jl,jb) = 0.0_wp
          field% lhflx     (jl,jb) = 0.0_wp
          field% shflx     (jl,jb) = 0.0_wp
          field%  evap     (jl,jb) = 0.0_wp
          field% rlus      (jl,jb) = 0.0_wp
          field% rsus      (jl,jb) = 0.0_wp
          field% csat      (jl,jb) = 0.0_wp
          field% cair      (jl,jb) = 0.0_wp
          field% z0h_lnd   (jl,jb) = 0.0_wp
          field% albvisdir (jl,jb) = 0.0_wp
          field% albnirdir (jl,jb) = 0.0_wp
          field% albvisdif (jl,jb) = 0.0_wp
          field% albnirdif (jl,jb) = 0.0_wp
          field% albedo    (jl,jb) = 0.0_wp
          field% ts        (jl,jb) = 0.0_wp
          field% ts_rad    (jl,jb) = 0.0_wp
        END DO
        IF (ASSOCIATED(field% q_snocpymlt)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% q_snocpymlt (jl,jb) = 0.0_wp
          END DO
        END IF
        IF (ASSOCIATED(field% q_rlw_impl)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% q_rlw_impl(jl,jb) = 0.0_wp
          END DO
        END IF
        IF (ASSOCIATED(tend% ta_sfc     )) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            tend% ta_sfc(jl,jb) = 0.0_wp
          END DO
        END IF
        IF (ASSOCIATED(tend% ta_rlw_impl)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            tend% ta_rlw_impl(jl,jb) = 0.0_wp
          END DO
        END IF
        !$ACC END PARALLEL
        !
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO jsfc = 1,nsfc_type
          DO jl = jcs,jce
            field% ts_tile        (jl,jb,jsfc) = 0.0_wp
            field% u_stress_tile  (jl,jb,jsfc) = 0.0_wp
            field% v_stress_tile  (jl,jb,jsfc) = 0.0_wp
            field% lhflx_tile     (jl,jb,jsfc) = 0.0_wp
            field% shflx_tile     (jl,jb,jsfc) = 0.0_wp
            field%  evap_tile     (jl,jb,jsfc) = 0.0_wp
            field% lwflxsfc_tile  (jl,jb,jsfc) = 0.0_wp
            field% swflxsfc_tile  (jl,jb,jsfc) = 0.0_wp
            field% albvisdir_tile (jl,jb,jsfc) = 0.0_wp
            field% albnirdir_tile (jl,jb,jsfc) = 0.0_wp
            field% albvisdif_tile (jl,jb,jsfc) = 0.0_wp
            field% albnirdif_tile (jl,jb,jsfc) = 0.0_wp
            field% albedo_tile    (jl,jb,jsfc) = 0.0_wp
            field% co2_flux_tile  (jl,jb,jsfc) = 0.0_wp
          END DO
        END DO
        !
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO jice=1,nice
          DO jl = jcs,jce
            field% Tsurf          (jl,jice,jb) = 0.0_wp
            field% T1             (jl,jice,jb) = 0.0_wp
            field% T2             (jl,jice,jb) = 0.0_wp
            field% Qtop           (jl,jice,jb) = 0.0_wp
            field% Qbot           (jl,jice,jb) = 0.0_wp
            field% albvisdir_ice  (jl,jice,jb) = 0.0_wp
            field% albnirdir_ice  (jl,jice,jb) = 0.0_wp
            field% albvisdif_ice  (jl,jice,jb) = 0.0_wp
            field% albnirdif_ice  (jl,jice,jb) = 0.0_wp
          END DO
        END DO
        !
        ! vdiff_up
        ! DA: fuse all the 1D copies in a single ACC kernel
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jl = jcs,jce
          field% z0m(jl,jb) = 0.0_wp
        END DO
        IF (ASSOCIATED(field% kedisp)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% kedisp(jl,jb) = 0.0_wp
          END DO
        END IF
        IF (ASSOCIATED(field% q_vdf_vi)) THEN
          !$ACC LOOP GANG VECTOR
          DO jl = jcs,jce
            field% q_vdf_vi(jl,jb) = 0.0_wp
          END DO
        END IF
        !$ACC END PARALLEL
        !
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO jk = 1,nlev
          DO jl = jcs,jce
            field% totte(jl,jk,jb) = 0.0_wp
          END DO
        END DO
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO jsfc = 1,nsfc_type
          DO jl = jcs,jce
            field% z0m_tile(jl,jb,jsfc) = 0.0_wp
          END DO
        END DO
        !
        IF (ASSOCIATED(field% q_vdf)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              field% q_vdf(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        !
        IF (ASSOCIATED(tend% ta_vdf)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              tend% ta_vdf(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(tend% ua_vdf)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              tend% ua_vdf(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(tend% va_vdf)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              tend% va_vdf(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        IF (ASSOCIATED(tend% wa_vdf)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev+1
            DO jl = jcs,jce
              tend% wa_vdf(jl,jk,jb) = 0.0_wp
            END DO
          END DO
        END IF
        !
        IF (ASSOCIATED(tend% qtrc_vdf)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 1,nlev
            DO jl = jcs,jce
              tend% qtrc_vdf(jl,jk,jb,iqv)  = 0.0_wp
              tend% qtrc_vdf(jl,jk,jb,iqc)  = 0.0_wp
              tend% qtrc_vdf(jl,jk,jb,iqi)  = 0.0_wp
              !$ACC LOOP SEQ
              DO jt = iqt,ntracer
                tend% qtrc_vdf(jl,jk,jb,jt) = 0.0_wp
              END DO
            END DO
          END DO
        END IF
        !$ACC WAIT
      !
      ENDDO !## jb loop5 END
!$OMP END PARALLEL DO
    !
    END IF ! if ( is_in_sd_ed_interval )

    !$ACC WAIT
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA

    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) CALL timer_stop(timer_vdf)

  END SUBROUTINE interface_aes_vdf

END MODULE mo_interface_aes_vdf
