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

MODULE mo_surface

  USE mo_kind,              ONLY: wp, i1






!#ifdef __NO_ICON_OCEAN__
  USE mo_exception,         ONLY: finish
!#endif

  USE mo_physical_constants,ONLY: grav, Tf, alf, albedoW, stbo, tmelt, rhos!!$, rhoi
  USE mo_physical_constants,ONLY: cvd, cpd
  USE mo_coupling_config,   ONLY: is_coupled_run
  USE mo_aes_phy_config,    ONLY: aes_phy_config
  USE mo_aes_phy_memory,    ONLY: cdimissval
  USE mo_aes_vdf_config,    ONLY: aes_vdf_config
  USE mo_turb_vdiff,        ONLY: ih, iqv, iu, iv, imh, imqv, imuv, &
                                & nmatrix, nvar_vdiff,              &
                                & matrix_to_richtmyer_coeff
  USE mo_turb_vdiff_params, ONLY: tpfac2
  USE mo_surface_diag,      ONLY: wind_stress, surface_fluxes
  USE mo_index_list,        ONLY: generate_index_list_batched
  USE mtime,                ONLY: t_datetime => datetime

  USE mo_jsb_interface,     ONLY: jsbach_interface

  USE mo_aes_sfc_indices,   ONLY: nsfc_type

  USE mo_ice_interface,     ONLY: ice_fast
  USE mo_ml_ocean,          ONLY: ml_ocean






  USE mo_physical_constants,ONLY: cpd
  USE mo_nh_testcases_nml,  ONLY: isrfc_type, shflx, lhflx
  USE mo_physical_constants,ONLY: alv

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_surface

  ! Shortcuts to components of aes_vdf_config
  !
  LOGICAL :: lsfc_mom_flux, lsfc_heat_flux

CONTAINS
  !>
  !!
  !!
  SUBROUTINE update_surface( jg,                                &! in
                           & jcs, kproma, kbdim,                &! in
                           & kice,                              &! in
                           & klev, ksfc_type,                   &! in
                           & idx_wtr, idx_ice, idx_lnd,         &! in
                           & datetime,                          &! in
                           & pdtime,                            &! in
                           & pfrc,                              &! in
                           & pcfh_tile, pcfm_tile,              &! in
                           & pfac_sfc, pocu, pocv,              &! in
                           & aa, aa_btm, bb, bb_btm,            &! inout
                           & pcpt_tile, pqsat_tile,             &! inout
                           & ptsfc_tile,                        &! inout
                           & pu_stress_gbm, pv_stress_gbm,      &! out
                           & plhflx_gbm, pshflx_gbm,            &! out
                           & pevap_gbm,                         &! out
                           & pu_stress_tile,   pv_stress_tile,  &! out
                           & plhflx_tile, pshflx_tile,          &! out
                           & pevap_tile,                        &! out
                           & pco2nat,                           &! out
                           !! optional
                           & nblock,                            &! in
                           & lsm,                               &! in
                           & alake,                             &! in
                           & pu,                                &! in
                           & pv,                                &! in
                           & ptemp,                             &! in
                           & pq,                                &! in
                           & pco2,                              &! in
                           & prsfl,                             &! in
                           & pssfl,                             &! in
                           & rlds,                              &! in
                           & rlus,                              &! inout
                           & rsds,                              &! in
                           & rsus,                              &! in
                           !
                           & rvds_dir,                          &! in
                           & rpds_dir,                          &! in
                           & rnds_dir,                          &! in
                           & rvds_dif,                          &! in
                           & rpds_dif,                          &! in
                           & rnds_dif,                          &! in
                           !
                           & pmair,                             &! in
                           & ps,                                &! in
                           & pcosmu0,                           &! in
                           & pch_tile,                          &! in
                           !! for JSBACH
                           & pcsat,                             &! inout
                           & pcair,                             &! inout
                           & q_snocpymlt,                       &! out
                           !
                           & z0m_tile,                          &! inout
                           & z0h_lnd,                           &! out
                           & albvisdir, albnirdir, albvisdif, albnirdif, &! out
                           & albvisdir_tile,                    &! out
                           & albnirdir_tile,                    &! out
                           & albvisdif_tile,                    &! out
                           & albnirdif_tile,                    &! out
                           & albedo, albedo_tile,               &! out
                           & emissivity,                        &! in
                           & pco2_flux_tile,                    &! inout
                           & ptsfc,                             &! out
                           & ptsfc_rad,                         &! out
                           & rsns_tile, rlns_tile,              &! out
                           & lake_ice_frc,                      &! out
                           !! Sea ice
                           & Tsurf,                             &! inout
                           & T1,                                &! inout
                           & T2,                                &! inout
                           & hi,                                &! in
                           & hs,                                &! inout
                           & Qtop,                              &! out
                           & Qbot,                              &! out
                           & conc,                              &! in
                           & albvisdir_ice, albvisdif_ice,      &! inout
                           & albnirdir_ice, albnirdif_ice)       ! inout

    TYPE(t_datetime), INTENT(IN), POINTER :: datetime ! date and time at the end of this time step
    REAL(wp),INTENT(IN) :: pdtime
    INTEGER, INTENT(IN) :: jg
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim
    INTEGER, INTENT(IN) :: klev, ksfc_type
    INTEGER, INTENT(IN) :: idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN) :: pfrc      (:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcfh_tile (:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcfm_tile (:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pfac_sfc  (:)   ! (kbdim)
    REAL(wp),INTENT(IN) :: pocu      (:)   ! (kbdim)
    REAL(wp),INTENT(IN) :: pocv      (:)   ! (kbdim)
    REAL(wp),INTENT(INOUT) :: aa     (:,:,:,:)    ! (kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT) :: aa_btm (:,:,:,imh:) ! (kbdim,3,ksfc_type,imh:imqv)
    REAL(wp),INTENT(INOUT) :: bb     (:,:,:)   ! (kbdim,klev,nvar_vdiff)
    REAL(wp),INTENT(INOUT) :: bb_btm (:,:,ih:) ! (kbdim,ksfc_type,ih:iqv)
    REAL(wp),INTENT(INOUT) :: pcpt_tile (:,:)  ! (kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: pqsat_tile(:,:)  ! (kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: ptsfc_tile (:,:) ! (kbdim,ksfc_type)

    REAL(wp),INTENT(OUT)   :: pu_stress_gbm (:) ! (kbdim)
    REAL(wp),INTENT(OUT)   :: pv_stress_gbm (:) ! (kbdim)
    REAL(wp),INTENT(OUT)   ::    plhflx_gbm (:) ! (kbdim)
    REAL(wp),INTENT(OUT)   ::    pshflx_gbm (:) ! (kbdim)
    REAL(wp),INTENT(OUT)   ::     pevap_gbm (:) ! (kbdim)

    REAL(wp),INTENT(OUT)   :: pu_stress_tile (:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pv_stress_tile (:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: plhflx_tile (:,:)    ! (kbdim,ksfc_type) see surface_fluxes
    REAL(wp),INTENT(OUT)   :: pshflx_tile (:,:)    ! (kbdim,ksfc_type) see surface_fluxes
    REAL(wp),INTENT(OUT)   :: pevap_tile (:,:)     ! (kbdim,ksfc_type) see surface_fluxes

    !! JSBACH input
    INTEGER, OPTIONAL,INTENT(IN) :: nblock
    REAL(wp),OPTIONAL,INTENT(IN) :: lsm(:)          ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(IN) :: alake(:)        ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(IN) :: pu        (:)   ! (kbdim) zonal wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: pv        (:)   ! (kbdim) meridional wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: ptemp     (:)   ! (kbdim) temperature of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: pq        (:)   ! (kbdim) humidity of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: pco2      (:)   ! (kbdim) co2 of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfl     (:)   ! (kbdim) rain large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfl     (:)   ! (kbdim) snow large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: rlds      (:)   ! (kbdim) downward surface  longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: rsds      (:)   ! (kbdim) downward surface shortwave flux [W/m2]
    
    REAL(wp),INTENT(IN) :: rvds_dir(:)        ! (kbdim) all-sky   vis. dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rpds_dir(:)        ! (kbdim) all-sky   par  dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rnds_dir(:)        ! (kbdim) all-sky   nir  dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rvds_dif(:)        ! (kbdim) all-sky   vis. dif. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rpds_dif(:)        ! (kbdim) all-sky   par  dif. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rnds_dif(:)        ! (kbdim) all-sky   nir  dif. downward flux at current   time [W/m2]

    REAL(wp),INTENT(IN) :: pmair(:,:)                   ! (kbdim,klev) air mass [kg/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: ps        (:)       ! (kbdim) surface pressure
    REAL(wp),OPTIONAL,INTENT(IN) :: pcosmu0   (:)       ! (kbdim) cos of zenith angle
    REAL(wp),OPTIONAL,INTENT(IN) :: pch_tile  (:,:)     ! (kbdim,ksfc_type)
    !! JSBACH output
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcsat(:)       ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcair(:)       ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: q_snocpymlt(:) ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: z0h_lnd(:)     ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: z0m_tile(:,:)  ! (kbdim,ksfc_type)
    !
    REAL(wp),OPTIONAL,INTENT(OUT)   :: albvisdir_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: albnirdir_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: albvisdif_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: albnirdif_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: albedo(:)           ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: albvisdir(:), albvisdif(:) ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: albnirdir(:), albnirdif(:) ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: albedo_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(IN)    :: emissivity(:) ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: pco2nat(:)    ! (kbdim)

    REAL(wp),OPTIONAL,INTENT(INOUT) :: pco2_flux_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: ptsfc(:)        ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: ptsfc_rad(:)    ! (kbdim)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rlus(:)         ! (kbdim)  INOUT upward surface  longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN)    :: rsus(:)         ! (kbdim) IN upward surface shortwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(OUT)   :: rsns_tile(:,:)  ! (kbdim,ksfc_type) shortwave net flux at surface on tiles
    REAL(wp),OPTIONAL,INTENT(OUT)   :: rlns_tile(:,:)  ! (kbdim,ksfc_type) longwave net flux at surface on tiles
    REAL(wp),OPTIONAL,INTENT(OUT)   :: lake_ice_frc(:) ! (kbdim) fraction of ice on lakes
    !! Sea ice
    INTEGER,          INTENT(IN)    :: kice ! Number of ice thickness classes
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Tsurf(:,:) ! (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T1   (:,:) ! (kbdim,kice) for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T2   (:,:) ! (kbdim,kice) for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(IN)    :: hi   (:,:) ! (kbdim,kice) for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hs   (:,:) ! (kbdim,kice) for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(OUT)   :: Qtop (:,:) ! (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: Qbot (:,:) ! (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(IN)    :: conc (:,:) ! (kbdim,kice) for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_ice(:,:) ! (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_ice(:,:) ! (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_ice(:,:) ! (kbdim,kice)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_ice(:,:) ! (kbdim,kice)

! locals

    INTEGER(i1) :: pfrc_test(kbdim,ksfc_type) !< integer mask to pass to CUB (can be removed later)
    INTEGER     :: loidx    (kbdim,ksfc_type) !< counter for masks
    INTEGER     :: is       (      ksfc_type) !< counter for masks

    INTEGER  :: jsfc, jk, jkm1, im, k, jl, jls, js

    REAL(wp) :: se_sum(kbdim), qv_sum(kbdim), wgt_sum(kbdim), wgt(kbdim)
    REAL(wp) :: zca(kbdim,ksfc_type), zcs(kbdim,ksfc_type)
    REAL(wp) :: zfrc_oce(kbdim)

    REAL(wp) :: zen_h (kbdim,ksfc_type)
    REAL(wp) :: zfn_h (kbdim,ksfc_type)
    REAL(wp) :: zen_qv(kbdim,ksfc_type)
    REAL(wp) :: zfn_qv(kbdim,ksfc_type)

    REAL(wp) ::                                                    &
      & zlhflx_lnd(kbdim), zlhflx_lwtr(kbdim), zlhflx_lice(kbdim), &
      & zshflx_lnd(kbdim), zshflx_lwtr(kbdim), zshflx_lice(kbdim), &
      & zevap_lnd(kbdim), zevap_lwtr(kbdim), zevap_lice(kbdim),    &
      & qsat_lnd(kbdim), qsat_lwtr(kbdim), qsat_lice(kbdim),       &
      & zcpt_lnd(kbdim), zcpt_lwtr(kbdim), zcpt_lice(kbdim),      &
      & ztsfc_lnd(kbdim), ztsfc_lnd_eff(kbdim),                    &
      & ztsfc_wtr(kbdim), ztsfc_lwtr(kbdim), ztsfc_lice(kbdim),    &
      & rvds(kbdim), rnds(kbdim), rpds(kbdim),                     &
      & rsns(kbdim), rlns(kbdim), fract_par_diffuse(kbdim),        &
      & zalbvis, zalbnir,                                          &
      & zalbedo_lwtr(kbdim), zalbedo_lice(kbdim),                  &
      & zwindspeed_lnd(kbdim), zwindspeed10m_lnd(kbdim)

    REAL(wp) :: rain_tmp(kbdim), snow_tmp(kbdim), drag_srf_tmp(kbdim), &
      & pch_tmp(kbdim), drag_wtr_tmp(kbdim), drag_ice_tmp(kbdim)

    REAL(wp) :: zgrnd_hflx(kbdim,ksfc_type), zgrnd_hcap(kbdim,ksfc_type)

    !REAL(wp) :: zt2s_conv(kbdim,ksfc_type)

    ! Sea ice
    REAL(wp) :: Tfw(kbdim)
    REAL(wp) :: swflx_ice(kbdim,kice), nonsolar_ice(kbdim,kice), dnonsolardT(kbdim,kice), conc_sum(kbdim)

    LOGICAL :: mask(kbdim)

    REAL(wp) :: delz(kbdim)

    !$ACC DATA &
    !$ACC   CREATE(loidx, is, se_sum, qv_sum, wgt_sum, wgt, zca, zcs) &
    !$ACC   CREATE(zfrc_oce, zen_h, zfn_h, zen_qv, zfn_qv, zlhflx_lnd) &
    !$ACC   CREATE(zlhflx_lwtr, zlhflx_lice, zshflx_lnd, zshflx_lwtr) &
    !$ACC   CREATE(zshflx_lice, pfrc_test) &
    !$ACC   CREATE(zevap_lnd, zevap_lwtr, zevap_lice) &
    !$ACC   CREATE(qsat_lnd, qsat_lwtr, qsat_lice) &
    !$ACC   CREATE(zcpt_lnd, zcpt_lwtr, zcpt_lice) &
    !$ACC   CREATE(ztsfc_lnd, ztsfc_lnd_eff, ztsfc_wtr, ztsfc_lwtr) &
    !$ACC   CREATE(ztsfc_lice, rvds, rnds, rpds, rsns, rlns) &
    !$ACC   CREATE(fract_par_diffuse, zalbedo_lwtr, zalbedo_lice) &
    !$ACC   CREATE(zgrnd_hflx, zgrnd_hcap, Tfw, swflx_ice, nonsolar_ice) &
    !$ACC   CREATE(dnonsolardT, conc_sum, mask, delz, zwindspeed_lnd) &
    !$ACC   CREATE(zwindspeed10m_lnd) &
    !$ACC   CREATE(rain_tmp, snow_tmp, drag_srf_tmp, pch_tmp, drag_wtr_tmp) &
    !$ACC   CREATE(drag_ice_tmp)

    ! Shortcuts to components of aes_vdf_config
    !
    lsfc_mom_flux  = aes_vdf_config(jg)% lsfc_mom_flux
    lsfc_heat_flux = aes_vdf_config(jg)% lsfc_heat_flux
  
    ! check for masks
    !
    ! DA: compute the index lists on the GPU
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jsfc = 1,ksfc_type
      DO jl = jcs,kproma
        pfrc_test(jl, jsfc) = MERGE(1_i1, 0_i1, pfrc(jl, jsfc) > 0.0_wp)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    CALL generate_index_list_batched(pfrc_test(:,:), loidx, jcs, kproma, is, 1)
    !$ACC UPDATE HOST(is) ASYNC(1)

    ! Compute factor for conversion temperature to dry static energy
    !DO jsfc=1,ksfc_type
    !  zt2s_conv(jcs:kproma,jsfc) = pcpt_tile(jcs:kproma,jsfc) / ptsfc_tile(jcs:kproma,jsfc)
    !END DO

    ! Compute downward shortwave surface fluxes
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs,kproma
      rvds(jl)      = rvds_dif(jl) + rvds_dir(jl)
      rnds(jl)      = rnds_dif(jl) + rnds_dir(jl)
      rpds(jl)      = rpds_dif(jl) + rpds_dir(jl)
    END DO

    !$ACC LOOP GANG VECTOR
    DO jl = jcs,kproma
      delz(jl) = (pmair(jl,klev) / pfac_sfc(jl) / tpfac2 * pdtime)
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT

    ! Turbulent transport of moisture:
    ! - finish matrix set up;
    ! - perform bottom level elimination;
    ! - convert matrix entries to Richtmyer-Morton coefficients
    IF (idx_lnd <= ksfc_type) THEN
      CALL matrix_to_richtmyer_coeff( jcs, kproma, klev, ksfc_type, idx_lnd, &! in
        & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
        & pdtime, delz,                            &! in
        & aa_btm, bb_btm,                          &! inout
        & zen_h, zfn_h, zen_qv, zfn_qv,            &! out
        & pcair = pcair(:),                        &! in
        & pcsat = pcsat(:))                         ! in
    ELSE
      CALL matrix_to_richtmyer_coeff( jcs, kproma, klev, ksfc_type, idx_lnd, &! in
        & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
        & pdtime, delz,                            &! in
        & aa_btm, bb_btm,                          &! inout
        & zen_h, zfn_h, zen_qv, zfn_qv             )! out
    END IF

    ! Set defaults
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jsfc = 1,ksfc_type
      DO jl = jcs,kproma
        zca(jl,jsfc) = 1._wp
        zcs(jl,jsfc) = 1._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !===========================================================================
    ! all surfaces
    !===========================================================================
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jk = 1,kbdim
        rlns_tile(jk,jsfc) = 0._wp
        rsns_tile(jk,jsfc) = 0._wp

        zgrnd_hflx(jk,jsfc)  = 0._wp
        zgrnd_hcap(jk,jsfc)  = 0._wp
      END DO
    END DO
    !===========================================================================
    ! Land surface
    !===========================================================================
    
    !$ACC LOOP GANG VECTOR
    DO jk = 1, kbdim
      zlhflx_lnd(jk)    = 0._wp
      zlhflx_lwtr(jk)   = 0._wp
      zlhflx_lice(jk)   = 0._wp

      zshflx_lnd(jk)    = 0._wp
      zshflx_lwtr(jk)   = 0._wp
      zshflx_lice(jk)   = 0._wp

      zevap_lnd(jk)     = 0._wp
      zevap_lwtr(jk)    = 0._wp
      zevap_lice(jk)    = 0._wp

      z0h_lnd(jk)       = 0._wp
      q_snocpymlt(jk)   = 0._wp
      lake_ice_frc(jk)  = 0._wp

      zalbedo_lwtr(jk)  = 0._wp
      zalbedo_lice(jk)  = 0._wp
    END DO
    !$ACC END PARALLEL

    IF (idx_lnd <= ksfc_type) THEN

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jl = 1, kbdim
        albvisdir_tile(jl,idx_lnd) = 0._wp
        albnirdir_tile(jl,idx_lnd) = 0._wp
        albvisdif_tile(jl,idx_lnd) = 0._wp
        albnirdif_tile(jl,idx_lnd) = 0._wp

        pco2_flux_tile(jl,idx_lnd) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

      ! If land is present, JSBACH is currently the only surface scheme supported by AES physcis package


      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jk = 1, kbdim
        qsat_lnd(jk)          = 0._wp
        qsat_lwtr(jk)         = 0._wp
        qsat_lice(jk)         = 0._wp
        zcpt_lnd(jk)          = 0._wp
        zcpt_lwtr(jk)         = 0._wp
        zcpt_lice(jk)         = 0._wp
        ztsfc_lnd(jk)         = 0._wp
        ztsfc_lnd_eff(jk)     = 0._wp
        ztsfc_lwtr(jk)        = 0._wp
        ztsfc_lice(jk)        = 0._wp
        z0m_tile(jk,idx_lnd)  = 0._wp
        zwindspeed_lnd(jk)    = 0._wp
        zwindspeed10m_lnd(jk) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs, kproma
        zwindspeed_lnd(jl) = SQRT(pu(jl)**2 + pv(jl)**2)
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs, kproma
        zwindspeed10m_lnd(jl)     = 0.8_wp * zwindspeed_lnd(jl)
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs, kproma
        IF (rpds(jl) > 0._wp) THEN
          fract_par_diffuse(jl) = rpds_dif(jl) / rpds(jl)
        ELSE
          fract_par_diffuse(jl) = 0._wp
        END IF
      END DO

      ! Prepare temporary fields to be passed to JSBACH since they cannot be
      ! computed in arguments.

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs, kproma
        rain_tmp(jl) = prsfl(jl)
        snow_tmp(jl) = pssfl(jl)
        drag_srf_tmp(jl) = grav*pfac_sfc(jl) * pcfh_tile(jl,idx_lnd)
        IF(lsm(jl)>0._wp) THEN
          pch_tmp(jl) = pch_tile(jl,idx_lnd) ! MERGE(pch_tile(jcs:kproma,idx_lnd),1._wp,lsm(jcs:kproma)>0._wp)
        ELSE
          pch_tmp(jl) = 1._wp
        END IF
        IF (aes_phy_config(jg)%llake) THEN
          drag_wtr_tmp(jl) = grav*pfac_sfc(jl) * pcfh_tile(jl,idx_wtr)
          drag_ice_tmp(jl) = grav*pfac_sfc(jl) * pcfh_tile(jl,idx_ice)
        END IF
      END DO
      !$ACC END PARALLEL

      !$ACC WAIT

      IF (aes_phy_config(jg)%ljsb ) THEN
      IF (aes_phy_config(jg)%llake) THEN
        CALL jsbach_interface ( jg, nblock, jcs, kproma,                                     & ! in
          & datetime, pdtime, pdtime,                                                        & ! in
          & t_air             = ptemp(jcs:kproma),                                           & ! in
          & q_air             = pq(jcs:kproma),                                              & ! in
          & rain              = rain_tmp(jcs:kproma),                                        & ! in
          & snow              = snow_tmp(jcs:kproma),                                        & ! in
          & wind_air          = zwindspeed_lnd(jcs:kproma),                                  & ! in
          & wind_10m          = zwindspeed10m_lnd(jcs:kproma),                               & ! in
          & lw_srf_down       = rlds(jcs:kproma),                                            & ! in
          & swvis_srf_down    = rvds(jcs:kproma),                                            & ! in
          & swnir_srf_down    = rnds(jcs:kproma),                                            & ! in
          & swpar_srf_down    = rpds(jcs:kproma),                                            & ! in
          & fract_par_diffuse = fract_par_diffuse(jcs:kproma),                               & ! in
          & press_srf         = ps(jcs:kproma),                                              & ! in
          & drag_srf          = drag_srf_tmp(jcs:kproma),                                    & ! in
          & t_acoef           = zen_h(jcs:kproma, idx_lnd),                                  & ! in
          & t_bcoef           = zfn_h(jcs:kproma, idx_lnd),                                  & ! in
          & q_acoef           = zen_qv(jcs:kproma, idx_lnd),                                 & ! in
          & q_bcoef           = zfn_qv(jcs:kproma, idx_lnd),                                 & ! in
          & pch               = pch_tmp(jcs:kproma),                                         & ! in
          & cos_zenith_angle  = pcosmu0(jcs:kproma),                                         & ! in
          & CO2_air           = pco2(jcs:kproma),                                            & ! in
          & t_srf             = ztsfc_lnd(jcs:kproma),                                       & ! out (T_s^(n+1)) surface temp
                                                                                               ! (filtered, if Asselin)
          & t_eff_srf         = ztsfc_lnd_eff(jcs:kproma),                                   & ! out (T_s^eff) surface temp
                                                                                               ! (effective, for longwave rad)
          & qsat_srf          = qsat_lnd(jcs:kproma),                                        & ! out
          & s_srf             = zcpt_lnd(jcs:kproma),                                        & ! out (s_s^star, for vdiff scheme)
          & fact_q_air        = pcair(jcs:kproma),                                           & ! out
          & fact_qsat_srf     = pcsat(jcs:kproma),                                           & ! out
          & evapotrans        = zevap_lnd(jcs:kproma),                                       & ! out
          & latent_hflx       = zlhflx_lnd(jcs:kproma),                                      & ! out
          & sensible_hflx     = zshflx_lnd(jcs:kproma),                                      & ! out
          & grnd_hflx         = zgrnd_hflx(jcs:kproma, idx_lnd),                             & ! out
          & grnd_hcap         = zgrnd_hcap(jcs:kproma, idx_lnd),                             & ! out
          & rough_h_srf       = z0h_lnd(jcs:kproma),                                         & ! out
          & rough_m_srf       = z0m_tile(jcs:kproma, idx_lnd),                               & ! out
          & q_snocpymlt       = q_snocpymlt(jcs:kproma),                                     & ! out
          & alb_vis_dir       = albvisdir_tile(jcs:kproma, idx_lnd),                         & ! out
          & alb_nir_dir       = albnirdir_tile(jcs:kproma, idx_lnd),                         & ! out
          & alb_vis_dif       = albvisdif_tile(jcs:kproma, idx_lnd),                         & ! out
          & alb_nir_dif       = albnirdif_tile(jcs:kproma, idx_lnd),                         & ! out
          & co2_flux          = pco2_flux_tile(jcs:kproma, idx_lnd),                         & ! out
          !
          & drag_wtr          = drag_wtr_tmp(jcs:kproma),                                    & ! in
          & drag_ice          = drag_ice_tmp(jcs:kproma),                                    & ! in
          & t_acoef_wtr       = zen_h(jcs:kproma, idx_wtr),                                  & ! in
          & t_bcoef_wtr       = zfn_h(jcs:kproma, idx_wtr),                                  & ! in
          & q_acoef_wtr       = zen_qv(jcs:kproma, idx_wtr),                                 & ! in
          & q_bcoef_wtr       = zfn_qv(jcs:kproma, idx_wtr),                                 & ! in
          & t_acoef_ice       = zen_h(jcs:kproma, idx_ice),                                  & ! in
          & t_bcoef_ice       = zfn_h(jcs:kproma, idx_ice),                                  & ! in
          & q_acoef_ice       = zen_qv(jcs:kproma, idx_ice),                                 & ! in
          & q_bcoef_ice       = zfn_qv(jcs:kproma, idx_ice),                                 & ! in
          & t_lwtr            = ztsfc_lwtr(jcs:kproma),                                      & ! out
          & qsat_lwtr         = qsat_lwtr(jcs:kproma),                                       & ! out
          & s_lwtr            = zcpt_lwtr(jcs:kproma),                                       & ! out
          & evapo_wtr         = zevap_lwtr(jcs:kproma),                                      & ! out
          & latent_hflx_wtr   = zlhflx_lwtr(jcs:kproma),                                     & ! out
          & sensible_hflx_wtr = zshflx_lwtr(jcs:kproma),                                     & ! out
          & albedo_lwtr       = zalbedo_lwtr(jcs:kproma),                                    & ! out
          & t_lice            = ztsfc_lice(jcs:kproma),                                      & ! out
          & qsat_lice         = qsat_lice(jcs:kproma),                                       & ! out
          & s_lice            = zcpt_lice(jcs:kproma),                                       & ! out
          & evapo_ice         = zevap_lice(jcs:kproma),                                      & ! out
          & latent_hflx_ice   = zlhflx_lice(jcs:kproma),                                     & ! out
          & sensible_hflx_ice = zshflx_lice(jcs:kproma),                                     & ! out
          & albedo_lice       = zalbedo_lice(jcs:kproma),                                    & ! out
          & ice_fract_lake    = lake_ice_frc(jcs:kproma)                                     & ! out
          )


      ELSE
        CALL jsbach_interface ( jg, nblock, jcs, kproma,                                     & ! in
          & datetime, pdtime, pdtime,                                                        & ! in
          & t_air             = ptemp(jcs:kproma),                                           & ! in
          & q_air             = pq(jcs:kproma),                                              & ! in
          & rain              = rain_tmp(jcs:kproma),                                        & ! in
          & snow              = snow_tmp(jcs:kproma),                                        & ! in
          & wind_air          = zwindspeed_lnd(jcs:kproma),                                  & ! in
          & wind_10m          = zwindspeed10m_lnd(jcs:kproma),                               & ! in
          & lw_srf_down       = rlds(jcs:kproma),                                            & ! in
          & swvis_srf_down    = rvds(jcs:kproma),                                            & ! in
          & swnir_srf_down    = rnds(jcs:kproma),                                            & ! in
          & swpar_srf_down    = rpds(jcs:kproma),                                            & ! in
          & fract_par_diffuse = fract_par_diffuse(jcs:kproma),                               & ! in
          & press_srf         = ps(jcs:kproma),                                              & ! in
          & drag_srf          = drag_srf_tmp(jcs:kproma),                                    & ! in
          & t_acoef           = zen_h(jcs:kproma, idx_lnd),                                  & ! in
          & t_bcoef           = zfn_h(jcs:kproma, idx_lnd),                                  & ! in
          & q_acoef           = zen_qv(jcs:kproma, idx_lnd),                                 & ! in
          & q_bcoef           = zfn_qv(jcs:kproma, idx_lnd),                                 & ! in
          & pch               = pch_tmp(jcs:kproma),                                         & ! in
          & cos_zenith_angle  = pcosmu0(jcs:kproma),                                         & ! in
          & CO2_air           = pco2(jcs:kproma),                                            & ! in
          & t_srf             = ztsfc_lnd(jcs:kproma),                                       & ! out (T_s^(n+1)) surface temp 
                                                                                               ! (filtered, if Asselin)
          & t_eff_srf         = ztsfc_lnd_eff(jcs:kproma),                                   & ! out (T_s^eff) surface temp 
                                                                                               ! (effective, for longwave rad)
          & qsat_srf          = qsat_lnd(jcs:kproma),                                        & ! out
          & s_srf             = zcpt_lnd(jcs:kproma),                                        & ! out (s_s^star, for vdiff scheme)
          & fact_q_air        = pcair(jcs:kproma),                                           & ! out
          & fact_qsat_srf     = pcsat(jcs:kproma),                                           & ! out
          & evapotrans        = zevap_lnd(jcs:kproma),                                       & ! out
          & latent_hflx       = zlhflx_lnd(jcs:kproma),                                      & ! out
          & sensible_hflx     = zshflx_lnd(jcs:kproma),                                      & ! out
          & grnd_hflx         = zgrnd_hflx(jcs:kproma, idx_lnd),                             & ! out
          & grnd_hcap         = zgrnd_hcap(jcs:kproma, idx_lnd),                             & ! out
          & rough_h_srf       = z0h_lnd(jcs:kproma),                                         & ! out
          & rough_m_srf       = z0m_tile(jcs:kproma, idx_lnd),                               & ! out
          & q_snocpymlt       = q_snocpymlt(jcs:kproma),                                     & ! out
          & alb_vis_dir       = albvisdir_tile(jcs:kproma, idx_lnd),                         & ! out
          & alb_nir_dir       = albnirdir_tile(jcs:kproma, idx_lnd),                         & ! out
          & alb_vis_dif       = albvisdif_tile(jcs:kproma, idx_lnd),                         & ! out
          & alb_nir_dif       = albnirdif_tile(jcs:kproma, idx_lnd),                         & ! out
          & co2_flux          = pco2_flux_tile(jcs:kproma, idx_lnd)                          & ! out
        )



      END IF ! llake
      END IF ! ljsb

      IF (aes_phy_config(jg)%ljsb ) THEN
        ptsfc_tile(jcs:kproma,idx_lnd) = ztsfc_lnd(jcs:kproma)
        pcpt_tile (jcs:kproma,idx_lnd) = zcpt_lnd(jcs:kproma)
        pqsat_tile(jcs:kproma,idx_lnd) = qsat_lnd(jcs:kproma)
      END IF
      IF (aes_phy_config(jg)%llake) THEN
        IF (idx_wtr <= ksfc_type) THEN
          WHERE (alake(jcs:kproma) > 0._wp)
            ptsfc_tile    (jcs:kproma, idx_wtr) = ztsfc_lwtr   (jcs:kproma)
            pcpt_tile     (jcs:kproma, idx_wtr) = zcpt_lwtr    (jcs:kproma)
            pqsat_tile    (jcs:kproma, idx_wtr) = qsat_lwtr    (jcs:kproma)
            albvisdir_tile(jcs:kproma, idx_wtr) = zalbedo_lwtr (jcs:kproma)
            albvisdif_tile(jcs:kproma, idx_wtr) = zalbedo_lwtr (jcs:kproma)
            albnirdir_tile(jcs:kproma, idx_wtr) = zalbedo_lwtr (jcs:kproma)
            albnirdif_tile(jcs:kproma, idx_wtr) = zalbedo_lwtr (jcs:kproma)
          ! security reasons
            pco2_flux_tile(jcs:kproma, idx_wtr) =  0._wp
          END WHERE
        END IF
        IF (idx_ice <= ksfc_type) THEN
          WHERE (alake(jcs:kproma) > 0._wp)
            ptsfc_tile    (jcs:kproma, idx_ice) = ztsfc_lice   (jcs:kproma)
            pcpt_tile     (jcs:kproma, idx_ice) = zcpt_lice    (jcs:kproma)
            pqsat_tile    (jcs:kproma, idx_ice) = qsat_lice    (jcs:kproma)
            albvisdir_tile(jcs:kproma, idx_ice) = zalbedo_lice (jcs:kproma)
            albvisdif_tile(jcs:kproma, idx_ice) = zalbedo_lice (jcs:kproma)
            albnirdir_tile(jcs:kproma, idx_ice) = zalbedo_lice (jcs:kproma)
            albnirdif_tile(jcs:kproma, idx_ice) = zalbedo_lice (jcs:kproma)
          ELSEWHERE
            lake_ice_frc(jcs:kproma) = 0._wp
          ENDWHERE
        END IF
      END IF

      ! Set the evapotranspiration coefficients, to be used later in
      ! blending and in diagnosing surface fluxes.
      !
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jl = jcs,kproma
        zca(jl,idx_lnd) = pcair(jl)
        zcs(jl,idx_lnd) = pcsat(jl)
      END DO
      !$ACC END PARALLEL LOOP
    END IF ! idx_lnd

    !===================================================================
    ! AFTER CALLING land/ocean/ice model:
    ! pco2_flux_tile set to zero for sea ice and ocean above...
    !===================================================================

    ! calculate grid box mean surface of co2
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jk = 1, kbdim
      pco2nat(jk) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jsfc=1,ksfc_type
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs,kproma
        pco2nat(jl) = pco2nat(jl) + pfrc(jl,jsfc) * pco2_flux_tile(jl,jsfc)
      END DO
    ENDDO

    !
    ! Turbulent transport of moisture and dry static energy:
    ! Get solution of the two variables on the lowest model level.
    !-------------------------------------------------------------------
    ! - Over individual tiles

    IF ( isrfc_type == 1) THEN
      !$ACC LOOP SEQ
      DO jsfc = 1,ksfc_type
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jl = jcs,kproma
          bb_btm(jl,jsfc,ih)  =  zen_h (jl,jsfc) + tpfac2*zfn_h (jl,jsfc)
          bb_btm(jl,jsfc,iqv) = zen_qv (jl,jsfc) + tpfac2*zfn_qv (jl,jsfc)
        END DO
      END DO
    ELSE
      !$ACC LOOP SEQ
      DO jsfc = 1,ksfc_type
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jl = jcs,kproma
          bb_btm(jl,jsfc,ih)  = tpfac2*(    zen_h (jl,jsfc)                      &
                              &         *pcpt_tile(jl,jsfc)                      &
                              &         +   zfn_h (jl,jsfc) )

          bb_btm(jl,jsfc,iqv) = tpfac2*(    zen_qv(jl,jsfc)                      &
                              &        *pqsat_tile(jl,jsfc)                      &
                              &        +    zfn_qv(jl,jsfc) )
        END DO
      END DO
    END IF

    ! - Grid box mean

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl = jcs,kproma
       se_sum(jl) = 0._wp    ! sum of weighted solution
       qv_sum(jl) = 0._wp    ! sum of weighted solution
      wgt_sum(jl) = 0._wp    ! sum of weights
    END DO

    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs,kproma
             wgt(jl) = pfrc(jl,jsfc)
         wgt_sum(jl) = wgt_sum(jl) + wgt(jl)
          se_sum(jl) = se_sum(jl) + bb_btm(jl,jsfc,ih ) * wgt(jl)
          qv_sum(jl) = qv_sum(jl) + bb_btm(jl,jsfc,iqv) * wgt(jl)
      ENDDO
    ENDDO

    IF (lsfc_heat_flux) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs,kproma
        bb(jl,klev,ih ) = se_sum(jl)/wgt_sum(jl)
        bb(jl,klev,iqv) = qv_sum(jl)/wgt_sum(jl)
      END DO
    ELSE
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jsfc)
      DO jl = jcs,kproma
        jsfc = 1
        bb(jl,klev,ih ) = bb_btm(jl,jsfc,ih )
        bb(jl,klev,iqv) = bb_btm(jl,jsfc,iqv)
      END DO
    END IF

    !-------------------------------------------------------------------
    ! Turbulent transport of u and v: adjust the right-hand side vector,
    ! then perform the bottom level elimination to get the solution
    !-------------------------------------------------------------------
    ! Add additional terms to the r.h.s. of the velocity equations
    ! to take into account ocean currents.
    ! Note that in subroutine rhs_setup the constant tpfac2 has been
    ! multiplied to the r.h.s. array bb. Thus the additional terms here
    ! need to be scaled by the same factor.

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl = jcs,kproma
      zfrc_oce(jl) = 0._wp 
    END DO
    IF (idx_wtr.LE.ksfc_type) THEN   ! Open water is considered
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs,kproma
        IF (idx_ice.LE.ksfc_type) THEN ! Sea ice is also considered
          zfrc_oce(jl) = pfrc(jl,idx_wtr)+pfrc(jl,idx_ice)
        ELSE ! only open water
          zfrc_oce(jl) = pfrc(jl,idx_wtr)
        ENDIF
        IF (idx_lnd.LE.ksfc_type) THEN
          IF (alake(jl) > 0._wp) THEN
            zfrc_oce(jl) = 0._wp
          END IF
        END IF
      END DO
    ENDIF

    ! Bottom level elimination

    im   = imuv
    jk   = klev    ! Bottom level index
    jkm1 = jk - 1

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl = jcs,kproma
      aa(jl,jk,2,im) =  aa(jl,jk,2,im) - aa(jl,jk,1,im)*aa(jl,jkm1,3,im)
      aa(jl,jk,3,im) =  aa(jl,jk,3,im)/aa(jl,jk,2,im)

      bb(jl,jk,iu) = - aa(jl,jk,3,im) * pocu(jl)*zfrc_oce(jl)*tpfac2 + (bb(jl,jk,iu) - aa(jl,jk,1,im)*bb(jl,jkm1,iu))/aa(jl,jk,2,im)
      bb(jl,jk,iv) = - aa(jl,jk,3,im) * pocv(jl)*zfrc_oce(jl)*tpfac2 + (bb(jl,jk,iv) - aa(jl,jk,1,im)*bb(jl,jkm1,iv))/aa(jl,jk,2,im)
    END DO
    !$ACC END PARALLEL

    ! Compute wind stress
    IF (lsfc_mom_flux) THEN
       CALL wind_stress( kbdim, ksfc_type,                    &! in
            &            pdtime,                              &! in
            &            loidx, is,                           &! in
            &            pfrc, pcfm_tile, pfac_sfc,           &! in
            &            bb(:,klev,iu), bb(:,klev,iv),        &! in
            &            pocu(:), pocv(:),                    &! in
            &            pu_stress_gbm,  pv_stress_gbm,       &! out
            &            pu_stress_tile, pv_stress_tile       )! out
    ELSE
       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       !$ACC LOOP SEQ
       DO jsfc = 1,ksfc_type
         !$ACC LOOP GANG VECTOR
         DO jk = 1, kbdim
           pu_stress_tile(jk, jsfc) = 0._wp
           pv_stress_tile(jk, jsfc) = 0._wp
         END DO
       END DO

       !$ACC LOOP GANG VECTOR
       DO jk = 1, kbdim
         pu_stress_gbm (jk)   = 0._wp
         pv_stress_gbm (jk)   = 0._wp
       END DO
       !$ACC END PARALLEL

    END IF

   !--------------------------------------------------------------------
   ! Various diagnostics
   ! Calculate turbulent heat fluxes before calling sea ice or mlo-ocean
   !--------------------------------------------------------------------


    IF (lsfc_heat_flux) THEN
       !$ACC WAIT
       CALL surface_fluxes( jcs, kproma, kbdim, ksfc_type,        &! in
            &               idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
            &               pdtime,                               &! in
            &               pfrc, lsm, alake,                     &! in
            &               pcfh_tile, pfac_sfc,                  &! in
            &               pcpt_tile, pqsat_tile,                &! in
            &               zca, zcs, bb_btm(:,:,ih:iqv),         &! in
            &               zlhflx_lnd, zlhflx_lwtr, zlhflx_lice, &! in
            &               zshflx_lnd, zshflx_lwtr, zshflx_lice, &! in
            &               zevap_lnd, zevap_lwtr, zevap_lice,    &! in
            &               plhflx_gbm, pshflx_gbm,               &! out
            &               pevap_gbm,                            &! out
            &               plhflx_tile, pshflx_tile,             &! out
            &               pevap_tile )                           ! out
      IF ( isrfc_type == 1 ) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jl = jcs,kproma
          pshflx_gbm (jl)   = -shflx*cpd*pfac_sfc(jl)*tpfac2/pdtime
          plhflx_gbm (jl)   = -lhflx*alv*pfac_sfc(jl)*tpfac2/pdtime
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    ELSE
       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       !$ACC LOOP SEQ
       DO jsfc = 1,ksfc_type
         !$ACC LOOP GANG VECTOR
         DO jl = 1, kbdim
           plhflx_tile(jl,jsfc) = 0._wp
           pshflx_tile(jl,jsfc) = 0._wp
           pevap_tile (jl,jsfc) = 0._wp
         END DO
       END DO
       !$ACC LOOP GANG VECTOR
       DO jk = 1,kbdim
         plhflx_gbm (jk)   = 0._wp
         pshflx_gbm (jk)   = 0._wp
         pevap_gbm  (jk)   = 0._wp
       END DO
       !$ACC END PARALLEL
    END IF


    !===========================================================================
    ! Ocean model
    !===========================================================================
    IF (idx_wtr <= ksfc_type) THEN

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jl = jcs,kproma
        rsns(jl)      = rsds(jl) - rsus(jl)
        rlns(jl)      = rlds(jl) - rlus(jl)
      END DO
      !$ACC END PARALLEL LOOP

      IF (aes_phy_config(jg)%lmlo) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jl = jcs,kproma
          ztsfc_wtr(jl)=ptsfc_tile(jl, idx_wtr)
        ENDDO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT
        CALL ml_ocean ( kbdim, jcs, kproma, pdtime, &
          & pahflw=plhflx_tile(:,idx_wtr),        & ! dependency on kproma has to be checked
          & pahfsw=pshflx_tile(:,idx_wtr),        & ! dependency on kproma has to be checked
          & ptrflw=rlns(:),                       &
          & psoflw=rsns(:),                       &
          & ptsw=ztsfc_wtr(:) )                     ! out
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jl = jcs,kproma
          IF (alake(jl) < EPSILON(1._wp)) THEN
            ptsfc_tile(jl, idx_wtr) = ztsfc_wtr(jl)
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Albedo model for the ocean
      ! TBD: This should be replaced by routine
      ! mo_surface_ocean:update_albedo_ocean from ECHAM6.2
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jl = jcs,kproma
        IF (alake(jl) < EPSILON(1._wp)) THEN
          albvisdir_tile(jl,idx_wtr) = albedoW
          albvisdif_tile(jl,idx_wtr) = albedoW
          albnirdir_tile(jl,idx_wtr) = albedoW
          albnirdif_tile(jl,idx_wtr) = albedoW
        END IF
      END DO
      !$ACC END PARALLEL LOOP

    END IF ! idx_wtr

    !===========================================================================
    ! Sea-ice model (thermodynamic)
    !===========================================================================

    IF (idx_ice <= ksfc_type .AND. aes_phy_config(jg)%lice) THEN

      ! DA: wait while the other kernels are not async
      !$ACC WAIT

      ! LL This is a temporary solution,
      ! we should restrcure ice thermodynamics in a more stand-alone way

      ! For explicit coupling to ice:

      ! Freezing point of sea-water
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jk = 1, kbdim
        Tfw(jk) = Tf
      END DO
      !$ACC END PARALLEL LOOP

      ! ECHAM has no tiles for SW & LW and this is how it's solved there
      ! Net shortwave on all bands.
      ! Net longwave - we don't have tiles yet
      ! First all ice classes

      IF (aes_phy_config(jg)%use_shflx_adjustment .AND. &
          .NOT. aes_phy_config(jg)%suppress_shflx_adjustment_over_ice) THEN
  

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO k=1,kice
          DO jl = jcs,kproma
            swflx_ice(jl,k) = rvds_dif(jl) * (1._wp - albvisdif_ice(jl,k)) +     &
                            & rvds_dir(jl) * (1._wp - albvisdir_ice(jl,k)) +     &
                            & rnds_dif(jl) * (1._wp - albnirdif_ice(jl,k)) +     &
                            & rnds_dir(jl) * (1._wp - albnirdir_ice(jl,k))

            nonsolar_ice(jl,k) = &
              emissivity(jl) * (rlds(jl) - stbo * (Tsurf(jl,k)+tmelt)**4) &  !  longwave net
              & + plhflx_tile(jl,idx_ice) + (cvd/cpd)*pshflx_tile(jl,idx_ice)

            dnonsolardT(jl,k) = -4._wp * emissivity(jl) * stbo * (Tsurf(jl,k)+tmelt)**3
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

      ELSE ! .NOT. aes_phy_config(jg)%use_shflx_adjustment .OR.
           ! aes_phy_config(jg)%suppress_shflx_adjustment_over_ice

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO k=1,kice
          DO jl = jcs,kproma
            swflx_ice(jl,k) = rvds_dif(jl) * (1._wp - albvisdif_ice(jl,k)) +     &
                            & rvds_dir(jl) * (1._wp - albvisdir_ice(jl,k)) +     &
                            & rnds_dif(jl) * (1._wp - albnirdif_ice(jl,k)) +     &
                            & rnds_dir(jl) * (1._wp - albnirdir_ice(jl,k))

            nonsolar_ice(jl,k) = &
              emissivity(jl) * (rlds(jl) - stbo * (Tsurf(jl,k)+tmelt)**4) &  !  longwave net
              & + plhflx_tile(jl,idx_ice) + pshflx_tile(jl,idx_ice)

            dnonsolardT(jl,k) = -4._wp * emissivity(jl) * stbo * (Tsurf(jl,k)+tmelt)**3
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

      ENDIF ! aes_phy_config(jg)%use_shflx_adjustment .AND.
            ! .NOT. aes_phy_config(jg)%suppress_shflx_adjustment_over_ice

      !$ACC WAIT
      CALL ice_fast(jcs, kproma, kbdim, kice, pdtime, &
        &   Tsurf,              & ! inout
        &   T1,                 & ! inout
        &   T2,                 & ! inout
        &   hi,                 & ! in
        &   hs,                 & ! in
        &   Qtop,               & ! out
        &   Qbot,               & ! out
        &   swflx_ice,          & ! in
        &   nonsolar_ice,       & ! in
        &   dnonsolardT,        & ! in
        &   Tfw,                & ! in
        &   albvisdir_ice,      & ! out
        &   albvisdif_ice,      & ! out
        &   albnirdir_ice,      & ! out
        &   albnirdif_ice,      & ! out
        &   lacc=.TRUE.)          ! in

      ! Update the thickness of snow on ice in atmosphere only simulation.
      ! In coupled experiments this is done by the ocean model in either
      ! ice_growth_zerolayer or ice_growth_winton.
      IF ( .NOT. is_coupled_run() ) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO k=1,kice
          DO jl = jcs,kproma
            ! Snowfall on ice - no ice => no snow
            IF ( hi(jl,k) > 0._wp ) THEN
              ! Snow only falls when it's below freezing
              IF ( Tsurf(jl,k) < 0._wp ) THEN
                hs(jl,k) = hs(jl,k) + pssfl(jl)*pdtime/rhos
              ENDIF
              ! Snow melt
              hs(jl,k) = hs(jl,k) - MIN( Qtop(jl,k)*pdtime/( alf*rhos ), hs(jl,k) )
            ELSE
              hs(jl,k) = 0._wp
            ENDIF
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      ! Average the albedo.
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jl = jcs,kproma
        conc_sum(jl) = SUM(conc(jl,:))
        IF (alake(jl) < EPSILON(1._wp)) THEN
          albvisdir_tile(jl,idx_ice) = 0._wp
          albvisdif_tile(jl,idx_ice) = 0._wp
          albnirdir_tile(jl,idx_ice) = 0._wp
          albnirdif_tile(jl,idx_ice) = 0._wp
          IF (conc_sum(jl) > 1.e-6_wp) THEN
            albvisdir_tile(jl,idx_ice) = SUM( conc(jl,:) * albvisdir_ice(jl,:)) / conc_sum(jl)
            albvisdif_tile(jl,idx_ice) = SUM( conc(jl,:) * albvisdif_ice(jl,:)) / conc_sum(jl)
            albnirdir_tile(jl,idx_ice) = SUM( conc(jl,:) * albnirdir_ice(jl,:)) / conc_sum(jl)
            albnirdif_tile(jl,idx_ice) = SUM( conc(jl,:) * albnirdif_ice(jl,:)) / conc_sum(jl)

            ! Set the tile temperature, convert back to K
            ptsfc_tile(jl,idx_ice) = Tsurf(jl,1) + tmelt
          END IF
        END IF
      END DO
      !$ACC END PARALLEL LOOP

      ! Compute new dry static energy
      ! (Switched off for now, should be used for implicit coupling)
      !pcpt_tile(jcs:kproma,idx_ice) = ptsfc_tile(jcs:kproma,idx_ice) * zt2s_conv(jcs:kproma,idx_ice)

    ENDIF ! lice

   !-------------------------------------------------------------------
   ! More diagnostics
   !-------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jsfc=1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jl = 1, kproma
        albedo_tile(jl,jsfc) = 0._wp
      END DO

      !$ACC LOOP GANG VECTOR PRIVATE(zalbvis, zalbnir)
      DO jl = jcs, kproma
        zalbvis = 0._wp
        IF(rvds(jl) > 0._wp) THEN
          zalbvis = &
            & (albvisdir_tile(jl,jsfc) * rvds_dir(jl) + albvisdif_tile(jl,jsfc) * rvds_dif(jl)) &
            & / rvds(jl)
        END IF
        zalbnir = 0._wp
        IF(rnds(jl) > 0._wp) THEN
          zalbnir = &
            & (albnirdir_tile(jl,jsfc) * rnds_dir(jl) + albnirdif_tile(jl,jsfc) * rnds_dif(jl)) &
            & / rnds(jl)
        END IF
        IF(rsds(jl) > 0._wp) THEN
          albedo_tile(jl,jsfc) = &
            & (zalbvis * rvds(jl) + zalbnir * rnds(jl)) &
            & / rsds(jl)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    ! calculate grid box mean surface temperature
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jk = 1, kbdim
      ptsfc(jk) = 0._wp
      ptsfc_rad(jk) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jsfc=1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jl= jcs, kproma
        ptsfc(jl) = ptsfc(jl) + pfrc(jl,jsfc) * ptsfc_tile(jl,jsfc)
      ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO jsfc=1,ksfc_type
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl= jcs, kproma
        ptsfc_rad(jl) = ptsfc_rad(jl) + pfrc(jl,jsfc) * ptsfc_tile(jl,jsfc)**4
      END DO
    ENDDO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl = jcs, kproma
      ptsfc_rad(jl) = ptsfc_rad(jl)**0.25_wp
    END DO
    !$ACC END PARALLEL

    ! Compute lw and sw surface radiation fluxes on tiles
    ! DA: can't collapse due to the jls bounds depending on jsfc
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jsfc=1,ksfc_type
      !$ACC LOOP GANG VECTOR PRIVATE(js)
      DO jls = 1,is(jsfc)
        ! set index
        js=loidx(jls,jsfc)
        IF (jsfc == idx_lnd) THEN
          rlns_tile(js,jsfc) = emissivity(js) * (rlds(js) - stbo * ztsfc_lnd_eff(js)**4)
        ELSE
          rlns_tile(js,jsfc) = emissivity(js) * (rlds(js) - stbo * ptsfc_tile(js,jsfc)**4)
        END IF

        rsns_tile(js,jsfc) = rvds_dif(js) * (1._wp - albvisdif_tile(js,jsfc)) + &
                           & rvds_dir(js) * (1._wp - albvisdir_tile(js,jsfc)) + &
                           & rnds_dif(js) * (1._wp - albnirdif_tile(js,jsfc)) + &
                           & rnds_dir(js) * (1._wp - albnirdir_tile(js,jsfc))
      END DO
    END DO
    !$ACC END PARALLEL

    ! Merge sw and lw surface fluxes
    ! This includes the update of the lw flux on land due to the new surface temperature where only part
    ! of the net radiation was used (due to the Taylor truncation in the surface energy balance)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jk = 1, kbdim
      rlns(jk) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !DA: can't move this loop into OpenACC w/o atomics
    DO jsfc=1,ksfc_type
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR PRIVATE(js) ASYNC(1)
      DO jls = 1,is(jsfc)
        ! set index
        js=loidx(jls,jsfc)
        rlns(js) = rlns(js) + pfrc(js,jsfc) * rlns_tile(js,jsfc)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jl = jcs, kproma
      rlus(jl) = rlds(jl) -rlns(jl)
    END DO
    !$ACC END PARALLEL LOOP

    ! Merge surface albedos
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jk = 1, kbdim
      albvisdir(jk) = 0._wp
      albvisdif(jk) = 0._wp
      albnirdir(jk) = 0._wp
      albnirdif(jk) = 0._wp
      albedo   (jk) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jsfc=1,nsfc_type
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs, kproma
        albvisdir(jl) = albvisdir(jl) + pfrc(jl,jsfc) * albvisdir_tile(jl,jsfc)
        albvisdif(jl) = albvisdif(jl) + pfrc(jl,jsfc) * albvisdif_tile(jl,jsfc)
        albnirdir(jl) = albnirdir(jl) + pfrc(jl,jsfc) * albnirdir_tile(jl,jsfc)
        albnirdif(jl) = albnirdif(jl) + pfrc(jl,jsfc) * albnirdif_tile(jl,jsfc)
        albedo   (jl) = albedo   (jl) + pfrc(jl,jsfc) * albedo_tile   (jl,jsfc)
      END DO
    END DO

    ! Mask out tiled variables
    !$ACC LOOP SEQ
    DO jsfc=1,ksfc_type
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs, kproma
        mask(jl) = pfrc(jl,jsfc) <= 0._wp
        !
        IF (mask(jl)) THEN
          pqsat_tile     (jl,jsfc) = cdimissval
          albedo_tile    (jl,jsfc) = cdimissval
          albvisdir_tile (jl,jsfc) = cdimissval
          albvisdif_tile (jl,jsfc) = cdimissval
          albnirdir_tile (jl,jsfc) = cdimissval
          albnirdif_tile (jl,jsfc) = cdimissval
        !  rsns_tile      (jl,jsfc) = cdimissval
        !  rlns_tile      (jl,jsfc) = cdimissval
          pevap_tile     (jl,jsfc) = cdimissval
        !  pshflx_tile    (jl,jsfc) = cdimissval
        !  plhflx_tile    (jl,jsfc) = cdimissval
          ptsfc_tile     (jl,jsfc) = cdimissval
        END IF
        ! land only
        IF (jsfc == idx_lnd) THEN
          IF (mask(jl)) THEN
            z0h_lnd        (jl)      = cdimissval
            z0m_tile       (jl,jsfc) = cdimissval
            rsns_tile      (jl,jsfc) = cdimissval
            rlns_tile      (jl,jsfc) = cdimissval
            pshflx_tile    (jl,jsfc) = cdimissval
            plhflx_tile    (jl,jsfc) = cdimissval
          END IF
        END IF
      END DO
    END DO

    !----------------------------------------------------------------------------
    ! For consistency z0m_tile for ice is masked out here
    !----------------------------------------------------------------------------
    IF (idx_ice<=ksfc_type) THEN  ! ice surface exists in the simulation
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl = jcs, kproma
        mask(jl) = pfrc(jl,idx_ice) == 0._wp
        IF (mask(jl)) THEN
          z0m_tile(jl,idx_ice) = cdimissval
        ELSE
          ! z0m for ice is not calculated yet, so in case the ice surface changes
          ! it is set to the initial value again
          z0m_tile(jl,idx_ice) = aes_vdf_config(jg)% z0m_ice
        ENDIF
      END DO
    ENDIF
    !$ACC END PARALLEL

    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE update_surface
  !-------------

END MODULE mo_surface
