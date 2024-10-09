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

! Module to provide updated radiative fluxes and heating rates
!
!   This module contains the "radheating" routine that diagnoses
!   SW and LW fluxes at TOA and at the surface and the heating
!   in the atmosphere for the current time.

MODULE mo_radheating

  USE mo_kind               , ONLY: wp
  USE mo_physical_constants , ONLY: stbo

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: radheating


CONTAINS

  !-----------------------------------------------------------------------------
  !>
  !! Compute shortwave and longwave heating rates
  !!
  !! The radheating subroutine computes the radiative heating rates in the
  !! shortwave (SW) and longwave (LW) part of the spectrum.
  !!
  !! SW and LW fluxes are computed only every radiation time step. The resulting
  !! upward and downward fluxes at all half-levels are used here to approximate
  !! fluxes at the current time, as follows:
  !!
  !! SW fluxes are rescaled by the SW indicent radiation at the top of the
  !! atmosphere. In order to make this work the radiative transfer uses a zenith
  !! angle limited to 84deg so that also at night time non-zero fluxes exist,
  !! as needed for the scaling.
  !!
  !! LW fluxes are kept constant except for the upward flux from the surface,
  !! which is corrected for the new radiative surface temperature.
  !!
  !! The corrected fluxes are used to diagnose various SW and LW flux components
  !! at the top of the atmosphere and at the surface, and to compute the SW and LW
  !! heating from the vertical flux convergences.
  !!
  !! The vertical ordering is assumed to be from the top of the atmosphere
  !! downward to the surface.
  !!

  SUBROUTINE radheating ( &
       !
       ! input
       ! -----
       !
       & jcs        ,&
       & jce        ,&
       & kbdim      ,&
       & klev       ,&
       & klevp1     ,&
       !
       & lclrsky_lw ,&
       & lclrsky_sw ,&
       !
       & rsdt0      ,&
       & cosmu0     ,&
       & daylght_frc,&
       !
       & emiss      ,&
       & tsr        ,&
       & tsr_rt     ,&
       !
       & rsd_rt     ,&
       & rsu_rt     ,&
       & rld_rt     ,&
       & rlu_rt     ,&
       !
       & rsdcs_rt   ,&
       & rsucs_rt   ,&
       & rldcs_rt   ,&
       & rlucs_rt   ,&
       !
       & rvds_dir_rt,&
       & rpds_dir_rt,&
       & rnds_dir_rt,&
       & rvds_dif_rt,&
       & rpds_dif_rt,&
       & rnds_dif_rt,&
       & rvus_rt    ,&
       & rpus_rt    ,&
       & rnus_rt    ,&
       !
       ! output
       ! ------
       !
       & rsdt       ,&
       & rsut       ,&
       & rsds       ,&
       & rsus       ,&
       !
       & rsutcs     ,&
       & rsdscs     ,&
       & rsuscs     ,&
       !
       & rvds_dir   ,&
       & rpds_dir   ,&
       & rnds_dir   ,&
       & rvds_dif   ,&
       & rpds_dif   ,&
       & rnds_dif   ,&
       & rvus       ,&
       & rpus       ,&
       & rnus       ,&
       !
       & rlut       ,&
       & rlds       ,&
       & rlus       ,&
       !
       & rlutcs     ,&
       & rldscs     ,&
       !
       & q_rsw      ,&
       & q_rlw      )


    INTEGER,  INTENT(in)  :: &
         &  jcs, jce, kbdim, &
         &  klev, klevp1

    LOGICAL,  INTENT(in)  ::        &
         &  lclrsky_lw             ,&! switch for longwave  clearsky computation
         &  lclrsky_sw               ! switch for shortwave clearsky computation

    REAL(wp), INTENT(in)  ::        &
         &  rsdt0                  ,&! indicent SW flux for sun in zenith
         &  cosmu0(:)              ,&! cosine of solar zenith angle at current time
         &  daylght_frc(:)         ,&! daylight fraction; with diurnal cycle 0 or 1, with zonal mean [0,1]
         &  emiss (:)              ,&! lw sfc emissivity
         &  tsr   (:)              ,&! radiative surface temperature at current   time [K]
         &  tsr_rt(:)              ,&! radiative surface temperature at radiation time [K]
         !
         &  rsd_rt(:,:)            ,&! all-sky   shortwave downward flux at radiation time [W/m2]
         &  rsu_rt(:,:)            ,&! all-sky   shortwave upward   flux at radiation time [W/m2]
         !
         &  rsdcs_rt(:,:)          ,&! clear-sky shortwave downward flux at radiation time [W/m2]
         &  rsucs_rt(:,:)          ,&! clear-sky shortwave upward   flux at radiation time [W/m2]
         !
         &  rld_rt(:,:)            ,&! all-sky   longwave  downward flux at radiation time [W/m2]
         &  rlu_rt(:,:)            ,&! all-sky   longwave  upward   flux at radiation time [W/m2]
         !
         &  rldcs_rt(:,:)          ,&! clear-sky longwave  downward flux at radiation time [W/m2]
         &  rlucs_rt(:,:)          ,&! clear-sky longwave  upward   flux at radiation time [W/m2]
         !
         &  rvds_dir_rt(:)         ,&! all-sky   vis. dir. downward flux at radiation time [W/m2]
         &  rpds_dir_rt(:)         ,&! all-sky   par  dir. downward flux at radiation time [W/m2]
         &  rnds_dir_rt(:)         ,&! all-sky   nir  dir. downward flux at radiation time [W/m2]
         &  rvds_dif_rt(:)         ,&! all-sky   vis. dif. downward flux at radiation time [W/m2]
         &  rpds_dif_rt(:)         ,&! all-sky   par  dif. downward flux at radiation time [W/m2]
         &  rnds_dif_rt(:)         ,&! all-sky   nir  dif. downward flux at radiation time [W/m2]
         &  rvus_rt    (:)         ,&! all-sky   visible   upward   flux at radiation time [W/m2]
         &  rpus_rt    (:)         ,&! all-sky   par       upward   flux at radiation time [W/m2]
         &  rnus_rt    (:)           ! all-sky   near-ir   upward   flux at radiation time [W/m2]

    REAL(wp), INTENT(out) ::        &
         &  rsdt  (:)              ,&! all-sky   shortwave downward flux at current   time [W/m2]
         &  rsut  (:)              ,&! all-sky   shortwave upward   flux at current   time [W/m2]
         &  rsds  (:)              ,&! all-sky   shortwave downward flux at current   time [W/m2]
         &  rsus  (:)              ,&! all-sky   shortwave upward   flux at current   time [W/m2]
         !
         &  rsutcs(:)              ,&! clear-sky shortwave upward   flux at current   time [W/m2]
         &  rsdscs(:)              ,&! clear-sky shortwave downward flux at current   time [W/m2]
         &  rsuscs(:)              ,&! clear-sky shortwave upward   flux at current   time [W/m2]
         !
         &  rvds_dir(:)            ,&! all-sky   vis. dir. downward flux at current   time [W/m2]
         &  rpds_dir(:)            ,&! all-sky   par  dir. downward flux at current   time [W/m2]
         &  rnds_dir(:)            ,&! all-sky   nir  dir. downward flux at current   time [W/m2]
         &  rvds_dif(:)            ,&! all-sky   vis. dif. downward flux at current   time [W/m2]
         &  rpds_dif(:)            ,&! all-sky   par  dif. downward flux at current   time [W/m2]
         &  rnds_dif(:)            ,&! all-sky   nir  dif. downward flux at current   time [W/m2]
         &  rvus    (:)            ,&! all-sky   visible   upward   flux at current   time [W/m2]
         &  rpus    (:)            ,&! all-sky   par       upward   flux at current   time [W/m2]
         &  rnus    (:)            ,&! all-sky   near-ir   upward   flux at current   time [W/m2]
         !
         &  rlut  (:)              ,&! all-sky   longwave  upward   flux at current   time [W/m2]
         &  rlds  (:)              ,&! all-sky   longwave  downward flux at current   time [W/m2]
         &  rlus  (:)              ,&! all-sky   longwave  upward   flux at current   time [W/m2]
         !
         &  rlutcs(:)              ,&! clear-sky longwave  upward   flux at current   time [W/m2]
         &  rldscs(:)              ,&! clear-sky longwave  downward flux at current   time [W/m2]
         !
         &  q_rsw (:,:)            ,&! radiative shortwave heating  [W/m2]
         &  q_rlw (:,:)              ! radiative longwave  heating  [W/m2]

    ! Local arrays
    REAL(wp) ::                     &
         &  xsdt  (kbdim)          ,&
         &  rsn   (kbdim,klevp1)   ,&
         &  rln   (kbdim,klevp1)

    REAL(wp) :: drlus_dtsr, dtsr
    INTEGER  :: jc, jk

    !$ACC DATA PRESENT(cosmu0, daylght_frc, emiss, tsr, tsr_rt, rsd_rt, rsu_rt, rsdcs_rt, rsucs_rt) &
    !$ACC   PRESENT(rld_rt, rlu_rt, rldcs_rt, rlucs_rt, rvds_dir_rt, rpds_dir_rt, rnds_dir_rt) &
    !$ACC   PRESENT(rvds_dif_rt, rpds_dif_rt, rnds_dif_rt, rvus_rt, rpus_rt, rnus_rt, rsdt, rsut) &
    !$ACC   PRESENT(rsds, rsus, rsutcs, rsdscs, rsuscs, rvds_dir, rpds_dir, rnds_dir, rvds_dif) &
    !$ACC   PRESENT(rpds_dif, rnds_dif, rvus, rpus, rnus, rlut, rlds, rlus, rlutcs, rldscs) &
    !$ACC   PRESENT(q_rsw, q_rlw) &
    !$ACC   CREATE(xsdt, rsn, rln)

    ! Shortwave fluxes
    ! ----------------
    !
    ! The original downward and upward fluxes form the radiative transfer (rt) calculation
    ! are scaled by the ratio of the incident solar fluxes of the current time and the
    ! rt-time. This assumes that the incident solar flux at rt-time is non-zero in all
    ! columns, where cosmu0 is currently positive.
    ! - incident solar radiation at rt-time     : rsdt_rt = rsd_rt(jk=1)
    ! - incident solar radiation at current time: rsdt    = rsdt0*MAX(0,cosmu0)
    ! - scaling ratio for fluxes at current time: xsdt    = rsdt / rsdt_rt
    !
    ! top of atmophere
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(jcs, jce, rsdt0, lclrsky_sw) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs, jce
      rsdt  (jc)   = rsdt0*cosmu0(jc)*daylght_frc(jc)
      IF (rsd_rt(jc,1) > 0.0_wp) THEN
         xsdt  (jc)   = rsdt  (jc) / rsd_rt(jc,1)
      ELSE
         xsdt  (jc)   = 0.0_wp
      END IF
      !
      ! all sky
      rsut  (jc)   = rsu_rt  (jc,1) * xsdt(jc)
      ! clear sky
      IF (lclrsky_sw) THEN
         rsutcs(jc)   = rsucs_rt(jc,1) * xsdt(jc)
      END IF
      !
    END DO
    !$ACC END PARALLEL
    !
    ! all half levels
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(klevp1, jcs, jce) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klevp1
      DO jc = jcs, jce
        rsn(jc,jk) = (rsd_rt(jc,jk) - rsu_rt(jc,jk)) * xsdt(jc)
      END DO
    END DO
    !$ACC END PARALLEL
    !
    ! surface
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(jcs, jce, klevp1, lclrsky_sw, lclrsky_lw) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs, jce
      ! all sky
      rsds  (jc)   = rsd_rt  (jc,klevp1) * xsdt(jc)
      rsus  (jc)   = rsu_rt  (jc,klevp1) * xsdt(jc)
      ! clear sky
      IF (lclrsky_sw) THEN
         rsdscs(jc)   = rsdcs_rt(jc,klevp1) * xsdt(jc)
         rsuscs(jc)   = rsucs_rt(jc,klevp1) * xsdt(jc)
      END IF
      !
      ! components
      rvds_dir(jc) = rvds_dir_rt(jc) * xsdt(jc)
      rpds_dir(jc) = rpds_dir_rt(jc) * xsdt(jc)
      rnds_dir(jc) = rnds_dir_rt(jc) * xsdt(jc)
      !
      rvds_dif(jc) = rvds_dif_rt(jc) * xsdt(jc)
      rpds_dif(jc) = rpds_dif_rt(jc) * xsdt(jc)
      rnds_dif(jc) = rnds_dif_rt(jc) * xsdt(jc)
      !
      rvus(jc)     = rvus_rt(jc) * xsdt(jc)
      rpus(jc)     = rpus_rt(jc) * xsdt(jc)
      rnus(jc)     = rnus_rt(jc) * xsdt(jc)

      ! Longwave fluxes
      ! ---------------
      !
      ! The original downward and upward fluxes form the radiative transfer (rt) calculation
      ! are kept constant, except for the upward flux from the surface, which is corrected
      ! for the change in radiative surface temperature using a 1st order Taylor expansion.
      ! - surface upward flux at rt-time      : rlus_rt = rlu(jk=klevp1)
      ! - rad. surface temp.  at rt-time      : tsr_rt
      ! - rad. surface temp.  at current time : tsr
      !
      ! top of atmophere
      ! all sky
      rlut  (jc)   = rlu_rt  (jc,1)
      ! clear sky
      IF (lclrsky_lw) THEN
         rlutcs(jc)   = rlucs_rt(jc,1)
      END IF
      !
    END DO
    !$ACC END PARALLEL
    ! all half levels
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(klevp1, jcs, jce) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klevp1
      DO jc = jcs, jce
        rln(jc,jk) = (rld_rt(jc,jk) - rlu_rt(jc,jk))
      END DO
    END DO
    !$ACC END PARALLEL
    !
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(klevp1, jcs, jce, lclrsky_lw) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(drlus_dtsr, dtsr)
    DO jc = jcs, jce
      ! surface
      ! all sky
      rlds  (jc)  = rld_rt  (jc,klevp1)
      ! clear sky
      IF (lclrsky_lw) THEN
         rldscs(jc)  = rldcs_rt(jc,klevp1)
      END IF
      !
      ! - correct upward flux for changed radiative surface temperature
      drlus_dtsr = emiss(jc)*4._wp*stbo*tsr(jc)**3 ! derivative
      dtsr       = tsr(jc) - tsr_rt(jc)            ! change in tsr
      rlus(jc)   = rlu_rt(jc,klevp1)             & ! rlus = rlus_rt
           &          +drlus_dtsr * dtsr           !       + correction
      !
      rln(jc,klevp1) = rlds(jc) - rlus(jc)
    END DO
    !$ACC END PARALLEL


    ! Heating rates in atmosphere
    !----------------------------
    !$ACC PARALLEL DEFAULT(NONE) FIRSTPRIVATE(klev, jcs, jce) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klev
      DO jc = jcs, jce
        q_rsw(jc,jk) = rsn(jc,jk)-rsn(jc,jk+1)
        q_rlw(jc,jk) = rln(jc,jk)-rln(jc,jk+1)
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE radheating

END MODULE mo_radheating
