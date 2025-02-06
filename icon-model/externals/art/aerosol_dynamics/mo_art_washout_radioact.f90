!
! mo_art_washout_radioact
! This module provides the washout routine of radioactive particles.
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

MODULE mo_art_washout_radioact
! ICON
  USE mo_kind,                          ONLY: wp
! ART
  USE mo_art_data,                      ONLY: t_art_data 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_washout_radioact

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_washout_radioact(rho, dz, qr, qs,                                     &
  &                             rr_gsp, rr_con, rr_con3d, ss_gsp, ss_con, ss_con3d,  &
  &                             itracer, wdf, wde, istart, iend, kstart, nlev, jb,   &
  &                             dtime, tracer, p_art_data)
!<
! SUBROUTINE art_washout_radioact
! Calculates washout for radioactive particles
! Based on: -
! Part of Module: mo_art_washout_radioact
! Author: Daniel Rieger, KIT
! Initial Release: 2013-02-28
! Modifications:
! 2016-11-07: Daniel Rieger, KIT
! - Adaption to ART physics coupling concept, cleanup
!>
  REAL(wp), INTENT(IN)              :: &
    &  rho(:,:),                       & !< Density (kg m-3)
    &  dz(:,:),                        & !< Vertical layer thickness (m)
    &  qr(:,:),                        & !< Mass mixing ratio of rain (kg kg-1)
    &  qs(:,:),                        & !< Mass mixing ratio of snow (kg kg-1)
    &  rr_gsp(:),                      & !< grid scale prec. rate at sfc
    &  rr_con(:),                      & !< conv. prec. rate at sfc
    &  rr_con3d(:,:),                  & !< 3D conv. prec. rate
    &  ss_gsp(:),                      & !< grid scale snow rate at sfc
    &  ss_con(:),                      & !< conv. snow rate at sfc
    &  ss_con3d(:,:),                  & !< 3D conv. snow prec. rate
    &  wdf, wde,                       & !< factor and exponent for wet deposition
    &  dtime                             !< model timestep (s)
  INTEGER, INTENT(IN)               :: &
    &  itracer,                        & !< Index in tracer container
    &  istart, iend,                   & !< Start and end index of jc loop
    &  kstart, nlev,                   & !< Start level and number of vertical levels
    &  jb                                !< Block index
  REAL(wp), INTENT(INOUT)           :: &
    &  tracer(:,:)                       !< tracer mixing ratio (Bq kg-1)
  TYPE(t_art_data),INTENT(INOUT)    :: &
    &  p_art_data                        !< ART data container
! Local variables
  REAL(wp)                          :: &
    &  tracer_before,                  & !< tracer concentration before washout
    &  washout_rate,                   & !< rate of radioactive washout
    &  factor_a,                       & !< factor for capturing the different deposition efficiency
    &  factor_b                          !< of liquid and solid precipitation
  INTEGER                           :: &
    &  jc, jk                            !<loop indizes

  factor_a = 1._wp
  factor_b = 5._wp


  !-----------------------------------------------------------------------------------------
  !--   Get the nuclide specific factors  --------------------------------------------------
  !-----------------------------------------------------------------------------------------

  DO jk=kstart,nlev
    DO jc=istart,iend

      tracer_before = tracer(jc,jk)

      ! Washout by rain
      IF ( qr(jc,jk)> 0.0_wp .OR. rr_con3d(jc,jk)> 0.0_wp) THEN
        IF ((rr_con(jc) + rr_gsp(jc)) > 0.0_wp) THEN
          washout_rate  = factor_a * wdf * ((rr_gsp(jc)+rr_con(jc)) ** wde)
          tracer(jc,jk) = tracer(jc,jk) - (tracer(jc,jk) * washout_rate * dtime)
        ENDIF
      ENDIF

      ! Washout by snow
      IF ( qs(jc,jk)> 0.0_wp .OR. ss_con3d(jc,jk)> 0.0_wp) THEN
        IF ((ss_con(jc) + ss_gsp(jc)) > 0.0_wp) THEN
          washout_rate  = factor_b * wdf * ((ss_gsp(jc)+ss_con(jc)) ** wde)
          tracer(jc,jk) = tracer(jc,jk) - (tracer(jc,jk) * washout_rate * dtime)
        ENDIF
      ENDIF

      ! Accumulated wet deposition [Bq/m2]
      p_art_data%diag%radioact(itracer)%wetdepo(jc,jb) =       &
        &    p_art_data%diag%radioact(itracer)%wetdepo(jc,jb)  &
        &  + ( tracer_before - tracer(jc,jk) )                   &
        &    * rho(jc,jk) * dz(jc,jk)

    ENDDO! jc
  ENDDO! jk

END SUBROUTINE art_washout_radioact
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_washout_radioact
