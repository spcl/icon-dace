!
! Computation of cloud cover and grid mean cloud liquid water and cloud ice
!
! This routine takes information from turbulence, convection and grid-scale
! to produce cloud properties used in radiation (and microphysics).
!
! Possible future options
! - simple diagnostic (from turbulence, convection and grid scale)
! - prognostic total water variance AND prognostic ice
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
#include "consistent_fma.inc"
!----------------------------

MODULE mo_art_cover_koe

  USE mo_kind,               ONLY: wp, i4

  USE mo_physical_constants, ONLY: rv         !! Rv

  USE mo_satad,              ONLY: sat_pres_ice   ! saturation pressure over ice

  USE mo_fortran_tools,      ONLY: set_acc_host_or_device, assert_acc_host_only


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_cover_dusty

  CHARACTER(len=*), PARAMETER :: modname = 'mo_art_cover_koe'

CONTAINS

!-------------------------------------------------------------------------
!
!  Parameterization of dusty cirrus cloud fraction
!  based on Seifert et al. 2023, ACP
!
!-------------------------------------------------------------------------

SUBROUTINE art_cover_dusty( &
  & kidia, kfdia, klon, kstart, klev, & ! in:    dimensions (turn off physics above kstart)
  & tt                              , & ! in:    temperature
  & deltaz                          , & ! in:    layer thickness
  & rho                             , & ! in:    density
  & qv                              , & ! inout: water vapor
  & dusta, dustb, dustc             , & ! in:    prognostic dust variables
  & dustyci_crit, dustyci_rhi       , & ! in:    thresholds
  & lacc                            , & ! in:    parameter to prevent openacc during init
  & cc_tot, qi_tot                    ) ! inout: cloud output diagnostic

!! Subroutine arguments:
!! --------------------

INTEGER(KIND=i4), INTENT(IN) ::  &
     & kidia            , & ! horizontal start index
     & kfdia            , & ! horizontal end   index
     & klon             , & ! horizontal dimension
     & kstart           , & ! vertical start index (turn off physics above)
     & klev                 ! vertical dimension

REAL(KIND=wp), INTENT(IN) ::  &
     & tt(:,:)          , & ! temperature                                   (  K  )
     & deltaz(:,:)      , & ! layer thickness                               (m)
     & rho(:,:)         , & ! density                                       (kg/m3)
     & qv(:,:)          , & ! specific water vapor content                  (kg/kg)
     & dusta(:,:)       , & ! dust mass concentration, mode A
     & dustb(:,:)       , & ! dust mass concentration, mode B
     & dustc(:,:)           ! dust mass concentration, mode C

REAL(KIND=wp), INTENT(IN) ::  &
     & dustyci_crit     , & ! dust threshold for dusty cirrus               (mug/kg)
     & dustyci_rhi          ! rhi threshold for dusty cirrus                (  -  )

LOGICAL, OPTIONAL, INTENT(IN) :: lacc ! parameter to prevent openacc during init

REAL(KIND=wp), INTENT(INOUT) ::   &
     & cc_tot(:,:)      , & ! cloud cover diagnostic
     & qi_tot(:,:)          ! specific cloud ice   content diagnostic       (kg/kg)

!! Local variables:
!! ----------------

LOGICAL ::  &
     & lzacc

INTEGER (KIND=i4) :: &
     & jl, jk, jkp1, jkp2, jkp3, jkp4, jkp5

REAL(KIND=wp) :: &
     & vap_pres

REAL(KIND=wp) :: &
     & zdust(klon,klev), zrhi(klon,klev), zdTdz(klon,klev), &
     & qi_dusty(klon,klev)

REAL(KIND=wp) :: mdust, mrhi, mgrad

!! Local parameters:
!! -----------------

CHARACTER(len=*), PARAMETER :: routine = modname//':art_cover_dusty'

REAL(KIND=wp), PARAMETER  :: &
     & dustyci_mdust = 100.0_wp,   & ! dust threshold for thick cirrus  [mug/kg]
     & dustyci_temp  = 240.0_wp,   & ! temp threshold for dusty cirrus  [K]
     & dustyci_lapse = -6.5e-3_wp, & ! max lapse rate for dusty cirrus  [K/m]
     & dustyci_mrhi  = 1.00_wp,    & ! rhi threshold for thick cirrus   [-]
     & dustyci_iwc   =  80e-6_wp     ! fixed IWC of dusty cirrus        [kg/m3]

!-----------------------------------------------------------------------

! JF:   CALL set_acc_host_or_device(lzacc, lacc)
  CALL assert_acc_host_only(routine, lacc)

!$ACC DATA CREATE(zdust, zrhi, zdTdz, qi_dusty) &
!$ACC   IF(lzacc)

!-----------------------------------------------------------------------
! Calculate water vapour saturation mixing ratio of over ice
! and prepare for dusty cirrus
!-----------------------------------------------------------------------

!$ACC PARALLEL IF(lzacc) DEFAULT(PRESENT) ASYNC(1)
!$ACC LOOP GANG VECTOR PRIVATE(vap_pres)
DO jk = kstart, klev
  jkp1 = MIN(jk+1,klev)
  DO jl = kidia, kfdia

    vap_pres = qv(jl,jk) * rho(jl,jk) * rv * tt(jl,jk)

    ! sum up dust with more weight on large modes
    IF (tt(jl,jk) < dustyci_temp) THEN
      zdust(jl,jk) = 1.0_wp * dustb(jl,jk) + 2.0_wp * dustc(jl,jk)
      zrhi(jl,jk)  = vap_pres / sat_pres_ice(tt(jl,jk))
      zdtdz(jl,jk) = (tt(jl,jk)-tt(jl,jkp1))/deltaz(jl,jk)
    ELSE
      zrhi(jl,jk)  = 0.0_wp
      zdust(jl,jk) = 0.0_wp
      zdtdz(jl,jk) = 0.0_wp
    END IF

  ENDDO
ENDDO
!$ACC END PARALLEL

!-----------------------------------------------------------------------
! Add dusty cirrus to cloud fraction
!-----------------------------------------------------------------------

!$ACC PARALLEL IF(lzacc) DEFAULT(PRESENT) ASYNC(1)
!$ACC LOOP GANG
DO jk = kstart,klev
  jkp1 = MIN(jk+1,klev)
  jkp2 = MIN(jk+2,klev)
  jkp3 = MIN(jk+3,klev)
  jkp4 = MIN(jk+4,klev)
  jkp5 = MIN(jk+5,klev)

  !$ACC LOOP VECTOR PRIVATE(mdust, mrhi, mgrad)
  DO jl = kidia,kfdia

    ! parameterization of dusty cirrus above a dust layer for ICON-ART
    IF (tt(jl,jk) < dustyci_temp) THEN

      mdust = MAX(zdust(jl,jkp1),zdust(jl,jkp2),zdust(jl,jkp3),zdust(jl,jkp4),zdust(jl,jkp5))
      mrhi  = MAX(zrhi(jl,jk),zrhi(jl,jkp1),zrhi(jl,jkp2),zrhi(jl,jkp3),zrhi(jl,jkp4))
      mgrad = MIN(zdtdz(jl,jk),zdtdz(jl,jkp1))

      IF ( mdust > dustyci_crit .AND.  mrhi  > dustyci_rhi .AND. mgrad < dustyci_lapse ) THEN
        qi_dusty(jl,jk)  = dustyci_iwc/rho(jl,jk)  &
          &              * MAX(MIN(1.0_wp,((mdust-dustyci_crit)/(dustyci_mdust-dustyci_crit))),0.1_wp) &
          &              * MAX(MIN(1.0_wp,((mrhi -dustyci_rhi )/(dustyci_mrhi -dustyci_rhi ))),0.1_wp)
        qi_tot(jl,jk)    = MAX( qi_tot(jl,jk), qi_dusty(jl,jk) )
        cc_tot(jl,jk)    = 1.0_wp
      END IF

    END IF
  ENDDO
ENDDO
!$ACC END PARALLEL

!$ACC WAIT
!$ACC END DATA

END SUBROUTINE art_cover_dusty

END MODULE mo_art_cover_koe
