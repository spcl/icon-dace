!
! mo_art_cond_aero
! Calculation of sulfate condensation
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

MODULE mo_art_cond_aero
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_physical_constants,            ONLY: argas
! ART

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_cond_aero'

  PUBLIC :: art_prepare_cond_so4
  PUBLIC :: art_finalize_cond_so4

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_cond_so4(temp, pres, rho, ch2so4, nconc, dg0, sigmag, & 
  &                             istart, iend, kstart, kend, dmdt)
!<
! SUBROUTINE art_prepare_cond_so4
! Calculation of sulfate condensation
! Based on: Riemer, 2002: Numerische Simulationen zur Wirkung des Aerosols auf die 
!                         troposphaerische Chemie und die Sichtweite
! Part of Module: mo_art_cond_aero
! Author: Daniel Rieger, KIT
! Initial Release: 2017-05-11
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)          :: &
    &  temp(:,:), pres(:,:),       & !< Temperature (K) and pressure (Pa)
    &  rho(:,:),                   & !< Air density (kg m-3)
    &  ch2so4(:,:),                & !< Current H2SO4 concentration (usually: mug kg-1)
    &  nconc(:,:),                 & !< Number concentration of tracer (kg-1)
    &  dg0(:,:),                   & !< Count median diameter of current mode
    &  sigmag                        !< Standard deviation
  INTEGER, INTENT(in)           :: &
    &  istart, iend,               & !< jc-loop start and end indices
    &  kstart, kend                  !< jk-loop start and end indices
  REAL(wp), INTENT(inout)       :: & !< inout as we might want to pass a pointer
    &  dmdt(:,:)                     !< Change in mass mixing ratio (same unit as ch2so4)
! Local variables
  REAL(wp)                      :: &
    &  ccofm,                      & !< Factor for free molecular condensation
    &  csqt,                       & !< Prefactor in eq. 3.47 Riemer Diss
    &  mom1,                       & !< First moment
    &  facmom1,                    & !< Conversion factor for 0th to 1st moment
    &  mom2,                       & !< Second moment
    &  facmom2,                    & !< Conversion factor for 0th to 2nd moment
    &  diffcorr,                   & !< diffsulf corrected for temperature and pressure
    &  int_nc,                     & !< Integral from 3.42 in Riemer Diss for near continuum (kg-1 s-1)
    &  int_fm,                     & !< Integral from 3.42 in Riemer Diss for free molecular (kg-1 s-1)
    &  int_harm,                   & !< Harmonic mean of int_nc and int_fm
    &  kg2mug                        !< Factor kg/kg to amug/kg
  REAL(wp), PARAMETER           :: &
    &  alphsulf = 1.0_wp,          & !< 
    &  diffsulf = 9.362223E-06_wp, & !< 
    &  mwh2so4  = 0.09807354_wp,   & !< Molecular weight of H2SO4 (kg mol-1)
    &  rhoso4   = 1.8E+3_wp,       & !< Bulk density of aerosol sulfate (kg m-3)
    &  cconc = 2._wp*pi * diffsulf   !< Prefactor
  INTEGER                       :: &
    &  jc, jk                        !< Loop counter


  kg2mug     = 1.E+09_wp
  ccofm      = alphsulf * SQRT( pi * argas / (2._wp * mwh2so4))
  dmdt(:,:)  = 0._wp

  facmom1  = EXP( 0.125_wp * (LOG(sigmag)**2) )**4
  facmom2  = EXP( 0.125_wp * (LOG(sigmag)**2) )**16

  DO jk = kstart, kend
!NEC$ ivdep
    DO jc = istart, iend
      IF (dg0(jc,jk) > 0._wp) THEN
        mom1        = nconc(jc,jk) * rho(jc,jk) * dg0(jc,jk) * facmom1
        diffcorr    = (101325._wp / pres(jc,jk)) * (temp(jc,jk) / 273.16_wp)**1.75_wp
        int_nc      = cconc * mom1 * diffcorr ! I_{l,nc} eq. 3.46

        mom2        = nconc(jc,jk) * rho(jc,jk) * dg0(jc,jk) * dg0(jc,jk) * facmom2
        csqt        = ccofm * SQRT(temp(jc,jk))
        int_fm      = csqt * mom2             ! I_{l,fm} eq. 3.47

        int_harm    = int_nc * int_fm / ( int_nc + int_fm )
        dmdt(jc,jk) = pi / 6._wp * ch2so4(jc,jk) * kg2mug * int_harm
      ENDIF
    ENDDO !jc
  ENDDO !jk

END SUBROUTINE art_prepare_cond_so4
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_finalize_cond_so4(dmdt_condso4, ch2so4, totmasscond, so4tracer, &
  &                              p_dtime, istart, iend, kstart, kend)
!<
! SUBROUTINE art_finalize_cond_so4
! Update of sulfate condensation
! Based on: Riemer, 2002: Numerische Simulationen zur Wirkung des Aerosols auf die 
!                         troposphaerische Chemie und die Sichtweite
! Part of Module: mo_art_cond_aero
! Author: Daniel Rieger, KIT
! Initial Release: 2017-05-15
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)    :: &
    &  dmdt_condso4(:,:),    & !< Change in sulfate mass mixing ratio (mug kg-1 s-1)
!   &  ch2so4(:,:),          & !< Sulfuric acid mass mixing ratio (mug kg-1)
    &  totmasscond(:,:),     & !< Total condensated mass (mug kg-1)
    &  p_dtime                 !< Model time step
  REAL(wp), INTENT(inout) :: &
    &  ch2so4(:,:),          & !< Sulfuric acid mass mixing ratio (mug kg-1)
    &  so4tracer(:,:)          !< Sulfate tracer mass mixing ratio (mug kg-1)
  INTEGER, INTENT(in)     :: &
    &  istart, iend,         & !< jc-loop start and end indice
    &  kstart, kend            !< jk-loop start and end indice
! Local variables
  REAL(wp)                :: &
    &  update,               & !< Update of concentration (mug kg-1)
    &  ratio,                & !< Ratio of total produced so4 mass to available h2so4 mass
    &  kg2mug                  !< Factor kg/kg to amug/kg
  INTEGER                 :: &
    &  jc, jk                  !< Loop counter

  kg2mug = 1.E+09_wp

  DO jk = kstart, kend
!NEC$ ivdep
    DO jc = istart, iend
      IF(totmasscond(jc,jk) > 0._wp .AND. ch2so4(jc,jk) > 0._wp) THEN

        ch2so4(jc,jk) = ch2so4(jc,jk) *  kg2mug

        ratio = totmasscond(jc,jk) / ch2so4(jc,jk)

        IF (ratio > 1._wp) THEN
          update = dmdt_condso4(jc,jk) / ratio * p_dtime
        ELSE
          update = dmdt_condso4(jc,jk) * p_dtime
        ENDIF

        so4tracer(jc,jk) = so4tracer(jc,jk) + update
!       so4tracer(jc,jk) = max(so4tracer(jc,jk),1.e-5)

        ch2so4(jc,jk) = ch2so4(jc,jk) / kg2mug
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE art_finalize_cond_so4
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_cond_aero
