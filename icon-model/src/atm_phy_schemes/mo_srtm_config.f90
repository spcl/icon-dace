!----------------------------------------------------------------------------
! Copyright (c) 2002-2020, Atmospheric & Environmental Research, Inc. (AER)
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!  * Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!  * Neither the name of Atmospheric & Environmental Research, Inc., nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL ATMOSPHERIC & ENVIRONMENTAL RESEARCH, INC.,
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
! THE POSSIBILITY OF SUCH DAMAGE.
!                        (http://www.rtweb.aer.com/)
!----------------------------------------------------------------------------

! SPDX-License-Identifier: BSD-3-Clause
! Code has been modified for the use in ICON (and formerly in ECHAM)

#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
#endif
  MODULE mo_srtm_config

    USE mo_kind,         ONLY: wp
    USE mo_srtm_kgb16,   ONLY: srtm_kgb16
    USE mo_srtm_kgb17,   ONLY: srtm_kgb17
    USE mo_srtm_kgb18,   ONLY: srtm_kgb18
    USE mo_srtm_kgb19,   ONLY: srtm_kgb19
    USE mo_srtm_kgb20,   ONLY: srtm_kgb20
    USE mo_srtm_kgb21,   ONLY: srtm_kgb21
    USE mo_srtm_kgb22,   ONLY: srtm_kgb22
    USE mo_srtm_kgb23,   ONLY: srtm_kgb23
    USE mo_srtm_kgb24,   ONLY: srtm_kgb24
    USE mo_srtm_kgb25,   ONLY: srtm_kgb25
    USE mo_srtm_kgb26,   ONLY: srtm_kgb26
    USE mo_srtm_kgb27,   ONLY: srtm_kgb27
    USE mo_srtm_kgb28,   ONLY: srtm_kgb28
    USE mo_srtm_kgb29,   ONLY: srtm_kgb29

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: setup_srtm, jpg, jpinpx, jpb1, jpb2, jpgpt, jpsw, ng16, ng17,&
      &       ng18, ng19, ng20, ng21, ng22, ng23, ng24, ng25, ng26, ng27,  &
      &       ng28, ng29, nspa, nspb, ngc, preflog, tref, repclc, replog,  &
      &       wavenum1, wavenum2, delwave, ssi_default, ssi_amip, ssi_rce, &
      &       ssi_preind, ssi_coddington


    SAVE

    !     ------------------------------------------------------------------
    !     Parameters relevant to AER's RRTM-SW radiation scheme

    !     030224  JJMorcrette

    !     Modified for g-point reduction from 224 to 112.
    !     Swap code below to restore 224 g-point set.
    !     Mar2004 MJIacono, AER
    !     ------------------------------------------------------------------

    INTEGER, PARAMETER :: jpg    = 16
    INTEGER, PARAMETER :: jpinpx = 35
    INTEGER, PARAMETER :: jpband = 29
    INTEGER, PARAMETER :: jpsw   = 14
    INTEGER, PARAMETER :: jpb1   = 16
    INTEGER, PARAMETER :: jpb2   = 29
    INTEGER, PARAMETER :: jpgpt  = 112 ! number of new g-points

    INTEGER, PARAMETER :: jmcmu  = 32
    INTEGER, PARAMETER :: jmumu  = 32
    INTEGER, PARAMETER :: jmphi  = 3
    INTEGER, PARAMETER :: jmxang = 4
    INTEGER, PARAMETER :: jmxstr = 16

    INTEGER, PARAMETER :: ng16 = 6
    INTEGER, PARAMETER :: ng17 = 12
    INTEGER, PARAMETER :: ng18 = 8
    INTEGER, PARAMETER :: ng19 = 8
    INTEGER, PARAMETER :: ng20 = 10
    INTEGER, PARAMETER :: ng21 = 10
    INTEGER, PARAMETER :: ng22 = 2
    INTEGER, PARAMETER :: ng23 = 10
    INTEGER, PARAMETER :: ng24 = 8
    INTEGER, PARAMETER :: ng25 = 6
    INTEGER, PARAMETER :: ng26 = 6
    INTEGER, PARAMETER :: ng27 = 8
    INTEGER, PARAMETER :: ng28 = 6
    INTEGER, PARAMETER :: ng29 = 12

    INTEGER, PARAMETER :: ngs16 = 0
    INTEGER, PARAMETER :: ngs17 = 16
    INTEGER, PARAMETER :: ngs18 = 32
    INTEGER, PARAMETER :: ngs19 = 48
    INTEGER, PARAMETER :: ngs20 = 64
    INTEGER, PARAMETER :: ngs21 = 80
    INTEGER, PARAMETER :: ngs22 = 96
    INTEGER, PARAMETER :: ngs23 = 112
    INTEGER, PARAMETER :: ngs24 = 128
    INTEGER, PARAMETER :: ngs25 = 144
    INTEGER, PARAMETER :: ngs26 = 160
    INTEGER, PARAMETER :: ngs27 = 176
    INTEGER, PARAMETER :: ngs28 = 192
    INTEGER, PARAMETER :: ngs29 = 208


    !=====================================================================
    ! Set arrays needed for the g-point reduction from 224 to 112 for the
    ! 14 SW bands:
    ! This mapping from 224 to 112 points has been carefully selected to
    ! minimize the effect on the resulting fluxes and cooling rates, and
    ! caution should be used if the mapping is modified.

    ! JPGPT   The total number of new g-points (NGPT)
    ! NGC
    ! NGS     The cumulative sum of new g-points for each band
    ! NGM     The index of each new g-point relative to the original
    !         16 g-points for each band.
    ! NGN     The number of original g-points that are combined to make
    !         each new g-point in each band.
    ! NGB     The band index for each new g-point.
    ! WT      RRTM weights for 16 g-points.


    INTEGER,  PARAMETER :: &
      &    ng(16:29) = (/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16 /), &
      &  nspa(16:29) = (/  9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1 /), &
      &  nspb(16:29) = (/  1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1 /), &
      &  nmpsrtm(14) = (/  6, 6, 5, 5, 5, 5, 5, 4, 4, 3, 2, 2, 1, 6 /)

    INTEGER, PARAMETER :: ngc(14) = (/& !< Number of new g-points in each band
      &  6, 12,  8,  8, 10, 10,  2, 10,  8,  6,  6,  8,  6, 12/)

    INTEGER, PARAMETER :: ngs(14) = (/& !< Cumulative sum of g-points in bands
      &  6, 18, 26, 34, 44, 54, 56, 66, 74, 80, 86, 94,100,112/)

    INTEGER, PARAMETER :: ngm(224) =  (/ & !< index of g-points relative to orig
      & 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, & ! Band 16
      & 1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9,10,10,11,12,12, & ! Band 17
      & 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, & ! Band 18
      & 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, & ! Band 19
      & 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10,10,10, & ! Band 20
      & 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10,10,10, & ! Band 21
      & 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2 ,2, & ! Band 22
      & 1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9,10,10,10,10, & ! Band 23
      & 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, & ! Band 24
      & 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! Band 25
      & 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! Band 26
      & 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, & ! Band 27
      & 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! Band 28
      & 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9,10,11,12 /) ! Band 29

    INTEGER, PARAMETER :: ngn(112) = (/         & !< number of combined g points
      & 2, 2, 2, 2, 4, 4,                    &            ! Band 16
      & 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2,  &            ! Band 17
      & 1, 1, 1, 1, 2, 2, 4, 4,              &            ! Band 18
      & 1, 1, 1, 1, 2, 2, 4, 4,              &            ! Band 19
      & 1, 1, 1, 1, 1, 1, 1, 1, 2, 6,        &            ! Band 20
      & 1, 1, 1, 1, 1, 1, 1, 1, 2, 6,        &            ! Band 21
      & 8, 8,                                &            ! Band 22
      & 2, 2, 1, 1, 1, 1, 1, 1, 2, 4,        &            ! Band 23
      & 2, 2, 2, 2, 2, 2, 2, 2,              &            ! Band 24
      & 1, 1, 2, 2, 4, 6,                    &            ! Band 25
      & 1, 1, 2, 2, 4, 6,                    &            ! Band 26
      & 1, 1, 1, 1, 1, 1, 4, 6,              &            ! Band 27
      & 1, 1, 2, 2, 4, 6,                    &            ! Band 28
      & 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1  /)            ! Band 29

    INTEGER, PARAMETER :: ngbsw(112)=(/         & !< band index for new points
      & 16,16,16,16,16,16,                   &            ! Band 16
      & 17,17,17,17,17,17,17,17,17,17,17,17, &            ! Band 17
      & 18,18,18,18,18,18,18,18,             &            ! Band 18
      & 19,19,19,19,19,19,19,19,             &            ! Band 19
      & 20,20,20,20,20,20,20,20,20,20,       &            ! Band 20
      & 21,21,21,21,21,21,21,21,21,21,       &            ! Band 21
      & 22,22,                               &            ! Band 22
      & 23,23,23,23,23,23,23,23,23,23,       &            ! Band 23
      & 24,24,24,24,24,24,24,24,             &            ! Band 24
      & 25,25,25,25,25,25,                   &            ! Band 25
      & 26,26,26,26,26,26,                   &            ! Band 26
      & 27,27,27,27,27,27,27,27,             &            ! Band 27
      & 28,28,28,28,28,28,                   &            ! Band 28
      & 29,29,29,29,29,29,29,29,29,29,29,29 /)            ! Band 29

    REAL(wp), PARAMETER :: repclc=1.e-12_wp ! epsilon for something
    REAL(wp), PARAMETER :: replog=1.e-12_wp ! epsilon for lapace transformation

    REAL(wp), PARAMETER :: wavenum1(16:29) = (/ & !< lower wave number limit
      &     2600._wp, 3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, &
      &     7700._wp, 8050._wp,12850._wp,16000._wp,22650._wp,29000._wp, &
      &    38000._wp, 820._wp  /)
    REAL(wp), PARAMETER :: wavenum2(16:29) = (/ & !< upper wave number limit
      &     3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, 7700._wp, &
      &     8050._wp,12850._wp,16000._wp,22650._wp,29000._wp,38000._wp, &
      &    50000._wp, 2600._wp /)

    REAL(wp), PARAMETER :: delwave(16:29) = wavenum2 - wavenum1 !< width of band

    REAL(wp), PARAMETER :: wt(16) =  (/ & !< rrtm weights for 16 g-points
      & 0.1527534276_wp, 0.1491729617_wp, 0.1420961469_wp, &
      & 0.1316886544_wp, 0.1181945205_wp, 0.1019300893_wp, &
      & 0.0832767040_wp, 0.0626720116_wp, 0.0424925000_wp, &
      & 0.0046269894_wp, 0.0038279891_wp, 0.0030260086_wp, &
      & 0.0022199750_wp, 0.0014140010_wp, 0.0005330000_wp, &
      & 0.0000750000_wp /)

    REAL(wp) ::  wtsm(16),  rwgt(224)

    !=======================================================================
    ! These pressures are chosen such that the ln of the first pressure
    ! has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
    !  each subsequent ln(pressure) differs from the previous one by 0.2.

    REAL(wp), PARAMETER :: pref(59) = (/   & !< reference pressure
      & 1.05363e+03_wp, 8.62642e+02_wp, 7.06272e+02_wp, 5.78246e+02_wp, &
      & 4.73428e+02_wp, 3.87610e+02_wp, 3.17348e+02_wp, 2.59823e+02_wp, &
      & 2.12725e+02_wp, 1.74164e+02_wp, 1.42594e+02_wp, 1.16746e+02_wp, &
      & 9.55835e+01_wp, 7.82571e+01_wp, 6.40715e+01_wp, 5.24573e+01_wp, &
      & 4.29484e+01_wp, 3.51632e+01_wp, 2.87892e+01_wp, 2.35706e+01_wp, &
      & 1.92980e+01_wp, 1.57998e+01_wp, 1.29358e+01_wp, 1.05910e+01_wp, &
      & 8.67114e+00_wp, 7.09933e+00_wp, 5.81244e+00_wp, 4.75882e+00_wp, &
      & 3.89619e+00_wp, 3.18993e+00_wp, 2.61170e+00_wp, 2.13828e+00_wp, &
      & 1.75067e+00_wp, 1.43333e+00_wp, 1.17351e+00_wp, 9.60789e-01_wp, &
      & 7.86628e-01_wp, 6.44036e-01_wp, 5.27292e-01_wp, 4.31710e-01_wp, &
      & 3.53455e-01_wp, 2.89384e-01_wp, 2.36928e-01_wp, 1.93980e-01_wp, &
      & 1.58817e-01_wp, 1.30029e-01_wp, 1.06458e-01_wp, 8.71608e-02_wp, &
      & 7.13612e-02_wp, 5.84256e-02_wp, 4.78349e-02_wp, 3.91639e-02_wp, &
      & 3.20647e-02_wp, 2.62523e-02_wp, 2.14936e-02_wp, 1.75975e-02_wp, &
      & 1.44076e-02_wp, 1.17959e-02_wp, 9.65769e-03_wp                 /)

    REAL(wp), PARAMETER :: preflog(59) = (/ & !< log of ref pressure
      & 6.9600e+00_wp, 6.7600e+00_wp, 6.5600e+00_wp, 6.3600e+00_wp, &
      & 6.1600e+00_wp, 5.9600e+00_wp, 5.7600e+00_wp, 5.5600e+00_wp, &
      & 5.3600e+00_wp, 5.1600e+00_wp, 4.9600e+00_wp, 4.7600e+00_wp, &
      & 4.5600e+00_wp, 4.3600e+00_wp, 4.1600e+00_wp, 3.9600e+00_wp, &
      & 3.7600e+00_wp, 3.5600e+00_wp, 3.3600e+00_wp, 3.1600e+00_wp, &
      & 2.9600e+00_wp, 2.7600e+00_wp, 2.5600e+00_wp, 2.3600e+00_wp, &
      & 2.1600e+00_wp, 1.9600e+00_wp, 1.7600e+00_wp, 1.5600e+00_wp, &
      & 1.3600e+00_wp, 1.1600e+00_wp, 9.6000e-01_wp, 7.6000e-01_wp, &
      & 5.6000e-01_wp, 3.6000e-01_wp, 1.6000e-01_wp,-4.0000e-02_wp, &
      &-2.4000e-01_wp,-4.4000e-01_wp,-6.4000e-01_wp,-8.4000e-01_wp, &
      &-1.0400e+00_wp,-1.2400e+00_wp,-1.4400e+00_wp,-1.6400e+00_wp, &
      &-1.8400e+00_wp,-2.0400e+00_wp,-2.2400e+00_wp,-2.4400e+00_wp, &
      &-2.6400e+00_wp,-2.8400e+00_wp,-3.0400e+00_wp,-3.2400e+00_wp, &
      &-3.4400e+00_wp,-3.6400e+00_wp,-3.8400e+00_wp,-4.0400e+00_wp, &
      -4.2400e+00_wp,-4.4400e+00_wp,-4.6400e+00_wp                 /)
    !
    ! These are the temperatures associated with the respective
    ! pressures for the MLS standard atmosphere.
    !
    REAL(wp), PARAMETER :: tref(59) = (/ & !< reference temperatures
      & 2.9420e+02_wp, 2.8799e+02_wp, 2.7894e+02_wp, 2.6925e+02_wp, &
      & 2.5983e+02_wp, 2.5017e+02_wp, 2.4077e+02_wp, 2.3179e+02_wp, &
      & 2.2306e+02_wp, 2.1578e+02_wp, 2.1570e+02_wp, 2.1570e+02_wp, &
      & 2.1570e+02_wp, 2.1706e+02_wp, 2.1858e+02_wp, 2.2018e+02_wp, &
      & 2.2174e+02_wp, 2.2328e+02_wp, 2.2479e+02_wp, 2.2655e+02_wp, &
      & 2.2834e+02_wp, 2.3113e+02_wp, 2.3401e+02_wp, 2.3703e+02_wp, &
      & 2.4022e+02_wp, 2.4371e+02_wp, 2.4726e+02_wp, 2.5085e+02_wp, &
      & 2.5457e+02_wp, 2.5832e+02_wp, 2.6216e+02_wp, 2.6606e+02_wp, &
      & 2.6999e+02_wp, 2.7340e+02_wp, 2.7536e+02_wp, 2.7568e+02_wp, &
      & 2.7372e+02_wp, 2.7163e+02_wp, 2.6955e+02_wp, 2.6593e+02_wp, &
      & 2.6211e+02_wp, 2.5828e+02_wp, 2.5360e+02_wp, 2.4854e+02_wp, &
      & 2.4348e+02_wp, 2.3809e+02_wp, 2.3206e+02_wp, 2.2603e+02_wp, &
      & 2.2000e+02_wp, 2.1435e+02_wp, 2.0887e+02_wp, 2.0340e+02_wp, &
      & 1.9792e+02_wp, 1.9290e+02_wp, 1.8809e+02_wp, 1.8329e+02_wp, &
      & 1.7849e+02_wp, 1.7394e+02_wp, 1.7212e+02_wp                /)

    INTEGER :: jn, jt, jp, igc, ipr, iprsm
    REAL(wp):: zsumk, zsumf, zsumf1, zsumf2, zsumf3, zsumf4

    !++hs
    REAL(wp), PARAMETER :: ssi_default(14) =  (/ & !< SRTM default solar flux (W/m2) in 14 SW bands
      & 12.1095682699999987_wp   , 20.3650825429849398_wp   , 23.7297328613475429_wp   , &
      & 22.4276934179066849_wp   , 55.6266126199999960_wp   , 1.0293153082385436E+02_wp, &
      & 24.2936128100154392_wp   , 3.4574251380000004E+02_wp, 2.1818712729860400E+02_wp, &
      & 3.4719231470000005E+02_wp, 1.2949501812000000E+02_wp, 50.1522503011060650_wp   , &
      & 3.0799387010047838_wp    , 12.8893773299999985_wp  /)
    ! sum of 14 bands is: 1.3682223735968237E+03

    ! for the cases that use a Radiative-Convective equilibrium setup
    ! assumes a cos(zenith angle) = pi/4
    REAL(wp), PARAMETER :: ssi_rce(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
      & 3.803964_wp,  6.413186_wp,  7.44969_wp,   7.032908_wp,  17.63879_wp, &
      & 32.63096_wp,  7.861645_wp,  110.624_wp,   69.16621_wp,  109.3144_wp, &
      & 41.19017_wp,  15.00594_wp, 1.009717_wp,   4.195554_wp  /)
      ! sum of 14 bands is: 433.3371

    REAL(wp), PARAMETER :: ssi_amip(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
                                ! AMIP-type CMIP5 simulation (average from 1979-1988)
      & 11.95053_wp, 20.14766_wp, 23.40394_wp, 22.09458_wp, 55.41401_wp, &
      & 102.5134_wp, 24.69814_wp, 347.5362_wp, 217.2925_wp, 343.4221_wp, &
      & 129.403_wp , 47.14264_wp, 3.172126_wp, 13.18075_wp /)
    ! sum of 14 bands is: 1361.371

    REAL(wp), PARAMETER :: ssi_preind(14) =  (/ & !< solar flux (W/m2) in 14 SW bands for
                                ! preindustrial CMIP5 simulation (average from 1844-1856)
      & 11.95005_wp, 20.14612_wp, 23.40302_wp, 22.09443_wp, 55.41679_wp, &
      & 102.512_wp , 24.69536_wp, 347.4719_wp, 217.2217_wp, 343.2816_wp, &
      & 129.3001_wp, 47.07624_wp, 3.130199_wp, 13.17521_wp /)
    ! sum of 14 bands is: 1360.875

    REAL(wp), PARAMETER :: ssi_coddington(14) = (/ & !< solar flux (W/m2) in 14 SW bands from Coddington et al. (2016)
      & 12.045647_wp  , 20.257584_wp  , 23.604472_wp  , 23.37569292_wp, 57.56843759_wp, &
      & 105.6339255_wp, 24.72360028_wp, 345.7746485_wp, 213.5909065_wp, 344.8864993_wp, &
      & 128.6916773_wp, 45.19260459_wp, 2.825112161_wp, 12.821340_wp  /)
    !--hs

  CONTAINS

    SUBROUTINE setup_srtm

      INTEGER :: igc, igcsm, ibnd, ig, ind, ipr, iprsm
      REAL(wp)    :: zwtsum

      !-- read in the molecular absorption coefficients

      CALL srtm_kgb16
      CALL srtm_kgb17
      CALL srtm_kgb18
      CALL srtm_kgb19
      CALL srtm_kgb20
      CALL srtm_kgb21
      CALL srtm_kgb22
      CALL srtm_kgb23
      CALL srtm_kgb24
      CALL srtm_kgb25
      CALL srtm_kgb26
      CALL srtm_kgb27
      CALL srtm_kgb28
      CALL srtm_kgb29

      !Mike Iacono 20050804
      !-- Perform g-point reduction from 16 per band (224 total points) to
      !-- a band dependent number (112 total points) for all absorption
      !-- coefficient input data and Planck fraction input data.
      !-- Compute relative weighting for new g-point combinations.

      igcsm = 0
      DO ibnd = 1,jpsw
        iprsm = 0
        IF (ngc(ibnd) < jpg) THEN
          DO igc = 1,ngc(ibnd)
            igcsm = igcsm + 1
            zwtsum = 0.0_wp
            DO ipr = 1, ngn(igcsm)
              iprsm = iprsm + 1
              zwtsum = zwtsum + wt(iprsm)
            ENDDO
            wtsm(igc) = zwtsum
          ENDDO

          DO ig = 1,ng(ibnd+15)
            ind = (ibnd-1)*jpg + ig
            rwgt(ind) = wt(ig)/wtsm(ngm(ind))
          ENDDO
        ELSE
          DO ig = 1,ng(ibnd+15)
            igcsm = igcsm + 1
            ind = (ibnd-1)*jpg + ig
            rwgt(ind) = 1.0_wp
          ENDDO
        ENDIF
      ENDDO

      CALL srtm_cmbgb16
      CALL srtm_cmbgb17
      CALL srtm_cmbgb18
      CALL srtm_cmbgb19
      CALL srtm_cmbgb20
      CALL srtm_cmbgb21
      CALL srtm_cmbgb22
      CALL srtm_cmbgb23
      CALL srtm_cmbgb24
      CALL srtm_cmbgb25
      CALL srtm_cmbgb26
      CALL srtm_cmbgb27
      CALL srtm_cmbgb28
      CALL srtm_cmbgb29

      !-----------------------------------------------------------------------
    END SUBROUTINE setup_srtm

    !  Original version:       Michael J. Iacono; July, 1998
    !  Revision for RRTM_SW:   Michael J. Iacono; November, 2002
    !  Revision for RRTMG_SW:  Michael J. Iacono; December, 2003

    !  The subroutines CMBGB16->CMBGB29 input the absorption coefficient
    !  data for each band, which are defined for 16 g-points and 14 spectral
    !  bands. The data are combined with appropriate weighting following the
    !  g-point mapping arrays specified in RRTMG_SW_INIT.  Solar source
    !  function data in array SFLUXREF are combined without weighting.  All
    !  g-point reduced data are put into new arrays for use in RRTMG_SW.


    !-----------------------------------------------------------------------------
    !>
    !! BAND 16:  2600-3250 cm-1 (low key- H2O,CH4; high key - CH4)
    !
    SUBROUTINE srtm_cmbgb16

      USE mo_yoesrta16, ONLY : ka , kb , selfref , forref , sfluxref , &
        &                      kac, kbc, selfrefc, forrefc, sfluxrefc

      DO jn = 1,9
        DO jt = 1,5
          DO jp = 1,13
            iprsm = 0
            DO igc = 1,ngc(1)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(igc)
                iprsm = iprsm + 1
                zsumk = zsumk + ka(jn,jt,jp,iprsm)*rwgt(iprsm)
              ENDDO
              kac(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,5
        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(1)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(igc)
              iprsm = iprsm + 1
              zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm)
            ENDDO
            kbc(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(1)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,3
        iprsm = 0
        DO igc = 1,ngc(1)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      iprsm = 0
      DO igc = 1,ngc(1)
        zsumf = 0.0_wp
        DO ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          zsumf = zsumf + sfluxref(iprsm)
        ENDDO
        sfluxrefc(igc) = zsumf
      ENDDO

    END SUBROUTINE srtm_cmbgb16
    !-----------------------------------------------------------------------------
    !>
    !! BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
    !
    SUBROUTINE srtm_cmbgb17

      USE mo_yoesrta17, ONLY : ka , kb , selfref , forref , sfluxref , &
        &                      kac, kbc, selfrefc, forrefc, sfluxrefc

      DO jn = 1,9
        DO jt = 1,5
          DO jp = 1,13
            iprsm = 0
            DO igc = 1,ngc(2)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(1)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + ka(jn,jt,jp,iprsm)*rwgt(iprsm+16)
              ENDDO
              kac(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jn = 1,5
        DO jt = 1,5
          DO jp = 13,59
            iprsm = 0
            DO igc = 1,ngc(2)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(1)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + kb(jn,jt,jp,iprsm)*rwgt(iprsm+16)
              ENDDO
              kbc(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(2)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+16)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,4
        iprsm = 0
        DO igc = 1,ngc(2)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+16)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jp = 1,5
        iprsm = 0
        DO igc = 1,ngc(2)
          zsumf = 0.0_wp
          DO ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            zsumf = zsumf + sfluxref(iprsm,jp)
          ENDDO
          sfluxrefc(igc,jp) = zsumf
        ENDDO
      ENDDO
    END SUBROUTINE srtm_cmbgb17
    !-----------------------------------------------------------------------------
    !>
    !! BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
    !
    SUBROUTINE srtm_cmbgb18

      USE mo_yoesrta18, ONLY : ka , kb , selfref , forref , sfluxref , &
        &                      kac, kbc, selfrefc, forrefc, sfluxrefc

      DO jn = 1,9
        DO jt = 1,5
          DO jp = 1,13
            iprsm = 0
            DO igc = 1,ngc(3)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(2)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + ka(jn,jt,jp,iprsm)*rwgt(iprsm+32)
              ENDDO
              kac(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,5
        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(3)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(2)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm+32)
            ENDDO
            kbc(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(3)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+32)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,3
        iprsm = 0
        DO igc = 1,ngc(3)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+32)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jp = 1,9
        iprsm = 0
        DO igc = 1,ngc(3)
          zsumf = 0.0_wp
          DO ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            zsumf = zsumf + sfluxref(iprsm,jp)
          ENDDO
          sfluxrefc(igc,jp) = zsumf
        ENDDO
      ENDDO

    END SUBROUTINE srtm_cmbgb18
    !-----------------------------------------------------------------------------
    !>
    !! BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
    !
    SUBROUTINE srtm_cmbgb19

      USE mo_yoesrta19, ONLY : ka , kb , selfref , forref , sfluxref , &
        &                      kac, kbc, selfrefc, forrefc, sfluxrefc

      DO jn = 1,9
        DO jt = 1,5
          DO jp = 1,13
            iprsm = 0
            DO igc = 1,ngc(4)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(3)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + ka(jn,jt,jp,iprsm)*rwgt(iprsm+48)
              ENDDO
              kac(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,5
        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(4)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(3)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm+48)
            ENDDO
            kbc(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(4)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(3)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+48)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,3
        iprsm = 0
        DO igc = 1,ngc(4)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(3)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+48)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jp = 1,9
        iprsm = 0
        DO igc = 1,ngc(4)
          zsumf = 0.0_wp
          DO ipr = 1, ngn(ngs(3)+igc)
            iprsm = iprsm + 1
            zsumf = zsumf + sfluxref(iprsm,jp)
          ENDDO
          sfluxrefc(igc,jp) = zsumf
        ENDDO
      ENDDO

    END SUBROUTINE srtm_cmbgb19
    !-----------------------------------------------------------------------------
    !>
    !! BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
    !
    SUBROUTINE srtm_cmbgb20

      USE mo_yoesrta20, ONLY : ka , kb , selfref , forref , absch4 , sfluxref , &
        &                      kac, kbc, selfrefc, forrefc, absch4c, sfluxrefc

      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(5)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(4)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + ka(jt,jp,iprsm)*rwgt(iprsm+64)
            ENDDO
            kac(jt,jp,igc) = zsumk
          ENDDO
        ENDDO

        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(5)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(4)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm+64)
            ENDDO
            kbc(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(5)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+64)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,4
        iprsm = 0
        DO igc = 1,ngc(5)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+64)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      iprsm = 0
      DO igc = 1,ngc(5)
        zsumf1 = 0.0_wp
        zsumf2 = 0.0_wp
        DO ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          zsumf1 = zsumf1 + sfluxref(iprsm)
          zsumf2 = zsumf2 + absch4(iprsm)*rwgt(iprsm+64)
        ENDDO
        sfluxrefc(igc) = zsumf1
        absch4c(igc) = zsumf2
      ENDDO

    END SUBROUTINE srtm_cmbgb20
    !-----------------------------------------------------------------------------
    !>
    !! BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
    !
    SUBROUTINE srtm_cmbgb21

      USE mo_yoesrta21, ONLY : ka , kb , selfref , forref , sfluxref , &
        &                      kac, kbc, selfrefc, forrefc, sfluxrefc

      DO jn = 1,9
        DO jt = 1,5
          DO jp = 1,13
            iprsm = 0
            DO igc = 1,ngc(6)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(5)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + ka(jn,jt,jp,iprsm)*rwgt(iprsm+80)
              ENDDO
              kac(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jn = 1,5
        DO jt = 1,5
          DO jp = 13,59
            iprsm = 0
            DO igc = 1,ngc(6)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(5)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + kb(jn,jt,jp,iprsm)*rwgt(iprsm+80)
              ENDDO
              kbc(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(6)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(5)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+80)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,4
        iprsm = 0
        DO igc = 1,ngc(6)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(5)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+80)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jp = 1,9
        iprsm = 0
        DO igc = 1,ngc(6)
          zsumf = 0.0_wp
          DO ipr = 1, ngn(ngs(5)+igc)
            iprsm = iprsm + 1
            zsumf = zsumf + sfluxref(iprsm,jp)
          ENDDO
          sfluxrefc(igc,jp) = zsumf
        ENDDO
      ENDDO

    END SUBROUTINE srtm_cmbgb21
    !-----------------------------------------------------------------------------
    !>
    !! BAND 22:  7700-8050 cm-1 (low - H2O,O2; high - O2)
    !
    SUBROUTINE srtm_cmbgb22

      USE mo_yoesrta22, ONLY : ka , kb , selfref , forref , sfluxref , &
        &                      kac, kbc, selfrefc, forrefc, sfluxrefc

      DO jn = 1,9
        DO jt = 1,5
          DO jp = 1,13
            iprsm = 0
            DO igc = 1,ngc(7)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(6)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + ka(jn,jt,jp,iprsm)*rwgt(iprsm+96)
              ENDDO
              kac(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,5
        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(7)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(6)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm+96)
            ENDDO
            kbc(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(7)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+96)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,3
        iprsm = 0
        DO igc = 1,ngc(7)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+96)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jp = 1,9
        iprsm = 0
        DO igc = 1,ngc(7)
          zsumf = 0.0_wp
          DO ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            zsumf = zsumf + sfluxref(iprsm,jp)
          ENDDO
          sfluxrefc(igc,jp) = zsumf
        ENDDO
      ENDDO

    END SUBROUTINE srtm_cmbgb22
    !-----------------------------------------------------------------------------
    !>
    !! BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
    !
    SUBROUTINE srtm_cmbgb23

      USE mo_yoesrta23, ONLY : ka , selfref , forref , sfluxref , rayl , &
        &                      kac, selfrefc, forrefc, sfluxrefc, raylc

      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(8)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(7)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + ka(jt,jp,iprsm)*rwgt(iprsm+112)
            ENDDO
            kac(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(8)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+112)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,3
        iprsm = 0
        DO igc = 1,ngc(8)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+112)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      iprsm = 0
      DO igc = 1,ngc(8)
        zsumf1 = 0.0_wp
        zsumf2 = 0.0_wp
        DO ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          zsumf1 = zsumf1 + sfluxref(iprsm)
          zsumf2 = zsumf2 + rayl(iprsm)*rwgt(iprsm+112)
        ENDDO
        sfluxrefc(igc) = zsumf1
        raylc(igc) = zsumf2
      ENDDO

    END SUBROUTINE srtm_cmbgb23
    !-----------------------------------------------------------------------------
    !>
    !! BAND 24:  12850-16000 cm-1 (low - H2O,O2; high - O2)
    !
    SUBROUTINE srtm_cmbgb24

      USE mo_yoesrta24, ONLY : ka     , kb     , selfref , forref , sfluxref , &
        &                      abso3a , abso3b , rayla   , raylb  ,            &
        &                      kac    , kbc    , selfrefc, forrefc, sfluxrefc, &
        &                      abso3ac, abso3bc, raylac  , raylbc

      DO jn = 1,9
        DO jt = 1,5
          DO jp = 1,13
            iprsm = 0
            DO igc = 1,ngc(9)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(8)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + ka(jn,jt,jp,iprsm)*rwgt(iprsm+128)
              ENDDO
              kac(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,5
        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(9)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(8)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm+128)
            ENDDO
            kbc(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(9)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+128)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,3
        iprsm = 0
        DO igc = 1,ngc(9)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+128)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      iprsm = 0
      DO igc = 1,ngc(9)
        zsumf1 = 0.0_wp
        zsumf2 = 0.0_wp
        zsumf3 = 0.0_wp
        DO ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          zsumf1 = zsumf1 + raylb(iprsm)*rwgt(iprsm+128)
          zsumf2 = zsumf2 + abso3a(iprsm)*rwgt(iprsm+128)
          zsumf3 = zsumf3 + abso3b(iprsm)*rwgt(iprsm+128)
        ENDDO
        raylbc(igc) = zsumf1
        abso3ac(igc) = zsumf2
        abso3bc(igc) = zsumf3
      ENDDO

      DO jp = 1,9
        iprsm = 0
        DO igc = 1,ngc(9)
          zsumf1 = 0.0_wp
          zsumf2 = 0.0_wp
          DO ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            zsumf1 = zsumf1 + sfluxref(iprsm,jp)
            zsumf2 = zsumf2 + rayla(iprsm,jp)*rwgt(iprsm+128)
          ENDDO
          sfluxrefc(igc,jp) = zsumf1
          raylac(igc,jp) = zsumf2
        ENDDO
      ENDDO

    END SUBROUTINE srtm_cmbgb24
    !-----------------------------------------------------------------------------
    !>
    !! BAND 25:  16000-22650 cm-1 (low - H2O; high - nothing)
    !
    SUBROUTINE srtm_cmbgb25

      USE mo_yoesrta25, ONLY : ka , sfluxref , abso3a , abso3b , rayl , &
        &                      kac, sfluxrefc, abso3ac, abso3bc, raylc

      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(10)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(9)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + ka(jt,jp,iprsm)*rwgt(iprsm+144)
            ENDDO
            kac(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      iprsm = 0
      DO igc = 1,ngc(10)
        zsumf1 = 0.0_wp
        zsumf2 = 0.0_wp
        zsumf3 = 0.0_wp
        zsumf4 = 0.0_wp
        DO ipr = 1, ngn(ngs(9)+igc)
          iprsm = iprsm + 1
          zsumf1 = zsumf1 + sfluxref(iprsm)
          zsumf2 = zsumf2 + abso3a(iprsm)*rwgt(iprsm+144)
          zsumf3 = zsumf3 + abso3b(iprsm)*rwgt(iprsm+144)
          zsumf4 = zsumf4 + rayl(iprsm)*rwgt(iprsm+144)
        ENDDO
        sfluxrefc(igc) = zsumf1
        abso3ac(igc) = zsumf2
        abso3bc(igc) = zsumf3
        raylc(igc) = zsumf4
      ENDDO

    END SUBROUTINE srtm_cmbgb25
    !-----------------------------------------------------------------------------
    !>
    !! BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)
    !
    SUBROUTINE srtm_cmbgb26

      USE mo_yoesrta26, ONLY : sfluxref , rayl , &
        &                      sfluxrefc, raylc

      iprsm = 0
      DO igc = 1,ngc(11)
        zsumf1 = 0.0_wp
        zsumf2 = 0.0_wp
        DO ipr = 1, ngn(ngs(10)+igc)
          iprsm = iprsm + 1
          zsumf1 = zsumf1 + rayl(iprsm)*rwgt(iprsm+160)
          zsumf2 = zsumf2 + sfluxref(iprsm)
        ENDDO
        raylc(igc) = zsumf1
        sfluxrefc(igc) = zsumf2
      ENDDO

    END SUBROUTINE srtm_cmbgb26
    !-----------------------------------------------------------------------------
    !>
    !! BAND 27:  29000-38000 cm-1 (low - O3; high - O3)
    !
    SUBROUTINE srtm_cmbgb27

      USE mo_yoesrta27, ONLY : ka , kb , sfluxref , rayl , &
        &                      kac, kbc, sfluxrefc, raylc

      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(12)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(11)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + ka(jt,jp,iprsm)*rwgt(iprsm+176)
            ENDDO
            kac(jt,jp,igc) = zsumk
          ENDDO
        ENDDO

        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(12)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(11)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm+176)
            ENDDO
            kbc(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      iprsm = 0
      DO igc = 1,ngc(12)
        zsumf1 = 0.0_wp
        zsumf2 = 0.0_wp
        DO ipr = 1, ngn(ngs(11)+igc)
          iprsm = iprsm + 1
          zsumf1 = zsumf1 + sfluxref(iprsm)
          zsumf2 = zsumf2 + rayl(iprsm)*rwgt(iprsm+176)
        ENDDO
        sfluxrefc(igc) = zsumf1
        raylc(igc) = zsumf2
      ENDDO

    END SUBROUTINE srtm_cmbgb27
    !-----------------------------------------------------------------------------
    !>
    !! BAND 28:  38000-50000 cm-1 (low - O3,O2; high - O3,O2)
    !
    SUBROUTINE srtm_cmbgb28

      USE mo_yoesrta28, ONLY : ka , kb , sfluxref , &
        &                      kac, kbc, sfluxrefc

      DO jn = 1,9
        DO jt = 1,5
          DO jp = 1,13
            iprsm = 0
            DO igc = 1,ngc(13)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(12)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + ka(jn,jt,jp,iprsm)*rwgt(iprsm+192)
              ENDDO
              kac(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jn = 1,5
        DO jt = 1,5
          DO jp = 13,59
            iprsm = 0
            DO igc = 1,ngc(13)
              zsumk = 0.0_wp
              DO ipr = 1, ngn(ngs(12)+igc)
                iprsm = iprsm + 1
                zsumk = zsumk + kb(jn,jt,jp,iprsm)*rwgt(iprsm+192)
              ENDDO
              kbc(jn,jt,jp,igc) = zsumk
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO jp = 1,5
        iprsm = 0
        DO igc = 1,ngc(13)
          zsumf = 0.0_wp
          DO ipr = 1, ngn(ngs(12)+igc)
            iprsm = iprsm + 1
            zsumf = zsumf + sfluxref(iprsm,jp)
          ENDDO
          sfluxrefc(igc,jp) = zsumf
        ENDDO
      ENDDO

    END SUBROUTINE srtm_cmbgb28
    !-----------------------------------------------------------------------------
    !>
    !! BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
    !
    SUBROUTINE srtm_cmbgb29

      USE mo_yoesrta29, ONLY : ka , kb , selfref , forref , sfluxref , absh2o , absco2 , &
        &                      kac, kbc, selfrefc, forrefc, sfluxrefc, absh2oc, absco2c

      DO jt = 1,5
        DO jp = 1,13
          iprsm = 0
          DO igc = 1,ngc(14)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(13)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + ka(jt,jp,iprsm)*rwgt(iprsm+208)
            ENDDO
            kac(jt,jp,igc) = zsumk
          ENDDO
        ENDDO

        DO jp = 13,59
          iprsm = 0
          DO igc = 1,ngc(14)
            zsumk = 0.0_wp
            DO ipr = 1, ngn(ngs(13)+igc)
              iprsm = iprsm + 1
              zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm+208)
            ENDDO
            kbc(jt,jp,igc) = zsumk
          ENDDO
        ENDDO
      ENDDO

      DO jt = 1,10
        iprsm = 0
        DO igc = 1,ngc(14)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+208)
          ENDDO
          selfrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      DO jt = 1,4
        iprsm = 0
        DO igc = 1,ngc(14)
          zsumk = 0.0_wp
          DO ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+208)
          ENDDO
          forrefc(jt,igc) = zsumk
        ENDDO
      ENDDO

      iprsm = 0
      DO igc = 1,ngc(14)
        zsumf1 = 0.0_wp
        zsumf2 = 0.0_wp
        zsumf3 = 0.0_wp
        DO ipr = 1, ngn(ngs(13)+igc)
          iprsm = iprsm + 1
          zsumf1 = zsumf1 + sfluxref(iprsm)
          zsumf2 = zsumf2 + absco2(iprsm)*rwgt(iprsm+208)
          zsumf3 = zsumf3 + absh2o(iprsm)*rwgt(iprsm+208)
        ENDDO
        sfluxrefc(igc) = zsumf1
        absco2c(igc) = zsumf2
        absh2oc(igc) = zsumf3
      ENDDO
    END SUBROUTINE srtm_cmbgb29

  END MODULE mo_srtm_config
