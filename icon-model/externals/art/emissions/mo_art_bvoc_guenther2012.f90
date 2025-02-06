!
! mo_art_bvoc
! Module including parameters for subroutine bvoc (biogenic VOC emissions)
!------------------------------------------------------------------------------
! parameters taken from BVOC emission parameterization:
! "The Model of Emissions of Gases and Aerosols from Nature version 2.1
! (MEGAN2.1): an extended and updated framework for modeling biogenic
! emissions"; Guenther et al. (2012), Geosci. model Dev., 5, 1471-1492
! or estimated if not available
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

MODULE mo_art_bvoc

! ICON
  USE mo_kind,                          ONLY: wp,i4
! ART

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bvoc_guenther2012, sel_bccnum, pftn

! ADAPT ACCORDING TO YOUR SIMULATION DOMAIN #####################
! climate zone parameter to choose mean temperatures 
! (24h mean temperature (T24) and 240h mean temperature (T240))
! and mean Photosynthetically active radiation 
! (24h mean PAR (PPFD24) and 240h mean PAR (PPFD240)) [see below]
! for 
! temperate zone: czone = 1 
! or 
! tropical zone : czone = 2
! [not a part of Guenther et al. (2012) but necessary because
! T24 and T240 are not available within the simulation].
! ***Please check whether the default values of t24, t240, ppfd24
! and ppfd240 fit to your simulation domain and time.***
  INTEGER(KIND=i4),PARAMETER  :: czone  = 2_i4
! ###############################################################

! number of czones
  INTEGER(KIND=i4), PARAMETER   :: nczone = 2_i4

! CO2 concentration [ppm]
! [not a constant but set to const owing to lack of data]
  REAL(KIND=wp),PARAMETER          :: co2     = 400._wp

! BVOC compound class number [Guenther et al. (2012), Tab.4]:
  INTEGER(KIND=i4), PARAMETER      :: bccn    = 19_i4
! 1  Isoprene
! 2  Myrcene 
! 3  Sabinene 
! 4  Limonene 
! 5  3-Carene
! 6  t-beta-Ocimene 
! 7  beta-Pinene 
! 8  alpha-Pinene 
! 9  Other Monoterpenes 
! 10 alpha-Farnesene 
! 11 beta-Caryophyllene 
! 12 Other Sesquiterpenes 
! 13 232-MBO
! 14 Methanol 
! 15 Acetone 
! 16 CO
! 17 Bidirectional VOC 
! 18 Stress VOC 
! 19 Other VOC 

! Select compound classes for the analysis in COSMO-ART:
  INTEGER(KIND=i4), PARAMETER      :: sel_bccnum = 19_i4
  INTEGER(KIND=i4)                 :: sel_bccn(sel_bccnum)
! data sel_bccn / 1,8,4 / ! Isoprene, alpha-Pinene, Limonene
  data sel_bccn / 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 /

! Plant functional type number [from CLM4 global land area; Guenther et al. (2012), Tab.3]:
  INTEGER(KIND=i4), PARAMETER      :: pftn    = 15_i4

! Emission factors EF (for BVOC compound class number and 
! vegetation type number) [Guenther et al. (2012), Tab.2]:
  REAL(KIND=wp)                           :: epsil(pftn,bccn)

! pftn          1       2/9     3/10      4/11       5/12      6/13       7/14       8/15   bccn
  DATA epsil / 600._wp, 3000._wp,   1._wp, 7000._wp, 10000._wp, 7000._wp, 10000._wp, 11000._wp,  &!1 
    &                   2000._wp,4000._wp, 4000._wp,  1600._wp,  800._wp,   200._wp,     1._wp,  &!1 
    &           70._wp,   70._wp,  60._wp,   80._wp,    30._wp,   80._wp,    30._wp,    30._wp,  &!2 
    &                     30._wp,  50._wp,   30._wp,    0.3_wp,   0.3_wp,    0.3_wp,    0.3_wp,  &!2 
    &           70._wp,   70._wp,  40._wp,   80._wp,    50._wp,   80._wp,    50._wp,    50._wp,  &!3 
    &                     50._wp,  70._wp,   50._wp,    0.7_wp,   0.7_wp,    0.7_wp,    0.7_wp,  &!3 
    &          100._wp,  100._wp, 130._wp,   80._wp,    80._wp,   80._wp,    80._wp,    80._wp,  &!4 
    &                     60._wp, 100._wp,   60._wp,    0.7_wp,   0.7_wp,    0.7_wp,    0.7_wp,  &!4 
    &          160._wp,  160._wp,  80._wp,   40._wp,    30._wp,   40._wp,    30._wp,    30._wp,  &!5 
    &                     30._wp, 100._wp,   30._wp,    0.3_wp,   0.3_wp,    0.3_wp,    0.3_wp,  &!5 
    &           70._wp,   70._wp,  60._wp,  150._wp,   120._wp,  150._wp,   120._wp,   120._wp,  &!6 
    &                     90._wp, 150._wp,   90._wp,     2._wp,    2._wp,     2._wp,     2._wp,  &!6 
    &          300._wp,  300._wp, 200._wp,  120._wp,   130._wp,  120._wp,   130._wp,   130._wp,  &!7 
    &                    100._wp, 150._wp,  100._wp,    1.5_wp,   1.5_wp,    1.5_wp,    1.5_wp,  &!7 
    &          500._wp,  500._wp, 510._wp,  600._wp,   400._wp,  600._wp,   400._wp,   400._wp,  &!8 
    &                    200._wp, 300._wp,  200._wp,     2._wp,    2._wp,     2._wp,     2._wp,  &!8 
    &          180._wp,  180._wp, 170._wp,  150._wp,   150._wp,  150._wp,   150._wp,   150._wp,  &!9 
    &                    110._wp, 200._wp,  110._wp,     5._wp,    5._wp,     5._wp,     5._wp,  &!9 
    &           40._wp,   40._wp,  40._wp,   60._wp,    40._wp,   60._wp,    40._wp,    40._wp,  &!10
    &                     40._wp,  40._wp,   40._wp,     3._wp,    3._wp,     3._wp,     4._wp,  &!10
    &           80._wp,   80._wp,  80._wp,   60._wp,    40._wp,   60._wp,    40._wp,    40._wp,  &!11
    &                     50._wp,  50._wp,   50._wp,     1._wp,    1._wp,     1._wp,     4._wp,  &!11
    &          120._wp,  120._wp, 120._wp,  120._wp,   100._wp,  120._wp,   100._wp,   100._wp,  &!12
    &                    100._wp, 100._wp,  100._wp,     2._wp,    2._wp,     2._wp,     2._wp,  &!12
    &          700._wp,   60._wp, 0.01_wp,  0.01_wp,   0.01_wp,  0.01_wp,   0.01_wp,     2._wp,  &!13
    &                    0.01_wp, 0.01_wp,  0.01_wp,   0.01_wp,  0.01_wp,   0.01_wp,   0.01_wp,  &!13
    &          900._wp,  900._wp, 900._wp,  500._wp,   900._wp,  500._wp,   900._wp,   900._wp,  &!14
    &                    900._wp, 900._wp,  900._wp,   500._wp,  500._wp,   500._wp,   900._wp,  &!14
    &          240._wp,  240._wp, 240._wp,  240._wp,   240._wp,  240._wp,   240._wp,   240._wp,  &!15
    &                    240._wp, 240._wp,  240._wp,    80._wp,   80._wp,    80._wp,    80._wp,  &!15
    &          600._wp,  600._wp, 600._wp,  600._wp,   600._wp,  600._wp,   600._wp,   600._wp,  &!16
    &                    600._wp, 600._wp,  600._wp,   600._wp,  600._wp,   600._wp,   600._wp,  &!16
    &          500._wp,  500._wp, 500._wp,  500._wp,   500._wp,  500._wp,   500._wp,   500._wp,  &!17
    &                    500._wp, 500._wp,  500._wp,    80._wp,   80._wp,    80._wp,    80._wp,  &!17
    &          300._wp,  300._wp, 300._wp,  300._wp,   300._wp,  300._wp,   300._wp,   300._wp,  &!18
    &                    300._wp, 300._wp,  300._wp,   300._wp,  300._wp,   300._wp,   300._wp,  &!18
    &          140._wp,  140._wp, 140._wp,  140._wp,   140._wp,  140._wp,   140._wp,   140._wp,  &!19
    &                    140._wp, 140._wp,  140._wp,   140._wp,  140._wp,   140._wp,   140._wp   /!19

! light dependent fraction of the bccn 1-19 LDF [Guenther et al. (2012), Tab.4]:
  REAL(KIND=wp)                       :: ldf(bccn)
! bccn     1/11   2/12    3/13    4/14   5/15   6/16    7/17     8/18   9/19      10
  DATA ldf /1._wp, 0.6_wp, 0.6_wp, 0.2_wp, 0.2_wp, 0.8_wp, 0.2_wp, 0.6_wp, 0.4_wp, 0.5_wp,  &
    &      0.5_wp, 0.5_wp,  1._wp, 0.8_wp, 0.2_wp,  1._wp, 0.8_wp, 0.8_wp, 0.2_wp/

! standard conditions for PPFD averaged over the past 24h P_S [umol m-2 s-1]
! (200 for sun leaves, 50 for shade leaves)
! [not a constant but set to const (avrg between sun and shade) owing to lack of data] 
  REAL(KIND=wp), PARAMETER            :: ppfds   = 125._wp 

! average ppfd of the past 24h P_24 [umol m-2 s-1]
! [not a constant but set to const owing to lack of data (COSMO simulation temporal mean 
! of June 1, 2014 and field mean 
! over SW-Africa delivered 421 umol m-2 s-1]
  REAL(KIND=wp)                       :: ppfd24(nczone)
!            temperate zone, tropical zone 
  DATA ppfd24 /400._wp, 400._wp/

! average ppfd of the past 240h P_240 [umol m-2 s-1]
! [not a constant but set to const owing to lack of data (COSMO simulation temporal mean 
! of June 1, 2014 and field mean 
! over SW-Africa delivered 421 umol m-2 s-1]
  REAL(KIND=wp)                       :: ppfd240(nczone)
!             temperate zone, tropical zone 
  DATA ppfd240 /400._wp, 400._wp/

! empirically determined coefficients of the bccn 1-19 C_T1 [Guenther et al. (2012), Tab.4]:
  REAL(KIND=wp)                       :: ct1(bccn)
! bccn     1/11   2/12    3/13    4/14   5/15   6/16    7/17     8/18   9/19      10
  DATA ct1 /95._wp, 80._wp, 80._wp, 80._wp, 80._wp, 80._wp, 80._wp, 80._wp, 80._wp, 130._wp,    &
    &      130._wp,130._wp, 95._wp, 60._wp, 80._wp, 60._wp, 95._wp, 80._wp, 80._wp/

! empirically determined coefficient C_T2 
  REAL(KIND=wp), PARAMETER            :: ct2     = 230._wp

! leave temperature at standard conditions T_S [K]
  REAL(KIND=wp), PARAMETER            :: ts      = 297._wp 

! average leave temperature of the past 24h T_24 [K]
! [not a constant but set to const (15 degC for temperate zone and 23.85 degC (equal to T_S) for
! tropical zone) owing to lack of data] 
  REAL(KIND=wp)                       :: t24(nczone)
!         temperate zone, tropical zone 
  DATA t24 /288.15_wp, 297._wp/

! average leave temperature of the past 240h T_240 [K]
! [not a constant but set to const (15 degC for temperate zone and 23.85 degC (equal to T_S) for
! tropical zone) owing to lack of data] 
  REAL(KIND=wp)                       :: t240(nczone)
!          temperate zone, tropical zone 
  DATA t240 /288.15_wp, 297._wp/

! empirically determined coefficients of the bccn 1-19 C_eo [Guenther et al. (2012), Tab.4]:
  REAL(KIND=wp)                       :: ceo(bccn)
! bccn     1/11   2/12    3/13    4/14   5/15   6/16    7/17     8/18   9/19      10
  DATA ceo /2._wp  ,1.83_wp,1.83_wp,1.83_wp,1.83_wp,1.83_wp,1.83_wp,1.83_wp,1.83_wp,2.37_wp,   &
    &       2.37_wp,2.37_wp,  2._wp, 1.6_wp,1.83_wp, 1.6_wp,  2._wp,1.83_wp,1.83_wp/

! empirically determined coefficients of the bccn 1-19 beta [Guenther et al. (2012), Tab.4]:
  REAL(KIND=wp)                       :: beta(bccn)
! bccn       1/11   2/12    3/13    4/14   5/15   6/16     7/17    8/18   9/19      10
  DATA beta /0.13_wp, 0.1_wp, 0.1_wp, 0.1_wp,0.1_wp, 0.1_wp, 0.1_wp,0.1_wp,0.1_wp,0.17_wp,   &
    &        0.17_wp,0.17_wp,0.13_wp,0.08_wp,0.1_wp,0.08_wp,0.13_wp,0.1_wp,0.1_wp/

! empirically determined coefficients of the bccn 1-19 A_new 
! (relative emission rates for new leaves) [Guenther et al. (2012), Tab.4]:
  REAL(KIND=wp)                       :: anew(bccn)
! bccn       1/11   2/12    3/13     4/14    5/15   6/16  7/17    8/18   9/19      10
  DATA anew /0.05_wp,  2._wp,   2._wp,  2._wp, 2._wp, 2._wp, 2._wp, 2._wp, 2._wp, 0.4_wp,  &
    &         0.4_wp, 0.4_wp, 0.05_wp, 3.5_wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp/

! empirically determined coefficients of the bccn 1-19 A_gro 
! (relative emission rates for growing leaves) [Guenther et al. (2012), Tab.4]:
  REAL(KIND=wp)                       :: agro(bccn)
! bccn       1/11   2/12    3/13    4/14   5/15   6/16     7/17    8/18   9/19      10
  DATA agro /0.6_wp, 1.8_wp, 1.8_wp, 1.8_wp, 1.8_wp, 1.8_wp, 1.8_wp, 1.8_wp, 1.8_wp, 0.6_wp,   &
   &         0.6_wp, 0.6_wp, 0.6_wp,  3._wp,  1._wp,  1._wp,  1._wp,  1._wp,  1._wp/

! empirically determined coefficients of the bccn 1-19 A_mat
! (relative emission rates for mature leaves) [Guenther et al. (2012), Tab.4]:
  REAL(KIND=wp)                       :: amat(bccn)
! bccn       1/11  2/12   3/13   4/14  5/15   6/16   7/17   8/18   9/19     10
  DATA amat /1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp,  &
    &        1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp, 1._wp/

! empirically determined coefficients of the bccn 1-19 A_sen
! (relative emission rates for senescing leaves) [Guenther et al. (2012), Tab.4]:
  REAL(KIND=wp)                       :: asen(bccn)
! bccn       1/11   2/12    3/13    4/14   5/15   6/16     7/17    8/18   9/19      10
  DATA asen /0.9_wp,1.05_wp,1.05_wp,1.05_wp,1.05_wp,1.05_wp,1.05_wp,1.05_wp,1.05_wp,0.95_wp,  &
   &        0.95_wp,0.95_wp, 0.9_wp, 1.2_wp,  1._wp,  1._wp,  1._wp,  1._wp,  1._wp/

! fraction of new foliage F_new [0 1]
! [not a constant but set to const owing to lack of data] 
  REAL(KIND=wp), PARAMETER            :: fnew    = 0.25_wp

! fraction of growing foliage F_gro [0 1]  
! [not a constant but set to const owing to lack of data] 
  REAL(KIND=wp), PARAMETER            :: fgro    = 0.25_wp

! fraction of mature foliage F_mat [0 1] 
! [not a constant but set to const owing to lack of data] 
  REAL(KIND=wp), PARAMETER            :: fmat    = 0.25_wp

! fraction of senescing foliage F_sen [0 1]
! [not a constant but set to const owing to lack of data] 
  REAL(KIND=wp), PARAMETER            :: fsen    = 0.25_wp

! empirical determined coefficient delta_theta_1 [m3 m-3]
  REAL(KIND=wp), PARAMETER            :: dtheta1 = 0.04_wp 

! estimated asymptote at which further decreases in intercellular CO2
! have a negligible effect on isoprene emission I_Smax [Heald et al.(2009), Tab. 2]
  REAL(KIND=wp), PARAMETER            :: ismax   = 1.344_wp

! scaling coefficient C_star [Heald et al.(2009), Tab. 2]
  REAL(KIND=wp), PARAMETER            :: cstar   = 585._wp

! exponential scalar h [Heald et al.(2009), Tab. 2]
  REAL(KIND=wp), PARAMETER            :: h       = 1.4614_wp    

! canopy environment coefficient C_CE
  REAL(KIND=wp), PARAMETER            :: cce     = 0.30_wp

! emperical parameter CHI related to leave angle distribution 
! (for calculation of LAI fraction lit by sun (Dai et al. 2004, Eq.2;
! 1=horizontal,-1=vertical,0=spherical leave distribution)
! [not a constant but set to const owing to lack of data] 

! 1  Needleleaf Evergreen Temperate Tree   
! 2  Needleleaf Evergreen Boreal Tree 
! 3  Needleleaf Deciduous Boreal Tree 
! 4  Broadleaf Evergreen Tropical Tree 
! 5  Broadleaf Evergreen Temperate Tree 
! 6  Broadleaf Deciduous Tropical Tree 
! 7  Broadleaf Deciduous Temperate Tree
! 8  Broadleaf Deciduous Boreal Tree 
! 9  Broadleaf Evergreen Temperate Shrub 
! 10 Broadleaf Deciduous Temperate Shrub 
! 11 Broadleaf Deciduous Boreal Shrub
! 12 Arctic C3 Grass 
! 13 Cool C3 Grass 
! 14 Warm C4 Grass 
! 15 Crop1
  REAL(KIND=wp)                       :: chi(pftn)
! assuming needleleaf to be spherical, broadleaf to be horizontal and grass to be vertical
! (leads not to realistic Isoprene values in Central Africa)
!  pftn      1     2     3     4     5     6     7     8     9    10    11    12   13   14   15 
!DATA chi /0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0, -1.0,-1.0,-1.0,-1.0/

! chi = 0.8 leads to realistic isoprene values in Central Africa
! pftn      1       2/9     3/10      4/11     5/12      6/13    7/14     8/15
  DATA chi /0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp,   &
    &                0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp,  0.8_wp/


! Subroutine BVOC for biogenic VOC Emission
!------------------------------------------------------------------------------
CONTAINS

SUBROUTINE bvoc_guenther2012(t,pabs,                                 &
!                &            wsoil,cpwp, soiltyp,                    &
                &            pft,lai,                                &
!                &            sun_sza,                                &
                &            istart,iend,cbio,laisun)
!<
! SUBROUTINE bvoc_guenther2012
! This subroutine calculates biogenic VOC emissions
! Based on: The Model of Emissions of Gases and Aerosols from Nature version 2.1
!           (MEGAN2.1): an extended and updated framework for modeling biogenic 
!           emissions"; Guenther et al. (2012), Geosci. model Dev., 5, 1471-1492
! Part of Module: mo_art_bvoc
! Author: Konrad Deetz, KIT
! Initial Release: 2015-04-01
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - ...
!>
 ! IN
  INTEGER(KIND=i4),INTENT(IN)        :: istart, iend
  REAL(KIND=wp),INTENT(IN)           :: t(:)        ! Temperature           (i)      [K]
  REAL(KIND=wp),INTENT(IN)           :: pabs(:)     ! Photosynthetic active radiation 
                                                    ! at the ground         (i)      [W m-2]]
!  REAL(KIND=wp),INTENT(IN)           :: wsoil(:)    ! Soil moisture         (i)      [m3 m-3] 
!  REAL(KIND=wp),INTENT(IN)           :: cpwp(:)     ! Plant wilting point   (soiltyp 1-10)
                                                    !   [fraction of volume]
!  REAL(KIND=wp),INTENT(IN)           :: soiltyp(:)  ! soiltyp               (i)      [-] 
  REAL(KIND=wp),INTENT(IN)           :: pft(:,:)    ! Plant functional type (i,npft) [-]
  REAL(KIND=wp),INTENT(IN)           :: lai(:)      ! Leaf area index       (i)      [m2 m-2] 
!  REAL(KIND=wp),INTENT(IN)           :: sun_sza(:)  ! Solar zenith angle    (i)      [deg] 
 
 ! INOUT
  REAL(KIND=wp),INTENT(INOUT)        :: cbio(:,:)   ! biogenic emissions     (i,spec) [kg m-2 h-1]
  REAL(KIND=wp),INTENT(INOUT)        :: laisun(:)   ! LAI that is lit by sun (i)      [m2 m-2] 
 
 ! local variables
  INTEGER(KIND=i4)                   :: iix, lbccnum, lpftn
  REAL(KIND=wp)                      :: par
!  REAL(KIND=wp)                      :: zpwp
  REAL(KIND=wp)                      :: lai_sun(pftn)!, phi_1, phi_2!, g_mu, kb
!  REAL(KIND=wp)                      :: sza
   
 !=======================================================================
 
 
 ! Call emission calculation function
 ! ######################################################################
 
  DO iix = istart,iend              ! lat 
    DO lbccnum = 1,sel_bccnum       ! BVOC compound class 
 
      ! Check LAI value
      ! LAI is sometimes negative in coastal areas (fill values), which could lead
      ! to problems in the BVOC routine. Owing the fact that no BVOC emissions
      ! with a LAI<=0 are possible, the calculation starts only with LAI>0.
      IF (lai(iix) > 0._wp) THEN
      
      ! PREPARATIONS
      
      ! 1. Calculate photosynthetic photon flux density PPFD [umol m-2 s-1]:
      ! Transform w m-2 s-1 to umol m-2 s-1, the factor 4.56 is defined for T=5800K (daylight)
        par = pabs(iix) * 4.56_wp
        IF (pabs(iix) < 1._wp) par = 0._wp
      
      ! 2. Calculate plant wilting point with soil type (by index)
      ! soil type: ice, rock, sand, sandy loam, loam, clay loam, clay, peat, sea water, sea ice
      
      ! ONLY a place holder for the moment. Soil type should be integer????
      ! In a first attempt we just make it to Integer
      !  zpwp = 1._wp !cpwp(INT(soiltyp(iix))) 
      
      ! 3.1 Calculate sun zenith angle
      !  sza = sun_sza(iix)
      
      ! 3.2 Calculate leaf area that is lit by sun (LAIsun)
      ! based on Dai et al. (2004): "A Two-Big-Leaf Model for Canopy Temperature,
      ! Photosynthesis, and Stomatal Conductance", Eq. 2,3 (also used in CLM4)
      ! (not part of the BVOC emission parameterization of Guenther et al. (2012),
      ! but neglecting the canopy shading would lead to emission overestimations)
      
      ! loop over plant functional types
        lai_sun(:) = 0._wp
        DO lpftn = 1,pftn
          !phi_1 = 0.5_wp - (0.633_wp * chi(lpftn)) - (0.33_wp *chi(lpftn)*chi(lpftn))
          !phi_2 = 0.877_wp * (1._wp - 2._wp*phi_1)
          ! SZA has very small influence on lai_sun (but can lead to unrealistic high
          ! values if sza < 2degree) ...
          !G_mu  = phi_1 + (phi_2*sza)
          !kb    = G_mu/sza
          ! ... use constant value sza=10.4degree (representation of mean conditions) instead. [KD 2015]
          !G_mu  = phi_1 + (phi_2 * COS(deg2rad * 10.4_wp))
          !kb    = G_mu/COS(deg2rad * 10.4_wp)
          lai_sun(lpftn) = lai(iix) !  (1._wp-EXP(-kb*lai(iix))) / kb  !  
        ENDDO ! lpftn
      
      ! 3.3 LAIsun as debug gribout variable LAI_SUN
        laisun(iix) = lai_sun(4_i4) ! 4 = Broadleaf Evergreen Tropical Tree
      
      ! ##############################
      ! call bvoc calculation function
        cbio(iix,lbccnum) = calc_bvoc(t(iix),par ,   & ! wsoil(iix),zpwp  ,    &
          &                           sel_bccn(lbccnum),pft(iix,:),lai_sun)
      ! parameterization function:    temp  ,ppfd,theta     ,thetaw,    &
      !   &                           lbccn            ,pft       ,leafai
      ! ##############################
      
      ! leaf area is zero - no BVOC emissions possible
      ELSE 
        cbio(iix,lbccnum) = 0._wp
      ENDIF
      
    ENDDO ! BVOC compound class
  ENDDO ! lon
 
 
 
 ! ######################################################################
 ! 4. Function definitions
 ! ######################################################################
 
 ! 4.1 Definition of emission calculation function
 ! ######################################################################
  CONTAINS
   
    REAL(KIND=wp) FUNCTION calc_bvoc(temp,ppfd,                    &
 !                            &    theta,thetaw,                 &
                              &    lbccn,pfts,leafai)
 
 ! loop variables
    INTEGER(KIND=i4)              :: lbccn
 
 ! 1. Data input (model fields and parameters)
 ! ######################################################################
 
 ! 1.1 Meteorological/vegetation input
 ! ----------------------------------------------------------------------
 ! 1.1.1 Air temperature [K]
 ! from data_fields: t 
    REAL(KIND=wp)                 :: temp
 
 ! 1.1.2 Photosynthetic active radiation at the ground [W m-2]
 ! from data_fields: pabs
 
 ! 1.1.3 Volumetric soil moisture [m3 m-3]
 ! from data_fields: w_g1, w_so
 !   REAL(KIND=wp)                 :: theta
 
 ! 1.1.4 Wilting point [m3 m-3]
 ! from src_soil_multlay: zpwp
 !   REAL(KIND=wp)                 :: thetaw
 
 ! 1.1.5 CO2 concentration [ppm]
 ! from mod_bvoc: co2
 
 ! 1.1.6 Leaf area index (lit by sun, depending on PFT!)
    REAL(KIND=wp)                 :: leafai(pftn)
 
 ! 1.1.7 Fractional area coverage of Plant functional type (PFT) 1-15
 ! from data_cosmo_art: pft
    REAL(KIND=wp)                 :: pfts(pftn)
 
 
 ! 1.2 Parameters for the BVOC parameterization
 ! ----------------------------------------------------------------------
 
 ! 1.2.1 Emission activity factor gamma [Guenther et al. (2012), Eq. 2]
    REAL(KIND=wp)                 :: gamma
 
 ! 1.2.2 Parameters for calculating total emission activity factor gamma
 !      [Guenther et al. (2012), Eq. 2] - canopy environment coefficient C_CE
 ! from mod_bvoc:  cce
 
 ! 1.2.3 BVOC compound class number [Guenther et al. (2012), Tab.4]:
 ! from mod_bvoc: bvoc=19
 
 ! 1  Isoprene
 ! 2  Myrcene 
 ! 3  Sabinene 
 ! 4  Limonene 
 ! 5  3-Carene
 ! 6  t-beta-Ocimene 
 ! 7  beta-Pinene 
 ! 8  alpha-Pinene 
 ! 9  Other Monoterpenes 
 ! 10 alpha-Farnesene 
 ! 11 beta-Caryophyllene 
 ! 12 Other Sesquiterpenes 
 ! 13 232-MBO
 ! 14 Methanol 
 ! 15 Acetone 
 ! 16 CO
 ! 17 Bidirectional VOC 
 ! 18 Stress VOC 
 ! 19 Other VOC 
 
 ! 1.2.4 Plant functional type number [from CLM4 global land area; Guenther et al. (2012), Tab.3]:
 ! from mod_bvoc: pftn = 15
 
 ! 1  Needleleaf Evergreen Temperate Tree   
 ! 2  Needleleaf Evergreen Boreal Tree 
 ! 3  Needleleaf Deciduous Boreal Tree 
 ! 4  Broadleaf Evergreen Tropical Tree 
 ! 5  Broadleaf Evergreen Temperate Tree 
 ! 6  Broadleaf Deciduous Tropical Tree 
 ! 7  Broadleaf Deciduous Temperate Tree
 ! 8  Broadleaf Deciduous Boreal Tree 
 ! 9  Broadleaf Evergreen Temperate Shrub 
 ! 10 Broadleaf Deciduous Temperate Shrub 
 ! 11 Broadleaf Deciduous Boreal Shrub
 ! 12 Arctic C3 Grass 
 ! 13 Cool C3 Grass 
 ! 14 Warm C4 Grass 
 ! 15 Crop1
 
 ! 1.2.5 Emission factors EF (for BVOC compound  class number and 
 !      egetation type number) [Guenther et al. (2012), Tab.2]:
 ! from mod_bvoc: epsil(pftn,bccn)
 
 ! 1.2.6 Parameters for calculating emission activity factor accounting 
 !       for light response gamma_P [Guenther et al. (2012), Eq. 3]:
 ! ---------------
 ! 1.2.6.1 Emission activity factor accounting for light response gamma_P 
 ! [Guenther et al. (2012), Eq. 3]
    REAL(KIND=wp)                  :: gamma_p
 
 ! 1.2.6.2 Light dependend fraction of the bccn 1-19 LDF [Guenther et al. (2012), Tab.4]:
 ! from mod_bvoc: ldf(bccn)
 
 ! 1.2.6.3 Light-dependent activity factor gamma_P_LDF [Guenther et al. (2012), Eq. 4]
    REAL(KIND=wp)                  :: gamma_pldf
 
 ! 1.2.6.4 Photosynthetic photon flux density PPFD
    REAL(KIND=wp)                  :: ppfd
 
 ! 1.2.6.5 Alpha [Guenther et al. (2012), Eq. 5]
    REAL(KIND=wp)                  :: alpha
 
 ! 1.2.6.6 C_P [Guenther et al. (2012), Eq. 6]
    REAL(KIND=wp)                  :: cp
 
 ! 1.2.6.7 Standard conditions for PPFD averaged over the past 24h P_S [umol m-2 s-1]
 ! (200 for sun leaves, 50 for shade leaves)
 ! [not a constant but set to const owing to lack of data] 
 ! from mod_bvoc: ppfds
 
 ! 1.2.6.8 Average ppfd of the past 24h P_24 [umol m-2 s-1]
 ! [not a constant but set to const owing to lack of data] 
 ! from mod_bvoc: ppfds24
 
 ! 1.2.6.9 Average ppfd of the past 240h P_240 [umol m-2 s-1]
 ! [not a constant but set to const owing to lack of data] 
 ! from mod_bvoc: ppfds240
 
 ! 1.2.7 Parameters for calculating emission activity factor accounting 
 !       for temperature gamma_T [Guenther et al. (2012), Eq. 7]:
 ! ---------------
 ! 1.2.7.1 Emission activity factor accounting for temperature gamma_T
 !        [Guenther et al. (2012), Eq. 7]
    REAL(KIND=wp)                  :: gamma_t
 
 ! 1.2.7.2 Gamma_T_LDF [Guenther et al. (2012), Eq. 8]
    REAL(KIND=wp)                  :: gamma_tldf
 
 ! 1.2.7.3 X
    REAL(KIND=wp)                  :: x
 
 ! 1.2.7.4 Empirically determined coefficients of the bccn 1-19 C_T1 [Guenther et al.(2012),Tab.4]:
 ! from mod_bvoc: ct1(bccn)
 
 ! 1.2.7.5 Empirically determined coefficient C_T2 
 ! from mod_bvoc: ct2
 
 ! 1.2.7.6 T_opt [Guenther et al. (2012), Eq. 9]
    REAL(KIND=wp)                  :: topt 
 
 ! 1.2.7.7 E_opt [Guenther et al. (2012), Eq. 10]
    REAL(KIND=wp)                  :: eopt 
 
 ! 1.2.7.8 Leaf temperature at standard conditions T_S [K]
 ! from mod_bvoc: ts
  
 ! 1.2.7.9 Average leave temperature of the past 24h T_24 [K]
 ! [not a constant but set to const owing to lack of data] 
 ! from mod_bvoc: t24
 
 ! 1.2.7.10 Average leave temperature of the past 240h T_240 [K]
 ! [not a constant but set to const owing to lack of data] 
 ! from mod_bvoc: t240
 
 ! 1.2.7.11 Empirically determined coefficients of the bccn 1-19 C_eo
 !          [Guenther et al. (2012), Tab.4]:
 ! from mod_bvoc: ceo(bccn)
 
 ! 1.2.7.12 Gamma_T_LIF [Guenther et al. (2012), Eq. 11]
    REAL(KIND=wp)                  :: gamma_tlif
 
 ! 1.2.7.13 Empirically determined coefficients of the bccn 1-19 beta
 !          [Guenther et al. (2012), Tab.4]:
 ! from mod_bvoc: beta(bccn)
 
 ! 1.2.8 Parameters for calculating emission activity factor accounting 
 !     for leaf age gamma_A [Guenther et al. (2012), Eq. 12]:
 ! ---------------
 ! 1.2.8.1 Emission activity factor accounting for leaf age gamma_A [Guenther et al. (2012),Eq. 12]
    REAL(KIND=wp)                  :: gamma_a
 
 ! 1.2.8.2 Empirically determined coefficients of the bccn 1-19 A_new 
 ! (relative emission rates for new leaves) [Guenther et al. (2012), Tab.4]:
 ! from mod_bvoc: anew(bccn)
 
 ! 1.2.8.3 Empirically determined coefficients of the bccn 1-19 A_gro 
 ! (relative emission rates for growing leaves) [Guenther et al. (2012), Tab.4]:
 ! from mod_bvoc: agro(bccn)
 
 ! 1.2.8.4 Empirically determined coefficients of the bccn 1-19 A_mat
 ! (relative emission rates for mature leaves) [Guenther et al. (2012), Tab.4]:
 ! from mod_bvoc: amat(bccn)
 
 ! 1.2.8.5 Empirically determined coefficients of the bccn 1-19 A_sen
 ! (relative emission rates for senescing leaves) [Guenther et al. (2012), Tab.4]:
 ! from mod_bvoc: asen(bccn)
 
 ! 1.2.8.6 Fraction of new foliage F_new [0 1]
 ! [not a constant but set to const owing to lack of data] 
 ! from mod_bvoc: fnew
 
 ! 1.2.8.7 fraction of growing foliage F_gro [0 1]  
 ! [not a constant but set to const owing to lack of data] 
 ! from mod_bvoc: fgro
 
 ! 1.2.8.8 fraction of mature foliage F_mat [0 1] 
 ! [not a constant but set to const owing to lack of data] 
 ! from mod_bvoc: fmat
 
 ! 1.2.8.9 fraction of senescing foliage F_sen [0 1]
 ! from mod_bvoc: fsen
 
 ! 1.2.9 Parameters for calculating emission activity factor accounting 
 !       for soil moisutre gamma_SM [Guenther et al. (2012), Eq. 13a-c]:
 ! ---------------
 ! 1.2.9.1 Emission activity factor accounting for soil moisture gamma_SM 
 !        [Guenther et al. (2012), Eq. 13]
 !   REAL(KIND=wp)                  :: gamma_sm

   INTEGER :: ipft  !< loop index
 
 ! Maximum rooting depth to estimate relevant soil moisture layer -
 ! estimated from: Canadell, J., Jackson, R. B., Ehleringer, J. R., 
 ! Mooney, H. A., Sala, O. E., Schulze,  E.-D., 1996: Maximum rooting 
 ! depth of vegetation types at the global scale, Oecologia, p. 583-595, 
 ! Appendix 1
 
 ! PFT                                      Estimated max. root depth (m)
 ! 1  Needleleaf Evergreen Temperate Tree   
 ! 2  Needleleaf Evergreen Boreal Tree      2 (boreal forest)
 ! 3  Needleleaf Deciduous Boreal Tree      2 (boreal forest)
 ! 4  Broadleaf Evergreen Tropical Tree     3 (tropical evergreen forest) 
 ! 5  Broadleaf Evergreen Temperate Tree    3 (temperate coniferous forest)
 ! 6  Broadleaf Deciduous Tropical Tree     4 (tropical deciduous forest)
 ! 7  Broadleaf Deciduous Temperate Tree    3 (temperate deciduous forest)
 ! 8  Broadleaf Deciduous Boreal Tree       2 (boreal forest)
 ! 9  Broadleaf Evergreen Temperate Shrub   3 (sclerophyllous shrubland and forest
 ! 10 Broadleaf Deciduous Temperate Shrub   3 (sclerophyllous shrubland and forest 
 ! 11 Broadleaf Deciduous Boreal Shrub      2 (boreal forest)
 ! 12 Arctic C3 Grass                       0.5 (tundra)
 ! 13 Cool C3 Grass                         2 (temperate grassland)
 ! 14 Warm C4 Grass                         2 (temperate grassland)
 ! 15 Crop1                                 2 (crops)
 ! -> Use mean soil moisture in the layer between 0 and approx. 2 m for all PFTs.
 
 ! 1.2.9.2 Empirical determined coefficient delta_theta_1 [m3 m-3]
 ! from mod_bvoc: dtheta1
 
 ! 1.2.9.3 Theta_1 [m3 m-3]
 !   REAL(KIND=wp)                  :: theta1 
 
 ! 1.2.10 Parameters for calculating emission activity factor accounting 
 !        for CO2 inhibition gamma_C [Guenther et al. (2012), Eq. 14]:
 ! ---------------
 ! 1.2.10.1 Emission activity factor accounting for CO2 inhibition gamma_C
 !          [Guenther et al. (2012), Eq. 14];
 ! based on Heald et al.(2009)]
 !   REAL(KIND=wp)                  :: gamma_C
 
 ! 1.2.10.2 Estimated asymptote at which further decreases in intercellular CO2
 ! have a negligible effect on isoprene emission I_Smax [Heald et al.(2009), Tab. 2]
 ! from mod_bvoc: ismax
 
 ! 1.2.10.3 Scaling coefficient C_star [Heald et al.(2009), Tab. 2]
 ! from mod_bvoc: cstar
 
 ! 1.2.10.4 Exponential scalar h [Heald et al.(2009), Tab. 2]
 ! from mod_bvoc: h 
 
 
 ! 2. Data output
 ! ######################################################################
 ! bvoc concentration
 ! in data_cosmo_art: cbio(ie,bccn)
 
 
 ! 3. General idea
 ! ######################################################################
 
 ! F(i) = gamma(i) * sum( epsil(i,j) * pft(j) ) [ug m-2 h-1]
 
 ! with F        emission of chemical species class i
 !      gamma    emission activity factor (~environmental and phenological conditions)
 !      epsil    emission factor at standard conditions
 !      pft      fractional grid box area coverage for plant functional type (PFT)
  
 ! i = 1,2,3,...,19 (bccn)  compound class
 ! j = 1,2,3,...,15 (pftn)  PFT
  
 
 ! 4.0 Initialization
    calc_bvoc=0._wp
 
 ! 4.1 Loop over all PFTs
    DO ipft = 1,pftn              ! plant functional type
 
 ! 4.2 Calculation of emission activity factor gamma
 ! 4.2.1 Emission activity factor accounting for light response (gamma_P)
      cp         = 0.0468_wp * EXP( 0.0005_wp * (ppfd24(czone)-ppfds)) * ((ppfd240(czone))**0.6_wp)
      alpha      = 0.004_wp - 0.0005_wp * LOG(ppfd240(czone))
      gamma_pldf = cp * ( (alpha * ppfd) / ( SQRT(1._wp + alpha*alpha * ppfd*ppfd)) )
 
      gamma_p    = (1._wp-ldf(lbccn)) + ldf(lbccn) * gamma_pldf
 
 ! 4.2.2 Emission activity factor accounting for temperature (gamma_T)
 ! LDF - light dependent fraction  
      topt       = 313._wp + (0.6_wp * (t240(czone)-ts))
      eopt       = ceo(lbccn) * EXP(0.05_wp * (t24(czone)-ts)) * EXP(0.05_wp * (t240(czone)-ts))
      x          = ((1._wp / topt) - (1._wp / temp)) / 0.00831_wp 
      gamma_tldf = eopt * (ct2 * EXP(ct1(lbccn) * x) / (ct2 - ct1(lbccn) * (1._wp - EXP(ct2 * x))))
 
 ! LIF - light independent fraction (follows monoterpene exponential temperature
 !       response function of Guenther et al. (1993))  
      gamma_tlif = EXP(beta(lbccn) * (temp - ts))
 
      gamma_t    = (1._wp - ldf(lbccn)) * gamma_tlif + ldf(lbccn) * gamma_tldf
 
 ! 4.2.3 Emission activity factor accounting for leaf age (gamma_A)
      gamma_a = fnew * anew(lbccn) + fgro * agro(lbccn) + fmat * amat(lbccn) + fsen * asen(lbccn)
 
 ! *** gamma for soil moisture and CO2 up to now not used *** 
 !!  If current bccn is isoprene  
 !      IF (lbccn.eq.1) THEN
 !
 !! 4.2.4 Emission activity factor accounting for soil moisture (only isoprene) (gamma_SM)
 !      theta1 = thetaw + dtheta1
 !
 !        IF (theta.gt.theta1) THEN
 !          gamma_sm = 1._wp
 !
 !        ELSE IF (theta.lt.theta1 .and. theta.gt.thetaw) THEN
 !          gamma_sm = (theta - thetaw) / dtheta1  
 !   
 !        ELSE IF (theta.lt.thetaw) THEN
 !          gamma_sm = 0._wp
 !
 !        ELSE IF (theta.eq.0._wp .and. thetaw.eq.0._wp) THEN
 !          write(*,*) "WARNING in subroutine art_bvoc: both soil moisture and plant &
 !            &         wilting point are zero."
 !          gamma_sm = 0._wp
 !        ENDIF
 !
 !! 4.2.5 Emission activity factor accounting for CO2 inhibition (only isoprene) (gamma_C)
 !! (CO2 inhibition parameterization is optimized for an ambient CO2 range of 365-717 ppm)
 !
 !        gamma_c = ismax - (ismax * ((co2*0.7_wp)**h)) / ((cstar**h) + ((co2*0.7_wp)**h))
 !
 !      ELSE  
 !!  if current bccn is not isoprene (all other BVOC compound classes) 
 !        gamma_sm = 1._wp 
 !        gamma_c  = 1._wp 
 !      ENDIF  
 ! ***
 
 ! 4.3 Calculate gamma by using 4.2
 ! (gamma_P,gamma_T,gamma_A,gamma_SM and gamma_C)
 !      gamma = cce * leafai * gamma_P * gamma_T * gamma_A * gamma_SM * gamma_C
 !                        ! gamma for soil moisture and CO2 up to now not used
      gamma = cce * leafai(ipft) * gamma_P * gamma_T * gamma_A
 
 ! 4.4 Calculation of biogenic emissions
      calc_bvoc = calc_bvoc + (gamma * epsil(ipft,lbccn) * pfts(ipft)) ! [ug m-2 h-1]
 
    ENDDO ! ipft
 
 ! 4.5 Unit conversion ug m-2 h-1 -> kg m-2 h-1
    calc_bvoc = calc_bvoc/1.E09_wp
 
  END FUNCTION calc_bvoc
 
END SUBROUTINE bvoc_guenther2012

END MODULE mo_art_bvoc
