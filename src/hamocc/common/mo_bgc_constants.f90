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

MODULE mo_bgc_constants

!- Description:
!
!  This module contains constants and parameters
!  

  USE mo_kind, ONLY: wp

 
  IMPLICIT NONE

  PUBLIC
   

   REAL(wp), PARAMETER :: PI = 3.141592654_wp 
   REAL(wp), PARAMETER :: tmelt = 273.15_wp        ! melting temperature of ice/snow
   !$ACC DECLARE COPYIN(tmelt)
   
   ! Conversion factors
   REAL(wp), PARAMETER:: kilo = 1.e3_wp
   REAL(wp), PARAMETER:: s2year = 3.1536e7_wp ! seconds to years
   REAL(wp), PARAMETER:: c2gtc = 12._wp*1.e-12_wp ! kmolC to GtC 
   REAL(wp), PARAMETER:: n2tgn = 14._wp*1e-9_wp ! kmolN to TgN
  !>  GRAVITATIONAL ACCELERATION (CONSTANT 9.81 M/S**2)
   REAL(wp), PARAMETER :: g = 9.81_wp
  !>  OCEAN REFERENCE DENSITY [kg m-3]
  REAL(wp), PARAMETER :: rhoref_water = 1025.0_wp
  REAL(wp), PARAMETER:: rmnit=14._wp
  REAL(wp), PARAMETER:: cmh2ms = 6.9722e-07_wp ! cm/hr to m/sec * a constant for piston velocity

  REAL(wp), PARAMETER :: molw_co2=44.011_wp

  ! Parameters used in MAGO sinking
  ! molecular dynamic viscosity
  REAL(wp),PARAMETER, DIMENSION(14) :: cdynv=(/1.79e-2_wp, -6.1299e-4_wp,1.4467e-5_wp, &
                           &      -1.6826e-7_wp, -1.8266e-7_wp, 9.8972e-12_wp, 2.4727e-5_wp,&
                           &      4.8429e-7_wp, -4.7172e-8_wp,7.5986e-10_wp,1.3817e-8_wp,&
                           &      -2.6363e-10_wp,6.3255e-13_wp,-1.2116e-14_wp/)


  !---------------PARAMETERS USED in CHEMCON

  !*         HALF PRECISION CONSTANTS
  REAL(wp), PARAMETER ::   SMICR = 1.E-6_wp, PERC  = 1.E-2_wp, THIRD = 1._wp/3._wp
  !*         CONVERSION FACTOR SALINITY -> CHLORINITY
  !             (AFTER WOOSTER ET AL., 1969)
  REAL(wp), PARAMETER ::  SALCHL = 1._wp/1.80655_wp
  !*         COEFFICIENTS FOR APPARENT SOLUBILITY EQUILIBRIUM
  !             OF CALCITE (Mucci 1983, in Z & W-G,2001)
  !          log10ksp = akcc1 + akcc2*t + akcc3/t + akcc4*LOG10(t)         &
  !               + (akcc5 +akcc6*t + akcc7/t) * sqrt(s)                  &
  !               + akcc8*s + akcc9*(s**(1.5))
  REAL(wp), PARAMETER ::  akcc1 = -171.9065_wp,   akcc2 = -0.077993_wp, &
       akcc3 = 2839.319_wp, akcc4 = 71.595_wp, akcc5 = -0.77712_wp,     &
       akcc6 = 0.0028426_wp, akcc7 = 178.34_wp, akcc8 = -0.07711_wp,    &
       akcc9 = 0.0041249_wp
  !*        FRACTION OF ARAGONITE IN BIOGENIC CACO3 PARTICLES
  !             --- -------- --------- -- -------- ----- ---------
  REAL(wp), PARAMETER :: arafra = 0._wp  ! simplification to use solubility product for calcite and
  !               ! adjust it in case of aragonit production
  !*         FRACTION OF CALCITE IN BIOGENIC CACO3 PARTICLES
  REAL(wp), PARAMETER ::   calfra = 1._wp - arafra
  !*          FACTOR TO GET APPARENT SOLUBILITY PRODUCT FOR
  !            ARAGONIT BY MULTIPLICATION WITH APPARENT SOLUBILITY
  !            PRODUCT FOR CALCIT (CF. BERNER, 1976,
  !            OR BROECKER ET AL., 1982)
  REAL(wp), PARAMETER :: aracal = arafra * 1.45_wp + calfra   ! simplification instead of calculating
                                       ! solubility product for aragonit
  !*           COEFFICIENTS FOR PRESSURE CORRECTION OF SOLUBILITY
  !            PRODUCT OF CACO3 CORRESPONDING TO ARAGONITE/CALCITE
  !            RATIO
  !             Millero (1995)
  !             this formulation is  only for molal volume changes,
  !             it neglects compressibility change
  !
  REAL(wp), PARAMETER ::  devkst = 0.5304_wp
  REAL(wp), PARAMETER ::  devks  = 46_wp * arafra + 48.76_wp * calfra

  !*        GAS CONSTANT  for pressure correction after Millero(1995)
  !             in cm**3 bar / mol K (from Zeebe and Wolf-Gladrow)
  REAL(wp), PARAMETER ::   rgas = 83.131_wp
  !*         BORON CONCENTRATION IN SEA WATER
  !*         IN G/KG PER O/OO CL ACCORDING
  !          TO RILEY AND SKIRROW, 1965 (P.250)
  REAL(wp), PARAMETER ::   bor1 = 0.00023_wp
  !*        INVERSE OF ATOMIC WEIGHT OF BORON [G**-1]
  !         (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
  REAL(wp), PARAMETER ::   bor2 = 1._wp/10.82_wp
  !*        COEFF. FOR 1. DISSOC. OF CARBONIC ACID
  !          Dickson, 2010
  REAL(wp), PARAMETER ::  c10 = -3633.86_wp, c11 = 61.2172_wp,&
    c12 = -9.67770_wp,  c13 = 0.011555_wp,   c14 = -0.0001152_wp
  !*         COEFF. FOR 2. DISSOC. OF CARBONIC ACID
  !          Dickson, 2010
  REAL(wp), PARAMETER :: c20 = -471.78_wp, c21 = -25.9290_wp,&
     c22 = 3.16967_wp,c23 = 0.0178_wp,c24 =-0.0001122_wp

 !       coefficients for pressure correction of
 !      aks, akf, ak1p, ak2p, ak3p, ksi, k1, k2, kb, kw, aksp 
  REAL(wp), PARAMETER,DIMENSION(11):: pa0= (/ -18.03_wp, -9.78_wp,  -14.51_wp, -23.12_wp, -26.57_wp, &
&          -29.48_wp, -25.5_wp, -15.82_wp, -29.48_wp, -20.02_wp, -45.96_wp/)
  REAL(wp), PARAMETER,DIMENSION(11):: pa1= (/ 0.0466_wp, -0.0090_wp, 0.1211_wp, 0.1758_wp, 0.2020_wp, &
&          0.1622_wp,  0.1271_wp, -0.0219_wp, 0.1622_wp, 0.1119_wp, 0.5304_wp/)
  REAL(wp), PARAMETER, DIMENSION(11)::pa2 =(/ 0.316e-3_wp, -0.942e-3_wp, -0.321e-3_wp, -2.647e-3_wp, &
&         -3.042e-3_wp, -2.6080e-3_wp,  0.0_wp,  0.0_wp, -2.608e-3_wp,-1.409e-3_wp, 0.0_wp /)
  REAL(wp), PARAMETER, DIMENSION(11)::pb0 = (/ -4.53e-3_wp,  -3.91e-3_wp,  -2.67e-3_wp, -5.15e-3_wp,-4.08e-3_wp, &
&           -2.84e-3_wp,  -3.08e-3_wp,   1.13e-3_wp, -2.84e-3_wp, -5.13e-3_wp,-11.76e-3_wp  /)
  REAL(wp), PARAMETER, DIMENSION(11)::pb1 = (/ 0.09e-3_wp,  0.054e-3_wp,  0.0427e-3_wp,  0.09e-3_wp, 0.0714e-3_wp,&
&           0.0_wp,     0.0877e-3_wp, -0.1475e-3_wp, 0.0_wp, 0.0794e-3_wp,-0.3692e-3_wp/)
  REAL(wp), PARAMETER,DIMENSION(11):: pb2 = 0.0_wp
 !        coeff. for silicic acid dissociation
 !        Millero p.671 (1995) using data from Yao and Millero (1995)
  REAL(wp), PARAMETER:: cksi1 =-8904.2_wp,  cksi2 = 117.385_wp, &
     cksi3  = -19.334_wp , cksi4  =   -458.79_wp,  cksi5  =    3.5913_wp,&
     cksi6  =    188.74_wp, cksi7  =  - 1.5998_wp, cksi8  =  -12.1652_wp,& 
     cksi9  =   0.07871_wp,cksi10 =  0.001005_wp
 !     coeff. for hydrogen sulfate dissociation
 !     Dickson (1990, J. chem. Thermodynamics 22, 113)
 REAL(wp),PARAMETER:: cks1  =   -4276.1_wp, cks2  =   141.328_wp,& 
   cks3  =  - 23.093_wp, cks4  =   -13856._wp,  cks5  =    324.57_wp, &
   cks6  =  - 47.986_wp, cks7  =    35474._wp,  cks8  =  - 771.54_wp,&
   cks9  =   114.723_wp, cks10 =   - 2698._wp,  cks11 =     1776._wp,&
   cks12 = -0.001005_wp  
  !      coeff. for hydrogen fluoride dissociation
  !      Dickson and Riley (1979) 
  REAL(wp), PARAMETER:: ckf1 = 874._wp, ckf2 = -9.68_wp,  ckf3 = 0.111_wp
  !      coeff. for phosphoric acid dissociation
  !      ak1p
  REAL(wp), PARAMETER:: ck1p1 = -4576.752_wp, ck1p2 =   115.525_wp, &
      ck1p3 =  - 18.453_wp, ck1p4 =  -106.736_wp, ck1p5 =   0.69171_wp ,& 
      ck1p6 =  -0.65643_wp,   ck1p7 =  -0.01844_wp 
  !      ak2p
  REAL(wp), PARAMETER:: ck2p1 = -8814.715_wp , ck2p2 =  172.0883_wp, &
   ck2p3 =  - 27.927_wp, ck2p4 =  -160.340_wp,  ck2p5 =    1.3566_wp, & 
   ck2p6 =   0.37335_wp,  ck2p7 = - 0.05778_wp 
  !      ak3p
  REAL(wp), PARAMETER:: ck3p1 =  -3070.75_wp, ck3p2 =  - 18.141_wp, &
   ck3p3 =  17.27039_wp, ck3p4 =   2.81197_wp, ck3p5 = -44.99486_wp, &
   ck3p6 = - 0.09984_wp
  !*        COEFF. FOR 1. DISSOC. OF BORIC ACID
  !         DOE (1994) (based on by Dickson 1990)
  REAL(wp), PARAMETER :: cb0 = -8966.90_wp,  cb1 = -2890.53_wp,  cb2 = -77.942_wp, &
       cb3 = 1.728_wp, cb4 = -0.0996_wp, cb5 = 148.0248_wp, cb6 = 137.1942_wp,      &
       cb7 = 1.62142_wp, cb8 = 24.4344_wp, cb9 = 25.085_wp, cb10 =0.2474_wp,        &
       cb11 =0.053105_wp
  !*     SET COEFF. FOR ION PRODUCT OF WATER
  !      DOE (1994)
  REAL(wp), PARAMETER ::  cw0 = 148.9802_wp, cw1 = -13847.26_wp, &
       cw2 = -23.6521_wp, cw3 = 118.67_wp, cw4 = -5.977_wp,       &
       cw5 = 1.0495_wp, cw6 = -0.01615_wp
  !*         VOLUMETRIC SOLUBILITY CONSTANTS FOR CO2 IN ML/L
  !          WEISS, R. F. (1974)
  !          CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A
  !           NON IDEAL GAS. MARINE CHEMISTRY, VOL. 2, 203-215.
  REAL(wp), PARAMETER :: c00 = 9345.17_wp,  c01 = -60.2409_wp, c02 = 23.3585_wp, &
       c03 = 0.023517_wp, c04 = -0.00023656_wp, c05 = 0.0047036_wp
  !*     VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L
  !         (WEISS, 1970)
  REAL(wp), PARAMETER ::  ox0 = -173.4292_wp, ox1 = 249.6339_wp, ox2 = 143.3483_wp, &
       ox3 = -21.8492_wp, ox4 = -0.033096_wp, ox5 = 0.014259_wp, ox6 = -0.0017_wp
  !*          VOLUMETRIC SOLUBILITY CONSTANTS FOR N2 IN ML/L
  !           WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN
  !           OXYGEN AND ARGON IN WATER AND SEAWATER.
  !           DEEP-SEA RESEARCH, VOL. 17, 721-735.
  REAL(wp), PARAMETER ::  an0 = -172.4965_wp, an1 = 248.4262_wp, an2 = 143.0738_wp, &
       an3 = -21.7120_wp, an4 = -0.049781_wp, an5 = 0.025018_wp, an6 = -0.0034861_wp
  !      Constants for laughing gas solubility
  !      (WEISS, 1974, MARINE CHEMISTRY)
  REAL, PARAMETER:: a1 = -165.8806_wp, a2 =  222.8743_wp, a3 = 92.0792_wp,&
      a4 = -1.48425_wp,  b1 = -0.056235_wp,  b2 = 0.031619_wp, b3 = -0.0048472_wp

  REAL(wp), PARAMETER :: atn2o = 3.e-7_wp

  REAL(wp), PARAMETER :: ten   = 10._wp
  !  ml/L to kmol/m3
  REAL(wp), PARAMETER :: oxyco = 1._wp/22391.6_wp
  !
  REAL(wp), PARAMETER :: invlog10 = 1.0_wp/log(10.0_wp)


 ! Parameters for oxygen saturation concentration req. in OMIP
 ! following Garcia& Gordon 1992, coefficients from
 ! Tab1, col2 Benson & Krause
  REAL(wp), PARAMETER :: oxya0=2.00907_wp, oxya1=3.22014_wp, oxya2=4.05010_wp,  &
  &     oxya3=4.94457_wp, oxya4=-2.56847E-1_wp, oxya5=3.88767_wp, &
  &     oxyb0=-6.24523E-3_wp, oxyb1=-7.37614E-3_wp, oxyb2=-1.03410E-2_wp, oxyb3=-8.17083E-3_wp, &
  &     oxyc0=-4.88682E-7_wp


END MODULE 
