!
!**** *CUCALCLPI*  ROUTINE FOR LPI COMPUTATION
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

MODULE mo_cucalclpi
  
  USE mo_kind   ,ONLY: JPRB=>wp     , &
    &                  jpim=>i4
  
  USE mo_cuparameters , ONLY :                                   &
    & rg, rd, rcpd
  
  USE mo_cufunctions, ONLY: foealfa, foeldcpm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cucalclpi, cucalcmlpi

CONTAINS

SUBROUTINE cucalclpi(klon, klev, ktype, ztu, zlu, zkineu, zmflxs        &
      &            , zten, pap, zdgeoh, ldland, lpi, lacc)

! Code Description:
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION

!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)

!     INPUT PARAMETERS (REAL)

!    *ZTU*          TEMPERATURE IN UPDRAFTS                         K
!    *ZLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *ZKINEU*       KINETIC ENERGY IN UPDRATFS                    M2/S2
!    *ZMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *ZTEN*         ENVIRONMENT TEMPERATURE (T+1)                 K
!    *PAP*          PRESSURE ON FULL LEVELS                       PA
!    *zdgeoh*       geopot thickness on full levels               M2/S2

!    OUTPUT PARAMETERS (REAL):

!    *LPI*          LIGHTNING POTENTIAL INDEX AS IN LYNN AND YAIR (2010) J/KG

IMPLICIT NONE

INTEGER(KIND=jpim),INTENT(in)  :: klon
INTEGER(KIND=jpim),INTENT(in)  :: klev
INTEGER(KIND=jpim),INTENT(in)  :: ktype(klon)
REAL(KIND=jprb)   ,INTENT(in)  :: ztu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)  :: zlu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)  :: zmflxs(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(in)  :: zten(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)  :: pap(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)  :: zdgeoh(klon,klev)
REAL(KIND=jprb)   ,INTENT(in)  :: zkineu(klon,klev)
LOGICAL           ,INTENT(in)  :: ldland(klon) 
LOGICAL           ,INTENT(in)  :: lacc
REAL(KIND=jprb)   ,INTENT(out) :: lpi(klon)

REAL(KIND=jprb)   :: zrho(klon,klev)
INTEGER (KIND=jpim)   :: kland(klon)

! Parameters as in Lopez 2016
! coefficient to split solid water mass flux in snow and graupel
! beta(1) for land, beta(2) for sea - Takahashi 2006 (mentioned
! in Lopez 2016)
REAL(KIND=jprb), PARAMETER :: beta(2)= [0.7_jprb , 0.45_jprb ]  
! Causes currently internal compiler error for Cray, ifdef should be removed when fixed
#ifndef _CRAYFTN
  !$ACC DECLARE COPYIN(beta)
#endif
REAL(KIND=jprb), PARAMETER :: Vgraup=3.0_jprb ! 3   m/s fall speed for graupel
REAL(KIND=jprb), PARAMETER :: Vsnow=0.5_jprb  ! 0.5 m/s fall speed for snow

! Auxiliary fields as in Lopez 2016 and Lynn and Yair 2010
REAL(KIND=jprb) :: zqIce, zqLiquid, zqGraup, zqSnow, zEps, zQi, zdz

REAL(KIND=jprb) :: unitVolume(klon) ! the volume used for the integral - region
                                  ! beween 0 and -20 celsius

INTEGER(KIND=jpim) :: jk, jl

  !$ACC DATA &
  !$ACC   PRESENT(ktype, ztu, zlu, zmflxs, zten, pap, zdgeoh, zkineu, ldland, lpi) &

  !$ACC   CREATE(unitVolume, zrho, kland) &
  !$ACC   IF(lacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

  ! Make land sea mask for beta
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl = 1, klon
    kland(jl)      = MERGE (1_jpim, 2_jpim, ldland(jl))
    LPI(jl)        = 0.0_jprb
    unitVolume(jl) = 0.0_jprb
  ENDDO

  !$ACC LOOP SEQ
  DO jk = 1, klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl = 1, klon
      zrho(JL,JK)=pap(JL,JK)/(zten(JL,JK)*Rd+1E-10_jprb)
    ENDDO
  ENDDO

  !$ACC LOOP SEQ
  DO jk = 1, klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqGraup, zqSnow, zqLiquid, zqIce, zQI, zeps, zdz)
    DO jl = 1, klon
      ! only for deep convection LPI is computed
      IF (ktype(JL) == 1) THEN
        IF (ztu(JL,JK) <=273.15_jprb .AND. ztu(JL,JK) >=273.15_jprb-20._jprb) THEN
        ! compute graupel and snow mixing ratio from the massflux
        ! of frozen precip - splitting it into snow and graupel with beta.
        ! See Eq. 1 and 2 in Lopez (2016) 
        ! "A lightning parameterization for the ECMWF integrated forecasting system"
          zqGraup=beta(kland(JL))*zmflxs(JL, JK)/zrho(JL,JK)/Vgraup
          zqSnow=(1-beta(kland(JL)))*zmflxs(JL,JK)/zrho(JL,JK)/Vsnow
          zqLiquid=zlu(JL,JK)*foealfa(ztu(JL,JK))
          zqIce=zlu(JL,JK)*(1._jprb-foealfa(ztu(JL,JK)))
        ! eq. 3 in Lynn and Yair (2010)
          zQI=zqGraup*(SQRT(zqSnow*zqGraup)/(zqSnow+zqGraup+1E-20_jprb)+&
                    SQRT(zqIce*zqGraup)/(zqIce+zqGraup+1E-20_jprb))
        ! eq. 2 in Lynn and Yair (2010)
          zeps=2._jprb*sqrt(zQI*zqLiquid)/(zQi+zqLiquid+1E-20_jprb)
          zdz=zdgeoh(jl,jk)/rg
          LPI(JL) = LPI(JL)+zeps*2._jprb*zkineu(JL, JK)*zdz
          unitVolume(JL) = unitVolume(JL)+zdz
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl = 1, klon
    lpi(JL)=lpi(JL)/(unitVolume(JL)+1E-20_jprb)
  ENDDO

!!  ! We use 15 km as vertial unit volume - after all it is just a scaling factor
!!  ! Note - In the paper the LPI is an average oer the unit volume. If we
!!  ! applied that here properly we would need to multiply with sigma
!!  ! (the updraft area proportion in the cloud). That would mean the numbers
!!  ! get smaller the larger the grid spacing. That is not necessarily
!!  ! desirable.
!!  ! In conclusion: The LPI here is representative for the updraft in the cloud
!!  ! only.
!!  lpi=lpi/15000._jprb 

  !$ACC END PARALLEL
  !$ACC WAIT(1)

  !$ACC END DATA

END SUBROUTINE cucalclpi
 
SUBROUTINE CUCALCMLPI(klon, klev, lpi, zten, zqen, pap, paph, koi, mlpi, lacc)
! Code Description:
! Computes a modified LPI using LPI as in Lynn and Yair and KOI.
!
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION


!     INPUT PARAMETERS (REAL)

!    *LPI*          LIGHTNING POTENTIAL INDEX AS IN LYNN AND YAIR (2010) J/KG
!    *ZTEN*         ENVIRONMENT TEMPERATURE (T+1)                 K
!    *ZQEN*         ENVIRONMENT SPEC. HUMIDITY (T+1)              KG/KG
!    *PAP*          PRESSURE ON FULL LEVELS                       PA
!    *PAPH*         PRESSURE ON HALF LEVELS                       PA

!    OUTPUT PARAMETERS (REAL):

!    *KOI*          CONVECTION INDEX (VERT. GRADIENT OF EQUIALENT POT.TEMP) K
!    *MLPI*         MODIFIED LIGHTNING POTENTIAL INDEX (COMBINATIN WITH KOI) J/KG

IMPLICIT NONE
  
INTEGER(KIND=jpim),INTENT(in) :: klon
INTEGER(KIND=jpim),INTENT(in) :: klev
REAL(KIND=jprb),INTENT(in)    :: zten(klon,klev) 
REAL(KIND=jprb),INTENT(in)    :: zqen(klon,klev) 
REAL(KIND=jprb),INTENT(in)    :: pap(klon,klev)
REAL(KIND=jprb),INTENT(in)    :: paph(klon,klev+1)
REAL(KIND=jprb),INTENT(in)    :: lpi(klon)
LOGICAL        ,INTENT(in)    :: lacc
REAL(KIND=jprb),INTENT(out)   :: mlpi(klon)
REAL(KIND=jprb),INTENT(out)   :: koi(klon)
! Equivalent temperatures in 600 and 900 hPa
! Note - the average of 500 to 700 hPa is computed for 600hPa
!        the average below 800 hPa is computed for 900 hPa
REAL(KIND=jprb) :: te
! Equivalent potential temperatures in 600 and 900 hPa
REAL(KIND=jprb) :: thetae600(klon), thetae900(klon), thetae(klon, klev)
! The total pressure difference for the integral
REAL(KIND=jprb) :: deltap600(klon), deltap900(klon)

!* Coefficients for the MLPI formula
!* They were computed in an optimization process
REAL(KIND=jprb) :: fa(klon),fb(klon)
REAL(KIND=jprb), PARAMETER :: fe=0.2960515_jprb
REAL(KIND=jprb), PARAMETER :: fd=4.548663_jprb
REAL(KIND=jprb), PARAMETER :: fg=15.52337_jprb 
REAL(KIND=jprb), PARAMETER :: fh=0.3845962_jprb
REAL(KIND=jprb), PARAMETER :: fi=0.04240491_jprb
REAL(KIND=jprb), PARAMETER :: fj=1.709239_jprb
REAL(KIND=jprb), PARAMETER :: LPIconst=86.15723_jprb

INTEGER(KIND=jpim) :: jk, jl

  !$ACC DATA &
  !$ACC   PRESENT(zten, zqen, pap, paph, lpi, mlpi, koi) &

  !$ACC   CREATE(thetae600, thetae900, thetae, deltap600, deltap900, fa, fb) &
  !$ACC   IF(lacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

  ! compute equivalent potential temperature
  !$ACC LOOP SEQ
  DO jk = 1_jpim, klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(te)
    DO jl = 1_jpim, klon
       te = zten(JL,JK)+foeldcpm(zten(JL,JK)+1E-20)*zqen(JL,JK)
       thetae(JL,JK)=te*(100000.D0/(pap(JL,JK)+1E-10))**(rd/rcpd)
    ENDDO
  ENDDO

  ! Now compute KOI
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl = 1, klon
    thetae900(jl)=0.0_jprb
    deltap900(jl)=0.0_jprb
    thetae600(jl)=0.0_jprb
    deltap600(jl)=0.0_jprb
  ENDDO

  !$ACC LOOP SEQ
  DO jk = 1_jpim, klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl = 1_jpim, klon
      IF (pap(JL,JK) > 80000.D0) THEN
        thetae900(JL) = thetae900(JL)                                   &
 &                       +(paph(JL,JK+1)-paph(JL,JK))*thetae(JL,JK)
        deltap900(JL) =  deltap900(JL)+(paph(JL,JK+1)-paph(JL,JK))
      ELSE IF (pap(JL,JK) < 70000.D0 .AND. pap(JL,JK) > 50000.D0)  THEN
        thetae600(JL) = thetae600(JL)                                   &
 &                       +(paph(JL,JK+1)-paph(JL,JK))*thetae(JL,JK)
        deltap600(JL) =  deltap600(JL)+(paph(JL,JK+1)-paph(JL,JK))
      ENDIF
    ENDDO
  ENDDO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl = 1, klon
    thetae900(jl) = thetae900(jl)/(deltap900(jl)+1E-20)
  ! Over mountains where the lowest pressure level is below 800hPa, 
  ! we use the surface value for thetae
    IF (pap(jl,klev) <= 80000.D0) THEN
      thetae900(jl)=thetae(jl,klev)
    ENDIF
    thetae600(jl) = thetae600(jl)/(deltap600(jl)+1E-20)
  ! same for the higher level thetae - above 500hPa height, KOI
  ! will be zero - and thetae600=thetae900=pap(:,klev).
    IF (pap(jl,klev) <= 50000.D0) THEN
      thetae600(jl)=thetae(jl,klev)
    ENDIF
    KOI(jl)=thetae600(jl)-thetae900(jl)
  ENDDO

  ! Compute the modified LPI
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl = 1, klon
    fa(jl) = fg*LPI(jl)**fh
    ! we require a >= b
    fb(jl) = min(fa(jl), fi*LPI(jl)**fj)
    MLPI(jl)= fb(jl)+(1_jprb+tanh(-fe*(KOI(jl)+fd)))/2_jprb*(fa(jl)-fb(jl))
  ! For LPI > LPIconst, we require MLPI=LPI.
  ! Note that the function is fit in such a way,
  ! that MLPI=fa=fb when LPI=LPIconst.
    IF (LPI(JL) >= LPIconst) THEN
      MLPI(JL)=LPI(JL)
    ENDIF
  ENDDO

  !$ACC END PARALLEL
  !$ACC WAIT(1)
                     
  !$ACC END DATA

END SUBROUTINE CUCALCMLPI

END MODULE mo_cucalclpi
