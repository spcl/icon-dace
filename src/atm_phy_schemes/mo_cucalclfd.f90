!
!**** *CUCALCLFD*  ROUTINE FOR LFD COMPUTATION
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

MODULE mo_cucalclfd
  
  USE mo_kind   ,ONLY: JPRB=>wp     , &
    &                  jpim=>i4
  
  USE mo_cuparameters , ONLY :                                   &
    & rg, rd, rcpd
  USE mo_cufunctions, ONLY: foealfa, foeldcpm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cucalclfd

CONTAINS

SUBROUTINE cucalclfd(klon, klev, ktype, ztu, zlu, kcbot, zcape, zmflxs &
      &            , zten, pap, zdgeoh, zgeoh, ldland, lfd, lacc)

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
!    *KCBOT*        CLOUD BASE LEVEL
!    *ZCAPE*        CONVECTVE AVAILABLE POTENTIAL ENERGY           J/KG
!    *ZMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *ZTEN*         ENVIRONMENT TEMPERATURE (T+1)                 K
!    *PAP*          PRESSURE ON FULL LEVELS                       PA
!    *zdgeoh*       geopot thickness on full levels               M2/S2
!    *zGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2

!    OUTPUT PARAMETERS (REAL):

!    *LFD*          LIGHTNING FLASH DENSITY AS IN LOPEZ (2016)    1/DAY/km2

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
REAL(KIND=jprb)   ,INTENT(in)  :: zgeoh(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(in)  :: zcape(klon) 
INTEGER(KIND=jpim),INTENT(in)  :: kcbot(klon)
LOGICAL           ,INTENT(in)  :: ldland(klon) 
LOGICAL           ,INTENT(in)  :: lacc
REAL(KIND=jprb)   ,INTENT(out) :: lfd(klon)

REAL(KIND=jprb)   :: zrho(klon,klev)
INTEGER (KIND=jpim)   :: kland(klon)
REAL(KIND=jprb)   :: zqr(klon)

! Some consants as in Lopez 2016
!--------------------------------
! coefficient to split solid water mass flux in snow and graupel
! beta(1) for land, beta(2) for sea - Takahashi 2006 (mentioned
! in Lopez 2016)
REAL(KIND=jprb), PARAMETER :: beta(2)= [0.7_jprb , 0.45_jprb ]
! Causes currently internal compiler error for Cray, ifdef should be removed when fixed
#ifndef _CRAYFTN 
  !$ACC DECLARE COPYIN(beta)
#endif
! some constants from Lopez 2016
REAL(KIND=jprb), PARAMETER :: Vgraup=3.0_jprb ! 3   m/s fall speed for graupel
REAL(KIND=jprb), PARAMETER :: Vsnow=0.5_jprb  ! 0.5 m/s fall speed for snow
REAL(KIND=jprb), PARAMETER :: alpha = 32.4_jprb

! Auxiliary variables as in Lopez 2016
REAL(KIND=jprb) :: zqIce, zqLiquid, zqGraup, zqSnow, zEps, zQi, zdz

INTEGER(KIND=jpim) :: jk, jl

  !$ACC DATA &
  !$ACC   PRESENT(ktype, ztu, zlu, zmflxs, zten, pap, zdgeoh, zgeoh, zcape, kcbot, ldland, lfd) &

  !$ACC   CREATE(zrho, zqr, kland) &
  !$ACC   IF(lacc)

  !US we could use default (none) here, but PGI complains about beta then!?
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

  ! Make land sea-mask for beta
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl = 1, klon
    kland(jl) = MERGE (1_jpim, 2_jpim, ldland(jl))
    zQR  (jl) = 0.0_jprb
    lfd  (jl) = 0.0_jprb
  ENDDO

  !$ACC LOOP SEQ
  DO jk = 1, klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO JL = 1, klon
      zrho(JL,JK)=pap(JL,JK)/(zten(JL,JK)*Rd+1E-10_jprb)
    ENDDO
  ENDDO

  ! Compute the LFD (Lopez 2016)
  !-----------------------------

  !$ACC LOOP SEQ
  DO jk = 1, klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqGraup, zqSnow, zdz)
    DO JL = 1, klon
      ! only for deep convection LPI is computed
      IF (ktype(JL) == 1) THEN
        IF (ztu(JL,JK) <=273.15_jprb .AND. ztu(JL,JK) >=273.15_jprb-25._jprb) THEN
        ! compute graupel and snow mixing ratio from the massflux
        ! of frozen precip - splitting it into snow and graupel with beta.
        ! See Eq. 1 and 2 in Lopez (2016) 
        ! "A lightning parameterization for the ECMWF integrated forecasting system"
          zqGraup=beta(kland(JL))*zmflxs(JL, JK)/zrho(JL,JK)/Vgraup
          zqSnow=(1_jprb-beta(kland(JL)))*zmflxs(JL,JK)/zrho(JL,JK)/Vsnow
          zdz=zdgeoh(jl,jk)/rg
        ! Eq. 3 in Lopez (2016) - integral
          zQR(JL)=zQR(JL) + zqGraup*(zlu(JL,JK)+zqSnow)*zrho(JL,JK)*zdz
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl = 1, klon
      ! Eq. 5 in Lopez (2016)
      ! Cloud base height is zgeoh(JL, kcbot(JL))/rg
    IF (ktype(JL) == 1_jpim .AND. kcbot(JL) > 0_jpim) THEN
        lfd(JL)=alpha*zQR(JL)*SQRT(MAX(zCAPE(JL),0._jprb))*             &
 &       MIN(zgeoh(JL, kcbot(JL))/rg/1000._jprb, 1.8_jprb)**2  
    ENDIF
  ENDDO

  !$ACC END PARALLEL
  !$ACC WAIT(1)

  !$ACC END DATA

END SUBROUTINE cucalclfd
 
END MODULE mo_cucalclfd

