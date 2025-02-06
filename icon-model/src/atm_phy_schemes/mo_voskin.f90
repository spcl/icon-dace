!
! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! This file has been modified for the use in ICON
!-------------------------------------------------------------------------------
! SPDX-License-Identifier: Apache-2.0
!-------------------------------------------------------------------------------
!
! COMPUTES WARM AND COLD SKIN EFFECTS OVER THE OCEAN

MODULE mo_voskin
 
  PUBLIC :: voskin

CONTAINS

SUBROUTINE VOSKIN(KIDIA,KFDIA,KLON,&
 & PTMST,&
 & PSSRFL ,PSLRFL ,PAHFS, PAHFL, PUSTR, PVSTR, &
 & PU10,PV10,PTSKM1M,PSST,&
 & PTSK )  
!     ------------------------------------------------------------------

!**   *VOSKIN* - COMPUTES WARM AND COLD SKIN EFFECTS OVER THE OCEAN

!     Original A. Beljaars       E.C.M.W.F.         02-01-2001

!     PURPOSE
!     -------

!     SOLVES FOR OCEAN SURFACE TEMPERATURE COLD SKIN AND WARM LAYER

!     INTERFACE
!     ---------

!     *VOSKIN* IS CALLED BY *CALLPAR*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PTMST*        PHYSICS TIME STEP
!     *PSSRFL*       NET SOLAR RADIATION AT THE SURFACE
!     *PSLRFL*       NET THERMAL RADIATION AT THE SURFACE
!     *PAHFS*        SURFACE SENSIBLE HEAT FLUX
!     *PAHFL*        SURFACE LATENT HEAT FLUX
!     *PUSTR*        SURFACE STRESS X-DIRECTION
!     *PVSTR*        SURFACE STRESS Y-DIRECTION
!     *PU10*         10 M WIND X-COMPONENT
!     *PV10*         10 M WIND Y-COMPONENT
!     *PTSKM1M*      SKIN TEMPERATURE AT PREVIOUS TIME LEVEL
!     *PSST*         SST

!     OUTPUT PARAMETERS (REAL):

!     *PTSK*         NEW SKIN TEMPERATURE 

!     METHOD
!     ------

!     The cool skin formulation follows Fairall et al. (1996) and  
!     depends ON surface energy balance and wind speed. The warm skin 
!     model uses the skin temperature as a prognostic variable for the 
!     top ocean layer. Two formulations exist:
!     Formulation A follows she diagnostic form by Webster et al. (1996)
!       cast into a empirical prognostic form. 
!     Formulation C is based on derivation by Xubin Zeng

!     For more details see Beljaars (1997): Air-sea interaction in the 
!     ECMWF model, in Seminar on Atmospher-surface interaction, 
!     8-12 September 1997. 

!     Both formulations can be activated independently by 
!     switches. 

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! 
! USE YOS_EXC  , ONLY : RKAP   ,LEOCWA   ,LEOCCO
! USE YOS_CST  , ONLY : RG     ,RCPD     ,RETV    ,RLVTT

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&           !yomcst  (& yos_exc)
      & RKAP     ,RG       ,RETV     ,RLVTT    ,&           !yoevdf  (& yos_exc)
      & RCPD                                                !yomcst  (& yos_cst)

!! There parameters turn on warm layer (WA) and cold skin (CO).
!! There could be added as namelist parameters but here are constants
!!      & LEOCWA   ,LEOCCO                                    !yoephy  (& yos_exc)

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSTR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSTR(:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU10(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV10(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSK(:) 

INTEGER(KIND=JPIM) :: JL
CHARACTER*1 CHVER

REAL(KIND=JPRB) :: ZNUW,ZROW,ZROA,ZCPW,ZKW,ZG,ZROC,ZCONM13,ZCON23,ZCON34,&
 & ZCON2,ZCON3,ZCON4,ZCON5,ZQ,ZQ2,ZLAMB,ZDELTA,ZFC,&
 & ZSRD,ZGU,ZDSST,ZZ,ZPARZI,ZEPDU2,&
 & ZEPUST,ZUST2,ZDZC,ZFI,ZDU2,ZDL,&
 & ZENHAN,ZA1,ZA2,ZD1,ZD2,ZROWT,ZROWQ,ZUSTW2,&
 & ZWST2,ZDZ,ZPHI,ZROADRW,ZAN

REAL(KIND=JPRB) :: ZBUO(KLON),ZU(KLON),ZALPHA(KLON),ZDCOOL(KLON),&
                 & ZDWARM(KLON),ZUST(KLON)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LEOCWA, LEOCCO

!     1. Initialize constants for ocean warm layer and cool skin

!     1.1 General

IF (LHOOK) CALL DR_HOOK('VOSKIN_MOD:VOSKIN',0,ZHOOK_HANDLE)

LEOCCO=.TRUE.      !    turn on cold skin  - if vokin ON
LEOCWA=.TRUE.      !    turn on warm layer - if vokin ON

CHVER="C"          !    formulation A,B, or C
ZNUW=1.E-6_JPRB    !    kinematic viscosity of water        (m2/s)
ZROW=1025._JPRB    !    density of water                    (kg/m3)
ZROA=1.2_JPRB      !    density of air (approximately)      (kg/m3)
ZCPW=4190._JPRB    !    heat capacity of water              (J/kgK)
ZKW=0.6_JPRB       !    thermal conductivity of water       (W/mK)
ZG=RG              !    gravitational constant              (m/s2)
ZPARZI=1000._JPRB  !    BL height for convective scaling    (m)
ZEPDU2=0.01_JPRB   !    security constant for velocity**2   (m2/s2)
ZEPUST=0.0001_JPRB !    security constant for velocity      (m/s)
ZROADRW=ZROA/ZROW  !    Density ratio                      (-)

ZROC=ZROW*ZCPW
ZCONM13=-1._JPRB/3._JPRB
ZCON23=2._JPRB/3._JPRB
ZCON34=3._JPRB/4._JPRB

!     1.2A Warm layer parametrization constants

ZDZ=2._JPRB        !    depth scale                         (m)
ZENHAN=1.4_JPRB    !    enhancement factor                  (-)

ZA1=0.002_JPRB     !    wind function constant (U<2)
ZD1=-0.000185_JPRB !    wind function constant (U<2)

ZA2=0.00265_JPRB   !    wind function constant (U>2)
ZD2=-0.00105_JPRB  !    wind function constant (U>2)

!     1.2C Warm layer parametrization constants
!
ZDZC=3._JPRB       !    depth scale                         (m)
                   !    ZFI=Fraction of solar radiation absorbed in warm layer (-)
ZFI=1._JPRB-0.28_JPRB*exp(-71.5_JPRB*ZDZC)-0.27_JPRB*exp(-2.8_JPRB*ZDZC) &
  &  - 0.45_JPRB*EXP(-0.07_JPRB*ZDZC)
! ZFI=1._JPRB
ZCON3=ZDZC*RKAP*RG/(ZROA/ZROW)**1.5_JPRB
ZAN=0.3_JPRB       !    Nu (exponent of temperature profile)
ZCON4=(ZAN+1.0_JPRB)*RKAP*SQRT(ZROA/ZROW)/ZDZC
ZCON5=(ZAN+1.0_JPRB)/(ZAN*ZDZC)
!     1.3 Cool skin parametrization constants

ZCON2=16._JPRB*ZG*ZROC*ZNUW**3/(ZKW**2)

!     2. General 

IF (LEOCWA .OR. LEOCCO) THEN
  ZDCOOL(KIDIA:KFDIA)=0.0_JPRB
  ZDWARM(KIDIA:KFDIA)=0.0_JPRB
  
  DO JL=KIDIA,KFDIA

!     Atmospheric buoyancy and wind
    ZROWQ=PAHFL(JL)/RLVTT
    ZROWT=PAHFS(JL)/RCPD
    ZBUO(JL)=RG*(-RETV*ZROWQ-ZROWT/PTSKM1M(JL))/ZROA
    IF (ZBUO(JL) < 0.0_JPRB) THEN
      ZWST2=0.0_JPRB
    ELSE
      ZWST2=(ZBUO(JL)*ZPARZI)**ZCON23
    ENDIF
    ZDU2=MAX(ZEPDU2,PU10(JL)**2+PV10(JL)**2)
    ZU(JL)=SQRT(ZDU2+ZWST2)

    ZUST2=SQRT(PUSTR(JL)**2+PVSTR(JL)**2)/ZROA
    ZUST2=ZUST2*((ZDU2+ZWST2)/ZDU2)
    ZUST(JL)=MAX(SQRT(ZUST2),ZEPUST)

!     Ocean buoyancy 
    ZALPHA(JL)=MAX(1.E-5_JPRB,1.E-5_JPRB*(PSST(JL)-273._JPRB))
  ENDDO
ENDIF

!     3. Cool skin (Fairall et al. 1996)

IF (LEOCCO) then
  DO JL=KIDIA,KFDIA

!      3.2 Apply empirical formulas

    ZUSTW2=ZROADRW*ZUST(JL)**2
    ZQ=MAX(1.0_JPRB,-PSLRFL(JL)-PAHFS(JL)-PAHFL(JL))
    ZLAMB=6._JPRB*(1.0_JPRB+(ZQ*ZALPHA(JL)*ZCON2/ZUSTW2**2)**ZCON34)**ZCONM13

    ZDELTA=ZLAMB*ZNUW/SQRT(ZUSTW2)

!          Solar absorption

    ZFC=0.065_JPRB+11._JPRB*ZDELTA&
     & -(6.6E-5_JPRB/ZDELTA)*(1.0_JPRB-EXP(-ZDELTA/8.E-4_JPRB))  
    ZFC=MAX(ZFC,0.01_JPRB)
    ZQ2=MAX(1.0_JPRB,-ZFC*PSSRFL(JL)+ZQ)
    ZDCOOL(JL)=-ZDELTA*ZQ2/ZKW
  ENDDO
ENDIF

IF (LEOCWA) then
  IF (CHVER == "A") THEN 

!     2.2 Warm layer; formulation A (empirical adapted from Webster al. 1996)

    DO JL=KIDIA,KFDIA
      ZDSST=MAX(PTSKM1M(JL)-PSST(JL)-ZDCOOL(JL),0.0_JPRB)
      IF (ZU(JL) < 2.0_JPRB) THEN
        ZGU=ZENHAN*(ZA1+ZD1*LOG(ZU(JL)))
      ELSE
        ZGU=ZENHAN*(ZA2+ZD2*LOG(ZU(JL)))
      ENDIF
      ZSRD=PSSRFL(JL)/0.93_JPRB
      ZZ=1.0_JPRB+PTMST/(ZGU*ZROC*ZDZ)
      ZDWARM(JL)=MAX(0.0_JPRB,(ZDSST+ZSRD*PTMST/(ZDZ*ZROC))/ZZ)
    ENDDO
  ELSEIF (CHVER == "C") THEN 
!
!     2.2 Warm layer; formulation C (Xubin Zeng)
!
    DO JL=KIDIA,KFDIA
        ZDSST=PTSKM1M(JL)-PSST(JL)-ZDCOOL(JL)

        ZSRD=(PSSRFL(JL)*ZFI+PSLRFL(JL)+PAHFS(JL)+PAHFL(JL))/ZROC

         IF (ZDSST > 0.0_JPRB .AND. ZSRD < 0.0_JPRB) THEN 
           ZDL=ZUST(JL)**2*(ZROA/ZROW)&
               &   *SQRT(ZDSST/(5._JPRB*ZDZC*RG*ZALPHA(JL)/ZAN))       
         ELSE
           ZDL=ZSRD
         ENDIF
        ZDL=ZCON3*ZALPHA(JL)*ZDL/ZUST(JL)**3

        IF (ZDL > 0.0_JPRB) THEN 
          ZPHI=1._JPRB+5._JPRB*ZDL
        ELSE
          ZPHI=1._JPRB/SQRT(1._JPRB-16._JPRB*ZDL)
        ENDIF 
        
        ZZ=1.0_JPRB+ZCON4*PTMST*ZUST(JL)/ZPHI
        ZDWARM(JL)=MAX(0.0_JPRB,(ZDSST+ZCON5*ZSRD*PTMST)/ZZ)
        
    ENDDO
  ENDIF
ENDIF

!     3. Apply warm layer and cool skin effects

PTSK(KIDIA:KFDIA)=PSST(KIDIA:KFDIA)+ZDWARM(KIDIA:KFDIA)+ZDCOOL(KIDIA:KFDIA)

IF (LHOOK) CALL DR_HOOK('VOSKIN_MOD:VOSKIN',1,ZHOOK_HANDLE)
END SUBROUTINE VOSKIN


END MODULE MO_VOSKIN
