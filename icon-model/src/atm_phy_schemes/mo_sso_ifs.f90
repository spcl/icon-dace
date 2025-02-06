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
! Source module for the IFS version of the Lott and Miller (1997) SSO scheme
!
!-----------------------------------------------------------------------------

MODULE mo_sso_ifs

  USE mo_kind   ,ONLY: JPRB=>wp     , &
    &                  JPIM=>i4     , &
    &                  vp

  USE mo_physical_constants , ONLY :   &
    & grav

  USE mo_cuparameters , ONLY :                                     &
    & rg       ,rd      ,rcpd                                     ,&
    & lhook,   dr_hook  ,lphylin, rlpdrag                         ,&
    & GRFPLM, GTENLIM, GSSEC, GTSEC, GVSEC

  USE mo_nwp_parameters,  ONLY: t_phy_params


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gwdrag

! Tunable parameters
! ------------------
REAL (KIND = JPRB) ::       &
  Gkdrag                  , &   ! gw drag constant (set in mo_nwp_tuning_nml)
  Gkwake                  , &   ! gw drag constant (set in mo_nwp_tuning_nml)
  Grcrit                  , &   ! critical Richardson number (set in mo_nwp_tuning_nml)
  Gfrcrit                       ! critical Froude number (determines depth of blocking layer; set in mo_nwp_tuning_nml)


CONTAINS

SUBROUTINE GWDRAG(KIDIA,KFDIA,KLON,KLEV,PARAMS,&
 & PAPHM1,PAPM1,PGEOM1,&
 & PTM1,PUM1,PVM1,&
 & PHSTD,PGAMMA,PTHETA,PSIG,&
 & PDT,IKENVH,&
 ! TENDENCY COEFFICIENTS (OUTPUT)
 & PSOTEU,PSOTEV,&
 & PUSTR_SSO,PVSTR_SSO) 
                               

!**** *GWDRAG* - DOES THE GRAVITY WAVE PARAMETRIZATION.

!     PURPOSE.
!     --------

!          COMPUTES THE COEFFICIENTS FOR COMPUTATION OF THE 
!          PHYSICAL TENDENCIES OF THE PROGNOSTIC VARIABLES U,V  AND T DUE TO  
!          VERTICAL TRANSPORTS BY SUBGRIDSCALE OROGRAPHICALLY EXCITED GRAVITY WAVES
!          THE TENDENCIES AND OTHER DIAGNOSTIC QUANTITIES ARE THEN COMPUTER 
!          IN CALLPAR.  

!          NOTE THAT THE VERTICAL INDEXING OF PAPHM1 IS FROM 0:KLEV TO MAKE IT 
!          COMPATIBLE WITH CALLPAR

!**   INTERFACE.
!     ----------
!          CALLED FROM *CALLPAR*.

!          THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
!          U,V,T AND P AT T-1.

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===

!        *KIDIA*       START POINT
!        *KFDIA*       END POINT
!        *KLON*        NUMBER OF GRID POINTS PER PACKET
!        *KLEV*        NUMBER OF VERTICAL LEVELS
!        *KGL1*        ROW NUMBER FOR PRINTINGS  (UNSUED)
!        *PHSTD*       STANDARD DEVIATION OF THE SUBGRID SCALE OROGRAPHY
!        *PTHETA*      GEOGRAPHICAL ORIENTATION 
!                      OF THE SUBGRID SCALE OROGRAPHY
!        *PGAMMA*      ANISOTROPY OF THE SUBGRID SCALE OROGRAPHY
!        *PSIG*        MEAN SLOPE OF THE SUBGRID SCALE OROGRAPHY
!        *PAPHM1*      PRESSURE AT HALF LEVELS AT T-1
!        *PAPM1*       PRESSURE AT FULL LEVELS AT T-1
!        *PGEOM1*      GEOPOTENTIAL AT T-1
!        *PTM1*        TEMPERATURE AT T-1
!        *PUM1*        X-VELOCITY COMPONENT AT T-1
!        *PVM1*        Y-VELOCITY COMPONENT AT T-1

!     ==== OUTPUTS ===

!        *PSOTEU*     EXPLICIT TENDENCY COEFFICIENT OF U-COMP OF WIND
!        *PSOTEV*     EXPLICIT TENDENCY COEFFICIENT OF V-COMP OF WIND
!        *PSOBETA*    IMPLICIT TENDENCY COEFFICIENT OF U,V-COMP OF WIND
!        *LLGWD*      TRUE IF PHSTD > 50 (IE GWDRAG SCHEME ACTIVE)

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------
!         THE SCHEME CONSISTS OF THREE PARTS,THE CALCULATION OF LOWLEVEL
!         DRAG DUE TO WAKE EFFECTS OF THE SUBGRID OROGRAPHY,THE GRAVITY 
!         WAVE STRESS AND THE STRESS PROFILE IN THE VERTICAL.
!         THE STRESS IS COMPUTED USING A LOW-LEVEL WIND,STATIC STABILITY
!         AND  SUBGRID OROGRAPHIC PARAMETERS.THESE ARE THE STANDARD 
!         DEVIATION OF THE HEIGHT,THE ANGLE RELATIVE TO EAST, THE
!         ANISOTROPY,AND THE SLOPE.
!         A WAV RICHARDSON NUMBER IS COMPUTED AT EACH LEVEL AND BY
!         REQUIRING THAT ITS VALUE IS NEVER LESS THAN A CRITICAL ONE
!         A VALUE OF STRESS IS DETERMINED AT EACH MODEL LEVEL.
!         THE CRITICAL FROUDE NUMBER OF THE FLOW DETERMINES THE DEPTH OF
!         THE LOWLEVEL DRAG.  

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!         SEE MODEL DOCUMENTATION

!        SEE ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S."

!     AUTHOR.
!     -------
!      M.MILLER + B.RITTER   E.C.M.W.F.     15/06/86.

!     MODIFICATIONS.
!     --------------
!      D.Salmond  15-Oct-01 FULLIMP mods
!      M.Hamrud   01-Oct-2003 CY28 Cleaning
!      A.Brown    17-Aug-2004 Introduction of cutoff moutain
!      A. Orr     01-Nov-2005 Subgrid oro and vdf integrated together
!      P.Bechtold 14-Jan-2009 move IKTOP computation to setup routine
!-----------------------------------------------------------------------

!x USE PARKIND1  ,ONLY : JPIM     ,JPRB
!x USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!x 
!x USE YOMCST   , ONLY : RG
!x USE YOEGWD   , ONLY : YREGWD
!x USE YOEPHLI  , ONLY : YREPHLI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
TYPE(t_phy_params),INTENT(in)    :: params
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) !note callpar indexing convention 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHSTD(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAMMA(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHETA(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIG(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
INTEGER(KIND=JPIM),INTENT(OUT)   :: IKENVH(KLON)
REAL(KIND=vp)     ,INTENT(OUT)   :: PSOTEU(KLON,KLEV)
REAL(KIND=vp)     ,INTENT(OUT)   :: PSOTEV(KLON,KLEV)
!REAL(KIND=JPRB)   ,INTENT(OUT), OPTIONAL   :: PSOBETA(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT), OPTIONAL   :: PUSTR_SSO(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT), OPTIONAL   :: PVSTR_SSO(KLON)

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) ::  &
 & ICRIT(KLON),&
 & IKCRITH(KLON),&
!& IKENVH(KLON),&
 & IKNU(KLON),&
 & IKNU2(KLON)  

REAL(KIND=JPRB) ::   ZTAU(KLON,KLEV+1),&
 & ZSTAB(KLON,KLEV+1),&
 & ZVPH(KLON,KLEV+1),&
 & ZRHO(KLON,KLEV+1),&
 & ZRI(KLON,KLEV+1),&
 & ZPSI(KLON,KLEV+1),&
 & ZZDEP(KLON,KLEV)  
REAL(KIND=JPRB) ::  ZULOW(KLON),&
 & ZVLOW(KLON),&
 & ZNU(KLON),&
 & ZD1(KLON),&
 & ZD2(KLON),&
 & ZDMOD(KLON)  

INTEGER(KIND=JPIM) :: ILEVP1, ILONS, JK, JL, IKTOP

REAL(KIND=JPRB) :: ZABSV, ZB, ZBLOCK, ZC, ZCONB,&
 & ZCS, ZDELP, ZDELP0, ZEFF, ZKDRAG, &
 & ZRATIO, ZSS, ZSTD, ZTEMP, ZZD1
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLGWD(KLON)

REAL(KIND=JPRB) :: zstrdu (klon,klev+1)    ! flux of u-momentum (GWD and Blocking)
REAL(KIND=JPRB) :: zstrdv (klon,klev+1)    ! flux of v-momentum (GWD and Blocking)
REAL(KIND=JPRB) :: zgdph

REAL(KIND=JPRB)   :: PSOBETA(KLON,KLEV)

!x #include "gwprofil.intfb.h"
!x #include "gwsetup.intfb.h"

!*************************************************************

!     ------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

IF (LHOOK) CALL DR_HOOK('GWDRAG',0,ZHOOK_HANDLE)

ILONS=KFDIA-KIDIA+1
!     ------------------------------------------------------------------

!       SELECT POINTS FOR SCHEME TO OPERATE(SGS OROG >50METRES)  ****

DO JL=KIDIA,KFDIA
  IF (PHSTD(JL) > 50._JPRB) THEN
    LLGWD(JL)=.TRUE.
  ELSE
    LLGWD(JL)=.FALSE.
  ENDIF
ENDDO

! Set tuning parameters
  Gkdrag  = params%Gkdrag
  Gkwake  = params%Gkwake
  Grcrit  = params%Grcrit
  Gfrcrit = params%Gfrcrit

!       Constant for surface stress

IF (LPHYLIN) THEN
  ZKDRAG=RLPDRAG
ELSE
  ZKDRAG=GKDRAG
ENDIF

!     ------------------------------------------------------------------

!*         2.     PRECOMPUTE BASIC STATE VARIABLES.
!*                ---------- ----- ----- ----------

CALL GWSETUP &
 & ( KIDIA , KFDIA , KLON  , KLEV , LLGWD,&
 & IKCRITH, ICRIT,&
 & IKENVH,IKNU,IKNU2,&
 & params%NKTOPG,&
 & PAPHM1, PAPM1 , PUM1   , PVM1 , PTM1 , PGEOM1, PHSTD,&
 & ZRHO  , ZRI   , ZSTAB  , ZTAU , ZVPH , ZPSI , ZZDEP,&
 & ZULOW , ZVLOW,&
 & PTHETA , PGAMMA,ZNU  ,ZD1  , ZD2, ZDMOD )  

!***********************************************************

!*         3.      COMPUTE THE SURFACE GRAVITY WAVE STRESS
!                  ---------------------------------------

ZBLOCK=0.0_JPRB

DO JL=KIDIA,KFDIA
  IF(LLGWD(JL)) THEN
 
    ZBLOCK=(PGEOM1(JL,IKENVH(JL))/RG)    
    ZSTD=PHSTD(JL)
    ZEFF=MAX(0.0_JPRB,3._JPRB*ZSTD-ZBLOCK)*2.0_JPRB
    ZTAU(JL,KLEV+1)=ZKDRAG*ZRHO(JL,KLEV+1)*PSIG(JL)*ZEFF**2 &
     &/9._JPRB/ZSTD*ZVPH(JL,KLEV+1)*ZDMOD(JL)*SQRT(ZSTAB(JL,KLEV+1))

  ELSE
    ZTAU(JL,KLEV+1)=0.0_JPRB
  ENDIF

ENDDO

!*         4.      COMPUTE GRAVITY WAVE STRESS PROFILE.
!*                 -----------------------------------

CALL GWPROFIL &
 & ( KIDIA ,  KFDIA , KLON  , KLEV,&
 & LLGWD,&
 & IKCRITH, ICRIT , IKENVH, IKNU,&
 & IKNU2  , PAPHM1, ZRHO  , ZSTAB , ZVPH,&
 & ZRI    , ZTAU  , ZDMOD , PSIG  , PHSTD)  


!*         5.      COMPUTE TENDENCY COEFFICIENTS.
!*                 -----------------------------

DO JK=1,KLEV
 DO JL=KIDIA,KFDIA
  PSOTEU(JL,JK)=0.0_JPRB
  PSOTEV(JL,JK)=0.0_JPRB
  PSOBETA(JL,JK)=0.0_JPRB
 ENDDO
ENDDO 


ILEVP1=KLEV+1

IKTOP=params%NGWDTOP

DO JK=1,KLEV
!OCL  VCT(CEX)
  DO JL=KIDIA,KFDIA
    IF(LLGWD(JL)) THEN
!      COMPUTE EXPLICIT COEFFICENTS FOR GRAVITY WAVE PORTION 
!      INITIALLY SET IMPLICIT BLOCKING COEFFICENTS TO ZERO 

      ZDELP =PAPHM1(JL,JK)-PAPHM1(JL,JK-1)
! adjust in pressure < 0.1hPa
      IF( PAPM1(JL,JK) < GRFPLM) THEN
        ZDELP0=PAPHM1(JL,IKTOP)-PAPHM1 (JL,0)
        ZTEMP=-RG*(ZTAU(JL,IKTOP+1)-ZTAU(JL,1))*ZDELP/ &
         &        (ZVPH(JL,ILEVP1)*ZDELP0**2)
      ELSE
        ZTEMP =-RG*(ZTAU(JL,JK+1)-ZTAU(JL,JK))/(ZDELP*ZVPH(JL,ILEVP1))
      ENDIF
      PSOTEU(JL,JK)=(ZULOW(JL)*ZD1(JL)-ZVLOW(JL)*ZD2(JL))*ZTEMP/ZDMOD(JL)
      PSOTEV(JL,JK)=(ZVLOW(JL)*ZD1(JL)+ZULOW(JL)*ZD2(JL))*ZTEMP/ZDMOD(JL)
      PSOBETA(JL,JK)=0.0_JPRB
      
!      LIMIT THE WIND TENDENCIES IN THE STRATOSPHERE (probably not necessary now ?)
       IF(JK < params%NGWDLIM) THEN
        PSOTEU(JL,JK)=SIGN(MIN(ABS(PSOTEU(JL,JK)),REAL(GTENLIM,vp)),PSOTEU(JL,JK))
        PSOTEV(JL,JK)=SIGN(MIN(ABS(PSOTEV(JL,JK)),REAL(GTENLIM,vp)),PSOTEV(JL,JK))
       ENDIF

!     FOR BLOCKING LEVELS RECALCULATE THE IMPLICIT COEFFICIENTS
!     AND SET EXPLICIT GRAVITY WAVE COEFFICIENTS TO ZERO      
      IF(JK >= IKENVH(JL)) THEN
        ZB    =1.0_JPRB-0.18_JPRB*PGAMMA(JL)-0.04_JPRB*PGAMMA(JL)**2
        ZC    =0.48_JPRB*PGAMMA(JL)+0.3_JPRB*PGAMMA(JL)**2
        ZCS   =COS(ZPSI(JL,JK))**2
        ZSS   =1.0_JPRB-ZCS
        ZCONB =GKWAKE*PSIG(JL)/(2.0_JPRB*PHSTD(JL)) !  ZCONB INDEPENDENT OF ZTMST
        ZABSV =0.5_JPRB*SQRT(PUM1(JL,JK)**2+PVM1(JL,JK)**2)
        ZZD1  =ZB*ZCS+ZC*ZSS
        ZRATIO=(ZCS+PGAMMA(JL)*ZSS)/(PGAMMA(JL)*ZCS+ZSS)      
!xmk    PSOBETA(JL,JK)=MAX(0.0_JPRB,2.0_JPRB-1.0_JPRB/ZRATIO)*ZCONB*ZZDEP(JL,JK)*ZZD1*ZABSV 
        PSOBETA(JL,JK)=MAX(0.0_JPRB,2.0_JPRB-1.0_JPRB/ZRATIO)*ZCONB*ZZDEP(JL,JK)*ZZD1*ZABSV * PDT
!xmk    PSOTEU(JL,JK)=0.0_JPRB
!xxx    PSOTEV(JL,JK)=0.0_JPRB
!amk
!        Partially implicit tendency calculation (taken from mo_sso_cosmo.f90)
!        ---------------------------------------------------------------------
        PSOTEU(JL,JK)=-pum1(jl,jk)/pdt * psobeta(jl,jk)/(1._jprb+(psobeta(jl,jk)))
        PSOTEV(JL,JK)=-pvm1(jl,jk)/pdt * psobeta(jl,jk)/(1._jprb+(psobeta(jl,jk)))
!xxx
      ENDIF

    ENDIF
  ENDDO
ENDDO


!     Surface stress calculation for output 
!     (taken from mo_sso_cosmo.f90)

!     Initialize flux at top
!     ----------------------
      DO jl=kidia,kfdia
        zstrdu(jl,1)=0.0_JPRB
        zstrdv(jl,1)=0.0_JPRB
      END DO

!     Increment flux based on tendency in each layer
!     ----------------------------------------------
      DO jk=1,klev
        DO jl=kidia,kfdia
          zgdph=-grav  /(paphm1(jl,jk+1)-paphm1(jl,jk))
          zstrdu(jl,jk+1)=psoteu(jl,jk)/zgdph + zstrdu(jl,jk)
          zstrdv(jl,jk+1)=psotev(jl,jk)/zgdph + zstrdv(jl,jk)
        END DO
      END DO

!     Store flux at surface
!     ---------------------
      IF(PRESENT(pustr_sso))THEN
        DO jl=kidia,kfdia
          pustr_sso(jl)=zstrdu(jl,klev+1)
          pvstr_sso(jl)=zstrdv(jl,klev+1)
        END DO
      ENDIF


IF (LHOOK) CALL DR_HOOK('GWDRAG',1,ZHOOK_HANDLE)
END SUBROUTINE GWDRAG

SUBROUTINE GWSETUP &
 & ( KIDIA , KFDIA , KLON , KLEV , LDGWD,&
 & KKCRITH, KCRIT,&
 & KKENVH, KKNU  , KKNU2,&
 & NKTOPG,&
 & PAPHM1, PAPM1 , PUM1   , PVM1 , PTM1  , PGEOM1, PHSTD,&
 & PRHO  , PRI   , PSTAB  , PTAU , PVPH  , PPSI,   PZDEP,&
 & PULOW , PVLOW,&
 & PTHETA, PGAMMA,  PNU  ,  PD1  ,  PD2  , PDMOD         )  

!**** *GWSETUP*

!     PURPOSE.

!**   INTERFACE.
!     ----------
!          FROM *GWDRAG*

!        EXPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===

!        *KIDIA*       START POINT
!        *KFDIA*       END POINT
!        *KLON*        NUMBER OF GRID POINTS PER PACKET
!        *KLEV*        NUMBER OF VERTICAL LEVELS
!        *LDGWD*       TRUE IF SSO ACTIVE
!        *PAPHM1*      PRESSURE AT HALF LEVELS AT T-1
!        *PAPM1*       PRESSURE AT FULL LEVELS AT T-1
!        *PUM1*        X-VELOCITY COMPONENT AT T-1
!        *PVM1*        Y-VELOCITY COMPONENT AT T-1
!        *PTM1*        TEMPERATURE AT T-1
!        *PGEOM1*      GEOPOTENTIAL AT T-1
!        *PHSTD*       STANDARD DEVIATION OF THE OROGRAPHY
!        *PTHETA*      GEOGRAPHICAL ORIENTATION OF THE OROGRAPHY
!        *PGAMMA*      ANISOTROPY OF THE OROGRAPHY

!     ==== OUTPUTS ===

!        *KKCRITH*     MAXIMUM LEVEL FOR WAVE BREAKING
!        *KCRIT*       CRITICAL LEVEL
!        *KKENVH*      BLOCKING LEVEL (TOP OF ENVELOPE LAYER)
!        *KKNU*        MODEL LEVEL AT 4*PHSTD HEIGHT 
!        *KKNU2*       MODEL LEVEL AT 3*PHSTD HEIGHT
!        *PRHO*        AIR DENSITY
!        *PRI*         MEAN FLOW RICHARDSON NUMBER
!        *PSTAB*       SQUARE OF THE BRUNT-VAISALA FREQUENCY
!        *PTAU*        STRESS PRODUCED BY GRAVITY WAVES
!        *PVPH*        WIND PROJECTION ALONG SURFACE STRESS AT HALF LEVELS
!        *PPSI*        ANGLE BETWEEN THE LOW-LEVEL WIND DIRECTION
!                      AND THE PRINCIPAL AXIS OF TOPOGRAPHY
!        *PZDEP*       VERTICAL LEAKINESS (LENGTH OF THE MOUNTAIN SEEN
!                      BY THE FLOW)
!        *PULOW*       LOW-LEVEL BLOCKED FLOW X-COMPONENT 
!        *PVLOW*       LOW-LEVEL BLOCKED FLOW Y-COMPONENT
!        *PNU*         NON-DIMENSIONAL HEIGHT 
!        *PD1*         X-DIRECTION OF THE SURFACE STRESS
!        *PD2*         Y-DIRECTION OF THE SURFACE STRESS
!        *PDMOD*       MODULUS OF THE SURFACE STRESS (PART)

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        SEE ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S."

!     AUTHOR.
!     -------
!      M.MILLER

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      A. Beljaars   15-Jun-2012 suppressing some resolution dependence
!-----------------------------------------------------------------------

!x USE PARKIND1  ,ONLY : JPIM     ,JPRB
!x USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!x 
!x USE YOEGWD   , ONLY : YREGWD
!x USE YOMCST   , ONLY : RG       ,RD       ,RCPD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
LOGICAL           ,INTENT(IN)    :: LDGWD(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KKCRITH(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCRIT(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KKENVH(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KKNU(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KKNU2(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: NKTOPG
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) !note callpar indexing convention  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHSTD(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRHO(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRI(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTAB(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPSI(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZDEP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PULOW(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVLOW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHETA(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAMMA(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD1(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMOD(KLON) 

!-----------------------------------------------------------------------

LOGICAL :: LL1(KLON,KLEV+1)
INTEGER(KIND=JPIM) :: IKNUB(KLON),IKNUL(KLON)

REAL(KIND=JPRB) :: ZHCRIT(KLON,KLEV),ZVPF(KLON,KLEV),ZDP(KLON,KLEV),ZSQST(KLON,KLEV)
REAL(KIND=JPRB) :: ZNORM(KLON),ZB(KLON),ZC(KLON),&
 & ZULOW(KLON),ZVLOW(KLON),ZNUP(KLON),ZNUM(KLON)  

INTEGER(KIND=JPIM) :: ILEVH, JK, JL

LOGICAL :: LLO

REAL(KIND=JPRB) :: ZCONS1, ZCONS2, ZDELP, ZDWIND, ZGGEENV, ZGGEOM1,&
 & ZGVAR, ZHGEO, ZPHI, ZRHOM, ZRHOP, ZST, ZSTABM, &
 & ZSTABP, ZU, ZVT1, ZVT2, ZWIND, Z1D2RG
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!     ------------------------------------------------------------------

!*         1.1   COMPUTATIONAL CONSTANTS
!                -----------------------

IF (LHOOK) CALL DR_HOOK('GWSETUP',0,ZHOOK_HANDLE)

ILEVH =KLEV/3

ZCONS1=1.0_JPRB/RD
ZCONS2=RG**2/RCPD
Z1D2RG=1.0_JPRB/(2.0_JPRB*RG)

!     ------------------------------------------------------------------

!*         2.
!                --------------

!     ------------------------------------------------------------------

!*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
!*                 LOW LEVEL WIND AND SET INDICATOR FOR CRITICAL LEVELS.

DO JL=KIDIA,KFDIA
  KKNU(JL)      =KLEV
  KKNU2(JL)     =KLEV
  IKNUB(JL)     =KLEV
  IKNUL(JL)     =KLEV
  PGAMMA(JL)    =MAX(PGAMMA(JL),GTSEC)
  LL1(JL,KLEV+1)=.FALSE.
ENDDO

!*      DEFINE VARIOUS LEVELS IN THE FLOW
!       ----------------------------
DO JK=KLEV,ILEVH,-1
  DO JL=KIDIA,KFDIA
    ZHCRIT(JL,JK)=4.4_JPRB*PHSTD(JL)
    ZHGEO=(PGEOM1(JL,JK)+PGEOM1(JL,JK-1))*Z1D2RG
    LL1(JL,JK)=(ZHGEO > ZHCRIT(JL,JK))
    IF(LL1(JL,JK).NEQV.LL1(JL,JK+1)) THEN
      KKNU(JL)=JK
    ENDIF
  ENDDO
ENDDO
DO JK=KLEV,ILEVH,-1
  DO JL=KIDIA,KFDIA
    ZHCRIT(JL,JK)=3.3_JPRB*PHSTD(JL)
    ZHGEO=(PGEOM1(JL,JK)+PGEOM1(JL,JK-1))*Z1D2RG
    LL1(JL,JK)=(ZHGEO > ZHCRIT(JL,JK))
    IF(LL1(JL,JK).NEQV.LL1(JL,JK+1)) THEN
      KKNU2(JL)=JK
    ENDIF
  ENDDO
ENDDO
DO JK=KLEV,ILEVH,-1
  DO JL=KIDIA,KFDIA
    ZHCRIT(JL,JK)=2.2_JPRB*PHSTD(JL)
    ZHGEO=(PGEOM1(JL,JK)+PGEOM1(JL,JK-1))*Z1D2RG
    LL1(JL,JK)=(ZHGEO > ZHCRIT(JL,JK))
    IF(LL1(JL,JK).NEQV.LL1(JL,JK+1)) THEN
      IKNUB(JL)=JK
    ENDIF
  ENDDO
ENDDO
DO JK=KLEV,ILEVH,-1
  DO JL=KIDIA,KFDIA
    ZHCRIT(JL,JK)=1.1_JPRB*PHSTD(JL)
    ZHGEO=(PGEOM1(JL,JK)+PGEOM1(JL,JK-1))*Z1D2RG
    LL1(JL,JK)=(ZHGEO > ZHCRIT(JL,JK))
    IF(LL1(JL,JK).NEQV.LL1(JL,JK+1)) THEN
!            IKNUL(JL)=KLEV   ! BUG
      IKNUL(JL)=JK
    ENDIF
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  KKNU(JL) =MIN(KKNU(JL),NKTOPG)
  IKNUB(JL)=MIN(IKNUB(JL),NKTOPG)
  IF(IKNUB(JL) == NKTOPG) IKNUL(JL)=KLEV
  IF(IKNUB(JL) == IKNUL(JL)) IKNUL(JL)=IKNUL(JL)+1
ENDDO

!C*     INITIALIZE VARIOUS ARRAYS

DO JL=KIDIA,KFDIA
  PTAU(JL,KLEV+1)  =0.0_JPRB
  PRHO(JL,KLEV+1)  =0.0_JPRB
  PSTAB(JL,KLEV+1) =0.0_JPRB
  PSTAB(JL,1)      =0.0_JPRB
  PRI(JL,KLEV+1)   =9999.0_JPRB
  PPSI(JL,KLEV+1)  =0.0_JPRB
  PRI(JL,1)        =0.0_JPRB
  PVPH(JL,1)       =0.0_JPRB
  PULOW(JL)        =0.0_JPRB
  PVLOW(JL)        =0.0_JPRB
  ZULOW(JL)        =0.0_JPRB
  ZVLOW(JL)        =0.0_JPRB
  KKCRITH(JL)      =KLEV
  KKENVH(JL)       =KLEV
  KCRIT(JL)        =1
  PNU (JL)         =0.0_JPRB
  ZNUM(JL)         =0.0_JPRB
  LL1(JL,KLEV+1)   =.FALSE.
ENDDO

!*       DEFINE STATIC STABILITY AND AIR DENSITY
!*       ---------------------------------------

DO JK=KLEV,2,-1
  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL)) THEN
      ZDP(JL,JK)=PAPM1(JL,JK)-PAPM1(JL,JK-1)
      PRHO(JL,JK)=2.0_JPRB*PAPHM1(JL,JK-1)*ZCONS1/(PTM1(JL,JK)+PTM1(JL,JK-1))
      PSTAB(JL,JK)=2.0_JPRB*ZCONS2/(PTM1(JL,JK)+PTM1(JL,JK-1))*&
       & (1.0_JPRB-RCPD*PRHO(JL,JK)*(PTM1(JL,JK)-PTM1(JL,JK-1))/ZDP(JL,&
       & JK))  
      PSTAB(JL,JK)=MAX(PSTAB(JL,JK),GSSEC)
      ZSQST(JL,JK)=SQRT(PSTAB(JL,JK))
    ENDIF
  ENDDO
ENDDO

!*         2.2     MEAN BRUNT-VAISALA FREQUENCY 
!                  AND DENSITY FOR WAVE STRESS

DO JK=ILEVH,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL)) THEN
      IF(JK >= (IKNUB(JL)+1).AND.JK <= IKNUL(JL)) THEN
        ZST=ZCONS2/PTM1(JL,JK)*(1.0_JPRB-RCPD*PRHO(JL,JK)*&
         & (PTM1(JL,JK)-PTM1(JL,JK-1))/ZDP(JL,JK))  
        PSTAB(JL,KLEV+1)=PSTAB(JL,KLEV+1)+ZST*ZDP(JL,JK)
        PSTAB(JL,KLEV+1)=MAX(PSTAB(JL,KLEV+1),GSSEC)
        PRHO(JL,KLEV+1)=PRHO(JL,KLEV+1)+PAPHM1(JL,JK-1)*2.0_JPRB*ZDP(JL,&
         & JK)&
         & *ZCONS1/(PTM1(JL,JK)+PTM1(JL,JK-1))  
      ENDIF
    ENDIF
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PSTAB(JL,KLEV+1)=PSTAB(JL,KLEV+1)/(PAPM1(JL,IKNUL(JL))-PAPM1(JL,IKNUB(JL)))
  PRHO(JL,KLEV+1)=PRHO(JL,KLEV+1)/(PAPM1(JL,IKNUL(JL))-PAPM1(JL,IKNUB(JL)))
ENDDO

!********************************************************************

!*     DEFINE LOW LEVEL FLOW FOR WAVE STRESS
!      -------------------------------------

DO JK=KLEV,ILEVH,-1
  DO JL=KIDIA,KFDIA
    IF(JK >= IKNUB(JL).AND.JK <= IKNUL(JL)) THEN
      PULOW(JL)=PULOW(JL)+PUM1(JL,JK)*(PAPHM1(JL,JK)-PAPHM1(JL,JK-1))
      PVLOW(JL)=PVLOW(JL)+PVM1(JL,JK)*(PAPHM1(JL,JK)-PAPHM1(JL,JK-1))
    ENDIF
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  PULOW(JL)=PULOW(JL)/(PAPHM1(JL,IKNUL(JL))-PAPHM1(JL,IKNUB(JL)-1))
  PVLOW(JL)=PVLOW(JL)/(PAPHM1(JL,IKNUL(JL))-PAPHM1(JL,IKNUB(JL)-1))
  ZNORM(JL)=MAX(SQRT(PULOW(JL)**2+PVLOW(JL)**2),GVSEC)
  PVPH(JL,KLEV+1)=ZNORM(JL)
ENDDO

!*******  SETUP OROGRAPHY AXES AND DEFINE PLANE OF PROFILES  *******

DO JL=KIDIA,KFDIA
  LLO=(PULOW(JL) < GVSEC).AND.(PULOW(JL) >= -GVSEC)
  IF(LLO) THEN
    ZU=PULOW(JL)+2.0_JPRB*GVSEC
  ELSE
    ZU=PULOW(JL)
  ENDIF
  ZPHI=ATAN(PVLOW(JL)/ZU)
  PPSI(JL,KLEV+1)=PTHETA(JL)-ZPHI
  ZB(JL)=1.0_JPRB-0.18_JPRB*PGAMMA(JL)-0.04_JPRB*PGAMMA(JL)**2
  ZC(JL)=0.48_JPRB*PGAMMA(JL)+0.3_JPRB*PGAMMA(JL)**2
  PD1(JL)=ZB(JL)-(ZB(JL)-ZC(JL))*(SIN(PPSI(JL,KLEV+1))**2)
  PD2(JL)=(ZB(JL)-ZC(JL))*SIN(PPSI(JL,KLEV+1))*COS(PPSI(JL,KLEV+1))
  PDMOD(JL)=SQRT(PD1(JL)**2+PD2(JL)**2)
ENDDO

!  ************ DEFINE FLOW IN PLANE OF LOWLEVEL STRESS *************

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL))  THEN
      ZVT1       =PULOW(JL)*PUM1(JL,JK)+PVLOW(JL)*PVM1(JL,JK)
      ZVT2       =-PVLOW(JL)*PUM1(JL,JK)+PULOW(JL)*PVM1(JL,JK)
      ZVPF(JL,JK)=(ZVT1*PD1(JL)+ZVT2*PD2(JL))/(ZNORM(JL)*PDMOD(JL))
    ENDIF

    PTAU(JL,JK)  =0.0_JPRB
    PZDEP(JL,JK) =0.0_JPRB
    PPSI(JL,JK)  =0.0_JPRB
    LL1(JL,JK)   =.FALSE.
  ENDDO
ENDDO

DO JK=2,KLEV
!DIR$ novector
  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL)) THEN
      ZDP(JL,JK)=PAPM1(JL,JK)-PAPM1(JL,JK-1)
      PVPH(JL,JK)=((PAPHM1(JL,JK-1)-PAPM1(JL,JK-1))*ZVPF(JL,JK)+&
       & (PAPM1(JL,JK)-PAPHM1(JL,JK-1))*ZVPF(JL,JK-1))&
       & /ZDP(JL,JK)  
      IF(PVPH(JL,JK) < GVSEC) THEN
        PVPH(JL,JK)=GVSEC
!U            KCRIT(JL)=JK     !BUGFIX
        IF (JK  <=  IKNUL(JL)) KCRIT(JL)=JK
      ENDIF
    ENDIF
  ENDDO
ENDDO
!DIR$ vector

!*         2.3     MEAN FLOW RICHARDSON NUMBER.

DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL)) THEN
      ZDWIND=MAX(ABS(ZVPF(JL,JK)-ZVPF(JL,JK-1)),GVSEC)
      PRI(JL,JK)=PSTAB(JL,JK)*(ZDP(JL,JK)/(RG*PRHO(JL,JK)*ZDWIND))**2
      PRI(JL,JK)=MAX(PRI(JL,JK),GRCRIT)
    ENDIF
  ENDDO
ENDDO

!*      DEFINE TOP OF BLOCKED LAYER
!       ----------------------------

DO JK=2,KLEV-1
  DO JL=KIDIA,KFDIA

    IF(LDGWD(JL)) THEN

      IF (JK == KKNU2(JL)) THEN

!       first level: integrate half layer only (below full level)
        ZNUM(JL)=PNU(JL)
        ZWIND=(PULOW(JL)*PUM1(JL,JK)+PVLOW(JL)*PVM1(JL,JK))/&
         & MAX(SQRT(PULOW(JL)**2+PVLOW(JL)**2),GVSEC)  
        ZWIND=MAX(ABS(ZWIND),GVSEC)
        ZDELP=PAPHM1(JL,JK)-PAPM1(JL,JK)
        ZSTABP=ZSQST(JL,JK+1)
        ZRHOP=PRHO(JL,JK+1)
        PNU(JL) = PNU(JL) + (ZDELP/RG)*&
         & (ZSTABP/ZRHOP)/ZWIND  
        IF((ZNUM(JL) <= GFRCRIT).AND.(PNU(JL) > GFRCRIT)&
         & .AND.(KKENVH(JL) == KLEV))&
         & KKENVH(JL)=JK  

      ELSEIF (JK > KKNU2(JL)) THEN

!       first level: integrate full layer (half above, half below full level)
        ZNUM(JL)=PNU(JL)
        ZWIND=(PULOW(JL)*PUM1(JL,JK)+PVLOW(JL)*PVM1(JL,JK))/&
         & MAX(SQRT(PULOW(JL)**2+PVLOW(JL)**2),GVSEC)  
        ZWIND=MAX(ABS(ZWIND),GVSEC)
        ZDELP=PAPHM1(JL,JK)-PAPHM1(JL,JK-1)
        ZSTABM=ZSQST(JL,JK)
        ZSTABP=ZSQST(JL,JK+1)
        ZRHOM=PRHO(JL,JK  )
        ZRHOP=PRHO(JL,JK+1)
        PNU(JL) = PNU(JL) + (ZDELP/RG)*&
         & ((ZSTABP/ZRHOP+ZSTABM/ZRHOM)/2.0_JPRB)/ZWIND  
        IF((ZNUM(JL) <= GFRCRIT).AND.(PNU(JL) > GFRCRIT)&
         & .AND.(KKENVH(JL) == KLEV))&
         & KKENVH(JL)=JK  

      ENDIF

    ENDIF
  ENDDO
ENDDO

!  CALCULATION OF A DYNAMICAL MIXING HEIGHT FOR THE BREAKING
!  OF GRAVITY WAVES:

DO JL=KIDIA,KFDIA
  ZNUP(JL)=0.0_JPRB
  ZNUM(JL)=0.0_JPRB
ENDDO

DO JK=KLEV-1,2,-1
  DO JL=KIDIA,KFDIA

    IF(LDGWD(JL)) THEN

      IF (JK < KKENVH(JL)) THEN

        ZNUM(JL)=ZNUP(JL)
        ZWIND=(PULOW(JL)*PUM1(JL,JK)+PVLOW(JL)*PVM1(JL,JK))/&
         & MAX(SQRT(PULOW(JL)**2+PVLOW(JL)**2),GVSEC)  
        ZWIND=MAX(ABS(ZWIND),GVSEC)
        ZDELP=PAPHM1(JL,JK)-PAPHM1(JL,JK-1)
        ZSTABM=ZSQST(JL,JK)
        ZSTABP=ZSQST(JL,JK+1)
        ZRHOM=PRHO(JL,JK  )
        ZRHOP=PRHO(JL,JK+1)
        ZNUP(JL) = ZNUP(JL) + (ZDELP/RG)*&
         & ((ZSTABP/ZRHOP+ZSTABM/ZRHOM)/2.0_JPRB)/ZWIND  
        IF((ZNUM(JL) <= 1.5_JPRB).AND.(ZNUP(JL) > 1.5_JPRB)&
         & .AND.(KKCRITH(JL) == KLEV))&
         & KKCRITH(JL)=JK  

      ENDIF

    ENDIF

  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  KKCRITH(JL)=MIN(KKCRITH(JL),KKNU(JL))
ENDDO

!     DIRECTIONAL INFO FOR FLOW BLOCKING *************************

DO JK=ILEVH,KLEV
  DO JL=KIDIA,KFDIA
    IF(JK >= KKENVH(JL)) THEN
      LLO=(PUM1(JL,JK) < GVSEC).AND.(PUM1(JL,JK) >= -GVSEC)
      IF(LLO) THEN
        ZU=PUM1(JL,JK)+2.0_JPRB*GVSEC
      ELSE
        ZU=PUM1(JL,JK)
      ENDIF
      ZPHI=ATAN(PVM1(JL,JK)/ZU)
      PPSI(JL,JK)=PTHETA(JL)-ZPHI
    ENDIF
  ENDDO
ENDDO
!      FORMS THE VERTICAL 'LEAKINESS' **************************

DO JK=ILEVH,KLEV
  DO JL=KIDIA,KFDIA
    IF(JK >= KKENVH(JL)) THEN
      ZGGEENV=MAX(1.0_JPRB,&
       & (PGEOM1(JL,KKENVH(JL))+PGEOM1(JL,KKENVH(JL)-1))/2.0_JPRB)  
      ZGGEOM1=MAX(PGEOM1(JL,JK),1.0_JPRB)
      ZGVAR  =MAX(PHSTD(JL)*RG,1.0_JPRB)
      PZDEP(JL,JK)=SQRT((ZGGEENV-ZGGEOM1)/(ZGGEOM1+ZGVAR))
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('GWSETUP',1,ZHOOK_HANDLE)
END SUBROUTINE GWSETUP

SUBROUTINE GWPROFIL &
 & ( KIDIA ,KFDIA   ,KLON   , KLEV,&
 & LDGWD,&
 & KKCRITH, KCRIT , KKENVH, KKNU,KKNU2,&
 & PAPHM1, PRHO   , PSTAB , PVPH , PRI , PTAU,&
 & PDMOD , PSIG   , PHSTD      )  

!**** *GWPROFIL*

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!          FROM *GWDRAG*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===

!        *KIDIA*       START POINT
!        *KFDIA*       END POINT
!        *KLON*        NUMBER OF GRID POINTS PER PACKET
!        *KLEV*        NUMBER OF VERTICAL LEVELS
!        *LDGWD*       TRUE IF SSO ACTIVE
!        *KKCRITH*     MAXIMUM HALF-LEVEL FOR WAVE BREAKING
!        *KCRIT*       CRITICAL LEVEL
!        *KKENVH*      HALF-LEVEL FOR BLOCKING (TOP OF ENVELOPE LAYER)
!        *KKNU*        MODEL LEVEL AT 4*PHSTD HEIGHT 
!        *KKNU2*       MODEL LEVEL AT 3*PHSTD HEIGHT
!        *PAPHM1*      PRESSURE AT HALF LEVELS AT T-1
!        *PRHO*        AIR DENSITY
!        *PSTAB*       SQUARE OF THE BRUNT-VAISALA FREQUENCY
!        *PVPH*        WIND PROJECTION ALONG SURFACE STRESS AT HALF LEVELS
!        *PRI*         MEAN FLOW RICHARDSON NUMBER
!        *PTAU*        STRESS PRODUCED BY GRAVITY WAVES
!        *PDMOD*       MODULUS OF THE SURFACE STRESS (PART)
!        *PSIG*        MEAN SLOPE OF THE SUBGRID SCALE OROGRAPHY
!        *PHSTD*       STANDARD DEVIATION OF THE SUBGRID SCALE OROGRAPHY

!     ==== OUTPUTS ===

!        *PTAU*        STRESS PRODUCED BY GRAVITY WAVES

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD:
!     -------
!     THE STRESS PROFILE FOR GRAVITY WAVES IS COMPUTED AS FOLLOWS:
!     IT IS CONSTANT (NO GWD) AT THE LEVELS BETWEEN THE GROUND
!     AND THE TOP OF THE BLOCKED LAYER (KKENVH).
!     IT DECREASES LINEARLY WITH HEIGHTS FROM THE TOP OF THE 
!     BLOCKED LAYER TO 3*VAROR (KKNU), TO SIMULATES LEE WAVES OR 
!     NONLINEAR GRAVITY WAVE BREAKING.
!     ABOVE IT IS CONSTANT, EXCEPT WHEN THE WAVE ENCOUNTERS A CRITICAL
!     LEVEL (KCRIT) OR WHEN IT BREAKS.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        SEE ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S."

!     AUTHOR.
!     -------

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!-----------------------------------------------------------------------

!x USE PARKIND1  ,ONLY : JPIM     ,JPRB
!x USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!x 
!x USE YOEGWD   , ONLY : YREGWD
!x USE YOEPHLI  , ONLY : YREPHLI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
LOGICAL           ,INTENT(IN)    :: LDGWD(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KKCRITH(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCRIT(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KKENVH(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KKNU(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KKNU2(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) !note callpar indexing convention
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHO(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTAB(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRI(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTAU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMOD(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIG(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHSTD(KLON) 

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDZ2 (KLON,KLEV) , ZNORM(KLON) , ZORO(KLON)
REAL(KIND=JPRB) :: ZTAU (KLON,KLEV+1)

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: ZALFA, ZALPHA, ZB, ZDEL, ZDELP, ZDELPT, ZDZ2N, ZKDRAG, ZRIW, ZSQR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------

!     Constant for surface stress

IF (LHOOK) CALL DR_HOOK('GWPROFIL',0,ZHOOK_HANDLE)

IF (LPHYLIN) THEN
  ZKDRAG=RLPDRAG
ELSE
  ZKDRAG=GKDRAG
ENDIF

DO JL=KIDIA,KFDIA
  IF(LDGWD(JL)) THEN
    ZORO(JL)=PSIG(JL)*PDMOD(JL)/4._JPRB/MAX(PHSTD(JL),1.0_JPRB)
    ZTAU(JL,KKNU(JL)+1)=PTAU(JL,KKNU(JL)+1)
    ZTAU(JL,KLEV+1)    =PTAU(JL,KLEV+1)
  ENDIF
ENDDO

DO JK=KLEV,2,-1

!*         1.1    CONSTANT WAVE STRESS UNTIL TOP OF THE
!                 BLOCKING LAYER.

  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL)) THEN
      IF(JK >= KKNU2(JL)) THEN
        PTAU(JL,JK)=ZTAU(JL,KLEV+1)
      ENDIF
    ENDIF
  ENDDO

!*         1.2    WAVE DISPLACEMENT AT NEXT LEVEL.

!OCL  VCT(CEX)
  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL)) THEN
      IF(JK < KKNU2(JL)) THEN
        ZNORM(JL)=ZKDRAG*PRHO(JL,JK)*SQRT(PSTAB(JL,JK))*PVPH(JL,JK)*ZORO(JL)
        ZDZ2(JL,JK)=PTAU(JL,JK+1)/MAX(ZNORM(JL),GSSEC)
      ENDIF
    ENDIF
  ENDDO

!*         1.3    WAVE RICHARDSON NUMBER, NEW WAVE DISPLACEMENT
!*                AND STRESS:  BREAKING EVALUATION AND CRITICAL 
!                 LEVEL

  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL)) THEN
      IF(JK < KKNU2(JL)) THEN
        IF((PTAU(JL,JK+1) < GTSEC).OR.(JK <= KCRIT(JL))) THEN
          PTAU(JL,JK)=0.0_JPRB
        ELSE
          ZSQR=SQRT(PRI(JL,JK))
          ZALFA=SQRT(PSTAB(JL,JK)*ZDZ2(JL,JK))/PVPH(JL,JK)
          ZRIW=PRI(JL,JK)*(1.0_JPRB-ZALFA)/(1+ZALFA*ZSQR)**2
          IF(ZRIW < GRCRIT) THEN
            ZDEL=4._JPRB/ZSQR/GRCRIT+1.0_JPRB/GRCRIT**2+4._JPRB/GRCRIT
            ZB=1.0_JPRB/GRCRIT+2.0_JPRB/ZSQR
            ZALPHA=0.5_JPRB*(-ZB+SQRT(ZDEL))
            ZDZ2N=(PVPH(JL,JK)*ZALPHA)**2/PSTAB(JL,JK)
            PTAU(JL,JK)=ZNORM(JL)*ZDZ2N
          ELSE
            PTAU(JL,JK)=ZNORM(JL)*ZDZ2(JL,JK)
          ENDIF
          PTAU(JL,JK)=MIN(PTAU(JL,JK),PTAU(JL,JK+1))
        ENDIF
      ENDIF
    ENDIF
  ENDDO

ENDDO

!  REORGANISATION OF THE STRESS PROFILE
!  IF BREAKING OCCURS AT LOW LEVEL:
DO JL=KIDIA,KFDIA
  IF(LDGWD(JL)) THEN
    ZTAU(JL,KKENVH(JL)) =PTAU(JL,KKENVH(JL))
    ZTAU(JL,KKCRITH(JL))=PTAU(JL,KKCRITH(JL))
  ENDIF
ENDDO

DO JK=1,KLEV

  DO JL=KIDIA,KFDIA
    IF(LDGWD(JL)) THEN
      IF(JK > KKCRITH(JL).AND.JK < KKENVH(JL))THEN

        ZDELP=PAPHM1(JL,JK-1)-PAPHM1(JL,KKENVH(JL)-1)
        ZDELPT=PAPHM1(JL,KKCRITH(JL)-1)-PAPHM1(JL,KKENVH(JL)-1)
        PTAU(JL,JK)=ZTAU(JL,KKENVH(JL)) +&
         & (ZTAU(JL,KKCRITH(JL))-ZTAU(JL,KKENVH(JL)) )*&
         & ZDELP/ZDELPT  

      ENDIF
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('GWPROFIL',1,ZHOOK_HANDLE)
END SUBROUTINE GWPROFIL


END MODULE mo_sso_ifs

