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
! Turbulent orographic form drag

MODULE mo_vdftofdc
 
  PUBLIC :: vdftofdc

CONTAINS

SUBROUTINE VDFTOFDC(KIDIA,KFDIA,KLON,KLEV,PTMST,&
 & PUM1,PVM1,PGEOM1,PSIGFLT,&
 & PTOFDC)  
!     ------------------------------------------------------------------

!**   *VDFTOFDC* - DETERMINES THE COEFFICIENTS FOR THE 
!                 TURBULENT OROGRAPHIC DRAG PARAMETRIZATION

!     Original  A. BELJAARS   ECMWF    17/11/2002.
!     Modified 

!     PURPOSE
!     -------

!     DETERMINE COEFFICIENTS FOR TURBULENT OROGRAPHIC DRAG

!     INTERFACE
!     ---------

!     *VDFTOFDC* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KLEV*         NUMBER OF LEVELS

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PSIGFLT*      FILTERED STANDARD DEVIATION OF SUBGRID OROGRAPHY

!     OUTPUT PARAMETERS (REAL):

!     *PTOFDC*        COEFFICIENTS IN DIAGONAL TO BE PASSED ON 
!                    TO IMPLICIT SOLVER. PTOFDC=alpha*DT*stressdiv/PSI

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

! USE PARKIND1  ,ONLY : JPIM     ,JPRB
! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE YOMCST   , ONLY : RG
! USE YOEVDF   , ONLY : RVDIFTS, REPDU2

!ICON definitions:
USE mo_kind         ,ONLY : JPRB=>wp ,JPIM=>i4
USE mo_cuparameters ,ONLY : lhook    ,dr_hook  ,&
                 & RG      ,&                     !yomcst
                 & RVDIFTS ,REPDU2                !yoevdf

IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGFLT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDC(KLON,KLEV) 
!*            LOCAL STORAGE
!             ----- -------

REAL(KIND=JPRB) ::    ZCOEF(KLON)

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) ::    ZZ,ZTOFDALPHA,ZTOFDBETA,ZTOFDCMD,ZTOFDCORR,ZTOFDK1,ZTOFDIH,&
 & ZTOFDKFLT,ZTOFDN1,ZTOFDN2,ZTOFDIC,ZTOFDIZ,ZTOFDIN,ZFACT1,&
 & ZFACT2,ZUABS,ZMAX  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VDFTOFDC',0,ZHOOK_HANDLE)
ZTOFDALPHA=27.
ZTOFDBETA=1.
ZTOFDCMD=0.005
ZTOFDCORR=0.6
ZTOFDK1=0.003
ZTOFDIH=0.00102
ZTOFDKFLT=0.00035
ZTOFDN1=-1.9
ZTOFDN2=-2.8

ZTOFDIC=2.109
ZTOFDIZ=1500.
ZTOFDIN=-1.2

ZFACT1=ZTOFDK1**(ZTOFDN1-ZTOFDN2)/(ZTOFDIH*ZTOFDKFLT**ZTOFDN1)
ZFACT2=ZTOFDALPHA*ZTOFDBETA*ZTOFDCMD*ZTOFDCORR*ZTOFDIC*ZFACT1&
 & *RVDIFTS*PTMST  

ZMAX=5000.

!        1. PREPARE ARRAY INDEPENDENT OF HEIGHT
!           ------- ----- ----------- -- ------
DO JL=KIDIA,KFDIA
  ZCOEF(JL)=ZFACT2*PSIGFLT(JL)**2
ENDDO

!        2. VERTICAL LOOP
!           -------- ----

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZZ=PGEOM1(JL,JK)/RG
    IF (ZZ > ZMAX) THEN
      PTOFDC(JL,JK)=0.
    ELSE
      ZUABS=SQRT(MAX(REPDU2,PUM1(JL,JK)**2+PVM1(JL,JK)**2))
      PTOFDC(JL,JK)=ZCOEF(JL)*ZUABS*EXP(-(ZZ/ZTOFDIZ)**1.5)*ZZ**ZTOFDIN
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('VDFTOFDC',1,ZHOOK_HANDLE)
END SUBROUTINE VDFTOFDC


END MODULE mo_vdftofdc
