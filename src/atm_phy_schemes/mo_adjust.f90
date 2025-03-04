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

MODULE mo_adjust

  USE mo_kind   ,ONLY: jprb=>wp     , &
    & jpim=>i4

!  USE yomhook   ,ONLY : lhook,   dr_hook

   USE mo_cufunctions, ONLY: foealfa, foeewmcu, foeewm,     &
     & foeldcpmcu, foeewmcu,foedem,   &
     & foealfcu,foedemcu, foeldcpm,   &
     & foeles_v,foeies_v, foeewmcu_v, &
     & foeewm_v , foeewl, foeewi

   !  USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT

   !  USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
   !&                       R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,&
   !&                       R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,&
   !&                       RTICE    ,RTICECU  ,RTWAT_RTICE_R      ,&
   !&                       RTWAT_RTICECU_R

   USE mo_cuparameters, ONLY : r2es     ,r3les    ,r3ies    ,r4les    ,&
     &                         r4ies    ,r5alvcp                      ,&
     &                         r5alscp  ,ralvdcp  ,ralsdcp            ,&
     &                         lphylin  ,rlptrc   ,rlpal1   ,rlpal2   ,&
     &                         retv     ,rtt                          ,&
     &                         lhook    ,dr_hook                      ,&
     &                         vdiv     ,vexp     ,vrec

   IMPLICIT NONE

   PRIVATE



   !     ------------------------------------------------------------

   PUBLIC :: satur, cuadjtq, cuadjtqs



 CONTAINS

   SUBROUTINE satur ( kidia , kfdia , klon, ktdia  , klev,&
     & paprsf, pt, pqv, pqsat , kflag)

     !$ACC ROUTINE GANG

     !>
     !! Description:

     !! **   *SATUR* -  COMPUTES SPECIFIC HUMIDITY AT SATURATION
     !!       J.F. MAHFOUF       E.C.M.W.F.     15/05/96
     !!       Modified J. HAGUE          13/01/03 MASS Vector Functions

     !!       PURPOSE.
     !!       --------

     !!       SPECIFIC HUMIDITY AT SATURATION IS USED BY THE
     !!       DIAGNOSTIC CLOUD SCHEME TO COMPUTE RELATIVE HUMIDITY
     !!       AND LIQUID WATER CONTENT

     !!
     !! Code Description:
     !!       PARAMETER     DESCRIPTION                                 UNITS
     !!       ---------     -----------                                 -----
     !!       INPUT PARAMETERS (INTEGER):

     !!      *KIDIA*        START POINT
     !!      *KFDIA*        END POINT
     !!      *KLON*         NUMBER OF GRID POINTS PER PACKET
     !!      *KTDIA*        START OF THE VERTICAL LOOP
     !!      *KLEV*         NUMBER OF LEVELS

     !!       INPUT PARAMETERS (REAL):

     !!      *PAPRSF*        PRESSURE ON FULL LEVELS                      PA
     !!      *PT*            TEMPERATURE AT T-DT                          K
     !!      *PQV*           specific humidity AT T-DT                          K

     !!       INPUT PARAMETERS (INTEGER):

     !!      *KFLAG*         FLAG TO DETECT CALL FROM

     !!                      CONVECTION  KFLAG=1
     !!                      OTHER       KFLAG=2

     !!       OUTPUT PARAMETER (REAL):

     !!      *PQSAT*         SATURATION SPECIFIC HUMIDITY                 KG/KG

     !-------------------------------------------------------------------------

     !USE PARKIND1  ,ONLY : JPIM     ,JPRB
     !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

     !USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT
     !USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
     ! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
     ! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
     ! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R

     !USE YOEPHLI  , ONLY : LPHYLIN
     !USE YOMJFH   , ONLY : N_VMASS

     IMPLICIT NONE

     INTEGER(KIND=jpim),INTENT(in)    :: klon
     INTEGER(KIND=jpim),INTENT(in)    :: klev
     INTEGER(KIND=jpim),INTENT(in)    :: kidia
     INTEGER(KIND=jpim),INTENT(in)    :: kfdia
     INTEGER(KIND=jpim),INTENT(in)    :: ktdia
     REAL(KIND=jprb)   ,INTENT(in)    :: paprsf(klon,klev)
     REAL(KIND=jprb)   ,INTENT(in)    :: pt(klon,klev)
     REAL(KIND=jprb)   ,INTENT(in)    :: pqv(klon,klev)
     REAL(KIND=jprb)   ,INTENT(inout) :: pqsat(klon,klev)
     INTEGER(KIND=jpim),INTENT(in)    :: kflag

     INTEGER(KIND=jpim) :: jk, jl, jlen

     REAL(KIND=jprb) :: zcor, zew, zfoeew, zqmax, zqs, ztarg
     REAL(KIND=jprb) :: zalfa !, zfoeewl, zfoeewi
     REAL(KIND=jprb) :: z_exparg1(kidia:kfdia)
     REAL(KIND=jprb) :: z_exparg2(kidia:kfdia)
     REAL(KIND=jprb) :: z_expout1(kidia:kfdia)
     REAL(KIND=jprb) :: z_expout2(kidia:kfdia)
     REAL(KIND=jprb) :: zhook_handle

     !#include "fcttre.h"

     !----------------------------------------------------------------------

     !>*    1.           DEFINE CONSTANTS
     !!                  ----------------


     IF (lhook) CALL dr_hook('SATUR',0,zhook_handle)

     zqmax=0.5_JPRB

     !     *
     !----------------------------------------------------------------------

     !!     *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
     !!                       --------------------------------------

     IF (lphylin) THEN     ! linear physics: set to .FALSE. in mo_cuparameters.f90

       !$ACC LOOP SEQ
       DO jk=ktdia,klev
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ztarg, zalfa, zfoeew, zqs, zcor)
         DO jl=kidia, kfdia
           ztarg = pt(jl,jk)
           zalfa = foealfa(ztarg)

           !> KF use functions
           !     ZFOEEWL = R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
           !      ZFOEEWI = R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
           !      ZFOEEW = ZALFA*ZFOEEWL+(1.0_JPRB-ZALFA)*ZFOEEWI
           zfoeew = zalfa*foeewl(ztarg)+(1.0_JPRB-zalfa)*foeewi(ztarg)
           !<KF

           zqs    = zfoeew/paprsf(jl,jk)
           IF (zqs > zqmax) THEN
             zqs=zqmax
           ENDIF
           zcor = 1.0_JPRB/(1.0_JPRB-retv*pqv(jl,jk))
           pqsat(jl,jk)=zqs*zcor
         ENDDO
       ENDDO

     ELSE

       !$ACC LOOP SEQ
       DO jk=ktdia,klev
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zew, zqs, zcor)
         DO jl=kidia, kfdia
           IF(kflag == 1) THEN
             zew  = foeewmcu(pt(jl,jk))
           ELSE
             zew  = foeewm(pt(jl,jk))
           ENDIF
           zqs  = zew/paprsf(jl,jk)
           zqs  = MIN(zqmax,zqs)
           ! Modification, GZ (2014-07-23): define qv_sat as qv/RH, implying that the qv_sat in the 
           ! denominator needs to be replaced with qv
           zcor = 1.0_JPRB/(1.0_JPRB-retv*pqv(jl,jk))
           pqsat(jl,jk)=zqs*zcor
         ENDDO
       ENDDO

     ENDIF


     IF (lhook) CALL dr_hook('SATUR',1,zhook_handle)


   END SUBROUTINE satur


   !!#ifdef RS6K
   !!@PROCESS HOT NOSTRICT
   !!#endif
   !
   SUBROUTINE cuadjtq &
     & (kidia,    kfdia,    klon,    klev,&
     & kk,&
     & psp,      pt,       pq,       ldflag,   kcall)

     !$ACC ROUTINE GANG

     !!
     !! Description:
     !!          M.TIEDTKE         E.C.M.W.F.     12/89

     !!          MODIFICATIONS
     !!          -------------
     !!          D.SALMOND         CRAY(UK))      12/8/91
     !!          J.J. MORCRETTE    ECMWF          92-09-18   Update to Cy44
     !!          J.F. MAHFOUF      ECMWF          96-06-11   Smoothing option
     !!          D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
     !!          J.HAGUE                          03-01-13   MASS Vector Functions
     !!          J.HAGUE                          03-07-07   More MASS V.F.
     !!        M.Hamrud              01-Oct-2003 CY28 Cleaning
     !!        J.Hague & D.Salmond   22-Nov-2005 Optimisations

     !!          PURPOSE.
     !!          --------
     !!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
     !!
     !!          INPUT ARE UNADJUSTED T AND Q VALUES,
     !!         IT RETURNS ADJUSTED VALUES OF T AND Q
     !!
     !! Code Description:
     !!     PARAMETER     DESCRIPTION                                   UNITS
     !!     ---------     -----------                                   -----
     !!     INPUT PARAMETERS (INTEGER):

     !!    *KIDIA*        START POINT
     !!    *KFDIA*        END POINT
     !!    *KLON*         NUMBER OF GRID POINTS PER PACKET
     !!   *KLEV*         NUMBER OF LEVELS
     !!    *KK*           LEVEL
     !!    *KCALL*        DEFINES CALCULATION AS
     !!                     KCALL=0  ENV. T AND QS IN*CUINI*
     !!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
     !!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)

     !!    INPUT PARAMETERS (LOGICAL):

     !!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)

     !!     INPUT PARAMETERS (REAL):

     !!    *PSP*          PRESSURE                                        PA

     !!     UPDATED PARAMETERS (REAL):

     !!    *PT*           TEMPERATURE                                     K
     !!    *PQ*           SPECIFIC HUMIDITY                             KG/KG

     !!          EXTERNALS
     !!          ---------
     !!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
     !!         FOR CONDENSATION CALCULATIONS.
     !!         THE TABLES ARE INITIALISED IN *SUPHEC*.
     !
     !----------------------------------------------------------------------

     !USE PARKIND1  ,ONLY : JPIM     ,JPRB
     !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

     !USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT
     !USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
     ! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
     ! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
     ! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R
     !USE YOEPHLI  , ONLY : LPHYLIN  ,RLPTRC   ,RLPAL1   ,RLPAL2
     !USE YOMJFH   , ONLY : N_VMASS

     IMPLICIT NONE

     INTEGER(KIND=jpim),INTENT(in)    :: klon
     INTEGER(KIND=jpim),INTENT(in)    :: klev
     INTEGER(KIND=jpim),INTENT(in)    :: kidia
     INTEGER(KIND=jpim),INTENT(in)    :: kfdia
     INTEGER(KIND=jpim),INTENT(in)    :: kk
     REAL(KIND=jprb)   ,INTENT(in)    :: psp(klon)
     REAL(KIND=jprb)   ,INTENT(inout) :: pt(klon,klev)
     REAL(KIND=jprb)   ,INTENT(inout) :: pq(klon,klev)
     LOGICAL           ,INTENT(in)    :: ldflag(klon)
     INTEGER(KIND=jpim),INTENT(in)    :: kcall
     INTEGER(KIND=jpim) :: jl, jlen

     REAL(KIND=jprb) :: z1s, z2s, zcond,zcond1, zcor,      &
       &                zoealfa, zqmax, zqsat, ztarg, zqp
     REAL(KIND=jprb) :: pt1, pq1
!     REAL(KIND=jprb) :: zfoeewi, zfoeewl

     REAL(KIND=jprb) :: zl, zi, zf

     !#include "fcttre.h"

     !     STATEMENT FUNCTIONS
     !REAL_B :: FOEALFAJ,FOEDEMJ,FOELDCPMJ,FOEEWMJ

     REAL(KIND=jprb) :: minj, x, y
     REAL(KIND=jprb) :: zhook_handle

     minj(x,y) = y - 0.5_JPRB*(ABS(x-y)-(x-y))

     !----------------------------------------------------------------------

     !     1.           DEFINE CONSTANTS
     !                  ----------------


     IF (lhook) CALL dr_hook('CUADJTQ',0,zhook_handle)


     zqmax=0.5_JPRB


     !*********************************************
     IF (.NOT.lphylin) THEN
     !*********************************************

       !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
       !                  -----------------------------------------------------


       IF (kcall == 1 ) THEN

!DIR$ IVDEP
!$NEC sparse
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqp, zl, zi, zqsat, zcor, zf, zcond, zcond1)
         DO jl=kidia,kfdia
           IF(ldflag(jl)) THEN
             zqp    =1.0_JPRB/psp(jl)
             !       ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
             ! FOEEWMCU ( PTARE ) = R2ES *&
             !  &(FOEALFCU(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
             !  &(1.0_JPRB-FOEALFCU(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
             zl=1.0_JPRB/(pt(jl,kk)-r4les)
             zi=1.0_JPRB/(pt(jl,kk)-r4ies)
             zqsat=r2es *(foealfcu(pt(jl,kk))*EXP(r3les*(pt(jl,kk)-rtt)*zl)+&
               & (1.0_JPRB-foealfcu(pt(jl,kk)))*EXP(r3ies*(pt(jl,kk)-rtt)*zi))
             zqsat=zqsat*zqp
             zqsat=MIN(0.5_JPRB,zqsat)
             zcor=1.0_JPRB-retv*zqsat
             !       ZCOND=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
             ! FOEDEMCU ( PTARE )=FOEALFCU(PTARE)*R5ALVCP*(1.0_JPRB/(PTARE-R4LES)**2)+&
             !   &(1.0_JPRB-FOEALFCU(PTARE))*R5ALSCP*(1.0_JPRB/(PTARE-R4IES)**2)
             zf=foealfcu(pt(jl,kk))*r5alvcp*zl**2 + &
               & (1.0_JPRB-foealfcu(pt(jl,kk)))*r5alscp*zi**2
             zcond=(pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
             !       ZCOND=MAX(ZCOND,0.0_JPRB)
             IF(zcond > 0.0_JPRB)THEN
               pt(jl,kk)=pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond
               pq(jl,kk)=pq(jl,kk)-zcond
               !         ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
               zl=1.0_JPRB/(pt(jl,kk)-r4les)
               zi=1.0_JPRB/(pt(jl,kk)-r4ies)
               zqsat=r2es *(foealfcu(pt(jl,kk))*EXP(r3les*(pt(jl,kk)-rtt)*zl)+&
                 & (1.0_JPRB-foealfcu(pt(jl,kk)))*EXP(r3ies*(pt(jl,kk)-rtt)*zi))
               zqsat=zqsat*zqp
               zqsat=minj(0.5_JPRB,zqsat)
               zcor=1.0_JPRB-retv*zqsat
               !         ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
               zf=foealfcu(pt(jl,kk))*r5alvcp*zl**2 + &
                 & (1.0_JPRB-foealfcu(pt(jl,kk)))*r5alscp*zi**2
               zcond1=(pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
               IF(zcond ==  0.0_JPRB)zcond1=0.0_JPRB
               pt(jl,kk)=pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
               pq(jl,kk)=pq(jl,kk)-zcond1
             ENDIF
           ENDIF
         ENDDO
       ENDIF  ! kcall == 1

       IF(kcall == 2) THEN

!DIR$ IVDEP
!OCL NOVREC
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqp, zqsat, zcor, zcond, zcond1)
         DO jl=kidia,kfdia
           IF(ldflag(jl)) THEN
             zqp    =1.0_JPRB/psp(jl)
             zqsat=foeewmcu(pt(jl,kk))*zqp
             zqsat=MIN(0.5_JPRB,zqsat)
             zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
             zqsat=zqsat*zcor
             zcond=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedemcu(pt(jl,kk)))
             zcond=MIN(zcond,0.0_JPRB)
             pt(jl,kk)=pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond
             pq(jl,kk)=pq(jl,kk)-zcond
             zqsat=foeewmcu(pt(jl,kk))*zqp
             zqsat=MIN(0.5_JPRB,zqsat)
             zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
             zqsat=zqsat*zcor
             zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedemcu(pt(jl,kk)))
             IF(zcond == 0.0_JPRB)zcond1=MIN(zcond1,0.0_JPRB)
             pt(jl,kk)=pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
             pq(jl,kk)=pq(jl,kk)-zcond1
           ENDIF
         ENDDO
       ENDIF  ! kcall == 2

       IF(kcall == 0) THEN

!DIR$ IVDEP
!OCL NOVREC
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqp, zqsat, zcor, zcond1)
         DO jl=kidia,kfdia
           zqp    =1.0_JPRB/psp(jl)
           zqsat=foeewm(pt(jl,kk))*zqp
           zqsat=MIN(0.5_JPRB,zqsat)
           zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
           zqsat=zqsat*zcor
           zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedem(pt(jl,kk)))
           pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
           zqsat=foeewm(pt(jl,kk))*zqp
           zqsat=MIN(0.5_JPRB,zqsat)
           zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
           zqsat=zqsat*zcor
           zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedem(pt(jl,kk)))
           pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
         ENDDO
       ENDIF  ! kcall == 0

       IF(kcall == 4 )THEN
        
!DIR$ IVDEP
!OCL NOVREC
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqp, zqsat, zcor, zcond, zcond1)
         DO jl=kidia,kfdia
           IF(ldflag(jl)) THEN
             zqp    =1.0_JPRB/psp(jl)
             zqsat=foeewm(pt(jl,kk))*zqp
             zqsat=MIN(0.5_JPRB,zqsat)
             zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
             zqsat=zqsat*zcor
             zcond=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedem(pt(jl,kk)))
             pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
             pq(jl,kk)=pq(jl,kk)-zcond
             zqsat=foeewm(pt(jl,kk))*zqp
             zqsat=MIN(0.5_JPRB,zqsat)
             zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
             zqsat=zqsat*zcor
             zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedem(pt(jl,kk)))
             pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
             pq(jl,kk)=pq(jl,kk)-zcond1
           ENDIF
         ENDDO
       ENDIF  ! kcall == 4
      
       IF(kcall == 5) THEN  ! Same as 4 but with LDFLAG all true
        
!DIR$ IVDEP
!OCL NOVREC
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqp, zqsat, zcor, zcond, zcond1)
         DO jl=kidia,kfdia
           zqp    =1.0_JPRB/psp(jl)
           zqsat=foeewm(pt(jl,kk))*zqp
           zqsat=MIN(0.5_JPRB,zqsat)
           zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
           zqsat=zqsat*zcor
           zcond=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedem(pt(jl,kk)))
           pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
           pq(jl,kk)=pq(jl,kk)-zcond
           zqsat=foeewm(pt(jl,kk))*zqp
           zqsat=MIN(0.5_JPRB,zqsat)
           zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
           zqsat=zqsat*zcor
           zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedem(pt(jl,kk)))
           pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
         ENDDO

       ENDIF     ! kcall == 5
      
       IF(kcall == 3) THEN

!DIR$ IVDEP
!OCL NOVREC
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqp, zqsat, zcor, zcond1, pt1, pq1)
         DO jl=kidia,kfdia
           zqp    =1.0_JPRB/psp(jl)
           zqsat=foeewmcu(pt(jl,kk))*zqp
           zqsat=MIN(0.5_JPRB,zqsat)
           zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
           zqsat=zqsat*zcor
           zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*foedemcu(pt(jl,kk)))
           pt1=pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
           pq1=pq(jl,kk)-zcond1
           zqsat=foeewmcu(pt1)*zqp
           zqsat=MIN(0.5_JPRB,zqsat)
           zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
           zqsat=zqsat*zcor
           zcond1=(pq1-zqsat)/(1.0_JPRB+zqsat*zcor*foedemcu(pt1))
           pt(jl,kk)=pt1+foeldcpmcu(pt1)*zcond1
           pq(jl,kk)=pq1-zcond1
         ENDDO
        
       ENDIF     ! kcall == 3

     !*********************************************
     ELSE   ! lphylin
     !*********************************************
     ! US not ported to GPUs, because lphylin is set to .FALSE. 

      !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
      !                  -----------------------------------------------------
      
      IF (kcall == 1 ) THEN
        
!DIR$ IVDEP
!OCL NOVREC
        DO jl=kidia,kfdia
          IF(ldflag(jl)) THEN
            zqp    =1.0_JPRB/psp(jl)
            ztarg=pt(jl,kk)
            zoealfa=0.5_JPRB*(TANH(rlpal1*(ztarg-rlptrc))+1.0_JPRB)
            !> KF use available functions
            !       ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
            !        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
            !       ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
            zqsat=zqp    *(zoealfa*foeewl(ztarg)+(1.0_JPRB-zoealfa)*foeewi(ztarg))
            !<KF
            
            z1s=TANH(rlpal2*(zqsat-zqmax))
            zqsat=0.5_JPRB*((1.0_JPRB-z1s)*zqsat+(1.0_JPRB+z1s)*zqmax)
            
            zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
            zqsat=zqsat*zcor
            
            z2s=    zoealfa *r5alvcp*(1.0_JPRB/(ztarg-r4les)**2)+&
              & (1.0_JPRB-zoealfa)*r5alscp*(1.0_JPRB/(ztarg-r4ies)**2)
            zcond=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*z2s)
            
            zcond=MAX(zcond,0.0_JPRB)
            
            IF(zcond /= 0.0_JPRB) THEN
              
              pt(jl,kk)=pt(jl,kk)+&
                & (zoealfa*ralvdcp+(1.0_JPRB-zoealfa)*ralsdcp)*zcond
              pq(jl,kk)=pq(jl,kk)-zcond
              ztarg=pt(jl,kk)
              zoealfa=0.5_JPRB*(TANH(rlpal1*(ztarg-rlptrc))+1.0_JPRB)
              !> KF use available functions
              !          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
              !          ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
              !          ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
              zqsat=zqp    *(zoealfa*foeewl(ztarg)+(1.0_JPRB-zoealfa)*foeewi(ztarg))
              !<KF
              z1s=TANH(rlpal2*(zqsat-zqmax))
              zqsat=0.5_JPRB*((1.0_JPRB-z1s)*zqsat+(1.0_JPRB+z1s)*zqmax)
              
              zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
              zqsat=zqsat*zcor
              
              z2s=    zoealfa *r5alvcp*(1.0_JPRB/(ztarg-r4les)**2)+&
                & (1.0_JPRB-zoealfa)*r5alscp*(1.0_JPRB/(ztarg-r4ies)**2)
              zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*z2s)
              
              pt(jl,kk)=pt(jl,kk)+(zoealfa*ralvdcp+(1.0_JPRB-zoealfa)*ralsdcp)*zcond1
              
              pq(jl,kk)=pq(jl,kk)-zcond1
            ENDIF
          ENDIF
        ENDDO
        
      ENDIF
      
      IF(kcall == 2) THEN
        
!DIR$ IVDEP
!OCL NOVREC
        DO jl=kidia,kfdia
          IF(ldflag(jl)) THEN
            zqp    =1.0_JPRB/psp(jl)
            
            ztarg=pt(jl,kk)
            zoealfa=0.5_JPRB*(TANH(rlpal1*(ztarg-rlptrc))+1.0_JPRB)
            !>KF
            !        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
            !        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
            !        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
            zqsat=zqp *(zoealfa*foeewl(ztarg)+(1.0_JPRB-zoealfa)*foeewi(ztarg))
            !<KF
            z1s=TANH(rlpal2*(zqsat-zqmax))
            zqsat=0.5_JPRB*((1.0_JPRB-z1s)*zqsat+(1.0_JPRB+z1s)*zqmax)
            
            zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
            zqsat=zqsat*zcor
            
            z2s=    zoealfa *r5alvcp*(1.0_JPRB/(ztarg-r4les)**2)+&
              & (1.0_JPRB-zoealfa)*r5alscp*(1.0_JPRB/(ztarg-r4ies)**2)
            zcond=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*z2s)
            
            zcond=MIN(zcond,0.0_JPRB)
            
            IF(zcond /= 0.0_JPRB) THEN
              
              pt(jl,kk)=pt(jl,kk)+&
                & (zoealfa*ralvdcp+(1.0_JPRB-zoealfa)*ralsdcp)*zcond
              pq(jl,kk)=pq(jl,kk)-zcond
              ztarg=pt(jl,kk)
              zoealfa=0.5_JPRB*(TANH(rlpal1*(ztarg-rlptrc))+1.0_JPRB)
              !>KF
              !          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
              !          ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
              !          ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
              zqsat=zqp*(zoealfa*foeewl(ztarg)+(1.0_JPRB-zoealfa)*foeewi(ztarg))
              !<KF
              z1s=TANH(rlpal2*(zqsat-zqmax))
              zqsat=0.5_JPRB*((1.0_JPRB-z1s)*zqsat+(1.0_JPRB+z1s)*zqmax)
              
              zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
              zqsat=zqsat*zcor
              
              z2s=    zoealfa *r5alvcp*(1.0_JPRB/(ztarg-r4les)**2)+&
                & (1.0_JPRB-zoealfa)*r5alscp*(1.0_JPRB/(ztarg-r4ies)**2)
              zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*z2s)
              
              pt(jl,kk)=pt(jl,kk)+(zoealfa*ralvdcp+(1.0_JPRB-zoealfa)*ralsdcp)*zcond1
              
              pq(jl,kk)=pq(jl,kk)-zcond1
            ENDIF
          ENDIF
        ENDDO
        
      ENDIF
      
      IF(kcall == 0) THEN
        
!DIR$ IVDEP
!OCL NOVREC
        DO jl=kidia,kfdia
          zqp    =1.0_JPRB/psp(jl)
          
          ztarg=pt(jl,kk)
          zoealfa=0.5_JPRB*(TANH(rlpal1*(ztarg-rlptrc))+1.0_JPRB)
          !>KF
          !      ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
          !      ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
          !      ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
          zqsat=zqp*(zoealfa*foeewl(ztarg)+(1.0_JPRB-zoealfa)*foeewi(ztarg))
          !>KF
          z1s=TANH(rlpal2*(zqsat-zqmax))
          zqsat=0.5_JPRB*((1.0_JPRB-z1s)*zqsat+(1.0_JPRB+z1s)*zqmax)
          
          zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
          zqsat=zqsat*zcor
          
          z2s=    zoealfa *r5alvcp*(1.0_JPRB/(ztarg-r4les)**2)+&
            & (1.0_JPRB-zoealfa)*r5alscp*(1.0_JPRB/(ztarg-r4ies)**2)
          zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*z2s)
          
          pt(jl,kk)=pt(jl,kk)+(zoealfa*ralvdcp+(1.0_JPRB-zoealfa)*ralsdcp)*zcond1
          
          pq(jl,kk)=pq(jl,kk)-zcond1
          
          ztarg=pt(jl,kk)
          zoealfa=0.5_JPRB*(TANH(rlpal1*(ztarg-rlptrc))+1.0_JPRB)
          !>KF
          !      ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
          !      ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
          !      ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
          zqsat=zqp*(zoealfa*foeewl(ztarg)+(1.0_JPRB-zoealfa)*foeewi(ztarg))
          !KF
          z1s=TANH(rlpal2*(zqsat-zqmax))
          zqsat=0.5_JPRB*((1.0_JPRB-z1s)*zqsat+(1.0_JPRB+z1s)*zqmax)
          
          zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
          zqsat=zqsat*zcor
          
          z2s=    zoealfa *r5alvcp*(1.0_JPRB/(ztarg-r4les)**2)+&
            & (1.0_JPRB-zoealfa)*r5alscp*(1.0_JPRB/(ztarg-r4ies)**2)
          zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*z2s)
          
          pt(jl,kk)=pt(jl,kk)+(zoealfa*ralvdcp+(1.0_JPRB-zoealfa)*ralsdcp)*zcond1
          
          pq(jl,kk)=pq(jl,kk)-zcond1
        ENDDO
        
      ENDIF
      
      IF(kcall == 4) THEN
        
!DIR$ IVDEP
!OCL NOVREC
        DO jl=kidia,kfdia
          IF(ldflag(jl)) THEN
            zqp    =1.0_JPRB/psp(jl)
            
            ztarg=pt(jl,kk)
            zoealfa=0.5_JPRB*(TANH(rlpal1*(ztarg-rlptrc))+1.0_JPRB)
            !>KF
            !        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
            !        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
            !        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
            zqsat=zqp *(zoealfa*foeewl(ztarg)+(1.0_JPRB-zoealfa)*foeewi(ztarg))
            !KF
            z1s=TANH(rlpal2*(zqsat-zqmax))
            zqsat=0.5_JPRB*((1.0_JPRB-z1s)*zqsat+(1.0_JPRB+z1s)*zqmax)
            
            zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
            zqsat=zqsat*zcor
            
            z2s=    zoealfa *r5alvcp*(1.0_JPRB/(ztarg-r4les)**2)+&
              & (1.0_JPRB-zoealfa)*r5alscp*(1.0_JPRB/(ztarg-r4ies)**2)
            zcond=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*z2s)
            
            pt(jl,kk)=pt(jl,kk)+(zoealfa*ralvdcp+(1.0_JPRB-zoealfa)*ralsdcp)*zcond
            
            pq(jl,kk)=pq(jl,kk)-zcond
            
            ztarg=pt(jl,kk)
            zoealfa=0.5_JPRB*(TANH(rlpal1*(ztarg-rlptrc))+1.0_JPRB)
            !>KF
            !        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
            !        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
            !        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_JPRB-ZOEALFA)*ZFOEEWI)
            zqsat=zqp *(zoealfa*foeewl(ztarg)+(1.0_JPRB-zoealfa)*foeewi(ztarg))
            !KF
            z1s=TANH(rlpal2*(zqsat-zqmax))
            zqsat=0.5_JPRB*((1.0_JPRB-z1s)*zqsat+(1.0_JPRB+z1s)*zqmax)
            
            zqsat=MIN(zqmax,zqsat)
            zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsat)
            zqsat=zqsat*zcor
            
            z2s=    zoealfa *r5alvcp*(1.0_JPRB/(ztarg-r4les)**2)+&
              & (1.0_JPRB-zoealfa)*r5alscp*(1.0_JPRB/(ztarg-r4ies)**2)
            zcond1=(pq(jl,kk)-zqsat)/(1.0_JPRB+zqsat*zcor*z2s)
            
            pt(jl,kk)=pt(jl,kk)+(zoealfa*ralvdcp+(1.0_JPRB-zoealfa)*ralsdcp)*zcond1
            
            pq(jl,kk)=pq(jl,kk)-zcond1
          ENDIF
        ENDDO
        
      ENDIF
      
     !*********************************************
     ENDIF  ! lphylin
     !*********************************************
    

     IF (lhook) CALL dr_hook('CUADJTQ',1,zhook_handle)


  END SUBROUTINE cuadjtq
  
  !!KF  NOTE: The following routine is only called in case of linearized Physics.
  !! Therefore no modifications for DWD purposes are made there
  !! And also no GPU port
  
  !
  SUBROUTINE cuadjtqs &
    & (kidia,    kfdia,    klon,    klev,&
    & kk,&
    & psp,      pt,       pq,       ldflag,   kcall)
    !
    !! Description:
    !!**   *CUADJTQS* - SIMPLIFIED VERSION OF MOIST ADJUSTMENT
    !!     J.F. MAHFOUF      ECMWF
    !!
    !!     PURPOSE.
    !!     --------
    !!     TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
    
    !! History:
    !!          MODIFICATIONS
    !!          -------------
    !!         D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
    !!        M.Hamrud      01-Oct-2003 CY28 Cleaning
    
    !!     INPUT ARE UNADJUSTED T AND Q VALUES,
    !!     IT RETURNS ADJUSTED VALUES OF T AND Q
    
    !!     PARAMETER     DESCRIPTION                                   UNITS
    !!     ---------     -----------                                   -----
    !!     INPUT PARAMETERS (INTEGER):
    
    !!    *KIDIA*        START POINT
    !!    *KFDIA*        END POINT
    !!    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!   *KLEV*         NUMBER OF LEVELS
    !!    *KK*           LEVEL
    !!    *KCALL*        DEFINES CALCULATION AS
    !!                      KCALL=0  ENV. T AND QS IN*CUINI*
    !!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
    !!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
    
    !!     INPUT PARAMETERS (LOGICAL):
    
    !!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)
    
    !!     INPUT PARAMETERS (REAL):
    
    !!    *PSP*          PRESSURE                                        PA
    
    !!     UPDATED PARAMETERS (REAL):
    
    !!    *PT*           TEMPERATURE                                     K
    !!    *PQ*           SPECIFIC HUMIDITY                             KG/KG
    
    !!          MODIFICATIONS
    !!          -------------
    !!         D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
    !!        M.Hamrud      01-Oct-2003 CY28 Cleaning
    
    !----------------------------------------------------------------------
    
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
    
    !USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT
    !USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
    ! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
    ! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
    ! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R
    
    IMPLICIT NONE
    
    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: kk
    REAL(KIND=jprb)   ,INTENT(in)    :: psp(klon)
    REAL(KIND=jprb)   ,INTENT(inout) :: pt(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pq(klon,klev)
    LOGICAL ,INTENT(in)    :: ldflag(klon)
    INTEGER(KIND=jpim),INTENT(in)    :: kcall
    REAL(KIND=jprb) ::     z3es(klon),             z4es(klon),&
      & z5alcp(klon),           zaldcp(klon)
    
    INTEGER(KIND=jpim) :: jl
    
    REAL(KIND=jprb) :: zfoeew
    REAL(KIND=jprb) :: zqmax, zqp, zcond, zcond1, ztarg, zcor, zqsat,  z2s
    REAL(KIND=jprb) :: zhook_handle
    
    !#include "fcttre.h"
    !----------------------------------------------------------------------
    
    !     1.           DEFINE CONSTANTS
    !                  ----------------
    
    IF (lhook) CALL dr_hook('CUADJTQS',0,zhook_handle)
    zqmax=0.5_JPRB
    
    !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
    !                  -----------------------------------------------------
    
    !*    ICE-WATER THERMODYNAMICAL FUNCTIONS
    
    DO jl=kidia,kfdia
      IF (pt(jl,kk) > rtt) THEN
        z3es(jl)=r3les
        z4es(jl)=r4les
        z5alcp(jl)=r5alvcp
        zaldcp(jl)=ralvdcp
      ELSE
        z3es(jl)=r3ies
        z4es(jl)=r4ies
        z5alcp(jl)=r5alscp
        zaldcp(jl)=ralsdcp
      ENDIF
    ENDDO
    
    IF (kcall == 1 ) THEN
      
!DIR$ IVDEP
!OCL NOVREC
      DO jl=kidia,kfdia
        IF(ldflag(jl)) THEN
          zqp    =1.0_JPRB/psp(jl)
          ztarg    =pt(jl,kk)
          zfoeew    =r2es*EXP(z3es(jl)*(ztarg    -rtt)/(ztarg    -z4es(jl)))
          zqsat    =zqp    *zfoeew
          IF (zqsat     > zqmax) THEN
            zqsat    =zqmax
          ENDIF
          zcor    =1.0_JPRB/(1.0_JPRB-retv*zqsat    )
          zqsat    =zqsat    *zcor
          z2s    =z5alcp(jl)/(ztarg    -z4es(jl))**2
          zcond    =(pq(jl,kk)-zqsat    )/(1.0_JPRB+zqsat    *zcor    *z2s    )
          zcond    =MAX(zcond    ,0.0_JPRB)
          !     IF(ZCOND /= _ZERO_) THEN
          pt(jl,kk)=pt(jl,kk)+zaldcp(jl)*zcond
          pq(jl,kk)=pq(jl,kk)-zcond
          ztarg    =pt(jl,kk)
          zfoeew    =r2es*EXP(z3es(jl)*(ztarg    -rtt)/(ztarg    -z4es(jl)))
          zqsat    =zqp    *zfoeew
          IF (zqsat     > zqmax) THEN
            zqsat    =zqmax
          ENDIF
          zcor    =1.0_JPRB/(1.0_JPRB-retv*zqsat    )
          zqsat    =zqsat    *zcor
          z2s    =z5alcp(jl)/(ztarg    -z4es(jl))**2
          zcond1    =(pq(jl,kk)-zqsat    )/(1.0_JPRB+zqsat    *zcor    *z2s    )
          IF(zcond ==  0.0_JPRB)zcond1=0.0_JPRB
          pt(jl,kk)=pt(jl,kk)+zaldcp(jl)*zcond1
          pq(jl,kk)=pq(jl,kk)-zcond1
          !     ENDIF
        ENDIF
      ENDDO
      
    ENDIF
    
    IF(kcall == 2) THEN
      
!DIR$ IVDEP
!OCL NOVREC
      DO jl=kidia,kfdia
        IF(ldflag(jl)) THEN
          zqp    =1.0_JPRB/psp(jl)
          ztarg    =pt(jl,kk)
          zfoeew    =r2es*EXP(z3es(jl)*(ztarg    -rtt)/(ztarg    -z4es(jl)))
          zqsat    =zqp    *zfoeew
          IF (zqsat     > zqmax) THEN
            zqsat    =zqmax
          ENDIF
          zcor    =1.0_JPRB/(1.0_JPRB-retv  *zqsat    )
          zqsat    =zqsat    *zcor
          z2s    =z5alcp(jl)/(ztarg    -z4es(jl))**2
          zcond    =(pq(jl,kk)-zqsat    )/(1.0_JPRB+zqsat    *zcor    *z2s    )
          zcond    =MIN(zcond    ,0.0_JPRB)
          !     IF(ZCOND /= _ZERO_) THEN
          pt(jl,kk)=pt(jl,kk)+zaldcp(jl)*zcond
          pq(jl,kk)=pq(jl,kk)-zcond
          ztarg    =pt(jl,kk)
          zfoeew    =r2es*EXP(z3es(jl)*(ztarg    -rtt)/(ztarg    -z4es(jl)))
          zqsat    =zqp    *zfoeew
          IF (zqsat     > zqmax) THEN
            zqsat    =zqmax
          ENDIF
          zcor    =1.0_JPRB/(1.0_JPRB-retv  *zqsat    )
          zqsat    =zqsat    *zcor
          z2s    =z5alcp(jl)/(ztarg    -z4es(jl))**2
          zcond1    =(pq(jl,kk)-zqsat    )/(1.0_JPRB+zqsat    *zcor    *z2s    )
          IF(zcond ==  0.0_JPRB)zcond1=0.0_JPRB
          pt(jl,kk)=pt(jl,kk)+zaldcp(jl)*zcond1
          pq(jl,kk)=pq(jl,kk)-zcond1
          !     ENDIF
        ENDIF
      ENDDO
      
    ENDIF
    
    IF(kcall == 0) THEN
      
!DIR$ IVDEP
!OCL NOVREC
      DO jl=kidia,kfdia
        zqp    =1.0_JPRB/psp(jl)
        ztarg    =pt(jl,kk)
        zfoeew    =r2es*EXP(z3es(jl)*(ztarg    -rtt)/(ztarg    -z4es(jl)))
        zqsat    =zqp    *zfoeew !(JL)
        IF (zqsat     > zqmax) THEN
          zqsat    =zqmax
        ENDIF
        zcor    =1.0_JPRB/(1.0_JPRB-retv  *zqsat    )
        zqsat    =zqsat    *zcor
        z2s    =z5alcp(jl)/(ztarg    -z4es(jl))**2
        zcond1    =(pq(jl,kk)-zqsat    )/(1.0_JPRB+zqsat    *zcor    *z2s    )
        pt(jl,kk)=pt(jl,kk)+zaldcp(jl)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
        ztarg    =pt(jl,kk)
        zfoeew    =r2es*EXP(z3es(jl)*(ztarg    -rtt)/(ztarg    -z4es(jl)))
        zqsat    =zqp    *zfoeew
        IF (zqsat     > zqmax) THEN
          zqsat    =zqmax
        ENDIF
        zcor    =1.0_JPRB/(1.0_JPRB-retv  *zqsat    )
        zqsat    =zqsat    *zcor
        z2s    =z5alcp(jl)/(ztarg    -z4es(jl))**2
        zcond1    =(pq(jl,kk)-zqsat    )/(1.0_JPRB+zqsat    *zcor    *z2s    )
        pt(jl,kk)=pt(jl,kk)+zaldcp(jl)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
      ENDDO
      
    ENDIF
    
    IF(kcall == 4) THEN
      
!DIR$ IVDEP
!OCL NOVREC
      DO jl=kidia,kfdia
        zqp    =1.0_JPRB/psp(jl)
        ztarg    =pt(jl,kk)
        zfoeew    =r2es*EXP(z3es(jl)*(ztarg    -rtt)/(ztarg    -z4es(jl)))
        zqsat    =zqp    *zfoeew!(JL)
        IF (zqsat     > zqmax) THEN
          zqsat    =zqmax
        ENDIF
        zcor    =1.0_JPRB/(1.0_JPRB-retv  *zqsat    )
        zqsat    =zqsat    *zcor
        z2s    =z5alcp(jl)/(ztarg    -z4es(jl))**2
        zcond    =(pq(jl,kk)-zqsat    )/(1.0_JPRB+zqsat    *zcor    *z2s    )
        pt(jl,kk)=pt(jl,kk)+zaldcp(jl)*zcond
        pq(jl,kk)=pq(jl,kk)-zcond
        ztarg    =pt(jl,kk)
        zqsat    =zqp    *zfoeew
        IF (zqsat     > zqmax) THEN
          zqsat    =zqmax
        ENDIF
        zcor    =1.0_JPRB/(1.0_JPRB-retv  *zqsat    )
        zqsat    =zqsat    *zcor
        z2s    =z5alcp(jl)/(ztarg    -z4es(jl))**2
        zcond1    =(pq(jl,kk)-zqsat    )/(1.0_JPRB+zqsat    *zcor    *z2s    )
        pt(jl,kk)=pt(jl,kk)+zaldcp(jl)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
      ENDDO
      
    ENDIF
    
    IF (lhook) CALL dr_hook('CUADJTQS',1,zhook_handle)
  END SUBROUTINE cuadjtqs
  
  
END MODULE mo_adjust

