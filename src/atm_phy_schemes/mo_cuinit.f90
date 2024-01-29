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

MODULE mo_cuinit

  USE mo_kind   ,ONLY: jprb=>wp     , &
    &                  jpim=>i4

!  USE yomhook   ,ONLY: lhook,   dr_hook
  
  USE mo_cuparameters,  ONLY : rkap,  r4les, r4ies   ,&
    &                          r5les ,r5ies ,ralfdcp ,&
    &                          lphylin               ,&
    &                          lmfdudv               ,&
    &                          rcpd   ,retv, rd, rg  ,&
    &                          rlmin                 ,&
    &                          lhook,   dr_hook

  USE mo_adjust ,ONLY: cuadjtq ,cuadjtqs

  !KF new - use module instead of include file
  USE mo_cufunctions, ONLY: foealfcu, foeewm, foealfa

#ifdef _OPENACC
  USE mo_cuparameters,  ONLY : r2es, r3les, r3ies, rtt, r5alvcp, r5alscp
  USE mo_cufunctions,   ONLY : foeldcpmcu
#endif

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: cuinin, cubasen

CONTAINS

  SUBROUTINE cuinin &
    & ( kidia,    kfdia,    klon,    ktdia,  klev, njkt2, &
    & pten,     pqen,     pqsen,    puen,     pven,&
    & pvervel,  pgeo,     paph,&
    & klwmin,   klab,&
    & ptenh,    pqenh,    pqsenh,   pgeoh,&
    & ptu,      pqu,      ptd,      pqd,&
    & puu,      pvu,      pud,      pvd,&
    & plu,      lacc  )

    !>
    !! Description:
    !!          M.TIEDTKE         E.C.M.W.F.     12/89

    !!          PURPOSE
    !!          -------

    !!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
    !!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
    !!          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
    !!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS

    !!          METHOD.
    !!          --------
    !!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)


    !!     PARAMETER     DESCRIPTION                                   UNITS
    !!     ---------     -----------                                   -----
    !!     INPUT PARAMETERS (INTEGER):

    !!    *KIDIA*        START POINT
    !!    *KFDIA*        END POINT
    !!    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!    *KLEV*         NUMBER OF LEVELS

    !!    INPUT PARAMETERS (REAL):

    !!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
    !!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
    !!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
    !!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
    !!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
    !!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
    !!    *PGEO*         GEOPOTENTIAL                                  M2/S2
    !!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
    !!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
    !!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA

    !!   OUTPUT PARAMETERS (INTEGER):

    !!    *KLWMIN*       LEVEL OF MAXIMUM VERTICAL VELOCITY
    !!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
    !!                        KLAB=2 FOR CONDENSATION LEVEL

    !!    OUTPUT PARAMETERS (REAL):

    !!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS         K
    !!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS    KG/KG
    !!    *PQSENH*       ENV. SPEC. SATURATION HUMIDITY (T+1)
    !!                   ON HALF LEVELS                              KG/KG
    !!    *PTU*          TEMPERATURE IN UPDRAFTS                       K
    !!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                  KG/KG
    !!    *PTD*          TEMPERATURE IN DOWNDRAFTS                     K
    !!    *PQU*          SPEC. HUMIDITY IN DOWNDRAFTS                KG/KG
    !!    *PUU*          U-VELOCITY IN UPDRAFTS                       M/S
    !!    *PVU*          V-VELOCITY IN UPDRAFTS                       M/S
    !!    *PUD*          U-VELOCITY IN DOWNDRAFTS                     M/S
    !!    *PVD*          V-VELOCITY IN DOWNDRAFTS                     M/S
    ! !   *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG

    !!          EXTERNALS
    !!          ---------
    !!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS

    !!          MODIFICATIONS
    !!          -------------
    !!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
    !!        M.Hamrud      01-Oct-2003 CY28 Cleaning
    !!             05-02-11 : Optimisation (NJKT2) P. BECHTOLD

    !!----------------------------------------------------------------------

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    !USE YOMCST   , ONLY : RCPD
    !USE YOECUMF  , ONLY : NJKT2
    !USE YOEPHLI  , ONLY : LPHYLIN

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktdia
    INTEGER(KIND=jpim),INTENT(in)    :: njkt2
    REAL(KIND=jprb)   ,INTENT(in)    :: pten(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqsen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: puen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pven(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pvervel(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    INTEGER(KIND=jpim),INTENT(inout)   :: klwmin(klon)
    INTEGER(KIND=jpim),INTENT(inout)   :: klab(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pqenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pqsenh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(inout)   :: ptu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pqu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: ptd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pqd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: puu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pvu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pud(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pvd(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: plu(klon,klev)
    LOGICAL           ,INTENT(in)      :: lacc

    REAL(KIND=jprb) ::     zwmax(klon)
    REAL(KIND=jprb) ::     zph(klon)
    LOGICAL ::  llflag(klon)

    INTEGER(KIND=jpim) :: icall, ik, jk, jl

    REAL(KIND=jprb) :: zalfa, zzs, zorcpd
    REAL(KIND=jprb) :: zhook_handle

    !#include "cuadjtq.intfb.h"
    !#include "cuadjtqs.intfb.h"

    !----------------------------------------------------------------------

    !*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
    !*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
    !*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
    !                  ----------------------------------------------

    IF (lhook) CALL dr_hook('CUININ',0,zhook_handle)

    !$ACC DATA &
    !$ACC   PRESENT(pten, pqen, pqsen, puen, pven, pvervel, pgeo, paph, klwmin) &
    !$ACC   PRESENT(klab, ptenh, pqenh, pqsenh, pgeoh, ptu, pqu, ptd, pqd, puu) &
    !$ACC   PRESENT(pvu, pud, pvd, plu) &

    !$ACC   CREATE(zwmax, zph, llflag) &
    !$ACC   IF(lacc)

    !is set below:  zph  (:) = 0.0_JPRB
   
!   zalfa=LOG(2.0_JPRB)
    zorcpd=1._jprb/rcpd

    !$ACC WAIT(1) ! ACCWA (cray 15.0.1): The WAIT(1) and missing ASYNC(1) are required to avoid the error
                  !                      "Invalid accelerator routine parallelism attribute: cuadjtq".
                  !                      Note that cuadjtq is also called in other regions of the code,
                  !                      but only the call here in cuinin raises the error above.
    !$ACC PARALLEL DEFAULT(PRESENT) IF(lacc)

    !$ACC LOOP SEQ
    DO jk=ktdia+1,klev

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        ptenh(jl,jk)=(MAX(rcpd*pten(jl,jk-1)+pgeo(jl,jk-1),&
          & rcpd*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk))*zorcpd
        pqenh(jl,jk)=pqen(jl,jk-1)
        pqsenh(jl,jk)=pqsen(jl,jk-1)
        zph(jl)=paph(jl,jk)
        llflag(jl)=.TRUE.
      ENDDO

      !orig   IF(JK.GE.KLEV-1) GO TO 130
      IF(jk >= klev-1 .OR. jk<njkt2) CYCLE
      ik=jk

      IF(lphylin)THEN
        icall=0
#ifndef _OPENACC
        CALL cuadjtqs &
          & ( kidia,    kfdia,    klon,  klev,&
          & ik,&
          & zph,      ptenh,    pqsenh,   llflag,   icall)
#endif
      ELSE
        icall=3
        CALL cuadjtq &
          & ( kidia,    kfdia,    klon,    klev,&
          & ik,&
          & zph,      ptenh,    pqsenh,   llflag,   icall)
      ENDIF

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        pqenh(jl,jk)=MIN(pqen(jl,jk-1),pqsen(jl,jk-1))&
          & +(pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk)=MAX(pqenh(jl,jk),0.0_JPRB)
      ENDDO

      !orig  130   continue
    ENDDO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      !KF next two statements are commented in cosmo-tiedtke to avoid drizzling 2010-03-10
      ptenh(jl,klev)=(rcpd*pten(jl,klev)+pgeo(jl,klev)-pgeoh(jl,klev))*zorcpd
      pqenh(jl,klev)=pqen(jl,klev)
      ptenh(jl,1)=pten(jl,1)
      pqenh(jl,1)=pqen(jl,1)
      klwmin(jl)=klev
      zwmax(jl)=0.0_JPRB
    ENDDO

    !$ACC LOOP SEQ
    DO jk=klev-1,ktdia+1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zzs)
      DO jl=kidia,kfdia
        zzs=MAX(rcpd*ptenh(jl,jk)+pgeoh(jl,jk),&
          & rcpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
        ptenh(jl,jk)=(zzs-pgeoh(jl,jk))*zorcpd
      ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO jk=klev,ktdia+2,-1
!DIR$ IVDEP
!OCL NOVREC
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        IF(pvervel(jl,jk) < zwmax(jl)) THEN
          zwmax(jl)=pvervel(jl,jk)
          klwmin(jl)=jk
        ENDIF
      ENDDO
    ENDDO

    !-----------------------------------------------------------------------

    !*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
    !*                 ---------------------------------------------

    !$ACC LOOP SEQ
    DO jk=ktdia,klev
      ik=jk-1
      IF(jk == ktdia) ik=ktdia
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        ptu(jl,jk)=ptenh(jl,jk)
        ptd(jl,jk)=ptenh(jl,jk)
        pqu(jl,jk)=pqenh(jl,jk)
        pqd(jl,jk)=pqenh(jl,jk)
        plu(jl,jk)=0.0_JPRB
        puu(jl,jk)=puen(jl,ik)
        pud(jl,jk)=puen(jl,ik)
        pvu(jl,jk)=pven(jl,ik)
        pvd(jl,jk)=pven(jl,ik)
        klab(jl,jk)=0
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC END DATA

    IF (lhook) CALL dr_hook('CUININ',1,zhook_handle)
  END SUBROUTINE cuinin


SUBROUTINE cubasen &
 & ( kidia,    kfdia,  klon,  ktdia, klev, njkt1, njkt2,  &
 & entrorg, entstpc1, entstpc2, rdepths, texc, qexc,      &
 & lgrz_deepconv, mtnmask, ldland, ldlake,                &
 & ptenh,  pqenh, pgeoh, paph,  pqhfl, pahfs,             &
!& PSSTRU,   PSSTRV,                                      &
 & pten,     pqen,     pqsen, pgeo,                       &
 & puen,     pven,                                        &
 & ptu,      pqu,      plu,      puu,      pvu ,  pwubase,&
 & klab,     ldcum,    kcbot,                             &
!& LDSC,     KBOTSC,
 & kctop,    kdpl,     pcape, pvervel650, lvv_shallow_deep, lacc )

!>
!! Description:
!!          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
!!          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT

!!          A. Pier Siebesma   KNMI ********
!!          modified C Jakob (ECMWF) (01/2001)
!!          modified P Bechtold (ECMWF) (08/2002)
!!         (include cycling over levels to find unstable departure/base level+
!!           mixed layer properties +w Trigger)

!!          PURPOSE.
!!          --------
!!          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR CU-PARAMETRIZATION

!!          INTERFACE
!!          ---------
!!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!!                 KLAB=0 FOR STABLE LAYERS
!!                 KLAB=1 FOR SUBCLOUD LEVELS
!!                 KLAB=2 FOR CLOUD LEVELS LEVEL

!!          METHOD.
!!          --------
!!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP
!!          (ENTRAINING PLUME, WITH ENTRAINMENT PROPORTIONAL TO (1/Z))

!! History:
!!          A. Pier Siebesma   KNMI ********
!!          modified C Jakob (ECMWF) (01/2001)
!!          modified P Bechtold (ECMWF) (08/2002)
!!          (include cycling over levels to find unstable departure/base level+
!!           mixed layer properties +w Trigger)
!!
!!          MODIFICATIONS
!!!          -------------
!!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!!             02-11-02 : Use fixed last possible departure level and
!!                        last updraft computation level for bit-reproducibility
!!                                            D.Salmond &  J. Hague
!!            03-07-03 : Tuning for p690     J. Hague
!!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!!
!!     PARAMETER     DESCRIPTION                                   UNITS
!!     ---------     -----------                                   -----
!!     INPUT PARAMETERS (INTEGER):

!!    *KIDIA*        START POINT
!!    *KFDIA*        END POINT
!!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!!    *KLEV*         NUMBER OF LEVELS

!!    INPUT PARAMETERS (REAL):

!! not used at the moment because we want to use linear intepolation
!! for fields on the half levels.

!!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG

!!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!!    *PSSTRU*       KINEMATIC surface U-MOMENTUM FLUX             (M/S)^2
!!    *PSSTRV*       KINEMATIC surface V-MOMENTUM FLUX             (M/S)^2
!!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!!    *zdph*         pressure thickness on full levels               PA
!!    *zdgeoh*       geopotential thickness on full levels          M2/S2
!!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!!    *PQSEN*        PROVISIONAL ENVIRONMENT SATU. HUMIDITY (T+1)  KG/KG
!!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!!    *PVERVEL650*   VERTICAL VELOCITY AT 650hPa                   hPa/s  

!!    UPDATED PARAMETERS (REAL):

!!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S

!!    UPDATED PARAMETERS (INTEGER):

!!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!!                        KLAB=2 FOR CLOUD LEVELS

!!    OUTPUT PARAMETERS (LOGICAL):

!!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

!!    OUTPUT PARAMETERS (INTEGER):

!!    *KCBOT*       CLOUD BASE LEVEL !
!!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL
!!                  WITH A NON-ZERO CLOUD UPDRAFT.
!!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS
!!    *KDPL*        DEPARTURE LEVEL
!!    *PCAPE*       PSEUDOADIABATIQUE max CAPE (J/KG)

!!    LOGICAL PARAMETERS
!!    *LVV_SHALLOW_DEEP* USE VERTICAL VELOCITY AT 650hPa TO DISTINGUISH
!!                       BETWEEN SHALLOW AND DEEP CONVECTION
!!          EXTERNALS
!!          ---------
!!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT

!!          MODIFICATIONS
!!          -------------
!!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!!             02-11-02 : Use fixed last possible departure level and
!!                        last updraft computation level for bit-reproducibility
!!                                            D.Salmond &  J. Hague
!!            03-07-03 : Tuning for p690     J. Hague
!!       M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------    

!
!USE parkind1  ,ONLY : jpim     ,jprb
!USE yomhook   ,ONLY : lhook,   dr_hook
!
!USE yomcst   , ONLY : rcpd     ,retv, rd, rg,&
! & rlvtt    ,rlstt    ,rtt
!USE yoevdf   , ONLY : rkap
!USE yoecumf  , ONLY : lmfdudv, ENTRORG, rdepths,  ENTSTPC1, ENTSTPC2, njkt1, njkt2
!USE yoecldp  , ONLY : rlmin
!USE yoethf   , ONLY : r2es     ,r3les    ,r3ies    ,r4les    ,&
! & r4ies    ,r5les    ,r5ies    ,r5alvcp  ,r5alscp  ,&
! & ralvdcp  ,ralsdcp  ,ralfdcp  ,rtwat    ,rtice    ,rticecu  ,&
! & rtwat_rticecu_r    ,rtwat_rtice_r
!


INTEGER(KIND=jpim),INTENT(in)    :: klon
INTEGER(KIND=jpim),INTENT(in)    :: klev
INTEGER(KIND=jpim),INTENT(in)    :: kidia
INTEGER(KIND=jpim),INTENT(in)    :: kfdia
INTEGER(KIND=jpim),INTENT(in)    :: ktdia
INTEGER(KIND=jpim),INTENT(in)    :: njkt1, njkt2
REAL(KIND=jprb)   ,INTENT(in)    :: entrorg, entstpc1, entstpc2
REAL(KIND=jprb)   ,INTENT(in)    :: rdepths
REAL(KIND=jprb)   ,INTENT(in)    :: texc, qexc
REAL(KIND=jprb)   ,INTENT(in)    :: mtnmask(klon)
LOGICAL           ,INTENT(in)    :: ldland(klon)
LOGICAL           ,INTENT(in)    :: ldlake(klon)
LOGICAL           ,INTENT(in)    :: lgrz_deepconv
REAL(KIND=jprb)   ,INTENT(in)    :: ptenh(klon,klev)
REAL(KIND=jprb)   ,INTENT(in)    :: pqenh(klon,klev)
REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1)
REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
REAL(KIND=jprb)   ,INTENT(in)    :: pqhfl(klon,klev+1)
REAL(KIND=jprb)   ,INTENT(in)    :: pahfs(klon,klev+1)
REAL(KIND=jprb)   ,INTENT(in)    :: pten(klon,klev)
REAL(KIND=jprb)   ,INTENT(in)    :: pqen(klon,klev)
REAL(KIND=JPRB)   ,INTENT(in)    :: pqsen(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: puen(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pven(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: ptu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: pqu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: plu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: puu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: pvu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: pwubase(klon) 
INTEGER(KIND=jpim),INTENT(inout) :: klab(klon,klev) 
LOGICAL           ,INTENT(inout) :: ldcum(klon) 
INTEGER(KIND=jpim),INTENT(inout) :: kcbot(klon) 
INTEGER(KIND=jpim),INTENT(out)   :: kctop(klon) 
INTEGER(KIND=jpim),INTENT(out)   :: kdpl(klon) 
REAL(KIND=jprb)   ,INTENT(out)   :: pcape(klon) 

REAL(KIND=jprb)   ,INTENT(in)    :: pvervel650(:) 
LOGICAL           ,INTENT(in)    :: lvv_shallow_deep
LOGICAL           ,INTENT(in)    :: lacc

! local variables

LOGICAL                          :: ldsc(klon) 
INTEGER(KIND=jpim)               :: kbotsc(klon) 

INTEGER(KIND=jpim) ::  ictop(klon),            icbot(klon),&
 & ibotsc(klon),           ilab(klon,klev),&
 & idpl(klon)  

LOGICAL ::         ll_ldbase(klon),&
 & llgo_on(klon),&
 & lldeep(klon),    lldcum(klon), &
 & lldsc(klon),     llfirst(klon)  
LOGICAL ::     llreset,        llresetjl(klon), ldocean(klon)

INTEGER(KIND=jpim) :: icall, ik, is, jk, jl, jkk, jkt1, jkt2, jkt, jkb ! ,IKB

REAL(KIND=jprb)    :: &
 & zsenh(klon,klev+1),&
 & zqenh(klon,klev+1),&
 & zsuh (klon,klev),&
 & zwu2h(klon,klev),&
 & zbuoh(klon,klev)  
REAL(KIND=jprb) :: zqold(klon),zph(klon)
REAL(KIND=jprb) :: zmix(klon)
REAL(KIND=jprb) :: zdz(klon)

REAL(KIND=jprb) ::    zlu(klon,klev),   zqu(klon,klev),&
 & ztu(klon,klev), &
 & zuu(klon,klev),   zvu(klon,klev)  

REAL(KIND=jprb) :: zcape(klon,klev) ! local for CAPE at every departure level

REAL(KIND=jprb) :: zbuof     ! BUOYANCY
REAL(KIND=jprb) :: zrho      ! DENSITY AT SURFACE (KG/M^3) 
REAL(KIND=jprb) :: zkhvfl    ! SURFACE BUOYANCY FLUX (K M/S)
REAL(KIND=jprb) :: zws       ! SIGMA_W AT LOWEST MODEL HALFLEVEL (M/S)
REAL(KIND=jprb) :: zqexc     ! HUMIDITY EXCESS AT LOWEST MODEL HALFLEVEL (KG/KG)
REAL(KIND=jprb) :: ztexc     ! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL (K)
REAL(KIND=jprb) :: ztex(klon), zqex(klon) ! Corresponding fields at lowest model level
REAL(KIND=jprb) :: zeps      ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
REAL(KIND=jprb) :: ztvenh    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
REAL(KIND=jprb) :: ztvuh     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
REAL(KIND=jprb) :: zlglac    ! UPDRAFT LIQUID WATER FROZEN IN ONE LAYER
REAL(KIND=jprb) :: zqsu, zcor, zdq, zalfaw, zfacw, zfaci, zfac,&
 & zesdp, zdqsdt, zdtdp, zdp,zpdifftop,zpdiffbot,zsf,zqf,zaw,zbw
!REAL(KIND=jprb) :: ztven1(klon), ztven2, ztvu1(klon), ztvu2 ! pseudoadiabatique T_v
REAL(KIND=jprb) :: ztven1(klon,klev),ztven2(klon,klev),ztvu1(klon,klev),ztvu2(klon,klev)
REAL(KIND=jprb) :: zdtvtrig(klon) ! virtual temperatures
REAL(KIND=jprb) :: zwork1, zwork2! work arrays for T and w perturbations
REAL(KIND=jprb) :: zrcpd, zrg, ztmp, zredfac, zcbase
REAL(KIND=jprb) :: zhook_handle

#ifdef _OPENACC
REAL(KIND=jprb) :: uzqp,uzl,uzi,uzqsat,uzcor,uzf,uzcond,uzcond1
#endif

!#include "cuadjtq.intfb.h"
!#include "fcttre.func.h"

!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!                  -------------------------------
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CUBASEN',0,zhook_handle)

!$ACC DATA &
!$ACC   PRESENT(mtnmask, ldland, ldlake, ptenh, pqenh, pgeoh, paph, pqhfl, pahfs) &
!$ACC   PRESENT(pten, pqen, pqsen, pgeo, puen, pven, ptu, pqu, plu, puu, pvu) &
!$ACC   PRESENT(pwubase, klab, ldcum, kcbot, kctop, kdpl, pcape, pvervel650) &

!$ACC   CREATE(ldsc, kbotsc, ictop, icbot, ibotsc, ilab, idpl, ll_ldbase, llgo_on) &
!$ACC   CREATE(lldeep, lldcum, lldsc, llfirst, llresetjl, ldocean, zsenh, zqenh) &
!$ACC   CREATE(zsuh, zwu2h, zbuoh, zqold, zph, zmix, zdz, zlu, zqu, ztu) &
!$ACC   CREATE(zuu, zvu, zcape, ztex, zqex, ztven1, ztven2, ztvu1, ztvu2, zdtvtrig) &
!$ACC   IF(lacc)

zaw    = 1.0_JPRB
zbw    = 1.0_JPRB
jkt1=njkt1
jkt2=njkt2
zrg=1.0_JPRB/rg
zrcpd=1.0_JPRB/rcpd

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

!$ACC LOOP SEQ
DO jk=1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    ztven1(jl,jk)  = 0._jprb
    ztven2(jl,jk)  = 0._jprb
    ztvu1 (jl,jk)  = 0._jprb
    ztvu2 (jl,jk)  = 0._jprb
    zsenh (jl,jk)  = 0._jprb
    zqenh (jl,jk)  = 0._jprb
    zsuh  (jl,jk)  = 0._jprb
    zwu2h (jl,jk)  = 0._jprb
    zbuoh (jl,jk)  = 0._jprb
  ENDDO
ENDDO

!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  zqold  (jl)         = 0._jprb
  zph    (jl)         = 0._jprb
  zdz    (jl)         = 0._jprb
  zsenh  (jl,klev+1)  = 0._jprb
  zqenh  (jl,klev+1)  = 0._jprb

  pwubase(jl) = 0.0_JPRB
  llgo_on(jl) = .TRUE.
  llfirst(jl) = .TRUE.
  ldocean(jl) = .NOT. (ldland(jl) .OR. ldlake(jl) .OR. lgrz_deepconv)
  kdpl   (jl) = klev
ENDDO

!$ACC LOOP SEQ
DO jk=ktdia,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    ztu(jl,jk) = ptu(jl,jk)
    zqu(jl,jk) = pqu(jl,jk)
    zlu(jl,jk) = plu(jl,jk)
    zuu(jl,jk) = puu(jl,jk)
    zvu(jl,jk) = pvu(jl,jk)
    ilab(jl,jk)= klab(jl,jk)
    zcape(jl,jk)= 0.0_JPRB
  ENDDO
ENDDO

!----------------------------------------------------------------------
!       -----------------------------------------------------------
!       1.1  PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!             OF SPECIFIC HUMIDITY AND STATIC ENERGY
!       -----------------------------------------------------------

!$ACC LOOP SEQ
DO jk=ktdia,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    ZWU2H(JL,JK)=0.0_JPRB
    zqenh(jl,jk) = pqenh(jl,jk)
    zsenh(jl,jk) = rcpd*ptenh(jl,jk)+pgeoh(jl,jk)
  ENDDO
ENDDO

!$ACC LOOP SEQ
DO jkk=klev,MAX(ktdia,jkt1),-1 ! Big external loop for level testing:
                               ! find first departure level that produces deepest cloud top
                               ! or take surface level for shallow convection and Sc
   !
   !        ---------------------------------------------------------
   !        1.2    INITIALISE FIELDS AT DEPARTURE HALF MODEL LEVEL
   !        ---------------------------------------------------------
   !

#ifndef _OPENACC
  is=0
#endif

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    IF (llgo_on(jl)) THEN
#ifndef _OPENACC
      is=is+1
#endif
      idpl(jl)    =jkk      ! departure level
      icbot  (jl) =jkk      ! cloud base level for convection, (-1 if not found)
      ibotsc (jl) =klev-1   ! sc    base level for sc-clouds , (-1 if not found)
      ictop(jl)   =klev-1   ! cloud top for convection (-1 if not found)
      lldcum(jl)  =.FALSE.  ! on exit: true if cloudbase=found
      lldsc (jl)  =.FALSE.  ! on exit: true if cloudbase=found
      ll_ldbase(jl)   =.FALSE. ! on exit: true if cloudbase=found
      lldeep(jl)  =.FALSE.  ! true if convection is deep
      zdtvtrig(jl) =0.0_JPRB
      zuu(jl,jkk) =puen(jl,jkk)*(paph(jl,jkk+1)-paph(jl,jkk))
      zvu(jl,jkk) =pven(jl,jkk)*(paph(jl,jkk+1)-paph(jl,jkk))
    ENDIF 
  ENDDO

#ifndef _OPENACC
  IF(is /= 0) THEN
#endif

    IF(jkk == klev) THEN

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrho, zkhvfl, zws, zredfac, ztvenh, ztvuh)
      DO jl=kidia,kfdia
        IF (llgo_on(jl)) THEN
          zrho  = paph(jl,jkk+1)/(rd*(pten(jl,jkk)*(1.0_JPRB+retv*pqen(jl,jkk))))
          zkhvfl= (pahfs(jl,jkk+1)*zrcpd+retv*pten(jl,jkk)*pqhfl(jl,jkk+1))/zrho
!          ZUST  = MAX(SQRT(PSSTRU(JL)**2 + PSSTRV(JL)**2),ZREPUST)     !u* (repust=10e-4)
!          ZWS=ZUST**3._JPRB- 1.5_JPRB*RKAP*ZKHVFL*PGEOH(JL,KLEV)/PTEN(JL,KLEV)
          zws=0.001_JPRB - 1.5_JPRB*rkap*zkhvfl*(pgeoh(jl,klev)-pgeoh(jl,klev+1))/pten(jl,klev)
          ztex(jl)= 0.0_JPRB
          zqex(jl)= 0.0_JPRB
          IF( zkhvfl < 0.0_JPRB ) THEN
            zws=1.2_JPRB*zws**.3333_JPRB
            ilab(jl,jkk)= 1
            zredfac = 1._jprb/(1._jprb+mtnmask(jl))
            ztex(jl)     = MAX(-1.5_JPRB*pahfs(jl,jkk+1)/(zrho*zws*rcpd),0.5_JPRB*texc)*zredfac
            zqex(jl)     = MAX(-1.5_JPRB*pqhfl(jl,jkk+1)/(zrho*zws),0.0_JPRB)*zredfac
            zqu (jl,jkk) = zqenh(jl,jkk) + zqex(jl)
            zsuh (jl,jkk) = zsenh(jl,jkk) + rcpd*ztex(jl)
            ztu (jl,jkk) = (zsenh(jl,jkk)-pgeoh(jl,jkk))*zrcpd + ztex(jl)
            zlu (jl,jkk) = 0.0_JPRB
            zwu2h(jl,jkk) = zws**2
        !
        !  determine buoyancy at lowest half level
        !
            ztvenh            = (1.0_JPRB+retv*zqenh(jl,jkk)) &
             & *(zsenh(jl,jkk)-pgeoh(jl,jkk))*zrcpd  
            ztvuh             = (1.0_JPRB+retv*zqu(jl,jkk))*ztu(jl,jkk)
            zbuoh(jl,jkk) = (ztvuh-ztvenh)*rg/ztvenh
          ELSE
            llgo_on(jl)=.FALSE.      ! non-convective point
          ENDIF
        ENDIF
      ENDDO

    ELSE

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zredfac, ztexc, zqexc, zwork1, jk, zwork2, ztvenh, ztvuh)
      DO jl=kidia,kfdia
        IF (llgo_on(jl)) THEN
          ilab(jl,jkk)= 1
          zredfac = 1._jprb/(1._jprb+mtnmask(jl))
          ztexc=texc*zredfac
          zqexc=qexc*pqenh(jl,jkk)*zredfac
          IF (jkk == klev-1 .AND. ldocean(jl) ) THEN
            ztexc = MAX(ztexc, ztex(jl))
            ztexc = MIN(ztexc, 3.0_JPRB)
            zqexc = MAX(zqexc, zqex(jl))
            zqexc = MIN(zqexc, 2.E-3_JPRB)
          ENDIF
          zqu (jl,jkk) = zqenh(jl,jkk) + zqexc
          zsuh (jl,jkk) = zsenh(jl,jkk) + rcpd*ztexc
          ztu (jl,jkk) = (zsenh(jl,jkk)-pgeoh(jl,jkk))*zrcpd + ztexc
          zlu (jl,jkk) = 0.0_JPRB
         ! construct mixed layer for parcels emanating in lowest 60 hPa
          IF (paph(jl,klev+1)-paph(jl,jkk-1)<60.e2_jprb) THEN
            zqu(jl,jkk) =0.0_JPRB
            zsuh(jl,jkk)=0.0_JPRB
            zwork1      =0.0_JPRB
            !original code does not vectorise!
            !DO jk=jkk+1,jkk-1,-1
            !  IF( zwork1 < 50.e2_jprb ) THEN
            !    zwork2=paph(jl,jk)-paph(jl,jk-1)
            !    zwork1      =zwork1+zwork2
            !    zqu(jl,jkk) =zqu(jl,jkk) +zqenh(jl,jk)*zwork2
            !    zsuh(jl,jkk)=zsuh(jl,jkk)+zsenh(jl,jk)*zwork2
            !  ENDIF
            !ENDDO
            !KF expand per hand
            jk = jkk+1
            !>KF
             ZWORK2=PAPH(JL,JK)-PAPH(JL,JK-1)
            !zwork2=zdph(jl,jk-1)
            zwork1      =zwork1+zwork2
            zqu(jl,jkk) =zqu(jl,jkk) +zqenh(jl,jk)*zwork2
            zsuh(jl,jkk)=zsuh(jl,jkk)+zsenh(jl,jk)*zwork2
            jk = jkk
            IF( zwork1 < 50.e2_jprb ) THEN
               ZWORK2=PAPH(JL,JK)-PAPH(JL,JK-1)
              !zwork2=zdph(jl,jk-1)
              zwork1      =zwork1+zwork2
              zqu(jl,jkk) =zqu(jl,jkk) +zqenh(jl,jk)*zwork2
              zsuh(jl,jkk)=zsuh(jl,jkk)+zsenh(jl,jk)*zwork2
            ENDIF
            jk = jkk-1
            IF( zwork1 < 50.e2_jprb ) THEN
               ZWORK2=PAPH(JL,JK)-PAPH(JL,JK-1)
              !zwork2=zdph(jl,jk-1)
               zwork1      =zwork1+zwork2
               zqu(jl,jkk) =zqu(jl,jkk) +zqenh(jl,jk)*zwork2
               zsuh(jl,jkk)=zsuh(jl,jkk)+zsenh(jl,jk)*zwork2
             ENDIF
             !<KF
            zqu(jl,jkk) =zqu(jl,jkk) /zwork1+zqexc
            zsuh(jl,jkk)=zsuh(jl,jkk)/zwork1+rcpd*ztexc
            ztu(jl,jkk) =(zsuh(jl,jkk)-pgeoh(jl,jkk))*zrcpd+ztexc
          ENDIF
          zwu2h(jl,jkk) = 1.0_JPRB
      !
      !  determine buoyancy at lowest half level
      !
          ztvenh            = (1.0_JPRB+retv*zqenh(jl,jkk)) &
           & *(zsenh(jl,jkk)-pgeoh(jl,jkk))*zrcpd  
          ztvuh             = (1.0_JPRB+retv*zqu(jl,jkk))*ztu(jl,jkk)
          zbuoh(jl,jkk) = (ztvuh-ztvenh)*rg/ztvenh
        ENDIF
      ENDDO
   
    ENDIF  ! jkk=klev
#ifndef _OPENACC
  ENDIF    ! is/=0  or ANY(llgo_on)
#endif
   
   !----------------------------------------------------------------------
   !     2.0          DO ASCENT IN SUBCLOUD AND LAYER,
   !                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
   !                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
   !                  CHECK FOR BUOYANCY AND SET FLAGS
   !                  -------------------------------------
   !       ------------------------------------------------------------
   !        1.2  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
   !       ------------------------------------------------------------

  !$ACC LOOP SEQ
  smalldoloop: DO jk=jkk-1,MAX(ktdia,jkt2),-1
#ifndef _OPENACC
    is=0
#endif

    IF(jkk==klev) THEN ! 1/z mixing for shallow

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zeps, zqf, zsf, ztmp)
      DO jl=kidia,kfdia
        IF (llgo_on(jl)) THEN
#ifndef _OPENACC
          is         = is+1
#endif
          zdz(jl)    = (pgeoh(jl,jk) - pgeoh(jl,jk+1))*zrg
          zeps       = entstpc1/((pgeoh(jl,jk)-pgeoh(jl,klev+1))*zrg) + entstpc2
  !        zmix(jl)   = 0.5_JPRB*zdz(jl)*zeps
          zmix(jl)   = (1.3_JPRB - MIN(1.0_JPRB,PQEN(JL,JK)/PQSEN(JL,JK)))*zdz(jl)*zeps
          zqf = (pqenh(jl,jk+1) + pqenh(jl,jk))*0.5_JPRB
          zsf = (zsenh(jl,jk+1) + zsenh(jl,jk))*0.5_JPRB
          ztmp = 1.0_JPRB/(1.0_JPRB+zmix(jl))
          zqu(jl,jk)= (zqu(jl,jk+1)*(1.0_JPRB-zmix(jl))&
         & +2.0_JPRB*zmix(jl)*zqf) * ztmp  
          zsuh (jl,jk)= (zsuh(jl,jk+1)*(1.0_JPRB-zmix(jl))&
         & +2.0_JPRB*zmix(jl)*zsf) * ztmp  
          zqold(jl)  = zqu(jl,jk)
          ztu (jl,jk) = (zsuh(jl,jk)-pgeoh(jl,jk))*zrcpd
          zph  (jl)    = paph(jl,jk)
        ENDIF
      ENDDO

    ELSE

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqf, zsf)
      DO jl=kidia,kfdia
        IF (llgo_on(jl)) THEN
#ifndef _OPENACC
          is         = is+1
#endif
          zdz(jl)    = (pgeoh(jl,jk) - pgeoh(jl,jk+1))*zrg
          zqf = (pqenh(jl,jk+1) + pqenh(jl,jk))*0.5_JPRB
          zsf = (zsenh(jl,jk+1) + zsenh(jl,jk))*0.5_JPRB
!         zmix(jl)=2.0_JPRB*0.8E-4_JPRB*zdz(jl)*(paph(jl,jk)/paph(jl,klev+1))**3
!         ZMIX(JL)=0.4_JPRB*ENTRORG*ZDZ(JL)*MIN(1.0_JPRB,(PQSEN(JL,JK)/PQSEN(JL,KLEV))**3)
          ZMIX(JL)=MERGE(0.3_JPRB, 1.3_JPRB - MIN(1.0_JPRB,PQEN(JL,JK)/PQSEN(JL,JK)), ldocean(jl)) &
         &  * ENTRORG*ZDZ(JL)*MAX(0.2_JPRB,MIN(1.0_JPRB,(PQSEN(JL,JK)/PQSEN(JL,KLEV))**2))
          ! Limitation to avoid trouble at very coarse vertical resolution
          zmix(jl) = MIN(1.0_jprb,zmix(jl))
          zqu(jl,jk)= zqu(jl,jk+1)*(1.0_JPRB-zmix(jl))+ zqf*zmix(jl)
          zsuh(jl,jk)= zsuh(jl,jk+1)*(1.0_JPRB-zmix(jl))+ zsf*zmix(jl)
          zqold(jl)  = zqu(jl,jk)
          ztu (jl,jk)= (zsuh(jl,jk)-pgeoh(jl,jk))*zrcpd
          zph  (jl)  = paph(jl,jk)
        ENDIF
      ENDDO

    ENDIF

#ifndef _OPENACC
    IF (is == 0) EXIT smalldoloop
#endif
 
    ik=jk
    icall=1
     
#if defined(_OPENACC) && (__NVCOMPILER_MAJOR__ <= 21)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    do jl=kidia,kfdia
    if (jl==-1292) print *, 'working on:  ', ik, jk, jkk
    enddo
#endif

    CALL cuadjtq &
     & ( kidia,    kfdia,    klon,    klev,      ik,&
     &   zph,      ztu,      zqu,     llgo_on,   icall)

    !DIR$ IVDEP
    !OCL NOVREC
   
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdq, zlglac, ztvuh, ztvenh, zbuof, ztmp, ik) &
    !$ACC   PRIVATE(zqsu, zcor, zalfaw, zfacw, zfaci, zfac, zesdp, zdqsdt, zdtdp) &
    !$ACC   PRIVATE(zdp, zpdifftop, zpdiffbot, zcbase, jkb)
    DO jl=kidia,kfdia
      IF(llgo_on(jl)) THEN
   
   ! add condensation to water
   
        zdq=MAX(zqold(jl)-zqu(jl,jk),0.0_JPRB)
        zlu(jl,jk)=zlu(jl,jk+1)+zdq

   ! freezing
   
        zlglac=zdq*((1.0_JPRB-foealfcu(ztu(jl,jk)))-&
         & (1.0_JPRB-foealfcu(ztu(jl,jk+1))))  
              
   
   ! pseudo-microphysics
   
        IF(jkk==klev) THEN  ! no precip for shallow
          zlu(jl,jk)=MIN(zlu(jl,jk),5.e-3_JPRB)
   !* chose a more pseudo-adiabatic formulation as original overestimates
   !* water loading efect and therefore strongly underestimates cloud thickness
        ELSE 
          zlu(jl,jk)=0.5_JPRB*zlu(jl,jk) 
        ENDIF
   
   ! update dry static energy after condensation + freezing
   
        zsuh(jl,jk)    = rcpd*(ztu(jl,jk)+ralfdcp*zlglac)+pgeoh(jl,jk)
         
   ! Buoyancy on half and full levels
            
        ztvuh           = (1.0_JPRB+retv*zqu(jl,jk)-zlu(jl,jk))*ztu(jl,jk)&
         & +ralfdcp*zlglac  
        ztvenh          = (1.0_JPRB+retv*zqenh(jl,jk)) &
         & *(zsenh(jl,jk)-pgeoh(jl,jk))*zrcpd  
        zbuoh(jl,jk)   = (ztvuh-ztvenh)*rg/ztvenh
        zbuof          = (zbuoh(jl,jk) + zbuoh(jl,jk+1))*0.5_JPRB
   
   ! solve kinetic energy equation
   
        ztmp=1.0_JPRB/(1.0_JPRB+2.0_JPRB*zbw*zmix(jl))
        zwu2h(jl,jk) = (zwu2h(jl,jk+1)*(1.0_JPRB-2.0_JPRB*zbw*zmix(jl))&
         & +2.0_JPRB*zaw*zbuof*zdz(jl)) * ztmp  
   
   ! compute pseudoadiabatique CAPE for diagnostics
   
        ztvu2(jl,jk) = ztu(jl,jk)  *(1.0_JPRB+retv*zqu(jl,jk))
        ztven2(jl,jk)= ptenh(jl,jk)*(1.0_JPRB+retv*pqenh(jl,jk))
        IF (jk + 1 == jkk) THEN
           ztvu1(jl,jk)  = ztvu2(jl,jk)
           ztven1(jl,jk) = ztven2(jl,jk)
        ENDIF
        zbuof = (ztvu2(jl,jk)+ztvu1(jl,jk)-ztven1(jl,jk)-ztven2(jl,jk))/ztven2(jl,jk)
        zbuof = zbuof*zdz(jl)*rg
        zcape(jl,jkk)  = zcape(jl,jkk) + MAX(0.0_JPRB,zbuof)
        ztvu1(jl,jk)=ztvu2(jl,jk)
        ztven1(jl,jk)=ztven2(jl,jk)
   
   ! first layer with liquid water - find exact cloud base
   
        IF(zlu(jl,jk) >0.0_JPRB.AND.ilab(jl,jk+1)==1) THEN
           
          ik=jk+1
          zqsu=foeewm(ztu(jl,ik))/paph(jl,ik)
          zqsu=MIN(0.5_JPRB,zqsu)
          zcor=1.0_JPRB/(1.0_JPRB-retv  *zqsu)
          zqsu=zqsu*zcor
          zdq=MIN(0._jprb,zqu(jl,ik)-zqsu)
          zalfaw=foealfa(ztu(jl,ik))
          zfacw=r5les/((ztu(jl,ik)-r4les)**2)
          zfaci=r5ies/((ztu(jl,ik)-r4ies)**2)
          zfac=zalfaw*zfacw+(1.0_JPRB-zalfaw)*zfaci
          zesdp=foeewm(ztu(jl,ik))/paph(jl,ik)
          zcor=1.0_JPRB/(1.0_JPRB-retv*zesdp)
          zdqsdt=zfac*zcor*zqsu
          zdtdp=rd*ztu(jl,ik)/(rcpd*paph(jl,ik))
          zdp=zdq/(zdqsdt*zdtdp)
          zcbase=paph(jl,ik)+zdp
           
   ! chose nearest half level as cloud base
   
          zpdifftop=zcbase-paph(jl,jk)
          zpdiffbot=paph(jl,jk+1)-zcbase
           
          IF(zpdifftop > zpdiffbot.AND.zwu2h(jl,jk+1)>0.0_JPRB) THEN
            jkb=MIN(klev-1,jk+1)
            ilab(jl,jkb)=2 
            ilab(jl,jk)=2
            ll_ldbase(jl) =.TRUE.
            lldsc(jl)   =.TRUE.
            ibotsc(jl) =jkb
            icbot(jl)  =jkb
            zlu(jl,jk+1) = rlmin
          ELSEIF(zpdifftop <= zpdiffbot.AND.zwu2h(jl,jk)>0.0_JPRB) THEN
            ilab(jl,jk)=2
            ll_ldbase(jl) =.TRUE.
            lldsc(jl)   =.TRUE.
            ibotsc(jl) =jk
            icbot(jl)  =jk
          ENDIF
          jkb=icbot(jl)
   
        ENDIF
   
   ! decide on presence of convection, cloud base and cloud top based on
   ! kinetic energy
   
        IF (zwu2h(jl,jk) < 0.0_JPRB) THEN
          llgo_on(jl) = .FALSE.             
          IF (zlu(jl,jk+1)>0.0_JPRB) THEN
            ictop(jl)   = jk
            lldcum(jl)   = .TRUE.
          ELSE
            lldcum(jl)   = .FALSE.
          ENDIF
        ELSE
          IF (zlu(jl,jk)>0.0_JPRB) THEN
            ilab(jl,jk) = 2
          ELSE
            ilab(jl,jk) = 1
          ENDIF
        ENDIF
      ENDIF
    ENDDO
   
    IF(lmfdudv.AND.jkk==klev) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        IF(.NOT.ll_ldbase(jl).AND.llgo_on(jl)) THEN
          zuu(jl,jkk)=zuu(jl,jkk)+puen(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
          zvu(jl,jkk)=zvu(jl,jkk)+pven(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
        ENDIF
      ENDDO
    ENDIF
   
!     IF (IS == 0) EXIT
  ENDDO smalldoloop  ! loop over jk
   
  IF( jkk==klev) THEN
    ! set values for departure level for PBL clouds = first model level

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jkt, jkb)
    DO jl=kidia,kfdia
      ldsc(jl)  = lldsc(jl)
      IF(ldsc(jl)) THEN
        kbotsc(jl)= ibotsc(jl)
      ELSE
        kbotsc(jl)=-1
      ENDIF
    
      llgo_on(jl) = .FALSE.
      jkt=ictop(jl)
      jkb=icbot(jl)
      
      ! use alternative criterion for "deep convection"
      IF (lvv_shallow_deep) THEN
         lldeep(jl)= (pvervel650(jl) < 0.0_jprb)
      ELSE
         lldeep(jl)=paph(jl,jkb)-paph(jl,jkt)>rdepths
      ENDIF
      IF(lldeep(jl)) lldcum(jl)=.FALSE. ! no deep allowed for KLEV
      lldeep(jl)=.FALSE. ! for deep convection start only at level KLEV-1
                            ! and form mixed layer, so go on
      ! test further for deep convective columns as not yet found
      IF ( lldeep(jl) ) llfirst(jl)=.FALSE.
      llgo_on(jl) = .NOT.lldeep(jl)
      IF(lldcum(jl)) THEN ! only write convective properties if cloud is shallow, and failed deep test
        kcbot(jl)= icbot(jl)
        kctop(jl)= ictop(jl)
        kdpl(jl)  = idpl(jl)
        ldcum(jl) = lldcum(jl)
        pwubase(jl)=SQRT(MAX(zwu2h(jl,jkb),0.0_JPRB))
      ELSE
        kctop(jl)=-1
        kcbot(jl)=-1
        kdpl(jl) =klev-1
        ldcum(jl)=.FALSE.
        pwubase(jl)=0.0_JPRB
      ENDIF
    ENDDO

    !$ACC LOOP SEQ
    DO jk=klev,ktdia,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jkt)
      DO jl=kidia,kfdia
        jkt=ictop(jl)
        IF ( jk>=jkt ) THEN
          klab(jl,jk)=ilab(jl,jk)
          ptu(jl,jk)=ztu(jl,jk)
          pqu(jl,jk)=zqu(jl,jk)
          plu(jl,jk)=zlu(jl,jk)
        ENDIF
      ENDDO
    ENDDO
  ENDIF !jkk==klev

  ! Modification by GZ to close vectorization gap
  IF( jkk < klev ) THEN
    !llreset=.FALSE.

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      IF ( .NOT.lldeep(jl) ) THEN
      !  jkt=ictop(jl)
      !  jkb=icbot(jl)
      ! test on cloud thickness and buoyancy
      !use alternative criterion for "deep convection"
         IF (lvv_shallow_deep) THEN
            lldeep(jl)= (pvervel650(jl) < 0.0_jprb)
         ELSE
            lldeep(jl)=paph(jl,icbot(jl))-paph(jl,ictop(jl))>=rdepths
         ENDIF
        !LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>=RDEPTHS &
        !   &.AND. ZDTVTRIG(JL)>0._JPRB
      ENDIF
      llresetjl(jl)=lldeep(jl).AND.llfirst(jl)
      !llreset=llreset.OR.llresetjl(jl)
    ENDDO

#ifndef _OPENACC
!RS this would be a reduction on GPUs
    llreset = ANY(llresetjl(kidia:kfdia))

    IF(llreset) THEN
#endif
      !$ACC LOOP SEQ
      DO jk=klev,ktdia,-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jkt, jkb)
        DO jl=kidia,kfdia
         ! keep first departure level that produces deep cloud
!          IF ( LLDEEP(JL) .AND. LLFIRST(JL) ) THEN 
          IF ( llresetjl(jl) ) THEN 
            jkt=ictop(jl)
            jkb=idpl(jl)
            IF ( jk<=jkb .AND. jk>=jkt ) THEN
              klab(jl,jk)=ilab(jl,jk)
              ptu(jl,jk)=ztu(jl,jk)
              pqu(jl,jk)=zqu(jl,jk)
              plu(jl,jk)=zlu(jl,jk)
            ELSE 
              klab(jl,jk)=1
              ptu(jl,jk)=ptenh(jl,jk)
              pqu(jl,jk)=pqenh(jl,jk)
              plu(jl,jk)=0.0_JPRB
            ENDIF
            IF ( jk<jkt ) klab(jl,jk)=0
          ENDIF
        ENDDO
      ENDDO
#ifndef _OPENACC
    ENDIF
#endif

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jkb)
    DO jl=kidia,kfdia
      IF ( lldeep(jl) .AND. llfirst(jl) ) THEN
        kdpl(jl)  = idpl(jl)
        kctop(jl) = ictop(jl)
        kcbot(jl) = icbot(jl)
        ldcum(jl) = lldcum(jl)
        ldsc(jl)  = .FALSE.
        kbotsc(jl)= -1
        jkb=kcbot(jl)
        pwubase(jl)=SQRT(MAX(zwu2h(jl,jkb),0.0_JPRB))
!  no initialization of wind for deep here, this is done in
!  CUINI and CUASCN
        llfirst(jl)=.FALSE.
      ENDIF
      llgo_on(jl) = .NOT.lldeep(jl)
    ENDDO
  ENDIF! jkk < klev

ENDDO ! end of big loop for search of departure level     

! compute maximum CAPE value
!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  pcape(jl) = 0._jprb
ENDDO

!$ACC LOOP SEQ
DO jk=ktdia,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    pcape(jl) = MAX(pcape(jl),zcape(jl,jk))
  ENDDO
ENDDO
!$ACC END PARALLEL

!$ACC WAIT(1)
!$ACC END DATA

IF (lhook) CALL dr_hook('CUBASEN',1,zhook_handle)
END SUBROUTINE cubasen



END MODULE mo_cuinit

