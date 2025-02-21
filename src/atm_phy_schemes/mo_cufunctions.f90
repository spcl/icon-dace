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
! functions for convection
!
!     ------------------------------------------------------------------
!     This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!       Two sets of functions are available. In the first set only the
!       cases water or ice are distinguished by temperature.  This set
!       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
!       The second set considers, besides the two cases water and ice
!       also a mix of both for the temperature range RTICE < T < RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.
!       FKOOP modifies the ice saturation mixing ratio for homogeneous
!       nucleation
!
!       Depending on the consideration of mixed phases either the first
!       set (e.g. surface, post-processing) or the second set
!       (e.g. clouds, condensation, convection) should be used.
!
! - Astronomical functions
!

MODULE mo_cufunctions

  USE mo_kind   ,ONLY: jprb=>wp     , &
    &                  jpim=>i4

!  USE yomhook   ,ONLY : lhook,   dr_hook
  
  
  USE mo_cuparameters, ONLY: r2es, r3les, r3ies, r4les, r4ies,   &
    & r5alvcp, r5alscp, &!r5les, r5ies,          &
    & ralvdcp, ralsdcp, rtice, rtwat,          &
    & rtwat_rtice_r,rticecu,rtwat_rticecu_r,   &
!    & rcdhalf, rcdhpi2, rcheta, rchetb, rchbb, &
!    & rchbcd,  rchbd, rchb23a, rchbbcd, rchba, &
!    & njkt1, njkt2,                            &
    & retv     ,rlvtt    ,rlstt    ,rtt,&
    & rd       ,rkappa   ,ratm ,               &
    & lhook,   dr_hook
  
  
  IMPLICIT NONE

  PRIVATE


  !! -  Thermodynamical functions (fcttre)
  !KF commented functions refer to obviously not used ones
  PUBLIC :: foedelta !! FOEDELTA = 1    water/ = 0    ice
  PUBLIC :: foeldcp
  PUBLIC :: foealfa
  PUBLIC :: foeewm
  PUBLIC :: foedem
  PUBLIC :: foeldcpm
  PUBLIC :: foedemcu
  PUBLIC :: foealfcu
  PUBLIC :: foeewmcu
  PUBLIC :: foeldcpmcu
  PUBLIC :: foelhmcu
  PUBLIC :: foeewm_v,foeewmcu_v,foeles_v,foeies_v
  PUBLIC :: foeewl, foeewi

  !! - Astronomical functions (fctast)
  !  PUBLIC:: RTETA,REL,REM,RRS,RLLS,RLLLS,RDS,RET
  !! - Time functions (fcttim)
  !  PUBLIC::  NDD,NMM,NCCAA,NAA,NAMD,NCTH,NZZAA,&
  !&           NZZMM,NCENT,NYEARC,&
  !&           NCONSTRUCT_DATE
  !  PUBLIC::  RJUDAT,RTIME
  !  PUBLIC::  RLV ! LATENT HEAT OF VAPOURISATION
  !  PUBLIC::  RLS ! LATENT HEAT OF SUBLIMATION
  !  PUBLIC::  RLF ! LATENT HEAT OF FUSION
  !  PUBLIC::  ESW ! SATURATION IN PRESENCE OF WATER
  !  PUBLIC::  ESS ! SATURATION IN PRESENCE OF ICE
  !  PUBLIC::  ES  ! SATURATION (IF T>RTT THEN WATER ; IF T<RTT THEN ICE)
  !  PUBLIC::  PHIMS,PHIHS !! - Stability functions (fcvdfs)
  !  PUBLIC:: CUBIDIAG

  PUBLIC:: cuancape2

  ! obviously not used functions
  !  PUBLIC :: FOEEWMO, FOEELIQ, FOEEICE
  !  PUBLIC :: FOELHM
  !  PUBLIC :: FOETB
  !  PUBLIC :: FOEEW
  !  PUBLIC :: FOEDE
  !  PUBLIC :: FOEDESU
  !  PUBLIC :: FOELH
  !  PUBLIC::  FOEW  ! FONCTION DE LA TENSION DE VAPEUR SATURANTE .
  !  PUBLIC::  FODLEW! FONCTION DERIVEE DU LOGARITHME NEPERIEN DE LA PRECEDENTE (FOEW)
  !  PUBLIC::  FOQS  ! FONCTION HUMIDITE SPECIFIQUE SATURANTE
  !  PUBLIC::  FODQS ! FONCTION DERIVEE EN TEMPERATURE DE LA PRECEDENTE (FOQS) .
  !  PUBLIC::  FOLH  ! FONCTION CHALEUR LATENTE .
  !  PUBLIC::  PHIHU,PHIMU,PSIHU,PSIMU
  !  PUBLIC::  PSIHS,PSIMS

CONTAINS

  !=======================================================================

  !!     ------------------------------------------------------------------
  !!     *****************************************************************

  !!                NO CONSIDERATION OF MIXED PHASES

  !!     *****************************************************************

  !     THERMODYNAMICAL FUNCTIONS .

  !     Pressure of water vapour at saturation
  !        INPUT : PTARE = TEMPERATURE


  ELEMENTAL FUNCTION foedelta (ptare)
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foedelta
    REAL(KIND=jprb), INTENT(in) :: ptare
    foedelta = MAX (0.0_JPRB,SIGN(1.0_JPRB,ptare-rtt))
  END FUNCTION foedelta

  ELEMENTAL FUNCTION foeldcp ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeldcp
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeldcp = foedelta(ptare)*ralvdcp + (1.0_JPRB-foedelta(ptare))*ralsdcp
  END FUNCTION foeldcp

  ELEMENTAL FUNCTION foeewl (ptare)
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeewl
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeewl=r2es*EXP(r3les*(ptare-rtt)/(ptare-r4les)) 
   END FUNCTION foeewl

  ELEMENTAL FUNCTION foeewi (ptare)
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeewi
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeewi=r2es*EXP(r3ies*(ptare-rtt)/(ptare-r4ies))
  END FUNCTION foeewi

  !FOEEW ( PTARE ) = R2ES*EXP (&
  !  &(R3LES*FOEDELTA(PTARE)+R3IES*(1.0_JPRB-FOEDELTA(PTARE)))*(PTARE-RTT)&
  !&/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE)))))

  !FOEDE ( PTARE ) = &
  !  &(FOEDELTA(PTARE)*R5ALVCP+(1.0_JPRB-FOEDELTA(PTARE))*R5ALSCP)&
  !&/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2

  !FOEDESU ( PTARE ) = &
  !  &(FOEDELTA(PTARE)*R5LES+(1.0_JPRB-FOEDELTA(PTARE))*R5IES)&
  !&/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(1.0_JPRB-FOEDELTA(PTARE))))**2

  !FOELH ( PTARE ) =&
  !         &FOEDELTA(PTARE)*RLVT + (1.0_JPRB-FOEDELTA(PTARE))*RLSTT

  !!     *****************************************************************

  !!           CONSIDERATION OF MIXED PHASES

  !!     *****************************************************************

  !!     FOEALFA is calculated to distinguish the three cases:

  !!                       FOEALFA=1            water phase
  !!                       FOEALFA=0            ice phase
  !!                       0 < FOEALFA < 1      mixed phase

  !!               INPUT : PTARE = TEMPERATURE


  ELEMENTAL FUNCTION foealfa (ptare)
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foealfa
    REAL(KIND=jprb), INTENT(in) :: ptare
    foealfa  = MIN(1.0_JPRB,((MAX(rtice,MIN(rtwat,ptare))-rtice)&
      & * rtwat_rtice_r)**2)
  END FUNCTION foealfa







  ELEMENTAL FUNCTION foeewm ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeewm
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeewm  = r2es                                * &
      & (foealfa(ptare)*EXP(r3les*(ptare-rtt)/(ptare-r4les))+ &
      & (1.0_JPRB-foealfa(ptare))*EXP(r3ies*(ptare-rtt)/(ptare-r4ies)))
  END FUNCTION foeewm

  ELEMENTAL FUNCTION foedem ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foedem
    REAL(KIND=jprb), INTENT(in) :: ptare
    foedem  = foealfa(ptare)*r5alvcp*(1.0_JPRB/(ptare-r4les)**2)+&
         & (1.0_JPRB-foealfa(ptare))*r5alscp*(1.0_JPRB/(ptare-r4ies)**2)
  END FUNCTION foedem

  ELEMENTAL FUNCTION foeldcpm ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeldcpm
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeldcpm  = foealfa(ptare)*ralvdcp+&
         & (1.0_JPRB-foealfa(ptare))*ralsdcp
  END FUNCTION foeldcpm

  !FOELHM ( PTARE ) =&
  !         &FOEALFA(PTARE)*RLVTT+(1.0_JPRB-FOEALFA(PTARE))*RLSTT







  !FOETB ( PTARE )=FOEALFA(PTARE)*R3LES*(RTT-R4LES)*(1.0_JPRB/(PTARE-R4LES)**2)+&
  !             &(1.0_JPRB-FOEALFA(PTARE))*R3IES*(RTT-R4IES)*(1.0_JPRB/(PTARE-R4IES)**2)


  ELEMENTAL FUNCTION foealfcu ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foealfcu
    REAL(KIND=jprb), INTENT(in) :: ptare
    foealfcu = MIN(1.0_JPRB,((MAX(rticecu,MIN(rtwat,ptare))&
         & -rticecu)*rtwat_rticecu_r)**2)
  END FUNCTION foealfcu


  ELEMENTAL FUNCTION foeewmcu ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeewmcu
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeewmcu  = r2es *                                              &
      & (foealfcu(ptare)*EXP(r3les*(ptare-rtt)/(ptare-r4les))+ &
      & (1.0_JPRB-foealfcu(ptare))*EXP(r3ies*(ptare-rtt)/(ptare-r4ies)))
  END FUNCTION foeewmcu

  !     Pressure of water vapour at saturation
  !        INPUT : PTARE = TEMPERATURE
  !

  ELEMENTAL FUNCTION foedemcu ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foedemcu
    REAL(KIND=jprb), INTENT(in) :: ptare
    foedemcu =foealfcu(ptare)*r5alvcp*(1.0_JPRB/(ptare-r4les)**2)+&
         & (1.0_JPRB-foealfcu(ptare))*r5alscp*(1.0_JPRB/(ptare-r4ies)**2)
  END FUNCTION foedemcu

  ELEMENTAL FUNCTION foeldcpmcu ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeldcpmcu
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeldcpmcu  = foealfcu(ptare)*ralvdcp+&
         & (1.0_JPRB-foealfcu(ptare))*ralsdcp
  END FUNCTION foeldcpmcu

  ELEMENTAL FUNCTION foelhmcu ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foelhmcu
    REAL(KIND=jprb), INTENT(in) :: ptare
    foelhmcu  =&
      & foealfcu(ptare) &
      & *rlvtt+(1.0_JPRB-foealfcu(ptare))*rlstt
  END FUNCTION foelhmcu

  !     ------------------------------------------------------------------

  !FOEEWMO( PTARE ) = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  !FOEELIQ( PTARE ) = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  !FOEEICE( PTARE ) = R2ES*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES))

  ELEMENTAL FUNCTION foeles_v ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeles_v
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeles_v =r3les*(ptare-rtt)/(ptare-r4les)
  END FUNCTION foeles_v

  ELEMENTAL FUNCTION foeies_v ( ptare )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeies_v
    REAL(KIND=jprb), INTENT(in) :: ptare
    foeies_v =r3ies*(ptare-rtt)/(ptare-r4ies)
  END FUNCTION foeies_v

  ELEMENTAL FUNCTION foeewm_v ( ptare,exp1,exp2 )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeewm_v
    REAL(KIND=jprb), INTENT(in) :: ptare,exp1,exp2
    foeewm_v =r2es*(foealfa(ptare)*exp1+ &
         & (1.0_JPRB-foealfa(ptare))*exp2)
  END FUNCTION foeewm_v

  ELEMENTAL FUNCTION foeewmcu_v ( ptare,exp1,exp2 )
    !$ACC ROUTINE SEQ
    REAL(KIND=jprb)             :: foeewmcu_v
    REAL(KIND=jprb), INTENT(in) :: ptare,exp1,exp2
    foeewmcu_v  = r2es*(foealfcu(ptare)*exp1+&
         & (1.0_JPRB-foealfcu(ptare))*exp2)
  END FUNCTION foeewmcu_v

  !     ------------------------------------------------------------------
  !     FONCTIONS THERMODYNAMIQUES : FONCTIONS DEFINIES DE LA PHYSIQUE .


  !     FONCTION DE LA TENSION DE VAPEUR SATURANTE .
  !        INPUT : PTARG = TEMPERATURE
  !                PDELARG = 0 SI EAU (QUELQUE SOIT PTARG)
  !                          1 SI GLACE (QUELQUE SOIT PTARG).

  !FOEW ( PTARG,PDELARG ) = EXP (&
  !    &( RALPW+PDELARG*RALPD )&
  !  &- ( RBETW+PDELARG*RBETD ) / PTARG &
  !  &- ( RGAMW+PDELARG*RGAMD ) * LOG(PTARG) )

  !     FONCTION DERIVEE DU LOGARITHME NEPERIEN DE LA PRECEDENTE (FOEW) .
  !        INPUT : PTARG = TEMPERATURE
  !                PDELARG = 0 SI EAU (QUELQUE SOIT PTARG)
  !                          1 SI GLACE (QUELQUE SOIT PTARG).

  !FODLEW ( PTARG,PDELARG ) = (&
  !      &( RBETW+PDELARG*RBETD )&
  !    &- ( RGAMW+PDELARG*RGAMD ) * PTARG )&
  !    &/ ( PTARG*PTARG )

  !     FONCTION HUMIDITE SPECIFIQUE SATURANTE .
  !        INPUT : PESPFAR = RAPPORT FOEW SUR PRESSION.

  !FOQS ( PESPFAR ) = PESPFAR / ( 1.0_JPRB+RETV*MAX(0.0_JPRB,&
  !    &(1.0_JPRB-PESPFAR)) )

  !     FONCTION DERIVEE EN TEMPERATURE DE LA PRECEDENTE (FOQS) .
  !        INPUT : PQSFARG = FOQS
  !                PESPFAR = RAPPORT FOEW SUR PRESSION
  !                PDLEFAR = FODLEW.

  !FODQS ( PQSFARG,PESPFAR,PDLEFAR ) = ( PQSFARG &
  !   &* (1.0_JPRB-PQSFARG)*PDLEFAR ) / (1.0_JPRB-PESPFAR)

  !     FONCTION CHALEUR LATENTE .
  !        INPUT : PTARG = TEMPERATURE
  !                PDELARG = 0 SI EAU (QUELQUE SOIT PTARG)
  !                          1 SI GLACE (QUELQUE SOIT PTARG).

  !FOLH ( PTARG,PDELARG ) =  RV * (&
  !    &( RBETW+PDELARG*RBETD )&
  !  &- ( RGAMW+PDELARG*RGAMD ) * PTARG )
  !     -------------------------------------


  !     ------------------------------------------------------------------

  !     *FCVDFS** CONTAINS STATEMENT FUNCTIONS DESCRIBING STAB. FUNCT.

  !     A.C.M. BELJAARS    E.C.M.W.F.      26/03/90.

  !     ------------------------------------------------------------------

  !          *THE STABILITY FUNCTIONS ARE THE SO-CALLED *PHI* AND
  !     *PSI*-FUNCTIONS. THE *PSI*-FUNCTIONS GIVE THE STABILITY
  !     CORRECTIONS IN THE LOGARITHMIC PROFILES FOR
  !     WIND, DRY STATIC ENERGY AND SPECIFIC HUMIDITY. THE FUNCTIONS
  !     DEPEND ON THE RATIO OF HEIGHT AND *OBUKHOV LENGTH (*ETA*).
  !          FOR THE UNSTABLE BOUNDARY LAYER, THE *DYER AND *HICKS
  !     FORMULATIONS ARE USED (CF. *DYER, 1974; *HOGSTROM, 1988). IN
  !     STABLE SITUATIONS, THE EMPIRICAL FORMS, PROPOSED BY *HOLTSLAG
  !     AND *DEBRUIN ARE USED WITH A MODIFICATION TO SATISFY A CRITICAL
  !     FLUX-*RICHARDSON NUMBER FOR LARGE *ETA*.
  !          THE *PHI* AND *PSI* FUNCTIONS ARE INTERRELATED. THE *PSI*
  !     FUNCTIONS CAN BE DERIVED FROM THE *PHI* FUNCTIONS BY INTEGRATION
  !     OF (1.-PHI)/ETA OR *PHI* FROM *PSI* BY COMPUTING
  !     (1.-ETA*DPSI/DETA) (SEE ALSO *HAUGEN, 1973; WORKSHOP ON
  !     MICROMETEOROLOGY, P. 77).

  !     ------------------------------------------------------------------

  !        *PHI AND *PSI FUNCTIONS FOR UNSTABLE SITUATIONS ACCORDING
  !        TO *DYER AND *HICKS
  !        3/2/03 Modification, J.Hague  X**1.5 -> X*SQRT(X)

  !PHIHU(PETA)= 1.0_JPRB/     SQRT(1.0_JPRB-RCDHALF*PETA)
  !PHIMU(PETA)= 1.0_JPRB/SQRT(SQRT(1.0_JPRB-RCDHALF*PETA))

  !PSIHU(PETA)= 2.0_JPRB*LOG((1.0_JPRB+     SQRT(1.0_JPRB-RCDHALF*PETA))*0.5_JPRB )
  !PSIMU(PETA)=    LOG((1.0_JPRB+SQRT(SQRT(1.0_JPRB-RCDHALF*PETA)))**2 &
  !                  &*(1.0_JPRB+     SQRT(1.0_JPRB-RCDHALF*PETA) ) *0.125_JPRB )&
  !                  &-2.0_JPRB*ATAN(SQRT(SQRT(1.0_JPRB-RCDHALF*PETA)))&
  !                  &+ RCDHPI2


  !>KF functions PHIHS and PHIMS are moved to mo_cuparameters to avoid circular
  !  dependenies
  !        *PHI AND *PSI FUNCTIONS FOR UNSTABLE SITUATIONS ACCORDING
  !        TO HOGSTROM FOR MOMENTUM AND DERIVED FROM THE ELLISON AND
  !        TURNER RELATION FOR THE RATIO OF PHIM AMD PHIH.

  ! ELEMENTAL FUNCTION PHIMS(PETA)
  !   REAL(KIND=JPRB)             :: PHIMS
  !   REAL(KIND=JPRB), INTENT(IN) :: PETA
  !   PHIMS = 1.0_JPRB+RCHETA*PETA
  ! END FUNCTION PHIMS
  !
  ! ELEMENTAL FUNCTION PHIHS(PETA)
  !   REAL(KIND=JPRB)             :: PHIHS
  !   REAL(KIND=JPRB), INTENT(IN) :: PETA
  !   PHIHS = (1.0_JPRB+RCHETB*PETA)**2
  ! END FUNCTION PHIHS

  !     PSI FUNCTIONS NOT COMPATIBLE WITH PHI FUNCTIONS IN THE STABLE CASE


  !PSIHS(PETA)= -RCHBB*(PETA-RCHBCD)*EXP  (-RCHBD*PETA)&
  !               &-(_ONE_+RCHB23A*PETA)**1.5_JPRB - RCHBBCD + _ONE_
  !PSIHS(PETA)= -RCHBB*(PETA-RCHBCD)*EXP  (-RCHBD*PETA)&
  !              &-(1.0_JPRB+RCHB23A*PETA)*SQRT(1.0_JPRB+RCHB23A*PETA)-RCHBBCD+1.0_JPRB
  !PSIMS(PETA)= -RCHBB*(PETA-RCHBCD)*EXP  (-RCHBD*PETA)&
  !               &-RCHBA*PETA - RCHBBCD

  !=======================================================================


  SUBROUTINE cuancape2(kidia,  kfdia,  klon,   klev, njkt1, njkt2, &
    & pap,    paph,   pt,    pq,     pcape)
    !
    ! Description:
    !***** CUANCAPE2 - COMPUTE APPROXIMATE CAPE USING THETAE AND THETAES
    !
    !     E. HOLM + P. BECHTOLD     E.C.M.W.F.     13/10/2005

    !     PURPOSE
    !     -------
    !                 ESTIMATE CAPE FIRST FOR A MIXED-LAYER PARCEL, THEN
    !                 LOOP OVER SUBSEQUENT DEPARTURE LAYERS IN LOWEST 350 hPa
    !                 Theta_e =Theta*exp[L q_v/(C_p T)]
    !                         = T*(P0/P)**(R_d/C_p) * exp[L q_v/(C_p T)]
    !                 -> THIS WILL BE THE UPDRAUGHT PARCEL (CONSERVING ITS
    !                 PROPERTIES)  (no entrainment)
    !                 CAPE    = Int ( g dTheta_v/Theta_v dz ) =
    !                   aprox = Int ( g (Theta_e_up-Theta_e_sat)/Theta_e_sat ) dz
    !                 WITH THIS FORMULATION THE ACTUAL CAPE IS OVERESTIMATED  BY
    !                 ROUGHLY 20%. DEEP CONVECTION CAN BE CONSIDERED FOR CAPE
    !                 VALUES ABOVE 200-500 J/KG

    !     PARAMETER     DESCRIPTION                              UNITS
    !     ---------     -----------                              -----
    !     INPUT PARAMETERS (INTEGER):

    !    *KIDIA*        START POINT
    !    *KFDIA*        END POINT
    !    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !    *KLEV*         NUMBER OF LEVELS

    !    INPUT PARAMETERS (REAL):

    !    *PAP*          PRESSURE ON FULL LEVELS                    PA
    !    *PAPH*         PRESSURE ON HALF LEVELS                    PA
    !    *PT*           TEMPERATURE ON FULL LEVELS                 K
    !    *PQ*           SPECIFIC HUMIDITY ON FULL LEVELS          KG/KG

    !    OUTPUT PARAMETERS (REAL):

    !    *CAPE*         CONVECTIVE AVAILABLE POT. ENERGY          J/KG

    !-------------------------------------------------------------------------------

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
    !USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
    ! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
    ! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
    ! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R
    !USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT,&
    !                     &RD       ,RKAPPA   ,RATM
    !USE YOECUMF  , ONLY : NJKT1, NJKT2

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: njkt1, njkt2
    REAL(KIND=jprb)   ,INTENT(in)    :: pap(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: pt(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(out)   :: pcape(klon)

    INTEGER(KIND=jpim) :: jl, jk, jkk

    REAL(KIND=jprb), DIMENSION(klon) :: zpmix, ztmix, zthmix, zqmix, ztheteu
    REAL(KIND=jprb)                  :: zcape(klon,klev)
    REAL(KIND=jprb) :: zdp, zthetes, zqs, zdz, ztemp, zrpap

    REAL(KIND=jprb) :: zhook_handle

    !#include "fcttre.h"
    !-------------------------------------------------------------------------------
    IF (lhook) CALL dr_hook('CUANCAPE2',0,zhook_handle)

    pcape(:)=0.0_JPRB

    DO jkk=klev-1,njkt1,-1

      DO jl=kidia,kfdia
        zcape(jl,jkk)=0.0_JPRB
        IF (paph(jl,klev+1)-paph(jl,jkk-1)<60.e2_jprb) THEN
          ztmix(jl)=0.0_JPRB
          zthmix(jl)=0.0_JPRB
          zqmix(jl)=0.0_JPRB
          zpmix(jl)=0.0_JPRB
          DO jk=jkk+1,jkk-1,-1
            IF(zpmix(jl)<30.e2_jprb) THEN
              zdp=paph(jl,jk+1)-paph(jl,jk)
              zpmix(jl)=zpmix(jl)+zdp
              zthmix(jl)=zthmix(jl)+pt(jl,jk)*zdp*(ratm/pap(jl,jk))**rkappa
              zqmix(jl)=zqmix(jl)+pq(jl,jk)*zdp
            ENDIF
          ENDDO
          zdp=1.0_JPRB/zpmix(jl)
          zqmix(jl)=zqmix(jl)*zdp
          zpmix(jl)=paph(jl,jkk+2)-0.5_JPRB*zpmix(jl)
          zthmix(jl)=zthmix(jl)*zdp
          ztmix(jl)=zthmix(jl)*(zpmix(jl)/ratm)**rkappa
        ELSE
          zqmix(jl)=pq(jl,jkk)
          zpmix(jl)=pap(jl,jkk)
          ztmix(jl)=pt(jl,jkk)
          zthmix(jl)=pt(jl,jkk)*(ratm/zpmix(jl))**rkappa
        ENDIF
        ztheteu(jl)=zthmix(jl)*&
          & EXP( foeldcp(ztmix(jl))*zqmix(jl)/ztmix(jl) )
      ENDDO

      DO jk=jkk,njkt2,-1
        DO jl=kidia,kfdia
          IF(pap(jl,jk)>80.e2_jprb.AND. &
            & (paph(jl,klev+1)-paph(jl,jkk))<350.e2_jprb) THEN
            zrpap=1.0_JPRB/pap(jl,jk)
            zqs = foeewm(pt(jl,jk))*zrpap
            zqs = MAX(1.e-8_JPRB,zqs)
            zqs = zqs/(1.0_JPRB-retv*zqs) ! small correction
            zthetes=pt(jl,jk)*(ratm*zrpap)**rkappa*&
              & EXP( foeldcp(pt(jl,jk))*zqs/pt(jl,jk) )
            ztemp=ztheteu(jl)/zthetes-1.0_JPRB
            IF(ztemp > 0.0_JPRB)THEN
              zdz=(paph(jl,jk+1)-paph(jl,jk))*zrpap*rd*pt(jl,jk)*&
                & (1.0_JPRB+retv*pq(jl,jk))
              zcape(jl,jkk)=zcape(jl,jkk)+ztemp*zdz
            ENDIF
          ENDIF
        ENDDO
      ENDDO

    ENDDO

    ! chose maximum CAPE value
    DO jl=kidia,kfdia
      pcape(jl)=MAXVAL(zcape(jl,njkt1:klev-1))
    ENDDO

    !-------------------------------------------------------------------------------

    IF (lhook) CALL dr_hook('CUANCAPE2',1,zhook_handle)
  END SUBROUTINE cuancape2



END MODULE mo_cufunctions

