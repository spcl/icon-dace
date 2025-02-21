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

!  cuentr: CALCULATES ENTRAINMENT/DETRAINMENT RATES FOR UPDRAFTS
!          IN CUMULUS PARAMETERIZATION

! The following subroutine is not used in the current setup
! *CUSTRAT* - COMPUTES T,Q TENDENCIES FOR STRATOCUMULUS
!                 CONVECTION

MODULE mo_cuascn

  USE mo_kind   ,ONLY: jprb=>wp     , &
    &                  jpim=>i4

!  USE yomhook   ,ONLY : lhook,   dr_hook

  USE mo_cuparameters, ONLY: rg   ,rcpd,  retv  ,  &
    &                        rtt   ,rd   , ralfdcp,&
    &                        rtber, rtbercu  ,rticecu   ,&
    &                        lphylin  ,rlptrc,           &
    &                        entshalp ,rmfcmin,          &
    &                        rmflic   ,rmflia ,rvdifts  ,&
    &                        rmfcmax, rlmin             ,&
    &                        lhook,   dr_hook, lmfglac


  USE mo_adjust ,ONLY: cuadjtq

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: cuascn

CONTAINS

  !
  !OPTIONS XOPT(HSFUN)
  !
  SUBROUTINE cuascn &
    & ( kidia,    kfdia,    klon,    ktdia,  klev, rmfcfl, &
    & entrorg, detrpen, rprcon, lmfmid, lgrz_deepconv, ptsphy,&
    & paer_ss,  ptenh,    pqenh,    ptenq,             &
    & pten,     pqen,     pqsen,    plitot,&
    & pgeo,     pgeoh,    pap,      paph,&
    & zdph,     zdgeoh,                  &
    & pvervel,  pwubase,  pcloudnum, pcore,     &
    & ldland,   ldlake,   ldcum,    ktype,    klab,&
    & ptu,      pqu,      plu,      plrain,        &
    & pmfu,     pmfub,    plglac,&
    & pmfus,    pmfuq,    pmful,    plude,    pdmfup,&
    & pdmfen,   pcape,    pcapethresh,  &
    & kcbot,    kctop,    kctop0,   kdpl,     pmfude_rate,    pkineu,  pwmean, lacc )

!>
!!          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!!          FOR CUMULUS PARAMETERIZATION

!!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89

!!          PURPOSE.
!!          --------
!!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!!           FLUXES AS WELL AS PRECIPITATION RATES)

!!          INTERFACE
!!          ---------

!!          THIS ROUTINE IS CALLED FROM *CUMASTR*.

!!          METHOD.
!!          --------
!!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!!          AND THEN CALCULATE MOIST ASCENT FOR
!!          ENTRAINING/DETRAINING PLUME.
!!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!!          SHALLOW AND DEEP CUMULUS CONVECTION.
!!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)

!!     PARAMETER     DESCRIPTION                                   UNITS
!!     ---------     -----------                                   -----
!!     INPUT PARAMETERS (INTEGER):

!!    *KIDIA*        START POINT
!!    *KFDIA*        END POINT
!!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!!    *KTDIA*        START OF THE VERTICAL LOOP
!!    *KLEV*         NUMBER OF LEVELS
!!    *KTYPE*        TYPE OF CONVECTION
!!                       1 = PENETRATIVE CONVECTION
!!                       2 = SHALLOW CONVECTION
!!                       3 = MIDLEVEL CONVECTION
!!    *KCBOT*        CLOUD BASE LEVEL
!!    *KDPL*         DEPARTURE LEVEL FOR CONVECTION

!!    INPUT PARAMETERS (REAL):

!!    *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
!!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
!!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
!!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)      M/S
!!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)      M/S
!!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)      K
!!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1) KG/KG
!!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)  KG/KG
!!    *PGEO*         GEOPOTENTIAL                                 M2/S2
!!    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT           KG/KG
!!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
!!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS           PA
!!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!!    *PTENQ*        MOISTURE TENDENCY                            KG/(KG S)
!!    *PVERVEL*      VERTICAL VELOCITY                            PA/S
!!    *zdgeoh*       geopot thickness on full levels              M2/S2
!!    *zdgeo*        geopot thickness on full levels               m2/s2
!!    *zdph*         pressure thickness on full levels             PA
!!    *pcape*        CAPE                                          J/kg
!!    *pcapethresh*  CAPE threshold beyond which entrainment parameter is reduced
!!    *pcore*        Updraft core fraction profile                 0-1    

!!    INPUT PARAMETERS (LOGICAL):

!!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!!    *LDLAKE*       LAKE MASK (.TRUE. FOR LAKE)
!!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

!!    UPDATED PARAMETERS (INTEGER):

!!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!!                        KLAB=2 FOR CLOUD LEVELS

!!    UPDATED PARAMETERS (REAL):

!!    *PTU*          TEMPERATURE IN UPDRAFTS                        K
!!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
!!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG
!     *PLRAIN*       RAIN WATER CONTENT IN UPDRAFTS               KG/KG

!!    OUTPUT PARAMETERS (INTEGER):

!!    *KCTOP*        CLOUD TOP LEVEL
!!    *KCTOP0*       FIRST GUESS OF CLOUD TOP LEVEL

!!    OUTPUT PARAMETERS (REAL):

!!    *PMFU*         MASSFLUX IN UPDRAFTS                         KG/(M2*S)
!!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)
!!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS         J/(M2*S)
!!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS           KG/(M2*S)
!!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS             KG/(M2*S)
!!    *PLUDE*        DETRAINED LIQUID WATER                       KG/(M2*S)
!!    *PLGLAC*       FROZEN CLOUD WATER CONTENT                   KG/KG
!!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS       KG/(M2*S)
!!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                     KG/(M2*S)
!!    *PKINEU*       UPDRAFT KINETIC ENERGY                       M2/S2
!!    *PWMEAN*       MEAN UPDRAUGHT VELOCITY                      M/S

!!          EXTERNALS
!!          ---------
!!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!!          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION

!!          REFERENCE
!!          ---------
!!          (TIEDTKE,1989)

!!          MODIFICATIONS
!!          -------------
!!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!!             99-06-14 : Optimisation        D.SALMOND
!!             01-05-22 : Modified flux limiter M.CULLEN
!!             02-08-14 : Allow for departure level =/ KLEV  P.BECHTOLD
!!             03-08-28 : Clean-up detrainment rates         P.BECHTOLD
!!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!!        J.Hague       08-Dec-2005 Tuning: LLFLAG indexing
!!             07-06-01 : Organized entrainment based on RH  P.BECHTOLD

!!----------------------------------------------------------------------

!USE PARKIND1  ,ONLY : JPIM     ,JPRB
!USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!USE YOMCST   , ONLY : RG       ,RCPD     ,RETV     ,RLVTT    ,RLSTT    ,RTT    
!USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
! & RALVDCP  ,RALSDCP  ,RALFDCP  ,RTWAT    ,RTBER    ,&
! & RTBERCU  ,RTICE    ,RTICECU  ,&
! & RTWAT_RTICECU_R    ,RTWAT_RTICE_R  
!USE YOECUMF  , ONLY : ENTRORG  ,ENTSHALP  ,RMFCMIN  ,RPRCON   ,RMFCFL   ,RMFLIC   ,RMFLIA
!USE YOEPHLI  , ONLY : LPHYLIN  ,RLPTRC
!USE YOECLDP  , ONLY : RLMIN    ,LAERLIQAUTOCP, LAERLIQAUTOCPB
!USE YOM_YGFL , ONLY : NACTAERO

  
USE mo_cufunctions       , ONLY: foealfcu
! GZ, 2013-09-13: tuning to reduce drizzle, and coupling of autoconversion to aerosols
USE mo_atm_phy_nwp_config, ONLY: ltuning_kessler, icpl_aero_conv


INTEGER(KIND=jpim),INTENT(in)    :: klon 
INTEGER(KIND=jpim),INTENT(in)    :: klev 
INTEGER(KIND=jpim),INTENT(in)    :: kidia 
INTEGER(KIND=jpim),INTENT(in)    :: kfdia 
INTEGER(KIND=jpim),INTENT(in)    :: ktdia
REAL(KIND=jprb)   ,INTENT(in)    :: rmfcfl 
REAL(KIND=jprb)   ,INTENT(in)    :: entrorg, rprcon, detrpen
LOGICAL           ,INTENT(in)    :: lmfmid, lgrz_deepconv
REAL(KIND=jprb)   ,INTENT(in)    :: ptsphy 
!KF
REAL(KIND=jprb)   ,INTENT(in), OPTIONAL:: paer_ss(klon)
!KF
REAL(KIND=jprb)   ,INTENT(inout) :: ptenh(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: pqenh(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pten(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pqen(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pqsen(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: plitot(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(in)    :: zdgeoh(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pap(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(in)    :: zdph(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: ptenq(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pvervel(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pwubase(klon) 
REAL(KIND=jprb)   ,INTENT(in)    :: pcloudnum(klon) 
REAL(KIND=jprb)   ,INTENT(in)    :: pcore(klon,klev) 
LOGICAL           ,INTENT(in)    :: ldland(klon) 
LOGICAL           ,INTENT(in)    :: ldlake(klon) 
LOGICAL           ,INTENT(inout) :: ldcum(klon) 
INTEGER(KIND=jpim),INTENT(inout) :: ktype(klon) 
INTEGER(KIND=jpim),INTENT(inout) :: klab(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: ptu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: pqu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: plu(klon,klev) 
REAL(KIND=JPRB)   ,INTENT(inout) :: plrain(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: pmfu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: pmfub(klon) 
REAL(KIND=jprb)   ,INTENT(out)   :: plglac(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: pmfus(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: pmfuq(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: pmful(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: plude(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: pdmfup(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: pdmfen(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pcape(klon)
REAL(KIND=jprb)   ,INTENT(in)    :: pcapethresh
INTEGER(KIND=jpim),INTENT(inout) :: kcbot(klon) 
INTEGER(KIND=jpim),INTENT(out)   :: kctop(klon) 
INTEGER(KIND=jpim),INTENT(inout) :: kctop0(klon) 
INTEGER(KIND=jpim),INTENT(in)    :: kdpl(klon) 
REAL(KIND=jprb)   ,INTENT(inout) :: pmfude_rate(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: pkineu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(out)   :: pwmean(klon) 
LOGICAL           ,INTENT(in)    :: lacc

REAL(KIND=jprb) ::     zdmfen(klon), zdmfde(klon),&
 & zqold(klon),  &
 & zbuo(klon,klev),    zluold(klon),&
 & zprecip(klon)  
REAL(KIND=jprb) ::     zdpmean(klon)
REAL(KIND=jprb) ::     zoentr(klon), zph(klon), zpbase(klon), zptop0(klon), zttop0(klon)
REAL(KIND=jprb) ::     zcrit(klon), zdrain(klon), zdnoprc(klon), zentrorg(klon)
LOGICAL ::  llflag(klon), llflaguv(klon), llo1(klon), llo3

INTEGER(KIND=jpim) :: icall, ik, is, jk, jl, ikb, kk
INTEGER(KIND=jpim) :: jll

REAL(KIND=jprb) :: z_cldmax, z_cprc2, z_cwdrag, z_cwifrac, zalfaw,&
 & zbc(klon), zbe, zbuoc, zc, zcbf, zcons2, zd, zdfi, &
 & zdkbuo, zdken, zrg, &
 & zdt, zfac, zfacbuo, zint, zkedke, zlcrit, &
 & zleen, zlnew, zmfmax, zmftest, zmfulk, zmfun, &
 & zmfuqk, zmfusk, zoealfa, zoealfap, zprcdgw, &
 & zprcon, zqeen, zqude, zrnew, zrold, zscde, &
 & zseen, ztglace, zvi, zvv, zvw, zwu, zzco, zzzmb, zdz, zmf, zglac, zentr_prof

REAL(KIND=jprb) ::  zchange,zxs,zxe
REAL(KIND=jprb) :: zhook_handle

!KF seasalt global mean, marine and continental
REAL(KIND=jprb) :: aer_ss_gm,aer_ss_m, aer_ss_c

LOGICAL :: llklab(klon)
! Set .true. if using updraft core fraction from explicit stochastic
! scheme to modify the lateral detrainment profile
LOGICAL, PARAMETER :: lstochdetr=.FALSE.

!#include "cuadjtq.intfb.h"
!#include "cubasmcn.intfb.h"
!#include "cuentr.intfb.h"

!#include "fcttre.func.h"
!----------------------------------------------------------------------

!*    1.           SPECIFY PARAMETERS
!                  ------------------

IF (lhook) CALL dr_hook('CUASCN',0,zhook_handle)

!$ACC DATA &
!$ACC   PRESENT(ptenh, pqenh, pten, pqen, pqsen, plitot, pgeo, pgeoh, zdgeoh, pap, paph) &
!$ACC   PRESENT(zdph, ptenq, pvervel, pwubase, pcloudnum, pcore, ldland, ldlake, ldcum, ktype) &
!$ACC   PRESENT(klab, ptu, pqu, plu, plrain, pmfu, pmfub, plglac, pmfus, pmfuq, pmful, plude) &
!$ACC   PRESENT(pdmfup, pdmfen, pcape, kcbot, kctop, kctop0, kdpl, pmfude_rate, pkineu, pwmean) &

!$ACC   CREATE(zdmfen, zdmfde, zqold, zbuo, zluold, zprecip, zdpmean, zoentr, zph, zpbase) &
!$ACC   CREATE(zptop0, zttop0, zcrit, zdrain, zdnoprc, zentrorg, llflag, llflaguv, llo1) &
!$ACC   CREATE(llklab, zbc) &
!$ACC   IF(lacc)


zcons2=rmfcfl/(rg*ptsphy)
ZRG=1.0_JPRB/RG
ztglace=rtt-13._jprb
zfacbuo=0.5_JPRB/(1.0_JPRB+0.5_JPRB)
zprcdgw=rprcon/rg
z_cldmax=5.e-3_JPRB
z_cwifrac=0.5_JPRB
z_cprc2=0.5_JPRB
z_cwdrag=(3._jprb/8._jprb)*0.506_JPRB/0.2_JPRB

IF(lmfglac) THEN
  zglac=0.5_JPRB
ELSE
  zglac=0.0_JPRB
ENDIF

!----------------------------------------------------------------------

!     2.           SET DEFAULT VALUES
!                  ------------------

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

!$ACC LOOP SEQ
DO jk=1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    !dimension: (klon,klev)
    plglac     (jl,jk) = 0.0_JPRB
    pmfus      (jl,jk) = 0.0_JPRB
    pmfuq      (jl,jk) = 0.0_JPRB
    pmful      (jl,jk) = 0.0_JPRB
    plude      (jl,jk) = 0.0_JPRB 
    pdmfup     (jl,jk) = 0.0_JPRB 
    pdmfen     (jl,jk) = 0.0_JPRB 
    pmfude_rate(jl,jk) = 0.0_JPRB 
    pkineu     (jl,jk) = 0.0_JPRB
    pmfu       (jl,jk) = 0.0_JPRB
    plrain     (jl,jk) = 0.0_JPRB
    zbuo       (jl,jk) = 0.0_JPRB
  ENDDO
ENDDO

!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  !dimension: (klon)
  pwmean       (jl) = 0.0_JPRB
  kctop        (jl) = 0
ENDDO


llo3=.FALSE.


!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  zluold(jl)=0.0_JPRB
  IF(.NOT.ldcum(jl)) THEN
    kcbot(jl)=-1
    pmfub(jl)=0.0_JPRB
    pqu(jl,klev)=0.0_JPRB
    ktype(jl)=0
  ENDIF
  ! pwmean(jl)=0.0_JPRB
  zdpmean(jl)=0.0_JPRB
  zoentr(jl)=0.0_JPRB
ENDDO

! initalize various quantities
! note that liquid water and kinetic energy at cloud base is 
! preserved from cubase

!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  llklab(jl)=.FALSE.
  IF(.NOT.ldcum(jl).OR.ktype(jl) == 3) llklab(jl)=.TRUE.
ENDDO

!$ACC LOOP SEQ
DO jk=ktdia,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    IF (jk /= kcbot(jl)) THEN 
      plu(jl,jk)=0.0_JPRB
    ENDIF
    IF( llklab(jl) ) klab(jl,jk)=0
    IF(.NOT.ldcum(jl).AND.paph(jl,jk) < 4.e4_jprb) kctop0(jl)=MAX(jk,2)
    ! pkineu(jl,jk)=0.0_JPRB
  ENDDO
ENDDO

!DIR$ IVDEP
!OCL NOVREC
!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  IF(ktype(jl) == 3) ldcum(jl)=.FALSE.
  ! Reduce entrainment in case of extreme CAPE in order to prevent numerical instabilities
  IF (pcape(jl) > pcapethresh) THEN
    zentrorg(jl) = entrorg*MAX(0.5_jprb,(pcapethresh/pcape(jl))**2)
  ELSE
    zentrorg(jl) = entrorg
  ENDIF
ENDDO

!----------------------------------------------------------------------

!IF(PRESENT (paer_ss)) THEN
!   
!      !!KF
!      !> define the seasonal dependend seasalt threshold instead of land/seamask
!      !!  derived from the annual mean
!      !!  aer_ss_c = seasalt_global and annMean*0.4
!      !!  aer_ss_m = seasalt_annMean*1.1
!
!      aer_ss_gm = 0.0066384_JPRB
!      aer_ss_c  = aer_ss_gm*0.4_JPRB
!      aer_ss_m  = aer_ss_gm*1.1_JPRB
!
!      DO jl=kidia,kfdia
!        ! to be used as a weight
!        !!
!        zcrit(jl) =(paer_ss(jl)-aer_ss_c)/(aer_ss_m - aer_ss_c)
!        zcrit(jl) = MAX(zcrit(jl),0._jprb)
!        zcrit(jl) = MIN(zcrit(jl),1._jprb)
!        !ZDRAIN(JL) = 8.e-4_JPRB - 4.e-4_JPRB*zcrit(jl) ! for PLU
!        zdrain(JL) = 2.8e4_JPRB - 1.4e4_JPRB*zcrit(jl)  ! for land
!
!        !KF
!      ENDDO
!ELSEIF (.NOT. PRESENT (paer_ss)) THEN

    IF (icpl_aero_conv == 1) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zd, zc)
      DO jl=kidia,kfdia
        zttop0(jl) = 0.5_jprb*(pten(jl,kctop0(jl))+pten(jl,kctop0(jl)-1)) ! cloud top temperature
        zdrain(jl)  = 0.25E4_JPRB + 3.75E-5_JPRB * pcloudnum(jl)  ! minimum cloud depth for generating precip: 30-150 hPa ...
        zdnoprc(jl) = 0.5E-4_JPRB + 8.75E-13_JPRB * pcloudnum(jl) ! QC autoconversion threshold: 0.065 - 0.35 g/kg ...
        !
        ! Distinction between warm clouds and mixed-phase clouds:
        ! Increase threshold for cloud depth and QC if cloud top temperature is above -19 C,
        ! particularly if the aerosol characteristics are continental
        zd = MIN(1._jprb,0.166_jprb*MAX(0._jprb,zttop0(jl)-254._jprb))     ! transition starts at about -13 C
        zc = MIN(2._jprb,1.e-8_jprb*MAX(0._jprb,pcloudnum(jl)-25.e6_jprb)) ! maximum is reached for a cloud number density of 225e6/m**3
        zdrain(jl)  = zdrain(jl)  + zc*zd*1.0e4_jprb    ! enhancement by at most 200 hPa
        zdnoprc(jl) = zdnoprc(jl) + zc*zd*1.25e-4_jprb  ! enhancement by at most 0.25 g/kg
      ENDDO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        IF(.NOT. ldland(jl) .AND. .NOT. ldlake(jl) .AND. .NOT. lgrz_deepconv) THEN
          zdrain(jl)  = MIN(0.6E4_JPRB, zdrain(jl) ) ! ... but over ocean at most 60 hPa
          zdnoprc(jl) = MIN(3.e-4_JPRB, zdnoprc(jl)) ! ... but over ocean at most 0.3 g/kg
        ENDIF
      ENDDO
    ELSE IF (ltuning_kessler) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        IF(ldland(jl) .OR. ldlake(jl)) THEN
          zdrain(jl)  = 1.0E4_JPRB  ! minimum cloud depth for generating precip: 100 hPa over land
          zdnoprc(jl) = 2.25e-4_JPRB ! autoconversion threshold: 0.225 g/kg over land
        ELSE
          zdrain(jl)  = 0.4E4_JPRB  ! minimum cloud depth for generating precip: 40 hPa over water
          zdnoprc(jl) = 1.0e-4_JPRB ! autoconversion threshold: 0.1 g/kg over water
        ENDIF
      ENDDO
    ELSE
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        IF(ldland(jl) .OR. ldlake(jl)) THEN
          zdnoprc(jl) = 5.e-4_JPRB ! 0.5 g/kg over land
        ELSE
          zdnoprc(jl) = 3.e-4_JPRB ! 0.3 g/kg over water
        ENDIF
        zdrain(jl)=-1.0E5_JPRB ! cloud depth criterion deactivated
      ENDDO
    ENDIF ! if over icpl_aero_conv == 1

!ENDIF !present aerosol
!

!----------------------------------------------------------------------

!     3.0          INITIALIZE VALUES AT cloud base LEVEL
!                  -------------------------------------

!$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb)
DO jl=kidia,kfdia
  kctop(jl)=kcbot(jl)
  IF(ldcum(jl)) THEN
    ikb=kcbot(jl)
    pkineu(jl,ikb)=0.5_JPRB*pwubase(jl)**2
    pmfu(jl,ikb)=pmfub(jl)
    pmfus(jl,ikb)=pmfub(jl)*(rcpd*ptu(jl,ikb)+pgeoh(jl,ikb))
    pmfuq(jl,ikb)=pmfub(jl)*pqu(jl,ikb)
    pmful(jl,ikb)=pmfub(jl)*plu(jl,ikb)
  ENDIF
ENDDO

!----------------------------------------------------------------------

!     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!                  -------------------------------------------------

!$ACC LOOP SEQ
DO jk=klev-1,ktdia+2,-1

!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------

!  ik=jk
!  CALL cubasmcn &
!   & ( kidia,    kfdia,    klon,    klev,&
!   & ik,&
!   & pten,     pqen,     pqsen,&
!   & pvervel,  pgeo,     pgeoh,    ldcum,    ktype,    klab,&
!   & kcbot,    pmfu,     pmfub,    plrain,&
!   & ptu,      pqu,      plu,&
!   & pmfus,    pmfuq,    pmful,    pdmfup)  

! cubasmcn is inlined for better efficiency
  kk=jk
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zzzmb)
  DO jl=kidia,kfdia

    IF(.NOT.ldcum(jl).AND.klab(jl,kk+1) == 0 ) THEN
      IF(lmfmid.AND.pgeo(jl,kk) > 5000.0_JPRB.AND.pgeo(jl,kk)<1.e5_jprb &
        & .AND.pqen(jl,kk) > 0.80_JPRB*pqsen(jl,kk)) THEN
        ptu(jl,kk+1)=(rcpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))/rcpd
        pqu(jl,kk+1)=pqen(jl,kk)
        plu(jl,kk+1)=0.0_JPRB
        zzzmb=MAX(rmfcmin,-pvervel(jl,kk)/rg)
        zzzmb=MIN(zzzmb,rmfcmax)
        pmfub(jl)=zzzmb
        pmfu(jl,kk+1)=pmfub(jl)
        pmfus(jl,kk+1)=pmfub(jl)*(rcpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
        pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
        pmful(jl,kk+1)=0.0_JPRB
        pdmfup(jl,kk+1)=0.0_JPRB
        plrain(jl,kk+1)=0.0_JPRB
        kcbot(jl)=kk
        klab(jl,kk+1)=1
        ktype(jl)=3
      ENDIF
    ENDIF
  ENDDO
! End of code inlined from cubasmcn


  is=0

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zmfmax, zfac)
  DO jl=kidia,kfdia
    llflag(jl)=.FALSE.
    zprecip(jl)=0.0_JPRB
    llo1(jl)=.FALSE.

    is=is+klab(jl,jk+1)

    IF(klab(jl,jk+1) == 0) klab(jl,jk)=0
    IF((ldcum(jl).AND.klab(jl,jk+1) == 2).OR.&
       & (ktype(jl) == 3 .AND. klab(jl,jk+1) == 1)) THEN  
      llflag(jl)=.TRUE.
    ENDIF
!   IF(klab(jl,jk+1) > 0) THEN
!     llflaguv(jl)=.TRUE.
!   ELSE
!     llflaguv(jl)=.FALSE.
!   ENDIF
    zph(jl)=paph(jl,jk)
    IF(ktype(jl) == 3.AND.jk == kcbot(jl)) THEN
      zmfmax=(paph(jl,jk)-paph(jl,jk-1))*zcons2*rmflic+rmflia
      IF(pmfub(jl) > zmfmax) THEN
        zfac=zmfmax/pmfub(jl)
        pmfu(jl,jk+1)=pmfu(jl,jk+1)*zfac
        pmfus(jl,jk+1)=pmfus(jl,jk+1)*zfac
        pmfuq(jl,jk+1)=pmfuq(jl,jk+1)*zfac
        pmfub(jl)=zmfmax
      ENDIF
    ENDIF
  ENDDO


  IF(is > 0) llo3=.TRUE.


! To avoid the reduction of is on GPUs, the variable llo3 is not set on GPUs and the
! following computations are performed on GPUs. But most of these computations are 
! enclosed by IF llfalg(jl) or IF (llo1(jl), which are .false.  if llo3 is .false.
!  (because all klab must be 0 then.)
! Only resetting of ptu, pqu, plu to environmental values if below departure level
! is not enclosed by any if-statement: what can go wrong here???


!*                  SPECIFY ENTRAINMENT RATES IN *CUENTR*
!                   -------------------------------------

!  ik=jk
!  CALL cuentr &
!   & ( kidia,    kfdia,    klon,     klev,&
!   & ik,       kcbot,&
!   & ldcum,    llo3,&
!   & paph,     pgeoh, zdgeoh,&
!   & pmfu,&
!   &  zpbase, zdmfen,   zdmfde )
  kk=jk

! Code inlined from cuentr for better efficiency

  IF(llo3) THEN

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ik)
    DO jl=kidia,kfdia
      zdmfen(jl)=0.0_JPRB
      zdmfde(jl)=0.0_JPRB
     !KF only needed for cloud depth threshold in cuascn.F90
      IF(ldcum(jl)) THEN
        ik=MAX(1,kcbot(jl))
        zpbase(jl)=paph(jl,ik)
        zptop0(jl) = paph(jl,kctop0(jl)) ! cloud top pressure
      ENDIF
    ENDDO

    !*    1.1          SPECIFY ENTRAINMENT RATES
    !                  -------------------------

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdz, zmf)
    DO jl=kidia,kfdia
      IF(ldcum(jl)) THEN
        ZDZ=(PGEOH(JL,KK)-PGEOH(JL,KK+1))*ZRG
        zmf=pmfu(jl,kk+1)*zdz
        IF(kk < kcbot(jl)) THEN
          zdmfen(jl)=0.0_JPRB*zmf
          zdmfde(jl)=detrpen*zmf
        ENDIF
      ENDIF
    ENDDO


  ENDIF

! End of code inlined from cuentr

!                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
!                  ---------------------------------------------------


  IF(llo3) THEN


    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      zqold(jl)=0.0_JPRB
    ENDDO

    !$NEC ivdep
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zmfmax, zxs, zmftest, zchange, zxe, zqeen, zseen) &
    !$ACC   PRIVATE(zleen, zscde, zqude, zmfusk, zmfuqk, zmfulk, zc)
    DO jl=kidia,kfdia
      IF (llflag(jl)) THEN
        zdmfde(jl)=MIN(zdmfde(jl),0.75_JPRB*pmfu(jl,jk+1))
        IF(jk==kcbot(jl)) THEN  ! bugfix 2014-10-06; was kcbot(jl)-1 before
          zoentr(jl)=-zentrorg(jl)*(MIN(1.0_JPRB,pqen(jl,jk)/pqsen(jl,jk))-1.0_JPRB)*&
          &(PGEOH(JL,JK)-PGEOH(JL,JK+1))*ZRG
          zoentr(jl)=MIN(0.4_JPRB,zoentr(jl))*pmfu(jl,jk+1)
        ENDIF
        IF(jk < kcbot(jl)) THEN
          zmfmax=(paph(jl,jk)-paph(jl,jk-1))*zcons2*rmflic+rmflia
          IF (ktype(jl)==2) zmfmax=zmfmax*0.5_jprb
          zxs=MAX(pmfu(jl,jk+1)-zmfmax,0.0_JPRB)
          pwmean(jl)=pwmean(jl)+pkineu(jl,jk+1)*(pap(jl,jk+1)-pap(jl,jk))
          zdpmean(jl)=zdpmean(jl)+pap(jl,jk+1)-pap(jl,jk)
          !CURRENT LEVEL'S ENTRAINMENT IS EQUAL TO ZOENTR VALUE OF LEVEL BELOW - ITERATIVE
          zdmfen(jl)=zoentr(jl)
          !OVERWRITE TURB DETR FOR DEEP CONV FROM CALL TO CUENTR WITH SHALLOW CU ENTRAINMENT
          IF(ktype(jl)>=2)THEN
            !FOR SHALLOW CONVECTION, DOUBLE THE ENTRAINMENT CALCULATED LATER FOR LEVEL BELOW
             zdmfen(jl)=ENTSHALP*zdmfen(jl)! Factor 2 on entrainment for shallow conv/ktype 2(initial ascent)
             !FOR SHALLOW CONVECTION, SET TURB DETRAINEMENT EQUAL TO ENTRAINMENT
             zdmfde(jl)=zdmfen(jl)
          ENDIF
          zc = 1.6_JPRB-MIN(1.0_JPRB,PQEN(JL,JK)/PQSEN(JL,JK))
          !DETRAINMENT out of cuentr (turb)*(1.6-RH) EQN 6.9
          ZDMFDE(JL)=ZDMFDE(JL)*MAX(0.9_jprb*zc,zc**4)
          ! safety check ???
          zmftest=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
          zchange=MAX(zmftest-zmfmax,0.0_JPRB)
          zxe=MAX(zchange-zxs,0.0_JPRB)
          zdmfen(jl)=zdmfen(jl)-zxe
          zchange=zchange-zxe
          zdmfde(jl)=zdmfde(jl)+zchange
        ENDIF
        ! extra diagnostic - net (entrainment-detrainment), or change of MF with height
        pdmfen(jl,jk) = zdmfen(jl)-zdmfde(jl)
        ! THIS LEVEL'S MASSFLUX CALCULATED HERE FROM LEVEL BELOW, PLUS ENTR, MINUS DETR
        pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
        zqeen=pqenh(jl,jk+1)*zdmfen(jl)
        zseen=(rcpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
        IF(plitot(jl,jk) > rlmin) THEN
          zleen=plitot(jl,jk)*zdmfen(jl)
        ELSE
          zleen=0.0_JPRB
        ENDIF
        zscde=(rcpd*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
        zqude=pqu(jl,jk+1)*zdmfde(jl)
        plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
        zmfusk=pmfus(jl,jk+1)+zseen-zscde
        zmfuqk=pmfuq(jl,jk+1)+zqeen-zqude
        zmfulk=pmful(jl,jk+1)+zleen-plude(jl,jk)
        plu(jl,jk)=zmfulk*(1.0_JPRB/MAX(rmfcmin,pmfu(jl,jk)))
        pqu(jl,jk)=zmfuqk*(1.0_JPRB/MAX(rmfcmin,pmfu(jl,jk)))
        ptu(jl,jk)=(zmfusk*(1.0_JPRB/MAX(rmfcmin,pmfu(jl,jk)))-&
         & pgeoh(jl,jk))/rcpd  
        ptu(jl,jk)=MAX(100._jprb,ptu(jl,jk))
        ptu(jl,jk)=MIN(400._jprb,ptu(jl,jk))
        zqold(jl)=pqu(jl,jk)
        plrain(jl,jk)=plrain(jl,jk+1)*(pmfu(jl,jk+1)-zdmfde(jl))*&
         & (1.0_JPRB/MAX(rmfcmin,pmfu(jl,jk)))  
        zluold(jl)=plu(jl,jk)
      ENDIF
    ENDDO

    ! reset to environmental values if below departure level
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      IF ( jk > kdpl(jl) ) THEN
        ptu(jl,jk)=ptenh(jl,jk)      
        pqu(jl,jk)=pqenh(jl,jk)      
        plu(jl,jk)=0.0_JPRB
        zluold(jl)=plu(jl,jk)
      ENDIF
    ENDDO

!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------

    ik=jk
    icall=1
    CALL cuadjtq &
     & ( kidia,    kfdia,    klon,     klev,    ik,&
     &   zph,      ptu,      pqu,      llflag,  icall)

    IF (lphylin) THEN
      ! lphylin=.TRUE. is not supported on GPUs
      !DIR$ IVDEP
      !$NEC ivdep
      DO jl=kidia,kfdia
        IF (llflag(jl)) THEN
          IF(pqu(jl,jk) /= zqold(jl)) THEN
            zoealfa   = 0.545_JPRB*(TANH(0.17_JPRB*(ptu(jl,jk  )-rlptrc))+1.0_JPRB)
            zoealfap  = 0.545_JPRB*(TANH(0.17_JPRB*(ptu(jl,jk+1)-rlptrc))+1.0_JPRB)
            plglac(jl,jk)=plu(jl,jk)*((1.0_JPRB-zoealfa)-(1.0_JPRB-zoealfap))
            ! add glaciation of rain
            ZFAC      = 0.545_JPRB*(TANH(0.17_JPRB*(PTEN(JL,JK  )-RLPTRC))+1.0_JPRB)
            PLGLAC(JL,JK)=PLGLAC(JL,JK)+ZFAC*PDMFUP(JL,JK+1)/MAX(RMFCMIN,PMFU(JL,JK+1))*&
                         &(0.5_JPRB+SIGN(0.5_JPRB,RTT-PTEN(JL,JK)))*zglac
            ptu(jl,jk)=ptu(jl,jk)+ralfdcp*plglac(jl,jk)
          ENDIF
        ENDIF
      ENDDO

    ELSE

      !DIR$ IVDEP
      !$NEC ivdep
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zfac)
      DO jl=kidia,kfdia
        IF (llflag(jl)) THEN
          IF(pqu(jl,jk) /= zqold(jl)) THEN
            plglac(jl,jk)=plu(jl,jk)*((1.0_JPRB-foealfcu(ptu(jl,jk)))-&
             & (1.0_JPRB-foealfcu(ptu(jl,jk+1))))  
            ! add glaciation of rain, only fraction added to updraught heat
            ZFAC=FOEALFCU(PTEN(JL,JK))
            PLGLAC(JL,JK)=PLGLAC(JL,JK)+ZFAC*PDMFUP(JL,JK+1)/MAX(RMFCMIN,PMFU(JL,JK+1))*&
                         &(0.5_JPRB+SIGN(0.5_JPRB,RTT-PTEN(JL,JK)))*zglac
            ptu(jl,jk)=ptu(jl,jk)+ralfdcp*plglac(jl,jk)
          ENDIF
        ENDIF
      ENDDO

    ENDIF

    !DIR$ IVDEP
    !$NEC sparse
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zbe, zbuoc, zdkbuo, zdken, zkedke, zmfun, zentr_prof, ikb)
    DO jl=kidia,kfdia
      IF (llflag(jl)) THEN
        IF(pqu(jl,jk) /= zqold(jl)) THEN
          klab(jl,jk)=2
          plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
          zbc(jl)=ptu(jl,jk)*(1.0_JPRB+retv*pqu(jl,jk)-plu(jl,jk+1)-plrain(jl,jk+1))
          zbe=ptenh(jl,jk)*(1.0_JPRB+retv*pqenh(jl,jk))
          zbuo(jl,jk)=zbc(jl)-zbe

! set flags in case of midlevel convection

          IF(ktype(jl) == 3 .AND. klab(jl,jk+1)== 1) THEN
            IF(zbuo(jl,jk) > -0.5_JPRB) THEN
              ldcum(jl)=.TRUE.
              kctop(jl)=jk
              pkineu(jl,jk)=0.5_JPRB
            ELSE
              klab(jl,jk)=0
              pmfu(jl,jk)=0.0_JPRB
              plude(jl,jk)=0.0_JPRB
              plu(jl,jk)=0.0_JPRB
            ENDIF
          ENDIF

          IF(klab(jl,jk+1) == 2) THEN! if layer below is inside cloud

            IF(zbuo(jl,jk) < 0.0_JPRB )THEN !.AND.klab(jl,jk+1) == 2) THEN
              ptenh(jl,jk)=0.5_JPRB*(pten(jl,jk)+pten(jl,jk-1))
              pqenh(jl,jk)=0.5_JPRB*(pqen(jl,jk)+pqen(jl,jk-1))
              zbuo(jl,jk)=zbc(jl)-ptenh(jl,jk)*(1.0_JPRB+retv*pqenh(jl,jk))
            ENDIF
            zbuoc=(zbuo(jl,jk)/(ptenh(jl,jk)*(1.0_JPRB+retv*pqenh(jl,jk)))&
             & +zbuo(jl,jk+1)/(ptenh(jl,jk+1)*(1.0_JPRB+retv*&
             & pqenh(jl,jk+1))))*0.5_JPRB  
            zdkbuo=(pgeoh(jl,jk)-pgeoh(jl,jk+1))*zfacbuo*zbuoc

! either use entrainment rate or if zero
! use detrainmnet rate as a subsitute for 
! mixing and "pressure" gradient term in upper
! troposphere

            IF(zdmfen(jl) > 0.0_JPRB)THEN
              zdken=MIN(1.0_JPRB,(1.0_JPRB + z_cwdrag)*&
               & zdmfen(jl)/MAX(rmfcmin,pmfu(jl,jk+1)))  
            ELSE
              zdken=MIN(1.0_JPRB,(1.0_JPRB + z_cwdrag)*&
               & zdmfde(jl)/MAX(rmfcmin,pmfu(jl,jk+1)))  
            ENDIF

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
            ! Alternative organised detrainment calcuation based on cloud area 
            ! profile.
            ! By what fraction is core area reduced from layer below?
            IF (lstochdetr) THEN              
              IF (pcore(jl,jk) .lt. pcore(jl,jk+1)) THEN
                 zmfun=(pcore(jl,jk+1)-pcore(jl,jk))/pcore(jl,jk+1)*.5_JPRB
                 ! Detrainment is corresponding fraction of MF in layer below,
                 ! if this is greater than the estimate previously calculated
                 zdmfde(jl)=MAX(zdmfde(jl),pmfu(jl,jk+1)*zmfun)
                 plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
                 pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)! Mass flux is same as layer below, minus detraiment
              ENDIF
            ENDIF 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

            pkineu(jl,jk)=(pkineu(jl,jk+1)*(1.0_JPRB-zdken)+zdkbuo)/(1.0_JPRB+zdken)
            ! CALCULATE ORGANISED DETRAINMENT WHEN UPDRAFT NEGATIVELY BUOYANT
            IF(zbuo(jl,jk) < 0.0_JPRB ) THEN ! .AND.klab(jl,jk+1) == 2) THEN
              ! Is stable -> no org. entrainment. This overwrites PMFU for current level
              ! which has been calculated above with organised entrainment
              zkedke=pkineu(jl,jk)/MAX(1.e-10_JPRB,pkineu(jl,jk+1))
              zkedke=MAX(1.e-30_JPRB,MIN(1.0_JPRB,zkedke))
              zmfun=SQRT(zkedke)
! ** suggestion by P. Bechtold (2013-11-21) - but degrades various scores in ICON **
!             zmfun = (1.6_JPRB-MIN(1.0_JPRB,pqen(JL,JK)/pqsen(JL,JK)))*zmfun
              zdmfde(jl)=MAX(zdmfde(jl),pmfu(jl,jk+1)*(1.0_JPRB-zmfun))
              plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
              pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)! Mass flux is same as layer below, minus detraiment
            ENDIF

            !SWITCHES *ON* CALCULATION OF ENTRAINMENT ONLY IN BUOYANT UPDRAFTS
            !ZOENTR IS OVERWRITTEN, BUT NOT USED UNTIL THE NEXT JK LEVEL IN THE LOOP
            IF(zbuo(jl,jk) > -0.2_JPRB) THEN !.AND.klab(jl,jk+1) == 2) THEN
              !when positively buoyant, have organised entrainment
              !which increases MF with height, while detrainment is pretty small and constant
              ikb=kcbot(jl)

              ! DOCUMENTATION EQN 6.7 - ENTRAINMENT TERM REPRESENTING BOTH ORGANISED AND TURBULENT
              ! Old entrainment profile shape
              !zoentr(jl)=zentrorg(jl)*(1.3_JPRB-MIN(1.0_JPRB,pqen(jl,jk-1)/pqsen(jl,jk-1)))*&
              !  &(pgeoh(jl,jk-1)-pgeoh(jl,jk))*zrg*MIN(1.0_JPRB,pqsen(jl,jk)/pqsen(jl,ikb))**3
              ! last term **3 is vertical scaling function, supposed to mimic the effects of a cloud ensemble and/
              ! or the effect of the inreasing Rup (cloud radius) with height. It strongly decreases with height,
              ! so that detrainment will eventually become greater than entrainment.
              !*********************************************************************************
              !Guenthers new entrainment profile shape 
              zentr_prof = MERGE((pqsen(jl,jk)/pqsen(jl,ikb))**2, (pqsen(jl,jk)/pqsen(jl,ikb))**3, lgrz_deepconv)
              zoentr(jl)=zentrorg(jl)*(1.3_JPRB-MIN(1.0_JPRB,pqen(jl,jk-1)/pqsen(jl,jk-1)))*&
                &(pgeoh(jl,jk-1)-pgeoh(jl,jk))*zrg*MIN(1.0_JPRB,zentr_prof)
              zoentr(jl)=MIN(0.4_JPRB,zoentr(jl))*pmfu(jl,jk)
            ELSE
              zoentr(jl)=0.0_JPRB
            ENDIF

            ! Erase values if below departure level
            IF ( jk > kdpl(jl) ) THEN
              pmfu(jl,jk)=pmfu(jl,jk+1)
              pkineu(jl,jk)=0.5_JPRB
            ENDIF
            ! determine convection top level;
            ! the last set of criteria serves to limit the overshooting of updrafts through the tropopause
            IF (pkineu(jl,jk) > 0.0_JPRB .AND. pmfu(jl,jk) > 0.0_JPRB .AND. (zbuo(jl,jk) > -2._JPRB .OR. &
               (pten(jl,jk-1)-pten(jl,jk))/(zrg*(pgeo(jl,jk-1)-pgeo(jl,jk))) < -3.e-3_jprb ) ) THEN
              kctop(jl)=jk
              llo1(jl)=.TRUE.
            ELSE
              klab(jl,jk)=0
              pmfu(jl,jk)=0.0_JPRB
              pkineu(jl,jk)=0.0_JPRB
              zdmfde(jl)=pmfu(jl,jk+1)
              plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
            ENDIF
          
! store detrainment rates for updraught

            IF ( pmfu(jl,jk+1) > 0.0_JPRB ) THEN
              pmfude_rate(jl,jk)=zdmfde(jl)
            ENDIF
          
          ENDIF ! klab=2
        ENDIF ! zqold
      ENDIF
    ENDDO

    !$NEC ivdep
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      IF (llflag(jl)) THEN
!       ELSEIF(LLFLAG(JL).AND.KTYPE(JL)==2.AND.PQU(JL,JK) == ZQOLD(JL)) THEN
!       ELSEIF(ktype(jl)==2.AND.pqu(jl,jk) == zqold(jl)) THEN
          IF(ktype(jl)==2.AND.pqu(jl,jk) == zqold(jl)) THEN
            klab(jl,jk)=0
            pmfu(jl,jk)=0.0_JPRB
            pkineu(jl,jk)=0.0_JPRB
            zdmfde(jl)=pmfu(jl,jk+1)
            plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
            pmfude_rate(jl,jk)=zdmfde(jl)
          ENDIF
      ENDIF
    ENDDO

!              CALCULATE PRECIPITATION RATE BY
!              ANALYTIC INTEGRATION OF EQUATION FOR L

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zwu, zprcon, zdt, zcbf, zzco, zlcrit, zdfi, zc, zd, zint, zlnew)
    DO jl=kidia,kfdia
      IF(llo1(jl)) THEN

        IF (plu(jl,jk) > zdnoprc(jl) .AND. (zpbase(jl)-MIN(zptop0(jl),paph(jl,jk))) > zdrain(jl)) THEN

          zwu=MIN(15._jprb,SQRT(2.0_JPRB*MAX(0.1_JPRB,pkineu(jl,jk+1))))
          zprcon=zprcdgw/(0.75_JPRB*zwu)

!           PARAMETERS FOR BERGERON-FINDEISEN PROCESS (T < -5C)

          zdt=MIN(rtbercu-rticecu,MAX(rtber-ptu(jl,jk),0.0_JPRB))
          zcbf=1.0_JPRB+z_cprc2*SQRT(zdt)
          zzco=zprcon*zcbf
          zlcrit=zdnoprc(jl)/zcbf

          zdfi=pgeoh(jl,jk)-pgeoh(jl,jk+1)
          zc=(plu(jl,jk)-zluold(jl))
          zd=zzco*(1.0_JPRB-EXP(-(plu(jl,jk)/zlcrit)**2))*zdfi
          zint=EXP(-zd)
          zlnew=zluold(jl)*zint+zc/zd*(1.0_JPRB-zint)
          zlnew=MAX(0.0_JPRB,MIN(plu(jl,jk),zlnew))
          zlnew=MIN(z_cldmax,zlnew)
          zprecip(jl)=MAX(0.0_JPRB,zluold(jl)+zc-zlnew)
          pdmfup(jl,jk)=zprecip(jl)*pmfu(jl,jk)
          plrain(jl,jk)=plrain(jl,jk)+zprecip(jl)
          plu(jl,jk)=zlnew
        ENDIF
      ENDIF
    ENDDO

    IF (lphylin) THEN
      ! lphylin=.TRUE. is not supported on GPUs
      DO jl=kidia,kfdia
        IF(llo1(jl)) THEN
          IF(plrain(jl,jk) > 0.0_JPRB) THEN
            zvw=21.18_JPRB*EXP(0.2_JPRB*LOG(plrain(jl,jk)))  ! optimization of plrain(JL,JK)**0.2_JPRB
            zvi=z_cwifrac*zvw
            zalfaw=0.545_JPRB*(TANH(0.17_JPRB*(ptu(jl,jk)-rlptrc))+1.0_JPRB)
            zvv=zalfaw*zvw+(1.0_JPRB-zalfaw)*zvi
            zrold=plrain(jl,jk)-zprecip(jl)
            zc=zprecip(jl)
            zwu=MIN(15._jprb,SQRT(2.0_JPRB*MAX(0.1_JPRB,pkineu(jl,jk))))
            zd=zvv/zwu
            zint=EXP(-zd)
            zrnew=zrold*zint+zc/zd*(1.0_JPRB-zint)
            zrnew=MAX(0.0_JPRB,MIN(plrain(jl,jk),zrnew))
            plrain(jl,jk)=zrnew
          ENDIF
        ENDIF
      ENDDO

    ELSE

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zvw, zvi, zalfaw, zvv, zrold, zc, zwu, zd, zint, zrnew)
      DO jl=kidia,kfdia
        IF(llo1(jl)) THEN
          IF(plrain(jl,jk) > 0.0_JPRB) THEN
            zvw=21.18_JPRB*EXP(0.2_JPRB*LOG(plrain(jl,jk)))
            zvi=z_cwifrac*zvw
            zalfaw=foealfcu(ptu(jl,jk))
            zvv=zalfaw*zvw+(1.0_JPRB-zalfaw)*zvi
            zrold=plrain(jl,jk)-zprecip(jl)
            zc=zprecip(jl)
            zwu=MIN(15._jprb,SQRT(2.0_JPRB*MAX(0.1_JPRB,pkineu(jl,jk))))
            zd=zvv/zwu
            zint=EXP(-zd)
            zrnew=zrold*zint+zc/zd*(1.0_JPRB-zint)
            zrnew=MAX(0.0_JPRB,MIN(plrain(jl,jk),zrnew))
            plrain(jl,jk)=zrnew
          ENDIF
        ENDIF
      ENDDO

    ENDIF

    !$NEC ivdep
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      IF (llflag(jl)) THEN
        pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
        pmfus(jl,jk)=(rcpd*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
      ENDIF
    ENDDO


  ENDIF  ! llo3

ENDDO    ! jk=klev-1,ktdia+2,-1

!----------------------------------------------------------------------

!     5.           FINAL CALCULATIONS 
!                  ------------------
!$NEC sparse
!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  IF(kctop(jl) == -1) ldcum(jl)=.FALSE.
  kcbot(jl)=MAX(kcbot(jl),kctop(jl))
  IF(ldcum(jl)) THEN
    pwmean(jl)=MAX(1.e-2_JPRB,pwmean(jl)/MAX(1.0_JPRB,zdpmean(jl)))
    pwmean(jl)=SQRT(2.0_JPRB*pwmean(jl))
  ENDIF
ENDDO

!$ACC END PARALLEL

!$ACC WAIT(1)
!$ACC END DATA

IF (lhook) CALL dr_hook('CUASCN',1,zhook_handle)
END SUBROUTINE cuascn


! This is inlined above for better efficiency
  SUBROUTINE cubasmcn &
    & (kidia,    kfdia,    klon,    klev,&
    & kk, lmfmid,                        &
    & pten,     pqen,     pqsen,    &
    & pvervel,  pgeo,     pgeoh,    ldcum,    ktype,    klab,&
    & kcbot,    pmfu,     pmfub,    plrain,&
    & ptu,      pqu,      plu,&
    & pmfus,    pmfuq,    pmful,    pdmfup )
    !>
    !! Description:
    !!          M.TIEDTKE         E.C.M.W.F.     12/89

    !!          PURPOSE.
    !!          --------
    !!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
    !!          FOR MIDLEVEL CONVECTION

    !!          INTERFACE
    !!          ---------

    !!          THIS ROUTINE IS CALLED FROM *CUASC*.
    !!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
    !!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION

    !!          METHOD.
    !!          --------
    !!          S. TIEDTKE (1989)

    !!
    !!     PARAMETER     DESCRIPTION                                   UNITS
    !!     ---------     -----------                                   -----
    !!     INPUT PARAMETERS (INTEGER):

    !!    *KIDIA*        START POINT
    !!    *KFDIA*        END POINT
    !!    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!    *KLEV*         NUMBER OF LEVELS
    !!    *KK*           ACTUAL LEVEL

    !!    INPUT PARAMETERS (REAL):

    !!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
    !!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
    !!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
    !!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
    !!    *PGEO*         GEOPOTENTIAL                                  M2/S2
    !!    *zdPGEO*       GEOPOTENTIAL thickness at full levels         M2/S2
    !!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
    !!    *PLRAIN*       RAIN WATER CONTENT IN UPDRAFTS                KG/KG

    !!    INPUT PARAMETERS (LOGICAL):

    !!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

    !!    UPDATED PARAMETERS (INTEGER):
    
    !!    *KTYPE*        TYPE OF CONVECTION
    !!                       1 = PENETRATIVE CONVECTION
    !!                       2 = SHALLOW CONVECTION
    !!                       3 = MIDLEVEL CONVECTION
    !!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
    !!                        KLAB=2 FOR CLOUD LEVELS
    !!    *KCBOT*        CLOUD BASE LEVEL
    
    !!    OUTPUT PARAMETERS (REAL):
    
    !!    *PMFU*         MASSFLUX IN UPDRAFTS                          KG/(M2*S)
    !!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE            KG/(M2*S)
    !!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
    !!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
    !!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
    !!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
    !!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
    !!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
    !!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS        KG/(M2*S)
    
    !!          EXTERNALS
    !!          ---------
    !!          NONE
    
    !!----------------------------------------------------------------------
    
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    !USE YOMCST   , ONLY : RG       ,RCPD
    !USE YOECUMF  , ONLY : ENTRMID  ,RMFCMAX  ,RMFCMIN  ,LMFMID

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: kk
    LOGICAL           ,INTENT(in)    :: lmfmid
    REAL(KIND=jprb)   ,INTENT(in)    :: pten(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqsen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pvervel(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1)
    LOGICAL           ,INTENT(in)    :: ldcum(klon)
    INTEGER(KIND=jpim),INTENT(inout) :: ktype(klon)
    INTEGER(KIND=jpim),INTENT(inout) :: klab(klon,klev)
    INTEGER(KIND=jpim),INTENT(inout) :: kcbot(klon)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pmfub(klon)
    REAL(KIND=jprb)   ,INTENT(inout) :: plrain(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: pqu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: plu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pmfus(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pmfuq(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pmful(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout)   :: pdmfup(klon,klev)
    INTEGER(KIND=jpim) :: jl

    !
    REAL(KIND=jprb) :: zzzmb
    REAL(KIND=jprb) :: zhook_handle

    !----------------------------------------------------------------------

    !!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    !!                  -------------------------------------------

!OCL NOVREC
    IF (lhook) CALL dr_hook('CUBASMCN',0,zhook_handle)
!DIR$ IVDEP
    DO jl=kidia,kfdia

      IF(.NOT.ldcum(jl).AND.klab(jl,kk+1) == 0 ) THEN
        IF(lmfmid.AND.pgeo(jl,kk) > 5000.0_JPRB.AND.pgeo(jl,kk)<1.e5_jprb &
        !! IF(LMFMID.AND.PGEO(JL,KK) > 4000.0_JPRB.AND.PGEO(JL,KK)<1.E5_JPRB &
          & .AND.pqen(jl,kk) > 0.80_JPRB*pqsen(jl,kk)) THEN
          !KF
          ptu(jl,kk+1)=(rcpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))/rcpd
          pqu(jl,kk+1)=pqen(jl,kk)
          plu(jl,kk+1)=0.0_JPRB
          zzzmb=MAX(rmfcmin,-pvervel(jl,kk)/rg)
          zzzmb=MIN(zzzmb,rmfcmax)
          pmfub(jl)=zzzmb
          pmfu(jl,kk+1)=pmfub(jl)
          pmfus(jl,kk+1)=pmfub(jl)*(rcpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
          pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
          pmful(jl,kk+1)=0.0_JPRB
          pdmfup(jl,kk+1)=0.0_JPRB
          plrain(jl,kk+1)=0.0_JPRB
          kcbot(jl)=kk
          klab(jl,kk+1)=1
          ktype(jl)=3
        ENDIF
      ENDIF
    ENDDO

    IF (lhook) CALL dr_hook('CUBASMCN',1,zhook_handle)
  END SUBROUTINE cubasmcn

! This is inlined above for better efficiency
  SUBROUTINE cuentr &
    & ( kidia,    kfdia,    klon,     klev,&
    & kk,       kcbot,&
    & ldcum,    ldwork, detrpen, &
    & paph,     pgeoh,  zdgeoh,      &
    & pmfu,&
    & pcbase,   pdmfen,   pdmfde )
    !>
    !! Description:
    !!          M.TIEDTKE         E.C.M.W.F.     12/89
    !!          P.BECHTOLD        E.C.M.W.F.     06/07

    !!          PURPOSE.
    !!          --------
    !!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
    !!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION

    !!          INTERFACE
    !!          ---------

    !!          THIS ROUTINE IS CALLED FROM *CUASC*.
    !!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
    !!          AND UPDRAFT VALUES T,Q ETC
    !!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES

    !!          METHOD.
    !!          --------
    !!          TURBULENT ENTRAINMENT IS SIMULATED BY A CONSTANT
    !!          MULTIPLIED BY A VERTICAL SCALING FUNCTION
    !!
    !!     PARAMETER     DESCRIPTION                                   UNITS
    !!     ---------     -----------                                   -----
    !!     INPUT PARAMETERS (INTEGER):

    !!    *KIDIA*        START POINT
    !!    *KFDIA*        END POINT
    !!    *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!    *KLEV*         NUMBER OF LEVELS
    !!    *KK*           CURRENT LEVEL
    !!    *KCBOT*        CLOUD BASE LEVEL

    !!    INPUT PARAMETERS (LOGICAL):

    !!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

    !!    INPUT PARAMETERS (REAL):

    !!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS          PA
    !!    *PGEOH*        PROVISIONAL GEOPOTENTIAL ON HALF LEVELS      m2/s2
    !!    *zdgeoh*       geopot thickness on full levels              M2/S2
    !!    *PMFU*         MASSFLUX IN UPDRAFTS                        KG/(M2*S)

    !    OUTPUT PARAMETERS (REAL):

    !!    *PCBASE*       PRESSURE AT CLOUD BASE                       PA
    !!    *PDMFEN*       ENTRAINMENT RATE                            KG/(M2*S)
    !!    *PDMFDE*       DETRAINMENT RATE                            KG/(M2*S)

    !          EXTERNALS
    !          ---------
    !          NONE

    !!----------------------------------------------------------------------

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    !USE YOMCST   , ONLY : RG
    !USE YOECUMF  , ONLY : DETRPEN

    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: kk
    INTEGER(KIND=jpim),INTENT(in)    :: kcbot(klon)
    LOGICAL           ,INTENT(in)    :: ldcum(klon)
    LOGICAL           ,INTENT(in)    :: ldwork
    REAL(KIND=jprb)   ,INTENT(in)    :: detrpen
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: zdgeoh(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pmfu(klon,klev)
    REAL(KIND=jprb)   ,INTENT(out)   :: pcbase(klon)
    REAL(KIND=jprb)   ,INTENT(out)   :: pdmfen(klon)
    REAL(KIND=jprb)   ,INTENT(out)   :: pdmfde(klon)

    LOGICAL ::  llo1

    INTEGER(KIND=jpim) :: jl, ik1

    REAL(KIND=jprb) :: zdz, zentr, zrg ! zentr2
    REAL(KIND=jprb) :: zhook_handle

    !----------------------------------------------------------------------

    !*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    !                  -------------------------------------------

    IF (lhook) CALL dr_hook('CUENTR',0,zhook_handle)
    IF(ldwork) THEN

      zrg=1.0_JPRB/rg
      DO jl=kidia,kfdia
        pdmfen(jl)=0.0_JPRB
        pdmfde(jl)=0.0_JPRB
        !KF only needed for cloud depth threshold in cuascn.F90
        IF(ldcum(jl)) THEN
          ik1=MAX(1,kcbot(jl))
          pcbase(jl)=paph(jl,ik1)
        ENDIF
      ENDDO

      !*    1.1          SPECIFY ENTRAINMENT RATES
      !!                  -------------------------

      DO jl=kidia,kfdia
        IF(ldcum(jl)) THEN
          !>KF
           ZDZ=(PGEOH(JL,KK)-PGEOH(JL,KK+1))*ZRG
          !zdz=zdgeoh(jl,kk+1)*zrg
          !<KF
          zentr=pmfu(jl,kk+1)*zdz
          llo1=kk < kcbot(jl)
          IF(llo1) THEN
            pdmfen(jl)=0.0_JPRB*zentr
            pdmfde(jl)=detrpen*zentr
          ENDIF
        ENDIF
      ENDDO

    ENDIF

    IF (lhook) CALL dr_hook('CUENTR',1,zhook_handle)
  END SUBROUTINE cuentr

  !=======================================================================

  !! KF NOTE: THIS ROUTINE iS NOT USED WITHIN THIS CODE BUT IS MAINTAINED HERE
  !!     FOR COMPLETENESS

! and this is never used
  SUBROUTINE custrat &
    & (  kidia,    kfdia,    klon,   ktdia,  klev,&
    & ldcum,    ptsphy,&
    & pap,      paph,     pgeo,&
    & pten,     pqen,     pqsat,    penth,&
    & ptent,    ptenq                              )
    !>
    !! Description:
    !!**** *CUSTRAT* - COMPUTES T,Q TENDENCIES FOR STRATOCUMULUS
    !!                 CONVECTION

    !!    M.TIEDTKE      E.C.M.W.F.    4/89 MODIF. 12/89

    !!    PURPOSE.
    !!    --------

    !!          THIS ROUTINE DOES THE PARAMETERIZATION OF BOUNDARY-LAYER
    !!    MIXING BY ENHANCED VERTICAL DIFFUSION OF SENSIBLE HEAT
    !!    AND MOISTURE FOR THE CASE OF STRATOCUMULUS CONVECTION.
    !!    THE SCHEME IS ONLY APPLIED IN THE BOUNDARY-LAYER AND ONLY
    !!    WHEN NEITHER PENETRATIVE NOR SHALLOW CONVECTION ARE ACTIVATED.

    !**   INTERFACE.
    !!    ----------

    !!          THIS ROUTINE IS CALLED FROM *CUCALL*:
    !!    IT TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
    !!    T,Q AT (T-1) AS WELL AS THE PROVISIONAL T,Q TENDENCIES AND
    !!    RETURNS ITS OUTPUT TO THE SAME SPACE:
    !!      MODIFIED TENDENCIES OF T AND Q

    !!    METHOD.
    !!    -------

    !!          ENHANCED VERTICAL DIFFUSION OF MOISTURE AND SENSIBLE
    !!    HEAT OCCURS, WHENEVER
    !!       1. LIFTED SURFACE AIR IS BUOYANT AND
    !!       2. CONDENSATION LEVEL EXISTS FOR FREE CONVECTION
    !!    THEN THE EXCHANGE COEFFICIENT IS AS FOLLOWS;
    !!        K=C1 IN CLOUD LAYER
    !!        K=C1*F(RH) AT CLOUD TOP (TOP ENTRAINMENT)

    !!    THE MATRIX INVERSION IS PERFORMED ANALOGOUSLY TO ROUTINE *VDIFF*

    !!    PARAMETER     DESCRIPTION                                   UNITS
    !!    ---------     -----------                                   -----
    !!    INPUT PARAMETERS (INTEGER):

    !!   *KIDIA*        START POINT
    !!   *KFDIA*        END POINT
    !!   *KLON*         NUMBER OF GRID POINTS PER PACKET
    !!   *KLEV*         NUMBER OF LEVELS

    !!   INPUT PARAMETERS (LOGICAL):

    !!   *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

    !C     INPUT PARAMETERS (REAL)

    !!   *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
    !!   *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS            PA
    !!   *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
    !!   *PGEO*         GEOPOTENTIAL                                  M2/S2
    !!   *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
    !!   *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
    !!   *PQSAT*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
    !!   *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)

    !!   UPDATED PARAMETERS (REAL):

    !!   *PTENT*        TEMPERATURE TENDENCY                           K/S
    !!   *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)

    !!    EXTERNALS.
    !!    ----------

    !!         *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT

    !!    Modifications.
    !!    --------------
    !!    G. Mozdzynski 2000-11-29: Corrections required for reproducibility
    !!       M.Hamrud      01-Oct-2003 CY28 Cleaning

    !----------------------------------------------------------------------

    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

    !USE YOMCST   , ONLY : RG       ,RD       ,RCPD     ,RETV
    !USE YOEVDF   , ONLY : RVDIFTS


    IMPLICIT NONE

    INTEGER(KIND=jpim),INTENT(in)    :: klon
    INTEGER(KIND=jpim),INTENT(in)    :: klev
    INTEGER(KIND=jpim),INTENT(in)    :: kidia
    INTEGER(KIND=jpim),INTENT(in)    :: kfdia
    INTEGER(KIND=jpim),INTENT(in)    :: ktdia
    LOGICAL ,INTENT(in)    :: ldcum(klon)
    REAL(KIND=jprb)   ,INTENT(in)    :: ptsphy
    REAL(KIND=jprb)   ,INTENT(in)    :: pap(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
    REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pten(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqen(klon,klev)
    REAL(KIND=jprb)   ,INTENT(in)    :: pqsat(klon,klev)
    REAL(KIND=jprb)   ,INTENT(out)   :: penth(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptent(klon,klev)
    REAL(KIND=jprb)   ,INTENT(inout) :: ptenq(klon,klev)
    REAL(KIND=jprb) ::     ztc(klon,klev),         zqc(klon,klev),&
      & zcf(klon,klev),         zcptgz(klon,klev),&
      & ztdif(klon,klev),       zqdif(klon,klev),&
      & zcpts(klon),            zqs(klon),&
      & zebs(klon,klev),        ztcoe(klon),&
      & zap(klon,klev+1),       zqold(klon)
    REAL(KIND=jprb) ::     zpp(klon)
    INTEGER(KIND=jpim) ::  ilab(klon,klev)
    LOGICAL ::  llflag(klon),           llo2(klon),            llbl(klon)
    INTEGER(KIND=jpim) :: icall, ik, ilevh, jk, jl

    REAL(KIND=jprb) :: zbuo, zcons1, zcons2, zcons3, zdisc,&
      & zdqdt, zdtdt, zfac, zkdiff1, zkdiff2, zqdp, ztmst, ztpfac1, ztpfac2
    REAL(KIND=jprb) :: zhook_handle

    !#include "cuadjtq.intfb.h"

    !-----------------------------------------------------------------------

    !!    2.           PHYSICAL CONSTANTS AND PARAMETERS.
    !!                 ---------------------------------

    IF (lhook) CALL dr_hook('CUSTRAT',0,zhook_handle)
    ztpfac1=rvdifts
    ztpfac2=1.0_JPRB/ztpfac1
    zkdiff1=10._jprb
    zkdiff2=2.5_JPRB
    ztmst=ptsphy
    zcons1=ztpfac1*ztmst*rg**2/(0.5_JPRB*rd)
    zcons2=1.0_JPRB/ztmst
    zcons3=ztmst*rcpd
    ilevh=klev/2

    !----------------------------------------------------------------------

    !*    3.           PRELIMINARY COMPUTATIONS.
    !!                 ------------------------

    DO jk=ktdia,klev
      DO jl=kidia,kfdia
        zcptgz(jl,jk)=pgeo(jl,jk)+pten(jl,jk)*rcpd
        zcf(jl,jk)=0.0_JPRB
        ilab(jl,jk)=0
      ENDDO
    ENDDO

    !-----------------------------------------------------------------

    !!    4.           DETERMINE EXCHANGE COEFFICIENTS THEREFORE
    !!                 (A) LIFT SURFACE AIR, CHECK FOR BUOYANCY AND SET FLAG
    !!                 (B) THEN DEFINE DIFFUSION COEFFICIENTS,I.E.
    !!                      K=C1 FOR CLOUD LAYER
    !!                      K=C1*F(RH) FOR CLOUD TOP (TOP ENTRAINMENT)
    !!                  ----------------------------------------------------

    DO jl=kidia,kfdia
      ztc(jl,klev)=pten(jl,klev)+0.25_JPRB
      zqc(jl,klev)=pqen(jl,klev)
      IF(.NOT.ldcum(jl)) THEN
        ilab(jl,klev)=1
      ELSE
        ilab(jl,klev)=0
      ENDIF
      llo2(jl)=.FALSE.
      llbl(jl)=.TRUE.
    ENDDO

    DO jk=klev-1,ktdia-1+ilevh,-1

      DO jl=kidia,kfdia
        IF(pap(jl,jk) < 0.9_JPRB*paph(jl,klev+1)) llbl(jl)=.FALSE.
      ENDDO

      DO jl=kidia,kfdia
        IF(llbl(jl)) THEN
          ztc(jl,jk)=(ztc(jl,jk+1)*rcpd+pgeo(jl,jk+1)-pgeo(jl,jk))/rcpd
          zqc(jl,jk)=zqc(jl,jk+1)
          IF(ilab(jl,jk+1) > 0) THEN
            llflag(jl)=.TRUE.
          ELSE
            llflag(jl)=.FALSE.
          ENDIF
          zap(jl,jk)=pap(jl,jk)
          zpp(jl)=pap(jl,jk)
          zqold(jl)=zqc(jl,jk)
        ENDIF
      ENDDO

      DO jl=kidia,kfdia
        IF(.NOT.llbl(jl)) llflag(jl)=.FALSE.
      ENDDO

      ik=jk
      icall=1
      CALL cuadjtq &
        & ( kidia,    kfdia,    klon,   klev,&
        & ik,&
        & zpp,      ztc,      zqc,      llflag,   icall)

      DO jl=kidia,kfdia
        IF(llbl(jl)) THEN
          IF(zqc(jl,jk) /= zqold(jl)) THEN
            ilab(jl,jk)=2
          ENDIF
        ENDIF
      ENDDO

!DIR$ IVDEP
!OCL NOVREC
      DO jl=kidia,kfdia
        IF(llbl(jl)) THEN
          zbuo=ztc(jl,jk)*(1.0_JPRB+retv  *zqc(jl,jk))-&
            & pten(jl,jk)*(1.0_JPRB+retv  *pqen(jl,jk))
          IF(zbuo < 0.0_JPRB) ilab(jl,jk)=0
          IF(zbuo > 0.0_JPRB.AND.ilab(jl,jk) == 0.AND.ilab(jl,jk+1) == 1)&
            & ilab(jl,jk)=1
          IF(ilab(jl,jk) == 2) llo2(jl)=.TRUE.
        ENDIF
      ENDDO

    ENDDO

    DO jl=kidia,kfdia
      llbl(jl)=.TRUE.
    ENDDO

    DO jk=klev-1,ktdia-1+ilevh,-1

      DO jl=kidia,kfdia
        IF(pap(jl,jk) < 0.9_JPRB*paph(jl,klev+1)) llbl(jl)=.FALSE.
      ENDDO

      DO jl=kidia,kfdia
        IF(llbl(jl)) THEN
          IF(ilab(jl,jk) == 2) THEN
            zcf(jl,jk)=zkdiff1
            IF(ilab(jl,klev-2) == 0) THEN
              zcf(jl,jk)=zkdiff2
            ENDIF
          ELSE
            zcf(jl,jk)=0.0_JPRB
          ENDIF
          IF(zcf(jl,jk+1) > 0.0_JPRB.AND.ilab(jl,jk) == 0) THEN
            zcf(jl,jk)=zcf(jl,jk+1)*5._jprb*&
              & MAX(pqen(jl,jk+1)/pqsat(jl,jk+1)-0.8_JPRB,0.0_JPRB)*&
              & MAX(pqen(jl,jk+1)/pqsat(jl,jk+1)-pqen(jl,jk)/&
              & pqsat(jl,jk),0.0_JPRB)
            llbl(jl)=.FALSE.
          ENDIF
        ENDIF
      ENDDO

    ENDDO

    !*    4.7          EXCHANGE COEFFICIENTS.

    DO jk=ktdia-1+ilevh,klev-1
      DO jl=kidia,kfdia
        zcf(jl,jk)=zcf(jl,jk)*zcons1*paph(jl,jk+1)/&
          & ((pgeo(jl,jk)-pgeo(jl,jk+1))*&
          & (pten(jl,jk)+pten(jl,jk+1)))
      ENDDO
    ENDDO

    !*    4.8          DUMMY SURFACE VALUES OF T AND Q AT SURFACE

    DO jl=kidia,kfdia
      zcpts(jl)=ztpfac2*zcptgz(jl,klev)
      zqs(jl)=ztpfac2*pqen(jl,klev)
    ENDDO

    !----------------------------------------------------------------------

    !!    5.           SOLUTION OF THE VERTICAL DIFFUSION EQUATION.
    !!                 --------------------------------------------

    !*    5.1          SETTING OF RIGHT HAND SIDES.

    DO jk=ktdia-1+ilevh,klev
      DO jl=kidia,kfdia
        ztdif(jl,jk)=ztpfac2*zcptgz(jl,jk)
        zqdif(jl,jk)=ztpfac2*pqen(jl,jk)
      ENDDO
    ENDDO

    !*    5.2          TOP LAYER ELIMINATION.

    DO jl=kidia,kfdia
      ztcoe(jl)=zcf(jl,ilevh)
      zqdp=1.0_JPRB/(paph(jl,ilevh+1)-paph(jl,ilevh))
      zdisc=1.0_JPRB/(1.0_JPRB+zcf(jl,ilevh)*zqdp)
      zebs(jl,ilevh)=zdisc*(zcf(jl,ilevh)*zqdp)
      zqdif(jl,ilevh)=zdisc*zqdif(jl,ilevh)
      ztdif(jl,ilevh)=zdisc*ztdif(jl,ilevh)
    ENDDO

    !*    5.3          ELIMINATION FOR LAYERS BELOW

    DO jk=ktdia+ilevh,klev
      DO jl=kidia,kfdia
        zqdp=1.0_JPRB/(paph(jl,jk+1)-paph(jl,jk))
        zfac=ztcoe(jl)*zqdp
        ztcoe(jl)=zcf(jl,jk)
        zdisc=1.0_JPRB/(1.0_JPRB+zfac*(1.0_JPRB-zebs(jl,jk-1))+zcf(jl,jk)*zqdp)
        zebs(jl,jk)=zdisc*(zcf(jl,jk)*zqdp)
        zqdif(jl,jk)=zdisc*(zqdif(jl,jk)+zfac*zqdif(jl,jk-1))
        ztdif(jl,jk)=zdisc*(ztdif(jl,jk)+zfac*ztdif(jl,jk-1))
      ENDDO
    ENDDO

    DO jl=kidia,kfdia
      zqdif(jl,klev)=zqdif(jl,klev)+(zebs(jl,klev)*zqs(jl))
      ztdif(jl,klev)=ztdif(jl,klev)+(zebs(jl,klev)*zcpts(jl))
    ENDDO

    !*    5.5          BACK-SUBSTITUTION.

    DO jk=klev-1,ktdia-1+ilevh,-1
      DO jl=kidia,kfdia
        zqdif(jl,jk)=zqdif(jl,jk)+(zebs(jl,jk)*zqdif(jl,jk+1))
        ztdif(jl,jk)=ztdif(jl,jk)+(zebs(jl,jk)*ztdif(jl,jk+1))
      ENDDO
    ENDDO

    !---------------------------------------------------------------------

    !*    6.           INCREMENTATION OF T AND Q TENDENCIES.
    !!                   -------------------------------------

    DO jk=ktdia-1+ilevh,klev
      DO jl=kidia,kfdia
        zdqdt=(zqdif(jl,jk)-ztpfac2*pqen(jl,jk))*zcons2
        ptenq(jl,jk)=ptenq(jl,jk)+zdqdt
        zdtdt=(ztdif(jl,jk)-ztpfac2*zcptgz(jl,jk))/zcons3
        ptent(jl,jk)=ptent(jl,jk)+zdtdt
        penth(jl,jk)=(ztdif(jl,jk)-ztpfac2*zcptgz(jl,jk))*zcons2
      ENDDO
    ENDDO

    !700 CONTINUE

    IF (lhook) CALL dr_hook('CUSTRAT',1,zhook_handle)
  END SUBROUTINE custrat


END MODULE mo_cuascn

