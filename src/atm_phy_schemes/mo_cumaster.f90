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
!**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME

MODULE mo_cumaster

  USE mo_kind   ,ONLY: JPRB=>wp, vp   , &
    &                  jpim=>i4

! USE yomhook   ,ONLY : lhook,   dr_hook
  USE mo_cuparameters , ONLY :                                   &
    & lmfdd    ,lmfdudv            ,lmfit                       ,&
    & rmflic  ,rmflia  ,rmflmax, rmfsoluv, rmfdef               ,&
    & ruvper    ,rmfsoltq,rmfsolct,rmfcmin  ,lmfsmooth,lmfwstar ,&
    & lmftrac   ,   LMFUVDIS                                    ,&
    & rg       ,rd, rv      ,rcpd  ,retv , rlvtt, rlstt, rvtmp2 ,&
    & lhook,   dr_hook, icapdcycl

  USE mo_adjust,      ONLY: satur
  USE mo_cufunctions, ONLY: foelhmcu, foealfcu
  USE mo_cuinit,      ONLY: cuinin, cubasen
  USE mo_cuascn,      ONLY: cuascn
  USE mo_cudescn,     ONLY: cudlfsn, cuddrafn
  USE mo_cuflxtends,  ONLY: cuflxn, cudtdqn,cududv,cuctracer
  USE mo_cucalclpi,   ONLY: cucalclpi, cucalcmlpi
  USE mo_cucalclfd,   ONLY: cucalclfd
  USE mo_nwp_parameters,  ONLY: t_phy_params
  USE mo_nwp_tuning_config, ONLY: tune_capdcfac_et, tune_capdcfac_tr, tune_lowcapefac, &
    &                             limit_negpblcape, tune_rcapqadv, tune_capethresh
  USE mo_fortran_tools,   ONLY: t_ptr_tracer
  USE mo_exception,   ONLY: finish
  USE mo_stoch_sde,            ONLY: shallow_stoch_sde, shallow_stoch_sde_passive
  USE mo_stoch_explicit,       ONLY: shallow_stoch_explicit
  USE mo_stoch_deep,           ONLY: deep_stoch_sde
  USE mo_nwp_phy_types, ONLY: t_ptr_cloud_ensemble
  
  IMPLICIT NONE

  PRIVATE


  PUBLIC :: cumastrn

CONTAINS

  !
SUBROUTINE cumastrn &
 & (  kidia,    kfdia,    klon,   ktdia,   klev, &
 & ldland, ldlake, ptsphy, phy_params, k950,     &
 & trop_mask, mtnmask,  paer_ss,                 &
 & pten,     pqen,     puen,     pven, plitot,   &
 & pvervel,                                      &
 & plen, pien, shfl_s, qhfl_s, pqhfl,    pahfs,  &
 & pap,      paph,     pgeo,     pgeoh,          &
 & zdph,               zdgeoh,   pcloudnum,      &
 & ptent,    ptenu,    ptenv,    ptenta, ptenqa, &
 & ptenq,    ptenrhoq, ptenrhol, ptenrhoi,       &
 & ptenrhor, ptenrhos,                           &
 & ldcum,      ktype , kcbot,    kctop,          &
 & LDSHCV,   fac_entrorg, fac_rmfdeps,           &
 & pmfu,     pmfd,                               &
 & pmfude_rate,        pmfdde_rate,              &
 & ptu,      pqu,      plu,  pcore,              &
 & pmflxr,   pmflxs,   prain, pdtke_con,         &
 & pcape,    pvddraf,                            &
 & pcen, ptenrhoc,                               &
 & l_lpi, l_lfd, lpi, mlpi, koi, lfd, peis,      &
! stochastic, extra diagnostics and logical switches
 & lspinup, k650,k700, temp_s,                   &
 & cell_area,iseed,                              &
 & mf_bulk,mf_perturb,mf_num,p_cloud_ensemble,   &
 & pclnum_a, pclmf_a, pclnum_p, pclmf_p,         &                  
 & pclnum_d, pclmf_d, lacc                       )


! Code Description:
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KSTART*       FIRST STEP OF MODEL

!    *k650*         LEVEL INDEX AT 650hPa
!    *iseed*        SEED FOR RANDOM NUMBER GENERATOR

!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDLAKE*       LAKE MASK (.TRUE. FOR LAKE)
!    *lspinup*      SPINUP CONVECTIVE CLOUD ENSEMBLE
!    *L_LPI*        COMPUTE LPI, MLPI, KOI
!    *L_LFD*        COMPUTE LFD

!     INPUT PARAMETERS (REAL)

!    paer_ss    monthly aerosol climatology sea salt (optical thickness)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATIONS KG/KG
!    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT            KG/KG
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *SHFL_S*       SENSIBLE HEAT FLUX (never halo-averaged)      W/M2
!    *QHFL_S*       MOISTURE FLUX (never halo-averaged)           KG/(SM2)
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *zdgeoh*       geopot thickness on full levels               M2/S2
!    *zdph*         pressure thickness on full levels               PA
!
!    *PSSTRU*       SURFACE MOMENTUM FLUX U                - not used presently
!    *PSSTRV*       SURFACE MOMENTUM FLUX V                - not used presently
!
!    *PTENTA*       TEMPERATURE TENDENCY DYNAMICS=TOT ADVECTION    K/S
!    *PTENQA*       MOISTURE    TENDENCY DYNAMICS=TOT ADVECTION    1/S

!    *temp_s*       TEMPERATURE IN LOWEST MODEL LEVEL                K
!    *cell_area*    GRID CELL AREA                                  M2? 
!!!  ALLOCATED ONLY IF lstoch_sde=.TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *pclnum_a*     ACTIVE CLOUD NUMBER (T)               M-2
!    *pclmf_a*      ACTIVE MASS FLUX (T)                KG/(M2*S)
!    *pclnum_p*     PASSIVE CLOUD NUMBER (T)              M-2
!    *pclmf_p*      PASSIVE  MASS FLUX (T)              KG/(M2*S)

!   !!!  ALLOCATED ONLY IF lstoch_deep=.TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *pclnum_d* PROGNOSTIC DEEP CLOUD NUMBER (T)               M-2
!    *pclmf_d*  PROGNOSTIC DEEP MASS FLUX (T)                KG/(M2*S)

!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
!    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2
!    *PTENRHOC*     TENDENCY OF CHEMICAL TRACERS                  KG/(M3*S)

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    *LDSC*         FLAG: .TRUE. FOR SC-POINTS

!    OUTPUT PARAMETERS (INTEGER):

!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    *KBOTSC*       CLOUD BASE LEVEL FOR SC-CLOUDS
!!!  ALLOCATED ONLY IF lstoch_sde=.TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *pclnum_a*     PROGNOSTIC ACTIVE CLOUD NUMBER (T+1)               M-2
!    *pclmf_a*      PROGNOSTIC ACTIVE MASS FLUX (T+1)                KG/(M2 S)
!    *pclnum_p*     PROGNOSTIC PASSIVE CLOUD NUMBER (T+1)              M-2
!    *pclmf_p*      PROGNOSTIC PASSIVE  MASS FLUX (T+1)              KG/(M2 S)
!!!  ALLOCATED ONLY IF lstoch_deep=.TRUE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *pclnum_d*     PROGNOSTIC DEEP CLOUD NUMBER (T+1)                 M-2
!    *pclmf_d*      PROGNOSTIC DEEP MASS FLUX (T+1)                  KG/(M2 S)

!    OUTPUT PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PCORE*        UPDRAFT CORE FRACTION                         0-1
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M2*S)
!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
!    *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
!    *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
!                   (NO EVAPORATION IN DOWNDRAFTS)
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PDTKE_CON     CONV. BUOYANT TKE-PRODUCTION AT HALF LEVELS   M2/S**3)
!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                      KG/(M3*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M3*S)
!    *PCAPE*        CONVECTVE AVAILABLE POTENTIAL ENERGY           J/KG
!    *PWMEAN*       VERTICALLY AVERAGED UPDRAUGHT VELOCITY         M/S
!    *pvddraf*      convective gust at surface                     M/S
!    *PTENRHOQ*     MOISTURE MASS DENSITY TENDENCY                KG/(M3*S)
!    *PTENRHOL*     LIQUID WATER MASS DENSITY TENDENCY            KG/(M3*S)
!    *PTENRHOI*     ICE CONDENSATE MASS DENSITY TENDENCY          KG/(M3*S)
!    *PTENRHOR*     DETRAINED RAIN MASS DENSITY TENDENCY          KG/(M3*S)
!    *PTENRHOS*     DETRAINED SNOW MASS DENSITY TENDENCY          KG/(M3*S)
!    *LPI*          LIGHTNING POTENTIAL INDEX AS IN LYNN AND YAIR (2010) J/KG
!    *MLPI*         MODIFIED LPI USING KOI
!    *KOI*          KONVEKTIONS INDEX                             K
!    *LFD*          LIGHTNING FLASH DENSITY AS IN LOPEZ(2016)     1/(KM2*DAY)
!    *mf_bulk*      CLOUD BASE MASS FLUX FROM T-B SCHEME          KG/(M2*S)
!    *mf_perturb*   CLOUD BASE MASS FLUX FROM STOCHASTIC SCHEME   KG/(M2*S)
!    *mf_num*       NUMBER OF SHALLOW CLOUDS STOCHASTIC SCHEME    1
  

!     EXTERNALS.
!     ----------

!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V

!     SWITCHES.
!     --------

!          LMFPEN=.TRUE.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.TRUE.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.TRUE.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFIT=.TRUE.    UPDRAUGHT ITERATION
!          LMFDD=.TRUE.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.TRUE.  CUMULUS FRICTION SWITCHED ON
!          LMFTRAC=.false. TRACER TRANSPORT

!     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
!     ------------------------------------------------
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     RMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     RMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     RMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     RPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN

!     REFERENCE.
!     ----------

!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!          DRAFT PAPER ON MASSFLUX SCHEME (NORDENG, 1995)
!          Bechtold et al. (2008 QJRMS 134,1337-1351), Rooy et al. (2012 QJRMS)  
!          Bechtold et al. (2013 JAS)

!     AUTHOR.
!     -------
!      M.TIEDTKE      E.C.M.W.F.     1986/1987/1989

!     MODIFICATIONS.
!     --------------
!      03-08-29 : Clean-up, deep/shallow switches  P.Bechtold
!      04-02-11 : Add tracer transport             P.Bechtold
!      05-02-11 : Positive scaling of total Mflux  P.Bechtold
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      04-12-03 : Turn off shallow convection over stratocu. M.Ko"hler
!      05-06-27 : Switch off ddraught if idtop<kctop  
!                 correction for detrainment rates P.Bechtold
!      05-11-22 : Mods for coarser/finer physics D.Salmond + M.Hortal
!      06-02-11 : Enable TQ implicit               P.Bechtold
!      07-06-01 : Only single updraught call with  P.Bechtold
!                 scaling, convective turnover time
!                 scale, contain momentum computations in cumastrn
!      07-10-09 : Added KE dissipation and convective scavenging   P. Bechtold
!      12-03-02 : remove all entrainment stuff     P. Bechtold
!      04-10-12 : Add RPLRG/RPLDARE for small planet  N.Semane+P.Bechtold  
!      13-02-23 : modif diurnal cycle CAPE closure P. Bechtold
!      13-10-01 : modified option RCAPDCYCL=1 for diurnal cycle N. Semane
!      16-01-27 : Introduced SPP scheme (LSPP)     M. Leutbecher & S.-J. Lock
!      20180303 : Gabor: Just a comment line to force recompilation due to 
!                        compiler wrapper optimization exception liat change
!      19-08-25 : Additional total moisture advection in CAPE closure RCAPQADV=1
!----------------------------------------------------------------------

!USE PARKIND1  ,ONLY : JPIM     ,JPRB
!USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!USE YOMCST   , ONLY : RG       ,RD       ,RCPD     ,RETV     ,&
! & RLVTT    ,RLSTT   ,RTT
!USE YOMDYNCORE, ONLY : RPLRG, RPLDARE 
!USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
! & RTWAT_RTICECU_R    ,RTWAT_RTICE_R  
!USE YOECUMF  , ONLY : LMFDD    ,LMFDUDV  ,&
! & RTAU      ,RTAU0   ,RCAPDCYCL,RDEPTHS ,LMFSCV  ,LMFPEN   ,&
! & NJKT3     ,NJKT2    ,&
! & RMFCFL    ,RMFLIC  ,RMFLIA  ,RMFSOLUV ,&
! & RUVPER    ,RMFSOLTQ,RMFSOLCT,RMFCMIN  ,LMFSMOOTH,LMFWSTAR ,&
! & LMFUVDIS  ,LMFCUCA
!USE YOEPHY   , ONLY : LMFTRAC, LMFSCAV

IMPLICIT NONE

INTEGER(KIND=jpim),INTENT(in)    :: klon
INTEGER(KIND=jpim),INTENT(in)    :: klev
INTEGER(KIND=jpim),INTENT(in)    :: kidia
INTEGER(KIND=jpim),INTENT(in)    :: kfdia
INTEGER(KIND=jpim),INTENT(in)    :: k950(klon)
INTEGER(KIND=jpim)               :: ktdia
LOGICAL           ,INTENT(in)    :: ldland(klon) 
LOGICAL           ,INTENT(in)    :: ldlake(klon)
REAL(KIND=jprb)   ,INTENT(in)    :: ptsphy
TYPE(t_phy_params),INTENT(in)    :: phy_params
REAL(KIND=jprb)   ,INTENT(in)    :: trop_mask(klon)
REAL(KIND=jprb)   ,INTENT(in)    :: mtnmask(klon)
!KF
REAL(KIND=jprb)   ,INTENT(in),OPTIONAL :: paer_ss(klon)
!KF
REAL(KIND=jprb)   ,INTENT(inout) :: pten(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: pqen(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: puen(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pven(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: plitot(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pvervel(klon,klev)
REAL(KIND=jprb)   ,INTENT(in)    :: plen(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pien(klon,klev) 
!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSEN(KLON,KLEV) 
REAL(KIND=jprb)   ,INTENT(in)    :: pqhfl(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(in)    :: pahfs(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(in)    :: shfl_s(klon) 
REAL(KIND=jprb)   ,INTENT(in)    :: qhfl_s(klon)
REAL(KIND=jprb)   ,INTENT(in)    :: fac_entrorg(klon), fac_rmfdeps(klon)
REAL(KIND=jprb)   ,INTENT(in)    :: pap(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: paph(klon,klev+1)
REAL(KIND=jprb)   ,INTENT(in)    :: zdph(klon,klev)
REAL(KIND=jprb)   ,INTENT(in)    :: pgeo(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: pgeoh(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(in)    :: zdgeoh(klon,klev)
REAL(KIND=jprb)   ,INTENT(in)    :: pcloudnum(klon)
TYPE(t_ptr_tracer),INTENT(in), POINTER :: pcen(:)
TYPE(t_ptr_tracer),INTENT(inout), POINTER :: ptenrhoc(:)
REAL(KIND=jprb)   ,INTENT(inout) :: ptent(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: ptenq(klon,klev)
REAL(KIND=vp)     ,INTENT(in)    :: ptenta(klon,klev) 
REAL(KIND=jprb)   ,INTENT(in)    :: ptenqa(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: ptenrhoq(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: ptenrhol(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: ptenrhoi(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: ptenrhor(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: ptenrhos(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: ptenu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: ptenv(klon,klev) 
LOGICAL           ,INTENT(inout) :: ldcum(klon)
INTEGER(KIND=jpim),INTENT(inout) :: ktype(klon)
INTEGER(KIND=jpim),INTENT(inout) :: kcbot(klon)
INTEGER(KIND=jpim),INTENT(inout) :: kctop(klon)
!INTEGER(KIND=JPIM),INTENT(OUT)   :: KBOTSC(KLON)
!LOGICAL           ,INTENT(OUT)   :: LDSC(KLON)
LOGICAL           ,INTENT(IN)    :: LDSHCV(KLON) 
REAL(KIND=jprb)   ,INTENT(inout) :: ptu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: pqu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: plu(klon,klev) 
REAL(KIND=jprb)   ,INTENT(inout) :: pcore(klon,klev)
!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV)
!REAL(KIND=JPRB)   ,INTENT(OUT)   :: PENTH(KLON,KLEV) 
REAL(KIND=jprb)   ,INTENT(inout) :: pmflxr(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(inout) :: pmflxs(klon,klev+1) 
REAL(KIND=jprb)   ,INTENT(out)   :: pdtke_con(klon,klev+1)
REAL(KIND=jprb)   ,INTENT(out)   :: prain(klon) 
REAL(KIND=jprb)   ,INTENT(inout) :: pmfu(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: pmfd(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: pmfude_rate(klon,klev)
REAL(KIND=jprb)   ,INTENT(inout) :: pmfdde_rate(klon,klev)
REAL(KIND=jprb)   ,INTENT(out)   :: pcape(klon) 
REAL(KIND=jprb)   ,INTENT(out)   :: pvddraf(klon)
! Stochastic convection diagnostics
REAL(KIND=jprb)   ,INTENT(out)   :: mf_bulk(:) 
REAL(KIND=jprb)   ,INTENT(out)   :: mf_perturb(:) 
REAL(KIND=jprb)   ,INTENT(out)   :: mf_num(:)
! Stochastic convection switches
LOGICAL           ,INTENT(in)    :: lspinup
! Variables needed to run stochastic convection
REAL(KIND=jprb)   ,INTENT(in)    :: temp_s(klon)
INTEGER(KIND=jpim),INTENT(in)    :: k650(klon)
INTEGER(KIND=jpim),INTENT(in)    :: k700(klon)
REAL(KIND=jprb)   ,INTENT(in)    :: cell_area(klon)
INTEGER(KIND=JPIM),INTENT(in)    :: iseed(klon)

! individual cloud's properties for explicit stochastic
! cloud ensemble.
TYPE(t_ptr_cloud_ensemble), INTENT(inout) :: p_cloud_ensemble

REAL(KIND=jprb)   ,INTENT(inout)   :: pclnum_a(:)     ! prognostic active cloud number
REAL(KIND=jprb)   ,INTENT(inout)   :: pclmf_a(:)      ! prognostic active mass flux
REAL(KIND=jprb)   ,INTENT(inout)   :: pclnum_p(:)     ! prognostic passive cloud number
REAL(KIND=jprb)   ,INTENT(inout)   :: pclmf_p(:)      ! prognostic passive mass flux
REAL(KIND=jprb)   ,INTENT(inout)   :: pclnum_d(:)     ! prognostic deep cloud number
REAL(KIND=jprb)   ,INTENT(inout)   :: pclmf_d(:)      ! prognostic deep mass flux

LOGICAL           ,OPTIONAL, INTENT(in)      :: l_lpi
LOGICAL           ,OPTIONAL, INTENT(in)      :: l_lfd
REAL(KIND=jprb)   ,OPTIONAL, INTENT(inout)   :: lpi(:)
REAL(KIND=jprb)   ,OPTIONAL, INTENT(inout)   :: mlpi(:)
REAL(KIND=jprb)   ,OPTIONAL, INTENT(inout)   :: koi(:)
REAL(KIND=jprb)   ,OPTIONAL, INTENT(inout)   :: lfd(:)
REAL(KIND=jprb)            , INTENT(inout)   :: peis(:)
LOGICAL                    , INTENT(in)      :: lacc

!*UPG change to operations
REAL(KIND=jprb) :: pwmean(klon)
REAL(KIND=jprb) :: plude(klon,klev) ! only local variable
REAL(KIND=jprb) :: penth(klon,klev) ! only local variable
REAL(KIND=jprb) :: pqsen(klon,klev) ! only local variable
REAL(KIND=JPRB) :: psnde(klon,klev,2)! only local variable
!*UPG

REAL(KIND=jprb) :: ztenq_sv(klon,klev) ! MOISTURE TENDENCY ON INPUT  KG/(KG S)

REAL(KIND=jprb) :: ztenh(klon,klev),       zqenh(klon,klev),&
  & zqsenh(klon,klev),&
  & ztd(klon,klev),         zqd(klon,klev),&
  & zmfus(klon,klev),       zmfds(klon,klev),&
  & zmfuq(klon,klev),       zmfdq(klon,klev),&
  & zdmfup(klon,klev),      zdmfdp(klon,klev),&
  & zmful(klon,klev),       zrfl(klon),&
  & zuu(klon,klev),         zvu(klon,klev),&
  & zud(klon,klev),         zvd(klon,klev),&
  & zkineu(klon,klev),      zkined(klon,klev), &
  & zvbuo(klon),            zlrain(klon,klev),&
  & zsatfr(klon)
REAL(KIND=jprb) :: &
  & zmfub(klon),            zmfub1(klon),  zkhvfl(klon), zkhfl(klon),&
  & zdqcv(klon)
REAL(KIND=jprb) :: zdpmel(klon,klev),zlglac(klon,klev)
REAL(KIND=jprb) :: zdhpbl(klon),     zwubase(klon)
REAL(KIND=jprb) :: pvervel650(klon)
REAL(KIND=jprb) :: zdmfen(klon,klev),zdmfde(klon,klev)

INTEGER(KIND=jpim) ::  ilab(klon,klev), idtop(klon), ictop0(klon), ilwmin(klon)
INTEGER(KIND=jpim) ::  idpl(klon) ! departure level for convection
REAL(KIND=jprb) ::     zcape(klon), zcape2(klon), zheat(klon), zcappbl(klon), zcapdcycl(klon)
LOGICAL ::             llddraf(klon), llddraf3(klon), lldcum(klon)
LOGICAL ::             llo1, llo2(klon), llo3

INTEGER(KIND=jpim) :: ikb, ikd, itopm2, jk, ik, jl

REAL(KIND=jprb) ::   zcons2, zcons, zdh,&
  & zdqmin, zdz, zeps, zfac, &
  & zmfmax, zpbmpt, zqumqe, zro, zmfa, zerate, zderate, zorcpd, zrdocpd,&
  & ZDUTEN, ZDVTEN, ZTDIS,&
  & zalv, zsfl(klon), zcapefac, zcapethr, zmaxkined, ztenh2, zqenh2

REAL(KIND=jprb) :: ztau(klon), ztaupbl(klon)  ! adjustment time

! scaling factor for momentum and tracer massflux
REAL(KIND=jprb) :: zmfs(klon),  zmfuus(klon,klev), zmfdus(klon,klev) ,&
  & zmfudr(klon,klev), zmfddr(klon,klev) ,&
  & ZTENU(KLON,KLEV),  ZTENV(KLON,KLEV)  ,&
  & zmfuub(klon), zmfuvb(klon),&
  & ZUV2(KLON,KLEV), ZSUM12(KLON), ZSUM22(KLON),&
  & zmf_shal(klon)
    
!   parameters to calculate near-surface gusts produced by convection
REAL(KIND=jprb), PARAMETER :: conv_gust_max  = 30.0_jprb ! max. speed of conv. gusts

REAL(KIND=jprb) :: zhook_handle


! total convective flux density or dry static energy and vater vapour
REAL(KIND=jprb) :: zcvfl_s, zcvfl_q, zcaplim

!*UPG change to operations
!    LOCALS FOR CONSERVATION CHECK
LOGICAL :: llconscheck=.FALSE.
INTEGER(KIND=jpim) :: jn
REAL(KIND=jprb), DIMENSION(:,:)  , ALLOCATABLE :: ztent, ztenq
REAL(KIND=jprb), DIMENSION(:,:)  , ALLOCATABLE :: zsumc         ! kg m-2 s-1
REAL(KIND=jprb), DIMENSION(:,:,:), ALLOCATABLE :: ztenrhoc      ! kg m-3 s-1

REAL(KIND=jprb) :: deprof(klon,klev)

REAL(KIND=jprb), DIMENSION(klon) :: dummy_mfp,dummy_mfa,dummy_clnum_p,dummy_clnum_a

REAL(KIND=jprb), DIMENSION(klon)  :: zdhout
LOGICAL, PARAMETER :: luse3d   = .FALSE. !write extra output from stoch scheme
LOGICAL, PARAMETER :: lpassive = .FALSE. !run stoch schemes in piggy-backing mode

INTEGER(KIND=jpim) :: ktrac  ! number of chemical tracers

REAL(KIND=jprb) :: msee(klon,klev)

!#include "cuascn.intfb.h"
!#include "cubasen.intfb.h"
!#include "cuddrafn.intfb.h"
!#include "cudlfsn.intfb.h"
!#include "cudtdqn.intfb.h"
!#include "cududv.intfb.h"
!#include "cuflxn.intfb.h"
!#include "cuinin.intfb.h"
!#include "cuctracer.intfb.h"

!#include "fcttre.func.h"

!$ACC DATA &
!$ACC   PRESENT(k950, ldland, ldlake, trop_mask, mtnmask, pten, pqen, puen, pven) &
!$ACC   PRESENT(plitot, pvervel, pqhfl, pahfs, pap, paph, pgeo, pgeoh, zdgeoh) &
!$ACC   PRESENT(pcloudnum, ptent, ptenq, ptenrhoq, ptenrhol, ptenrhoi, ptenrhor) &
!$ACC   PRESENT(ptenrhos, ptenu, ptenv, ldcum, ktype, kcbot, kctop, ldshcv, ptu, pqu) &
!$ACC   PRESENT(plu, pmflxr, pmflxs, pdtke_con, prain, pmfu, pmfd, pmfude_rate) &
!$ACC   PRESENT(pmfdde_rate, pcape, pvddraf, phy_params, zdph, shfl_s, qhfl_s, pcore) &
!$ACC   PRESENT(ptenta, ptenqa, k700, plen, pien, peis, fac_entrorg, fac_rmfdeps) &

!$ACC   CREATE(pwmean, plude, penth, pqsen, psnde, ztenq_sv, ztenh, zqenh, zqsenh) &
!$ACC   CREATE(ztd, zqd, zmfus, zmfds, zmfuq, zmfdq, zdmfup, zdmfdp, zmful, zrfl) &
!$ACC   CREATE(zuu, zvu, zud, zvd, zkineu, zkined, zvbuo, zlrain, zmfub, zmfub1) &
!$ACC   CREATE(zkhvfl, zkhfl, zdqcv, zdpmel, zlglac, zdhpbl, zwubase, zdmfen, zdmfde) &
!$ACC   CREATE(ilab, idtop, ictop0, ilwmin, idpl, zcape, zheat, zcappbl, zcapdcycl) &
!$ACC   CREATE(llddraf, llddraf3, lldcum, llo2, zsfl, ztau, ztaupbl, zmfs, zmfuus) &
!$ACC   CREATE(zmfdus, zmfudr, zmfddr, ZTENU, ZTENV, zmfuub, zmfuvb, ZUV2, ZSUM12) &
!$ACC   CREATE(ZSUM22, zmf_shal, pvervel650, deprof, zdhout, zsatfr, zcape2, msee) &
!$ACC   IF(lacc)
    
!$ACC DATA &
!$ACC   PRESENT(lpi, mlpi, koi) &
!$ACC   IF(lacc .and. l_lpi)

!$ACC DATA &
!$ACC   PRESENT(lfd) &
!$ACC   IF(lacc .and. l_lfd)

! Special creates for the ALLOCATABLE arrays, when they are used within the code
! ztent, ztenq, zsumc, ztenrhoc

!---------------------------------------------------------------------

!     0.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------

IF (lhook) CALL dr_hook('CUMASTRN',0,zhook_handle)
zcons2  = phy_params%mfcfl/(rg*ptsphy)
zcons   = 1.0_JPRB/(rg*ptsphy)
zorcpd  = 1.0_JPRB/rcpd
zrdocpd = rd*zorcpd
zeps    = 0.0_JPRB
zcaplim = MERGE(0.0_jprb, 0.05_jprb, phy_params%lgrayzone_deepconv)


IF (ASSOCIATED(pcen) .AND. ASSOCIATED(ptenrhoc)) THEN
  ktrac = SIZE(pcen)
  !
  ! sanity check
  IF (SIZE(pcen) /= SIZE(ptenrhoc)) THEN
    CALL finish('mo_cumaster:', 'Size of pcen and ptenrhoc does not match')
  ENDIF
ELSE
  ktrac = 0
ENDIF


!---------------------------------------------------------------------
!*UPG Change to operations call SATUR routine here

!     1.           Compute Saturation specific humidity
!                  ------------------------------------


!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

!$ACC LOOP SEQ
DO jk=1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    pqsen(jl,jk)=pqen(jl,jk)
  ENDDO
ENDDO

! GZ, 2015-02-17: Change start level of qsat computation from 60 hPa (phy_params%kcon2) to ktdia,
! i.e. the general start level of moisture physics, because the previous implementation led to a crash
! in a pathological case of anomalously deep convection
CALL satur (kidia, kfdia, klon, ktdia, klev, pap, pten, pqen, pqsen, 1)


!*UPG
!set local fields zero at begin!
!dimensions: (kidia:kfdia, ktdia:klev)
!CDIR BEGIN COLLAPSE

!$ACC LOOP SEQ
DO jk=1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    ztenh (jl,jk)=0.0_JPRB
    zqsenh(jl,jk)=0.0_JPRB
    zqenh (jl,jk)=0.0_JPRB
    ztd   (jl,jk)=0.0_JPRB
    zqd   (jl,jk)=0.0_JPRB
    zmfus (jl,jk)=0.0_JPRB
    zmfds (jl,jk)=0.0_JPRB
    zmfuq (jl,jk)=0.0_JPRB
    zmfdq (jl,jk)=0.0_JPRB
    zdmfup(jl,jk)=0.0_JPRB
    zdmfdp(jl,jk)=0.0_JPRB
    zmful (jl,jk)=0.0_JPRB
    zuu   (jl,jk)=0.0_JPRB
    zvu   (jl,jk)=0.0_JPRB
    zud   (jl,jk)=0.0_JPRB
    zvd   (jl,jk)=0.0_JPRB
    zkineu(jl,jk)=0.0_JPRB
    zkined(jl,jk)=0.0_JPRB
    !    zdpmel(jl,jk)=0.0_JPRB
    zlglac(jl,jk)=0.0_JPRB
    zdmfen(jl,jk)=0.0_JPRB
    zdmfde(jl,jk)=0.0_JPRB
    ilab  (jl,jk)=0
    zmfuus(jl,jk)=0.0_JPRB
    zmfdus(jl,jk)=0.0_JPRB
    zmfudr(jl,jk)=0.0_JPRB
    zmfddr(jl,jk)=0.0_JPRB
  ENDDO
ENDDO

!dimensions: (kidia:kfdia, ktdia,klev+1)
!$ACC LOOP SEQ
DO jk=1,klev+1
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    pdtke_con(jl,jk)=0.0_JPRB
  ENDDO
ENDDO
!CDIR END

!dimensions: (kidia:kfdia)
!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  pvddraf (jl)=0.0_JPRB ! in case that it is not actually calculated !
  zmfuub  (jl)=0.0_JPRB
  zmfuvb  (jl)=0.0_JPRB
  zmf_shal(jl)=0.0_JPRB
  zmfs    (jl)=0.0_JPRB
  zdhpbl  (jl)=0.0_JPRB
  zwubase (jl)=0.0_JPRB
  zrfl    (jl)=0.0_JPRB
  !    zhcbase (jl)=0.0_JPRB
  zmfub   (jl)=0.0_JPRB
  zmfub1  (jl)=0.0_JPRB
  zdqcv   (jl)=0.0_JPRB
  idtop   (jl)=0
  ictop0  (jl)=klev ! originally 0; change needed to avoid segfaults in cuascn (but otherwise no impact on results)
  ilwmin  (jl)=0
  idpl    (jl)=0
  zcape   (jl)=0.0_JPRB
  zheat   (jl)=0.0_JPRB
  ldcum   (jl)=.FALSE.
  llddraf (jl)=.FALSE.
  llddraf3(jl)=.FALSE.
  lldcum  (jl)=.FALSE.
  llo2    (jl)=.FALSE.
  ! from new stochastic convection
  zdhout  (jl)=0.0_JPRB
ENDDO
    
!$ACC END PARALLEL

!----------------------------------------------------------------------

!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------

CALL cuinin &
  & ( kidia,    kfdia,    klon,   ktdia,    klev, phy_params%kcon2, &
  & pten,     pqen,     pqsen,    puen,     pven,&
  & pvervel,  pgeo,     paph,&
  & ilwmin,   ilab,&
  & ztenh,    zqenh,    zqsenh,   pgeoh,&
  & ptu,      pqu,      ztd,      zqd,&
  & zuu,      zvu,      zud,      zvd,&
  & plu,      lacc     )

!---------------------------------------------------------------------

!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------

!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------

IF (phy_params%lvv_shallow_deep .OR. phy_params%lvvcouple) THEN
  DO jl=kidia,kfdia
    ! vertical velocity on the model level nearest 650hPa
    ! needed for shallow/deep distinction by vertical velocity
    pvervel650(jl)=pvervel(jl,k650(jl))
  ENDDO
ENDIF

CALL cubasen &
  & ( kidia,    kfdia,    klon,   ktdia,    klev, &
  & phy_params%kcon1, phy_params%kcon2, phy_params%entrorg,fac_entrorg, &
  & phy_params%entstpc1, phy_params%entstpc2, phy_params%rdepths, &
  & phy_params%texc, phy_params%qexc, phy_params%lgrayzone_deepconv, mtnmask, ldland, ldlake, &
  & ztenh,    zqenh,    pgeoh,    paph,&
  & pqhfl,    pahfs,    &
  & pten,     pqen,     pqsen,    pgeo,&
  & puen,     pven,&
  & ptu,      pqu,      plu,      zuu,      zvu,    zwubase,&
  & ilab,     ldcum,       kcbot,   &
!&  LDSC,  KBOTSC,
  & ictop0,   idpl,     pcape, pvervel650, phy_params%lvv_shallow_deep, lacc )

!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 DECIDE ON TYPE OF CUMULUS CONVECTION
!*                 ONE THE BASIS OF THE DEPTH OF THE CONVECTION
!*                 DEEP IF CLOUD DEPTH > 200MB
!*                 SHALLOW IF CLOUD DEPTH <200MB
!                  -----------------------------------------

! CALCULATE COLUMN AND SUB CLOUD LAYER MOISTURE CONVERGENCE
! AND SUB CLOUD LAYER MOIST STATIC ENERGY CONVERGENCE

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  zdhpbl(jl)=0.0_JPRB
  idtop(jl)=0
  zcappbl(jl)=0.
  zkhvfl(jl)= -pahfs(jl,klev+1) * zorcpd - retv * pten(jl,klev) * pqhfl(jl,klev+1)
  zkhfl(JL) = -pahfs(jl,klev+1) - rlvtt * pqhfl(jl,klev+1)
ENDDO

!$ACC LOOP SEQ
DO jk=MAX(ktdia,phy_params%kcon2),klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdz)
  DO jl=kidia,kfdia
    IF(ldcum(jl).AND.jk >= kcbot(jl)) THEN
      ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))
      zdhpbl(jl)=zdhpbl(jl)+(rlvtt*ptenq(jl,jk)+rcpd*ptent(jl,jk))*zdz
      zcappbl(jl)=zcappbl(jl)+(ptent(jl,jk)+retv*pten(jl,jk)*ptenq(jl,jk))*zdz
    ENDIF
! calculate EIS diagnostic
    msee(JL,JK) = pgeo(JL,JK)+pten(JL,JK)*rcpd*(1.0_JPRB+rvtmp2*pqen(JL,JK)) - rlvtt * plen(JL,JK) - rlstt * pien(JL,JK)+ &
         &        rcpd*pten(JL,JK)*5.87_JPRB*(pqen(JL,JK)+plitot(JL,JK))
!     S=cpd*(1+5.87*q_tot)*T-Lv*ql-Lx*qi+gz
  ENDDO
ENDDO

!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
   peis(JL)=MAX(msee(JL,k700(jl))-msee(JL,k950(jl)),msee(JL,k950(jl))-msee(JL,KLEV))/rcpd
   !MAX(S700-S950; S(950)-S(surf))
ENDDO
 
!*                 ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC AND INITIAL DETERMINATION OF
!*                 CLOUD TYPE
!*                 (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------

!DIR$ IVDEP
!OCL NOVREC
    
!*                 SPECIFY INITIAL CLOUD TYPE
!*

!$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, itopm2, zpbmpt)
DO jl=kidia,kfdia
  IF (ldcum(jl)) THEN
    ikb=kcbot(jl)
    itopm2=ictop0(jl)
    zpbmpt=paph(jl,ikb)-paph(jl,itopm2)
    ! Use vertical velocity criterion to distinguish between
    ! shallow and deep cloud types. Ascending motion at 650hPa
    ! -> deep convection
    IF (phy_params%lvv_shallow_deep) THEN
       IF (pvervel650(jl).LT.0.0) THEN
          ktype(jl)=1
       ELSE
          ktype(jl)=2
       ENDIF
    ELSE
       ! Default criterion using cloud depth (in Pa) to
       ! distinguish between shallow and deep convection
       ! Clouds thicker than rdepths are deep
       IF (zpbmpt >= phy_params%rdepths) THEN
          ktype(jl)=1
       ELSE
          ktype(jl)=2
       ENDIF
    ENDIF
  ELSE
    ktype(jl)=0
  ENDIF
ENDDO

!*             (C) calculate initial updraught mass flux
!*                 and set lateral mixing rates
!*
!*                 for deep convection assume it is 10% of
!*                 maximum value which is determined by the
!*                 thickness of the layer and timestep
!*
!*                 for shallow convection calculated assuming
!*                 a balance of moist static energy in the
!*                 sub-cloud layer (ignores present of downdraughts)
!                  ------------------------------------------

!OCL NOVREC
! Use Grant closure
IF (lmfwstar) THEN
!DIR$ IVDEP
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, zdz, zmfmax)
  DO jl=kidia,kfdia
    IF (ldcum(jl)) THEN
      ikb=kcbot(jl)
      zdz=MAX(0.0_JPRB,MIN(1.5E3_JPRB,(pgeoh(jl,ikb)-pgeoh(jl,klev+1))/rg))
      zmf_shal(jl)=0.07_JPRB*(rg/pten(jl,klev)*zdz*&
        & MAX(0.0_JPRB,zkhvfl(jl)))**.3333_jprb
      ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2*RMFLIC+RMFLIA
      zmf_shal(jl)=MIN(zmf_shal(jl),zmfmax)
    ENDIF
  ENDDO
ENDIF

!$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, zmfmax, zqumqe, zdqmin, zdh)
DO jl=kidia,kfdia
  IF (ldcum(jl)) THEN
    ikb=kcbot(jl)
    ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2*RMFLIC+RMFLIA

    ! deep convection
    
    IF (ktype(jl) == 1) THEN
      ! first guess mass flux for deep convection based on arbitrary fixed
      ! parameter rmfdef
      zmfub(jl)=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*rmfdef/(rg*ptsphy)
    ELSEIF (ktype(jl) == 2) THEN

      ! shallow convection

      IF (zdhpbl(jl) > 0.0_JPRB) THEN
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
        zdqmin=MAX(0.01_JPRB*zqenh(jl,ikb),1.e-10_JPRB)
        zdh=rcpd*(ptu(jl,ikb)-ztenh(jl,ikb))+rlvtt*zqumqe
        zdh=rg*MAX(zdh,1.e5_jprb*zdqmin)
        zmfub(jl)=zdhpbl(jl)/zdh
        zmfub(jl)=MIN(zmfub(jl),0.5_jprb*zmfmax)
      ELSE
        ldcum(jl)=.FALSE.
      ENDIF!zdhpbl check
      IF(lmfwstar) zmfub(jl)=zmf_shal(jl)
    ENDIF ! if ktype=2 shallow
  ELSE !ldcum=F

   ! no buoyancy cloud base from surface
   ! set cloud base mass flux and mixing rate
   ! to default value for safety

    zmfub(jl)=0.0_JPRB
  ENDIF !ldcum=F

ENDDO

!$ACC END PARALLEL

!-----------------------------------------------------------------------

!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------

!*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------

! CALCULATIONS NOW DONE IS SECTION 3 ABOVE SO THAT
! INITIAL CLOUD DEPTH CAN BE USED TO SPECIFY
! THE TYPE OF CONVECTION

!*             (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
!                  --------------------------------------------

CALL cuascn &
  & ( kidia,    kfdia,    klon,   ktdia,   klev, phy_params%mfcfl, &
  & phy_params%entrorg, fac_entrorg, phy_params%detrpen, phy_params%rprcon, phy_params%lmfmid,      &
  & phy_params%lgrayzone_deepconv, ptsphy, paer_ss,                &
  & ztenh,    zqenh,&
  & ptenq, &
  & pten,     pqen,     pqsen,    plitot,&
  & pgeo,     pgeoh,    pap,      paph,&
  & zdph,     zdgeoh,                  &
  & pvervel,  zwubase, pcloudnum,deprof,      &
  & ldland,   ldlake,  ldcum,    ktype,    ilab,&
  & ptu,      pqu,      plu,     zlrain,        &
  & pmfu,     zmfub,    zlglac,&
  & zmfus,    zmfuq,    zmful,    plude,    zdmfup,&
  & zdmfen,   pcape,    tune_capethresh, &
  & kcbot,    kctop,    ictop0,   idpl,     pmfude_rate,   zkineu,   pwmean,    lacc )

!*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!              -----------------------------------------------------

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

!DIR$ IVDEP
!OCL NOVREC
!$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, itopm2, zpbmpt)
DO jl=kidia,kfdia
  IF (ldcum(jl)) THEN
    ikb=kcbot(jl)
    itopm2=kctop(jl)
    zpbmpt=paph(jl,ikb)-paph(jl,itopm2)
    ! Use vertical velocity criterion to distinguish between
    ! shallow and deep cloud types. Ascending motion at 650hPa
    ! -> deep convection
    IF (phy_params%lvv_shallow_deep) THEN
       ! Option: delay redefining deep cloud to shallow if no m.s.e.
       ! convergence exists, so that deep CAPE closure can be applied
       IF(ktype(jl) == 1.AND.pvervel650(jl).ge.0.0) ktype(jl)=2
       IF(ktype(jl) == 2.AND.pvervel650(jl).LT.0.0) ktype(jl)=1
    ELSE
       ! Default criterion using cloud depth (in Pa) to
       ! distinguish between shallow and deep convection
       ! Clouds thicker than rdepths are deep
       ! Option: delay redefining deep cloud to shallow if no m.s.e.
       ! convergence exists, so that deep CAPE closure can be applied
       IF(ktype(jl) == 1.AND.zpbmpt < phy_params%rdepths) ktype(jl)=2
       IF(ktype(jl) == 2.AND.zpbmpt >= phy_params%rdepths) ktype(jl)=1
    ENDIF
    ! Reset to deep convection for extreme CAPE values
    IF(pcape(jl) > tune_capethresh) ktype(jl) = 1
    ictop0(jl)=kctop(jl)
  ENDIF
  zrfl(jl)=zdmfup(jl,1)
ENDDO

!$ACC LOOP SEQ
DO jk=ktdia+1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
  ENDDO
ENDDO

!$ACC LOOP SEQ
DO jk=ktdia,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    pmfd(jl,jk)=0.0_JPRB
    zmfds(jl,jk)=0.0_JPRB
    zmfdq(jl,jk)=0.0_JPRB
    zdmfdp(jl,jk)=0.0_JPRB
    zdpmel(jl,jk)=0.0_JPRB
  ENDDO
ENDDO

IF(LMFUVDIS) THEN
  !$ACC LOOP SEQ
  DO JK=ktdia,KLEV
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO JL=KIDIA,KFDIA
      ZTENU(JL,JK)=PTENU(JL,JK)
      ZTENV(JL,JK)=PTENV(JL,JK)
    ENDDO
  ENDDO
ENDIF

!$ACC END PARALLEL

!-----------------------------------------------------------------------

!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!                  ------------------------------

IF(lmfdd) THEN

!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------

  CALL cudlfsn &
    & ( kidia,    kfdia,    klon,   ktdia,    klev,&
    & kcbot,    kctop,    ldcum, fac_rmfdeps, &
    & ztenh,    zqenh,         &
    & pten,     pqsen,    pgeo,&
    & pgeoh,    paph,     ptu,      pqu, &
    & zmfub,    zrfl,&
    & ztd,      zqd,&
    & pmfd,     zmfds,    zmfdq,    zdmfdp,&
    & idtop,    llddraf, ldland,   ldlake, lacc )

!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------

  CALL cuddrafn &
    & ( kidia,    kfdia,    klon,   ktdia,  klev,&
    & k950, llddraf, phy_params%entrdd, ztenh,    zqenh,&
    & pgeo,     pgeoh,    paph,     zrfl,&
    & zdph,     zdgeoh,                  &
    & ztd,      zqd,      pmfu,&
    & pmfd,     zmfds,    zmfdq,    zdmfdp,&
    & zdmfde,   pmfdde_rate,        zkined, zvbuo, lacc )

ENDIF

!-----------------------------------------------------------------------
!*    6.0          CLOSURE
!                  ------

!*                  RECALCULATE CLOUD BASE MASSFLUX FROM A
!*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)          
!                  --------------------------------------------

!   DEEP CONVECTION

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  zheat(jl)=0.0_JPRB
  zcape(jl)=0.0_JPRB
  zcape2(jl)=0.0_JPRB
  zmfub1(jl)=zmfub(jl)
  zdqcv(jl)=0.0_JPRB
  zsatfr(jl)=0.0_JPRB
ENDDO

!$ACC LOOP SEQ
DO jk=ktdia,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(llo1, llo3, zdz, ztenh2, zqenh2)
  DO jl=kidia,kfdia
    llo1=ldcum(jl).AND.ktype(jl) == 1
    llo3 = llo1 .AND. jk <= kcbot(jl) .AND. jk > kctop(jl)
    IF (llo3) THEN
      ZDZ=(PGEO(JL,JK-1)-PGEO(JL,JK))
      zheat(jl)=zheat(jl) +&
       & (  (pten(jl,jk-1)-pten(jl,jk) + zdz*zorcpd)/ztenh(jl,jk)&
       & +  retv*(pqen(jl,jk-1)-pqen(jl,jk))  ) *&
       & (rg*(pmfu(jl,jk)+pmfd(jl,jk)))
      zdz=(pap(JL,JK)-pap(JL,JK-1))
      zcape(jl)=zcape(jl) +&
       & ((ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)&
       & +retv*(pqu(jl,jk)-zqenh(jl,jk))&
       & -plu(jl,jk) ) * zdz
      IF (tune_rcapqadv > 0.0_JPRB) THEN
        ztenh2=ztenh(JL,JK)-0.5_JPRB*(ptenta(JL,JK)+ptenta(JL,JK-1))*ptsphy
        zqenh2=zqenh(JL,JK)-0.5_JPRB*(ptenqa(JL,JK)+ptenqa(JL,JK-1))*ptsphy
        zcape2(JL)=zcape2(JL) + ( (ptu(JL,JK)-ztenh2)/ztenh2&
         & +RETV*(pqu(JL,JK)-zqenh2) -plu(JL,JK) ) * ZDZ
      ENDIF
    ENDIF
    IF (llo1) THEN
      ZDZ=(paph(JL,JK+1)-paph(JL,JK))
      zdqcv(jl)=zdqcv(jl)+ptenqa(jl,jk)*zdz*(pqen(jl,jk)/pqsen(jl,jk))
      IF (jk >= kctop(jl)) zsatfr(jl)=zsatfr(jl)+(pqen(jl,jk)/pqsen(jl,jk))*zdz
    ENDIF
  ENDDO
ENDDO

! time scale and subcloud contribution to CAPE to be subtracted for better diurnal cycle over land

!$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikd, ikb, ik, llo1, zdz, zduten, zcapefac, zcapethr)
DO jl = kidia, kfdia
  zcapdcycl(jl) = 0.0_jprb
  IF (ldcum(jl) .AND. ktype(jl) == 1 .AND. icapdcycl == 3) THEN
    ikd = idpl(jl)
    ikb = kcbot(jl)
    ik  = kctop(jl)
    ztau(jl) = (pgeoh(jl,ik)-pgeoh(jl,ikb))/((2.0_jprb+MIN(15.0_jprb,pwmean(jl)))*rg)*phy_params%tau/ &
               MAX(1._jprb,fac_entrorg(jl)**2)
    llo1 = (paph(jl,klev+1)-paph(jl,ikd)) < 50.e2_jprb
    IF (llo1 .AND. ldland(jl)) THEN
      ! Use PBL CAPE for diurnal cycle correction, including a reduction term for low-CAPE situations
      ! and an increased correction over small-scale mountain peaks to reduce excessive precipitation maxima
      zcapethr = MERGE(100._jprb,-100._jprb,phy_params%lgrayzone_deepconv)
      zcapefac = (tune_capdcfac_et*(1._jprb-trop_mask(jl)) + tune_capdcfac_tr*trop_mask(jl)) *      &
                 MIN(1._jprb,MAX(0._jprb,(tune_lowcapefac*pcape(jl)+zcapethr)/300._jprb))
      zcapdcycl(jl) = (zcapefac+mtnmask(jl))*MAX(limit_negpblcape,zcappbl(jl))*ztau(jl)*phy_params%tau0
    ENDIF
    ! This largely suppresses convective drizzle over mountain ridges in grayzone deep convection mode
    IF (phy_params%lgrayzone_deepconv) zcapdcycl(jl) = MAX(zcapdcycl(jl),                     &
      MAX(0._jprb,mtnmask(jl)-0.2_jprb)*MERGE(10._jprb,0.1_jprb,llo1)*ztau(jl)*phy_params%tau0)
    ! Reduce adjustment time scale for extreme CAPE values
    IF (pcape(jl) > tune_capethresh) ztau(jl) = ztau(jl)*phy_params%tau0
    ! dynamic contribution to cape correction
    zdqcv(jl)=zdqcv(jl)*rlvtt/pgeoh(jl,ik)*ztau(jl)*phy_params%tau0*tune_rcapqadv
    zsatfr(jl)=zsatfr(jl)/(paph(jl,klev+1)-paph(jl,ik)) 
  ELSE IF (ldcum(jl) .AND. ktype(jl) == 1) THEN
    ikd = idpl(jl)
    ikb = kcbot(jl)
    ik  = kctop(jl)
    ztau(jl) = (pgeoh(jl,ik)-pgeoh(jl,ikb))/((2.0_jprb+MIN(15.0_jprb,pwmean(jl)))*rg)*phy_params%tau
    llo1 = (paph(jl,klev+1)-paph(jl,ikd)) < 50.e2_jprb
    IF (llo1 .AND. ldland(jl) .AND. icapdcycl==1) THEN
       zdz = MIN(1.e4_jprb,pgeoh(jl,ikb)-pgeoh(jl,klev+1))/rg
       zcapdcycl(jl) = ztau(jl)*MAX(0.0_jprb,zkhvfl(jl))*rcpd/zdz
    ENDIF
    IF (llo1 .AND. icapdcycl==2) THEN
      IF (ldland(jl)) THEN
        zcapdcycl(jl) = zcappbl(jl)*ztau(jl)*phy_params%tau0
      ELSE
        zduten = 2.0_jprb + SQRT(0.5*(puen(jl,ikb)**2 + pven(jl,ikb)**2 + &
          puen(jl,k950(jl))**2 + pven(jl,k950(jl))**2))
        ztaupbl(jl) = MIN(1.e4_jprb, pgeoh(jl,ikb)-pgeoh(jl,klev+1))/(rg*zduten)
        zcapdcycl(jl) = zcappbl(jl)*ztaupbl(jl)
      ENDIF
    ENDIF
  ENDIF
ENDDO

!$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, zmfmax)
DO jl=kidia,kfdia
  IF(ldcum(jl).AND.ktype(jl) == 1) THEN
    ikb=kcbot(jl)
    zcape2(jl)=tune_rcapqadv*zcape2(jl)+(1.0_JPRB-tune_rcapqadv)*zcape(jl)
    zcapdcycl(jl)=MAX(zcapdcycl(jl),-2*zcape2(jl))
    IF(zsatfr(jl)<=0.94_JPRB) THEN
!GZ: the additional limitation to 2.5*zcape2 is needed in ICON
      zcape(jl)=MIN(2.5_jprb*zcape2(jl),MAX(zcaplim*zcape(jl),zcape2(jl)-zcapdcycl(jl)+zdqcv(jl)))
    ELSE
      zcape(jl)=MAX(zcaplim*zcape(jl),zcape(jl)-zcapdcycl(jl))
    ENDIF
    zcape(jl)=MAX(0.0_jprb,MIN(zcape(jl),5000.0_jprb))
    zheat(jl)=MAX(1.e-4_jprb,zheat(jl))
    ztau(jl)=MAX(720._jprb,ztau(jl))
    zmfub1(jl)=(zcape(jl)*zmfub(jl))/(zheat(jl)*ztau(jl))
    zmfub1(jl)=MAX(zmfub1(jl),0.001_jprb)
    zmfmax=(paph(jl,ikb)-paph(jl,ikb-1))*zcons2*rmflic+rmflia
    zmfub1(jl)=MIN(zmfub1(jl),zmfmax)
  ENDIF
ENDDO

!  SHALLOW CONVECTION AND MID_LEVEL

!DIR$ IVDEP
!OCL NOVREC

!$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, zeps, zqumqe, zdqmin, zmfmax, zdh)
DO jl=kidia,kfdia
  IF ( ldcum(jl) .AND. (ktype(jl) == 2.OR. ktype(jl) == 3) ) THEN
    ikb=kcbot(jl)
    IF(pmfd(jl,ikb) < 0.0_JPRB) THEN
      zeps=-pmfd(jl,ikb)/MAX(zmfub(jl),1.e-10_JPRB)
    ELSE
      zeps=0.0_JPRB
   ENDIF
    zdqmin=MAX(0.01_JPRB*zqenh(jl,ikb),1.e-10_JPRB)
    ! maximum permisable value of ud base mass flux
    zmfmax=(paph(jl,ikb)-paph(jl,ikb-1))*zcons2*rmflic+rmflia

    ! shallow convection

    IF(ktype(jl) == 2) THEN

      IF(zdhpbl(jl) > 0.0_JPRB) THEN
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-&
             & zeps*zqd(jl,ikb)-(1.0_JPRB-zeps)*zqenh(jl,ikb)
        zdh=rcpd*(ptu(jl,ikb)-zeps*ztd(jl,ikb)-&
             & (1.0_JPRB-zeps)*ztenh(jl,ikb))+rlvtt*zqumqe
        zdh=rg*MAX(zdh,1.e5_jprb*zdqmin)
        zmfub1(jl)=zdhpbl(jl)/zdh
        ! m.s.e. difference between updraft and environment
        ! is required in shallow stochastic scheme. Save!
        zdhout(jl)=zdh
      ELSE
        !MA: cleanup convection types and set default values
        zmfub1(jl)=0.0_JPRB
        ldcum(jl)=.FALSE.
        ktype(jl)=0
      ENDIF !zdhpbl check
      ! Why is shallow mass flux limited by half the zmfmax, but deep
      ! mass flux limited by zmfmax?
      zmfub1(jl)=MIN(zmfub1(jl),0.5_jprb*zmfmax)

      IF(lmfwstar) zmfub1(jl)=zmf_shal(jl)
    ENDIF ! if shallow convection

! mid-level convection

    IF(ktype(jl) == 3)THEN
      zmfub1(jl)=zmfub(jl)*(1.0_JPRB+zeps)
      zmfub1(jl)=MAX(zmfub(jl),zkhfl(jl)/pgeoh(jl,ikb))&
     &                    *(1.0_JPRB+zeps)
      zmfub1(jl)=MIN(zmfub1(jl),zmfmax)
   ENDIF
   
  ENDIF
ENDDO
!$ACC END PARALLEL

!!!BEGINNING OF STOCHASTIC ROUTINES

IF (phy_params%lstoch_expl .or. phy_params%lstoch_sde .or. phy_params%lstoch_deep) THEN

  ! Make sure zdhout is calculated - is required for stoch routines,
  ! and may not have been calculated for points redefined from ktype
  ! 1 to ktype 2
  DO jl=kidia,kfdia
    IF(ldcum(jl)) THEN
      ! last check on convection type
      ikb=kcbot(jl)
      itopm2=kctop(jl)
      zpbmpt=paph(jl,ikb)-paph(jl,itopm2)
      ! In case m.s.e. difference between updraft and environment
      ! was not calculated as part of the shallow closure,
      ! calculate it here. Is required in shallow stochastic scheme.
      ! This can be the case if Boeing or CAPE closures were used
      ! instead of m.s.e. closure, e.g. because point was redefined
      ! from deep to shallow, based on cloud depth.
      IF (zdhout(jl).eq.0._JPRB) THEN
        zdqmin=MAX(0.01_JPRB*zqenh(jl,ikb),1.e-10_JPRB)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-&
              & zeps*zqd(jl,ikb)-(1.0_JPRB-zeps)*zqenh(jl,ikb)
        zdh=rcpd*(ptu(jl,ikb)-zeps*ztd(jl,ikb)-&
              & (1.0_JPRB-zeps)*ztenh(jl,ikb))+rlvtt*zqumqe
        zdh=rg*MAX(zdh,1.e5_jprb*zdqmin)!
        zdhout(jl) = zdh
      ENDIF
    ENDIF
  ENDDO
   
  ! Save mass flux calculated by convectional T-B scheme 
  ! into mf_bulk diagnostic, and initialise the stochastically
  ! perturbed mass flux diagnostic mf_perturb.
  ! This applies to both versions of the stochastic scheme
  DO jl=kidia,kfdia
    mf_bulk(jl) = zmfub1(jl)
    mf_perturb(jl) = 0._JPRB
  ENDDO
ENDIF

! Note the order of the routines! The following four calls cover
! the options to run the shallow stochastic scheme explicitly or
! as SDE, with or without "piggy-backing" mode enabled.
! In piggy backing mode, it is crucial that the passive SDE routine
! is called before the explicit routine such that it can act on the
! convection state of the previous time step.


IF( phy_params%lstoch_expl .and. lpassive ) THEN
  ! Call SDE stochastic scheme (passively, in piggy-backing mode)
  CALL shallow_stoch_sde_passive(                                  &
&                         i_startidx   = kidia,                    & !IN
&                         i_endidx     = kfdia,                    & !IN
&                         klon         = klon,                     & !IN
&                         klev         = klev,                     & !IN
&                         ptsphy       = ptsphy,                   & !IN
&                         pgeoh        = pgeoh,                    & !IN
&                         mbas_con     = kcbot,                    & !IN
&                         mfb          = zmfub1,                   & !IN
&                         shfl         = pahfs(:,klev+1),          & !IN
&                         lhfl         = pqhfl(:,klev+1)*rlvtt,    & !IN
&                         temp_s       = temp_s,                   & !IN
&                         dh           = zdhout,                   & !IN
&                         ktype        = ktype,                    & !IN
!&                         extra_3d     = extra_3d,                 & !INOUT
&                         mfp          = pclmf_p,                  & !IN
&                         mfa          = pclmf_a,                  & !IN
&                         clnum_p      = pclnum_p,                 & !IN
&                         clnum_a      = pclnum_a,                 & !IN
&                         lseed        = iseed,                    & !IN
&                         luse3d       = luse3d,                   & !IN
&                         cell_area    = cell_area,                & !IN
&                         lgrayzone    = phy_params%lgrayzone_deepconv)
ENDIF

IF( phy_params%lstoch_expl ) THEN
  ! Call explicit stochastic scheme (interactively)
  CALL  shallow_stoch_explicit(                                       &
&                             i_startidx = kidia,                     & !IN
&                             i_endidx   = kfdia,                     & !IN
&                             klon       = klon,                      & !IN
&                             klev       = klev,                      & !IN
&                             dt         = ptsphy,                    & !IN
&                             pgeoh      = pgeoh,                     & !IN
&                             mbas_con   = kcbot,                     & !IN
&                             mtop_con   = kctop,                     & !IN
&                             mf_bulk    = zmfub1,                    & !IN
&                             shfl       = shfl_s,                    & !IN
&                             lhfl       = qhfl_s*rlvtt,              & !IN
&                             temp_s     = temp_s,                    & !IN
&                             mf_perturb = mf_perturb,                & !OUT 
&                             mfp        = pclmf_p,                   & !OUT
&                             mfa        = pclmf_a,                   & !OUT
&                             dh         = zdhout,                    & !IN
&                             ktype      = ktype,                     & !IN
&                             clnum      = mf_num,                    & !OUT
&                             lseed      = iseed,                     & !IN
&                             cell_area  = cell_area,                 & !IN
&                             core       = pcore,                     & !OUT
&                             deprof     = deprof,                    & !OUT
!&                             extra_3d   = extra_3d,                  & !INOUT
&                             ncloudspout= pclnum_p,                  & !OUT
&                             ncloudsaout= pclnum_a,                  & !OUT
&                             time_i     = p_cloud_ensemble%time_i,   & !INOUT
&                             life_i     = p_cloud_ensemble%life_i,   & !INOUT
&                             mf_i       = p_cloud_ensemble%mf_i,     & !INOUT
&                             type_i     = p_cloud_ensemble%type_i,   & !INOUT
&                             ktype_i    = p_cloud_ensemble%ktype_i,  & !INOUT
&                             area_i     = p_cloud_ensemble%area_i,   & !INOUT
&                             depth_i    = p_cloud_ensemble%depth_i,  & !INOUT
&                             base_i     = p_cloud_ensemble%base_i,   & !INOUT
&                             used_cell  = p_cloud_ensemble%used_cell,& !INOUT
&                             lpassive   = .FALSE.,                   & !IN
&                             luse3d     = luse3d,                    & !IN
&                             lspinup    = lspinup,                   & !IN
&                             lgrayzone  = phy_params%lgrayzone_deepconv)!IN

  ! For low bulk MF, no cloud may be generated in the stochastic scheme,
  ! leading to zero perturbed MF. Point should be switched off in this
  ! case. Otherwise, overwrite "final" mass flux after closure from
  ! conventional T-B scheme with stochastically perturbed value.
  DO jl=kidia,kfdia
    IF ((ktype(jl) .EQ. 2 .OR. (phy_params%lgrayzone_deepconv .and. (ktype(jl) .eq. 1))) &
         & .AND. (mf_perturb(jl) .GT. 0._JPRB) ) THEN
      zmfub1(jl) = mf_perturb(jl)
    ELSE
      zmfub1(jl) = 0._JPRB
      ktype(jl)  = 0
      ldcum(jl)  = .FALSE.
      kcbot(jl)  = 0
      kctop(jl)  = 0
    ENDIF
  ENDDO
   
ENDIF !lstoch_expl

IF( phy_params%lstoch_sde.and.lpassive ) THEN
  ! Call explicit stochastic scheme (passively in piggy-backing mode)
  CALL  shallow_stoch_explicit(                                       &
&                             i_startidx = kidia,                     & !IN
&                             i_endidx   = kfdia,                     & !IN
&                             klon       = klon,                      & !IN
&                             klev       = klev,                      & !IN
&                             dt         = ptsphy,                    & !IN
&                             pgeoh      = pgeoh,                     & !IN
&                             mbas_con   = kcbot,                     & !IN
&                             mtop_con   = kctop,                     & !IN
&                             mf_bulk    = zmfub1,                    & !IN
&                             shfl       = shfl_s,                    & !IN
&                             lhfl       = qhfl_s*rlvtt,              & !IN
&                             temp_s     = temp_s,                    & !IN
&                             mf_perturb = mf_perturb,                & !OUT
&                             mfp        = dummy_mfp,                 & !OUT
&                             mfa        = dummy_mfa,                 & !OUT
&                             dh         = zdhout,                    & !IN
&                             ktype      = ktype,                     & !IN
&                             clnum      = mf_num,                    & !OUT
&                             lseed      = iseed,                     & !IN
&                             cell_area  = cell_area,                 & !IN
&                             core       = pcore,                     & !OUT
&                             deprof     = deprof,                    & !OUT
!&                             extra_3d   = extra_3d,                  & !INOUT
&                             ncloudspout= dummy_clnum_p,             & !OUT
&                             ncloudsaout= dummy_clnum_a,             & !OUT
&                             time_i     = p_cloud_ensemble%time_i,   & !INOUT
&                             life_i     = p_cloud_ensemble%life_i,   & !INOUT
&                             mf_i       = p_cloud_ensemble%mf_i,     & !INOUT
&                             type_i     = p_cloud_ensemble%type_i,   & !INOUT
&                             ktype_i    = p_cloud_ensemble%ktype_i,  & !INOUT
&                             area_i     = p_cloud_ensemble%area_i,   & !INOUT
&                             depth_i    = p_cloud_ensemble%depth_i,  & !INOUT
&                             base_i     = p_cloud_ensemble%base_i,   & !INOUT
&                             used_cell  = p_cloud_ensemble%used_cell,& !INOUT
&                             lpassive   = .TRUE.,                    & !IN
&                             luse3d     = luse3d,                    & !IN
&                             lspinup    = lspinup,                   & !IN
&                             ncloudspin = pclnum_p,                  & !IN, OPT
&                             ncloudsain = pclnum_a,                  & !IN, OPT
&                             lgrayzone  = phy_params%lgrayzone_deepconv)!IN

ENDIF

IF(phy_params%lstoch_sde) THEN
  ! Call SDE stochastic scheme (interactively)
  CALL shallow_stoch_sde(                                          &
&                         i_startidx   = kidia,                    & !IN
&                         i_endidx     = kfdia,                    & !IN
&                         klon         = klon,                     & !IN
&                         klev         = klev,                     & !IN
&                         ptsphy       = ptsphy,                   & !IN
&                         pgeoh        = pgeoh,                    & !IN
&                         mbas_con     = kcbot,                    & !IN
&                         mfb          = zmfub1,                   & !IN
&                         shfl         = pahfs(:,klev+1),          & !IN
&                         lhfl         = pqhfl(:,klev+1)*rlvtt,    & !IN
&                         temp_s       = temp_s,                   & !IN
&                         mfp          = mf_perturb,               & !OUT
&                         dh           = zdhout,                   & !IN
&                         ktype        = ktype,                    & !IN
&                         clnum        = mf_num,                   & !OUT
&                         lseed        = iseed,                    & !IN
&                         luse3d       = luse3d,                   & !IN
&                         cell_area    = cell_area,                & !IN
&                         pclnum_a     = pclnum_a,                 & !OUT
&                         pclmf_a      = pclmf_a,                  & !OUT
&                         pclnum_p     = pclnum_p,                 & !OUT
&                         pclmf_p      = pclmf_p,                  & !OUT
&                         lspinup      = lspinup,                  & !IN
!&                         extra_3d     = extra_3d,                 & !INOUT
&                         lgrayzone    = phy_params%lgrayzone_deepconv)

  ! For low bulk MF, no cloud may be generated in the stochastic scheme,
  ! leading to zero perturbed MF. Point should be switched off in this
  ! case. Otherwise, overwrite "final" mass flux after closure from
  ! conventional T-B scheme with stochastically perturbed value.
  DO jl=kidia,kfdia
     IF ((ktype(jl) .EQ. 2 .OR. (phy_params%lgrayzone_deepconv .and. (ktype(jl) .eq. 1))) &
          & .AND. (mf_perturb(jl) .GT. 0._JPRB) ) THEN
      zmfub1(jl) = mf_perturb(jl)
    ELSE
      zmfub1(jl) = 0._JPRB
      ktype(jl)=0
      ldcum(jl)=.FALSE.
      kcbot(jl)=0
      kctop(jl)=0
    ENDIF
  ENDDO
ENDIF !lstoch_sde

IF(phy_params%lstoch_deep) THEN
  ! Call SDE stochastic scheme for deep convection. This should only be used
  ! at appropriately coarse (global) resolution, and not be used in conjunction
  ! with the shallow stochastic scheme. If resolution is coarse enough to
  ! parameterise deep convection, the shallow scheme converges back to the
  ! conventional T-B scheme, so there is no added benefit.
  CALL deep_stoch_sde( &
&                         i_startidx   = kidia,                    & !IN
&                         i_endidx     = kfdia,                    & !IN
&                         klon         = klon,                     & !IN
&                         ptsphy       = ptsphy,                   & !IN
&                         mfb          = zmfub1,                   & !IN
&                         mfp          = mf_perturb,               & !OUT
&                         ktype        = ktype,                    & !IN
&                         clnum        = mf_num,                   & !OUT
&                         lseed        = iseed,                    & !IN
&                         luse3d       = luse3d,                   & !IN
&                         cell_area    = cell_area,                & !IN
&                         pclnum_d     = pclnum_d,                 & !OUT
&                         pclmf_d      = pclmf_d                   ) !OUT
!&                         pclmf_d      = pclmf_d,                  & !OUT
!&                         extra_3d     = extra_3d)                   !INOUT

  ! For low bulk MF, no cloud may be generated in the stochastic scheme,
  ! leading to zero perturbed MF. Point should be switched off in this
  ! case. Otherwise, overwrite "final" mass flux after closure from
  ! conventional T-B scheme with stochastically perturbed value.
  DO jl=kidia,kfdia
    ! Only modify grid points with deep convection active
    IF (ktype(jl) .EQ. 1) THEN
      if (mf_perturb(jl) .gt. 0._JPRB) THEN
        zmfub1(jl) = mf_perturb(jl)
      ELSE
        ! If stochastic scheme didn't generate clouds despite test parcel ascent
        ! indicating deep convection, switch off point
        zmfub1(jl) = 0._JPRB
        ktype(jl)=0
        ldcum(jl)=.FALSE.
        kcbot(jl)=0
        kctop(jl)=0
      ENDIF
    ENDIF
  ENDDO
ENDIF !lstoch_deep

!!!END OF STOCHASTIC ROUTINES

! rescale DD fluxes if deep and shallow convection

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
!$ACC LOOP SEQ
DO jk=ktdia,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zfac)
  DO jl=kidia,kfdia
    IF ( llddraf(jl) .AND.( ktype(jl) == 1.OR. ktype(jl) == 2 ) ) THEN
      zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_JPRB)
      pmfd(jl,jk)=pmfd(jl,jk)*zfac
      zmfds(jl,jk)=zmfds(jl,jk)*zfac
      zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
      zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
!  also rescale detrainment flux for ERA pp
      pmfdde_rate(jl,jk)=pmfdde_rate(jl,jk)*zfac
    ENDIF
  ENDDO
ENDDO
!$ACC END PARALLEL

! Updraft iteration is .FALSE. by default
IF(lmfit) THEN

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
  !$ACC LOOP SEQ
  DO jk=ktdia+1,klev-1
    !$ACC LOOP GANG VECTOR
    DO jl=kidia,kfdia
      zuu(jl,jk)=puen(jl,jk-1)
      zvu(jl,jk)=pven(jl,jk-1)
    ENDDO
  ENDDO

  ! reset updraught mass flux at cloud base

  !$ACC LOOP GANG VECTOR
  DO jl=kidia,kfdia
    zmfub(jl)=zmfub1(jl)
  ENDDO
  !$ACC END PARALLEL

  !-----------------------------------------------------------------------

  !*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
  !*                 FOR PENETRATIVE CONVECTION (TYPE=1),
  !*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
  !*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
  !                  -------------------------------------------------

  CALL cuascn &
    & ( kidia,    kfdia,    klon,   ktdia,   klev, phy_params%mfcfl, &
    & phy_params%entrorg, fac_entrorg, phy_params%detrpen,phy_params%rprcon, phy_params%lmfmid, &
    & phy_params%lgrayzone_deepconv, ptsphy, paer_ss,&
    & ztenh,    zqenh,    &
    & ptenq,            &
    & pten,     pqen,     pqsen,    plitot,&
    & pgeo,     pgeoh,    pap,      paph,&
    & zdph,     zdgeoh,         &
    & pvervel,  zwubase,  pcloudnum,deprof,     &
    & ldland,   ldlake,   ldcum,    ktype,    ilab,&
    & ptu,      pqu,      plu,      zlrain,        &
    & pmfu,     zmfub,    zlglac,&
    & zmfus,    zmfuq,    zmful,    plude,    zdmfup,&
    & zdmfen,   pcape,    tune_capethresh, &
    & kcbot,    kctop,    ictop0,   idpl,     pmfude_rate,    zkineu,   pwmean, lacc )

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
  !$ACC LOOP GANG VECTOR PRIVATE(ikb, itopm2, zpbmpt)
  DO jl=kidia,kfdia
    IF (ldcum(jl)) THEN
      ikb=kcbot(jl)
      itopm2=kctop(jl)
      zpbmpt=paph(jl,ikb)-paph(jl,itopm2)
      IF (phy_params%lvv_shallow_deep) THEN 
         IF(ktype(jl) == 1.AND.pvervel650(jl).ge.0.0) ktype(jl)=2
         IF(ktype(jl) == 2.AND.pvervel650(jl).LT.0.0) ktype(jl)=1
      ELSE
         IF(ktype(jl) == 1.AND.zpbmpt < phy_params%rdepths) ktype(jl)=2
         IF(ktype(jl) == 2.AND.zpbmpt >= phy_params%rdepths) ktype(jl)=1
      ENDIF
      ! Reset to deep convection for extreme CAPE values
      IF(pcape(jl) > tune_capethresh) ktype(jl) = 1
    ENDIF
  ENDDO
  !$ACC END PARALLEL

ELSE

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    IF(ldcum(jl)) THEN
      zmfs(jl)=zmfub1(jl)/MAX(rmfcmin,zmfub(jl))
    ENDIF
  ENDDO

  !$ACC LOOP SEQ
  DO jk=ktdia+1,klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, zdz, zmfmax)
    DO jl=kidia,kfdia
      IF(ldcum(jl).AND.jk>=kctop(jl)-1) THEN
        ikb=kcbot(jl)
        IF(jk>ikb) THEN
          zdz=((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
          pmfu(jl,jk)=pmfu(jl,ikb)*zdz
        ENDIF
        zmfmax=MIN(rmflmax,(paph(jl,jk)-paph(jl,jk-1))*zcons2*rmflic+rmflia)
        IF(pmfu(jl,jk)*zmfs(jl)>zmfmax) &
          & zmfs(jl)=MIN(zmfs(jl),zmfmax/pmfu(jl,jk))
      ENDIF
    ENDDO                       !
  ENDDO

  !$ACC LOOP SEQ
  DO jk=ktdia+1,klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      IF(ldcum(jl).AND.jk<=kcbot(jl).AND.jk>=kctop(jl)-1) THEN
        pmfu(jl,jk)=pmfu(jl,jk)*zmfs(jl)
        zmfus(jl,jk)=zmfus(jl,jk)*zmfs(jl)
        zmfuq(jl,jk)=zmfuq(jl,jk)*zmfs(jl)
        zmful(jl,jk)=zmful(jl,jk)*zmfs(jl)
        zdmfup(jl,jk)=zdmfup(jl,jk)*zmfs(jl)
        zdmfen(jl,jk)=zdmfen(jl,jk)*zmfs(jl)
        plude(jl,jk)=plude(jl,jk)*zmfs(jl)
        pmfude_rate(jl,jk)=pmfude_rate(jl,jk)*zmfs(jl)
        ! Calculate a simple estimate of the updraft core fraction
        ! by assuming a density of 1kg/kg, and an updraft velocity
        ! that is the square root of the updraft kinetic energy.
        ! Only perform this calculation if the updraft fraction is
        ! not being calculated by the explicit stochastic scheme and
        ! on layers where the updraft kinetic energy exceeds a minimum
        ! threshold of 0.25 (to avoid excessively large updraft area 
        ! near cloud base/top. The factor 2 is a tuning factor?
        IF (.NOT.(phy_params%lstoch_expl) .AND. &
             zkineu(jl,jk) > 0.25_JPRB .AND. ktype(jl) == 2) THEN 
           pcore(jl,jk)=pmfu(jl,jk)/(2._jprb*SQRT(zkineu(jl,jk)))
           pcore(jl,jk)=MAX(0._JPRB,pcore(jl,jk))
        ENDIF
      endif
    ENDDO
  ENDDO
  !$ACC END PARALLEL

ENDIF

!-----------------------------------------------------------------------

!*    6.5          IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
!                  RESET LDCUM TO FALSE-> FLUXES SET TO ZERO IN CUFLXN
!                  ---------------------------------------------------

!                 exclude pathological KTYPE=2 KCBOT=KCTOP=KLEV-1

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  IF(ktype(jl)==2.AND.kcbot(jl)==kctop(jl).AND.kcbot(jl)>=klev-1) THEN
    ldcum(jl)=.FALSE.
    ktype(jl)=0
  ENDIF
ENDDO

!USE EIS to switch off convection in stable conditions
!$ACC LOOP GANG(STATIC: 1) VECTOR
DO JL=KIDIA,KFDIA
   IF (peis(JL)>phy_params%eiscrit .AND. ktype(jl)==2) THEN
      llo2(jl)=.TRUE.
      ldcum(jl)=.FALSE.
   ENDIF
ENDDO
!                  turn off shallow convection if stratocumulus PBL type
!$ACC LOOP GANG(STATIC: 1) VECTOR
DO JL=KIDIA,KFDIA
!xmk IF((.NOT.LDSHCV(JL) .AND. KTYPE(JL)==2)) THEN
!RN added condition: maintain shallow cumulus starting above lowest level (KLEV)
  IF((.NOT.LDSHCV(JL) .AND. KTYPE(JL)==2 .AND. IDPL(JL)==KLEV)) THEN
!xxx
    LLO2(JL)=.TRUE.
    LDCUM(JL)=.FALSE.
  ENDIF
ENDDO


IF (.NOT.phy_params%lmfscv .OR. .NOT.phy_params%lmfpen) THEN
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    IF((.NOT.phy_params%lmfscv .AND. ktype(jl)==2).OR.(.NOT.phy_params%lmfpen .AND. ktype(jl)==1))THEN
      llo2(jl)=.TRUE.
      ldcum(jl)=.FALSE.
    ENDIF
  ENDDO
ENDIF

IF (phy_params%lvvcouple) THEN
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    ! Use vertical velocity as a criterion to decide when to turn off parameterization
    ! at a grid point and resolve convection.
    ! For rising motion at 650hPa, turn off parameterized convection.
    IF(.NOT.phy_params%lmfpen .AND. ktype(jl)==2 .AND. pvervel650(jl) < 0.0_jprb)THEN
      llo2(jl)=.TRUE.
      ldcum(jl)=.FALSE.
      ktype(jl)=0
      zmfub(jl)=0.0_jprb
    ENDIF
  ENDDO
ENDIF

!-----------------------------------------------------------------------

!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------

!- set DD mass fluxes to zero above cloud top
!  (because of inconsistency with second updraught)

!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  IF(llddraf(jl).AND.idtop(jl)<=kctop(jl)) THEN
    idtop(jl)=kctop(jl)+1
  ENDIF
ENDDO

!$ACC LOOP SEQ
DO jk=ktdia+1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    IF (llddraf(jl)) THEN
      IF (jk<idtop(jl)) THEN
        pmfd(jl,jk)=0.0_JPRB
        zmfds(jl,jk)=0.0_JPRB
        zmfdq(jl,jk)=0.0_JPRB
        pmfdde_rate(jl,jk)=0.0_JPRB
        zdmfdp(jl,jk)=0.0_JPRB
      ELSEIF (jk==idtop(jl)) THEN
        pmfdde_rate(jl,jk)=0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
ENDDO
!$ACC END PARALLEL

CALL cuflxn &
  & ( kidia,    kfdia,    klon,   ktdia,    klev, phy_params%mfcfl, &
  & phy_params%rhebc_land, phy_params%rhebc_ocean, phy_params%rcucov, &
  & phy_params%rhebc_land_trop, phy_params%rhebc_ocean_trop, &
  & phy_params%rcucov_trop, phy_params%lmfdsnow, trop_mask,   &
  & ptsphy,  pten,     pqen,     pqsen,    ztenh,    zqenh,&
  & paph,     pap,      pgeoh,    ldland,   ldlake, ldcum,&
  & kcbot,    kctop,    idtop,    itopm2,&
  & ktype,    llddraf,&
  & pmfu,     pmfd,     zmfus,    zmfds,&
  & zmfuq,    zmfdq,    zmful,    plude,  zlrain,  psnde, &
  & zdmfup,   zdmfdp,   zdpmel,   zlglac,&
  & pmflxr,   pmflxs,   prain,    pmfude_rate,  pmfdde_rate, lacc )

!- rescale DD fluxes if total mass flux becomes negative
!- correct DD detrainment rates if entrainment becomes negative
!- correct UD detrainment rates if entrainment becomes negative
!- conservation correction for precip

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

!$ACC LOOP GANG(STATIC: 1) VECTOR
DO jl=kidia,kfdia
  zmfs(jl)=1.0_JPRB
ENDDO

!DO JK=2,KLEV-1
!$ACC LOOP SEQ
DO jk=ktdia+1,klev ! change for stability
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zmfmax)
  DO jl=kidia,kfdia
    IF ( llddraf(jl) .AND. jk>=idtop(jl)-1 ) THEN
      zmfmax=pmfu(jl,jk)*0.98_JPRB
      IF(pmfd(jl,jk)+zmfmax+1.e-15_JPRB<0.0_JPRB) THEN
        zmfs(jl)=MIN(zmfs(jl),-zmfmax/pmfd(jl,jk))
      ENDIF
    ENDIF
  ENDDO
ENDDO

! done above:  zmfuub(:)=0.0_JPRB

!$ACC LOOP SEQ
DO jk=ktdia+1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    IF ( zmfs(jl)<1.0_JPRB .AND. jk>=idtop(jl)-1 ) THEN
      pmfd(jl,jk)=pmfd(jl,jk)*zmfs(jl)
      zmfds(jl,jk)=zmfds(jl,jk)*zmfs(jl)
      zmfdq(jl,jk)=zmfdq(jl,jk)*zmfs(jl)
      pmfdde_rate(jl,jk)=pmfdde_rate(jl,jk)*zmfs(jl)
      zmfuub(jl)=zmfuub(jl)-(1.0_JPRB-zmfs(jl))*zdmfdp(jl,jk)
      pmflxr(jl,jk+1)=pmflxr(jl,jk+1)+zmfuub(jl)
      zdmfdp(jl,jk)=zdmfdp(jl,jk)*zmfs(jl)
    ENDIF
  ! ZDMFUPC(JL,JK)=ZDMFUP(JL,JK)
  ! ZDMFDPC(JL,JK)=ZDMFDP(JL,JK)
  ENDDO
ENDDO

!$ACC LOOP SEQ
DO jk=ktdia+1,klev-1
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zerate)
  DO jl=kidia,kfdia
    IF ( llddraf(jl) .AND. jk>=idtop(jl)-1 ) THEN
      zerate=-pmfd(jl,jk)+pmfd(jl,jk-1)+pmfdde_rate(jl,jk)
      IF(zerate<0.0_JPRB) THEN
        pmfdde_rate(jl,jk)=pmfdde_rate(jl,jk)-zerate
      ENDIF
    ENDIF
    IF ( ldcum(jl) .AND. jk>=kctop(jl)-1 ) THEN
      zerate=pmfu(jl,jk)-pmfu(jl,jk+1)+pmfude_rate(jl,jk)
      IF(zerate<0.0_JPRB) THEN
        pmfude_rate(jl,jk)=pmfude_rate(jl,jk)-zerate
      ENDIF
      ! ZDMFUP(JL,JK)=ZDMFUP(JL,JK)+ZDMFDP(JL,JK)
      zdmfup(jl,jk)=pmflxr(jl,jk+1)+pmflxs(jl,jk+1)&
        & -pmflxr(jl,jk)-pmflxs(jl,jk)
      zdmfdp(jl,jk)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

! avoid negative humidities at ddraught top
!$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ik, jk)
DO jl=kidia,kfdia
  IF(llddraf(jl)) THEN
    jk=idtop(jl)
    ik=MIN(jk+1,klev)
    IF(zmfdq(jl,jk)<0.3_JPRB*zmfdq(jl,ik)) THEN
      IF(rmfsoltq==0.0_JPRB) THEN
        zmfdq(jl,jk)=0.3_JPRB*zmfdq(jl,ik)
      ELSE
        pmfd(jl,jk)=0.3_JPRB*pmfd(jl,ik)
      ENDIF
    ENDIF
  ENDIF
ENDDO

! avoid negative humidities near cloud top because gradient of precip flux
! and detrainment / liquid water flux too large
!$ACC LOOP SEQ
DO jk=ktdia+1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdz, zmfa)
  DO jl=kidia,kfdia
    IF(ldcum(jl).AND.jk>=kctop(jl)-1.AND.jk<kcbot(jl)) THEN
      ZDZ=PTSPHY*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
      zmfa=zmfuq(jl,jk+1)+zmfdq(jl,jk+1)-zmfuq(jl,jk)-zmfdq(jl,jk)+&
      &zmful(jl,jk+1)-zmful(jl,jk)+zdmfup(jl,jk)
      zmfa=(zmfa-plude(jl,jk))*zdz
      IF(pqen(jl,jk)+zmfa<0.0_JPRB) THEN
        plude(jl,jk)=plude(jl,jk)+2.0_JPRB*(pqen(jl,jk)+zmfa)/zdz
      ENDIF
      IF(plude(jl,jk)<0.0_JPRB) THEN
        plude(jl,jk)=0.0_JPRB
      ENDIF
    ENDIF
    IF(.NOT.LDCUM(JL)) THEN
      PMFUDE_RATE(JL,JK)=0.0_JPRB
    ENDIF
    IF(PMFD(JL,JK-1)==0.0_JPRB) THEN
      PMFDDE_RATE(JL,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

#ifndef _OPENACC
!*UPG change to operations
IF ( llconscheck ) THEN
#ifdef _OPENACC
  CALL finish('mo_cumaster:', 'llconscheck=.TRUE. not available on GPU')
#endif
  ALLOCATE(ztent(klon,klev))
  ALLOCATE(ztenq(klon,klev))
  DO jk=ktdia+1,klev
    DO jl=kidia,kfdia
      IF ( ldcum(jl) ) THEN
        ztent(jl,jk)=ptent(jl,jk)
        ztenq(jl,jk)=ptenq(jl,jk)
        ztenu(jl,jk)=ptenu(jl,jk)
        ztenv(jl,jk)=ptenv(jl,jk)
      ENDIF
    ENDDO
  ENDDO

   IF ( lmftrac .AND. ktrac>0 ) THEN
     ! this should only be possible, if PRESENT(ptenrhoc),
     ! which is not the case for .NOT. lart, so this kernel is not ported to GPU
     ALLOCATE(ztenrhoc(klon,klev,ktrac))
     ALLOCATE(zsumc(klon,4+ktrac))
     DO jn=1,ktrac
       DO jk=ktdia+1,klev
         DO jl=kidia,kfdia
           IF ( ldcum(jl) ) THEN
             ztenrhoc(jl,jk,jn)=ptenrhoc(jn)%ptr(jl,jk)
           ENDIF
         ENDDO
       ENDDO
     ENDDO
   ELSE
     ALLOCATE(zsumc(klon,4))
   ENDIF
ENDIF
!*UPG change to operations
#endif

! Calculation of kinetic energy production by the convective buoyant heat flux:
!$ACC LOOP SEQ
DO jk=ktdia+1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zcvfl_s, zcvfl_q)
  DO jl=kidia,kfdia
    IF ( ldcum(jl) ) THEN
      zcvfl_s  =  zmfus (jl,jk) + zmfds (jl,jk)
      zcvfl_q  =  zmfuq (jl,jk) + zmfdq (jl,jk)

      pdtke_con(jl,jk) = MAX( 0.0_JPRB, rg*rd/paph(jl,jk) * ( (1.0_JPRB+retv*zqenh(jl,jk))*zcvfl_s/rcpd + &
                                                              retv*ztenh(jl,jk)*zcvfl_q ) )
    ENDIF
  ENDDO
ENDDO

!----------------------------------------------------------------------

!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------

! save moisture tendency prior to update
!$ACC LOOP SEQ
DO jk=1,klev
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    ztenq_sv(jl,jk) = ptenq(jl,jk)
  ENDDO
ENDDO

IF( rmfsoltq>0.0_JPRB) THEN
! derive draught properties for implicit

  !$ACC LOOP SEQ
  DO jk=klev,ktdia+1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zmfa)
    DO jl=kidia,kfdia
      IF(ldcum(jl)) THEN
        IF(jk>kcbot(jl)) THEN
          zmfa=1.0_JPRB/MAX(1.e-15_JPRB,pmfu(jl,jk))
          pqu(jl,jk)=zqenh(jl,jk)+zmfuq(jl,jk)*zmfa
          ptu(jl,jk)=ztenh(jl,jk)+zmfus(jl,jk)*zmfa*zorcpd
          zmfus(jl,jk)=pmfu(jl,jk)*(rcpd*ptu(jl,jk)+pgeoh(jl,jk))
          zmfuq(jl,jk)=pmfu(jl,jk)*pqu(jl,jk)
          IF(llddraf(jl)) THEN
            zmfa=1.0_JPRB/MIN(-1.e-15_JPRB,pmfd(jl,jk))
            zqd(jl,jk)=zqenh(jl,jk)+zmfdq(jl,jk)*zmfa
            ztd(jl,jk)=ztenh(jl,jk)+zmfds(jl,jk)*zmfa*zorcpd
            zmfdq(jl,jk)=pmfd(jl,jk)*zqd(jl,jk)
            zmfds(jl,jk)=pmfd(jl,jk)*(rcpd*ztd(jl,jk)+pgeoh(jl,jk))
          ENDIF
        ELSEIF(jk<=kcbot(jl).AND.jk>=kctop(jl)) THEN
          zmfus(jl,jk)=pmfu(jl,jk)*(rcpd*ptu(jl,jk)+pgeoh(jl,jk))
          zmfuq(jl,jk)=pmfu(jl,jk)*pqu(jl,jk)
          zmfds(jl,jk)=pmfd(jl,jk)*(rcpd*ztd(jl,jk)+pgeoh(jl,jk))
          zmfdq(jl,jk)=pmfd(jl,jk)*zqd(jl,jk)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  
ENDIF

!$ACC END PARALLEL

CALL cudtdqn &
  & ( kidia,    kfdia,    klon,   ktdia,    klev,&
  & itopm2,   ktype,    kctop,    idtop,    ldcum,    llddraf,   ptsphy,&
  & paph,     pgeoh,    pgeo,&
  & zdph,                    &
  & pten,     ztenh,    pqen,     zqenh,    pqsen,&
  & zlglac,   plude,    psnde,    pmfu,     pmfd,&
  & zmfus,    zmfds,    zmfuq,    zmfdq,&
  & zmful,    zdmfup,   zdpmel,&
  & ptent,    ptenq,    penth,    lacc )


!----------------------------------------------------------------------

!*    9.0          COMPUTE MOMENTUM IN UPDRAUGHT AND DOWNDRAUGHT
!                  ---------------------------------------------

IF(lmfdudv) THEN

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

  !$ACC LOOP SEQ
  DO jk=klev-1,ktdia+1,-1
    ik=jk+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, zfac, zerate, zderate, zmfa)
    DO jl=kidia,kfdia
      IF(ldcum(jl)) THEN
        IF(jk==kcbot(jl).AND.ktype(jl)<3) THEN
          ikb=idpl(jl)
          zuu(jl,jk)=puen(jl,ikb-1)
          zvu(jl,jk)=pven(jl,ikb-1)
        ELSEIF(jk==kcbot(jl).AND.ktype(jl)==3) THEN
          zuu(jl,jk)=puen(jl,jk-1)
          zvu(jl,jk)=pven(jl,jk-1)
        ENDIF
        IF( jk<kcbot(jl).AND.jk>=kctop(jl)) THEN
          zfac=0.0_JPRB
! ** suggestion by P. Bechtold (40r3) - delete following to lines to improve momentum transport **
          IF(ktype(jl)==1.OR.ktype(jl)==3) zfac=2.0_JPRB
          IF(ktype(jl)==1.AND.jk<=kctop(jl)+2) zfac=3.0_jprb
          zerate=pmfu(jl,jk)-pmfu(jl,ik)+(1.0_JPRB+zfac)*pmfude_rate(jl,jk)
          zderate=(1.0_JPRB+zfac)*pmfude_rate(jl,jk)
          zmfa=1.0_JPRB/MAX(rmfcmin,pmfu(jl,jk))
          zuu(jl,jk)=(zuu(jl,ik)*pmfu(jl,ik)+zerate*puen(jl,jk)-zderate*zuu(jl,ik))*zmfa
          zvu(jl,jk)=(zvu(jl,ik)*pmfu(jl,ik)+zerate*pven(jl,jk)-zderate*zvu(jl,ik))*zmfa
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  !$ACC LOOP SEQ
  DO jk=ktdia+2,klev
    ik=jk-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zerate, zmfa)
    DO jl=kidia,kfdia
      IF( ldcum(jl)) THEN
        IF(jk==idtop(jl)) THEN
          zud(jl,jk)=0.5_JPRB*(zuu(jl,jk)+puen(jl,ik))
          zvd(jl,jk)=0.5_JPRB*(zvu(jl,jk)+pven(jl,ik))
        ELSEIF(jk>idtop(jl)) THEN
          zerate=-pmfd(jl,jk)+pmfd(jl,ik)+pmfdde_rate(jl,jk)
          zmfa=1.0_JPRB/MIN(-rmfcmin,pmfd(jl,jk))
          zud(jl,jk)=(zud(jl,ik)*pmfd(jl,ik)-zerate*puen(jl,ik)&
            & + pmfdde_rate(jl,jk)*zud(jl,ik))*zmfa
          zvd(jl,jk)=(zvd(jl,ik)*pmfd(jl,ik)-zerate*pven(jl,ik)&
            & +pmfdde_rate(jl,jk)*zvd(jl,ik))*zmfa
        ENDIF
      ENDIF
    ENDDO
  ENDDO

!*    9.1          UPDATE TENDENCIES FOR U AND V IN SUBROUTINE CUDUDV
!                  --------------------------------------------------

! for explicit/semi-implicit rescale massfluxes for stability in Momentum
!------------------------------------------------------------------------

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    zmfs(jl)=1.0_JPRB
  ENDDO

! IF(RMFSOLUV<=0.5_JPRB) THEN
  IF(rmfsoluv<=1.0_JPRB) THEN
    !$ACC LOOP SEQ
    DO jk=ktdia+1,klev
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zmfmax)
      DO jl=kidia,kfdia
        IF(ldcum(jl).AND.jk>=kctop(jl)-1) THEN
          zmfmax=(paph(jl,jk)-paph(jl,jk-1))*zcons
          IF(pmfu(jl,jk)>zmfmax.AND.jk>=kctop(jl)) &
           & zmfs(jl)=MIN(zmfs(jl),zmfmax/pmfu(jl,jk))  
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  !$ACC LOOP SEQ
  DO jk=ktdia,klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      zmfuus(jl,jk)=pmfu(jl,jk)
      zmfdus(jl,jk)=pmfd(jl,jk)
      IF(ldcum(jl).AND.jk>=kctop(jl)-1) THEN
        zmfuus(jl,jk)=pmfu(jl,jk)*zmfs(jl)
        zmfdus(jl,jk)=pmfd(jl,jk)*zmfs(jl)
      ENDIF
    ENDDO
  ENDDO

! recompute Draught properties below for Implicit
! based on linear flux profiles

  IF(rmfsoluv>0.0_JPRB) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jk, ik)
    DO jl=kidia,kfdia
      IF(ldcum(jl)) THEN
        jk=kcbot(jl)
        ik=jk-1
        zmfuub(jl)=zmfuus(jl,jk)*(zuu(jl,jk)-puen(jl,ik))
        zmfuvb(jl)=zmfuus(jl,jk)*(zvu(jl,jk)-pven(jl,ik))
      ENDIF
    ENDDO

    !$ACC LOOP SEQ
    DO jk=ktdia+1,klev
      ik=jk-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ikb, zdz, zmfa)
      DO jl=kidia,kfdia
        IF ( ldcum(jl).AND.jk>kcbot(jl) ) THEN
          ikb=kcbot(jl)
          zdz=((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
          IF(ktype(jl) == 3) THEN
            zdz=zdz*zdz
          ENDIF
          zmfa=1.0_JPRB/MAX(rmfcmin,zmfuus(jl,jk))
          zuu(jl,jk)=puen(jl,ik)+zmfuub(jl)*zdz*zmfa
          zvu(jl,jk)=pven(jl,ik)+zmfuvb(jl)*zdz*zmfa

          zmfdus(jl,jk)=zmfdus(jl,ikb)*zdz
          zud(jl,jk)=puen(jl,ik)+zud(jl,ikb)-puen(jl,ikb-1)
          zvd(jl,jk)=pven(jl,ik)+zvd(jl,ikb)-pven(jl,ikb-1)

        ENDIF
    ! add UV perturb to correct wind bias
        IF ( ldcum(jl).AND.jk>=kctop(jl) ) THEN
          zuu(jl,jk)=zuu(jl,jk)-MIN(ABS(zuu(jl,jk)),ruvper)*SIGN(1.0_JPRB,zuu(jl,jk))
          zvu(jl,jk)=zvu(jl,jk)-MIN(ABS(zvu(jl,jk)),ruvper)*SIGN(1.0_JPRB,zvu(jl,jk))
        ENDIF
      ENDDO
    ENDDO
  ENDIF


! Maximum possible convective gust
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zmaxkined)
  DO jl = kidia, kfdia
    zmaxkined   = MAX (zkined(jl,klev-2), zkined(jl,klev-1), zkined(jl,klev))
    pvddraf(jl) = SQRT( 2._jprb*zmaxkined)
    pvddraf(jl) = MIN( pvddraf(jl), conv_gust_max)
  ENDDO
  !$ACC END PARALLEL

!-------------------------------------------------------------------
! End
! Intermediate Solution for stability in EPS: 
! For original code replace line
!  &, PUEN,     PVEN,     ZMFUUS,   ZMFDUS &
!by
!  &, PUEN,     PVEN,     PMFU,     PMFD

  CALL cududv &
    & ( kidia,  kfdia,    klon,     ktdia,    klev,&
    & itopm2,   ktype,    kcbot,    kctop,    ldcum,    ptsphy,&
    & zdph,                                          &
    & paph,     puen,     pven,     zmfuus,   zmfdus,&
    & zuu,      zud,      zvu,      zvd,&
    & ptenu,    ptenv,    lacc     )
  
  IF(LMFUVDIS) THEN
! add KE dissipation

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO JL=KIDIA,KFDIA
      ZSUM12(JL)=0.0_JPRB
      ZSUM22(JL)=0.0_JPRB
    ENDDO

    !$ACC LOOP SEQ
    DO JK=ktdia,KLEV
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdz, zduten, zdvten)
      DO JL=KIDIA,KFDIA
        ZUV2(JL,JK)=0.0_JPRB
        IF (LDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
          ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))
          ZDUTEN=PTENU(JL,JK)-ZTENU(JL,JK)
          ZDVTEN=PTENV(JL,JK)-ZTENV(JL,JK)
          ZUV2(JL,JK)=SQRT(ZDUTEN**2+ZDVTEN**2)
          ZSUM22(JL)=ZSUM22(JL)+ZUV2(JL,JK)*ZDZ
          ZSUM12(JL)=ZSUM12(JL)-(PUEN(JL,JK)*ZDUTEN+PVEN(JL,JK)*ZDVTEN)*ZDZ
        ENDIF
      ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO JK=ktdia,KLEV
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdz, ztdis)
      DO JL=KIDIA,KFDIA
        IF (LDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
          ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))
          ZTDIS=ZORCPD*ZSUM12(JL)*&
               &             ZUV2(JL,JK)/MAX(1.E-15_JPRB,ZSUM22(JL))
          PTENT(JL,JK)=PTENT(JL,JK)+ZTDIS
        ENDIF
      ENDDO
    ENDDO

    !$ACC END PARALLEL
  ENDIF

ENDIF

!----------------------------------------------------------------------

!*   10.           IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
!                  NEED TO SET SOME VARIABLES A POSTERIORI TO ZERO
!                  ---------------------------------------------------

IF (.NOT.phy_params%lmfscv .OR. .NOT.phy_params%lmfpen) THEN
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

  !$ACC LOOP SEQ
  DO jk=ktdia+1,klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      IF(llo2(jl).AND.jk>=kctop(jl)-1) THEN
        ptu(jl,jk)  =pten(jl,jk)
        pqu(jl,jk)  =pqen(jl,jk)
        plu(jl,jk)  =0.0_JPRB
        penth(jl,jk) =0.0_JPRB
        pmfude_rate(jl,jk) =0.0_JPRB
        pmfdde_rate(jl,jk) =0.0_JPRB
      ENDIF
    ENDDO
  ENDDO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    IF(llo2(jl)) THEN
      kctop(jl)=klev-1
      kcbot(jl)=klev-1
    ENDIF
  ENDDO

  !$ACC END PARALLEL
ENDIF

!----------------------------------------------------------------------

!*   11.0          CHEMICAL TRACER TRANSPORT
!                  -------------------------


IF ( lmftrac .AND. ktrac>0 ) THEN

!US this is only the case for lart, which is not considered yet for GPUs

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
  !$ACC LOOP GANG VECTOR
  DO jl = 1, klon
    zmfs(jl)=1.0_JPRB
  END DO
  !$ACC END PARALLEL

  ! transport switched off for mid-level convection
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO jl=kidia,kfdia
    !IF( LDCUM(JL).AND.KTYPE(JL)/=3 ) THEN
    IF( ldcum(jl).AND.ktype(jl)/=3.AND.kcbot(jl)-kctop(jl)>=1 ) THEN
      lldcum(jl)=.TRUE.
      llddraf3(jl)=llddraf(jl)
    ELSE
      lldcum(jl)=.FALSE.
      llddraf3(jl)=.FALSE.
    ENDIF
  ENDDO

  ! check and correct mass fluxes for CFL criterium
  IF(rmfsolct<=3.0_JPRB) THEN
    !$ACC LOOP SEQ
    DO jk=ktdia+1,klev
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zmfmax)
      DO jl=kidia,kfdia
        IF(lldcum(jl).AND.jk>=kctop(jl)) THEN
          zmfmax=(paph(jl,jk)-paph(jl,jk-1))*0.8_JPRB*zcons
          IF(pmfu(jl,jk)>zmfmax) &
            & zmfs(jl)=MIN(zmfs(jl),zmfmax/pmfu(jl,jk))
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !$ACC LOOP SEQ
  DO jk=ktdia,klev
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      IF(lldcum(jl).AND.jk>=kctop(jl)-1) THEN
        zmfuus(jl,jk)=pmfu(jl,jk)*zmfs(jl)
        zmfudr(jl,jk)=pmfude_rate(jl,jk)*zmfs(jl)
      ELSE
        zmfuus(jl,jk)=0._jprb
        zmfudr(jl,jk)=0._jprb
      ENDIF
      IF ( llddraf3(jl) .AND. jk>=idtop(jl)-1) THEN
        zmfdus(jl,jk)=pmfd(jl,jk)*zmfs(jl)
        zmfddr(jl,jk)=pmfdde_rate(jl,jk)*zmfs(jl)
      ELSE
        zmfdus(jl,jk)=0._jprb
        zmfddr(jl,jk)=0._jprb
      ENDIF
    ENDDO
  ENDDO

  IF( lmfsmooth ) THEN
    ! smmoothing of mass fluxes (gradients) at top and bottom of draughts
    !$ACC LOOP SEQ
    DO jk=ktdia+1,klev-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zerate)
      DO jl=kidia,kfdia
        IF(llddraf3(jl).AND.zmfdus(jl,jk)<0.0_JPRB .AND. zmfdus(jl,jk+1)==0.0_JPRB) THEN
          zerate=MIN(0._jprb,zmfdus(jl,jk)-0.5_JPRB*zmfdus(jl,jk-1))
          zmfdus(jl,jk)=zmfdus(jl,jk)-zerate
          zmfddr(jl,jk)=zmfddr(jl,jk)-zerate
          zmfddr(jl,jk+1)=-zmfdus(jl,jk)
        ENDIF
        IF(lldcum(jl).AND.jk==kctop(jl)) THEN
          zerate=MAX(0.0_JPRB,zmfuus(jl,jk)-0.5_JPRB*zmfuus(jl,jk+1))
          zmfuus(jl,jk)=zmfuus(jl,jk)-zerate
          zmfudr(jl,jk)=zmfudr(jl,jk)+zerate
          zmfudr(jl,jk-1)=zmfuus(jl,jk)
        ENDIF
      ENDDO
    ENDDO
    !$ACC LOOP SEQ
    DO jk=klev-1,ktdia+1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        IF(lldcum(jl)) THEN
          IF(zmfudr(jl,jk)==0.0_JPRB.AND.zmfudr(jl,jk-1)>0.0_JPRB) THEN
            zmfudr(jl,jk)=0.5_JPRB*zmfudr(jl,jk-1)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  !$ACC END PARALLEL

  IF ( ktrac > 0 ) THEN
    CALL cuctracer &
      & ( kidia,    kfdia,    klon,  ktdia,  klev,     ktrac,&
      &   kctop,     idtop,&
      & lldcum,   llddraf3,  ptsphy,  &
      & paph,     zdph,     zdgeoh,          &
      & zmfuus,   zmfdus,   zmfudr,   zmfddr,&
      & pcen,     ptenrhoc, lacc = .TRUE. )
  ENDIF
ENDIF

!----------------------------------------------------------------------

!*   12.           PUT DETRAINMENT RATES FROM MFLX UNITS IN UNITS MFLX/M 
!                  FOR ERA40
!                  ---------------------------------------------------

!$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
!$ACC LOOP GANG VECTOR PRIVATE(zro) COLLAPSE(2)
DO jk=ktdia+1,klev
  DO jl=kidia,kfdia
    IF ( ldcum(jl) ) THEN
      zro=rg/(pgeoh(jl,jk)-pgeoh(jl,jk+1))  ! 1/dz
      pmfude_rate(jl,jk)=pmfude_rate(jl,jk)*zro
      pmfdde_rate(jl,jk)=pmfdde_rate(jl,jk)*zro
      IF(jk<kctop(jl)) THEN
        plu(jl,jk)=0.0_JPRB
        ptu(jl,jk)=pten(jl,jk)
        pqu(jl,jk)=pqen(jl,jk)
      ENDIF
    ENDIF
  ENDDO
ENDDO
!$ACC END PARALLEL

!----------------------------------------------------------------------
!*UPG change to operations

IF ( llconscheck ) THEN

!US is set to .FALSE. above, so is also not considered for GPUs

!*   13.0          CONSERVATION CHECK and CORRECTION
!                  ---------------------------------

  DO jl=kidia,kfdia
    zsumc(jl,:)=0._jprb
  ENDDO
  DO jk=klev,ktdia+1,-1
    DO jl=kidia,kfdia
      IF ( ldcum(jl) .AND. jk>=kctop(jl)-1 ) THEN
        ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))/RG
        zsumc(jl,1)=zsumc(jl,1)+(ptenq(jl,jk)-ztenq(jl,jk))*zdz+plude(jl,jk)
        zalv=foelhmcu(pten(jl,jk))
        zsumc(jl,2)=zsumc(jl,2)+rcpd*(ptent(jl,jk)-ztent(jl,jk))*zdz-zalv*plude(jl,jk)
        zsumc(jl,3)=zsumc(jl,3)+(ptenu(jl,jk)-ztenu(jl,jk))*zdz
        zsumc(jl,4)=zsumc(jl,4)+(ptenv(jl,jk)-ztenv(jl,jk))*zdz
      ENDIF
    ENDDO
  ENDDO
  IF ( lmftrac .AND. ktrac>0 ) THEN
    DO jn=1,ktrac
      DO jk=klev,ktdia+1,-1
        DO jl=kidia,kfdia
          IF ( ldcum(jl) .AND. jk>=kctop(jl)-1) THEN
       !DR     zdz=(paph(jl,jk+1)-paph(jl,jk))/rg
       !DR     zsumc(jl,4+jn)=zsumc(jl,4+jn)+(ptenc(jl,jk,jn)-ztenc(jl,jk,jn))*zdz
            zdz=zdgeoh(jl,jk)/rg  ! dz
            zsumc(jl,4+jn)=zsumc(jl,4+jn)+(ptenrhoc(jn)%ptr(jl,jk)-ztenrhoc(jl,jk,jn))*zdz
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  DO jl=kidia,kfdia
    IF ( ldcum(jl) ) THEN
      zalv=foelhmcu(pten(jl,klev))
      zsfl(jl)=pmflxr(jl,klev+1)+pmflxs(jl,klev+1)

      WRITE(61,'(i4,a9,2f15.8,i4,a9,f15.8,a10,2f15.8)')jl,' CONS q: ',&
       & -zsumc(jl,1)*zalv,zsfl(jl)*zalv,ktype(jl),&
       & ' CONS h: ',zsumc(jl,2),' CONS uv: ',zsumc(jl,3),zsumc(jl,4)

      ikb=kctop(jl)
      zdz=(paph(jl,klev+1)-paph(jl,ikb-1))/rg
      zsumc(jl,1)=(zsumc(jl,1)+zsfl(jl))/zdz
      zsumc(jl,2)=(zsumc(jl,2)-zalv*zsfl(jl))/(zdz*rcpd)
    ENDIF
  END DO

  DEALLOCATE(zsumc)
  IF ( lmftrac .AND. ktrac>0 ) THEN
    DEALLOCATE(ztenrhoc)
  ENDIF
  DEALLOCATE(ztenq)
  DEALLOCATE(ztent)

ENDIF

!----------------------------------------------------------------------

!*    14.0         COMPUTE CONVECTIVE TENDENCIES FOR LIQUID AND SOLID
!                  CLOUD CONDENSATE, CHANGE PRECIP UNITS IN M/S (if wanted, K. Froehlich, 19.01.2009)
!                  --------------------------------------------------


  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO jk=ktdia,klev
     DO jl=kidia,kfdia
       ptenrhoq(jl,jk)= (ptenq(jl,jk) - ztenq_sv(jl,jk)) &
         &            * (zdph(jl,jk)/zdgeoh(jl,jk))
       ptenrhol(jl,jk)=plude(jl,jk)*rg/zdgeoh(jl,jk)
       ptenrhoi(jl,jk)=(1.0_JPRB-foealfcu(pten(jl,jk)))*ptenrhol(jl,jk)
       ptenrhol(jl,jk)= ptenrhol(jl,jk)-ptenrhoi(jl,jk)
    ENDDO
  ENDDO
  !$ACC END PARALLEL

  IF (phy_params%lmfdsnow) THEN
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
    !$ACC LOOP GANG VECTOR PRIVATE(zdz) COLLAPSE(2)
    DO jk=ktdia,klev
      DO jl=kidia,kfdia
        zdz = rg/zdgeoh(jl,jk)
        ptenrhos(jl,jk)=psnde(jl,jk,1)*zdz
        ptenrhor(jl,jk)=psnde(jl,jk,2)*zdz
      ENDDO
    ENDDO
    !$ACC END PARALLEL
  ENDIF

!----------------------------------------------------------------------

!*    15.0         COMPUTE LIGHTENING POTENTIAL INDEX AND LIGHTENING 
!                  FLASH DENSITY BASED ON UPDRAFT PROFILE 
!                  --------------------------------------------------

  IF (PRESENT(l_lpi)) THEN
    IF (l_lpi) THEN
      CALL cucalclpi(klon, klev, ktype, ptu, plu, zkineu, pmflxs        &
      &          , pten, pap, zdgeoh, ldland, lpi, lacc)

!! alternative commented out - in the alternative the mass flux
!! is used to computed the LPI - which is tricky, as the area fraction
!! of the updraft is unknown
!!  CALL cucalclpi(klon, klev, ktype, ptu, plu,                           &
!!  &              2*(pmfu*pten*Rd/(pap+1E-10_jprb))**2, pmflxs        &
!!  &          , pten, pap, zdgeoh, ldland, lpi, lacc)
      CALL cucalcmlpi(klon, klev, lpi, pten, pqen, pap, paph, koi, mlpi, lacc)
    ENDIF
  ENDIF

  IF (PRESENT(l_lfd)) THEN
    IF (l_lfd) THEN
      CALL cucalclfd(klon, klev, ktype, ptu, plu, kcbot, pcape, pmflxs &
      &          , pten, pap, zdgeoh, pgeoh, ldland, lfd, lacc)
    ENDIF
  ENDIF

!  DO JL=KIDIA,KFDIA
!   PMFLXR(JL,KLEV+1)=PMFLXR(JL,KLEV+1)*1.E-3
!   PMFLXS(JL,KLEV+1)=PMFLXS(JL,KLEV+1)*1.E-3
!  ENDDO
!----------------------------------------------------------------------
!*UPG Change to operations

IF (lhook) CALL dr_hook('CUMASTRN',1,zhook_handle)

!end for l_lfd
!$ACC END DATA

!end for l_lpi
!$ACC END DATA

!$ACC WAIT(1)
!$ACC END DATA

END SUBROUTINE cumastrn

END MODULE mo_cumaster

