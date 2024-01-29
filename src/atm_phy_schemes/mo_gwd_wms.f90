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
!     THIS ROUTINE COMPUTES NON-OROGRAPHIC GRAVITY WAVE DRAG
!     AFTER SCINOCCA (2003) AND Mc LANDRESS AND SCINOCCIA (JAS 2005)
!     HYDROSTATIC NON-ROTATIONAL SIMPLIFIED VERSION OF THE
!     WARNER AND MCINTYRE (1996) NON-OROGRAPHIC GRAVITY WAVE PARAMETERIZATION
!     CONSTANTS HAVE BEEN OPTIMIZED FOLLOWING M. ERN ET AL. (ATMOS. CHEM. PHYS. 2006)
!
!     REFERENCE: Orr, A., P. Bechtold, J. Scinoccia, M. Ern, M. Janiskova, 2010:
!                Improved middle atmosphere climate and analysis in the ECMWF forecasting system
!                through a non-orographic gravity wave parametrization. J.  Climate., 23, 5905-5926.
!
!     LAUNCH SPECTRUM - GENERALIZED DESAUBIES
!     INCLUDES A CRITICAL-LEVEL CORRECTION THAT PREVENTS THE
!     MOMEMTUM DEPOSITION IN EACH AZIMUTH FROM DRIVING THE FLOW TO SPEEDS FASTER
!     THAN THE PHASE SPEED OF THE WAVES, I.E. WHEN WAVES BREAK THEY DRAG THE MEAN
!     FLOW TOWARDS THEIR PHASE SPEED - NOT PAST IT.
!

MODULE mo_gwd_wms
  
  USE mo_kind,               ONLY: jprb=>wp, vp, jpim =>i4
  !------------------------------------------------------------------------------
  USE mo_math_constants , ONLY : rpi => pi
  
  USE mo_physical_constants , ONLY :   &
    & rd           , & ! gas constant for dry air
    & rcpd => cpd  , & ! specific heat capacity at constant pressure
    & rg   => grav     ! acceleration due to gravity
  
  USE data_gwd,    ONLY : gcstar, gptwo, nslope, gfluxlaun, & ! nlaunch, &
    & ggaussa, ggaussb, ngauss, gcoeff, lozpr
  
  IMPLICIT NONE
  
  PRIVATE
  


  PUBLIC :: gwdrag_wms
  
 
CONTAINS
  
  SUBROUTINE gwdrag_wms(kidia,  kfdia,   klon,  klev, klaunch, ptstep, &
    & ptm1 ,  pum1,    pvm1,  papm1,  paphm1, pgeo1, &
    & pgelat, pprecip, &
    & ptenu, ptenv,   pfluxu, pfluxv, &
    & lacc)
    
   !<**** *GWDRAG_WMS*
   !!
   !!     Original Fortran Code by   J. SCINOCCIA
   !!     Rewritten in IFS format by A. ORR          E.C.M.W.F.     August 2008
   !! 
   !!     MODIFICATIONS
   !!     -------------
   !!           October 2008 : Cleaning and full adaptation   P. Bechtold/JJMorcrette
   !!                          Optimisation+bug corrections   P. Bechtold
   !! ---------------------------------------------------------------------------------
    
    
    !USE PARKIND1,    ONLY : JPIM, JPRB
    !USE YOEGWWMS,    ONLY : GCSTAR, GPTWO, NSLOPE, GFLUXLAUN, klaunch, &
    !                      & GGAUSSA, GGAUSSB, NGAUSS, GCOEFF, LOZPR
    !USE YOMCST,      ONLY : RG, RD, RCPD, RPI
    !USE YOMHOOK,     ONLY : LHOOK,   DR_HOOK
    !USE YOEGWD,      ONLY : GSSEC
    
    IMPLICIT NONE
    
    !in
    INTEGER(KIND=jpim),INTENT(in) :: kidia, kfdia
    INTEGER(KIND=jpim),INTENT(in) :: klon         ! horizontal dimension
    INTEGER(KIND=jpim),INTENT(in) :: klev         ! vertical levels
    INTEGER(KIND=jpim),INTENT(in) :: klaunch      ! launch level
    REAL(KIND=jprb),INTENT(in) :: ptstep      !model time step
    REAL(KIND=jprb),INTENT(in) :: pum1(:,:)   ! (klon,klev)   full model level zonal velocity (t-dt)
    REAL(KIND=jprb),INTENT(in) :: pvm1(:,:)   ! (klon,klev)   full model level meridional velocity (t-dt)
    REAL(KIND=jprb),INTENT(in) :: ptm1(:,:)   ! (klon,klev)   full model level temperature (t-dt)
    REAL(KIND=jprb),INTENT(in) :: papm1(:,:)  ! (klon,klev)   full model level pressure (t-dt)
    REAL(KIND=jprb),INTENT(in) :: paphm1(:,:) ! (klon,klev+1) half-model level pressure (t-dt)
    REAL(KIND=jprb),INTENT(in) :: pgeo1(:,:)  ! (klon,klev)   full model level geopotential
    REAL(KIND=jprb),INTENT(in) :: pgelat(:)   ! (klon)        latitude
    REAL(KIND=jprb),INTENT(in) :: pprecip(:)  ! (klon)        total surface precipitation

    !inout
    REAL(KIND=vp),  INTENT(out):: ptenu(:,:)  ! (klon,klev)   full-model level zonal momentum tendency
    REAL(KIND=vp),  INTENT(out):: ptenv(:,:)  ! (klon,klev)   full-model level meridional momentum tendency
    REAL(KIND=jprb),INTENT(out):: pfluxu(:,:) ! (klon,klev+1) zonal component of vertical momentum flux (Pa)
    REAL(KIND=jprb),INTENT(out):: pfluxv(:,:) ! (klon,klev+1) meridional component of vertical momentum flux (Pa)
    
    LOGICAL, INTENT(in) :: lacc               ! If .TRUE. use openACC. (.FALSE. is used during initialization)

    !work
    INTEGER(KIND=jpim), PARAMETER :: iazidim=4      !number of azimuths
    INTEGER(KIND=jpim), PARAMETER :: incdim=20      !number of discretized c spectral elements in launch spectrum
    REAL(KIND=jprb) :: zuhm1(klon,klev)             !half-model level zonal velocity
    REAL(KIND=jprb) :: zvhm1(klon,klev)             !half-model level meridional velocity
    REAL(KIND=jprb) :: zbvfhm1(klon,klev)           !half-model level Brunt-Vaisalla frequency
    REAL(KIND=jprb) :: zrhohm1(klon,klev)           !half-model level density
    REAL(KIND=jprb) :: zx(incdim)                   !coordinate transformation
    REAL(KIND=jprb) :: zci(incdim)                  !phase speed element
    REAL(KIND=jprb) :: zdci(incdim)
    REAL(KIND=jprb) :: zui(klon,klev,iazidim)       !intrinsic velocity
    REAL(KIND=jprb) :: zul(klon,iazidim)            !velocity in azimuthal direction at launch level
    REAL(KIND=jprb) :: zbvfl(klon)                  !buoyancy at launch level
    REAL(KIND=jprb) :: zcosang(iazidim)             !cos of azimuth angle
    REAL(KIND=jprb) :: zsinang(iazidim)             !sin of azimuth angle
    REAL(KIND=jprb) :: zfct(klon,klev)
    REAL(KIND=jprb) :: zfnorm(klon)                 !normalisation factor (A)
    REAL(KIND=jprb) :: zci_min(klon,iazidim)
    REAL(KIND=jprb) :: zthm1(klon,klev)             !temperature on half-model levels
    REAL(KIND=jprb) :: zflux(klon,incdim,iazidim)   !momentum flux at each vertical level and azimuth
    REAL(KIND=jprb) :: zpu(klon,klev,iazidim)       !momentum flux
    REAL(KIND=jprb) :: zdfl(klon,klev,iazidim)
    REAL(KIND=jprb) :: zact(klon,incdim,iazidim)    !if =1 then critical level encountered
    REAL(KIND=jprb) :: zacc(klon,incdim,iazidim)
    REAL(KIND=jprb) :: zcrt(klon,klev,iazidim)
    
    REAL(KIND=jprb) :: gssec=1.e-12_JPRB
    
    INTEGER(KIND=jpim) :: inc, jk, jl, iazi
    REAL(KIND=jprb) :: zradtodeg, zgelatdeg
    REAL(KIND=jprb) :: zcimin, zcimax
    REAL(KIND=jprb) :: zgam, zpexp, zxmax, zxmin, zxran, zdx, zx1, zx2
    REAL(KIND=jprb) :: zang, zaz_fct, znorm, zang1, ztx
    REAL(KIND=jprb) :: zu, zcin
    REAL(KIND=jprb) :: zcin4, zbvfl4, zcin2, zbvfl2, zcin3, zbvfl3, zcinc
    REAL(KIND=jprb) :: zatmp, zfluxs, zdep, zfluxsq, zulm, zdft, ze1, ze2
    REAL(KIND=jprb) :: zms_l,zms, z0p5, z0p0, z50s
    REAL(KIND=jprb) :: zgauss(klon), zfluxlaun(klon), zcngl(klon)
    REAL(KIND=jprb) :: zcons1,zcons2,zdelp,zrgpts
    
!     REAL(KIND=jprb) :: zhook_handle
    
    !--------------------------------------------------------------------------
    
    !IF (LHOOK) CALL DR_HOOK('GWDRAG_WMS',0,ZHOOK_HANDLE)
    
    !--------------------------------------------------------------------------
    
    !$ACC DATA PRESENT(pum1, pvm1, ptm1, papm1, paphm1, pgeo1, pgelat, pprecip, ptenu, ptenv, pfluxu, pfluxv) &
    !$ACC   CREATE(zfnorm, zfct, zacc, zui, zcngl, zrhohm1, zsinang, zthm1, zul, zgauss, zcosang, zdci, zci, zcrt) &
    !$ACC   CREATE(zbvfl, zx, zvhm1, zact, zbvfhm1, zuhm1, zdfl, zpu, zflux, zfluxlaun, zci_min) &
    !$ACC   IF(lacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)

    !*       INPUT PARAMETERS
    !*       ----------------
    
    zradtodeg=57.29577951_JPRB
    
    !m_star
    zms_l=2000.0_JPRB
    zms=2._jprb*rpi/zms_l
    
    
    !*       INITIALIZE FIELDS TO ZERO
    !*       -------------------------
    
    !$ACC LOOP SEQ
    DO jk=1,klev
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        ptenu(jl,jk)=0.0_JPRB
        ptenv(jl,jk)=0.0_JPRB
      ENDDO
    ENDDO
    
    !$ACC LOOP SEQ
    DO iazi=1,iazidim
      !$ACC LOOP SEQ
      DO jk=1,klev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jl=kidia,kfdia
          zpu(jl,jk,iazi)=0.0_JPRB
          zcrt(jl,jk,iazi)=0.0_JPRB
          zdfl(jl,jk,iazi)=0.0_JPRB
        ENDDO
      ENDDO
    ENDDO
    
    !$ACC LOOP SEQ
    DO jk=1,klev+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        pfluxu(jl,jk)=0.0_JPRB
        pfluxv(jl,jk)=0.0_JPRB
      ENDDO
    ENDDO
    
    
    !*       INITIALIZE PARAMETERS FOR COORDINATE TRANSFORM
    !*       ----------------------------------------------
    
    ! ZCIMIN,ZCIMAX - min,max intrinsic launch-level phase speed (c-U_o) (m/s)
    ! ZGAM - half=width of coordinate stretch
    
    zcimin=0.25_JPRB
    zcimax=100.0_JPRB
    zgam=0.25_JPRB
    
    zpexp=gptwo/2.0_JPRB
    
    ! set initial min ci in each column and azimuth (used for critical levels)
    
    !$ACC LOOP SEQ
    DO iazi=1,iazidim
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zci_min(jl,iazi)=zcimin
      ENDDO
    ENDDO
    
    
    !*       DEFINE HALF MODEL LEVEL WINDS AND TEMPERATURE
    !*       -----------------------------------
    
    !$ACC LOOP SEQ
    DO jk=2,klev
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zthm1(jl,jk) =0.5_JPRB*(ptm1(jl,jk-1)+ptm1(jl,jk))
        zuhm1(jl,jk) =0.5_JPRB*(pum1(jl,jk-1)+pum1(jl,jk))
        zvhm1(jl,jk) =0.5_JPRB*(pvm1(jl,jk-1)+pvm1(jl,jk))
      ENDDO
    ENDDO
    jk=1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      zthm1(jl,jk)=ptm1(jl,jk)
      zuhm1(jl,jk)=pum1(jl,jk)
      zvhm1(jl,jk)=pvm1(jl,jk)
    ENDDO
    
    
    !*       DEFINE STATIC STABILITY AND AIR DENSITY ON HALF MODEL LEVELS
    !*       ------------------------------------------------------------
    
    zcons1=1.0_JPRB/rd
    zcons2=rg**2/rcpd
    !$ACC LOOP SEQ
    DO jk=klev,2,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zrhohm1(jl,jk)=paphm1(jl,jk)*zcons1/zthm1(jl,jk)
        zbvfhm1(jl,jk)=zcons2/zthm1(jl,jk)*(1.0_JPRB+rcpd           &
          *(ptm1(jl,jk)-ptm1(jl,jk-1))/(pgeo1(jl,jk)-pgeo1(jl,jk-1)))
        zbvfhm1(jl,jk)=MAX(zbvfhm1(jl,jk),gssec)
        zbvfhm1(jl,jk)=SQRT(zbvfhm1(jl,jk))
      ENDDO
    ENDDO
    
    !*       SET UP AZIMUTH DIRECTIONS AND SOME TRIG FACTORS
    !*       -----------------------------------------------
    
    zang=2.0_JPRB*rpi/REAL(iazidim,jprb)
    zaz_fct=1.0_JPRB
    
    ! get normalization factor to ensure that the same amount of momentum
    ! flux is directed (n,s,e,w) no mater how many azimuths are selected.
    ! note, however, the code below assumes a symmetric distribution of
    ! of azimuthal directions (ie 4,8,16,32,...)
    
    znorm=0.0_JPRB
    !$ACC LOOP SEQ
    DO iazi=1,iazidim
      zang1=(REAL(iazi,jprb)-1.0_jprb)*zang
      zcosang(iazi)=COS(zang1)
      zsinang(iazi)=SIN(zang1)
      znorm=znorm+ABS(zcosang(iazi))
    ENDDO
    zaz_fct=2._jprb*zaz_fct/znorm
    
    
    !*       DEFINE COORDINATE TRANSFORM
    !*       -----------------------------------------------
    
    ! note that this is expresed in terms of the intrinsic phase speed
    ! at launch ci=c-u_o so that the transformation is identical at every
    ! launch site.
    ! See Eq. 28-30 of Scinocca 2003.
    
    zxmax=1.0_JPRB/zcimin
    zxmin=1.0_JPRB/zcimax
    
    zxran=zxmax-zxmin
    zdx=zxran/REAL((incdim-1),KIND(zdx))
    zx1=zxran/(EXP(zxran/zgam)-1.0_JPRB)
    zx2=zxmin-zx1
    
    !$ACC LOOP SEQ
    DO inc=1,incdim
     !ZTX=REAL(INC-1)*ZDX+ZXMIN
      ztx=REAL((inc-1),KIND(1._jprb))*zdx+zxmin
      zx(inc)=zx1*EXP((ztx-zxmin)/zgam)+zx2                            !Eq. 29 of Scinocca 2003
      zci(inc)=1.0_JPRB/zx(inc)                                        !Eq. 28 of Scinocca 2003
      zdci(inc)=zci(inc)**2*(zx1/zgam)*EXP((ztx-zxmin)/zgam)*zdx       !Eq. 30 of Scinocca 2003
    ENDDO
    
    
    !*       DEFINE INTRINSIC VELOCITY (RELATIVE TO LAUNCH LEVEL VELOCITY) U(Z)-U(Zo), AND COEFFICINETS
    !*       ------------------------------------------------------------------------------------------
    
    !$ACC LOOP SEQ
    DO iazi=1,iazidim
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zul(jl,iazi)=zcosang(iazi)*zuhm1(jl,klaunch)+zsinang(iazi)*zvhm1(jl,klaunch)
      ENDDO
    ENDDO
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      zbvfl(jl)=zbvfhm1(jl,klaunch)
    ENDDO
    
    !$ACC LOOP SEQ
    DO jk=2,klaunch
      !$ACC LOOP SEQ
      DO iazi=1,iazidim
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zu)
        DO jl=kidia,kfdia
          zu=zcosang(iazi)*zuhm1(jl,jk)+zsinang(iazi)*zvhm1(jl,jk)
          zui(jl,jk,iazi)=zu-zul(jl,iazi)
        ENDDO
      ENDDO
    ENDDO
    
    !*       DEFINE RHO(Zo)/N(Zo)
    !*       -------------------
    
    !$ACC LOOP SEQ
    DO jk=2,klaunch
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zfct(jl,jk)=zrhohm1(jl,jk)/zbvfhm1(jl,jk)
      ENDDO
    ENDDO
    
    !*       SET LAUNCH MOMENTUM FLUX SPECTRAL DENSITY
    !*       -----------------------------------------
    
    ! Eq. (25) of Scinocca 2003 (not including the 'A' component), and with U-Uo=0
    ! do this for only one azimuth since it is identical to all azimuths, and it will be renormalized
    
    IF(nslope==1) THEN
      ! s=1 case
      !$ACC LOOP SEQ
      DO inc=1,incdim
        zcin=zci(inc)
        zcin4=(zms*zcin)**4
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zbvfl4)
        DO jl=kidia,kfdia
          zbvfl4=zbvfl(jl)**4
          zflux(jl,inc,1)=zfct(jl,klaunch)*zbvfl4*zcin/(zbvfl4+zcin4)
          zact(jl,inc,1)=1.0_JPRB
        ENDDO
      ENDDO
    ELSEIF(nslope==-1) THEN
      ! s=-1 case
      !$ACC LOOP SEQ
      DO inc=1,incdim
        zcin=zci(inc)
        zcin2=(zms*zcin)**2
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zbvfl2)
        DO jl=kidia,kfdia
          zbvfl2=zbvfl(jl)**2
          zflux(jl,inc,1)=zfct(jl,klaunch)*zbvfl2*zcin/(zbvfl2+zcin2)
          zact(jl,inc,1)=1.0_JPRB
        ENDDO
      ENDDO
    ELSEIF(nslope==0) THEN
      ! s=0 case
      !$ACC LOOP SEQ
      DO inc=1,incdim
        zcin=zci(inc)
        zcin3=(zms*zcin)**3
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zbvfl3)
        DO jl=kidia,kfdia
          zbvfl3=zbvfl(jl)**3
          zflux(jl,inc,1)=zfct(jl,klaunch)*zbvfl3*zcin/(zbvfl3+zcin3)
          zact(jl,inc,1)=1.0_JPRB
          zacc(jl,inc,1)=1.0_JPRB
        ENDDO
      ENDDO
    ENDIF
    
    !*       NORMALIZE LAUNCH MOMENTUM FLUX
    !*       ------------------------------
    
    ! (rho x F^H = rho_o x F_p^total)
    
    ! integrate (ZFLUX x dX)
    !$ACC LOOP SEQ
    DO inc=1,incdim
      zcinc=zdci(inc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zpu(jl,klaunch,1)=zpu(jl,klaunch,1)+zflux(jl,inc,1)*zcinc
      ENDDO
    ENDDO
    
    !*       NORMALIZE GFLUXLAUN TO INCLUDE SENSITIVITY TO PRECIPITATION
    !*       -----------------------------------------------------------
    
    ! Also other options to alter tropical values
    
    ! A=ZFNORM in Scinocca 2003.  A is independent of height.
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jl=kidia,kfdia
      zfluxlaun(jl)=gfluxlaun
      zfnorm(jl)=zfluxlaun(jl)/zpu(jl,klaunch,1)
    ENDDO
    
    ! If LOZPR=TRUR then increase EPLAUNCH over tropics
    IF (lozpr) THEN
      IF (ngauss==1) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jl=kidia,kfdia
          zfluxlaun(jl)=gfluxlaun*(1.0_JPRB+MIN(0.5_JPRB,gcoeff*pprecip(jl)))     !precip
          !  ZFLUXLAUN(JL)=GFLUXLAUN*(1.0_JPRB+MIN(0.5_JPRB,1.0E-3_JPRB*PPRECIP(JL)))!cape
          zfnorm(jl)=zfluxlaun(jl)/zpu(jl,klaunch,1)
        ENDDO
      ELSEIF (ngauss==2) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zgelatdeg)
        DO jl=kidia,kfdia
          zgelatdeg=pgelat(jl)*zradtodeg
          zgauss(jl)=ggaussb*EXP((-zgelatdeg*zgelatdeg)/(2._jprb*ggaussa*ggaussa))
          zfluxlaun(jl)=(1.0_JPRB+zgauss(jl))*gfluxlaun
          zfnorm(jl)=zfluxlaun(jl)/zpu(jl,klaunch,1)
        ENDDO
      ELSEIF (ngauss==4) THEN
        ! Set latitudinal dependence to optimize stratospheric winds for 36r1
        z50s=-50.0_JPRB
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zgelatdeg)
        DO jl=kidia,kfdia
          zgelatdeg=pgelat(jl)*zradtodeg-z50s
          zgauss(jl)=ggaussb*EXP((-zgelatdeg*zgelatdeg)/(2._jprb*ggaussa*ggaussa))
          zfluxlaun(jl)=(1.0_JPRB+zgauss(jl))*gfluxlaun
          zfnorm(jl)=zfluxlaun(jl)/zpu(jl,klaunch,1)
        ENDDO
      ENDIF
    ENDIF
    
    !$ACC LOOP SEQ
    DO iazi=1,iazidim
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zpu(jl,klaunch,iazi)=zfluxlaun(jl)
      ENDDO
    ENDDO
    
    !*       ADJUST CONSTANT ZFCT
    !*       --------------------
    !$ACC LOOP SEQ
    DO jk=2,klaunch
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zfct(jl,jk)=zfnorm(jl)*zfct(jl,jk)
      ENDDO
    ENDDO
    
    !*       RENORMALIZE EACH SPECTRAL ELEMENT IN FIRST AZIMUTH
    !*       --------------------------------------------------
    !$ACC LOOP SEQ
    DO inc=1,incdim
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zflux(jl,inc,1)=zfnorm(jl)*zflux(jl,inc,1)
      ENDDO
    ENDDO
    
    !*       COPY ZFLUX INTO ALL OTHER AZIMUTHS
    !*       --------------------------------
    
    ! ZACT=1 then no critical level
    ! ZACT=0 then critical level
    
    !$ACC LOOP SEQ
    DO iazi=2,iazidim
      !$ACC LOOP SEQ
      DO inc=1,incdim
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jl=kidia,kfdia
          zflux(jl,inc,iazi)=zflux(jl,inc,1)
          zact(jl,inc,iazi)=1.0_JPRB
          zacc(jl,inc,iazi)=1.0_JPRB
        ENDDO
      ENDDO
    ENDDO
    
    ! -----------------------------------------------------------------------------
    
    !*       BEGIN MAIN LOOP OVER LEVELS
    !*       ---------------------------
    
    !* begin IAZIDIM do-loop
    !* --------------------
    
    !$ACC LOOP SEQ
    DO iazi=1,iazidim
      
      !* begin JK do-loop
      !* ----------------
      !$ACC LOOP SEQ
      DO jk=klaunch-1,2,-1
        
        
        !* first do critical levels
        !* ------------------------
        
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jl=kidia,kfdia
          zci_min(jl,iazi)=MAX(zci_min(jl,iazi),zui(jl,jk,iazi))
        ENDDO
        
        !* set ZACT to zero if critical level encountered
        !* ----------------------------------------------
        
        z0p5=0.5_JPRB
        !$ACC LOOP SEQ
        DO inc=1,incdim
          zcin=zci(inc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zatmp)
          DO jl=kidia,kfdia
            zatmp=z0p5+SIGN(z0p5,zcin-zci_min(jl,iazi))
            zacc(jl,inc,iazi)=zact(jl,inc,iazi)-zatmp
            zact(jl,inc,iazi)=zatmp
          ENDDO
        ENDDO
        
        !* integrate to get critical-level contribution to mom deposition on this level, i.e. ZACC=1
        !* ----------------------------------------------------------------------------------------
        
        !$ACC LOOP SEQ
        DO inc=1,incdim
          zcinc=zdci(inc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jl=kidia,kfdia
            zdfl(jl,jk,iazi)=zdfl(jl,jk,iazi)+&
              & zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zcinc
          ENDDO
        ENDDO
        
        !* get weighted average of phase speed in layer
        !* --------------------------------------------
        
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zatmp)
        DO jl=kidia,kfdia
          IF(zdfl(jl,jk,iazi)>0.0_JPRB) THEN
            zatmp=zcrt(jl,jk,iazi)
!$NEC unroll(incdim)
            !$ACC LOOP SEQ
            DO inc=1,incdim
              zatmp=zatmp+zci(inc)*&
                & zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zdci(inc)
            ENDDO
            zcrt(jl,jk,iazi)=zatmp/zdfl(jl,jk,iazi)
          ELSE
            zcrt(jl,jk,iazi)=zcrt(jl,jk+1,iazi)
          ENDIF
        ENDDO
        
        !* do saturation (Eq. (26) and (27) of Scinocca 2003)
        !* -------------------------------------------------
        
        IF(gptwo==3.0_JPRB) THEN
          !$ACC LOOP SEQ
          DO inc=1,incdim
            zcin=zci(inc)
            zcinc=1.0_JPRB/zcin
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ze1, ze2, zfluxsq, zdep)
            DO jl=kidia,kfdia
              ze1=zcin-zui(jl,jk,iazi)
              ze2=gcstar*zfct(jl,jk)*ze1
              zfluxsq=ze2*ze2*ze1*zcinc
              !  ZFLUXSQ=ZE2*ZE2*ZE1/ZCIN
              zdep=zact(jl,inc,iazi)*(zflux(jl,inc,iazi)**2-zfluxsq)
              IF(zdep>0.0_JPRB) THEN
                zflux(jl,inc,iazi)=SQRT(zfluxsq)
              ENDIF
            ENDDO
          ENDDO
        ELSEIF(gptwo==2.0_JPRB) THEN
          !$ACC LOOP SEQ
          DO inc=1,incdim
            zcin=zci(inc)
            zcinc=1.0_JPRB/zcin
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zfluxs, zdep)
            DO jl=kidia,kfdia
              zfluxs=gcstar*zfct(jl,jk)*&
                & (zcin-zui(jl,jk,iazi))**2*zcinc
              !  ZFLUXS=GCSTAR*ZFCT(JL,JK)*(ZCIN-ZUI(JL,JK,IAZI))**2/ZCIN
              zdep=zact(jl,inc,iazi)*(zflux(jl,inc,iazi)-zfluxs)
              IF(zdep>0.0_JPRB) THEN
                zflux(jl,inc,iazi)=zfluxs
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        
        !* integrate spectrum
        !* ------------------
        
        !$ACC LOOP SEQ
        DO inc=1,incdim
          zcinc=zdci(inc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jl=kidia,kfdia
            zpu(jl,jk,iazi)=zpu(jl,jk,iazi)+&
              & zact(jl,inc,iazi)*zflux(jl,inc,iazi)*zcinc
          ENDDO
        ENDDO
        
        !* end JK do-loop
        !* --------------
        
      ENDDO
      
      !* end IAZIDIM do-loop
      !* ---------------
      
    ENDDO
    
    ! -----------------------------------------------------------------------------
    
    !*       MAKE CORRECTION FOR CRITICAL-LEVEL MOMENTUM DEPOSITION
    !*       ------------------------------------------------------
    
    z0p0=0._jprb
    zrgpts=1.0_JPRB/(rg*ptstep)
    !$ACC LOOP SEQ
    DO iazi=1,iazidim
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=kidia,kfdia
        zcngl(jl)=0.0_JPRB
      ENDDO
      !$ACC LOOP SEQ
      DO jk=2,klaunch
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zulm, zdft)
        DO jl=kidia,kfdia
          zulm=zcosang(iazi)*pum1(jl,jk)+zsinang(iazi)*pvm1(jl,jk)-zul(jl,iazi)
          zdfl(jl,jk-1,iazi)=zdfl(jl,jk-1,iazi)+zcngl(jl)
          zdft=MIN(zdfl(jl,jk-1,iazi),2.0_JPRB*(papm1(jl,jk-1)-papm1(jl,jk))*&
            & (zcrt(jl,jk-1,iazi)-zulm)*zrgpts)
          zdft=MAX(zdft,z0p0)
          zcngl(jl)=(zdfl(jl,jk-1,iazi)-zdft)
          zpu(jl,jk,iazi)=zpu(jl,jk,iazi)-zcngl(jl)
        ENDDO
      ENDDO
    ENDDO
    
    
    !*       SUM CONTRIBUTION FOR TOTAL ZONAL AND MERIDIONAL FLUX
    !*       ---------------------------------------------------
    
    !$ACC LOOP SEQ
    DO iazi=1,iazidim
      !$ACC LOOP SEQ
      DO jk=klaunch,2,-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jl=kidia,kfdia
          pfluxu(jl,jk)=pfluxu(jl,jk)+zpu(jl,jk,iazi)*zaz_fct*zcosang(iazi)
          pfluxv(jl,jk)=pfluxv(jl,jk)+zpu(jl,jk,iazi)*zaz_fct*zsinang(iazi)
        ENDDO
      ENDDO
    ENDDO
    
    
    !*    UPDATE U AND V TENDENCIES
    !*    ----------------------------
    
    zcons1=1.0_JPRB/rcpd
    !$ACC LOOP SEQ
    DO jk=1,klaunch
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdelp, ze1, ze2)
      DO jl=kidia,kfdia
        zdelp= rg/(paphm1(jl,jk+1)-paphm1(jl,jk))
        ze1=(pfluxu(jl,jk+1)-pfluxu(jl,jk))*zdelp
        ze2=(pfluxv(jl,jk+1)-pfluxv(jl,jk))*zdelp
        ptenu(jl,jk)=ze1
        ptenv(jl,jk)=ze2
      ENDDO
    ENDDO
    
    !---------------------------------------------------------------------------
    !$ACC END PARALLEL
    
    !$ACC WAIT IF(lacc)
    !$ACC END DATA
    !IF (LHOOK) CALL DR_HOOK('GWDRAG_WMS',1,ZHOOK_HANDLE)
    
  END SUBROUTINE gwdrag_wms
  
  
END MODULE mo_gwd_wms

