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
! YOEGWWMS --  parameters/switches for the Warner McIntyre GW parameterization

MODULE data_gwd
  
  ! Modules used:
  
  USE mo_kind,               ONLY: jprb=>wp     , &
    & jpim =>i4
  
  USE mo_cuparameters,       ONLY: lhook,   dr_hook
  USE mo_nwp_tuning_config,  ONLY: tune_gfluxlaun, tune_gcstar
  !==============================================================================
  
  IMPLICIT NONE

  PRIVATE
  

  !==============================================================================
  
  !SAVE
  !-----------------------------------------------------------------
  
  LOGICAL :: lozpr     !If .TRUE. then enhancement of mom. flux over tropics
  ! INTEGER(KIND=jpim) :: nlaunch   !launch height of gravity wave spectrum (Pa)
  INTEGER(KIND=jpim) :: nslope    !slope at small-m end of spectrum
  INTEGER(KIND=jpim) :: ngauss    !if LOZPR=TRUE and GGAUSS=2 then gaussian
  ! distribution of GFLUXLAUN based on GGAUSSA and GGAUSSB
  !if LOZPR=TRUE and GGAUSS=1 then FLUXLAUN=GCOEFF*PPRECIP
  REAL(KIND=jprb)    :: gfluxlaun !total launch momentum flux in each azimuth
  REAL(KIND=jprb)    :: gcstar    !C* (see McLandress and Scinocca 2005)
  REAL (KIND=jprb)   :: gptwo     !2*p where p is the exponent of omega for the expression
  ! of the launch energy density
  REAL(KIND=jprb)    :: gtphygwwms!Time frequency (s) of call of scheme
  
  REAL(KIND=jprb)    :: ggaussa   !gaussian distribution half-width (used if GGAUSS=1)
  REAL(KIND=jprb)    :: ggaussb   !height of gaussian distribution (ie amplification factor)
  REAL(KIND=jprb)    :: gcoeff    !if GGAUSS=1 then FLUXLAUN=GCOEFF*PPRECIP
  
  !==============================================================================
  
  
  PUBLIC :: sugwwms          ! All constants and variables in this module are public
  PUBLIC :: gfluxlaun, gcstar, gptwo, nslope, &
    & lozpr, ggaussa, ggaussb, ngauss, gcoeff !, GTPHYGWWMS
  
  
 
CONTAINS
  
  
  SUBROUTINE sugwwms(ksmax,nflevg,ppref,klaunch)
    
    !
    ! INITIALZE YOEGWWMS, THE MODULE THAT CONTROLS THE WARNER MCINTYRE GW PARAMETRIZATION
    
    !          A.ORR            ECMWF     August 2008
    
    !          INTERFACE
    !          ---------
    !          CALLED FROM *SUPHEC*
    
    !          MODIFICATIONS
    !          -------------
    !          P. Bechtold, ECMWF (October 2008) Redefine computation of launch level
    
    !
    !--------------------------------------------------------------------------------
    
    !USE PARKIND1  ,ONLY : JPIM     ,JPRB
    !!USE YOMDIM    ,ONLY : NFLEVG
    !!USE YOMSTA    ,ONLY : STPRE
    !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
    
    !USE YOEGWWMS   ,ONLY : NLAUNCH, GFLUXLAUN, GCSTAR, GPTWO, NSLOPE, &
    !                     & LOZPR, GGAUSSA, GGAUSSB, NGAUSS, GCOEFF, GTPHYGWWMS
    
    !IMPLICIT NONE
    
    INTEGER(KIND=jpim),INTENT(in),OPTIONAL:: ksmax  ! horizontal (spectral) resolution of host model
    
    INTEGER(KIND=jpim),INTENT(in):: nflevg ! number of vertical levels
    REAL(KIND=jprb), INTENT(in)  :: ppref(nflevg) ! reference pressure to determine departure level
    INTEGER(KIND=jpim),INTENT(out):: klaunch      ! launch height of gravity wave spectrum (Pa)
    
    ! internal
    INTEGER(KIND=jpim) :: jk
    REAL(KIND=jprb)    :: zlaunchp        ! launch height of gw spectrum in Pa
    
    REAL(KIND=jprb) :: zhook_handle
    
    !             Set values of parameters/switches
    !              --------------------------------
    
    
    IF (lhook) CALL dr_hook('SUGWWMS',0,zhook_handle)
    
    !*  SPECIFICATION OF SOURCE SPECTRUM
    !*  --------------------------------
    
    zlaunchp=45000.0_JPRB   !launch height (Pa)
    gfluxlaun=tune_gfluxlaun !total launch momentum flux in each azimuth (rho_o x F_o)
                             ! (set in mo_nwp_tuning_nml)
    nslope=1                !s (1,0,-1 are valid values) s is the slope at small-m
    !end of the launch spectrum
    gptwo=2.0_JPRB          !2*p (3 or 2 are valid values) p is the exponent of omega
    !in the expression of the launch spectrum
    gcstar=tune_gcstar         !C* (see McLandress and Scinocca 2005)
    
    !* Extra parameters to introduce some latitudinal/seasonal variation of the launch spectrum.
    !* If LOZPR=TRUE and GGAUSS=1 then launch momentum flux is proportional to total precipitation
    !* If LOZPR=TRUE and GGAUSS=2 then launch momentum flux has a gaussian distribution
    
    lozpr=.FALSE.           !If .TRUE. then variable launch momemtum flux
    ngauss=1
    ggaussa=10.0_JPRB       !gaussian distribution half-width in degrees (used if GGAUSS=2)
    ggaussb=0.5_JPRB        !height of gaussian distribution (used if GGAUSS=2)
    gcoeff=2000.0_JPRB      !multiplicative factor (used if GGAUSS=1)
    
    
    !*  COMPUTE MODEL LEVEL LAUNCH HEIGHT OF GW SPECTRUM
    !*  ------------------------------------------------
    
    klaunch=nflevg-1
    DO jk=nflevg,2,-1
      ! IF(STPRE(JK) > ZLAUNCHP)klaunch=JK
      IF(ppref(jk) > zlaunchp)klaunch=jk
    ENDDO
    
    !*  SETUP TIME FREQUENCY CALL OF SCHEME AS FUNCTION OF MODEL RESOLUTION
    !*  -------------------------------------------------------------------
    
    IF(PRESENT(ksmax))THEN
      gtphygwwms=3600.0_JPRB
      IF (ksmax == 255) THEN
        gtphygwwms=5400.0_JPRB
      ELSEIF (ksmax < 255) THEN
        gtphygwwms=7200.0_JPRB
      ENDIF
    ENDIF
    !---------------------------------------------------------------------------
    
    IF (lhook) CALL dr_hook('SUGWWMS',1,zhook_handle)
    
  END SUBROUTINE sugwwms
  
  
END MODULE data_gwd
