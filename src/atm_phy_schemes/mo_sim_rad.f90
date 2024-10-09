! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!  Description:
!  This module computes simplified radiation for stratocumulus case - DYCOMS-II
!  see Stevens et al. 2005
!
! Routines (module procedure)
!
!     - sim_rad!

MODULE mo_sim_rad

USE mo_kind,               ONLY: ireals=>wp     , &
                                 iintegers=>i4

USE mo_physical_constants, ONLY: cpd

  PRIVATE

  INTEGER, PARAMETER :: DKIND = SELECTED_REAL_KIND(13,300)
  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  REAL(KIND=DKIND)   :: F0_RAD=70.0_DKIND !W/m2
  REAL(KIND=DKIND)   :: F1_RAD=22.0_DKIND !W/m2
  REAL(KIND=DKIND)   :: KAPPA_RAD=85.0_DKIND !m2/kg
  REAL(KIND=DKIND)   :: ALFA_Z=1.0_DKIND !m^{-4/3}
  REAL(KIND=DKIND)   :: D_RAD=3.75E-6_DKIND !s^{-1}
  REAL(KIND=DKIND)   :: QT_ZI=8E-3_DKIND !g/kg


  PUBLIC  :: sim_rad

CONTAINS

SUBROUTINE SIM_RAD(KIDIA,KFDIA,KLON,KTDIA,KLEV,&
&               MZF,MZH,QV,QL,QI,RHO,&
&               DT_RAD)


IMPLICIT NONE
INTEGER, PARAMETER :: DKIND = SELECTED_REAL_KIND(13,300)
INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)

INTEGER(KIND=IKIND)   ,INTENT(IN)    :: KIDIA
INTEGER(KIND=IKIND)   ,INTENT(IN)    :: KFDIA
INTEGER(KIND=IKIND)   ,INTENT(IN)    :: KLON
INTEGER(KIND=IKIND)   ,INTENT(IN)    :: KTDIA
INTEGER(KIND=IKIND)   ,INTENT(IN)    :: KLEV
REAL(KIND=DKIND)   ,INTENT(IN)    :: MZF(KLON,KLEV)
REAL(KIND=DKIND)   ,INTENT(IN)    :: MZH(KLON,0:KLEV)
REAL(KIND=DKIND)   ,INTENT(IN)    :: QV(KLON,KLEV)
REAL(KIND=DKIND)   ,INTENT(IN)    :: QL(KLON,KLEV)
REAL(KIND=DKIND)   ,INTENT(IN)    :: QI(KLON,KLEV)
REAL(KIND=DKIND)   ,INTENT(IN)    :: RHO(KLON,KLEV)
REAL(KIND=DKIND)   ,INTENT(OUT)   :: DT_RAD(KLON,KLEV)

REAL(KIND=DKIND)  :: ZQT(KLON,KLEV)       !total specific moisture
REAL(KIND=DKIND)  :: ZF_RAD(KLON,KLEV)    !flux
REAL(KIND=DKIND)  :: ZQ_RAD(KLON,KLEV)    !kappa_rad . liquid water path
REAL(KIND=DKIND)  :: ZRHOH(KLON,0:KLEV)   !density on half levels
REAL(KIND=DKIND)  :: ZSIG_ZI(KLON,0:KLEV) !index indicating whethe we are above ZI
REAL(KIND=DKIND)  :: ZI(KLON)             !inversion layer height
REAL(KIND=IKIND)  :: ZI_FOUND(KLON)       !index for ZI detection


REAL(KIND=DKIND)  ::  Z4THIRDS=4.0_DKIND/3.0_DKIND
REAL(KIND=DKIND)  ::  ZTHIRD=1.0_DKIND/3.0_DKIND
INTEGER(KIND=IKIND)  :: JLON, JLEV
  
!liquid water path and zi
DO JLON=KIDIA,KFDIA
  ZQ_RAD(JLON,KLEV) = 0.0_DKIND
  ZI(JLON) = 0.0_DKIND
  ZI_FOUND(JLON)=0_IKIND
ENDDO
DO JLEV=KLEV-1,KTDIA,-1
  DO JLON=KIDIA,KFDIA
    ZQT(JLON,JLEV)=QV(JLON,JLEV)+QL(JLON,JLEV)+QI(JLON,JLEV)
    ZQ_RAD(JLON,JLEV) = ZQ_RAD(JLON,JLEV+1)+KAPPA_RAD*RHO(JLON,JLEV)*&
    &QL(JLON,JLEV)/(1.0_DKIND-ZQT(JLON,JLEV))*(MZF(JLON,JLEV)-MZF(JLON,JLEV+1))
    IF (ZQT(JLON,JLEV).LE.QT_ZI.AND.(ZI_FOUND(JLON).NE.1_IKIND)) THEN
      ZI(JLON)=MZF(JLON,JLEV)
      ZI_FOUND(JLON)=1_IKIND
    ENDIF
    !if(jlon.eq.1) THEN
    !        write(*,*) "sim_radm1:",jlev,ZQ_RAD(JLON,JLEV),&
    !        &ZQT(JLON,JLEV),RHO(JLON,JLEV),QL(JLON,JLEV),ZI(JLON)
    !endif
  ENDDO
ENDDO

!rho on half levels
DO JLEV=KTDIA,KLEV-1
  DO JLON=KIDIA,KFDIA
    ZRHOH(JLON,JLEV)=0.5_DKIND*(RHO(JLON,JLEV)+RHO(JLON,JLEV+1))
    IF(MZH(JLON,JLEV)>ZI(JLON)) THEN
      ZSIG_ZI(JLON,JLEV)=1.0_DKIND
    ELSE IF (MZH(JLON,JLEV)==ZI(JLON)) THEN
      ZSIG_ZI(JLON,JLEV)=0.5_DKIND
    ELSE
      ZSIG_ZI(JLON,JLEV)=0.0_DKIND
    ENDIF
  ENDDO
ENDDO
DO JLON=KIDIA,KFDIA
    ZRHOH(JLON,KLEV)=ZRHOH(JLON,KLEV-1)
    ZSIG_ZI(JLON,KLEV)=0.0_DKIND
ENDDO

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
     ZF_RAD(JLON,JLEV)=F0_RAD*EXP(-1.0_DKIND*(ZQ_RAD(JLON,KTDIA)-ZQ_RAD(JLON,JLEV)))+&
     &F1_RAD*EXP(-1.0_DKIND*(ZQ_RAD(JLON,JLEV)))
     IF(ZSIG_ZI(JLON,JLEV)>0 )THEN
       ZF_RAD(JLON,JLEV)=ZF_RAD(JLON,JLEV)+ZSIG_ZI(JLON,JLEV)*&
       &ZRHOH(JLON,JLEV)*cpd*D_RAD*ALFA_Z*((MZH(JLON,JLEV)-ZI(JLON))**Z4THIRDS*0.25_DKIND+&
       &ZI(JLON)*(MZH(JLON,JLEV)-ZI(JLON))**ZTHIRD)
     ENDIF
     
    !if(jlon.eq.1) THEN
    !        write(*,*) "sim_rad0:",jlev,ZF_RAD(JLON,JLEV),&
    !        F0_RAD,F1_RAD,ZQ_RAD(JLON,JLEV),ZSIG_ZI(JLON,JLEV),&
    !        ZRHOH(JLON,JLEV),cpd*D_RAD*ALFA_Z,(MZH(JLON,JLEV)-ZI(JLON))**Z4THIRDS,&
    !        ZI(JLON)*(MZH(JLON,JLEV)-ZI(JLON))**ZTHIRD,MZH(JLON,JLEV)-ZI(JLON)
    !endif
  ENDDO
ENDDO

!update tendencies with contribution from radiation
DO JLEV=KTDIA+1,KLEV
  DO JLON=KIDIA,KFDIA
    DT_RAD(JLON,JLEV)=-1.0_DKIND*(ZF_RAD(JLON,JLEV-1)-ZF_RAD(JLON,JLEV))/&
    &(MZH(JLON,JLEV-1)-MZH(JLON,JLEV))/(cpd*RHO(JLON,JLEV))
    !if(jlon.eq.1) THEN
    !        write(*,*) "sim_rad:",jlev,DT_RAD(JLON,JLEV),ZF_RAD(JLON,JLEV-1),&
    !        ZF_RAD(JLON,JLEV-1),MZH(JLON,JLEV-1)-MZH(JLON,JLEV),RHO(JLON,JLEV),&
    !        cpd
    !ENDIF
  ENDDO
ENDDO

DO JLON=KIDIA,KFDIA
    DT_RAD(JLON,KTDIA)=DT_RAD(JLON,KTDIA+1)
ENDDO

END SUBROUTINE SIM_RAD


END MODULE mo_sim_rad
