!
! mo_art_emission_seas
!   This module calculates the emission fluxes of sea salt aerosol. This
!   module is splitted into three subroutines, one for each mode.
!   For further description of the emission flux parameterizations,
!   see the corresponding subroutines.
!   Based on COSMO-ART code by K. LUNDGREN and H. VOGEL
!   Based on Lundgren (2010) - Direct Radiative Effects of Sea Salt on the Regional Scale
!         Dissertation an der Fakultaet fuer Physik des Karlsruher Instituts fuer Technologie (KIT)
!
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

MODULE mo_art_emission_seas
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_seas_emiss_martensson
  PUBLIC :: art_seas_emiss_monahan
  PUBLIC :: art_seas_emiss_smith
  PUBLIC :: art_seas_emiss_mode1
  PUBLIC :: art_seas_emiss_mode2
  PUBLIC :: art_seas_emiss_mode3

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_seas_emiss_martensson(u10m,v10m,dz,sst,fr_land,fr_ice,fr_lake, &
  &                                  istart,iend,emiss_rate)
!<
! SUBROUTINE art_seas_emiss_martensson
! Subroutine for the smallest particle sizes. 
! Based on: Martensson et al. (2003) - Laboratory simulations and parameterization 
!                                      of the primary marine aerosol production
!           COSMO-ART code by K. LUNDGREN and H. VOGEL
! Part of Module: mo_art_emission_seas
! Author: Daniel Rieger, KIT
! Initial Release: 2013-05-16
! Modifications:
! 2014-11-10: Daniel Rieger, KIT
! - Unified ICON/COSMO-ART emission routine.
!>
  REAL(wp),INTENT(IN) :: &
    &  u10m(:),v10m(:),  & !< 10m windspeed (both components) [m s-1]
    &  dz(:),            & !< layer thickness of bottommost layer [m]
    &  sst(:),           & !< Sea surface temperature [K]
    &  fr_land(:),       & !< Fraction of land in grid cell
    &  fr_ice(:),        & !< Fraction of sea ice in grid cell
    &  fr_lake(:)          !< Fraction of lakes in grid cell
  INTEGER,INTENT(IN)  :: &
    &  istart,iend         !< Loop indizes (start and end)
  REAL(wp),INTENT(INOUT)::    &
    &  emiss_rate(istart:iend) !< Output emission rate [ug m-3 s-1]
  ! Local Variables
  REAL(wp)            :: &
    &   c41,c31,c21,     &
    &   c11,c01,c42,     &
    &   c32,c22,c12,     &
    &   c02,c43,c33,     &
    &   c23,c13,c03
    
  REAL(wp)            :: &
    &   d41,d31,d21,     &
    &   d11,d01,d42,     &
    &   d32,d22,d12,     &
    &   d02,d43,d33,     &
    &   d23,d13,d03
    
  INTEGER             :: &
    &   n,nmax, i
    
  PARAMETER (nmax=18)
  
  REAL(wp)            :: &
    &   w,fi,            &
    &   rau,             &       !< density of dry sea salt particles [kgm^-3]
    &   d1,              &       !< median diameter in [m] for the smallest particle mode according
                                 !    to measurements of O'Dowd et al. [1997]
    &   sigma1,          &
    &   mcomp,           & 
    &   afak(nmax),      &
    &   bfak(nmax),      &
    &   dp(nmax),        &
    &   summe, vb
  
  PARAMETER (c41 = -2.576e035_wp,c31 =  5.932e028_wp,c21 = -2.867e021_wp,  &
    &        c11 = -3.003e013_wp,c01 = -2.881e006_wp,c42 = -2.452e033_wp,  &
    &        c32 =  2.404e027_wp,c22 = -8.148e020_wp,c12 =  1.183e014_wp,  &
    &        c02 = -6.743e006_wp,c43 =  1.085e029_wp,c33 = -9.841e023_wp,  &
    &        c23 =  3.132e018_wp,c13 = -4.165e012_wp,c03 =  2.181e006_wp)

  PARAMETER (d41 =  7.188e037_wp,d31 = -1.616e031_wp,d21 =  6.791e023_wp,  &
    &        d11 =  1.829e016_wp,d01 =  7.609e008_wp,d42 =  7.368e035_wp,  &
    &        d32 = -7.310e029_wp,d22 =  2.528e023_wp,d12 = -3.787e016_wp,  &
    &        d02 =  2.279e009_wp,d43 = -2.859e031_wp,d33 =  2.601e026_wp,  &
    &        d23 = -8.297e020_wp,d13 =  1.105e015_wp,d03 = -5.800e008_wp)
     
  PARAMETER (rau=2.2e3_wp,d1=0.2e-6_wp,sigma1=1.9_wp)

  dp =(/2.24000E-08_wp,2.81999E-08_wp,3.55016E-08_wp,4.46939E-08_wp,         &
    &   5.62663E-08_wp,7.08350E-08_wp,8.91760E-08_wp,1.12266E-07_wp,         &
    &   1.41335E-07_wp,1.77930E-07_wp,2.24000E-07_wp,2.81999E-07_wp,         &
    &   3.55016E-07_wp,4.46939E-07_wp,5.62663E-07_wp,7.08350E-07_wp,         &
    &   8.91760E-07_wp,1.12266E-06_wp/)
  
  DO i = istart, iend
    ! ----------------------------------
    ! --- Check: Is there open water within this gridpoint?
    ! ----------------------------------
    IF ((fr_land(i)+fr_ice(i)+fr_lake(i)) < 1.0_wp) THEN
    
      ! Get horizontal 10m wind speed
      vb=SQRT(u10m(i)*u10m(i) + v10m(i)*v10m(i))
    
      ! ----------------------------------
      ! --- Calculate the flux
      ! ----------------------------------
   
      summe = 0.0_wp
      w     = (3.84e-6_wp)*(vb**(3.41_wp))
   
      DO n=1,9
        afak(n)=c41*dp(n)**4+c31*dp(n)**3+c21*dp(n)**2   &
          &     + c11*dp(n)+c01
            
        bfak(n)=d41*dp(n)**4+d31*dp(n)**3+d21*dp(n)**2   &
          &     + d11*dp(n)+d01
        fi=sst(i)*afak(n)+bfak(n)
        summe=summe+w*fi
      ENDDO
   
      DO n=10, 13
        afak(n)=c42*dp(n)**4+c32*dp(n)**3+c22*dp(n)**2   &
          &     + c12*dp(n)+c02
        bfak(n)=d42*dp(n)**4+d32*dp(n)**3+d22*dp(n)**2   &
          &     + d12*dp(n)+d02
        fi=sst(i)*afak(n)+bfak(n)
        summe=summe+w*fi
      ENDDO
   
   
      DO n=14,nmax
        afak(n)=c43*dp(n)**4+c33*dp(n)**3+c23*dp(n)**2   &
          &     + c13*dp(n)+c03
        bfak(n)=d43*dp(n)**4+d33*dp(n)**3+d23*dp(n)**2   &
          &     + d13*dp(n)+d03
        fi=sst(i)*afak(n)+bfak(n)
        summe=summe+w*fi
      ENDDO
   
      mcomp=(1._wp/6._wp)*pi*rau*d1**3*EXP((9._wp/2._wp*(LOG(sigma1))**2))
      
      ! Calculate emission rate
      emiss_rate(i) = mcomp*summe*0.1_wp/dz(i)*1.e09_wp
      ! And consider fraction of grid cell which is not open water
      emiss_rate(i) = emiss_rate(i)*(1.0_wp - (fr_land(i)+fr_ice(i)+fr_lake(i)))
    ELSE !(fr_land+fr_ice+fr_lake) .GE. 1.0_wp
      emiss_rate(i) = 0.0_wp
    ENDIF
  ENDDO

END SUBROUTINE art_seas_emiss_martensson
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_seas_emiss_monahan(u10m,v10m,dz,fr_land,fr_ice,fr_lake,istart,iend,emiss_rate)
!<
! SUBROUTINE art_seas_emiss_monahan
! Subroutine for the medium sized particles.
! Based on: Monahan et al. (1986) - A model of marine aerosol generation
!                                   via whitecaps and wave disruption
!           COSMO-ART code by K. LUNDGREN and H. VOGEL
! Part of Module: mo_art_emission_seas
! Author: Daniel Rieger, KIT
! Initial Release: 2013-05-16
! Modifications:
! 2014-11-10: Daniel Rieger, KIT
! - Unified ICON/COSMO-ART emission routine.
!>
  REAL(wp),INTENT(IN) :: &
    &  u10m(:),v10m(:),  & !< 10m windspeed (both components) [m s-1]
    &  dz(:),            & !< layer thickness of bottommost layer [m]
    &  fr_land(:),       & !< Fraction of land in grid cell
    &  fr_ice(:),        & !< Fraction of sea ice in grid cell
    &  fr_lake(:)          !< Fraction of lakes in grid cell
  INTEGER,INTENT(IN)  :: &
    &  istart,iend         !< Loop indizes (start and end)
  REAL(wp),INTENT(INOUT)   :: &
    &  emiss_rate(istart:iend) !< Output emission rate [ug m-3 s-1]
  ! Local Variables
  INTEGER             ::  &
    &   n,nmax, i
  
  PARAMETER (nmax=9)
  
  REAL(wp)            ::  &
    &   rh,r80,aa,bb,     &
    &   m,frac,           &
    &   A,B,rd,cfak,      &
    &   dFddp(nmax),      &
    &   dFdr80,           &
    &   dFdlogdp(nmax),   &
    &   summe,            &
    &   mcompM,           &
    &   rau,d2,sigma2,    &
    &   dp(nmax),         &
    &   vb
    
    
  PARAMETER (rau=2.2e3_wp,d2=2.e-6_wp,sigma2=2._wp)
  
  dp =(/1.41335_wp,1.77930_wp,2.24000_wp,2.81999_wp,3.55016_wp,  &
    &   4.46939_wp,5.62663_wp,7.08351_wp,8.91761_wp/)

  rh=0.8_wp
  m=10._wp
  
  DO i = istart, iend
    summe=0.0_wp
    ! ----------------------------------
    ! --- Check: Is there open water within this gridpoint?
    ! ----------------------------------
    IF ((fr_land(i)+fr_ice(i)+fr_lake(i)) < 1.0_wp) THEN
    
      ! Get horizontal 10m wind speed
      vb=SQRT(u10m(i)*u10m(i) + v10m(i)*v10m(i))
    
      ! ----------------------------------
      ! --- Calculate the flux
      ! ----------------------------------
      
      DO n=1,nmax
        rd=0.5_wp*dp(n)
        frac=(2.0_wp-rh)/(1.0_wp-rh)
        r80=rd*(4.0_wp/3.7_wp)*frac**(1._wp/3._wp)
        B=(0.380_wp-LOG10(r80))/0.650_wp
        aa=-B**2
        bb=EXP(aa)
        A=10_wp**(1.19_wp*bb)
     
        cfak=1.373_wp*r80**(-3)*(1._wp+0.057_wp*r80**1.05_wp)*A
     
        dFdr80=vb**(3.41_wp)*cfak
        dFddp=2.0_wp/3.7_wp*frac**(1._wp/3._wp)*dFdr80
        dFdlogdp=dp(n)*LOG(m)*dFddp
     
        summe=summe+dFdlogdp(n)
      ENDDO
     
      mcompM=(1._wp/6._wp)*pi*rau*d2**3*EXP((9._wp/2._wp*(LOG(sigma2))**2))
     
      emiss_rate(i) = mcompM*summe*0.1_wp/dz(i)*1.e09_wp
      ! And consider fraction of grid cell which is not open water
      emiss_rate(i) = emiss_rate(i)*(1.0_wp - (fr_land(i)+fr_ice(i)+fr_lake(i)))
    ELSE !(fr_land+fr_ice+fr_lake) .GE. 1.0_wp
      emiss_rate(i) = 0.0_wp
    ENDIF
  ENDDO
END SUBROUTINE art_seas_emiss_monahan
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_seas_emiss_smith(u10m,v10m,dz,fr_land,fr_ice,fr_lake,istart,iend,emiss_rate)
!<
! SUBROUTINE art_seas_emiss_smith
! Subroutine for the large particles.
! Based on: Smith et al. (1993) - Marine aerosol concentrations and
!                                 estimated fluxes over the sea
!           COSMO-ART code by K. LUNDGREN and H. VOGEL
! Part of Module: mo_art_emission_seas
! Author: Daniel Rieger, KIT
! Initial Release: 2013-05-16
! Modifications:
! 2014-11-10: Daniel Rieger, KIT
! - Unified ICON/COSMO-ART emission routine.
!>
  REAL(wp),INTENT(IN) :: &
    &  u10m(:),v10m(:),  & !< 10m windspeed (both components) [m s-1]
    &  dz(:),            & !< layer thickness of bottommost layer [m]
    &  fr_land(:),       & !< Fraction of land in grid cell
    &  fr_ice(:),        & !< Fraction of sea ice in grid cell
    &  fr_lake(:)          !< Fraction of lakes in grid cell
  INTEGER,INTENT(IN)  :: &
    &  istart,iend         !< Loop indizes (start and end)
  REAL(wp),INTENT(INOUT)   :: &
    &  emiss_rate(istart:iend)  !< Output emission rate [ug m-3 s-1]
  ! Local Variables
  INTEGER             :: &
    &   n,nmax, i
    
  PARAMETER(nmax=5)
  
  REAL(wp)            ::  &
    &   r01, r02,         &
    &   f1 , f2,          &
    &   a1,a2,m,          &
    &   dfak,efak,        &
    &   dp(nmax),         &
    &   dfdlogdp,dr80ddp, &
    &   ddpdlogdp,        &
    &   dfdr801,dfdr802,  &
    &   rau,d3,sigma3,    &
    &   mcomps,           &
    &   summe, vb
          
  PARAMETER(r01 = 2.1_wp,r02 = 9.2_wp,f1  = 3.1_wp,f2  = 3.3_wp)
  
  PARAMETER(rau=2.2e3_wp,d3=12e-6_wp,sigma3=1.7_wp)

  m=10_wp

  dp =(/11.2266_wp,14.1335_wp,17.7930_wp,22.4000_wp,    &
    &   28.2000_wp/)
  DO i = istart, iend
    ! ----------------------------------
    ! --- Check: Is there open water within this gridpoint?
    ! ----------------------------------
    IF ((fr_land(i)+fr_ice(i)+fr_lake(i)) < 1.0_wp) THEN
    
      ! Get horizontal 10m wind speed
      vb=SQRT(u10m(i)*u10m(i) + v10m(i)*v10m(i))
      
      ! ----------------------------------
      ! --- Calculate the flux
      ! ----------------------------------
    
      summe=0.0_wp
      
      ! Derivative of PARAMETERization of Lewis&Schwartz
      dr80ddp=2.0_wp/3.7_wp*(6._wp**(1._wp/3._wp))
    
    
      DO n=1,nmax
        ddpdlogdp=dp(n)*LOG(m)
    
        dfak=EXP(-f1*(LOG((dp(n)*2.0_wp/3.7_wp*6._wp**(1._wp/3._wp))/(r01)))**2)
        efak=EXP(-f2*(LOG((dp(n)*2.0_wp/3.7_wp*6._wp**(1._wp/3._wp))/(r02)))**2)
    
        a1=10._wp**(0.0676_wp*vb+2.43_wp)
        a2=10._wp**(0.959_wp*(vb**(1._wp/2._wp))-1.476_wp)
    
        dfdr801=a1*dfak
        dfdr802=a2*efak
        dfdlogdp=ddpdlogdp*dr80ddp*(dfdr801+dfdr802)
    
        summe=summe+dfdlogdp
      END DO
    
    
      mcomps = (1._wp/6._wp)*pi*rau*d3**3*EXP((9._wp/2._wp*(LOG(sigma3)**2)))
    
      emiss_rate(i) = mcomps*summe*0.1_wp/dz(i)*1.e09_wp
      ! And consider fraction of grid cell which is not open water
      emiss_rate(i) = emiss_rate(i)*(1.0_wp - (fr_land(i)+fr_ice(i)+fr_lake(i)))
    ELSE !(fr_land+fr_ice+fr_lake) .GE. 1.0_wp
      emiss_rate(i) = 0.0_wp
    ENDIF
  ENDDO
END SUBROUTINE art_seas_emiss_smith
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_seas_emiss_mode1(u10m,v10m,dz,sst,fr_land,fr_ice,fr_lake, &
  &                                  istart,iend,emiss_rate)
!<
! SUBROUTINE art_seas_emiss_mode1
! Subroutine for the smallest particle sizes. 
! Size range: (0.01 - 0.3)*10^-6 m
! Based on: Grythe et al. (2014) - A review of sea-spray aerosol source
!                                    functions using a large global set of sea
!                                    salt aerosol concentration measurements
! Part of Module: mo_art_emission_seas
! Author: Lisa Muth, KIT
! Initial Release: 2021-06-01
!>
  REAL(wp),INTENT(IN) :: &
    &  u10m(:),v10m(:),  & !< 10m windspeed (both components) [m s-1]
    &  dz(:),            & !< layer thickness of bottommost layer [m]
    &  sst(:),           & !< Sea surface temperature [K]
    &  fr_land(:),       & !< Fraction of land in grid cell
    &  fr_ice(:),        & !< Fraction of sea ice in grid cell
    &  fr_lake(:)          !< Fraction of lakes in grid cell
  INTEGER,INTENT(IN)  :: &
    &  istart,iend         !< Loop indizes (start and end)
  REAL(wp),INTENT(INOUT):: &
    &  emiss_rate(:)       !< Output emission rate [ug m-3 s-1]

  INTEGER             :: &
    &   n,nmax, i

  PARAMETER (nmax=25)

  REAL(wp)            :: &
    &   rau,             &       !< density of dry sea salt particles [kgm^-3]
    &   d1,              &       !< median diameter in [m] for the smallest
                                 !    particle mode according
                                 !    to measurements of O'Dowd et al. [1997]
    &   sigma1,          &
    &   mcomp,           &
    &   dFddp(nmax),     &
    &   dFdlogdp(nmax),  &
    &   dp(nmax),        &
    &   summe, vb,       &
    &   T_w, m

  PARAMETER (rau=2.2e3_wp,d1=0.1e-6_wp,sigma1=1.9_wp)

  m =10._wp
  ! dp [10^-6 m], dlogdp = 0.1
  dp =(/0.01122_wp,0.01413_wp,0.01778_wp,0.02240_wp,         &
    &   0.02819_wp,0.03550_wp,0.04469_wp,0.05626_wp,         &
    &   0.07083_wp,0.08917_wp,0.11226_wp,0.14133_wp,         &
    &   0.17793_wp,0.22400_wp,0.28199_wp,0.35516_wp,         &
    &   0.44694_wp,0.56266_wp,0.70835_wp,0.89176_wp,         &
    &   1.12266_wp,1.41335_wp,1.77930_wp,2.24000_wp,         &
    &   2.28199_wp/)

  DO i = istart, iend
    ! ----------------------------------
    ! --- Check: Is there open water within this gridpoint?
    ! ----------------------------------
    IF ((fr_land(i)+fr_ice(i)+fr_lake(i)) < 1.0_wp &
      & .AND. sst(i) > 273.15_wp) THEN

      ! Get horizontal 10m wind speed
      vb=SQRT(u10m(i)*u10m(i) + v10m(i)*v10m(i))

      ! Get temperature fit function
      T_w = 0.3_wp+0.1_wp*(sst(i)-273.15_wp)      &
          & -0.0076_wp*(sst(i)-273.15_wp)**2      &
          & +0.00021_wp*(sst(i)-273.15_wp)**3 

      ! ----------------------------------
      ! --- Calculate the flux
      ! ----------------------------------

      summe = 0.0_wp


      DO n=1,nmax
        dFddp=T_w*235._wp*vb**3.5_wp*EXP(-0.55_wp*(LOG(dp(n)/0.1_wp))**2)
        dFdlogdp=dp(n)*LOG(m)*dFddp

        summe=summe+dFdlogdp(n)
      ENDDO

      mcomp=(1._wp/6._wp)*pi*rau*d1**3*EXP((9._wp/2._wp*(LOG(sigma1))**2))

      ! Calculate emission rate
      emiss_rate(i) = mcomp*summe*0.1_wp/dz(i)*1.e09_wp
      ! And consider fraction of grid cell which is not open water
      emiss_rate(i) = emiss_rate(i)*(1.0_wp - (fr_land(i)+fr_ice(i)+fr_lake(i)))
    ELSE !(fr_land+fr_ice+fr_lake) .GE. 1.0_wp
      emiss_rate(i) = 0.0_wp
    ENDIF
  ENDDO
END SUBROUTINE art_seas_emiss_mode1
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_seas_emiss_mode2(u10m,v10m,dz,sst,fr_land,fr_ice,fr_lake,istart,iend,emiss_rate)
!<
! SUBROUTINE art_seas_emiss_mode3
! Subroutine for the medium sized particles.
! Size range: (0.3 - 10)*10^-6 m
! Based on: Grythe et al. (2014) - A review of sea-spray aerosol source
!                                    functions using a large global set of sea
!                                    salt aerosol concentration measurements
! Part of Module: mo_art_emission_seas
! Author: Lisa Muth, KIT
! Initial Release: 2021-06-01
!>
  REAL(wp),INTENT(IN) :: &
    &  u10m(:),v10m(:),  & !< 10m windspeed (both components) [m s-1]
    &  dz(:),            & !< layer thickness of bottommost layer [m]
    &  sst(:),           & !< Sea surface temperature in K
    &  fr_land(:),       & !< Fraction of land in grid cell
    &  fr_ice(:),        & !< Fraction of sea ice in grid cell
    &  fr_lake(:)          !< Fraction of lakes in grid cell

  INTEGER,INTENT(IN)  :: &
    &  istart,iend         !< Loop indizes (start and end)
  REAL(wp),INTENT(INOUT):: &
    &  emiss_rate(:)       !< Output emission rate [ug m-3 s-1]
  ! Local Variables
  INTEGER             ::  &
    &   n,nmax, i

  PARAMETER (nmax=5)

  REAL(wp)            ::  &
    &   dFddp(nmax),      &
    &   dFdlogdp(nmax),   &
    &   summe,            &
    &   nrfluxM(nmax),    &
    &   mcompM,           &
    &   rau,d2,sigma2,    &
    &   dp(nmax),         &
    &   vb,T_w, m


  PARAMETER (rau=2.2e3_wp,d2=3.0e-6_wp,sigma2=2.0_wp)

  ! dp [10^-6 m], dlogdp = 0.1
  dp =(/ 3.55016_wp,4.46939_wp,5.62663_wp,7.08351_wp,        &
     &   8.91761_wp/)

  m=10._wp

  DO i = istart, iend
    ! ----------------------------------
    ! --- Check: Is there open water within this gridpoint?
    ! ----------------------------------
   IF ((fr_land(i)+fr_ice(i)+fr_lake(i)) < 1.0_wp &
     & .AND. sst(i) > 273.15_wp) THEN

      ! Get horizontal 10m wind speed
      vb=SQRT(u10m(i)*u10m(i) + v10m(i)*v10m(i))

      ! Get temperature fit function
      T_w = 0.3_wp+0.1_wp*(sst(i)-273.15_wp)      &
          & -0.0076_wp*(sst(i)-273.15_wp)**2  &
          & +0.00021_wp*(sst(i)-273.15_wp)**3 

      ! ----------------------------------
      ! --- Calculate the flux
      ! ----------------------------------

      summe = 0.0_wp


      DO n=1,nmax
        dFddp=T_w*0.2_wp*vb**3.5_wp*EXP(-1.5_wp*(LOG(dp(n)/3._wp))**2)
        dFdlogdp=dp(n)*LOG(m)*dFddp

        summe=summe+dFdlogdp(n)
      ENDDO

      mcompM=(1._wp/6._wp)*pi*rau*d2**3*EXP((9._wp/2._wp*(LOG(sigma2))**2))

      emiss_rate(i) = mcompM*summe*0.1_wp/dz(i)*1.e09_wp
      ! And consider fraction of grid cell which is not open water
      emiss_rate(i) = emiss_rate(i)*(1.0_wp - (fr_land(i)+fr_ice(i)+fr_lake(i)))
    ELSE !(fr_land+fr_ice+fr_lake) .GE. 1.0_wp
      emiss_rate(i) = 0.0_wp
    ENDIF
  ENDDO
END SUBROUTINE art_seas_emiss_mode2
!!
!!------------------------------------------------------------------------
!!
  SUBROUTINE art_seas_emiss_mode3(u10m,v10m,dz,sst,fr_land,fr_ice,fr_lake,istart,iend,emiss_rate)
!<
! SUBROUTINE art_seas_emiss_mode3
! Subroutine for the large particles.
! Size range: (10 - 28)*10^-6 m
! Based on: Grythe et al. (2014) - A review of sea-spray aerosol source
!                                    functions using a large global set of sea
!                                    salt aerosol concentration measurements
! Part of Module: mo_art_emission_seas
! Author: Lisa Muth, KIT
! Initial Release: 2021-06-01
!>
  REAL(wp),INTENT(IN) :: &
    &  u10m(:),v10m(:),  & !< 10m windspeed (both components) [m s-1]
    &  dz(:),            & !< layer thickness of bottommost layer [m]
    &  sst(:),           & !< Sea surface temperature in K
    &  fr_land(:),       & !< Fraction of land in grid cell
    &  fr_ice(:),        & !< Fraction of sea ice in grid cell
    &  fr_lake(:)          !< Fraction of lakes in grid cell
  INTEGER,INTENT(IN)  :: &
    &  istart,iend         !< Loop indizes (start and end)
  REAL(wp),INTENT(INOUT):: &
    &  emiss_rate(:)       !< Output emission rate [ug m-3 s-1]
  ! Local Variables
  INTEGER             ::  &
    &   n,nmax, i

  PARAMETER (nmax=5)

  REAL(wp)            ::  &
    &   dFddp(nmax),      &
    &   dFdlogdp(nmax),   &
    &   summe,            &
    &   nrfluxM(nmax),    &
    &   mcompM,           &
    &   rau,d3,sigma3,    &
    &   dp(nmax),         &
    &   vb,T_w, m


  PARAMETER (rau=2.2e3_wp,d3=30.0e-6_wp,sigma3=1.7_wp)

  ! dp [10^-6 m], dlogdp = 0.1
  dp =(/11.2266_wp,14.1335_wp,17.7930_wp,22.4000_wp,    &
     &  28.2000_wp/)

  m=10._wp

  DO i = istart, iend
    ! ----------------------------------
    ! --- Check: Is there open water within this gridpoint?
    ! ----------------------------------
    IF ((fr_land(i)+fr_ice(i)+fr_lake(i)) < 1.0_wp &
      & .AND. sst(i) > 273.15_wp) THEN

      ! Get horizontal 10m wind speed
      vb=SQRT(u10m(i)*u10m(i) + v10m(i)*v10m(i))

      ! Get temperature fit function
      T_w = 0.3_wp+0.1_wp*(sst(i)-273.15_wp)      &
          & -0.0076_wp*(sst(i)-273.15_wp)**2      &
          & +0.00021_wp*(sst(i)-273.15_wp)**3 


      ! ----------------------------------
      ! --- Calculate the flux
      ! ----------------------------------

      summe = 0.0_wp

      DO n=1,nmax
        dFddp=T_w*6.8e-3_wp*vb**3*EXP(-1._wp*(LOG(dp(n)/30._wp))**2)
        dFdlogdp=dp(n)*LOG(m)*dFddp

        summe=summe+dFdlogdp(n)
      ENDDO

      mcompM=(1._wp/6._wp)*pi*rau*d3**3*EXP((9._wp/2._wp*(LOG(sigma3))**2))

      emiss_rate(i) = mcompM*summe*0.1_wp/dz(i)*1.e09_wp
      ! And consider fraction of grid cell which is not open water
      emiss_rate(i) = emiss_rate(i)*(1.0_wp - (fr_land(i)+fr_ice(i)+fr_lake(i)))
    ELSE !(fr_land+fr_ice+fr_lake) .GE. 1.0_wp
      emiss_rate(i) = 0.0_wp
    ENDIF
  ENDDO
END SUBROUTINE art_seas_emiss_mode3
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_seas
