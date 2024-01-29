!
! mo_art_nucleation_cold
! This module calculates the number cocentration of activated aerosol particles for ice clouds,
! according to the model Barahona and Nenes (2008, 2009).
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

MODULE mo_art_nucleation_cold
  USE mo_kind,                          ONLY: wp

  IMPLICIT NONE

  PRIVATE  
  PUBLIC :: IceParam

!==================================================================

! VARIABLES USED IN THE ICE NUCLEATION PARAMETERIZATION 
!
!=================================================================

REAL(wp) :: amw_ice, wmw_ice, rgas_ice, dhs_ice, cpa_ice,                    &
  &         rv_ice, denice_ice, vpresw_ice, vpresi_ice,                      &
  &         denair_ice, depcoef_ice, diff_ice, aircond_ice, thaccom_ice,     &
  &         ddry_ice, np_ice, nin_ice, sh_ice,                               &
  &         grav_ice, pi_ice ,alfa_ice, beta_ice, shom_ice, koft_ice,        &
  &         dliq_ice, normv_ice, g1_ice, g2_ice, gdoin_ice, z_ice, doin_ice, &
  &         vmin_ice, vmax_ice, miuv_ice, sigmav_ice,                        &
  &         norg_ice, sigorg_ice, dorg_ice, dbc_ice, sigbc_ice, S_grid 

! kbc_ice and kdust_ice are set, but never given a value!!!
REAL(wp):: lambda_ice, kdust_ice, kbc_ice, shdust_ice,                       &
  &        shbc_ice, effdust_ice, effbc_ice, del1dust_ice, si0dust_ice,      &
  &        del1bc_ice, si0bc_ice, nbc_ice

REAL(wp), ALLOCATABLE :: ndust_ice(:), sigdust_ice(:), ddust_ice(:) 

INTEGER :: typeofspec_ice, nbindust_ice 
LOGICAL :: purehet_ice, purehom_ice
REAL(wp)    :: dH1smooth

DATA wmw_ice /.018_wp/  ! Water molecular weight
DATA amw_ice  /0.029_wp/  ! Air molecular Weight 
DATA rgas_ice /8.314_wp/  ! Universal gas constant
DATA grav_ice /9.81_wp/  ! gravity
DATA cpa_ice  /1005.1_wp/   !Air thermal gas capacity
DATA pi_ice /3.1415927_wp/   
DATA depcoef_ice /0.1_wp/    ! Default deposition coefficient
DATA thaccom_ice /0.7_wp/ !Default thermal accommodation coefficient

TYPE t_ullrich_param
  REAL(wp) :: a
  REAL(wp) :: b
  REAL(wp) :: c
  REAL(wp) :: d
  REAL(wp) :: e
END TYPE t_ullrich_param

TYPE(t_ullrich_param) :: param_soot
TYPE(t_ullrich_param) :: param_dust

!
!=======END of decalarations================================================================

CONTAINS

!*************************************************************************
!
! ICE FREEZING PARAMETERIZATION FILES START HERE
!
! ************************************************************************
!
!
!======================================================================
!
!       Code Developer
!       Donifan Barahona, GA TECH
!       donifan@gatech.edu
!    -------------------------------
!     DESCRIPTION
!
!***********************************************************
!** Parameterization  of ICE crystal number concentration 
!** for large scale models. 
!** Donifan Barahona, Athanasios Nenes
!   JGR, 111, D11211,  2008
!   ACP, 9, 369-381,   2009 
!   ACP  9, 5933-5948, 2009
!   Homogeneoeus and heterogeneous nucleation considered   
!*** SI units unless otherwise specified. 
! 
!
! *** WRITTEN BY DONIFAN BARAHONA
!
!=======================================================================

SUBROUTINE IceParam (lscalew,sigma_w, T, P,s_gridin, np_hom,ddry, ndust, ddust, &
     &                nbc, dia_bc, norg, dia_org,                               &
     &                lhet,lhom,sigdust,sigbc,sigorg,                           &
     &                nhet, nice, smax, nlim, iaci_cold) 

  REAL(wp), INTENT(in) :: T, P, s_gridin, np_hom, nbc, sigma_w, norg, dia_org,  &
    &                      dia_bc,lscalew,sigbc,sigorg,ddry

  REAL(wp)             :: ndust(:), ddust(:),sigdust(:)

  LOGICAL, INTENT(in)  :: lhet,lhom

  REAL(wp), INTENT(out) :: nice, nhet, smax, nlim  

  INTEGER, INTENT(in) :: iaci_cold
  
  REAL(wp) :: wpar_ice, T_ice, P_ice

  LOGICAL :: use_av_v, use_grids

  !Aerosol properties 
  ! Number concentration

  np_ice       = np_hom 
  nbindust_ice = SIZE(ndust) 

  ALLOCATE (sigdust_ice(nbindust_ice))
  ALLOCATE (ndust_ice(nbindust_ice))
  ALLOCATE (ddust_ice(nbindust_ice))

  ndust_ice = ndust
  nbc_ice   = nbc
  norg_ice  = norg

  ! Mean size  
  ddry_ice  = ddry      !Homogeneous aerosol mean geometric diameter, m, by now it is constant but
                        ! can be passed
  dbc_ice   = dia_bc
  dorg_ice  = dia_org
  ddust_ice = ddust

! Std deviation 

  sigdust_ice = DLOG(sigdust)  ! calculate log directly here     
  sigbc_ice   = DLOG(sigbc)
  sigorg_ice  = DLOG(sigorg)

! Other conditions
  P_ice       = P
  T_ice       = T
  S_grid      = s_gridin
  depcoef_ice = 0.1_wp  !Deposition coefficient 

!Heteroegeneous effects =============================== 

  typeofspec_ice = iaci_cold
  !-1 - Monodisperse
  ! 1- Meyers 
  ! 2- BKG, Phillips 2007 
  ! 3- Barahona (2011, in prep) Formulation based shifted water activity method
  ! 4- PDA08, using fixed size distributions.
  ! 5- Ullrich (2014, in prep) ice parameterization

  purehet_ice= lhet  !logical True supresses homog nucleation
  purehom_ice= lhom   ! True supresses het nucleation

!Updraft Velocity Distribution =================================
  use_av_v= .FALSE. !set to false to integrate over a pdf of updraft
  use_grids = .TRUE. !set to true to use avoid pdf calculations

!parameters used in the CNT spectrum 
  shdust_ice  = 0.2_wp  !maximum freezing threeshold dust only used in CNT
  effdust_ice = 0.05_wp !maximum freezing efficiency dust  
  shbc_ice    = 0.3_wp  !maximum freezing threeshold bc only used in CNT
  effbc_ice   = 0.05_wp !maximum freezing efficiency bc

!parameters used in the Monodisperse approximation 
  sh_ice  = 0.3_wp  !assume freeizng threeshold
  nin_ice = (ndust_ice(1)+nbc_ice)*0.05_wp

  CALL prop_ice(T_ice, P_ice)

  IF (use_av_v) THEN

    wpar_ice = MAX(0.01_wp, sigma_w*0.8_wp) !m/s !characteristic velocity. Minimum 1 cm/s
    CALL nice_param(wpar_ice, T_ice, nice, smax, nhet, nlim) 

  ELSEIF(use_grids .AND. (T_ice >= 235._wp .OR. purehet_ice) ) THEN

    CALL nice_param_sgrid(T_ice, nice, smax, nhet, nlim)

  ELSE

    vmin_ice=0.00001_wp
    vmax_ice=5._wp

    sigmav_ice = sigma_w  !standard deviation of the V distribution m/
    miuv_ice   = lscalew  !mean velocity for the V distribution 
                          !(can be set to the large scale updraft)

    vmin_ice   = MAX(miuv_ice-(4._wp*sigmav_ice), vmin_ice)
    vmax_ice   = MIN(miuv_ice+(4._wp*sigmav_ice), vmax_ice)

    CALL nice_Vdist(T_ice, nice, smax, nhet, nlim)

  END IF          

  DEALLOCATE (sigdust_ice)
  DEALLOCATE (ndust_ice)
  DEALLOCATE (ddust_ice)      

  RETURN

END SUBROUTINE IceParam


!*************************************************************
!    SUBROUTINE nice_Vdist. Calculates the ice crystal number concentration
!    at the maximum supersaturation using a PDF of updraft using a 
!    sixth order Gauss-LegENDre quadrature  
!     Inputs:  T, and P all SI units)
!    Output NC (m-3)
!    Barahona and Nenes, JGR, D11211 (2008) and ACPD, 15665-15698, (2008) 
!   Written by Donifan Barahona
!************************************************************ 


SUBROUTINE nice_Vdist(T_ice, nice, smax, nhet, nlim) 

  REAL(wp) :: quadx(6), wpar, sum1, quadw(6), dp, &
    &        sum2, sum3, sum4, x1, x2 
  REAL(wp), INTENT(in)  :: T_ice
  REAL(wp), INTENT(out) :: nice, smax, nhet, nlim 
  INTEGER :: index

  DATA quadx/0.23861918_wp, -0.23861918_wp, 0.66120939_wp, &
    &   -0.66120939_wp, 0.93246951_wp, -0.93246951_wp/

  DATA quadw/0.46791393_wp, 0.46791393_wp, 0.36076157_wp, &
    &   0.36076157_wp, 0.17132449_wp,  0.17132449_wp/

!calculate the integral in the denominator

  x1=(vmin_ice-miuv_ice)/(SQRT(2._wp)*sigmav_ice)
  x2=(vmax_ice-miuv_ice)/(SQRT(2._wp)*sigmav_ice)

  normv_ice=(ERF(x2)-ERF(x1))*0.5_wp!  !cummulative width of the distribution of velocities 

  IF(normv_ice==0.0_wp) THEN
    PRINT*, 'NORMVICE:',normv_ice,miuv_ice,sigmav_ice,vmin_ice,vmax_ice, &
            'N=',ERF(x2),ERF(x1),x2,x1
    !^ Transformation into "CALL message(..,..)"
    miuv_ice=0.5_wp
    sigmav_ice=0.05_wp 
    x1=(vmin_ice-miuv_ice)/(SQRT(2._wp)*sigmav_ice)
    x2=(vmax_ice-miuv_ice)/(SQRT(2._wp)*sigmav_ice)
    normv_ice=(ERF(x2)-ERF(x1))*0.5_wp!  !cummulative width of the distribution of velocities
  ENDIF

  sum1 = 0._wp
  sum2 = 0._wp
  sum3 = 0._wp
  sum4 = 0._wp

!use a Gauss-LegENDre Quadrature

  DO index =1, 6
    wpar=0.5_wp*(((vmax_ice-vmin_ice)*quadx(index)) &
          +(vmax_ice+vmin_ice))
    CALL nice_param(wpar, T_ice, nice, smax, nhet, nlim)

    CALL gausspdf(wpar, dp)
    sum1=sum1+(nice*dp*quadw(index))
    sum2=sum2+(smax*dp*quadw(index))
    sum3=sum3+(nhet*dp*quadw(index))
    sum4=sum4+(nlim*dp*quadw(index))

  END DO

  nice=sum1*(vmax_ice-vmin_ice)*0.5_wp
  smax=sum2*(vmax_ice-vmin_ice)*0.5_wp
  nhet=sum3*(vmax_ice-vmin_ice)*0.5_wp
  nlim=sum4*(vmax_ice-vmin_ice)*0.5_wp

  RETURN

END SUBROUTINE nice_Vdist


!*************************************************************     
         !Approximation to the error FUNCTION
!*************************************************************


REAL(wp) FUNCTION erfapp(x)

  REAL(wp), INTENT(in) :: x 
  REAL(wp) :: ax  
  ax=x*x*(1.27324_wp+(0.147_wp*x*x))/(1._wp+(0.147_wp*x*x))
  erfapp=SQRT(1._wp-EXP(-ax))

END FUNCTION erfapp


!*************************************************************
!    SUBROUTINE nice_param. Calculates the ice crystal number concentration
!    at the maximum supersaturation . Inputs: Wpar, T, and P all SI units)
!    Output NC (m-3)
!    Barahona and Nenes, JGR, D11211 (2008) and ACPD, 15665-15698, (2008) 
!    Written by Donifan Barahona
!************************************************************ 


SUBROUTINE nice_param(wpar_ice,T_ice,          &
    &                  nice, smax, nhet, nlim_)     

  REAL(wp), INTENT(in)  :: T_ice, wpar_ice
  REAL(wp), INTENT(out) :: nice, smax, nhet, nlim_

  REAL(wp) :: findsmax
  REAL(wp) :: tao
  REAL(wp) :: s_temp

  REAL(wp) :: aux1, aux2, g,                                         &
    &        dpmax, miu, monvol, fds, nlim, dlim,  dstar,  nstar,    &
    &        nhom, fc, phido, auxnc, sizecorr, dsh, nhet_, f1, f2,   &
    &        saux, sup, slow, shalf, fhalf, dpmin, gam

  INTEGER :: index 

  monvol=np_ice*1.0_wp-6._wp*ddry_ice*ddry_ice*ddry_ice
  aux1=1.6397e-14_wp*T_ice-3.1769e-12_wp
  dpmax=aux1*(monvol**(-0.373_wp))*(wpar_ice**(-0.05_wp))

  IF (dpmax > 1.0e-4_wp) THEN
    dpmax=1.0e-4_wp
  END IF

  dpmin=dliq_ice+(0.02_wp/SQRT(alfa_ice*wpar_ice*g1_ice))  !MInimum size for DPmax (added 09/11/09)
  dpmax=MAX(dpmax,dpmin)

  aux1=dpmax-dliq_ice
  aux2=DLOG((g2_ice+(g1_ice*dpmax))/(g2_ice+(g1_ice*dliq_ice)))     
  g=1.3346_wp*((g1_ice*aux1)-(g2_ice*aux2))/(aux1*g1_ice*g1_ice)
  lambda_ice=lambda_ice/SQRT(wpar_ice)
  nstar=((g1_ice*alfa_ice*wpar_ice)**1.5_wp)/beta_ice/z_ice/SQRT(2._wp) 

  gam=g2_ice/g1_ice

!   WRITE(*, *) NSTAR, g1_ice, alfa_ice, beta_ice, z_ice

!*********IS HOMOGENEOUS FREEZING HAPPENING?********************

  fds  = 1._wp !CORRECTION TO Nc FROM HET FREEZING
  nhom = 0._wp

  IF (T_ice >= 235._wp) THEN  !no homogeneous above 235 K
    nhom=0._wp
    GOTO 685
  END IF

  IF (purehet_ice) THEN
    NHOM=0._wp
    GOTO 685
  END IF

!calculate limiting NIN for combined hom_het  

  IF (typeofspec_ice >= 0._wp) THEN!polydisperse expressions
    CALL inspec_ice(shom_ice, T_ice, nhet_, dsh)
    sizecorr=EXP(-2._wp/lambda_ice/shom_ice)
    dstar=((4._wp*dsh*dsh/3._wp)+(2._wp*dsh*(shom_ice-dsh))) &
      &    /(shom_ice-dsh+1._wp)
    nlim=nstar*(shom_ice+1._wp)/shom_ice/SQRT(dstar)/sizecorr

  ELSE !monodisperse approximation
    dsh=shom_ice-sh_ice
    dstar=((4._wp*dsh*dsh/3._wp)+(2._wp*dsh*(shom_ice-dsh))) &
      &     /(shom_ice-dsh+1._wp)
    dlim=-gam+SQRT((gam*gam)+(2._wp*dstar/g1_ice/alfa_ice/wpar_ice))
    nlim=alfa_ice*wpar_ice*(shom_ice+1._wp)/z_ice/beta_ice/shom_ice
    nlim=nlim*((g1_ice*dlim)+g2_ice)/dlim/dlim
    nhet_=nin_ice
  END IF

  nlim_=MIN(nlim, 1e+10_wp)
  nlim_=MAX(nlim, 1e-6_wp)  !avoid overflow errors

  fds=1._wp-((nhet_/nlim)**1.5_wp)

  IF (fds/=fds) THEN
    fds=1._wp
  END IF

  IF (purehom_ice) THEN
    nhet_=0._wp
    fds=1.0_wp
  END IF

  IF (FDS <= 0._wp) THEN  !Homogeneous nucleation completely inhibited by IN
    nhom=0._wp
    GOTO 685
  END IF

!********FRACTION OF FROZEN DROPLETS********************

  miu=fds*alfa_ice*(shom_ice+1._wp)*wpar_ice*koft_ice/shom_ice

  phido=((pi_ice*g/miu/2._wp)**0.5_wp)*(g/miu)
  auxnc=2._wp*denair_ice/beta_ice/koft_ice/denice_ice/pi_ice/np_ice
  fc=auxnc/phido

!calculating hom Nc

  IF (np_ice > 0._wp) THEN

    IF (fc <= 0.6_wp) THEN
      nhom=np_ice*EXP(-fc)*(1.0_wp-EXP(-fc)) 
    ELSE
      nhom=np_ice/(1._wp+EXP((9._wp-2._wp*fc)/7._wp))  !correction needed for convective clouds 
                                                !(very high updraft) (Barahona et al. JGR, 2010)
    END IF

  ELSE  
    nhom=0._wp
  END IF

  smax=shom_ice
  nhet=nhet_

  GOTO 686  !finish	

!********PURE HETEROEGENEOUS FREEZING********************

!find interval for bisection

685     smax=0_wp
        nhet=0_wp
        saux=0.01_wp 

  IF (typeofspec_ice < 0._wp) THEN
    saux=sh_ice+0.00000000001_wp !minimun smax in monodisperse case
  END IF

  s_temp=saux
  IF (typeofspec_ice >= 0._wp) THEN!polydisperse expressions

    CALL inspec_ice(s_temp, T_ice, nhet_, dsh)
    sizecorr=EXP(-2._wp/lambda_ice/s_temp)
    dstar=((4._wp*dsh*dsh/3._wp)+(2._wp*dsh*(s_temp-dsh)))/(s_temp-dsh+1._wp)
    dstar=dstar+(gdoin_ice*alfa_ice*wpar_ice)
    tao=nhet_*sizecorr*s_temp*SQRT(dstar)/(s_temp+1._wp)/nstar

  ELSE !monodisperse approximation

    dsh=s_temp-sh_ice
    dstar=((4._wp*dsh*dsh/3._wp)+(2._wp*dsh*(s_temp-dsh))) &
      &     /(s_temp-dsh+1._wp)
    dlim=-gam+SQRT((gam*gam)+(2._wp*dstar/g1_ice/alfa_ice/wpar_ice))
    tao=alfa_ice*wpar_ice*(s_temp+1._wp)/z_ice/beta_ice/s_temp
    tao=tao*((g1_ice*dlim)+g2_ice)/dlim/dlim/nin_ice

  END IF

  findsmax=1._wp-tao

  f1=findsmax

  DO index =1, 20

    saux=saux+0.1_wp
    s_temp=saux
    IF (typeofspec_ice >= 0._wp) THEN !polydisperse expressions

      CALL inspec_ice(s_temp, t_ice, nhet_, dsh)
      sizecorr=EXP(-2._wp/lambda_ice/s_temp)
      dstar=((4._wp*dsh*dsh/3._wp)+(2._wp*dsh*(s_temp-dsh)))/(s_temp-dsh+1._wp)
      dstar=dstar+(gdoin_ice*alfa_ice*wpar_ice)
      tao=nhet_*sizecorr*s_temp*SQRT(dstar)/(s_temp+1._wp)/nstar

    ELSE !monodisperse approximation
      dsh=s_temp-sh_ice
      dstar=((4._wp*dsh*dsh/3._wp)+(2._wp*dsh*(s_temp-dsh))) &
        &     /(s_temp-dsh+1._wp)
      dlim=-gam+SQRT((gam*gam)+(2._wp*dstar/g1_ice/alfa_ice/wpar_ice))
      tao=alfa_ice*wpar_ice*(s_temp+1._wp)/z_ice/beta_ice/s_temp
      tao=tao*((g1_ice*dlim)+g2_ice)/dlim/dlim/nin_ice

    END IF

    findsmax=1_wp-tao

    f2=findsmax

    IF (f2*f1 < 0._wp) GOTO 677
    f2=f1
  END DO 

  IF (f2*f1 > 0._wp) THEN
    nhet=0._wp
    smax=saux    !No NIN present in pure heterogeneous mode smax>200%
    GOTO 686
  END IF

!Perform bisection

677    sup=saux
       slow=saux-0.1_wp

  DO index=1,100
    shalf=0.5_wp*(sup+slow) 
    s_temp=shalf
    IF (typeofspec_ice >= 0._wp) THEN!polydisperse expressions
      CALL inspec_ice(s_temp, t_ice, nhet_, dsh)
      sizecorr=EXP(-2._wp/lambda_ice/s_temp)
      dstar=((4._wp*dsh*dsh/3._wp)+(2._wp*dsh*(s_temp-dsh)))/(s_temp-dsh+1._wp)
      dstar=dstar+(gdoin_ice*alfa_ice*wpar_ice)
      tao=nhet_*sizecorr*s_temp*SQRT(dstar)/(s_temp+1._wp)/nstar

    ELSE !monodisperse approximation
      dsh=s_temp-sh_ice
      dstar=((4._wp*dsh*dsh/3._wp)+(2._wp*dsh*(s_temp-dsh))) &
        &     /(s_temp-dsh+1._wp)
      dlim=-gam+SQRT((gam*gam)+(2._wp*dstar/g1_ice/alfa_ice/wpar_ice))
      tao=alfa_ice*wpar_ice*(s_temp+1._wp)/z_ice/beta_ice/s_temp
      tao=tao*((g1_ice*dlim)+g2_ice)/dlim/dlim/nin_ice

    END IF

    findsmax=1._wp-tao

    fhalf=findsmax

    IF (SIGN(1._wp,f1)*SIGN(1._wp,fhalf) <= 0._wp) THEN  
      f2    = fhalf
      sup   = shalf
    ELSE
      f1    = fhalf
      slow   = shalf
    ENDIF

    IF (ABS(slow-sup) <= 1.e-3_wp) GOTO 678
  
  END DO

678    smax=shalf

  IF (typeofspec_ice >= 0._wp) THEN
    CALL inspec_ice(smax, T_ice, nhet, dsh)
  ELSE
    nhet=nin_ice !monodisperse approximation
  END IF

! Adding Heterogeneous 

686   IF (nhet/=nhet) THEN !avoiding errors
        nhet=0._wp
      END IF

  IF (nhom/=nhom) THEN !avoiding errors
    nhom=0._wp
  END IF

  nice=nhom+nhet

END SUBROUTINE nice_param


!*************************************************************
!    SUBROUTINE nice_param. Calculates the ice crystal number concentration
!    at the maximum supersaturation . Inputs: Wpar, T, and P all SI units)
!    Output NC (m-3)
!    Barahona and Nenes, JGR, D11211 (2008) and ACPD, 15665-15698, (2008)
!    Written by Donifan Barahona
!************************************************************


SUBROUTINE nice_param_sgrid(T_ice,                    &
    &                        nice, smax, nhet, nlim_)

  REAL(wp), INTENT(in)  :: T_ice
  REAL(wp), INTENT(out) :: nice, smax, nhet, nlim_

  REAL(wp) ::   dsh  ! not needed but kept to use only one INSPEC routine

! USE COSMO GRID-SCALE ICE SATURATION (PDF updraft not necessary)
  IF (typeofspec_ice >= 0._wp) THEN
    CALL inspec_ice(S_grid, T_ice, nhet, dsh)
  ELSE
    nhet=nin_ice !monodisperse approximation
  END IF

  IF (nhet/=nhet) THEN !avoiding errors
    nhet=0._wp
  END IF

! SET OUTPUTS 
  nice  = nhet     ! only heterogeneous
  nlim_ = 0.0_wp
  smax  = S_grid

  RETURN

END SUBROUTINE nice_param_sgrid

!*************************************************************
!    FUNCTION VPRESWATER. Calculates the saturated vapor pressure
!    of water (Pa) according to Murphy & Koop (2005)
!    T in K (173.15-373.15)
!************************************************************


REAL(wp) FUNCTION vpreswater_ice(t)

  real(wp), INTENT(in) :: t
  real(wp)             :: a(0:9) 

  DATA a/54.842763_wp, -6763.22_wp,  -4.21_wp, 0.000367_wp,            &
    &       0.0415_wp, 218.8_wp, 53.878_wp, -1331.22_wp, -9.44523_wp,  &
    &       0.014025_wp/

  vpreswater_ice = a(0)+(a(1)/t)+(a(2)*LOG(t))+(a(3)*t)    &
    &            + (TANH(a(4)*(t-a(5)))*((a(6)+(a(7)/t))   &
    &            + (a(8)*LOG(t))+ (a(9)*t))) 

  vpreswater_ice=EXP(vpreswater_ice)

  RETURN

END FUNCTION vpreswater_ice


!*************************************************************
!    FUNCTION VPRESICE. Calculates the saturated vapor pressure
!    of ice (pa) according to Murphy & Koop (2005)
!    T in K (>110)
!************************************************************         

REAL(wp) FUNCTION vpresice(t)

  REAL(wp), INTENT(in) :: t
  REAL(wp)             :: a(0:3)

  DATA a/9.550426_wp, -5723.265_wp, 3.53068_wp, -0.00728332_wp/

  vpresice = a(0)+(a(1)/t)+(a(2)*LOG(t))+(a(3)*t)
  vpresice = EXP(vpresice)

  RETURN

END FUNCTION vpresice


!*************************************************************
!    FUNCTION DHSUB. Calculates the latent heat of sublimation
!    of ice (J/Kg) according to Murphy & Koop (2005)
!    T in K (>30)
!*************************************************************


REAL(wp) FUNCTION dhsub_ice(t)

  REAL(wp), INTENT(in) :: t
  REAL(wp)             :: a(0:4)

  DATA a/46782.5_wp, 35.8925_wp, -0.07414_wp, 541.5_wp, 123.75_wp/

  dhsub_ice = a(0) + (a(1) * t) + (a(2)*t*t) + (a(3) &
    &       * EXP(-((t/ a(4))**2)))

  dhsub_ice=1000_wp*dhsub_ice/18_wp

  RETURN

END FUNCTION dhsub_ice


!*************************************************************
!    FUNCTION ICEDENSITY. Calculates the DENSITY OF ICE
!    of ice (Kg/m3) according to PK97 
!    T in K (>30)
!************************************************************ 

        
REAL(wp) FUNCTION densityice(t)

  REAL(wp), INTENT(in) :: t
  REAL(wp)             :: a(0:2), ttemp

  DATA a/0.9167_wp, -1.75e-4_wp, -5.0e-7_wp/

  ttemp=t-273_wp

  densityice= 1000_wp*(a(0)+(a(1)*ttemp)+(a(2)*ttemp*ttemp))

  RETURN

END FUNCTION densityice


!*************************************************************
!    FUNCTION WATDENSITY. Calculates the DENSITY OF ICE
!    of liquid water (Kg/m3) according to PK97 
!    T in K (>240)
!************************************************************        


REAL(wp) FUNCTION watdensity_ice(t)

  REAL(wp), INTENT(in) :: t
  REAL(wp)             :: a(0:6),  ttemp, watdensity
  INTEGER              :: i

  DATA a/0.99986_wp, 6.690e-5_wp, -8.486e-6_wp, 1.518e-7_wp, & 
    &      -6.9984e-9_wp, -3.6449e-10_wp, -7.497e-12_wp /

  ttemp=t-273._wp

  IF (ttemp <= -40._wp) THEN
    ttemp=-40._wp
  END IF

  watdensity=a(6)*ttemp 

    IF (t >= 240._wp) THEN 
      DO i=5,1, -1
        watdensity= (watdensity+a(i))*(ttemp)
      ENDDO
      watdensity=watdensity + a(0)
    ELSE
      watdensity=0.979_wp
    END IF 

  watdensity=watdensity*1000._wp
  watdensity_ice=watdensity

  RETURN

END FUNCTION watdensity_ice


!*************************************************************
!    SUBROUTINE PROPERTIES. Set physical an thermodynamic 
!    properties at T and P 
!************************************************************   

SUBROUTINE prop_ice(T_ice_in, P_ice)

  REAL(wp), INTENT(in)  :: T_ice_in, P_ice
  REAL(wp)              :: aux1, aux2, sw, Tc,T_ice, hdust, hbc, &
    &                      b0, b1, b2, b3, x, T0bc, T0dust, gam, gamma

  ! Sanity check

  IF (T_ice_in >= 273._wp) THEN ! Use this properties only for ice
    T_ice=273._wp   ! MB cahnged to 273
  ELSE
    T_ice=T_ice_in
  END IF 

  rv_ice     = rgas_ice/wmw_ice
  dhs_ice    = dhsub_ice(T_ice)
  vpresw_ice = vpreswater_ice(T_ice)
  vpresi_ice = vpresice(T_ice)
  denice_ice = densityice(T_ice)
  denair_ice = P_ice*amw_ice/rgas_ice/T_ice
! Kinetic properties of the bulk vapor (SI UNITS, Seinfel and Pandis, 1997)

  diff_ice=(0.211_wp*101325._wp/P_ice)*((T_ice/273._wp)**1.94_wp)*1.0e-4_wp !m^2/s
  aux1=1.0e-3_wp*(4.39_wp+0.071_wp*T_ice) !W/m

!correcting Kair for size assuming D=1e-6_wp m

  aux2=(2._wp*aux1/(thaccom_ice*1.0e-6_wp*denair_ice*cpa_ice)) &  
        *((58.0e-3_wp*pi_ice/(rgas_ice*T_ice))**0.5_wp)

  aircond_ice=aux1/(1._wp+aux2)

!Physical constants

  aux1     = grav_ice*dhs_ice/rv_ice/T_ice/T_ice/cpa_ice
  aux2     = grav_ice*amw_ice/rgas_ice/T_ice
  alfa_ice = aux1-aux2
  beta_ice = amw_ice*P_ice/wmw_ice/vpresi_ice
  gamma    = 1.5_wp*dhs_ice*dhs_ice/rv_ice/T_ice/T_ice/cpa_ice  !Correction for T>250 K

  beta_ice = beta_ice+gamma  !only for high T (>250 K, Barahona et al. JGR 2010)

!Homogeneous freezing only 
  shom_ice=2.349_wp-(T_ice/259._wp) !hom threeshold Si according to Ren & McKenzie, 2005
  sw=shom_ice*vpresi_ice/vpresw_ice
  shom_ice=shom_ice-1._wp
  koft_ice=(0.0240_wp*T_ice*T_ice)-(8.035_wp*T_ice)+934.0_wp ! constant related to Jmax, Barahona & Nenes JGR 2008

!Calculate Dliq using an approximation derived from the equilbrium calculations and the
!Approximation proposed by Lewis (2008), 13, D03205, JGR 

  IF (sw < 0.99_wp) THEN  !only subsaturated regime (Haze Aerosols)
    aux1=(1._wp/(1._wp-sw))-1.1764_wp
  ELSE
    aux1=(1._wp/0.01_wp)-1.1764_wp
  END IF 
  dliq_ice=ddry_ice*0.9344_wp*(aux1**(1._wp/3._wp))

! calculate average G for homogeneous freezing   
  aux1   = denice_ice*rv_ice*T_ice/vpresi_ice/diff_ice
  aux2   = dhs_ice*denice_ice/aircond_ice/T_ice
  aux2   = aux2*((dhs_ice/rv_ice/T_ice)-1.0_wp)
  g1_ice = (aux1+aux2)/4.0_wp

  g2_ice = denice_ice*rv_ice*T_ice/2.0_wp/vpresi_ice/depcoef_ice
  g2_ice = g2_ice*((2.0_wp*pi_ice/rv_ice/T_ice)**0.5_wp)      

  doin_ice  = 1.e-6_wp !assumed IN diameter
  gdoin_ice = (g1_ice*0.5_wp*doin_ice*doin_ice)+(g2_ice*doin_ice)
  z_ice     = denice_ice*pi_ice/2.0_wp/denair_ice

  gam=g2_ice/g1_ice
  lambda_ice=1._wp/SQRT(alfa_ice*g1_ice*gam*gam)  !divided by sqrt(wparcel) in niceparam


!!============Parameters needed for IN spectra=========

!For Phillips, et. al. 2008 spectra PDA08!!!!!!!!!!!!!!!!!!!!!

  Tc          = T_ice-273._wp
  hdust       = 0.15_wp
  T0dust      = -40._wp
  b0          = -1.0261_wp
  b1          = 3.1656e-3_wp
  b2          = 5.3938e-4_wp
  b3          = 8.2584e-6_wp
  x           = b0+(b1*Tc)+(b2*Tc*Tc)+(b3*Tc*Tc*Tc)
  si0dust_ice = 1._wp+(10._wp**x)

  del1dust_ice=cubicint_ice(Tc, T0dust, T0dust+5._wp, 1._wp, hdust) !bug corrected

  hbc        = 0._wp
  T0bc       = -50._wp
  b0         = 0.5652_wp
  b1         = 1.085e-2_wp
  b2         = -3.118e-5_wp
  si0bc_ice  = b0+(b1*T_ice)+(b2*T_ice*T_ice)-0.1_wp  !bug corrected C to K
  del1bc_ice = cubicint_ice(Tc, T0bc, T0bc+5._wp, 1._wp, hbc)

  RETURN

END SUBROUTINE prop_ice


!*************************************************************
!   SUBROUTINE gauspdf (normalized over the width of the distribution).  
!************************************************************  

SUBROUTINE gausspdf(x, dp)

  REAL(wp), INTENT(in)  :: x
  REAL(wp), INTENT(out) :: dp
  IF(sigmav_ice==0.0_wp .OR. SQRT(2._wp*pi_ice)==0._wp .OR. normv_ice==0._wp) THEN
    PRINT*, 'divzero:',miuv_ice,x, dp,sigmav_ice,pi_ice,normv_ice
    !^ Transformation into "CALL message(..,..)"
  ENDIF

  dp=EXP(-0.5_wp*(x-miuv_ice)*(x-miuv_ice)/sigmav_ice/sigmav_ice) & 
    &   /sigmav_ice/sqrt(2._wp*pi_ice)/normv_ice 

  RETURN

END SUBROUTINE gausspdf


!*************************************************************
!   FUNCTION cubicint_ice (cubic interpolation between y1 and y2 within a and b).  
!************************************************************  

REAL(wp) FUNCTION cubicint_ice(y, y1, y2, a, b)

  REAL(wp), INTENT(in) :: y, y1, y2, a, b   
  REAL(wp)             :: a_, b_, a0, a1, a2, a3, d, aux

  IF (y <= y1) THEN
    d=a
    goto 5065
  END IF 

  IF (y >= y2) THEN
    d=b
    goto 5065
  END IF 

  aux = y2-y1      
  a_  = 6._wp*(a-b)/(aux*aux*aux)
  b_  = a+(a_*(y1*y1*y1)/6._wp)-(a_*(y1*y1)*y2*0.5_wp)

  a0 = b_
  a1 = a_*y1*y2
  a2 = -a_*(y1+y2)*0.5_wp
  a3 = a_/3._wp
  d  = a0+(a1*y)+(a2*y*y)+(a3*y*y*y)

5065  cubicint_ice = d

END FUNCTION cubicint_ice


!*************************************************************
!   FUNCTION dcubicint_ice (used in the PDA08 spectrum).  
!************************************************************  


REAL(wp) FUNCTION dcubicint_ice(y, y1, y2, a, b)

  REAL(wp), INTENT(in) :: y, y1, y2, a, b   
  REAL(wp)             :: a_, a1, a2, a3, d, aux

  IF (y <= y1) THEN
    d=0
    goto 5065
    END IF 

  IF (y >= y2) THEN
    d=0
    goto 5065
  END IF 

  aux = y2-y1      
  a_  = 6._wp*(a-b)/(aux*aux*aux)

  a1 = a_*y1*y2
  a2 = -a_*(y1+y2)*0.5_wp
  a3 = a_/3._wp
  d  = (a1)+(2._wp*a2*y)+(3._wp*a3*y*y)

5065  dcubicint_ice = d

END FUNCTION dcubicint_ice


!*************************************************************
! FUNCTION PDG07 (simplIFied ice nucleation 
!                     spectra according to Phillips et. al. 2007).  
! si is supersaturation wrt ice and T is in K 
!************************************************************  

REAL(wp) FUNCTION pdg07_ice(si, t)     

  REAL(wp), INTENT(in) :: si, t
  REAL(wp)             :: n 

  IF (t <= 243._wp)THEN
    n=1000._wp*EXP(-0.388_wp)*(EXP(3.88_wp*si)-1._wp)
  ELSE
    n=60._wp*EXP(-0.639_wp)*(EXP(12.96_wp*si)-1._wp)
  END IF

  pdg07_ice=n

END FUNCTION pdg07_ice     


!*************************************************************
! SUBROUTINE INSPEC_ice
!  Provides the Ice Nuclei concentration (m-3) 
! and the chracteristic freezing threeshold, DSh (Barahona & Nenes 2009), at given 
! si and T. The variable typeofspec_ice (integer) has the values
! 1 Meyers et. al. 1992
! 2  Phillips et. al. 2007
! 3  Barahona 2011
! 4  Phillips et. al. 2008 (simplIFed) 
! 5  Ullrich et al. 2014 (in prep.)
! si is supersaturation wrt ice and T is in K 

!      Written by Donifan Barahona 
!      donifanb@umbc.edu

!************************************************************  


SUBROUTINE inspec_ice(si, T, n, dsh)

  REAL(wp), INTENT(in)  :: si, T
  REAL(wp), INTENT(out) :: n, dsh
  REAL(wp)              :: nd, nbc, aux, si_, sw, del0, ddel0,         &
    &                      fc, delw0, ddelw0, sw0, hdust, hbc,         &
    &                      nbase, dnd, dnbc, dnbase, dh,               &
    &                      dfc, ndaux, dndaux, dnorg, norg, ndustaux,  &
    &                      frac, aux2, dx2, xi_t, D_grid_bio, n_grid_bio, SIW

  REAL(wp)              ::  surface_soot               !< Only necessary for typeofspec_ice = 5
  REAL(wp), ALLOCATABLE :: surface_dust(:)             !< Only necessary for typeofspec_ice = 5
  REAL, dimension  (nbindust_ice) ::  ndust_s, ddust_s
  REAL :: n_iw, DSh_s ,  nbc_s, dbc_s, Asolo
  INTEGER               :: index

  SELECT CASE  (typeofspec_ice)

    CASE(1) !Meyers1992 
      n=1000._wp*EXP(-0.639_wp)*(EXP(12.96_wp*si)-1._wp)
      dsh=1._wp/12.96_wp

    CASE(2) !Phillips2007
      n=pdg07_ice(si, T)
      IF (T <= 243._wp)THEN
        dsh=1._wp/3.88_wp
      ELSE
        dsh=1._wp/12.96_wp
      END IF

    CASE(3) ! Barahona2010

!dust contribution
      ndustaux=0.0_wp
      DO index=1, nbindust_ice
        ndustaux=ndustaux+ndust_ice(index)
      END DO

      IF (si <= shdust_ice) THEN 
        nd=(si/shdust_ice)*ndustaux*effdust_ice &
          & * EXP(-kdust_ice*(shdust_ice-si))
        dnd=nd*((1._wp/si)+kdust_ice)

      ELSE
        nd=ndustaux*effdust_ice
        dnd=0._wp
      END IF

!soot contribution
      IF (si <= shbc_ice) THEN   
        nbc=(si/shbc_ice)*nbc_ice*effbc_ice &
          &  * EXP(-kbc_ice*(shbc_ice-si))
        dnbc=nbc*((1._wp/si)+kbc_ice)
      ELSE
        nbc=nbc_ice*effbc_ice
        dnbc=0._wp
      END IF

      n=nd+nbc
      IF ((dnd+dnbc) > 0._wp) THEN
        dsh=n/(dnd+dnbc)
      ELSE
        dsh=si
      END IF 

    CASE(4) !PDA2008. Allows multiple lognormal modes for dust. Single mode lognormal distributions
            !         are assumed for bc and organics

      si_=si+1._wp
      sw=si_*vpresi_ice/vpresw_ice

      sw0=0.97_wp
      delw0=cubicint_ice(sw, sw0, 1._wp, 0._wp, 1._wp)
      ddelw0=dcubicint_ice(sw, sw0, 1._wp, 0._wp, 1._wp)

      nbase=pdg07_ice(si, T)

      IF (T <= 243._wp) THEN
        dnbase=3.88_wp*nbase
      ELSE
        dnbase=12.96_wp*nbase
      END IF

      ! for high temperature:
      xi_t = cubicint_ice(T, -2.0_wp, -5.0_wp,  0._wp, 1._wp)

      !dust contribution
      del0=cubicint_ice(si_, si0dust_ice, si0dust_ice+0.1_wp, 0._wp, 1._wp)
      ddel0=dcubicint_ice(si_, si0dust_ice, si0dust_ice+0.1_wp, 0._wp, 1._wp)

      fc=0.5_wp*del1dust_ice*del0
      dfc=0.5_wp*del1dust_ice*ddel0

      hdust=fc+((1._wp-fc)*delw0) 
      dh=(dfc*(1._wp-delw0))+(ddelw0*(1._wp-fc))

      IF (hdust > 1._wp) THEN 
        hdust=1._wp
        dh=0._wp
      END IF

      aux=(2._wp/3._wp)*hdust*xi_t*(nbase/0.76_wp)*pi_ice/5.e-7_wp/4._wp
      aux2=(2._wp/3._wp)*pi_ice/0.76_wp/5.e-7_wp/4._wp !The 4d0 was introduced as recommnedation of
                                                      !V Phillips

      nd=0._wp
      dnd=0._wp

      DO index =1, nbindust_ice
        dx2= ddust_ice(index)*ddust_ice(index)

        frac=0.5_wp*(1._wp+erfapp(-LOG(ddust_ice(index)/0.1e-6_wp) & !fraction above 0.1 microns
          &    /sigdust_ice(index)/SQRT(2._wp)))                     !sigma_dust=log(sigma_bc_g)

        ndaux=frac*ndust_ice(index)*(1_wp-EXP(-aux*dx2))
      
        nd=nd+ndaux
        ndaux=(frac*ndust_ice(index)-ndaux)
        dndaux=ndaux*((dh*nbase)+(hdust*xi_t*dnbase))*aux2*dx2

        dnd=dnd+dndaux

      END DO

!soot contribution

      del0=cubicint_ice(si_, si0bc_ice, si0bc_ice+0.1_wp, 0._wp, 1._wp)
      ddel0=dcubicint_ice(si_, si0bc_ice, si0bc_ice+0.1_wp, 0._wp, 1._wp)

      fc=0.5_wp*del1bc_ice*del0
      hbc=fc+((1._wp-fc)*delw0)
      dfc=0.5_wp*del1bc_ice*ddel0
      dh=(dfc*(1._wp-delw0))+(ddelw0*(1_wp-fc))

      IF (hbc > 1._wp) THEN 
        hbc=1._wp
        dh=0._wp
      END IF

      frac=0.5_wp*(1._wp +erfapp(-LOG(dbc_ice/0.1e-6_wp) & 
        &   /sigbc_ice/SQRT(2._wp)))                       !sigbc=log(sigma_bc_g)
      dx2=dbc_ice*dbc_ice       
      aux=((1._wp/3._wp)-0.06_wp)*hbc*xi_t*(nbase/0.76_wp)*pi_ice/2.7e-7_wp 
      aux2=((1._wp/3._wp)-0.06_wp)*pi_ice/0.76_wp/2.7e-7_wp

      nbc=nbc_ice*frac*(1._wp-EXP(-aux*dx2))
      dnbc=(nbc_ice*frac-nbc)*((dh*nbase)+(hbc*xi_t*dnbase))*aux2*dx2

      !Organics contribution

      frac=0.5_wp*(1._wp+erfapp(-LOG(dorg_ice/0.1e-6_wp) & 
        &   /sigorg_ice/SQRT(2._wp)))         !sigorg=log(sigma_org_g)

      dx2=dorg_ice*dorg_ice
      aux=0.06_wp*hbc*xi_t*(nbase/0.76_wp)*pi_ice/9.1e-7_wp 
      aux2=0.06_wp*pi_ice/0.76_wp/9.1e-7_wp

      norg=norg_ice*frac*(1._wp-EXP(-aux*dx2))
      dnorg=(norg_ice*frac-norg)*((dh*nbase)+(hbc*xi_t*dnbase))*aux2*dx2

      n=nd+nbc+norg

      IF ((dnd+dnbc+dnorg) > 0._wp) THEN
        dsh=n/(dnd+dnbc+dnorg)
      ELSE
        dsh=si
      END IF 

    CASE(5)

      ALLOCATE (surface_dust(nbindust_ice))
    
      ! total surface area for the given number-size-distributions in m^2 for each mode
      DO index =1, nbindust_ice
        surface_dust(index) = pi_ice * ddust_ice(index)**2 * EXP( 2._wp * sigdust_ice(index)**2 )
      ENDDO
      surface_soot = pi_ice * dbc_ice**2 * EXP( 2._wp * sigbc_ice**2 )

      sw=(si+1._wp)*vpresi_ice/vpresw_ice               !< saturation ratio wrt water

      IF (sw >= 1._wp) THEN           !< freezing above water saturation
        nd = 0._wp
        DO index =1, nbindust_ice
          nd = nd + ndust_ice(index)*( 1._wp - EXP( EXP( 151.548_wp - 0.521_wp*T ) &
            &                           * surface_dust(index) ) )
        ENDDO
        !nbc = exp( 127.936 - 0.444*t ) * surface_soot * nbc_ice
        !norg = ...
        n = nd! + nbc + norg                          !< total IN number concentration
      ELSE                                            !< nucleation below water saturation
        param_soot%a = 58._wp
        param_soot%b = 0.01_wp
        param_soot%c = 200._wp
        param_soot%d = 0.015_wp
        param_soot%e = 240.0_wp

        param_dust%a = 44.0_wp
        param_dust%b = 0.01_wp
        param_dust%c = 248.0_wp
        param_dust%d = 0.2_wp
        param_dust%e = 238.0_wp

        nbc = nbc_ice*( 1._wp - EXP( het_icenuc_ullrich(T,si,param_soot) * surface_soot ) )
        nd  = 0.0_wp
        DO index =1, nbindust_ice
          nd  =nd + ndust_ice(index)*( 1._wp - EXP( het_icenuc_ullrich(T,si,param_dust) * surface_dust(index) ) )
        ENDDO
        !norg  = het_icenuc_ullrich(t,si,param_org) * surface_org

        n = nbc + nd! + norg                    !< total IN number concentration
      END IF

      dsh=si
      DEALLOCATE(surface_dust)
      
    case (6) !Phillips et al 2013.
      
      D_grid_bio =0.e0_wp !1.0e-8_wp
      n_grid_bio = 1.e-12_wp
      
      Si_=si+1._wp
      SW=Si_*vpresi_ice/vpresw_ice
      SIW=vpresw_ice/vpresi_ice
      !dust
      DO index =1, nbindust_ice

        IF(ddust_ice(index).GT.0._wp) THEN        
          frac=0.5_wp*(1._wp-erfapp(log(0.1e-6_wp/ddust_ice(index)) & !fraction above 0.1 microns
              /sigdust_ice(index)/sqrt(2._wp)))
        ELSE 
          frac=0._wp  !TESTESTEST
        ENDIF
        
        ndust_s(index) = SNGL(frac*ndust_ice(index))
        ddust_s=SNGL(ddust_ice(index))
        
      end do
      !black carbon
      
      frac=0.5_wp*(1._wp-erfapp(log(0.1e-6_wp/dbc_ice) & !fraction above 0.1 microns
          /sigbc_ice/sqrt(2._wp)))
      nbc_s = SNGL(frac*nbc_ice)
      
      dbc_ice=min(dbc_ice, 2.e-6_wp)
      dbc_s = SNGL(dbc_ice*1._wp)
         
      !Soluble organics (spherical)
      frac=0.25_wp*0.5_wp*(1._wp-erfapp(log(0.1e-6_wp/dorg_ice) & 
          /sigorg_ice/sqrt(2._wp)))     !sigorg=log(sigma_org_g) !20-30% of organics are soluble (Saxena and HIdeelman 1996). 
      
      Asolo = SNGL(frac*norg_ice*3.1415*dorg_ice*dorg_ice) !
      
      call      EMPIRICAL_PARAM_PHILLIPS(SNGL(T), SNGL(Si_), SNGL(SIW), SNGL(SW), &
                                         ddust_s, ndust_s, 3, & !5
                                         (/dbc_s/), (/nbc_s/), 1, &
                                         (/SNGL(D_grid_bio)/), (/SNGL(n_grid_bio)/), 1, Asolo, &
                                          n_iw, DSh_s)
      N=DBLE(n_iw)
      DSh=DBLE(DSh_s)
      
    CASE DEFAULT

    n=0._wp
    dsh=si

  END SELECT

  IF (dsh >= si) THEN
    dsh=si
  END IF 

END SUBROUTINE inspec_ice
!
!-------------------------------------------------------------------------------------
!
FUNCTION het_icenuc_ullrich(temp,si,param)RESULT(inas_density)
  REAL(wp),INTENT(in)                :: temp         !< temperature in K
  REAL(wp),INTENT(in)                :: si           !< supersaturation over ice
  TYPE(t_ullrich_param),INTENT(in)   :: param        !< storage for constant parameters 
                                                     !    (fct(aerosol class))
  REAL(wp)                           :: inas_density !< ice nucleating active site density in m^-2
  REAL(wp)                           :: cotang       !< Arc cotangent
  
  cotang = pi_ice / 2._wp - ATAN(param%d * (temp - param%e))
    !< conversion formula see Bronstein,Semendjajew,Musiol,Muehlig:
    !< Taschenbuch der Mathematik, 6. Auflage, Harri Deutsch Verlag Ch.2.8, Eqn. 2.149
  
  inas_density = EXP(param%a*(si)**(0.25_wp) * COS(param%b * ((temp - param%c)))**2 * cotang &
    &                / pi_ice)

END FUNCTION het_icenuc_ullrich

!=======================================================================================
!=======================================================================================
!=======================================================================================
!!====================================================================================
! EMPIRICAL PARAMETERISATION (Phillips et al. 2013, JAS)
! contributed by Vaughan Phillips, 2012
! University of Leeds
! Implementation:   Donifan Barahona donifan.o.barahona@nasa.gov 
!====================================================================================


SUBROUTINE EMPIRICAL_PARAM_PHILLIPS(temperature_K, SI, SIW, SW,          &
                                    D_grid_dust, n_grid_dust, ijstop_dust,      &
                                    D_grid_soot, n_grid_soot, ijstop_soot,      &
                                    D_grid_bio, n_grid_bio, ijstop_bio, A_solo, &
                                    n_iw, DSH)
  implicit none
  real, intent(IN):: temperature_K, SI, SIW, SW, A_solo
  real, dimension(:), intent(IN):: D_grid_dust, n_grid_dust, &
                      D_grid_soot, n_grid_soot, D_grid_bio, n_grid_bio
  integer, intent(IN):: ijstop_dust, ijstop_soot, ijstop_bio
  real :: nin_a_nuc_dust, nin_a_nuc_soot, nin_a_nuc_bio, nin_a_nuc_solo, &
   num_ic_dust_imm, num_ic_soot_imm, num_ic_bio_imm, num_ic_solo_imm 
   
  real, intent (inout)  :: DSH, n_iw!DONIF
  
  real ::  dn_in_dust, dn_in_soot, dn_in_bio, dNall, dNaux, naux, SS_w, &
           dH_frac_dust, dH_frac_soot,  aux, dfdep !DONIF
   
  
  REAL ::  RHO_CFDC, &
           BASE_DUST_OMEGA, BASE_SOOT_PHILIC_OMEGA, BASE_BIO_OMEGA, &
           ALPHA_DUST, &
           ALPHA_SOOT, ALPHA_bio, FRACTION_DEPNUCL_WARM_DUST, PIE, BASE_SOLO_OMEGA, &
           TEMP_MAX_DUST_DEGC, TEMP_MAX_SOOT_DEGC, TEMP_MAX_bio_DEGC, GLASS_FRAC
  PARAMETER( BASE_DUST_OMEGA = 2.0e-6_wp, &
             BASE_SOOT_PHILIC_OMEGA = 1.e-7_wp, &  
             BASE_BIO_OMEGA = 0.89e-6_wp, BASE_SOLO_OMEGA = 5.6e-5_wp, GLASS_FRAC = 0.5,&  
             ALPHA_DUST = 2./3., ALPHA_SOOT = 1./3. - 0.03, ALPHA_bio = 0.03, &
             RHO_CFDC = 50000./(287.*228.15), FRACTION_DEPNUCL_WARM_DUST = 0.15, PIE = 3.1415926, &
             TEMP_MAX_DUST_DEGC = -10., TEMP_MAX_SOOT_DEGC = -15., TEMP_MAX_bio_DEGC = -2.)
  
  real :: FAC_CORRECT_RH = 2.,  rho_AIDA
  real::   H_frac_dust, n_in, n_in_dust, n_in_ultra, n_in_dust_ultra,  &
           CIHENC_dust, SS_i, &
           n_in_soot_ultra, &
           H_frac_soot, H_frac_bio, n_in_soot, n_in_bio, n_in_bio_ultra, &
           CIHENC_soot, CIHENC_bio,   &
           n_in_max, SS_iw
  
  real ::  H_frac_solO, RHI, n_in_solO, n_in_solO_star, CIHENC_solO, &
           Psi_solO
  
  real ::  mu, S_i_0, tc_HM_degC
  real ::  S_w_0, dep_frac, n_in_hat, n_in_tilde
  integer :: ij
  !intrinsic :: exp, DEXP, SIZE, DBLE
  
  !print *, SIZE(n_grid_dust(:))
  if(ijstop_dust .ne. SIZE(n_grid_dust)) stop 6366
  if(ijstop_soot .ne. SIZE(n_grid_soot)) stop 6366
  if(ijstop_bio .ne. SIZE(n_grid_bio)) stop 6366
  
  
  dNaux=12.96 !default
  naux =0.0
  dNall =dNaux
  DSh=0.0
  n_iw=0.0
  nin_a_nuc_dust=0.0; nin_a_nuc_soot=0.0; nin_a_nuc_bio=0.0; nin_a_nuc_solo=0.0
  num_ic_dust_imm=0.0; num_ic_soot_imm=0.0; num_ic_bio_imm=0.0; num_ic_solo_imm=0.0 
  dn_in_dust=0.0; dn_in_soot=0.0; dn_in_bio=0.0
  !n_in_dust=0.0;n_in_soot=0.0;n_in_bio=0.0
  
  dH_frac_dust = 0.0
  dH_frac_soot = 0.0
  dH1smooth=0.0
  aux=0.0
  
  
  
  !====================================================================================
  ! COMPUTATION BLOCK 
  !
  !====================================================================================
  !
  rho_AIDA = 90000./(287.*205.)	
  
  Psi_solO = A_solO/BASE_SOLO_OMEGA
  SS_i = SI-1.0 !everything is based on supersaturation 
  SS_w = SW-1.0
  SS_iw = SIW - 1.0
  
  
  
  if(SS_i >  0.0) then
    if(temperature_K < 273.15 .and. temperature_K > 273.15 - 90. ) then
      !SS_iw = QSW/QSI - 1.
      
      if(SS_w   > 0.) then
        SS_i = SS_iw
        SS_w = 0.0
      end if 
!      S_i_zero = 1.15 !this is taken care of
      
      
      tc_HM_degC = temperature_K - 273.15
      
      
      S_i_0 = 1. + 10.**(8.2584e-6_wp*tc_HM_degC*tc_HM_degC*tc_HM_degC + 5.3938E-4_wp*tc_HM_degC*tc_HM_degC &
            + 3.1656E-3_wp*tc_HM_degC - 1.0261)
      
      
      S_w_0 = 0.97
      
      S_w_0 = 0.97
      
      aux =H_1_smooth(-(temperature_K-273.15), 35., 40., FRACTION_DEPNUCL_WARM_DUST, 1.)/FAC_CORRECT_RH
      
      dep_frac = H_1_smooth(SS_i + 1, S_i_0,  S_i_0 + 0.1, 0.,1.)* aux
      dfdep=dH1smooth*aux
      
      aux= H_1_smooth(SS_w + 1.0, S_w_0, 1., 0.,1.)
      
      H_frac_dust = dep_frac  + (1. - dep_frac)*aux
      
      dH_frac_dust = dfdep + (SIW*(1. - dep_frac)*dH1smooth)- aux*dfdep
      
      if(H_frac_dust > 1.) H_frac_dust = 1.
      
      if ((H_frac_dust .gt. 1.0e-6_wp) .and. (H_frac_dust .lt. 1.)) then 
        dH_frac_dust = dH_frac_dust/H_frac_dust
      else
        dH_frac_dust =0.0
      end if 
      
      !soluble organics
      S_i_0 = 1.2
      
      aux =H_1_smooth(-(temperature_K-273.15), 65., 75., 0.,1.)
      dep_frac = H_1_smooth(SS_i + 1, S_i_0, S_i_0+0.1, 0.,1.)*aux   
      H_frac_solO = dep_frac  
      
      if(H_frac_solO > 1.) H_frac_solO = 1.
      
      S_w_0 = 0.97
      
      S_i_0 = 1.3
          
      aux = H_1_smooth(-(temperature_K-273.15), 40., 50., 0.,1.) &
      	 /FAC_CORRECT_RH
      dep_frac = H_1_smooth(SS_i + 1, S_i_0, S_i_0+0.1, 0.,1.)* aux 
      
      dfdep= dH1smooth*aux
      
      aux = H_1_smooth(SS_w + 1.0, S_w_0, 1., 0.,1.)
      H_frac_soot = dep_frac  + (1. - dep_frac)*aux
      if(H_frac_soot > 1.) H_frac_soot = 1.
      
      dH_frac_soot = dfdep + (SIW*(1. - dep_frac)*dH1smooth)- aux*dfdep
      if ((H_frac_soot .gt. 1.0e-6_wp) .and. (H_frac_soot .lt. 1.)) then 
         dH_frac_soot = dH_frac_soot/H_frac_soot
      else
         dH_frac_soot =0.0
      end if 
      
      
      
      H_frac_bio = H_frac_soot 
      
      if(temperature_K < 273.15 .and. temperature_K >= 273.15 - 35.) then
        n_in = 1.E3* (exp(12.96*SS_i - 0.639)/RHO_CFDC) *0.0587*FAC_CORRECT_RH
        if( temperature_K > 273.15 -5. .and. temperature_K < 273.15 - 2. ) then
          n_in = n_in*H_1_smooth(-(temperature_K-273.15), 2., 5., 0., 1.)
        endif
        if(temperature_K >= 273.15 - 2. ) n_in = 0.
        
        
        if(temperature_K < 273.15 -25. ) then
          n_in_tilde = 1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC  
          n_in_hat = n_in
          
          if(temperature_K >= 273.15 - 30.) n_in_max = 1.E3* (exp(12.96*SS_iw - 0.639)/RHO_CFDC) *0.0587*FAC_CORRECT_RH
          if(temperature_K < 273.15 - 30.) n_in_max = 1000.*(exp(0.1296*(SS_iw*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC 
          
          if(n_in_hat > n_in_max) n_in_hat = n_in_max
          if(n_in_tilde > n_in_max) n_in_tilde = n_in_max
          
          
          n_in = n_in_hat * ((n_in_tilde/n_in_hat)**(H_1_smooth(-(temperature_K-273.15), 25., 35., 0., 1.)))
          
          if(n_in > n_in_max) n_in = n_in_max
          
        endif
        n_in_dust = 0.
        dn_in_dust = 0.
        
        
        if(temperature_K < 273.15 - 30.) then   !DONIF
          dnaux = 3.88 !this is a simplified derivative of dNds
        else
          dnaux = 12.96
        end if 
        
        naux=0.0
        
        do ij = 1, ijstop_dust
          mu = n_in*ALPHA_DUST*H_frac_dust*PIE*D_grid_dust(ij)*D_grid_dust(ij)&
             /BASE_DUST_OMEGA
          naux = (1. - exp(-mu))*n_grid_dust(ij)
          n_in_dust = n_in_dust + naux 
          dn_in_dust = max(mu*(n_grid_dust(ij)-naux)*(dnaux + dH_frac_dust), 0.0) + dn_in_dust
        enddo
        
        if( temperature_K > 273.15 +TEMP_MAX_DUST_DEGC - 20. .and. temperature_K < 273.15 + TEMP_MAX_DUST_DEGC) then
          n_in_dust = n_in_dust*H_1_smooth(-(temperature_K-273.15),-TEMP_MAX_DUST_DEGC,-TEMP_MAX_DUST_DEGC+20., 0., 1.)
        endif
        if(temperature_K >= 273.15 + TEMP_MAX_DUST_DEGC) n_in_dust = 0.
        
        
        n_in_soot = 0.
        dn_in_soot = 0.
        do ij = 1, ijstop_soot
          mu = n_in*ALPHA_SOOT*H_frac_soot*PIE*(D_grid_soot(ij)**2.) &
          /BASE_SOOT_PHILIC_OMEGA
          naux = (1. - exp(-mu))*n_grid_soot(ij)
          n_in_soot = n_in_soot + naux 
          dn_in_soot = max(mu*(n_grid_soot(ij)-naux)*(dnaux+dH_frac_soot), 0.0) + dn_in_soot
        enddo
        
        if( temperature_K > 273.15 + TEMP_MAX_SOOT_DEGC - 10. .and. temperature_K < 273.15 + TEMP_MAX_SOOT_DEGC) then
          n_in_soot = n_in_soot*H_1_smooth(-(temperature_K-273.15),-TEMP_MAX_SOOT_DEGC,-TEMP_MAX_SOOT_DEGC+10., 0., 1.)
        endif
        if(temperature_K >= 273.15 + TEMP_MAX_SOOT_DEGC) n_in_soot = 0.
        
        n_in_bio = 0.
        dn_in_bio = 0.
        do ij = 1, ijstop_bio
          mu = n_in*ALPHA_bio*H_frac_bio*PIE*(D_grid_bio(ij)**2.) &
             /BASE_BIO_OMEGA
          
          mu = n_in*ALPHA_bio*H_frac_bio
          naux  =  (1. - exp(-mu))*n_grid_bio(ij)
          !!!!!!!!!!!!!!remember!!!!!!!!!!!!!!!!
          !naux =  n_in*ALPHA_bio*H_frac_bio
          !!!!!!!!!!!!!!remember!!!!!!!!!!!!!!!!
          n_in_bio = n_in_bio + naux 
          dn_in_bio = max(mu*(n_grid_bio(ij)-naux)*dnaux, 0.0) + dn_in_bio
        enddo
        
        
        if( temperature_K > 273.15 + TEMP_MAX_bio_DEGC - 3. .and. temperature_K < 273.15 + TEMP_MAX_bio_DEGC) then
          n_in_bio = n_in_bio*H_1_smooth(-(temperature_K-273.15),-TEMP_MAX_bio_DEGC,-TEMP_MAX_bio_DEGC+3., 0., 1.)
        endif
        
        if(temperature_K >= 273.15 + TEMP_MAX_bio_DEGC ) n_in_bio = 0.
        
        
        
      else
        n_in = 0.; n_in_ultra = 0.; n_in_dust = 0.;  n_in_soot  = 0.; n_in_bio = 0.;
      endif
      
      if(temperature_K < 273.15 - 35.) then
        n_in_ultra = 1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC
        dnaux = 3.88 !DONIF simplified treatment of derivative asusming dH small 
        naux=0.0
        
        
        RHI = (SS_i+1.)*100.
        if(RHI < 0.) RHI = 0.
        n_in_solO_star = 1000.e6_wp*(7.7211e-5_wp * RHI - 9.2688e-3_wp)/rho_AIDA
        
        n_in_dust_ultra = 0.; 
        dn_in_dust = 0.0
        do ij = 1, ijstop_dust
          mu = n_in_ultra*ALPHA_DUST*H_frac_dust*PIE*(D_grid_dust(ij)**2.) &
             /BASE_DUST_OMEGA
          naux = (1. - exp(-mu))*n_grid_dust(ij)
          n_in_dust_ultra = n_in_dust_ultra + naux 
          dn_in_dust = max(mu*(n_grid_dust(ij)-naux)*(dnaux +dH_frac_dust), 0.0) + dn_in_dust
        
        enddo
        
        
        n_in_soot_ultra = 0.0
        dn_in_soot = 0.0
        do ij = 1, ijstop_soot
          mu = n_in_ultra*ALPHA_SOOT*H_frac_soot*PIE*(D_grid_soot(ij)**2.) &
             /BASE_SOOT_PHILIC_OMEGA
          naux =   (1. - exp(-mu))*n_grid_soot(ij)						
          n_in_soot_ultra = n_in_soot_ultra + naux 
          dn_in_soot = max(mu*(n_grid_soot(ij)-naux)*(dnaux +dH_frac_soot), 0.0) + dn_in_soot
        enddo
        
        
        n_in_bio_ultra = 0.
        dn_in_bio = 0.0
        do ij = 1, ijstop_bio
          mu = n_in_ultra*ALPHA_bio*H_frac_bio*PIE*(D_grid_bio(ij)**2.) &
             /BASE_BIO_OMEGA
          naux  =  (1. - exp(-mu))*n_grid_bio(ij)
          n_in_bio_ultra = n_in_bio_ultra + naux 
          dn_in_bio = max(mu*(n_grid_bio(ij)-naux)*dnaux, 0.0) + dn_in_bio
        enddo
        
        
        n_in_solO = Psi_solO*glass_frac*H_frac_solO*n_in_solO_star
        
      else
        n_in_ultra = 0.; n_in_dust_ultra = 0.; n_in_soot_ultra = 0.; n_in_solO = 0.; n_in_bio_ultra = 0.; 
      endif
      
      
      n_in_dust = n_in_dust + n_in_dust_ultra;
      n_in_soot = n_in_soot + n_in_soot_ultra;
      n_in_bio = n_in_bio + n_in_bio_ultra;
      
      
  
  ! PROBLEM:  how to ensure that the frozen fraction does not exceed 1 ?
  
      if(n_in_dust + n_in_bio + n_in_soot + n_in_solO > 0.) then
        
        CIHENC_dust = n_in_dust - nin_a_nuc_dust 
        if(CIHENC_dust < 0.) CIHENC_dust = 0.
        
        CIHENC_soot = n_in_soot - nin_a_nuc_soot 
        if(CIHENC_soot < 0.) CIHENC_soot = 0.
        
        CIHENC_bio = n_in_bio - nin_a_nuc_bio 
        if(CIHENC_bio < 0.) CIHENC_bio = 0.
        
        CIHENC_solO = n_in_solO - nin_a_nuc_solO 
        if(CIHENC_solO < 0.) CIHENC_solO = 0.
        
        n_iw =  n_iw + CIHENC_dust
        nin_a_nuc_dust = nin_a_nuc_dust + CIHENC_dust
        num_ic_dust_imm = num_ic_dust_imm + CIHENC_dust
        
        n_iw =  n_iw + CIHENC_soot
        nin_a_nuc_soot = nin_a_nuc_soot + CIHENC_soot
        num_ic_soot_imm = num_ic_soot_imm + CIHENC_soot
  
        n_iw =  n_iw + CIHENC_bio
        nin_a_nuc_bio = nin_a_nuc_bio + CIHENC_bio
        num_ic_bio_imm = num_ic_bio_imm + CIHENC_bio
        
        n_iw =  n_iw + CIHENC_solO
        nin_a_nuc_solO = nin_a_nuc_solO + CIHENC_solO
        num_ic_solO_imm = num_ic_solO_imm + CIHENC_solO
      endif
    endif
  endif
  
  dNall = dn_in_dust + dn_in_bio + n_in_soot + n_in_solO  !DONIF
  
  if  (( dNall > 0.) .and. (n_iw .gt. 0.0))   then
    Dsh=max(min(n_iw/dNall, SS_i), 0.005)
  else
    Dsh=0.005
  end if
  
END SUBROUTINE EMPIRICAL_PARAM_PHILLIPS

real function H_1(X, X_1, X_2, Hlo)
  real, intent(in) :: Hlo, X, X_1, X_2
  
  if(X >= X_2) H_1 = 1
  if(X <= X_1) H_1 = Hlo 
  if(X > X_1 .and. X < X_2) H_1 = (X - X_1)/(X_2 - X_1) 
  
  if( X_2 <= X_1) stop 91919
  
  return 
end function 


real function H_1_smooth(X, X_1, X_2, Hlo, Hhi)
  real, intent(in) :: Hlo, Hhi, X, X_1, X_2
  real :: a_0, a_1, a_2, a_3, A, B
  
  if(X >= X_2) H_1_smooth = Hhi
  if(X <= X_1) H_1_smooth = Hlo 
  
  if(X >= X_2) dH1smooth = 0.0
  if(X <= X_1) dH1smooth = 0.0 
  
  if(X > X_1 .and. X < X_2) then
    A = 6.*(Hlo - Hhi)/(X_2**3. - X_1**3. + 3.*(X_2*X_1*X_1 - X_1*X_2*X_2) )
    a_3 = (A/3.)
    a_2 = -(A/2.)*(X_1 + X_2)
    a_1 = A*X_2*X_1
    B = Hlo + A*(X_1**3.)/6. - A*X_1*X_1*X_2/2.	
    a_0 = B
    H_1_smooth = a_0 + a_1*X + a_2*X*X + a_3*X*X*X
    dH1smooth =  a_1 + 2.0*a_2*X + 3.0*a_3*X*X
  endif
  
  dH1smooth =min(dH1smooth , 1.0e6_wp)
  dH1smooth =max(dH1smooth , 1.0e-12_wp)
  
  if( X_2 <= X_1) stop 91919
  
  return 
end function 
END MODULE mo_art_nucleation_cold
