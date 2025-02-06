!
! mo_art_emission_plumerise
! This module provides the plume rise model for biomass burning parameterization.
! All changes made after initial creation are marked in comments. Main parts equals to the
! original version of S. Freitas.
!
! Plume rise model for vegetation fires (CPTEC/INPE 2005-2006,2009)                         !
! Refs.:                                                                                    !
! Freitas, S. R., K. M. Longo, J. Trentmann, D. Latham. Technical Note: Sensitivity         !
! of 1D smoke plume rise models to the inclusion of environmental wind drag.                !
! Atmospheric Chemistry  and Physics, 2010.                                                 !
!                                                                                           !
! Freitas, S. R., K. M. Longo, R. Chatfield, D. Latham, M. A. F. Silva Dias, M. O. Andreae, !
! E. Prins, J. C. Santos, R. Gielow and J. A. Carvalho Jr.: Including the sub-grid scale    !
! plume rise of vegetation fires in low resolution atmospheric transport models.            !
!  Atmospheric Chemistry and Physics,2007.                                                  !
!-                                                                                          !
! Freitas, S. R.; Longo, K. M.; M. Andreae. Impact of including the plume rise of vegetation!
! fires in numerical simulations of associated atmospheric pollutants. Geophys. Res. Lett., !
! 33, L17808, DOi:10.1029/2006GL026608, 2006.                                               !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
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

MODULE mo_art_emission_plumerise
! ICON
  USE mo_kind,                          ONLY: wp
  USE mtime,                            ONLY: datetime
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
! ART


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine = 'mo_art_emission_plumerise'

  PUBLIC :: smk_pr_driver


!!
!!-------------------------------------------------------------------------
!!

!JS!! First Compile Test
!JS!USE data_modelconfig,     ONLY: ke


  INTEGER, PARAMETER :: n_biome = 4
  INTEGER, PARAMETER :: tropical_forest = 1
  INTEGER, PARAMETER :: boreal_forest = 2
  INTEGER, PARAMETER :: savanna = 3
  INTEGER, PARAMETER :: grassland = 4
  INTEGER, PARAMETER :: nkp = 200, ntime = 200                !cw nkp war 200
 
  REAL(wp),DIMENSION(nkp) ::  w,t,theq,qv,qc,qh,qi,sc,vth,vti,rho,txs,qt,           &
                          &   est,qsat,qpas,qtotal,td,vel_p,rad_p
  REAL(wp),DIMENSION(nkp) ::  wc,wt,tt,qvt,qct,qht,qit,sct,vel_t,rad_t
  REAL(wp),DIMENSION(nkp) ::  dzm,dzt,vctr1,vctr2,vt3dc,vt3df,vt3dk,vt3dg,scr1,     &
                          &   vt3dj,vt3dn,vt3DO,vt3da,scr2,vt3db,vt3dd,vt3dl,vt3dm, &
                          &   vt3di,vt3dh,vt3de, rbuoy,dwdt_entr
  REAL(wp),DIMENSION(nkp) ::  the, thee,rhe,sce,tde ! environment at plume grid
  REAL(wp), DIMENSION(nkp) :: te, pe, upe, vpe, qvenv, vel_e, zt, zm !wt !cw:davor eine zeile drueber
  !                                                         andere auch noch allocatable machen?
  ! REAL(wp), DIMENSION(40) :: ucon, zcon, vcon, tmpcon, prcon, rvcon   !be careful with 40, should be ke, 
  !                                                                  JS: ported to smk_pr_driver() 
  REAL(wp),DIMENSION(1200) :: wcon,thtcon,zzcon,scon,urcon,tdcon ! environment at RAMS  grid
  REAL(wp) :: dz,dqsdz,visc(nkp),viscosity,tstpf   
  INTEGER :: n,nm1,l
  REAL(wp) :: advw,advt,advv,advc,advh,advi,cvh(nkp),cvi(nkp),adiabat,&
              wbar,alast(10),vhrel,virel  ! advection
  REAL(wp) :: zsurf,zbase,ztop
  INTEGER :: lbase
  REAL(wp) :: area,rsurf,alpha,radius(nkp)  ! entrain
  REAL(wp) :: heating(ntime),fmoist,bload,heat_fluxw   ! heating
  REAL(wp) :: dt,time,tdur
  INTEGER :: mintime,mdur,maxtime
  REAL(wp) :: ztop_(ntime)
  !turn ON (=1), OFF (=0) the env wind effect on plume rise
  INTEGER, PARAMETER :: wind_eff = 1
  INTEGER :: kmt, nk
 
!---------------------------------------------------------------------------
  !Module rconstants
  REAL(wp), PARAMETER ::                   &
       rgas     = 287.05_wp,               &
       cp       = 1004._wp,                &
       cv       = 717._wp,                 &
       rm       = 461._wp,                 &
       p00      = 1.e5_wp,                 &
       t00      = 273.16_wp,               &
       g        = 9.80_wp,                 &
       pi180    = 3.1415927_wp / 180._wp,  &
       pi4      = 3.1415927_wp * 4._wp,    &
       spcon    = 111120._wp,              &
       erad     = 6367000._wp,             &
       vonk     = 0.40_wp,                 &
       tkmin    = 5.e-4_wp,                &
       alvl     = 2.50e6_wp,               &
       alvi     = 2.834e6_wp,              &
       alli     = 0.334e6_wp,              &
       alvl2    = 6.25e12_wp,              &
       alvi2    = 8.032e12_wp,             &
       solar    = 1.3533e3_wp,             &
       stefan   = 5.6696e-8_wp,            &
       cww      = 4218._wp,                &
       c0       = 752.55_wp * 4.18684e4_wp,&
       viscos   = .15e-4_wp,               &
       rowt     = 1.e3_wp,                 &
       dlat     = 111120._wp,              &
       omega    = 7.292e-5_wp,             &
       rocp     = rgas / cp,               &
       p00i     = 1._wp / p00,             &
       cpor     = cp / rgas,               &
       rocv     = rgas / cv,               &
       cpi      = 1._wp / cp,              &
       cpi4     = 4._wp * cpi,             &
       cp253i   = cpi / 253._wp,           & 
       allii    = 1._wp / alli,            &
       aklv     = alvl / cp,               &
       akiv     = alvi / cp,               &
       gama     = cp / cv,                 &
       gg       = .5_wp * g,               &
       ep       = rgas / rm,               & 
       p00k     = 26.870941_wp,            &  !  = p00 ** rocp  
       p00ki    = 1._wp / p00k
!---------------------------------------------------------------------------

!JS!:
CONTAINS

SUBROUTINE  smk_pr_driver(ke,heat_flux_max,heat_flux_min,burnt_area,temp_in,pres_in,u_in, &
                            v_in,hml_in,qv_in,topo_in,ztop_out,zbot_out)
  IMPLICIT NONE

  !JS: Instead of USE ke: in ICON nlev is passed by SUBROUTINE PARAMETER, not as a global variable 
  INTEGER, INTENT(IN) :: ke
  REAL(wp), INTENT(IN) :: heat_flux_min, heat_flux_max, burnt_area
  REAL(wp),INTENT(IN) :: temp_in(:),pres_in(:),u_in(:),v_in(:),hml_in(:),qv_in(:),topo_in
  REAL(wp),INTENT(INOUT) :: ztop_out, zbot_out
  INTEGER :: m1,m2,m3,ia,iz,ja,jz,ibcon,mynum,i,j,k,&
          !   k_CO_smold,k_PM25_smold,
      &      imm,k2,kk  !kmt
   
  REAL(wp) :: dz_flam,rhodzi
  REAL(wp), DIMENSION(2) :: ztopmax
  !JS: for ICON they have to be defined here with correct DIMENSION 
  REAL(wp), POINTER :: ucon(:), zcon(:), vcon(:), tmpcon(:), prcon(:), rvcon(:)   
  
  !JS 
  ALLOCATE( ucon(ke) )
  ALLOCATE( zcon(ke) )
  ALLOCATE( vcon(ke) )
  ALLOCATE( tmpcon(ke) )
  ALLOCATE( prcon(ke) )
  ALLOCATE( rvcon(ke) )

 !- initial settings (should be changed for the coupling with 3d host model)
  ia=1;iz=1;ja=1;jz=1; m1=ke
 
 !-initialize several PARAMETERs (only need to be DOne at the 1st time)
  CALL zero_plumegen_coms(ucon,zcon,vcon,tmpcon,prcon,rvcon) !JS: passing PARAMETER

!*********** READ VARIABLES OF HOST MODEL **************************************
  DO k = 1,m1 ! z-direction of the host model 
    rvcon (k) = qv_in(k)   !fill with the water vapor mixing ratio (kg/kg) of the host model 
    zcon  (k) = hml_in(k)-topo_in  !fill with the cartesian height (meters) of the termodynamic grid point  
                           !(half-levels) above local surface--> muss full level sein, siehe Aufruf htint
    tmpcon(k) = temp_in(k)
    prcon (k) = pres_in(k) 
    ucon  (k) = u_in(k)    !fill with the zonal wind (m/s) of the host model
    vcon  (k) = v_in(k)    !fill with the meridional wind (m/s) of the host model
  ENDDO
!***********


  ! get envinronmental state (temp, water vapor mix ratio, ...)
  CALL get_env_condition(ke,kmt,ucon,zcon,vcon,tmpcon,prcon,rvcon) !JS ucon,... as well
  DO imm=1,2  !  min and max power allowed 
    !- get fire properties (burned area, plume radius, heating rates ...)
    CALL get_fire_properties(burnt_area,heat_flux_max,heat_flux_min,imm) 
    !------  generate plume rise    ------       
    CALL makeplume (kmt,ztopmax(imm),imm)    
    ! print*,burnt_area,ztopmax(imm)/1000.,heat_fluxW
    ztop_out = ztopmax(2) 
    zbot_out = ztopmax(1)        
  ENDDO ! ENDDO DO loop em imm (min and max power allowed

  ! JS Dealloc
  DEALLOCATE( ucon,zcon,vcon,tmpcon,prcon,rvcon )

END SUBROUTINE smk_pr_driver
!-------------------------------------------------------------------------

SUBROUTINE get_env_condition(ke,kmt,ucon,zcon,vcon,tmpcon,prcon,rvcon)

  IMPLICIT NONE
 
  REAL(wp),INTENT(IN) :: ucon(:),zcon(:),vcon(:),tmpcon(:),prcon(:),rvcon(:) !JS as passed PARAMETER
  INTEGER :: ke,k,kmt,nk,i,kk
  REAL(wp) :: znz,dummy
  REAL(wp) :: topo,es!,esat

  CALL set_grid ! define vertical grid of plume model
                ! zt(k) =  thermo and water levels
                ! zm(k) =  dynamical levels 

  !k2=40                         !cw k2 war 200
  znz=zcon(ke)
  DO k=nkp,1,-1
    IF(zt(k)<znz) EXIT !go to 13
  ENDDO
  IF (k==1) stop ' envir stop 12'
!  13 continue
  kmt=k


  kmt=nkp !cwalter
  nk = ke
  !nk=40  ! number of vertical levels of the sounding

  CALL htint(nk,  ucon,zcon,kmt,upe  ,zt)
  CALL htint(nk,  vcon,zcon,kmt,vpe  ,zt)
  CALL htint(nk, rvcon,zcon,kmt,qvenv,zt)
  CALL htint(nk,tmpcon,zcon,kmt,te   ,zt)
  CALL htint(nk, prcon,zcon,kmt,pe   ,zt)

  DO k=1,kmt
    !  PE esta em kPa  - ESAT DO RAMS esta em mbar = 100 Pa = 0.1 kPa
    !  ES       = 0.1*ESAT (TE(k)) !blob saturation vapor pressure, em kPa
    !  QSAT (k) = (.622 * ES) / (PE (k)*1.e-3 - ES) !saturation lwc g/g
    qvenv(k)=MAX(qvenv(k),1e-8_wp)
  ENDDO

  DO k=1,kmt
    vel_e(k) = SQRT(upe(k)**2+vpe(k)**2)
  ENDDO

  !ewe - env wind effect
  IF(wind_eff < 1)  vel_e = 0._wp

  !--------- converte press de Pa para kPa para uso modelo de plumerise
  DO k=1,kmt
   pe(k) = pe(k)*1.e-3_wp
  ENDDO 

END SUBROUTINE get_env_condition
!-------------------------------------------------------------------------

SUBROUTINE set_grid()
 
  IMPLICIT NONE
  INTEGER :: k,mzp

  dz=100._wp ! set constant grid spacing of plume grid model(meters)

  mzp=nkp
  zt(1) = zsurf
  zm(1) = zsurf
  zt(2) = zt(1) + 0.5_wp*dz
  zm(2) = zm(1) + dz
  DO k=3,mzp
    zt(k) = zt(k-1) + dz ! thermo and water levels
    zm(k) = zm(k-1) + dz ! dynamical levels    
  ENDDO
  !print*,zsurf
  !Print*,zt(:)

  DO k = 1,mzp-1
    dzm(k) = 1._wp / (zt(k+1) - zt(k))
  ENDDO 
  dzm(mzp)=dzm(mzp-1)

  DO k = 2,mzp
    dzt(k) = 1._wp / (zm(k) - zm(k-1))
  ENDDO
  dzt(1) = dzt(2) * dzt(2) / dzt(3)
  !   dzm(1) = 0.5/dz
  !   dzm(2:mzp) = 1./dz

END SUBROUTINE set_grid
!-------------------------------------------------------------------------

SUBROUTINE get_fire_properties(burnt_area,heat_flux_max,heat_flux_min,imm)
  
  IMPLICIT NONE

  INTEGER :: moist, i, icount, imm
  REAL(wp) :: bfract, effload, heat, hinc, burnt_area, heat_flux_min, heat_flux_max
  !REAL(wp),    DIMENSION(2,4) :: heat_flux
  !data heat_flux/  &
  !---------------------------------------------------------------------
  !  heat flux      !IGBP Land Cover      ! 
  ! min  ! max      !Legend and           ! reference
  !    kW/m^2       !description          ! 
  !--------------------------------------------------------------------
  ! 30.0,  80.0,   &! Tropical Forest         ! igbp 2 & 4
  ! 30.0,   80.0,   &! Boreal forest           ! igbp 1 & 3
  !  4.4,  23.0,   &! cerrado/woody savanna   | igbp  5 thru 9
  !  3.3,   3.3    /! Grassland/cropland      ! igbp 10 thru 17
  !--------------------------------------------------------------------


  !-- fire at the surface
  area = burnt_area   ! area of burn, m^2
  !- heat flux

  IF (imm==1) THEN
    heat_fluxW = heat_flux_min  !cw
  ELSE
    heat_fluxW = heat_flux_max
  ENDIF

  mdur = 200       ! duration of burn, minutes
  bload = 10._wp      ! total loading, kg/m**2 
  moist = 10       ! fuel moisture, %. average fuel moisture,percent dry

  maxtime =mdur-1  ! model time, min

  !heat = 21.e6    !- joules per kg of fuel consumed                   
  !heat = 15.5e6   !- joules/kg - cerraDO/savannah
  heat = 19.3e6_wp    !- joules/kg - tropical forest (mt)
  !alpha = 0.1      !- entrainment constant
  alpha = 0.05_wp      !- entrainment constant

  ! ******************** fix up inputs *********************************
                                             
  IF (MOD (maxtime, 2) /=0) maxtime = maxtime+1  !make maxtime even
                                                  
  maxtime = maxtime * 60  ! and put in seconds

  rsurf = SQRT (area / 3.14159_wp) !- entrainment surface radius (m)

  fmoist   = moist / 100._wp       !- fuel moisture fraction

  ! calculate the energy flux and water content at lboundary.
  ! fills heating() on a minute basis. could ask for a file at this po
  ! in the program. whatever is input has to be adjusted to a one
  ! minute timescale.
                        
  DO i = 1, ntime         !- make sure of energy release
    heating (i) = 0.0001_wp  !- avoid possible divide by 0
  ENDDO  
                                  
  tdur = mdur * 60._wp       !- number of seconds in the burn

  bfract = 1._wp             !- combustion factor

  effload = bload * bfract  !- patchy burning
  
  !     spread the burning evenly over the interval
  !     except for the first few minutes for stability
  icount = 1  
  !
  IF(mdur > ntime) STOP 'Increase time duration (ntime) in min - see file "plumerise_mod.f90"'

  DO WHILE (icount<=mdur)                             
 !  HEATING (ICOUNT) = HEAT * EFFLOAD / TDUR  ! W/m**2 
 !  HEATING (ICOUNT) = 80000.  * 0.55         ! W/m**2 
    heating (icount) = heat_fluxW  * 0.55_wp     ! W/m**2 (0.55 converte para energia convectiva)
    icount = icount + 1  
  ENDDO  
  !     ramp for 5 minutes
  hinc = heating (1) / 4._wp  
  heating (1) = 0.1_wp 
  heating (2) = hinc  
  heating (3) = 2._wp * hinc  
  heating (4) = 3._wp * hinc 

END SUBROUTINE get_fire_properties
!-------------------------------------------------------------------------------

SUBROUTINE makeplume (kmt,ztopmax,imm)  

! *********************************************************************
!
!   EQUATION SOURCE--Kessler Met.Monograph No. 32 V.10 (K)
!    Alan Weinstein, JAS V.27 pp 246-255. (W),
!    Ogura and Takahashi, Monthly Weather Review V.99,pp895-911 (OT)
!    Roger Pielke,Mesoscale Meteorological Modeling,Academic Press,1984
!    Originally developed by: Don Latham (USFS)
!
!
! ************************ VARIABLE ID ********************************
!
!     DT=COMPUTING TIME INCREMENT (SEC)
!     DZ=VERTICAL INCREMENT (M)
!     LBASE=LEVEL ,CLOUD BASE
!
!     CONSTANTS:
!       G = GRAVITATIONAL ACCELERATION 9.80796 (M/SEC/SEC).
!       R = DRY AIR GAS CONSTANT (287.04E6 JOULE/KG/DEG K)
!       CP = SPECIFIC HT. (1004 JOULE/KG/DEG K)
!       HEATCOND = HEAT OF CONDENSATION (2.5E6 JOULE/KG)
!       HEATFUS = HEAT OF FUSION (3.336E5 JOULE/KG)
!       HEATSUBL = HEAT OF SUBLIMATION (2.83396E6 JOULE/KG)
!       EPS = RATIO OF MOL.WT. OF WATER VAPOR TO THAT OF DRY AIR (0.622)
!       DES = DIFFERENCE BETWEEN VAPOR PRESSURE OVER WATER AND ICE (MB)
!       TFREEZE = FREEZING TEMPERATURE (K)
!
!
!     PARCEL VALUES:
!       T = TEMPERATURE (K)
!       TXS = TEMPERATURE EXCESS (K)
!       QH = HYDROMETEOR WATER CONTENT (G/G DRY AIR)
!       QHI = HYDROMETEOR ICE CONTENT (G/G DRY AIR)
!       QC = WATER CONTENT (G/G DRY AIR)
!       QVAP = WATER VAPOR MIXING RATIO (G/G DRY AIR)
!       QSAT = SATURATION MIXING RATIO (G/G DRY AIR)
!       RHO = DRY AIR DENSITY (G/M**3) MASSES = RHO*Q'S IN G/M**3
!       ES = SATURATION VAPOR PRESSURE (kPa)
!
!     ENVIRONMENT VALUES:
!       TE = TEMPERATURE (K)
!       PE = PRESSURE (kPa)
!       QVENV = WATER VAPOR (G/G)
!       RHE = RELATIVE HUMIDITY FRACTION (e/esat)
!       DNE = dry air density (kg/m^3)
!
!     HEAT VALUES:
!       HEATING = HEAT OUTPUT OF FIRE (WATTS/M**2)
!       MDUR = DURATION OF BURN, MINUTES
!
!       W = VERTICAL VELOCITY (M/S)
!       RADIUS=ENTRAINMENT RADIUS (FCN OF Z)
!   RSURF = ENTRAINMENT RADIUS AT GROUND (SIMPLE PLUME, TURNER)
!   ALPHA = ENTRAINMENT CONSTANT
!       MAXTIME = TERMINATION TIME (MIN)
!
!
!**********************************************************************
!**********************************************************************               
 
  IMPLICIT NONE   

  CHARACTER (LEN=10) :: varn
  INTEGER :: iconv, itime, k, kk, kkmax, deltak, kmt, &
           & nrectotal, i_micro, n_sub_step, imm
  REAL(wp) :: ztopmax, wmax, rmaxtime, es, heat, dt_save !esat
  REAL(wp), PARAMETER :: vc = 5.107387_wp  
  REAL(wp), PARAMETER :: g = 9.80796_wp,r = 287.04_wp,cp = 1004._wp,eps = 0.622_wp,tmelt = 273.3_wp
  REAL(wp), PARAMETER :: heatsubl = 2.834e6_wp, heatfus = 3.34e5_wp, heatcond =2.501e6_wp
  REAL(wp), PARAMETER :: tfreeze = 269.3_wp
  ! ******************* SOME CONSTANTS **********************************
  !      XNO=10.0E06 median volume diameter raindrop (K table 4)
  !      VC = 38.3/(XNO**.125) mean volume fallspeed eqn. (K)

  tstpf = 2.0_wp     !- timestep factor
  nrectotal=150

  !*************** PROBLEM SETUP AND INITIAL CONDITIONS *****************
  mintime = 1  
  ztopmax = 0._wp 
  ztop    = 0._wp 
  time    = 0._wp  
  dt      = 1._wp
  wmax    = 1._wp 
  kkmax   = 10
  deltaK  = 20
  l       = 1   ! l initialization
  viscosity = 500._wp!- viscosity constant (original value: 0.001)

  !--- initialization
  CALL initial(kmt)      

  ! ******************* model evolution ******************************
  rmaxtime = REAL(maxtime,wp)

  DO WHILE (time<=rmaxtime)  !beginning of time loop

    !-- set model top integration
    nm1 = MIN(kmt, kkmax + deltak)
!-- set timestep
    !dt = (zm(2)-zm(1)) / (tstpf * wmax) 
    dt = MIN(5._wp,((zm(2)-zm(1)) / (tstpf * wmax))) 
!-- elapsed time, sec
    time = time+dt 
!-- elapsed time, minutes                                      
    mintime = 1 + INT (time) / 60     
    wmax = 1._wp  !no zeroes allowed.

!************************** BEGIN SPACE LOOP **************************

!-- zerout all model tendencies
    CALL tend0_plumerise

!-- surface bounday conditions (k=1)
    l=1
    CALL lbound_mtt()

!-- dynamics for the level k>1 

!-- W advection 
    CALL vel_advectc_plumerise(nm1,wc,wt,rho,dzm)
  
!-- scalars advection
    CALL scl_advectc_plumerise('SC',nm1)
    !CALL scl_advectc_plumerise2('SC',NM1)

!-- scalars entrainment, adiabatic
    CALL scl_misc(nm1)

!-- scalars dinamic entrainment
    CALL  scl_dyn_entrain(nm1)
    
!-- gravity wave damping using Rayleigh friction layer fot T
    CALL damp_grav_wave(1,nm1,deltak,zt,w,t,tt,qv,qh,qi,qc,te,pe,qvenv, &
                      & vel_p,vel_t,vel_e)

!-- microphysics
   !goto 101 ! bypass microphysics
    dt_save=dt
    n_sub_step=3
    dt=dt/REAL(n_sub_step,wp)

    DO i_micro=1,n_sub_step
!-- sedim ?
      CALL fallpart(nm1)
!-- microphysics
      DO l=2,nm1-1
        wbar    = 0.5_wp*(w(l)+w(l-1))
        es      = 0.1_wp*esat(t(l)) !blob saturation vapor pressure, em kpa
!   print*,'blob741',t(l),l
        qsat(l) = (eps * es) / (pe(l) - es)  !blob saturation lwc g/g dry air
        est (l) = es  
        rho (l) = 3483.8_wp * pe (l) / t (l) ! air parcel density , g/m**3
!srf18jun2005
!   if (w(l) >= 0.) dqsdz = (qsat(l  ) - qsat(l-1)) / (zt(l  ) -zt(l-1))
!   if (w(l) < 0.) dqsdz = (qsat(l+1) - qsat(l  )) / (zt(l+1) -zt(l  ))
        IF (w(l) >= 0._wp) THEN 
          dqsdz = (qsat(l+1) - qsat(l-1)) / (zt(l+1 )-zt(l-1))
        ELSE
          dqsdz = (qsat(l+1) - qsat(l-1)) / (zt(l+1) -zt(l-1))
        ENDIF 
    
        CALL waterbal  
      ENDDO
    ENDDO
    dt=dt_save

!    101 continue
    
!-- W-viscosity for stability 
    CALL visc_w(nm1,deltak,kmt)

!-- update scalars
    CALL update_plumerise(nm1,'S')
       !print*,'wi apos update=',w(1:nm1)
       !print*,'Ti apos update=',T(1:nm1)

    CALL hadvance_plumerise(1,nm1,wc,wt,w,mintime) 

!-- Buoyancy
    CALL buoyancy_plumerise(nm1, t, te, qv, qvenv, qh, qi, qc, wt, scr1)
 
!-- Entrainment 
    CALL entrainment(nm1,w,wt,radius,alpha)

!-- update W
    CALL update_plumerise(nm1,'W')

    CALL hadvance_plumerise(2,nm1,wc,wt,w,mintime) 


!-- misc
    DO k=2,nm1
!    pe esta em kpa  - esat do rams esta em mbar = 100 Pa = 0.1 kpa
      es       = 0.1_wp*esat(t(k)) !blob saturation vapor pressure, em kPa
!    rotina do plumegen calcula em kPa
      qsat(k) = (eps * es) / (pe(k) - es)  !blob saturation lwc g/g dry air
      est (k) = es  
      txs (k) = t(k) - te(k)
      rho (k) = 3483.8_wp * pe (k) / t (k) ! air parcel density , g/m**3
                                       ! no pressure diff with radius
                       
      IF((ABS(wc(k)))>wmax) wmax = ABS(wc(k)) ! keep wmax largest w

!srf-27082005
!     IF((ABS(wt(k)))>wtmax) wtmax = ABS(wt(k)) ! keep wmax largest w
    ENDDO  

! Gravity wave damping using Rayleigh friction layer for W
    CALL damp_grav_wave(2,nm1,deltak,zt,w,t,tt,qv,qh,qi,qc,te,pe,qvenv,vel_p,vel_t,vel_e)

!- update radius (para versao original DO modelo, comente as 3 linhas abaixo
    DO k=2,nm1
      radius(k) = rad_p(k)
    ENDDO

!-- try to find the plume top (above surface height)
    kk = 1
    DO WHILE (w (kk) > 1._wp)  
      kk = kk + 1  
      ztop =  zm(kk) 
      !print*,'W=',w (kk)
    ENDDO  
        
    ztop_(mintime) = ztop
    ztopmax = MAX (ztop, ztopmax) 
    kkmax   = MAX (kk  , kkmax  ) 
!    print * ,'ztopmax=', mintime,'mn ',ztop_(mintime), ztopmax!,wtmax

! IF the solution is going to a stationary phase, exit
    IF(mintime > 40) THEN
      IF( ABS(ztop_(mintime)-ztop_(mintime-10)) < dz ) EXIT   !??????
    ENDIF
 
  ENDDO   !DO next timestep

END SUBROUTINE makeplume

!-------------------------------------------------------------------------------
SUBROUTINE burn(eflux, water)  
    
  IMPLICIT NONE
  !calculates the energy flux and water content at lboundary                               
  !REAL(wp), PARAMETER :: HEAT = 21.E6 !Joules/kg
  !REAL(wp), PARAMETER :: HEAT = 15.5E6 !Joules/kg - cerraDO
  REAL(wp), PARAMETER :: heat = 19.3E6_wp !Joules/kg - floresta em Alta Floresta (MT)
  REAL(wp) :: eflux,water

  ! The emission factor for water is 0.5. The water produced, in kg,
  ! is then  fuel mass*0.5 + (moist/100)*mass per square meter.
  ! The fire burns for DT out of TDUR seconds, the total amount of
  ! fuel burned is AREA*BLOAD*(DT/TDUR) kg. this amount of fuel is
  ! considered to be spread over area AREA and so the mass burned per
  ! unit area is BLOAD*(DT/TDUR), and the rate is BLOAD/TDUR.

  IF (time>tdur) THEN !is the burn over?   
    eflux = 0.000001_wp    !prevent a potential divide by zero
    water = 0._wp  
    RETURN  
  ELSE                   
    eflux = heating (mintime)                                ! watts/m**2                                                   
    water = eflux * (dt / heat) * (0.5 + fmoist)             ! kg/m**2 
    water = eflux * (dt / heat) * (0.5_wp + fmoist) /0.55_wp ! kg/m**2 
    water = water * 1000._wp                                 ! g/m**2
  ENDIF  
    
END SUBROUTINE burn
!-------------------------------------------------------------------------------
SUBROUTINE lbound_mtt ()  

  ! ********** BOUNDARY CONDITIONS AT ZSURF FOR PLUME AND CLOUD ********
  ! source of equations: J.S. Turner Buoyancy Effects in Fluids
  !                      Cambridge U.P. 1973 p.172,
  !                      G.A. Briggs Plume Rise, USAtomic Energy Commissio
  !                      TID-25075, 1969, P.28
  ! fundamentally a point source below ground. at surface, this produces
  ! a velocity w(1) and temperature T(1) which vary with time. There is
  ! also a water load which will first saturate, then remainder go into
  ! QC(1).
  ! EFLUX = energy flux at ground,watt/m**2 for the last DT
    
  IMPLICIT NONE

  REAL(wp), PARAMETER :: g = 9.80796_wp,r = 287.04_wp,cp = 1004.6_wp,eps = 0.622_wp,tmelt = 273.3_wp
  REAL(wp), PARAMETER :: tfreeze = 269.3_wp, pi = 3.14159_wp, e1 = 1._wp/3._wp, e2 = 5._wp/3._wp
  REAL(wp) :: es, eflux, water,  pres, c1, c2, f, zv,  denscor, xwater
             
  qh (1) = qh (2)   !soak up hydrometeors
  qi (1) = qi (2)              
  qc (1) = 0._wp       !no cloud here

  CALL burn (eflux, water)  

!  calculate PARAMETERs at boundary from a virtual buoyancy point source

  pres = pe (1) * 1000._wp   !need pressure in n/m**2
                              
  c1 = 5._wp / (6._wp * alpha)  !alpha is entrainment constant

  c2 = 0.9_wp * alpha  

  f = eflux / (pres * cp * pi)  
                            
  f = g * r * f * area  !buoyancy flux
                
  zv = c1 * rsurf  !virtual boundary height
                                  
  w (1) = c1 * ( (c2 * f) **e1) / zv**e1  !boundary velocity
                                        
  denscor = c1 * f / g / (c2 * f) **e1 / zv**e2   !density correction

  t (1) = te (1) / (1._wp - denscor)    !temperature of virtual plume at zsurf

  wc(1) = w(1)
  
  vel_p(1) = 0._wp
  rad_p(1) = rsurf

  sc(1) = 1000._wp!sce(1)+f/1000.*dt  ! gas/particle (g/g)

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !     match dw/dz,dt/dz at the boundary. f is conserved.
  !wbar = w (1) * (1. - 1. / (6. * zv) )  
  !advw = wbar * w (1) / (3. * zv)  
  !advt = wbar * (5. / (3. * zv) ) * (denscor / (1. - denscor) )  
  !advc = 0.  
  !advh = 0.  
  !advi = 0.  
  !adiabat = - wbar * g / cp  
   vth (1) = - 4._wp  
   vti (1) = - 3._wp  
   txs (1) = t (1) - te (1)  

   visc (1) = viscosity  

   rho (1) = 3483.8_wp * pe (1) / t (1)   !air density at level 1, g/m**3

   xwater = water / (w (1) * dt * rho (1) )   !firewater mixing ratio
                                            
   qv (1) = xwater + qvenv (1)  !plus what's already there 

!  pe esta em kpa  - esat do rams esta em mbar = 100 pa = 0.1 kpa
   es       = 0.1_wp*esat(t(1)) !blob saturation vapor pressure, em kpa

   est  (1)  = es                                  
   qsat (1) = (eps * es) / (pe (1) - es)   !blob saturation lwc g/g dry air
  
   IF (qv (1) > qsat (1) ) THEN  
       qc (1) = qv   (1) - qsat (1) + qc (1)  !remainder goes into cloud drops
       qv (1) = qsat (1)  
   ENDIF  

   CALL waterbal  

END SUBROUTINE lbound_mtt
!-------------------------------------------------------------------------------
SUBROUTINE initial ( kmt)  

! ************* SETS UP INITIAL CONDITIONS FOR THE PROBLEM ************
  IMPLICIT NONE 

  REAL(wp), PARAMETER :: tfreeze = 269.3_wp
  INTEGER :: isub, k, n1, n2, n3, lbuoy, itmp, isubm1, kmt
  REAL(wp) :: xn1, xi, es
 
  n=kmt
! initialize temperature structure,to the END of equal spaced sounding,
  DO k = 1, n             
    txs (k) = 0.0_wp  
    w (k) = 0.0_wp             
    t (k) = te(k) !blob set to environment
!print*,'t(',k,')=',t(k)          
    wc(k) = 0.0_wp
    wt(k) = 0.0_wp
    qv(k) = qvenv (k)                      
    vth(k) = 0._wp      !initial rain velocity = 0                       
    vti(k) = 0._wp      !initial ice  velocity = 0                       
    qh(k) = 0._wp           !no rain                 
    qi(k) = 0._wp           !no ice                  
    qc(k) = 0._wp      !no cloud drops                      
!  pe esta em kpa  - esat do rams esta em mbar = 100 pa = 0.1 kpa
    es       = 0.1_wp*esat(t(k)) !blob saturation vapor pressure, em kpa
    est  (k) = es  
    qsat (k) = (.622_wp * es) / (pe (k) - es) !saturation lwc g/g
    rho  (k) = 3483.8_wp * pe (k) / t (k)   !dry air density g/m**3     
    vel_p(k) = 0._wp
    rad_p(k) = 0._wp
  ENDDO  

! Initialize the entrainment radius, Turner-style plume
  radius(1) = rsurf
  DO k=2,N
    radius(k) = radius(k-1)+(6._wp/5._wp)*alpha*(zt(k)-zt(k-1))
    rad_p(k)  = radius(k)
  ENDDO
    
  rad_p(1) = rsurf

!  Initialize the viscosity
  visc (1) = viscosity
  DO k=2,N
    visc (k) = MAX(1.e-3_wp,visc(k-1) - 1._wp* viscosity/REAL(nkp,wp))
  ENDDO

  CALL lbound_mtt()

END SUBROUTINE initial
!-------------------------------------------------------------------------------

SUBROUTINE damp_grav_wave(ifrom,nm1,deltak,zt,w,t,tt,qv,qh,qi,qc,te,pe,qvenv,vel_p,vel_t,vel_e)
  IMPLICIT NONE

  INTEGER :: nm1, ifrom, deltak

  REAL(wp), DIMENSION(nm1) :: w,t,tt,qv,qh,qi,qc,te,pe,qvenv,dummy,zt,vel_p,vel_t,vel_e

  IF(ifrom==1) THEN
    CALL friction(ifrom,nm1,deltak,zt,t,tt,te)
! CALL friction(ifrom,nm1,deltak,dt,zt,zm,vel_p,vel_t,vel_e)
! CALL friction(ifrom,nm1,dt,zt,zm,qv,qvt,qvenv)
    RETURN
  ENDIF 

  dummy(:) = 0._wp
  IF(ifrom==2) CALL friction(ifrom,nm1,deltak,zt,w,dummy ,dummy)
!CALL friction(ifrom,nm1,dt,zt,zm,qi,qit ,dummy)
!CALL friction(ifrom,nm1,dt,zt,zm,qh,qht ,dummy)
!CALL friction(ifrom,nm1,dt,zt,zm,qc,qct ,dummy)
END SUBROUTINE damp_grav_wave
!-------------------------------------------------------------------------------
SUBROUTINE friction(ifrom,nm1,deltak,zt,var1,vart,var2)
  IMPLICIT NONE

  INTEGER :: k, nfpt, kf, nm1, ifrom, deltak
  REAL(wp), DIMENSION(nm1) :: var1, var2, vart, zt
  REAL(wp) zmkf,ztop,distim,c1,c2
  !nfpt=50
  !kf = nm1 - nfpt
  !kf = nm1 - int(deltak/2) ! orig
  kf = nm1 - INT(deltak)
  !IF( kf < 10) return !ver necessidade

  zmkf = zm(kf) !old: float(kf )*dz
  ztop = zm(nm1)

  !distim = 60. ! orig
  distim = MIN(3._wp*dt,60._wp)

  c1 = 1._wp / (distim * (ztop - zmkf))
  c2 = dt * c1

  IF(ifrom == 1) THEN  
    DO k = nm1,2,-1
      IF (zt(k) <= zmkf) CYCLE  ! exit ???
        vart(k) = vart(k)   + c1 * (zt(k) - zmkf)*(var2(k) - var1(k))
    ENDDO
  ELSEIF(ifrom == 2) THEN
    DO k = nm1,2,-1
      IF (zt(k) <= zmkf) CYCLE  ! exit ???
        var1(k) =  var1(k) + c2 * (zt(k) - zmkf)*(var2(k) - var1(k))
    ENDDO
  ENDIF
END SUBROUTINE friction
!-------------------------------------------------------------------------------

SUBROUTINE vel_advectc_plumerise(m1,wc,wt,rho,dzm)

  IMPLICIT NONE

  INTEGER :: k,m1
  REAL(wp), DIMENSION(m1) :: wc,wt,flxw,dzm,rho
  REAL(wp), DIMENSION(m1) :: dn0 ! var local
  REAL(wp) :: c1z

  !dzm(:)= 1./dz

  dn0(1:m1)=rho(1:m1)*1.e-3_wp ! converte de cgs para mks

  flxw(1) = wc(1) * dn0(1) 

  DO k = 2,m1-1
    flxw(k) = wc(k) * .5_wp * (dn0(k) + dn0(k+1))
  ENDDO

  ! Compute advection contribution to W tendency

  c1z = .5_wp 

  DO k = 2,m1-2
    wt(k) = wt(k)  &
      &   + c1z * dzm(k) / (dn0(k) + dn0(k+1))             &
      &   * ( (flxw(k) + flxw(k-1))  * (wc(k) + wc(k-1))   &
      &   - (flxw(k) + flxw(k+1))  * (wc(k) + wc(k+1))     &
      &   + (flxw(k+1) - flxw(k-1)) * 2._wp * wc(k)       )
  ENDDO

END SUBROUTINE vel_advectc_plumerise
!-------------------------------------------------------------------------------

SUBROUTINE hadvance_plumerise(iac,m1,wc,wt,w,mintime)

  IMPLICIT NONE

  INTEGER :: k,iac
  INTEGER :: m1,mintime
  REAL(wp), DIMENSION(m1) :: dummy, wc,wt,w
  REAL(wp) :: eps
  !     It is here that the Asselin filter is applied.  For the velocities
  !     and pressure, this must be done in two stages, the first when
  !     IAC=1 and the second when IAC=2.
  eps = .2_wp
  IF(mintime == 1) eps=0.5_wp

  !     For both IAC=1 and IAC=2, CALL PREDICT for U, V, W, and P.

  CALL predict_plumerise(m1,wc,w,wt,dummy,iac,2._wp*dt,eps)

END SUBROUTINE hadvance_plumerise
!-------------------------------------------------------------------------------

SUBROUTINE predict_plumerise(npts,ac,ap,fa,af,iac,dtlp,epsu)
  IMPLICIT NONE

  INTEGER :: npts,iac,m
  REAL(wp) :: epsu
  REAL(wp),INTENT(IN) :: dtlp
  REAL(wp), DIMENSION(*) :: ac,ap,fa,af

  !     For IAC=3, this routine moves the arrays AC and AP forward by
  !     1 time level by adding in the prescribed tendency. It also
  !     applies the Asselin filter given by:

  !              {AC} = AC + EPS * (AP - 2 * AC + AF)

  !     where AP,AC,AF are the past, current and future time levels of A.
  !     All IAC=1 DOes is to perform the {AC} calculation without the AF
  !     term present.  IAC=2 completes the calculation of {AC} by adding
  !     the AF term only, and advances AC by filling it with input AP
  !     values which were already updated in ACOUSTC.

  IF (iac == 1) THEN
    DO m = 1,npts
      ac(m) = ac(m) + epsu * (ap(m) - 2._wp * ac(m))
    ENDDO
    RETURN
  ELSEIF (iac == 2) THEN
    DO m = 1,npts
      af(m) = ap(m)
      ap(m) = ac(m) + epsu * af(m)
    ENDDO
!ELSEIF (iac == 3) then
!   DO m = 1,npts
!      af(m) = ap(m) + dtlp * fa(m)
!   ENDDO
!   IF (ngrid == 1 .and. ipara == 0) CALL cyclic(nzp,nxp,nyp,af,'T')
!   DO m = 1,npts
!      ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
!   ENDDO
  ENDIF

  DO m = 1,npts
    ac(m) = af(m)
  ENDDO
END SUBROUTINE predict_plumerise
!-------------------------------------------------------------------------------
!
SUBROUTINE  buoyancy_plumerise(m1, t, te, qv, qvenv, qh, qi, qc, wt, scr1)

  IMPLICIT NONE

  INTEGER :: k,m1
  REAL(wp), PARAMETER :: g = 9.8_wp, eps = 0.622_wp, gama = 0.5_wp ! mass virtual coeff.
  REAL(wp), DIMENSION(m1) :: t, te, qv, qvenv, qh, qi, qc, wt, scr1
  REAL(wp) :: tv,tve,qwtotl,umgamai
  REAL(wp), PARAMETER :: mu = 0.15_wp 

  !- orig
  umgamai = 1._wp/(1._wp+gama) ! compensa a falta DO termo de aceleracao associaDO `as
                               ! das pertubacoes nao-hidrostaticas no campo de pressao

  !- new                 ! Siesbema et al, 2004
  !umgamai = 1./(1.-2.*mu)

  DO k = 2,m1-1

    tv =   t(k) * (1._wp + (qv(k)   /eps))/(1._wp + qv(k)   )  !blob virtual temp.                                               
    tve = te(k) * (1._wp + (qvenv(k)/eps))/(1._wp + qvenv(k))  !and environment

    qwtotl = qh(k) + qi(k) + qc(k)                       ! qwtotl*g is drag
!- orig
   !scr1(k)= g*( umgamai*(  tv - tve) / tve   - qwtotl) 
    scr1(k)= g*  umgamai*( (tv - tve) / tve   - qwtotl) 

    !IF(k < 10)print*,'BT',k,TV,TVE,TVE,QWTOTL
  ENDDO

  DO k = 2,m1-2
!srf- just for output
    rbuoy(k)=0.5_wp*(scr1(k)+scr1(k+1))

    wt(k) = wt(k)+0.5_wp*(scr1(k)+scr1(k+1))
!   print*,'W-BUO',k,wt(k),scr1(k),scr1(k+1)
  ENDDO
!srf- just for output
  rbuoy(1)=rbuoy(2)
    
END SUBROUTINE  buoyancy_plumerise
!-------------------------------------------------------------------------------

SUBROUTINE entrainment(m1,w,wt,radius,alpha)

  IMPLICIT NONE

  INTEGER :: k,m1
  REAL(wp), DIMENSION(m1) :: w,wt,radius
  REAL(wp) :: dmdtm,alpha,wbar,radius_bar,umgamai,dyn_entr
  REAL(wp), PARAMETER :: mu = 0.15_wp ,gama = 0.5_wp ! mass virtual coeff.

  !- new - Siesbema et al, 2004
  !umgamai = 1./(1.-2.*mu)

  !- orig
  umgamai = 1._wp/(1._wp+gama) ! compensa a falta DO termo de aceleracao associaDO `as
                         !  pertubacoes nao-hidrostaticas no campo de pressao

!-- ALPHA/RADIUS(L) = (1/M)DM/DZ  (W 14a)
  DO k=2,m1-1
!-- for W: WBAR is only W(k)
!    WBAR=0.5*(W(k)+W(k-1))           
    wbar=w(k)          
    radius_bar = 0.5_wp*(radius(k) + radius(k-1))

! orig plump model
   ! dmdtm =           2. * alpha * abs (wbar) / radius_bar  != (1/m)dm/dt
    dmdtm = umgamai * 2._wp * alpha * ABS (wbar) / radius_bar  != (1/m)dm/dt

!--  dmdtm*w(l) entrainment,
    wt(k) = wt(k)  - dmdtm*ABS (wbar) !- dmdtm*abs (wbar)*1.875*0.5

     !if(vel_p (k) - vel_e (k) > 0.) cycle

!-   dynamic entrainment
    dyn_entr =  (2._wp/3.1416_wp)*0.5_wp*ABS (vel_p(k)-vel_e(k)+vel_p(k-1)-vel_e(k-1)) /radius_bar
    wt(k) = wt(k)  - dyn_entr*ABS (wbar)

!- entraiment acceleration for output only
    dwdt_entr(k) =  - dmdtm*ABS (wbar)- dyn_entr*ABS (wbar)
  ENDDO
!- entraiment acceleration for output only
  dwdt_entr(1) = dwdt_entr(2)


END SUBROUTINE  entrainment
!     ****************************************************************

SUBROUTINE scl_misc(m1)

  IMPLICIT NONE

  REAL(wp), PARAMETER :: g = 9.81_wp, cp=1004._wp
  INTEGER :: m1,k
  REAL(wp) :: dmdtm

  DO k=2,m1-1
    wbar    = 0.5_wp*(w(k)+w(k-1))  
!-- dry adiabat
    adiabat = - wbar * g / cp 
      
!-- entrainment     
    dmdtm = 2._wp * alpha * ABS (wbar) / radius (k)  != (1/m)dm/dt
      
!-- tendency temperature = adv + adiab + entrainment
    tt(k) = tt(k) + adiabat - dmdtm * ( t  (k) -    te (k) ) 

!-- tendency water vapor = adv  + entrainment
    qvt(k) = qvt(k)         - dmdtm * ( qv (k) - qvenv (k) )

    qct(k) = qct(k)         - dmdtm * ( qc (k)  )
    qht(k) = qht(k)         - dmdtm * ( qh (k)  )
    qit(k) = qit(k)         - dmdtm * ( qi (k)  )

!-- tendency horizontal speed = adv  + entrainment
    vel_t(k) = vel_t(k)     - dmdtm * ( vel_p (k) - vel_e (k) )

!-- tendency horizontal speed = adv  + entrainment
    rad_t(k) = rad_t(k)     + 0.5_wp*dmdtm*(6._wp/5._wp)*radius (k)

!-- tendency gas/particle = adv  + entrainment
    sct(k) = sct(k)         - dmdtm * ( sc (k) -   sce (k) )

  ENDDO
END SUBROUTINE scl_misc

!     ****************************************************************

SUBROUTINE scl_dyn_entrain(m1)

  IMPLICIT NONE

  REAL(wp), PARAMETER :: g = 9.81_wp, cp=1004._wp, pi=3.1416_wp
  INTEGER :: m1,k
  REAL(wp) :: dmdtm

  DO k=2,m1-1
   
!-- tendency horizontal radius from dyn entrainment
    !rad_t(K) = rad_t(K)   +     (vel_e(k)-vel_p(k)) /pi
    rad_t(K) = rad_t(K)   + ABS((vel_e(k)-vel_p(k)))/pi

!-- entrainment     
    !DMDTM = (2./3.1416)  *     (VEL_E (k) - VEL_P (k)) / RADIUS (k)  
    dmdtm = (20._wp/3.1416_wp)  *  ABS(vel_e (k) - vel_p (k)) / radius (k)  

!-- tendency horizontal speed  from dyn entrainment
    vel_t(k) = vel_t(k)     - dmdtm * ( vel_p (k) - vel_e (k) )

!     if(vel_p (k) - vel_e (k) > 0.) cycle

!-- tendency temperature  from dyn entrainment
    tt(k) = tt(k)           - dmdtm * ( t (k) - te  (k) ) 

!-- tendency water vapor  from dyn entrainment
    qvt(k) = qvt(k)         - dmdtm * ( qv (k) - qvenv (k) )

    qct(k) = qct(k)         - dmdtm * ( qc (k)  )
    qht(k) = qht(k)         - dmdtm * ( qh (k)  )
    qit(k) = qit(k)         - dmdtm * ( qi (k)  )

!-- tendency gas/particle  from dyn entrainment
    sct(k) = sct(k)         - dmdtm * ( sc (k) - sce (k) )

  ENDDO
END SUBROUTINE scl_dyn_entrain

!-------------------------------------------------------------------------------

SUBROUTINE scl_advectc_plumerise(varn,mzp)

  IMPLICIT NONE

  INTEGER :: mzp
  CHARACTER(LEN=*) :: varn
  REAL(wp) :: dtlto2
  INTEGER :: k

!  wp => w
!- Advect  scalars
    dtlto2   = .5_wp * dt
    vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * rho(1)*1.e-3_wp!converte de CGS p/ MKS
    vt3df(1) = .5_wp * (w(1) + wc(1)) * dtlto2 * dzm(1)

    DO k = 2,mzp
      vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5_wp * (rho(k) + rho(k+1))*1.e-3_wp
      vt3df(k) =  (w(k) + wc(k)) * dtlto2 *.5_wp *  dzm(k)
    ENDDO
 
!  DO k = 1,mzp-1
    DO k = 1,mzp
      vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
      vctr2(k) = (zm(k)   - zt(k)) * dzm(k)
      vt3dk(k) = dzt(k) /(rho(k)*1.e-3_wp)
    ENDDO

!      scalarp => scalar_tab(n,ngrid)%var_p
!      scalart => scalar_tab(n,ngrid)%var_t

!- temp advection tendency (tt)
    scr1=T
    CALL fa_zc_plumerise(mzp,t,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)

    CALL advtndc_plumerise(mzp,t,scr1(1),tt)

!- water vapor advection tendency (qvt)
    scr1=qv
    CALL fa_zc_plumerise(mzp,qv,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk (1),vctr1,vctr2)

    CALL advtndc_plumerise(mzp,qv,scr1(1),qvt)

!- liquid advection tendency (qct)
    scr1=qc
    CALL fa_zc_plumerise(mzp,qc,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)

    CALL advtndc_plumerise(mzp,qc,scr1(1),qct)

!- ice advection tendency (qit)
    scr1=qi
    CALL fa_zc_plumerise(mzp,qi,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)

    CALL advtndc_plumerise(mzp,qi,scr1(1),qit)

!- hail/rain advection tendency (qht)
!   IF(ak1 > 0. .or. ak2 > 0.) then

    scr1=qh
    CALL fa_zc_plumerise(mzp,qh,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)

    CALL advtndc_plumerise(mzp,qh,scr1(1),qht)
!   ENDIF


!- horizontal wind advection tendency (vel_t)

    scr1=vel_p
    CALL fa_zc_plumerise(mzp,vel_p,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)

    CALL advtndc_plumerise(mzp,vel_p,scr1(1),vel_t)

!- vertical radius transport

    scr1=rad_p
    CALL fa_zc_plumerise(mzp,rad_p,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)

    CALL advtndc_plumerise(mzp,rad_p,scr1(1),rad_t)
!   return

!- gas/particle advection tendency (sct)
!    IF(varn == 'SC')return
   scr1=sc
   CALL fa_zc_plumerise(mzp,sc,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)
   
   CALL advtndc_plumerise(mzp,sc,scr1(1),sct)


END SUBROUTINE scl_advectc_plumerise
!-------------------------------------------------------------------------------

SUBROUTINE fa_zc_plumerise(m1,scp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2)

  IMPLICIT NONE

  INTEGER :: m1,k
  REAL(wp) :: dfact
  REAL(wp), DIMENSION(m1) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dk
  REAL(wp), DIMENSION(m1) :: vctr1,vctr2

  dfact = .5_wp

! Compute scalar flux VT3DG
  DO k = 1,m1-1
    vt3dg(k) = vt3dc(k)                   &
             * (vctr1(k) * scr1(k)        &
             +  vctr2(k) * scr1(k+1)      &
             +  vt3df(k) * (scr1(k) - scr1(k+1)))
  ENDDO
      
! Modify fluxes to retain positive-definiteness on scalar quantities.
!    If a flux will remove 1/2 quantity during a timestep,
!    reduce to first order flux. This will remain positive-definite
!    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 IF
!    both fluxes are evacuating the box.

  DO k = 1,m1-1
   IF (vt3dc(k) > 0._wp) THEN
     IF (vt3dg(k) * vt3dk(k)    > dfact * scr1(k)) THEN
       vt3dg(k) = vt3dc(k) * scr1(k)
     ENDIF
   ELSEIF (vt3dc(k) < 0._wp) THEN
     IF (-vt3dg(k) * vt3dk(k+1) > dfact * scr1(k+1)) THEN
       vt3dg(k) = vt3dc(k) * scr1(k+1)
     ENDIF
   ENDIF

  ENDDO

! Compute flux divergence
  DO k = 2,m1-1
      scr1(k) = scr1(k)  &
              + vt3dk(k) * ( vt3dg(k-1) - vt3dg(k) &
              + scp  (k) * ( vt3dc(k)   - vt3dc(k-1)))
  ENDDO

END SUBROUTINE fa_zc_plumerise
!-------------------------------------------------------------------------------

SUBROUTINE advtndc_plumerise(m1,scp,sca,sct)
  IMPLICIT NONE

  INTEGER :: m1,k
  REAL(wp) :: dtli
  REAL(wp), DIMENSION(m1) :: scp,sca,sct

  dtli = 1._wp / dt
  DO k = 2,m1-1
    sct(k) = sct(k) + (sca(k)-scp(k)) * dtli
  ENDDO
END SUBROUTINE advtndc_plumerise
!-------------------------------------------------------------------------------

SUBROUTINE tend0_plumerise

 wt  (1:nm1)  = 0._wp
 tt  (1:nm1)  = 0._wp
qvt  (1:nm1)  = 0._wp
qct  (1:nm1)  = 0._wp
qht  (1:nm1)  = 0._wp
qit  (1:nm1)  = 0._wp
vel_t(1:nm1)  = 0._wp
rad_t(1:nm1)  = 0._wp
sct  (1:nm1)  = 0._wp

END SUBROUTINE tend0_plumerise
!-------------------------------------------------------------------------------
SUBROUTINE visc_w(m1,deltak,kmt)

  IMPLICIT NONE
  INTEGER :: m1,k,deltak,kmt,m2
  REAL(wp) :: dz1t,dz1m,dz2t,dz2m,d2wdz,d2tdz,d2qvdz ,d2qhdz ,d2qcdz ,d2qidz ,d2scdz,  &
    &         d2vel_pdz,d2rad_dz
  INTEGER :: mi,mf,deltam

  mi = 2
  mf = MIN(m1,kmt-1)
  deltam = 1

  DO k=mi,mf,deltam !v2
 !DO k=2,m2-1 !orig
    dz1t   = 0.5_wp*(zt(k+1)-zt(k-1))
    dz2t   = visc (k) / (dz1t * dz1t)  
    dz1m   = 0.5_wp*(zm(k+1)-zm(k-1))
    dz2m   = visc (k) / (dz1m * dz1m)  
    d2wdz  = (w  (k + 1) - 2 * w  (k) + w  (k - 1) ) * dz2m  
    d2tdz  = (t  (k + 1) - 2 * t  (k) + t  (k - 1) ) * dz2t  
    d2qvdz = (qv (k + 1) - 2 * qv (k) + qv (k - 1) ) * dz2t  
    d2qhdz = (qh (k + 1) - 2 * qh (k) + qh (k - 1) ) * dz2t 
    d2qcdz = (qc (k + 1) - 2 * qc (k) + qc (k - 1) ) * dz2t  
    d2qidz = (qi (k + 1) - 2 * qi (k) + qi (k - 1) ) * dz2t  
    d2scdz = (sc (k + 1) - 2 * sc (k) + sc (k - 1) ) * dz2t 
    d2vel_pdz=(vel_p  (k + 1) - 2 * vel_p  (k) + vel_p  (k - 1) ) * dz2t
    d2rad_dz =(rad_p  (k + 1) - 2 * rad_p  (k) + rad_p  (k - 1) ) * dz2t
    
    wt(k) =   wt(k) + d2wdz
    tt(k) =   tt(k) + d2tdz
               
! print*,'v=',k,d2tdz,t  (k + 1) - 2 * t  (k) + t  (k - 1),dz2t        
               
    qvt(k) =  qvt(k) + d2qvdz 
    qct(k) =  qct(k) + d2qcdz
    qht(k) =  qht(k) + d2qhdz 
    qit(k) =  qit(k) + d2qidz     

    vel_t(k) =   vel_t(k) + d2vel_pdz

    rad_t(k) =   rad_t(k) + d2rad_dz

    sct(k) =  sct(k) + d2scdz
    !print*,'w-visc=',k,d2wdz
  ENDDO  

END SUBROUTINE visc_w
!     ****************************************************************

SUBROUTINE update_plumerise(m1,varn)

  IMPLICIT NONE

  INTEGER :: m1,k
  CHARACTER(LEN=*) :: varn
   
  IF(varn == 'W') THEN
    DO k=2,m1-1
      w(k) =  w(k) +  wt(k) * dt  
    ENDDO
    RETURN
  ELSE 
    DO k=2,m1-1
      t(k)  =  t(k) +  tt(k) * dt  
      qv(k) = qv(k) + qvt(k) * dt  
      qc(k) = qc(k) + qct(k) * dt !cloud drops travel with air 
      qh(k) = qh(k) + qht(k) * dt  
      qi(k) = qi(k) + qit(k) * dt 
      qv(k) = MAX(0._wp, qv(k))
      qc(k) = MAX(0._wp, qc(k))
      qh(k) = MAX(0._wp, qh(k))
      qi(k) = MAX(0._wp, qi(k))
      vel_p(k) =  vel_p(k) + vel_t(k) * dt  
      rad_p(k) =  rad_p(k) + rad_t(k) * dt  
      sc(k)    =  sc(k)    + sct(k)  * dt 
    ENDDO
  ENDIF

END SUBROUTINE update_plumerise
!-------------------------------------------------------------------------------

SUBROUTINE fallpart(m1)

  IMPLICIT NONE

  INTEGER :: m1,k
  REAL(wp) :: vtc, dfhz,dfiz,dz1
  !srf==================================
  !   verificar se o gradiente esta correto  
  !srf==================================
  !     XNO=1.E7  [m**-4] median volume diameter raindrop,Kessler
  !     VC = 38.3/(XNO**.125), median volume fallspeed eqn., Kessler
  !     for ice, see (OT18), use F0=0.75 per argument there. rho*q
  !     values are in g/m**3, velocities in m/s
  REAL(wp), PARAMETER :: vconst = 5.107387_wp, eps = 0.622_wp, f0 = 0.75_wp
  REAL(wp), PARAMETER :: g = 9.81_wp, cp = 1004._wp

  DO k=2,m1-1
    vtc = vconst * rho (k) **.125_wp   ! median volume fallspeed (ktable4)                              
!  hydrometeor assembly velocity calculations (k table4)
!  vth(k)=-vtc*qh(k)**.125  !median volume fallspeed, water            
    vth (k) = - 4._wp      !small variation with qh
    vhrel = w (k) + vth (k)  !relative to surrounding cloud
!  rain ventilation coefficient for evaporation
    cvh(k) = 1.6_wp + 0.57e-3_wp * (ABS (vhrel) ) **1.5_wp
!  vti(k)=-vtc*f0*qi(k)**.125    !median volume fallspeed,ice             
    vti (k) = - 3._wp                !small variation with qi
    virel = w (k) + vti (k)       !relative to surrounding cloud
!  ice ventilation coefficient for sublimation
    cvi(k) = 1.6_wp + 0.57e-3_wp * (ABS (virel) ) **1.5_wp / f0  

    IF (vhrel>=0.0_wp) THEN  
      dfhz=qh(k)*(rho(k  )*vth(k  )-rho(k-1)*vth(k-1))/rho(k-1)
    ELSE  
      dfhz=qh(k)*(rho(k+1)*vth(k+1)-rho(k  )*vth(k  ))/rho(k)
    ENDIF  

    IF (virel>=0.0_wp) THEN  
      dfiz=qi(k)*(rho(k  )*vti(k  )-rho(k-1)*vti(k-1))/rho(k-1)
    ELSE  
      dfiz=qi(k)*(rho(k+1)*vti(k+1)-rho(k  )*vti(k  ))/rho(k)
    ENDIF
   
    dz1=zm(k)-zm(k-1)
   
    qht(k) = qht(k) - dfhz / dz1 !hydrometeors don't
          
    qit(k) = qit(k) - dfiz / dz1  !nor does ice? hail, what about
  ENDDO

END SUBROUTINE fallpart
! *********************************************************************

SUBROUTINE waterbal  
                                        
  IMPLICIT NONE

  IF (qc (l) <=1.0e-10_wp) qc (l) = 0._wp  !DEFEAT UNDERFLOW PROBLEM
  IF (qh (l) <=1.0e-10_wp) qh (l) = 0._wp  
  IF (qi (l) <=1.0e-10_wp) qi (l) = 0._wp  

  CALL evaporate    !vapor to cloud,cloud to vapor  
                             
  CALL sublimate    !vapor to ice  
                            
  CALL glaciate     !rain to ice 
                         
  CALL melt         !ice to rain
       
!IF(ak1 > 0. .or. ak2 > 0.) &
  CALL convert () !(auto)conversion and accretion 
!CALL convert2 () !(auto)conversion and accretion 
END SUBROUTINE waterbal
! *********************************************************************

SUBROUTINE evaporate  
!- evaporates cloud,rain and ice to saturation
  
  IMPLICIT NONE
  !     XNO=10.0E06
  !     HERC = 1.93*1.E-6*XN035        !evaporation constant
  REAL(wp), PARAMETER :: herc = 5.44e-4_wp, cp = 1.004_wp, heatcond = 2.5E3_wp  
  REAL(wp), PARAMETER :: heatsubl = 2834._wp, tmelt = 273._wp, tfreeze = 269.3_wp
  REAL(wp), PARAMETER :: frc = heatcond / cp, src = heatsubl / cp
  REAL(wp) :: evhdt, evidt, evrate, evap, sd, quant, dividend, divisor, devidt

  sd = qsat (l) - qv (l)  !vapor deficit
  IF (sd==0.0_wp)  RETURN  
  !IF (ABS(SD)<1.e-7)  RETURN  

  evhdt = 0._wp
  evidt = 0._wp  
  !evrate =0.; evap=0.; sd=0.0; quant=0.0; dividend=0.0; divisor=0.0; devidt=0.0                            
  evrate = ABS(wbar * dqsdz)   ! evaporation rate (Kessler 8.32)
  evap = evrate * dt   ! what we can get in dt
                                  
  IF (sd<=0.0_wp) THEN  ! condense. sd is negative 
    IF (evap>=ABS (sd) ) THEN    !we get it all                               
      qc (l) = qc  (l) - sd  !deficit,remember?
      qv (l) = qsat(l)       !set the vapor to saturation  
      t  (l) = t   (l) - sd * frc  !heat gained through condensation per gram of dry air
      RETURN  
    ELSE  
      qc (l) = qc (l) + evap         !get what we can in dt 
      qv (l) = qv (l) - evap         !remove it from the vapor
      t  (l) = t  (l) + evap * frc   !get some heat
      RETURN  
    ENDIF  
  ELSE                               !sd is positive, need some water
! not saturated. saturate if possible. use everything in order
! cloud, rain, ice. sd is positive                                      
    IF (evap<=qc (l) ) THEN        !enough cloud to last dt  
      IF (sd<=evap) THEN          !enough time to saturate
        qc (l) = qc (l) - sd       !remove cloud                                          
        qv (l) = qsat (l)          !saturate
        t (l) = t (l) - sd * frc   !cool the parcel                                          
        RETURN  !done                    
      ELSE   !not enough time                                 
        sd = sd-evap               !use what there is
        qv (l) = qv (l) + evap     !add vapor
        t (l) = t (l) - evap * frc !lose heat
        qc (l) = qc (l) - evap     !lose cloud
                                   !go on to rain.                                      
      ENDIF     
    ELSE                !not enough cloud to last dt     
      IF (sd<=qc (l) ) THEN   !but there is enough to sat                                      
        qv (l) = qsat (l)  !use it
        qc (l) = qc (l) - sd  
        t  (l) = t (l) - sd * frc  
        RETURN                             
      ELSE            !not enough to sat
        sd = sd-qc (l)  
        qv (l) = qv (l) + qc (l)  
        t  (l) = t (l) - qc (l) * frc         
        qc (l) = 0.0_wp  !all gone                                  
      ENDIF       !on to rain                          
    ENDIF          !finished with cloud
!  but still not saturated, so try to use some rain
!  this is tricky, because we only have time dt to evaporate. If there
!  is enough rain, we can evaporate it for dt. ice can also sublimate
!  at the same time. there is a compromise here.....use rain first, then
!  ice. saturation may not be possible in one dt time.
!  rain evaporation rate (W12),(OT25),(K Table 4). evaporate rain first
!  sd is still positive or we wouldn't be here.
  
    IF (qh (l) >1.E-10_wp) THEN                                  
!srf-25082005
!  quant = (qc (l) + qv (l) - qsat (l) ) * rho (l)   !g/m**3
      quant = ( qsat (l)-qc (l) - qv (l)  ) * rho (l)   !g/m**3
      evhdt = (dt * herc * (quant) * (qh (l) * rho (l) ) **.65_wp) / rho (l)
!             rain evaporation in time DT                                   
      IF (evhdt<=qh (l) ) THEN           !enough rain to last DT
        IF (sd<=evhdt) THEN         !enough time to saturate      
          qh (l) = qh (l) - sd    !remove rain      
          qv (l) = qsat (l)       !saturate     
          t (l) = t (l) - sd * frc    !cool the parcel          
          !if(mintime>40) print*,'1',l,t(l)-273.15,qv(l)*1000,qh(l)*1000    
          RETURN              !DOne                      
        ELSE                               !not enough time
          sd = sd-evhdt        !use what there is
          qv (l) = qv (l) + evhdt      !add vapor
          t (l) = t (l) - evhdt * frc      !lose heat
          qh (l) = qh (l) - evhdt      !lose rain
          !if(mintime>40.and. l<40) print*,'2',l,t(l)-273.15,qv(l)*1000,qh(l)*1000
          !if(mintime>40.and. l<40) print*,'3',l,evhdt,quant
          !if(mintime>40.and. l<40) print*,'4',l,qc (l)*1000. , qv (l)*1000. , qsat (l)*1000.                                             
        ENDIF                   !go on to ice.
                                  
      ELSE  !not enough rain to last dt
        IF (sd<=qh (l) ) THEN           !but there is enough to sat
          qv (l) = qsat (l)               !use it
          qh (l) = qh (l) - sd  
          t (l) = t (l) - sd * frc  
          RETURN                            
        ELSE                              !not enough to sat
          sd = sd-qh (l)  
          qv (l) = qv (l) + qh (l)  
          t (l) = t (l) - qh (l) * frc    
          qh (l) = 0.0_wp                    !all gone                                          
        ENDIF                             !on to ice                                       
      ENDIF                                !finished with rain
    ENDIF
!  now for ice
!  equation from (OT); correction factors for units applied

!   33    continue

    IF (qi (l) <=1.E-10_wp) RETURN            !no ice there

    dividend = ( (1.e6_wp / rho (l) ) **0.475_wp) * (sd / qsat (l) &
            &  - 1) * (qi (l) **0.525_wp) * 1.13_wp
    divisor = 7.e5_wp + 4.1e6_wp / (10._wp * est (l) )  
                                                 
    devidt = - cvi(l) * dividend / divisor   !rate of change
                                                  
    evidt = devidt * dt                      !what we could get

! logic here is identical to rain. could get fancy and make SUBROUTINE
! but duplication of code is easier. God bless the screen editor.
                                         
    IF (evidt<=qi (l) ) THEN               !enough ice to last dt                                    
      IF (sd<=evidt) THEN                  !enough time to saturate
        qi (l) = qi (l) - sd                 !remove ice
        qv (l) = qsat (l)                    !saturate
        t  (l) = t (l) - sd * src            !cool the parcel
      RETURN                                 !done                                          
      ELSE                                   !not enough time                                
        sd = sd-evidt                        !use what there is
        qv (l) = qv (l) + evidt              !add vapor
        t  (l) =  t (l) - evidt * src        !lose heat
        qi (l) = qi (l) - evidt              !lose ice                                    
      ENDIF                                  !go on,unsatisfied                                    
    ELSE                                     !not enough ice to last dt         
      IF (sd<=qi (l) ) THEN                !but there is enough to sat     
        qv (l) = qsat (l)                    !use it
        qi (l) = qi   (l) - sd  
        t  (l) =  t   (l) - sd * src  
        RETURN  
      ELSE                                   !not enough to sat
        sd = sd-qi (l)  
        qv (l) = qv (l) + qi (l)  
        t  (l) = t (l) - qi (l) * src             
        qi (l) = 0.0_wp                         !all gone
      ENDIF                                  !on to better things
                                             !finished with ice
    ENDIF                              
  ENDIF                                      !finished with the sd decision

END SUBROUTINE evaporate

! *********************************************************************
SUBROUTINE convert ()  

  IMPLICIT NONE

  !- ACCRETION AND AUTOCONVERSION
  REAL(wp),      PARAMETER ::  ak1 = 0.001_wp    !conversion rate constant
  REAL(wp),      PARAMETER ::  ak2 = 0.0052_wp   !collection (accretion) rate
  REAL(wp),      PARAMETER ::  th  = 0.5_wp      !Kessler threshold
  INTEGER,   PARAMETER ::iconv = 1        !Kessler conversion                                     
  !REAL(wp), PARAMETER :: ANBASE =  50.!*1.e+6 !Berry-number at cloud base #/m^3(maritime)
  REAL(wp), PARAMETER :: anbase =100000._wp!*1.e+6 !Berry-number at cloud base #/m^3(continental)
  ! na formulacao abaixo use o valor em #/cm^3  
  !REAL(wp), PARAMETER :: bdisp = 0.366       !Berry--size dispersion (maritime)
  REAL(wp), PARAMETER :: bdisp = 0.146_wp       !Berry--size dispersion (continental)
  REAL(wp), PARAMETER :: tfreeze = 269.3_wp  !ice formation temperature
  REAL(wp) :: accrete, con, q, h, bc1, bc2, total

  !     selection rules
  IF (t (l)  <= tfreeze) RETURN  !process not allowed above ice
  IF (qc (l) == 0._wp     ) RETURN  
  !srf IF (qc (l) < 1.e-3     ) RETURN  

  accrete = 0._wp  
  con = 0._wp  
  q = rho (l) * qc (l)  
  h = rho (l) * qh (l)  
              
  IF (qh (l) > 0._wp     ) accrete = ak2 * q * (h**.875_wp)  !accretion, Kessler

  IF (iconv/=0) THEN   !select Berry or Kessler
    !old   bc1 = 120.  
    !old   bc2 = .0266 * anbase * 60.  
    !old   con = bdisp * q * q * q / (bc1 * q * bdisp + bc2)    
    con = q*q*q*bdisp/(60._wp*(5._wp*q*bdisp+0.0366_wp*anbase))
  ELSE  
  !   con = ak1 * (q - th)   !kessler autoconversion rate   
  !   if (con<0.0) con = 0.0   !havent reached threshold
    con = max(0._wp,ak1 * (q - th)) ! versao otimizada
  ENDIF  

  total = (con + accrete) * dt / rho (l)  

  IF (total<qc (l) ) THEN  
    qc (l) = qc (l) - total  
    qh (l) = qh (l) + total    !no phase change involved
    RETURN  
  ELSE       
    qh (l) = qh (l) + qc (l)    !uses all there is
    qc (l) = 0.0_wp  
  ENDIF  


END SUBROUTINE convert
!**********************************************************************
SUBROUTINE convert2 ()  
  
  IMPLICIT NONE

  LOGICAL,PARAMETER ::  aerosol=.true.
  REAL(wp), PARAMETER :: tnull=273.16_wp, lat=2.5008E6_wp,  & 
    &                    epsi=0.622_wp ,db=1._wp ,nb=1500._wp !alpha=0.2 
  REAL(wp) :: ka,keins,kzwei,kdrei,vt 
  REAL(wp) :: a,b,c,d, con,accrete,total   
  REAL(wp) :: y(6),roh
        
  a=0.
  b=0.
  y(1) = t(l)
  y(4) = w(l)
  y(2) = qc(l)
  y(3) = qh(l)
  y(5) = radius(l)
  roh =  rho(l)*1.e-3_wp ! dens (mks) ??

  ! autoconversion
  ka = 0.0005_wp 
  IF( y(1) < 258.15_wp )THEN
  !   keins=0.00075
    keins=0.0009_wp 
    kzwei=0.0052_wp
    kdrei=15.39_wp
  ELSE
    keins=0.0015_wp
    kzwei=0.00696_wp
    kdrei=11.58_wp
  ENDIF
        
  !   roh=pe/rd/te
  vt=-kdrei* (y(3)/roh)**0.125_wp

  IF (y(4)>0.0_wp ) THEN
    IF (aerosol) THEN
      a = 1/y(4)  *  y(2)*y(2)*1000._wp/( 60._wp *( 5._wp + 0.0366_wp*nb/(y(2)*1000._wp*db) )  )
    ELSE
      IF (y(2)>(ka*roh)) THEN
        a = keins/y(4) *(y(2) - ka*roh )
      ENDIF
    ENDIF
  ELSE
    a = 0.0_wp
  ENDIF

  ! accretion
  IF(y(4)>0.0_wp) THEN
    b = kzwei/(y(4) - vt) * MAX(0.0_wp,y(2)) *   &
        MAX(0.001_wp,roh)**(-0.875_wp)*(MAX(0.0_wp,y(3)))**(0.875_wp)
  ELSE
     b = 0.0_wp
  ENDIF

      !psatw=610.7*exp( 17.25 *( y(1) - tnull )/( y(1)-36. ) )
      !psate=610.7*exp( 22.33 *( y(1) - tnull )/( y(1)- 2. ) )
      !qsatw=epsi*psatw/( pe-(1.-epsi)*psatw )
      !qsate=epsi*psate/( pe-(1.-epsi)*psate )
      !mu=2.*alpha/y(5)
      !c = mu*( roh*qsatw - roh*qve + y(2) )
      !d = roh*lat*qsatw*epsi/y1/y1/rd *dydx1
      !dydx(2) = - a - b - c - d  ! d rc/dz
      !dydx(3) = a + b            ! d rh/dz
      ! rc=rc+dydx(2)*dz
      ! rh=rh+dydx(3)*dz
   
  con      = a
  accrete  = b
   
  total = (con + accrete) *(1._wp/dzm(l))    ! dt / rho (l)  

  IF (total<qc (l) ) THEN  
    qc (l) = qc (l) - total  
    qh (l) = qh (l) + total    !no phase change involved
    RETURN  
  ELSE               
    qh (l) = qh (l) + qc (l)    !uses all there is
    qc (l) = 0.0_wp
  ENDIF  

END SUBROUTINE convert2
! ice - effect on temperature
!      TTD = 0.0 
!      TTE = 0.0  
!       CALL ICE(QSATW,QSATE,Y(1),Y(2),Y(3), &
!               TTA,TTB,TTC,DZ,ROH,D,C,TTD,TTE)
!       DYDX(1) = DYDX(1) + TTD  + TTE ! DT/DZ on Temp

!**********************************************************************

SUBROUTINE sublimate  

  ! ********************* VAPOR TO ICE (USE EQUATION OT22)***************
  IMPLICIT NONE

  REAL(wp), PARAMETER :: eps = 0.622_wp, heatfus = 334._wp, heatsubl = 2834._wp, cp = 1.004_wp
  REAL(wp), PARAMETER :: src = heatsubl / cp, frc = heatfus / cp, tmelt = 273.3_wp
  REAL(wp), PARAMETER :: tfreeze = 269.3_wp
  REAL(wp) :: dtsubh, dividend,divisor, subl

  dtsubh = 0._wp  
  !selection criteria for sublimation
  IF (t (l)  > tfreeze  ) RETURN  
  IF (qv (l) <= qsat (l) ) RETURN  
  !     from (OT); correction factors for units applied
  dividend = ( (1.e6_wp / rho (l) ) **0.475_wp) * (qv (l) / qsat (l) &
              - 1._wp) * (qi (l) **0.525_wp) * 1.13_wp
  divisor = 7.e5_wp + 4.1e6_wp / (10._wp * est (l) )                                         
  dtsubh = ABS (dividend / divisor)   !sublimation rate
  subl = dtsubh * dt                  !and amount possible

  !     again check the possibilities

  IF (subl<qv (l) ) THEN  
    qv (l) = qv (l) - subl             !lose vapor
    qi (l) = qi (l) + subl             !gain ice
    t  (l) = t  (l) + subl * src       !energy change, warms air
    RETURN  
  ELSE                                      
    qi (l) = qv (l)                    !use what there is
    t  (l) = t (l) + qv (l) * src      !warm the air
    qv (l) = 0.0_wp  
  ENDIF  

END SUBROUTINE sublimate
! *********************************************************************

SUBROUTINE glaciate  

! *********************** CONVERSION OF RAIN TO ICE *******************
!     uses equation OT 16, simplest. correction from W not applied, but
!     vapor pressure differences are supplied.
  IMPLICIT NONE

  REAL(wp), PARAMETER :: heatfus = 334._wp, cp = 1.004_wp, eps = 0.622_wp, heatsubl = 2834._wp
  REAL(wp), PARAMETER :: frc = heatfus / cp, frs = heatsubl / cp, tfreeze = 269.3_wp
  REAL(wp), PARAMETER :: glconst = 0.025_wp   !glaciation time constant, 1/sec
  REAL(wp) :: dfrzh
                                      
  dfrzh = 0._wp    !rate of mass gain in ice
  !selection rules for glaciation
  IF (qh (l) <= 0._wp       ) RETURN  
  IF (qv (l) < qsat (l) ) RETURN                                        
  IF (t  (l) > tfreeze  ) RETURN  
  !      NT=TMELT-T(L)
  !      IF (NT>50) NT=50                                  
  dfrzh = dt * glconst * qh (l)    ! from OT(16)
  IF (dfrzh<qh (l) ) THEN  
    qi (l) = qi (l) + dfrzh  
    qh (l) = qh (l) - dfrzh  
    t  (l) = t  (l) + frc * dfrzh  !warms air   
    RETURN  
  ELSE  
    qi (l) = qi (l) + qh (l)  
    t  (l) = t  (l) + frc * qh (l)  
    qh (l) = 0.0_wp
  ENDIF  

END SUBROUTINE glaciate
! *********************************************************************

SUBROUTINE melt  
! ******************* MAKES WATER OUT OF ICE **************************                                             

  IMPLICIT NONE

  REAL(wp), PARAMETER :: frc = 332.27_wp, tmelt = 273._wp, f0 = 0.75_wp   !ice velocity factor
  REAL(wp) :: dtmelt                                  
  dtmelt = 0._wp   !conversion,ice to rain
  !selection rules
  IF (qi (l) <= 0.0_wp  ) RETURN  
  IF (t (l)  < tmelt) RETURN  
                                                        !OT(23,24)
  dtmelt = dt * (2.27_wp / rho (l) ) * cvi(l) * (t (l) - tmelt) * ( (rho(l)  &
           * qi (l) * 1.e-6_wp) **0.525_wp) * (f0** ( - 0.42_wp) )
                                                        !after mason,1956
  !     check the possibilities
  IF (dtmelt<qi (l) ) THEN  
    qh (l) = qh (l) + dtmelt  
    qi (l) = qi (l) - dtmelt  
    t  (l) = t  (l) - frc * dtmelt     !cools air
    RETURN  
  ELSE  
    qh (l) = qh (l) + qi (l)   !get all there is to get
    t  (l) = t (l) - frc * qi (l)  
    qi (l) = 0.0_wp  
  ENDIF  

END SUBROUTINE melt
!-----------------------------------------------------------------------------
SUBROUTINE zero_plumegen_coms(ucon,zcon,vcon,tmpcon,prcon,rvcon)

  IMPLICIT NONE
  ! JS: passed variables as inout, for init
  REAL(wp),INTENT(INOUT) :: ucon(:),zcon(:),vcon(:),tmpcon(:),prcon(:),rvcon(:)

  w=0._wp;t=0._wp;td=0._wp;theq=0._wp
  qv=0._wp;qc=0._wp;qh=0._wp;qi=0._wp;sc=0._wp;vel_p=0._wp;rad_p=0._wp;rad_t=0._wp
  vth=0._wp;vti=0._wp;rho=0._wp;txs=0._wp
  vt3dc=0._wp;vt3df=0._wp;vt3dk=0._wp;vt3dg=0._wp;scr1=0._wp;vt3dj=0._wp;vt3dn=0._wp;vt3DO=0._wp
  vt3da=0._wp;scr2=0._wp;vt3db=0._wp;vt3dd=0._wp;vt3dl=0._wp;vt3dm=0._wp;vt3di=0._wp;vt3dh=0._wp;vt3de=0._wp
  wc=0._wp;wt=0._wp;tt=0._wp;qvt=0._wp;qct=0._wp;qht=0._wp;qit=0._wp;sct=0._wp;vel_t=0.
  est=0._wp;qsat=0._wp;qpas=0._wp;qtotal=0._wp
  dzm=0._wp;dzt=0._wp;zm=0._wp;zt=0._wp;vctr1=0._wp;vctr2=0._wp
  the=0._wp;thee=0._wp;rhe=0._wp;sce=0._wp
  ucon=0._wp;vcon=0._wp;wcon=0._wp;thtcon =0._wp;rvcon=0._wp;tmpcon=0._wp;prcon=0._wp 
  zcon=0._wp;zzcon=0._wp;scon=0._wp 
  dz=0._wp;dqsdz=0._wp;visc=0._wp;viscosity=0._wp;tstpf=0._wp
  advw=0._wp;advt=0._wp;advv=0._wp;advc=0._wp;advh=0._wp;advi=0._wp;cvh=0._wp;cvi=0._wp;adiabat=0._wp
  wbar=0._wp;alast=0._wp;vhrel=0._wp;virel=0._wp  
  zsurf=0._wp;zbase=0._wp;ztop=0._wp;area=0._wp;rsurf=0._wp;alpha=0._wp;radius=0._wp;heating=0._wp
  fmoist=0._wp;bload=0._wp;dt=0._wp;time=0._wp;tdur=0._wp
  ztop_=0._wp
  n=0;nm1=0;l=0;lbase=0;mintime=0;mdur=0;maxtime=0
  dwdt_entr=0._wp
  kmt = 0
  vel_e(:)=0._wp;pe(:)=0._wp;te(:)=0._wp;qvenv(:)=0._wp;upe(:)=0._wp;vpe(:)=0._wp

END SUBROUTINE zero_plumegen_coms
!     ******************************************************************

SUBROUTINE thetae(p,t,rv,the,tdd)

  IMPLICIT NONE

  REAL(wp) :: p,t,rv,the,tdd
  REAL(wp), PARAMETER :: cp=1004._wp,g=9.8_wp,r=287._wp,alvl=2.35e6_wp,cpg=cp/g
  REAL(wp) :: pit,tupo,ttd,dz,tupn,tmn
  INTEGER :: itter

  pit=p
  tupo=t
  ttd=ftd(p,rv)
  tdd=ttd
  dz=cpg*(t-ttd)
  !print*,'t-ttd:',t-ttd
  IF(dz>0._wp) THEN
    DO itter=1,50
      tupn=t-g/cp*dz
      tmn=(tupn+t)*.5_wp*(1._wp+.61_wp*rv)
!print*,'tmn:',tmn, 'tupn:', tupn
      pit=p*EXP(-g*dz/(r*tmn))
      IF(ABS(tupn-tupo)<0.001_wp) EXIT !GOTO 20
      ttd=ftd(pit,rv)
      tupo=tupn
      dz=dz+cpg*(tupn-ttd)
    ENDDO
    IF (itter==50) print *, 'stop: problems with thetae calculation'
  ENDIF
!  20 continue
  !print*,'pit:',pit,'cp*tupo:',cp*tupo
  the=tupo*(1e5_wp/pit)**.286_wp*EXP(alvl*rv/(cp*tupo))

END SUBROUTINE thetae
!     ******************************************************************

REAL(wp) FUNCTION ftd(p,rs)

  IMPLICIT NONE
  REAL(wp) :: rr,rs,es,esln,p

  rr=rs+1.0e-8_wp
  es=p*rr/(.622_wp+rr)
  esln=LOG(es)
  ftd=(35.86_wp*esln-4947.2325_wp)/(esln-23.6837_wp)

END FUNCTION ftd
!     ******************************************************************
 
SUBROUTINE htint (nzz1, vctra, eleva, nzz2, vctrb, elevb)

  IMPLICIT NONE
  INTEGER, INTENT(IN ) :: nzz1
  INTEGER, INTENT(IN ) :: nzz2
  REAL(wp),    INTENT(IN ) :: vctra(nzz1)
  REAL(wp),    INTENT(OUT) :: vctrb(nzz2)
  REAL(wp),    INTENT(IN ) :: eleva(nzz1)
  REAL(wp),    INTENT(IN ) :: elevb(nzz2)
  INTEGER :: l
  INTEGER :: k
  INTEGER :: kk
  REAL(wp)    :: wt

  l=1

  DO k=1,nzz2
    DO
      IF ( (elevb(k) <  eleva(1)) .OR. &
             ((elevb(k) >= eleva(l)) .AND. (elevb(k) <= eleva(l+1))) ) THEN
        wt       = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
        vctrb(k) = vctra(l)+(vctra(l+1)-vctra(l))*wt
        EXIT
      ELSE IF ( elevb(k) >  eleva(nzz1))  THEN
        wt       = (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
        vctrb(k) = vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
      EXIT
      END IF

      l=l+1
      IF(l == nzz1) THEN
        PRINT *,'htint:nzz1',nzz1
        DO kk=1,l
          PRINT*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
        END DO
        STOP 'htint'
      END IF
    END DO
  END DO
  
END SUBROUTINE htint
!     ******************************************************************
!     ******************************************************************
REAL(wp) FUNCTION  esat(t)

  IMPLICIT NONE
  REAL(wp) :: t
  !     esat(millibars),t(kelvin)
  REAL(wp), PARAMETER :: abz=273.15_wp
  REAL(wp) :: tc

  tc=t-abz
  esat=6.1078_wp*EXP((17.2693882_wp*tc)/(tc+237.3_wp))
  RETURN

END FUNCTION  esat

!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_plumerise

