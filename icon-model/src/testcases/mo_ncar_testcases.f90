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

MODULE mo_ncar_testcases

  !=======================================================================
  !
  !  Functions for setting up initial conditions for the Jablonowski-Williamson test case.
  !
  !  Given longitude (radians), latitude (radians), eta (pressure) and rotation_angle (degrees)
  !  the functions will return temperature, surface geopotential, zonal and meridional wind
  !  components, respectively.
  !
  !  lperturb=.FALSE. result in initial conditions for the steady-state test case.
  !  lperturb=.TRUE.  result in initial conditions for the baroclinic wave test case.
  !
  !     T   : FUNCTION temperature         (lon,lat,eta,rotation_angle)
  !     PHIS: FUNCTION surface_geopotential(lon,lat,rotation_angle)
  !     U   : FUNCTION u_wind              (lon,lat,eta,lperturb,rotation_angle)
  !     V   : FUNCTION v_wind              (lon,lat,eta,lperturb,rotation_angle)
  !     PS  : set to the constant p0
  !
  !  The non-rotated (rotation_angle=0) version of the test cases is described in:
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A baroclinic instability
  !                 test case for atmospheric model dynamical cores.
  !                 Quart. J. Roy. Meteor. Soc., 132, 2943-2975.
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A Baroclinic Wave Test Case
  !                 for Dynamical Cores of General Circulation Models: Model Intercomparisons,
  !                 NCAR Technical Note, NCAR/TN-469+STR, 89 pp.
  !
  !  The rotated version simply rotates the initial conditions so that the spherical coordinate
  !  poles do not conicide with the earth's rotation axis. Thereby the Coriolis parameter is
  !  a function of latitude and longitude:
  !
  !      f = 2*Omega*(-cos(lon)*cos(lat)*sin(rotation_angle)+sin(lat)*cos(rotation_angle))
  !
  !  where Omega = 7.292 x 10E-5/s and rotation_angle is the angle between the flow direction
  !  and equator.
  !
  !  Author: Peter Hjort Lauritzen (NCAR, pel@ucar.edu)
  !          Christiane Jablonowski (University of Michigan, cjablono@umich.edu)
  !
  !=======================================================================

  USE mo_kind,           ONLY: wp
  USE mo_math_constants, ONLY: pi, pi_2, deg2rad
  USE mo_model_domain,   ONLY: t_patch
  USE mo_grid_config,    ONLY: grid_sphere_radius, grid_angular_velocity
  USE mo_exception,      ONLY: finish
  
  IMPLICIT NONE

!=======================================================================
!  physical constants
!=======================================================================

  PUBLIC

  REAL(wp), PARAMETER ::                       &
       Rd         = 287.04_wp,                 & ! gas constant J/(K kg)
       cp         = 1004.64_wp,                & ! specific heat at constant pressure J/(K kg)
       kappa      = Rd/cp,                     & ! kappa = 2/7
       g          = 9.80616_wp                   ! gravitational acceleration (m/s^2)

!-----------------------------------------------------------------------
! steady-state and baroclinic wave tuning parameter
!-----------------------------------------------------------------------
  REAL(wp), PARAMETER ::                             &
       eta_tropo  = 0.2_wp     ,                     & ! tropopause level
       T0         = 288._wp    ,                     & ! horizontal mean T at surface
       eta0       = 0.252_wp   ,                     & ! center of jets (hybrid)
       u0         = 35._wp     ,                     & ! 35 m/s  
       !
       radius                 = 10._wp,             & ! reciprocal radius of the
                                                      ! perturbation without 'a'
       perturbation_amplitude =  1._wp,             & ! amplitude of u perturbation 1 m/s
       perturbation_longitude = 20._wp,             & ! longitudinal position, 20E
       perturbation_latitude  = 40._wp,             & ! latitudinal position, 40N
       eta_sfc                = 1._wp,              & ! hybrid value at surface
       delta_T                = 480000._wp,         & ! in K, for T mean calculation
       gamma                  = 0.005_wp,           & ! lapse rate
       !
       perturbation_latitude_tracer = 55._wp,        &
       exponent               = Rd*gamma/g

   REAL(wp) :: a       ! Domain radius in m
   REAL(wp) :: omega   ! angular velocity 1/s
   REAL(wp) :: a_omega ! periphery velocity =a*omega

CONTAINS

  !-------------------------------------------------------------------
  !>
  SUBROUTINE init_ncar_testcases_domain()
    a = grid_sphere_radius
    omega = grid_angular_velocity
    a_omega = a*omega
  END SUBROUTINE init_ncar_testcases_domain
  !-------------------------------------------------------------------

  
!********************************************************************
!
! Temperature (equation (6) in Jablonowski and Williamson, 2006)
!
!********************************************************************
  REAL(wp) FUNCTION temperature(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(wp)             :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8_wp) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,1)
    ENDIF

    temperature  = t_mean(eta) + t_deviation(rot_lat,eta)
  END FUNCTION temperature
  !
  ! Horizontally averaged temperature (equation (4) and (5) in Jablonowski and Williamson (2006))
  !
  REAL(wp) FUNCTION t_mean(eta)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: eta

    IF (eta.GT.(eta_tropo)) THEN
       t_mean = T0*eta**exponent     ! mean temperature at each level (troposphere)
    ELSE
       t_mean = T0*eta**exponent + delta_T*(eta_tropo-eta)**5  ! ...  (stratosphere)
    ENDIF
  END FUNCTION t_mean
  !
  ! Temperature deviation from the horizontal mean
  ! (equation (6) minus horizontally averaged temperature)
  !
  REAL(wp) FUNCTION t_deviation(lat,eta)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: eta, lat
    REAL(wp)             :: factor, phi_vertical, rot_lat

    factor       = eta*pi*u0/Rd
    phi_vertical = (eta - eta0) * 0.5_wp*pi

    rot_lat = lat

    t_deviation = factor * 1.5_wp * SIN(phi_vertical) * (COS(phi_vertical))**0.5_wp               &
         &      * ((-2._wp*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1._wp/3._wp) + 10._wp/63._wp) &
         &      * u0 * (COS(phi_vertical))**1.5_wp + (8._wp/5._wp*(COS(rot_lat))**3               &
         &      * ((SIN(rot_lat))**2 + 2._wp/3._wp) - pi/4._wp)*a_omega*0.5_wp )

  END FUNCTION t_deviation

!**************************************************************************
!
! Surface geopotential (equaiton (7) in Jablonowski and Williamson, 2006)
!
!**************************************************************************
  REAL(wp) FUNCTION surface_geopotential(lon,lat,rotation_angle)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: lon, lat, rotation_angle
    REAL(wp)             :: cos_tmp, rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8_wp) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,1)
    ENDIF

    cos_tmp    = u0 * (COS((eta_sfc-eta0)*pi*0.5_wp))**1.5_wp

    surface_geopotential = ((-2._wp*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1._wp/3._wp) &
         &               + 10._wp/63._wp)*COS_tmp + (8._wp/5._wp*(COS(rot_lat))**3        &
         &               * ((SIN(rot_lat))**2 + 2._wp/3._wp) - pi/4._wp)*a_omega)*COS_tmp

  END FUNCTION surface_geopotential

!********************************************************************
!
! wind components (equation 2 in Jablonowski and Williamson, 2006)
!
!********************************************************************
  REAL(wp) FUNCTION u_wind(lon,lat,eta,lperturb,rotation_angle)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: lon,lat,eta,rotation_angle
    LOGICAL, INTENT(IN)  :: lperturb
    REAL(wp) :: u_lat, phi_vertical, rot_lon, rot_lat, sin_tmp, cos_tmp, r, u_perturb, v_lat
    REAL(wp) :: perturb_lon, perturb_lat, v_tmp

    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    IF (ABS(rotation_angle)<1.0E-8_wp) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,1)
    ENDIF

    phi_vertical = (eta - eta0) *0.5_wp*pi
    u_lat = (COS(phi_vertical))**1.5_wp * 4._wp * u0 * (SIN(rot_lat))**2 * (COS(rot_lat))**2
    u_wind = u_lat

    IF (lperturb) THEN

       sin_tmp = SIN(perturb_lat)*SIN(rot_lat)
       cos_tmp = COS(perturb_lat)*COS(rot_lat)

       r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-perturb_lon) )    ! great circle distance
                                                                 ! without radius 'a'
       u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )
       IF (u_perturb <= 1.e-6_wp) u_perturb = 0._wp
       u_lat     = u_perturb + u_lat                             ! zonal wind
    ENDIF
    IF (ABS(rotation_angle)<1.0E-8_wp) THEN
       u_wind = u_lat
    ELSE
       v_lat = 0.0_wp
       !
       ! rotate wind components
       !
       CALL turnwi(u_lat,v_lat, u_wind,v_tmp,lon,lat,rot_lon,rot_lat, &
            &      0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,-1)
       IF (ABS(u_wind)<1.0E-10_wp) u_wind=0.0_wp
    ENDIF
  END FUNCTION u_wind

  REAL(wp) FUNCTION v_wind(lon,lat,eta,lperturb,rotation_angle)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: lon,lat,eta,rotation_angle
    LOGICAL, INTENT(IN)  :: lperturb
    REAL(wp) :: u_lat, phi_vertical, rot_lon, rot_lat, sin_tmp, cos_tmp, r, u_perturb, v_lat
    REAL(wp) :: perturb_lon, perturb_lat, u_tmp

    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    IF (ABS(rotation_angle)<1.0E-8_wp) THEN
       v_wind = 0.0_wp
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,1)


       phi_vertical = (eta - eta0) *0.5_wp*pi
       u_lat = (COS(phi_vertical))**1.5_wp * 4._wp * u0 * (SIN(rot_lat))**2 * (COS(rot_lat))**2

       IF (lperturb) THEN

          sin_tmp = SIN(perturb_lat)*SIN(rot_lat)
          cos_tmp = COS(perturb_lat)*COS(rot_lat)

          r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-perturb_lon) )    ! great circle distance
                                                                    ! withour radius 'a'
          u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )
          IF (u_perturb <= 1.e-6_wp) u_perturb = 0._wp
          u_lat     = u_perturb + u_lat
       ENDIF

       v_lat = 0.0_wp
       !
       ! pole point velocities are not well-defined
       !
       IF (ABS(pi*0.5_wp-lat)<1.0E-8_wp.OR.ABS(pi*0.5_wp+lat)<1.0E-8_wp) THEN
          v_wind = 0.0_wp
       ELSE
          !
          ! rotate wind components
          !
          CALL turnwi(u_lat,v_lat, u_tmp,v_wind,lon,lat,rot_lon,rot_lat, &
               &      0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,-1)
       ENDIF
    ENDIF
  END FUNCTION v_wind

!******************************************************************************
!
! Subroutines for rotation
!
!******************************************************************************
  SUBROUTINE regrot(pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)
    IMPLICIT NONE
!
!----------------------------------------------------------------------
!
!*    conversion between regular and rotated spherical coordinates.
!*
!*    pxreg     longitudes of the regular coordinates
!*    pyreg     latitudes of the regular coordinates
!*    pxrot     longitudes of the rotated coordinates
!*    pyrot     latitudes of the rotated coordinates
!*              all coordinates given in degrees n (negative for s)
!*              and degrees e (negative values for w)
!*    pxcen     regular longitude of the south pole of the rotated grid
!*    pycen     regular latitude of the south pole of the rotated grid
!*
!*    kcall=-1: find regular as functions of rotated coordinates.
!*    kcall= 1: find rotated as functions of regular coordinates.
!
!-----------------------------------------------------------------------
!
      INTEGER  ::kcall
      REAL(wp) :: pxreg,pyreg,&
                  pxrot,pyrot,&
                  pxcen,pycen
!
!-----------------------------------------------------------------------
!
      REAL(wp) zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
               zsyrot,zcyrot,zcxrot,zsxrot
!
!----------------------------------------------------------------------
!
      zsycen = SIN((pycen+pi_2))
      zcycen = COS((pycen+pi_2))
!
      IF (kcall.EQ.1) THEN
!
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
         zsyrot = MAX(zsyrot,-1._wp)
         zsyrot = MIN(zsyrot,+1._wp)
         !
         pyrot = ASIN(zsyrot)
         !
         zcyrot = COS(pyrot)
         zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
         zcxrot = MAX(zcxrot,-1._wp)
         zcxrot = MIN(zcxrot,+1._wp)
         zsxrot = zcyreg*zsxmxc/zcyrot
         !
         pxrot = ACOS(zcxrot)
         !
         IF (zsxrot<0.0_wp) pxrot = -pxrot
               !
      ELSEIF (kcall.EQ.-1) THEN
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
         zsyreg = MAX(zsyreg,-1._wp)
         zsyreg = MIN(zsyreg,+1._wp)
         !
         pyreg = ASIN(zsyreg)
         !
         zcyreg = COS(pyreg)
         zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
         zcxmxc = MAX(zcxmxc,-1._wp)
         zcxmxc = MIN(zcxmxc,+1._wp)
         zsxmxc = zcyrot*zsxrot/zcyreg
         zxmxc  = ACOS(zcxmxc)
         IF (zsxmxc<0.0_wp) zxmxc = -zxmxc
         !
         pxreg = zxmxc + pxcen
         !
      ELSE
         WRITE(6,'(1x,''invalid kcall in regrot'')')
         CALL finish('','')
      ENDIF
    END SUBROUTINE regrot

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    SUBROUTINE turnwi(puarg,pvarg,pures,pvres,   &
                      pxreg,pyreg,pxrot,pyrot,   &
                      pxcen,pycen,kcall)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!*    turn horizontal velocity components between regular and
!*    rotated spherical coordinates.
!
!*    puarg : input u components
!*    pvarg : input v components
!*    pures : output u components
!*    pvres : output v components
!*    pa    : transformation coefficients
!*    pb    :    -"-
!*    pc    :    -"-
!*    pd    :    -"-
!*    pxreg : regular longitudes
!*    pyreg : regular latitudes
!*    pxrot : rotated longitudes
!*    pyrot : rotated latitudes
!*    kxdim              : dimension in the x (longitude) direction
!*    kydim              : dimension in the y (latitude) direction
!*    kx                 : number of gridpoints in the x direction
!*    ky                 : number of gridpoints in the y direction
!*    pxcen              : regular longitude of the south pole of the
!*                         transformed grid
!*    pycen              : regular latitude of the south pole of the
!*                         transformed grid
!*
!*    kcall < 0          : find wind components in regular coordinates
!*                         from wind components in rotated coordinates
!*    kcall > 0          : find wind components in rotated coordinates
!*                         from wind components in regular coordinates
!*    note that all coordinates are given in degrees n and degrees e.
!*       (negative values for s and w)
!
!-----------------------------------------------------------------------

      INTEGER  kcall
      REAL(wp) puarg,pvarg,    &
               pures,pvres,    &
               pa,   pb,       &
               pc,   pd,       &
               pxreg,pyreg,    &
               pxrot,pyrot
      REAL(wp) pxcen,pycen
!-----------------------------------------------------------------------
      REAL(wp) zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,&
               zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot
!-----------------------------------------------------------------------
      IF (kcall.EQ.1) THEN
         zsyc = SIN(pycen+pi_2)
         zcyc = COS(pycen+pi_2)
         !
         zsxreg = SIN(pxreg)
         zcxreg = COS(pxreg)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         !
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         !
         pa = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
         pb = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - &
              zsxmxc*zsyreg*zcxrot
         pc = zsyc*zsxmxc/zcyrot
         pd = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
         !
         pures = pa*puarg + pb*pvarg
         pvres = pc*puarg + pd*pvarg
      ELSEIF (kcall.EQ.-1) THEN
         zsyc = SIN(pycen+pi_2)
         zcyc = COS(pycen+pi_2)
         !
         zsxreg = SIN(pxreg)
         zcxreg = COS(pxreg)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         !
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         !
         pa = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
         pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot - &
              zcxmxc*zsxrot*zsyrot
         pc =-zsyc*zsxrot/zcyreg
         pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
         !
         pures = pa*puarg + pb*pvarg
         pvres = pc*puarg + pd*pvarg
      ELSE
         WRITE(6,'(1x,''invalid kcall in turnwi'')')
         CALL finish('','')
      ENDIF
    END SUBROUTINE turnwi

!********************************************************************
!
! Tracers
!
!********************************************************************

!-----------------------------------------------------------------------
! Tracer q1 and q2
!-----------------------------------------------------------------------
  REAL(wp) FUNCTION tracer_q1_q2(lon,lat,eta,rotation_angle, eta_c)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: eta, lon, lat, rotation_angle, eta_c
    REAL(wp) :: rot_lon, rot_lat, sin_tmp, cos_tmp, r
    REAL(wp) :: rot_perturb_lon, rot_perturb_lat, tmp

    rot_perturb_lon = perturbation_longitude*deg2rad
    rot_perturb_lat = perturbation_latitude_tracer *deg2rad

    IF (ABS(rotation_angle)<1.0E-8_wp) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,1)
    ENDIF
    sin_tmp = SIN(rot_perturb_lat)*SIN(rot_lat)
    cos_tmp = COS(rot_perturb_lat)*COS(rot_lat)
    r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-rot_perturb_lon) )    ! great circle distance

    tmp = EXP(- ((r*radius)**2 + ((eta-eta_c)/0.1_wp)**2))
    IF (ABS(tmp)<1.0E-8_wp) tmp = 0.0_wp
    tracer_q1_q2 = tmp
  END FUNCTION tracer_q1_q2

!-----------------------------------------------------------------------
! Tracer q3
!-----------------------------------------------------------------------
  REAL(wp) FUNCTION tracer_q3(lon,lat,rotation_angle)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: lon, lat, rotation_angle
    REAL(wp) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8_wp) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q3 = 0.5_wp * ( TANH( 3._wp*ABS(rot_lat)-pi ) + 1._wp)

  END FUNCTION tracer_q3

!-----------------------------------------------------------------------
! Tracer q, absolute value of the relative vorticity of the unperturbed initial state
!           multiplied by 10^5
!-----------------------------------------------------------------------
  REAL(wp) FUNCTION tracer_q(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(wp) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8_wp) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q = ABS(-4._wp * u0/a * (COS((eta-eta0)*pi*0.5_wp))**1.5_wp * SIN(rot_lat) * &
               COS(rot_lat) * (2._wp-5._wp*(SIN(rot_lat))**2)) * 1.e5_wp
    IF (tracer_q < 1.e-9_wp) tracer_q = 0._wp  !  otherwise error in netcdf file

  END FUNCTION tracer_q

!==========================================================================================
! pure 3D advection, time-dependent
!==========================================================================================
  SUBROUTINE init_pure_adv_wind (lon, lat, rotation_angle,  &
    &                            u_wind, v_wind)
  !-----------------------------------------------------------------------
  !     input parameters
  !-----------------------------------------------------------------------
  REAL(wp), INTENT(in)  :: lon,              & ! longitude in radians
                           lat,              & ! latitude in radians
                           rotation_angle      ! alpha in degrees
  !-----------------------------------------------------------------------
  !     output parameters
  !-----------------------------------------------------------------------
  REAL(wp), INTENT(out) :: u_wind,           & ! zonal wind in m/s
                           v_wind              ! meridional wind in m/s
  !-----------------------------------------------------------------------
  !     test case parameters
  !-----------------------------------------------------------------------
  REAL(wp) :: this_u0   ! circumference / 12 days
  !-----------------------------------------------------------------------
  !     local variables
  !-----------------------------------------------------------------------
  REAL(wp) :: alpha, rot_lon, rot_lat, u_lat, v_lat, u_tmp, v_tmp

  this_u0 = (2._wp*pi*a)/(12._wp*86400._wp)  ! circumference / 12 days
  alpha   = rotation_angle*deg2rad
  !-----------------------------------------------------------------------
  !    initialize the wind components
  !-----------------------------------------------------------------------
  IF (ABS(rotation_angle)<1.0E-8_wp) THEN
    rot_lon = lon
    rot_lat = lat
  ELSE
    CALL regrot(lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+alpha,1)
  ENDIF

  u_lat  =  this_u0 *  COS(rot_lat)

  IF (ABS(rotation_angle)<1.0E-8_wp) THEN
    u_wind = u_lat
  ELSE
    v_lat = 0.0_wp
    !
    ! rotate wind components
    !
    CALL turnwi(u_lat,v_lat,u_wind,v_tmp,lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+alpha,-1)
    IF (ABS(u_wind)<1.0E-10_wp) u_wind=0.0_wp
  ENDIF
  IF (ABS(rotation_angle)<1.0E-8_wp) THEN
    v_wind = 0.0_wp
  ELSE
    v_lat = 0.0_wp
    !
    ! pole point velocities are not well-defined
    !
    IF (ABS(pi*0.5_wp-lat)<1.0E-8_wp.OR.ABS(pi*0.5_wp+lat)<1.0E-8_wp) THEN
      v_wind = 0.0_wp
    ELSE
      !
      ! rotate wind components
      !
      CALL turnwi(u_lat,v_lat,u_tmp,v_wind,lon,lat,rot_lon,rot_lat,0.0_wp,-0.5_wp*pi+alpha,-1)
    ENDIF
    !
  ENDIF
  !-----------------------------------------------------------------------
  !     initialization of the vertical velocity:
  !     must be implemented in the dynamical core
  !-----------------------------------------------------------------------

  END SUBROUTINE init_pure_adv_wind

  SUBROUTINE init_pure_adv_tracers (tracer_variant, lon, lat, height,  &
    &                               q4, q5, q6, q7, q8)
  !-----------------------------------------------------------------------
  !     input parameters
  !-----------------------------------------------------------------------
!!$  CHARACTER*5, INTENT(in) :: tracer_variant    ! identifies test variant 'yy', here tracers
  INTEGER, INTENT(in) :: tracer_variant(:)    ! identifies test variant 'yy', here tracers
  ! e.g. 0 : no tracer, set to zero
  !      5 : tracer q5 only
  !     56 : both tracers q5 and q6
  !    456 : tracers q4, q5 and q6
  REAL(wp), INTENT(in)  :: lon,              & ! longitude in radians
    lat,                                     & ! latitude in radians
    height
  !-----------------------------------------------------------------------
  !     output parameters
  !-----------------------------------------------------------------------
  REAL(wp), INTENT(out) :: q4,      & ! tracer q4
   &  q5,                           & ! tracer q5
   &  q6,                           & ! tracer q6
   &  q7,                           & ! tracer q7
   &  q8                              ! tracer q8
  !-----------------------------------------------------------------------
  !     test case parameters
  !-----------------------------------------------------------------------
  REAL(wp), PARAMETER ::                       &
    RR         = 1._wp/3._wp,                  & ! horizontal half width divided by 'a'
    ZZ         = 1000._wp,                     & ! vertical half width
    z0         = 4500._wp,                     & ! center point in z
    lambda0    = 1.5_wp*pi,                    & ! center point in longitudes
    RR_q7      = 1._wp/2._wp,                  & ! half width divided by 'a' for q7
    lambda0_q7 = 1.5_wp*pi-pi,                 & ! center point in longitudes for q7
    phi0       = 0._wp,                        & ! center point in latitudes
    slot       = 1._wp/8._wp                     ! half width of the slot in radians
  !-----------------------------------------------------------------------
  !     local variables
  !-----------------------------------------------------------------------
  REAL(wp) :: sin_tmp, cos_tmp
  REAL(wp) :: d1, d2, r

  !-----------------------------------------------------------------------
  !     Tracer variables
  !-----------------------------------------------------------------------
  q4 = 0._wp   ! default
  q5 = 0._wp   ! default
  q6 = 0._wp   ! default
  q7 = 0._wp   ! default
  q8 = 0._wp   ! default
  !-----------------------------------------------------------------------
  !     tracer q4 (2D field after Williamson et. al. (1992))
  !-----------------------------------------------------------------------
!!$  IF (tracer_variant(1:1) == '4' .OR. tracer_variant(2:2) == '4'  &
!!$    &                            .OR. tracer_variant(3:3) == '4'  &
!!$    &                            .OR. tracer_variant(4:4) == '4'  &
!!$    &                            .OR. tracer_variant(5:5) == '4')    THEN
  IF (ANY(tracer_variant == 4)) THEN
    sin_tmp = SIN(lat) * SIN(phi0)
    cos_tmp = COS(lat) * COS(phi0)
    r  = ACOS (sin_tmp + cos_tmp*COS(lon-lambda0))       ! great circle distance without 'a'
    d1 = MIN( 1._wp, (r/RR) )
    q4 = 0.5_wp  * (1._wp + COS(pi*d1))
!DR    q4 = 0.5_wp * ZZ * (1._wp + COS(pi*d1)) ! Almut: for comparison with Bill
  ENDIF
  !-----------------------------------------------------------------------
  !     tracer q5 (3D-field after Jablonowski et al (2008))
  !-----------------------------------------------------------------------
!!$  IF (tracer_variant(1:1) == '5' .OR. tracer_variant(2:2) == '5'  &
!!$    &                            .OR. tracer_variant(3:3) == '5'  &
!!$    &                            .OR. tracer_variant(4:4) == '5'  &
!!$    &                            .OR. tracer_variant(5:5) == '5') THEN
  IF (ANY(tracer_variant == 5)) THEN
    sin_tmp = SIN(lat) * SIN(phi0)
    cos_tmp = COS(lat) * COS(phi0)
    r  = ACOS (sin_tmp + cos_tmp*COS(lon-lambda0))       ! great circle distance without 'a'
    d1 = MIN( 1._wp, (r/RR)**2 + ((height-z0)/ZZ)**2 )
    q5 = 0.5_wp * (1._wp + COS(pi*d1))
  ENDIF
  !-----------------------------------------------------------------------
  !     tracer q6 (3D slotted cylinder after Jablonowski et al (2008))
  !-----------------------------------------------------------------------
!!$  IF (tracer_variant(1:1) == '6' .OR. tracer_variant(2:2) == '6'  &
!!$    &                            .OR. tracer_variant(3:3) == '6'  &
!!$    &                            .OR. tracer_variant(4:4) == '6'  &
!!$    &                            .OR. tracer_variant(5:5) == '6') THEN
  IF (ANY(tracer_variant == 6)) THEN
    sin_tmp = SIN(lat) * SIN(phi0)
    cos_tmp = COS(lat) * COS(phi0)
    r  = ACOS (sin_tmp + cos_tmp*COS(lon-lambda0))       ! great circle distance without 'a'
    d2 = (r/RR)**2 + ((height-z0)/ZZ)**2
    IF (d2 <= 1._wp) THEN
      q6 = 1._wp
    ELSE
      q6 = 0._wp
    ENDIF
    IF ((height > z0) .AND. ((phi0-slot) < lat .AND. lat < (phi0+slot)) ) q6 = 0._wp   ! slotted
                                                                                       ! ellipse
  ENDIF
  !-----------------------------------------------------------------------
  !     tracer q7 (2D slotted cylinder)
  !     See e.g.: Nair and Lauritzen (2010), JCP
  ! Please note: -pi < lon < pi.
  ! Therefore pi has been added below for the computation of q7 and
  ! substracted from the original lambda0.
  !
  ! Constants according to Lipscomb and Ringler, MWR (2005)
  ! radius RR*   : 1/2
  ! slot width*  : (1/6)*RR
  ! slot length* : (5/3)*RR
  ! *divided by earth radius
  !-----------------------------------------------------------------------
!!$  IF (tracer_variant(1:1) == '7' .OR. tracer_variant(2:2) == '7'  &
!!$    &                            .OR. tracer_variant(3:3) == '7'  &
!!$    &                            .OR. tracer_variant(4:4) == '7'  &
!!$    &                            .OR. tracer_variant(5:5) == '7') THEN
  IF (ANY(tracer_variant == 7)) THEN
    sin_tmp = SIN(lat) * SIN(phi0)
    cos_tmp = COS(lat) * COS(phi0)
    r  = ACOS (sin_tmp + cos_tmp*COS(lon+pi-lambda0_q7)) ! great circle distance without 'a'
    d2 = r/RR_q7

    IF (d2 <= 1._wp .AND. ABS(lon+pi-lambda0_q7) >= RR_q7/6._wp) THEN
      q7 = 1._wp
     ELSE IF (d2 <= 1._wp .AND. ABS(lon+pi-lambda0_q7) < RR_q7/6._wp &
      &                   .AND. (phi0-lat) > (2._wp/3._wp)*RR_q7)THEN
      q7 = 1._wp
    ELSE
      q7 = 0._wp
    ENDIF

  ENDIF
  !-----------------------------------------------------------------------
  ! tracer q8 (2D field after Harris and Lauritzen (2010))
  !
  ! This initial condition has continuous second and third derivatives.
  ! It is thus can be used to check for third order convergence
  !-----------------------------------------------------------------------
!!$  IF (tracer_variant(1:1) == '8' .OR. tracer_variant(2:2) == '8'  &
!!$    &                            .OR. tracer_variant(3:3) == '8'  &
!!$    &                            .OR. tracer_variant(4:4) == '8'  &
!!$    &                            .OR. tracer_variant(5:5) == '8') THEN
  IF (ANY(tracer_variant == 8)) THEN
    sin_tmp = SIN(lat) * SIN(phi0)
    cos_tmp = COS(lat) * COS(phi0)
    r  = ACOS (sin_tmp + cos_tmp*COS(lon-lambda0))       ! great circle distance without 'a'
    d1 = MIN( 1._wp, (r/RR) )
    q8 = 0.25_wp * (1._wp + COS(pi*d1))**2
  ENDIF

  END SUBROUTINE init_pure_adv_tracers

!==========================================================================================
! Rossby_Haurwitz wave, wavenumber 4
!==========================================================================================
  SUBROUTINE Rossby_Haurwitz (lon, lat, pressure,                                &
                              u_wind, v_wind, temperature, surface_geopotential, &
                              surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      REAL(wp), INTENT(in)  :: lon,                  & ! longitude in radians
                               lat,                  & ! latitude in radians
                               pressure                ! pressure at full model level in Pa
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      REAL(wp), INTENT(out) :: u_wind,               & ! zonal wind in m/s
                               v_wind,               & ! meridional wind in m/s
                               temperature,          & ! temperature in K
                               surface_geopotential, & ! surface geopotential in m^2/s^2
                               surface_pressure        ! surface pressure in Pa

!-----------------------------------------------------------------------
!     test case parameters
!-----------------------------------------------------------------------
      REAL(wp),PARAMETER :: u0      = 50._wp,        &   ! reference wind
                            T0      = 288._wp,       &   ! reference temperature
                            n       = 4._wp,         &   ! wavenumber
                            gamma   = 0.0065_wp,     &   ! lapse rate in K/m
                            p_ref   = 95500._wp          ! reference pressure

!-----------------------------------------------------------------------
!     local
!-----------------------------------------------------------------------
      REAL(wp) ::   this_MM,      &   ! parameter M and p_ref=95500.
                    this_KK           ! parameter K
      REAL(wp) :: tmp1, tmp2, tmp3
      REAL(wp) :: cos_lat, sin_lat
      REAL(wp) :: exponent_1, exponent_2
      REAL(wp) :: AA, BB, CC
      REAL(wp) :: phis_perturb

      this_MM      = u0/(n*a)   ! parameter M and p_ref=95500.
      this_KK      = u0/(n*a)   ! parameter K
!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      cos_lat = COS(lat)
      sin_lat = SIN(lat)
      tmp1 = a * this_MM * cos_lat
      tmp2 = a * this_KK * cos_lat**(n-1._wp)*(n*sin_lat**2 - cos_lat**2)
      tmp3 = -a * this_KK * n * cos_lat**(n-1._wp) * sin_lat
      u_wind = tmp1 + tmp2 * COS(n*lon)
      v_wind = tmp3 * SIN(n*lon)
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0._wp
!-----------------------------------------------------------------------
!     initialize surface pressure and temperature
!-----------------------------------------------------------------------
      tmp1       = gamma/(g*T0)
      tmp2       = a*a
      exponent_1 = g/(gamma*Rd)
      exponent_2 = (gamma*Rd)/g

      cos_lat = COS(lat)
      AA = tmp2 * (0.5_wp * this_MM*(2._wp*omega+this_MM) * cos_lat**2 + 0.25_wp * this_KK**2  &
           &    * cos_lat**(2._wp*n) * ( (n+1._wp)*cos_lat**2 + (2._wp*n*n - n - 2._wp)) &
           &    - 0.5_wp*n*n*this_KK**2 * cos_lat**(2._wp*(n-1._wp)))
      BB = tmp2 * (2._wp*(omega+this_MM)*this_KK/((n+1._wp)*(n+2._wp)) * cos_lat**n &
           &    * ( (n*n + 2._wp*n +2._wp) - (n+1._wp)**2 * cos_lat**2 ))
      CC = tmp2 * (0.25_wp * this_KK**2 * cos_lat**(2._wp*n) * ( (n+1._wp)*cos_lat**2 - (n+2._wp)))
      phis_perturb = AA + BB * COS(n*lon) + CC * COS(2._wp*n*lon)
      surface_pressure = p_ref * (1._wp + tmp1*phis_perturb)**exponent_1   ! surface pressure
      temperature      = T0 * (pressure/p_ref)**exponent_2                 ! temperature

  END SUBROUTINE Rossby_Haurwitz

!==========================================================================================
! Mountain induced Rossby wave
!==========================================================================================
  SUBROUTINE mountain_Rossby (lon, lat,                                          &
                              u_wind, v_wind, temperature, surface_geopotential, &
                              surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      REAL(wp), INTENT(in)  :: lon,                  & ! longitude in radians
                               lat                     ! latitude in radians
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      REAL(wp), INTENT(out) :: u_wind,               & ! zonal wind in m/s
                               v_wind,               & ! meridional wind in m/s
                               temperature,          & ! temperature in K
                               surface_geopotential, & ! surface geopotential in m^2/s^2
                               surface_pressure        ! surface pressure in Pa
!-----------------------------------------------------------------------
!     test case parameters
!-----------------------------------------------------------------------
      REAL(wp),PARAMETER :: u0      = 20._wp,        & ! 20 m/s
                            T0      = 288._wp,       & ! temperature
                            N2      = g*g/(cp*T0),   & ! squared Brunt Vaisala frequency N^2
                            h0      = 2000._wp,      & ! amplitude of the mountain, 2km
                            d       = 1500.e3_wp,    & ! half width 1500 km
                            lambda0 = 0.5_wp*pi,     & ! center point in longitudes
                            phi0    = pi/6._wp,      & ! center point in latitudes
                            p_sp    = 93000._wp        ! pressure at the South Pole in Pa
!-----------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------
      REAL(wp) :: sin_tmp, cos_tmp
      REAL(wp) :: tmp1, tmp2, tmp3
      REAL(wp) :: r

!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * COS(lat)
      v_wind = 0._wp
!-----------------------------------------------------------------------
!     initialize T (temperature)
!-----------------------------------------------------------------------
      temperature = T0
!-----------------------------------------------------------------------
!     initialize surface geopotential and surface pressure
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2._wp * g*g * kappa) * (u0/a + 2._wp * omega)
      tmp2 = N2 / (g*g * kappa)
      sin_tmp = SIN(lat) * SIN(phi0)
      cos_tmp = COS(lat) * COS(phi0)
      tmp3 = tmp1*((SIN(lat))**2 - 1._wp)
      r = a * ACOS (sin_tmp + cos_tmp*COS(lon-lambda0))   ! great circle distance with 'a'
      surface_geopotential = g*h0 * EXP(-(r/d)**2)        ! Gaussian profile of the mountain
      surface_pressure     = p_sp * EXP( -tmp3 - tmp2*surface_geopotential)

  END SUBROUTINE mountain_Rossby

!==========================================================================================
! gravity waves
!==========================================================================================
  SUBROUTINE gravity_wave (choice, lon, lat, height,                          &
                           u_wind, v_wind, temperature, surface_geopotential, &
                           surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      INTEGER, INTENT(in)   :: choice                  ! identifies test variant 'x', e.g:
                                                       ! 0 : no Coriolis, N=0.01 1/s, u0=0  m/s
                                                       ! 1 : no Coriolis, N=0.01 1/s, u0=40 m/s
      REAL(wp), INTENT(in)  :: lon,                  & ! longitude in radians
                               lat,                  & ! latitude in radians
                               height                  ! height of the level in m
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      REAL(wp), INTENT(out) :: u_wind,               & ! zonal wind in m/s
                               v_wind,               & ! meridional wind in m/s
                               temperature,          & ! temperature in K
                               surface_geopotential, & ! surface geopotential in m^2/s^2
                               surface_pressure

!-----------------------------------------------------------------------
!     test case parameters
!-----------------------------------------------------------------------
      REAL(wp),PARAMETER :: p0      = 100000._wp,    & ! reference pressure
                            T0      = 300._wp,       & ! reference temperature
                            Lz      = 20.e3_wp,      & ! vertical wave length, 20 km
                            delta_theta = 10._wp       ! potential temp. perturbation amplitude

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
      REAL(wp) :: this_RR                                   ! half width
      REAL(wp) :: sin_tmp, cos_tmp
      REAL(wp) :: tmp1
      REAL(wp) :: theta                                ! potential temperature
      REAL(wp) :: theta_mean
      REAL(wp) :: pres                                 ! pressure
      REAL(wp) :: r                                    ! great circle distance

!-----------------------------------------------------------------------
!     more test case parameters
!-----------------------------------------------------------------------
      REAL(wp) :: N2,                              &   ! squared Brunt-Vaisala frequency
                  S,                               &   ! parameter
                  u0,                              &   ! background wind speed
                  lambda0,                         &   ! center point in longitudes
                  phi0,                            &   ! center point in latitudes
                  gw_omega,                        &   ! rotation
                  p_eq

      this_RR      = a/3._wp   ! half width
!-----------------------------------------------------------------------
!    initialize parameters
!-----------------------------------------------------------------------
      lambda0 = pi                                     ! center point in longitudes
      p_eq    = p0                                     ! surface pressure at the equator

      SELECT CASE (choice)
      CASE (0)
        N2 = 1.e-4_wp                                  ! squared Brunt Vaisala frequency N^2
        u0 = 0._wp                                     ! background wind speed
        phi0 = 0._wp                                   ! center point in latitudes (0 deg)
        gw_omega = 0._wp                               ! no rotation
      CASE (1)
        N2 = (g*g)/(cp*T0)                             ! squared Brunt Vaisala frequency N^2
        u0 = 0._wp                                     ! background wind speed
        phi0 = 0._wp                                   ! center point in latitudes (0 deg)
        gw_omega = 0._wp                               ! no rotation
      CASE (2)
        N2 = (g*g)/(cp*T0)                             ! squared Brunt Vaisala frequency N^2
        u0 = 40._wp                                    ! background wind speed
        phi0 = 0._wp                                   ! center point in latitudes (0 deg)
        gw_omega = 0._wp                               ! no rotation
      CASE (3)
        N2 = (g*g)/(cp*T0)                             ! squared  Brunt Vaisala frequency N^2
        u0 = 0._wp                                     ! background wind speed
        phi0 = pi/4._wp                                ! center point in latitudes (45 deg)
        gw_omega = omega                               ! Earth's rotation
      CASE DEFAULT
        N2 = 0._wp                                     ! squared  Brunt Vaisala frequency N^2
        u0 = 0._wp                                     ! background wind speed
        phi0 = 0._wp                                   ! center point in latitudes (45 deg)
        gw_omega = 0._wp                               ! Earth's rotation
        CALL finish('','')
      END SELECT
      S = g*g/(cp*N2)                                  ! parameter

!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * COS(lat)
      v_wind = 0._wp
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0._wp
!-----------------------------------------------------------------------
!     initialize surface pressure
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2._wp * g*g * kappa) * (u0/a + 2._wp * gw_omega)
      surface_pressure = p_eq * EXP( - tmp1*((SIN(lat))**2) )
!-----------------------------------------------------------------------
!     initialize temperature
!-----------------------------------------------------------------------
      pres    = p0 * ( (1._wp - S/T0) + S/T0 * EXP(- (N2*height)/g) )**(cp/Rd)
      sin_tmp = SIN(lat) * SIN(phi0)
      cos_tmp = COS(lat) * COS(phi0)
      theta_mean = T0 /( T0/S * ((pres/p0)**kappa - 1._wp) + 1._wp )
      r = a * ACOS (sin_tmp + cos_tmp*COS(lon-lambda0))     ! great circle distance with radius
      IF (r < this_RR) THEN
        tmp1 = 0.5_wp * (1._wp + COS(pi*r/this_RR))
      ELSE
        tmp1 = 0._wp
      ENDIF
      theta = theta_mean + delta_theta * tmp1 * SIN(2._wp*pi*height/Lz)
      temperature = theta * (pres/p0)**kappa

  END SUBROUTINE gravity_wave

END MODULE mo_ncar_testcases
