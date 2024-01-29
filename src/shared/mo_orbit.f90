! Computation of orbital parameters for use in radiative transfer
! calculations (among other things)
!
! The code in this module provides routines to calculate the
! distance to, the right ascension, and the declination of the sun.
!
! Two orbital models are provided:
!   - orbit_vsop87 : standard and accurate model
!   - orbit_kepler : simple model, appropriate for idealized work as
!                    well as an inquiry function for the  declination
!                    (if it has been calculated).
!
! The code is based on the following published articles and books and
! do neither describe patented, copyrighted, or licensed algorithms or
! code.
!
! P. Bretagnon and G. Francou. Planetary theories in rectangular and
! spherical variables - VSOP 87 solutions. Astronomy and Astrophysics,
! 202:309-315, 1988.
!
! Simon, J.L., P. Bretagnon, et al., "Numerical expressions for precession 
! formulae and mean elements for the Moon and the planets", Astron. and 
! Astrophys. 282 (1994) 663-683.
!
! J. Laskar, F. Joutel, and F. Boudin. Orbital, precessional, and
! insolation quantities for the earth from -20 myr to +10 myr. Astronomy
! and Astrophysics, 270:522-533, 1993.
!
! Laskar, J.: 1986, "Secular Terms of Classical Planetary Theories
! Using the Results of General Theory". Astron. and Astrophys. 157,59
!
! A. Berger. Obliquity and precession for the last 5,000,000
! years. Astronomy and Astrophysics, 51:127-135, 1976.
!
! Monin, A. S.: An Introduction to the Theory of Climate  D. Reidel 
! Publishing Company, Dordrecht, 1986 (pp 10-12).
!
! J. Meeus. Astronomical Algorithms. Willmann-Bell, Richmond, 2 edition,
! 1998.
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

MODULE mo_orbit

  USE mo_kind,           ONLY : wp, i8
  USE mo_math_constants, ONLY : pi           ,& ! pi
       &                        twopi => pi2 ,& ! pi*2
       &                        deg2rad         ! pi/180
  USE mo_exception,      ONLY : finish, message, message_text,warning, print_value
  USE mtime,             ONLY : julianday, newJulianday, deallocateJulianday, getJulianDayFromDatetime, &
       &                        newDateTime, deallocateDateTime, datetime, no_of_ms_in_a_day, &
       &                        getNoOfDaysInYearDateTime

  IMPLICIT NONE
  PRIVATE 
  PUBLIC :: orbit_kepler, orbit_vsop87, inquire_declination, &
            get_orbit_times

  TYPE terms 
    REAL(wp) :: A
    REAL(wp) :: B
    REAL(wp) :: C
  END TYPE terms

  LOGICAL,  SAVE      :: initialized = .FALSE.
  REAL(wp), PARAMETER :: sec2rad     = deg2rad/3600.0_wp

  REAL(wp), SAVE      :: declination  = 0.0_wp
  LOGICAL , SAVE      :: decl_warning = .TRUE.

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief Simple model for sun-earth geometry versus time of year
  !!
  !! @remarks 
  !!   This routine computes three orbital parameters depending on the time of
  !!   the day as well as of the year (both in radians). The main parameters 
  !!   are the eccentricity (cecc); the obliquity (cobld) and the longitude of
  !!   of perihelion (clonp).
  !
  SUBROUTINE orbit_kepler (cecc, cobld, clonp, time, rasc_sun, decl_sun, dist_sun)

    REAL(wp), PARAMETER   :: ceps = 1.0e-9_wp        

    REAL(wp), INTENT(in)  :: cecc  !< Eccentricity of the Kepler orbit
    REAL(wp), INTENT(in)  :: cobld !< Obliquity of the Earth axis [Deg]
    REAL(wp), INTENT(in)  :: clonp !< Longitude of perihelion [Deg]
    REAL(wp), INTENT(in)  :: time  !< Time of year rel. to vernal equinox [rad]
    REAL(wp), INTENT(out) :: &
         rasc_sun,           & !< Right Ascension of the Sun
         decl_sun,           & !< Declination of the Sun
         dist_sun              !< Distance to the Sun in Astronomical Units

    INTEGER  :: iter
    REAL(wp) ::    &
         obl_rad,  & ! oblquity [radians]
         phl_rad     ! Longitude of perihelion [radians]
    REAL(wp) :: sq_ecc, lmbd, delt, guess, a, b, z1, z2, z3, diff, cos_e, big_e

    !
    ! set local variables for eccentricity and obliquity
    obl_rad = cobld*deg2rad
    phl_rad = clonp*deg2rad
    WRITE(message_text, '(a14,f9.3,a23,f9.3)') &
         & ' eccentricity=',cecc,'  obliquity in degrees=',cobld
    CALL message('orbit_kepler (mo_radiation_orbit):',message_text)
    !
    ! Calculation of eccentric anomaly (big_e) of vernal equinox using
    ! Lacaille's formula.
    ! --------------------------------
    sq_ecc = SQRT((1.0_wp+cecc)/(1.0_wp-cecc))
    big_e  = 2.0_wp*ATAN(TAN(0.5_wp*phl_rad)/sq_ecc)
    !
    ! Calculation of true anomaly of vernal equinox (Kepler) 
    ! --------------------------------
    !
    ! -- make first guess of the eccentric anomaly (big_e) a correction is
    ! applied for special cases where the first guest will not lead to 
    ! convergence.  
    !
    delt = time-(big_e-cecc*SIN(big_e))
    guess = delt/(1.0_wp-cecc)
    big_e = delt
    IF ( cecc > 0.975_wp .AND. ABS(delt) <  0.52359_wp) THEN  
      a     = (1.0_wp-cecc)/(4.0_wp*cecc+0.5_wp)
      b     = delt/(8.0_wp*cecc+1.0_wp)
      z1    = SIGN(SQRT(b*b+a*a*a), b)
      z2    = SIGN((ABS(b+z1))**1.5_wp,b+z1) - 0.5_wp*a
      z3    = z2-(0.078_wp*z2**5)/(1.0_wp+cecc)        
      big_e = delt+cecc*(3.0_wp*z3 - 4.0_wp*z3*z3*z3)
    END IF
    !
    ! --- Newton-Raphson Iteration
    !
    iter = 0
    diff  = guess - big_e
    DO WHILE (iter < 25 .AND. ABS(diff) >= ceps)
      guess = big_e
      cos_e = COS(big_e)
      big_e = (delt+cecc*(SIN(big_e)-big_e*cos_e)) / (1.0_wp-cecc*cos_e)
      diff  = guess - big_e
      iter  = iter+1      
    END DO
    IF (iter >= 25 ) CALL finish('orbit_kepler','Eccentric anomaly not found!')
    !
    ! finalize output values and make declination available for inquiry
    ! --------------------------------
    !
    ! --- lamda defines the true longitude of the earth (ccw from vernl equinox)
    !
    lmbd    = 2.0_wp*ATAN(sq_ecc*TAN(big_e*0.5_wp)) + phl_rad
    rasc_sun  = ATAN2(COS(obl_rad)*SIN(lmbd), COS(lmbd))
    decl_sun  = ASIN (SIN(obl_rad)*SIN(lmbd))
    dist_sun  = 1.0_wp-cecc*COS(big_e)

    declination = decl_sun
    initialized = .TRUE.

  END SUBROUTINE orbit_kepler
  !-----------------------------------------------------------------------------
  !>
  !! @brief Advanced model for calculating orbital parameters
  !!
  !! @remarks 
  !!   This routine computes orbital parameters depending on the Julian Day.
  !!   The model is based on the Variations Seculaires des Orbites Plan
  !!   Planetaires (VSOP) method, as implemented in the VSOP87 model. The VSOP 
  !!   model provides the sun-earth distance, right ascention and declination
  !!   of the sun and the hour angle as output.  This is the standard orbital 
  !!   model used by ECHAM/ICON.
  !
  SUBROUTINE orbit_vsop87 (julian_day, rasc_sun, decl_sun, dist_sun)

    REAL(wp), INTENT(in)  :: julian_day
    REAL(wp), INTENT(out) :: &
         rasc_sun,           & !< Right Ascension of the Sun
         decl_sun,           & !< Declination of the Sun
         dist_sun              !< Distance to the Sun in Astronomical Units

    REAL(wp), PARAMETER     :: sun_helio(3) = (/ 0.0_wp, 0.0_wp, 0.0_wp/)
    REAL(wp)                :: t, obl, coords(2), sun_geo(3)
    !
    ! --- Preliminary calculations 
    !
    IF (ABS(julian_day) >= 1.e+07_wp) CALL finish('orbit_vsop87',             &
         'Orbital model not valid for extreme times')
    t       = jd2t(julian_day) ! time in centuries since J2000.0
    obl     = obliquity(t)
    coords  = nutation (t)     ! returns nutation longitude and obliquity
    sun_geo = helio2geo(sun_helio,earth_position (t))
    !
    ! --- Finalize output values and make declination available for inquiry
    !     ecl2equ returns right ascension and declination
    !
    coords   = ecl2equ (sun_geo(1)+coords(1), sun_geo(2), obl+coords(2))
    rasc_sun = coords(1)
    decl_sun = coords(2)
    dist_sun = sun_geo(3)
    CALL aberration (t, obl, rasc_sun, decl_sun)

    declination = decl_sun
    initialized = .TRUE.

  END SUBROUTINE orbit_vsop87
  !-----------------------------------------------------------------------------
  !>
  !! @brief corrects the right ascension and declination
  !
  SUBROUTINE aberration (t, obl, rasc, decl)

    REAL(wp), INTENT(in)    :: t, obl
    REAL(wp), INTENT(inout) :: rasc, decl

    REAL(wp), PARAMETER :: abconst = 20.49552_wp*sec2rad
    REAL(wp)            :: cosobl, sinobl, cosdecl, sindecl, cosrasc, sinrasc &
         &                ,coslon, sinlon, pir, tmp, coords(2)
    !
    cosobl  = COS(obl)
    sinobl  = SIN(obl)
    cosrasc = COS(rasc)
    sinrasc = SIN(rasc)
    cosdecl = COS(decl)
    sindecl = SIN(decl)
    !
    coords = sun_position (t) ! returns longitude and eccentric
    coslon = COS (coords(1))
    sinlon = SIN (coords(1))
    ! 
    ! --- the next three terms are corrections for the FK5 system 
    !
    pir    = (102.93735_wp+t*(1.71954_wp+t*0.00046_wp))*deg2rad
    coslon = coslon-coords(2)*COS(pir)
    sinlon = sinlon-coords(2)*SIN(pir)
    !
    rasc = rasc-abconst*(cosrasc*coslon*cosobl+sinrasc*sinlon)/cosdecl
    tmp  = coslon*(sinobl*cosdecl-sinrasc*sindecl*cosobl)+cosrasc*sindecl*sinlon
    decl = decl-abconst*tmp

  END SUBROUTINE aberration
  !-----------------------------------------------------------------------------
  !>
  !! @brief Ecliptic to Equatorial transformation of coordinates
  !
  PURE FUNCTION ecl2equ(l, b, obl)

    REAL(wp)              :: ecl2equ(2) !< right ascension and delcination
    REAL(wp), INTENT(in)  :: l, b, obl  !< longitude, latitude, obliquity

    REAL(wp) :: sinobl, cosobl, sinl, cosl, sinb, cosb

    sinobl = SIN(obl)
    cosobl = COS(obl)
    sinl = SIN(l)
    cosl = COS(l)
    sinb = SIN(b)
    cosb = COS(b)

    ecl2equ(1) = ATAN2(cosb*sinl*cosobl-sinb*sinobl,cosb*cosl)
    IF (ecl2equ(1) < 0) ecl2equ(1) = ecl2equ(1)+twopi
    ecl2equ(2) = ASIN(sinb*cosobl+cosb*sinobl*sinl)

  END FUNCTION ecl2equ
  !-----------------------------------------------------------------------------
  !>
  !! @brief Converts from spherical to rectangular coordinates
  !
  PURE FUNCTION sph2rect(s)

    REAL(wp) :: sph2rect(3)
    REAL(wp), INTENT(in)  :: s(3)

    sph2rect(1) = s(3)*COS(s(1))*COS(s(2))
    sph2rect(2) = s(3)*SIN(s(1))*COS(s(2))
    sph2rect(3) = s(3)*SIN(s(2))

  END FUNCTION sph2rect
  !-----------------------------------------------------------------------------
  !>
  !! @brief Converts from rectangular to spherical coordinates
  !
  PURE FUNCTION rect2sph (r)

    REAL(wp) :: rect2sph(3)
    REAL(wp), INTENT(in)  :: r(3)

    rect2sph(1) = ATAN2(r(2),r(1))
    IF (rect2sph(1) < 0.0_wp) rect2sph(1) = rect2sph(1)+twopi
    rect2sph(3) = SQRT(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
    rect2sph(2) = ASIN(r(3)/rect2sph(3))

  END FUNCTION rect2sph
  !-----------------------------------------------------------------------------
  !>
  !! @brief Converts Heliocentric to Geocentric coordinates
  !
  PURE FUNCTION helio2geo(shelio, searth)

    REAL(wp) :: helio2geo(3) 
    REAL(wp), INTENT(in)  :: shelio(3), searth(3)

    helio2geo(1:3)   = rect2sph (sph2rect (shelio) - sph2rect (searth) )

  END FUNCTION helio2geo
  !-----------------------------------------------------------------------------
  !>
  !! @brief Calculations longitude and obliquity of nutation constant
  !
  PURE FUNCTION nutation(t)

    REAL(wp) :: nutation(2)    !< Longitude and Obliquity of Nutation Constant
    REAL(wp), INTENT(in)  :: t !< Time in Centuries since J2000.0

    REAL(wp), PARAMETER :: po(0:3) = (/ &
         &  125.0445222_wp,  -1934.1362608_wp,  0.00207833_wp,  2.220e-6_wp /)
    REAL(wp), PARAMETER :: ps(0:3) = (/ &
         &  357.5277233_wp,  35999.0503400_wp, -0.00016030_wp, -3.330e-6_wp /)
    REAL(wp), PARAMETER :: pm(0:3) = (/ &
         &  134.9629814_wp, 477198.8673981_wp,  0.00869720_wp,  1.778e-5_wp /)
    REAL(wp), PARAMETER :: pf(0:3) = (/ &
         &  93.2719103_wp, 483202.0175381_wp, -0.00368250_wp,   3.056e-6_wp /)
    REAL(wp), PARAMETER :: pd(0:3) = (/ &
         &  297.8503631_wp, 445267.1114800_wp, -0.00191420_wp,  5.278e-6_wp /)

    REAL(wp) :: &
         xo,  & !< longitude of mean ascending node of the moons mean orbit
         xs,  & !< mean anomaly of the sun
         xm,  & !< mean anomoaly of the moon
         xf,  & !< moon argument of latitude
         xd     !< mean elongation of the moon from the sun

    REAL(wp) :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, nutlon, nutobl

    xo = MODULO ((po(0) + t*(po(1) + t*(po(2) + t*po(3))))*deg2rad,twopi)
    xs = MODULO ((ps(0) + t*(ps(1) + t*(ps(2) + t*ps(3))))*deg2rad,twopi)
    xm = MODULO ((pm(0) + t*(pm(1) + t*(pm(2) + t*pm(3))))*deg2rad,twopi)
    xf = MODULO ((pf(0) + t*(pf(1) + t*(pf(2) + t*pf(3))))*deg2rad,twopi)
    xd = MODULO ((pd(0) + t*(pd(1) + t*(pd(2) + t*pd(3))))*deg2rad,twopi)

    n1 = (-0.01742_wp * t - 17.1996_wp)*SIN(xo)                               &
         + (-1.3187_wp - 0.00016_wp * t)*SIN(2*xf-2*xd+2*xo)                  &
         + (-0.2274_wp - 0.00002_wp * t)*SIN(2*xf+2*xo)                       &
         + (0.20620_wp + 0.00002_wp * t)*SIN(2*xo)                            &
         + (0.14260_wp - 0.00034_wp * t)*SIN(xs)                              &
         + (0.07120_wp + 0.00001_wp*t)*SIN(xm)

    n2 = + (-0.0517_wp + 0.00012_wp*t)*SIN(xs+2*xf-2*xd+2*xo)                 &
         + (-0.0386_wp-0.00004_wp*t)*SIN(2*xf+xo)-0.0301_wp*SIN(xm+2*xf+2*xo) &
         + (0.0217_wp - 0.00005_wp*t)*SIN(-xs+2*xf-2*xd+2*xo)                 &
         -0.0158_wp*SIN(xm-2*xd) + (0.0129_wp + 0.00001_wp*t)*SIN(2*xf-2*xd+xo)

    n3 = +0.0123_wp*SIN(-xm+2*xf+2*xo)                                        &
         +0.0063_wp*SIN(2*xd) + (0.0063_wp + 0.00001_wp*t)*SIN(xm+xo)         &
         -0.0059_wp*SIN(-xm+2*xf+2*xd+2*xo)                                   &
         + (-0.0058_wp - 0.00001_wp*t)*SIN(-xm+xo) -0.0051_wp*SIN(xm+2*xf+xo)

    n4 = +0.0048_wp*SIN(2*xm-2*xd)      + 0.0046_wp*SIN(-2*xm+2*xf+xo)        &
         -0.0038_wp*SIN(2*xf+2*xd+2*xo) - 0.0031_wp*SIN(2*xm+2*xf+2*xo)       &
         +0.0029_wp*SIN(2*xm)           + 0.0029_wp*SIN(xm+2*xf-2*xd+2*xo)

    n5 = +0.0026_wp*SIN(2*xf)        - 0.0022_wp*SIN(2*xf-2*xd)               &
         +0.0021_wp*SIN(-xm+2*xf+xo) + (0.0017_wp - 0.00001_wp*t)*SIN(2*xs)   &
         +0.0016_wp*SIN(-xm+2*xd+xo)                                          &
         +(-0.0016_wp + 0.00001_wp*t)*SIN(2*xs+2*xf-2*xd+2*xo)

    n6 = -0.0015_wp*SIN(xs+xo)           - 0.0013_wp*SIN(xm-2*xd+xo)          &
         -0.0012_wp*SIN(-xs+xo)          + 0.0011_wp*SIN(2*xm-2*xf)           &
         -0.0010_wp*SIN(-xm+2*xf+2*xd+xo)- 0.0008_wp*SIN(xm+2*xf+2*xd+2*xo)

    n7 = +0.0007_wp*SIN(xs+2*xf+2*xo)    - 0.0007_wp*SIN(xm+xs-2*xd)          &
         -0.0007_wp*SIN(-xs+2*xf+2*xo)   - 0.0007_wp*SIN(2*xf+2*xd+xo)        &
         +0.0006_wp*SIN(xm+2*xd)         + 0.0006_wp*SIN(2*xm+2*xf-2*xd+2*xo)

    n8 = +0.0006_wp*SIN(xm+2*xf-2*xd+xo) - 0.0006_wp*SIN(-2*xm+2*xd+xo)       &
         -0.0006_wp*SIN(2*xd+xo)         + 0.0005_wp*SIN(xm-xs)               &
         -0.0005_wp*SIN(-xs+2*xf-2*xd+xo)- 0.0005_wp*SIN(-2*xd+xo)

    n9 = -0.0005_wp*SIN(2*xm+2*xf+xo)    - 0.0003_wp*SIN(-2*xm+2*xf+2*xo)     &
         +0.0004_wp*SIN(2*xm-2*xd+xo)    + 0.0004_wp*SIN(xs+2*xf-2*xd+xo)     &
         -0.0003_wp*SIN(xm-xs+2*xf+2*xo) - 0.0003_wp*SIN(-xm-xs+2*xf+2*xd+2*xo)

    n10 =-0.0003_wp*SIN(3*xm+2*xf+2*xo) - 0.0003_wp*SIN(-xs+2*xf+2*xd+2*xo)   &
         -0.0003_wp*SIN(xm-xs-xd)       - 0.0004_wp*SIN(xm-xd)                &
         -0.0004_wp*SIN(xs-2*xd)        + 0.0004_wp*SIN(xm-2*xf)

    n11 =-0.0004_wp*SIN(xd) - 0.0003_wp*SIN(xm+xs) + 0.0003_wp*SIN(xm+2*xf)

    nutlon = n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11

    n1 =   (0.00089_wp*t + 9.2025_wp)*COS(xo)                                 &
         + (0.5736_wp - 0.00031_wp*t)*COS(2*xf-2*xd+2*xo)                     &
         + (0.0977_wp - 0.00005_wp*t)*COS(2*xf+2*xo)                          &
         + (-0.0895_wp + 0.00005_wp*t)*COS(2*xo)                              &
         + (0.0054_wp - 0.00001_wp*t)*COS(xs) - 0.0007_wp*COS(xm)

    n2 = + (0.0224_wp-0.00006_wp*t)*COS(xs+2*xf-2*xd+2*xo)                    &
         + 0.0200_wp*COS(2*xf+xo) + (0.0129_wp-0.00001_wp*t)*COS(xm+2*xf+2*xo)&
         + (-0.0095_wp + 0.00003_wp*t)*COS(-xs+2*xf-2*xd+2*xo)                &
         - 0.0070_wp*COS(2*xf-2*xd+xo) - 0.0053_wp*COS(-xm+2*xf+2*xo)

    n3 =  -0.0033_wp*COS(xm+xo)         + 0.0026_wp*COS(-xm+2*xf+2*xd+2*xo)   &
         + 0.0032_wp*COS(-xm+xo)        + 0.0027_wp*COS(xm+2*xf+xo)           &
         - 0.0024_wp*COS(-2*xm+2*xf+xo) + 0.0016_wp*COS(2*xf+2*xd+2*xo)

    n4 = +0.0013_wp*COS(2*xm+2*xf+2*xo)     - 0.0012_wp*COS(xm+2*xf-2*xd+2*xo)&
         -0.0010_wp*COS(-xm+2*xf+xo)        - 0.0008_wp*COS(-xm+2*xd+xo)      &
         +0.0007_wp*COS(2*xs+2*xf-2*xd+2*xo)+ 0.0009_wp*COS(xs+xo)

    n5 = +0.0007_wp*COS(xm-2*xd+xo)       + 0.0006_wp*COS(-xs+xo)             &
         +0.0005_wp*COS(-xm+2*xf+2*xd+xo) + 0.0003_wp*COS(xm+2*xf+2*xd+2*xo)  &
         -0.0003_wp*COS(xs+2*xf+2*xo)     + 0.0003_wp*COS(-xs+2*xf+2*xo)

    n6 = +0.0003_wp*COS(2*xf+2*xd+xo)    - 0.0003_wp*COS(2*xm+2*xf-2*xd+2*xo) &
         -0.0003_wp*COS(xm+2*xf-2*xd+xo) + 0.0003_wp*COS(-2*xm+2*xd+xo)       &
         +0.0003_wp*COS(2*xd+xo)         + 0.0003_wp*COS(-xs+2*xf-2*xd+xo)

    n7 = +0.0003_wp*COS(-2*xd+xo)+0.0003_wp*COS(2*xm+2*xf+xo)

    nutobl = n1+n2+n3+n4+n5+n6+n7

    nutation(:) = (/nutlon,nutobl/)*sec2rad

  END FUNCTION nutation
  !-----------------------------------------------------------------------------
  !>
  !! @brief Calculates obliquity given Julian time
  !! 
  !
  PURE FUNCTION obliquity(t) 

    REAL(wp) :: obliquity     !< obliquity
    REAL(wp), INTENT(in) :: t !< Centuries since J2000.0

    REAL(wp), PARAMETER ::  c(0:10) = (/                                &
         &  84381.448_wp, -4680.93_wp, -1.55_wp, 1999.25_wp, -51.38_wp, &
         &    -249.67_wp,  -39.050_wp,  7.12_wp,   27.87_wp,   5.79_wp, &
         &       2.45_wp /)

    REAL(wp)             :: x

    x         = 0.01_wp*t ! convert to years since J2000.
    obliquity = sec2rad*(c(0)+x*(c(1) + x*(c(2) + x*(c(3) + x*(c(4) + x*(c(5) &
         &       + x*(c(6) + x*(c(7) + x*(c(8) + x*(c(9) + x*c(10)))))))))))

  END FUNCTION obliquity
  !-----------------------------------------------------------------------------
  !>
  !! @brief Julian Day to Time (in centuries since J2000.0)
  !! 
  !! @note 2451545 is the Julian Day of 2000-01-01.
  !
  PURE FUNCTION jd2t(jd)

    REAL(wp) :: jd2t           !< centuries since J2000.0
    REAL(wp), INTENT(in) :: jd !< Julian Day

    jd2t = (jd-2451545.0_wp)/36525.0_wp 

  END FUNCTION jd2t
  !-----------------------------------------------------------------------------
  !>
  !! @brief calculate the longitude and eccentricity of the sun
  !
  PURE FUNCTION sun_position (t)

    REAL(wp) :: sun_position(2)
    REAL(wp), INTENT(in)  :: t ! Time in centuries since J2000.0

    REAL(wp) :: l, m, c, e2, e

    l  = (280.46646_wp+t*(36000.76983_wp+t*0.0003032_wp))*deg2rad;
    m  = (357.52910_wp+t*(35999.05028_wp-t*0.0001561_wp))*deg2rad;
    e  = 0.016708617_wp-t*(0.000042040_wp+t*0.0000001236_wp);
    e2 = e*e
    !
    ! --- equation of the center, in terms of e and m
    !
    c = e*(2-0.25_wp*e2)*SIN(m)+1.25_wp*e2*SIN(2*m) &
         +1.0833333333_wp*e*e2*SIN(3*m)

    sun_position(1) = l+c ! longitude
    sun_position(2) = e   ! eccentricity

  END FUNCTION sun_position
  !-----------------------------------------------------------------------------
  !>
  !! @brief calculate the heliocentric earth position at time t
  !! 
  !! @remarks
  !!   Incorporated from mo_vosp87 and rewritten as a function. Trailing zeros
  !!   removed from data, and continuation lines limited to no more than 40 
  !!   (per FORTRAN standard).
  !
  PURE FUNCTION earth_position (t)

    REAL(wp) :: earth_position(3) !< long, lat, radius
    REAL(wp), INTENT(in)  :: t    !< number of centuries since J2000
    !
    ! --- Data given are based on VSOP87D:  Some terms are broken up so that the
    !     number of continuation linds conform with the F95 standard, hence 
    !     L1 = (/L1A,L1B/) and R1 = (/R1A,R1B/)
    !
    TYPE (terms), PARAMETER :: L1A(32) = (/ &
         terms(1.75347045673e+00_wp,0.0000000e+00_wp, 0.00000000000000e+00_wp),&
         terms(3.341656456e-02_wp,4.66925680417e+00_wp,6.2830758499914e+03_wp),&
         terms(3.4894275e-04_wp,4.62610241759e+00_wp, 1.25661516999828e+04_wp),&
         terms(3.497056e-05_wp, 2.74411800971e+00_wp, 5.75338488489680e+03_wp),&
         terms(3.417571e-05_wp, 2.82886579606e+00_wp, 3.52311834900000e+00_wp),&
         terms(3.135896e-05_wp, 3.62767041758e+00_wp, 7.77137714681205e+04_wp),&
         terms(2.676218e-05_wp, 4.41808351397e+00_wp, 7.86041939243920e+03_wp),&
         terms(2.342687e-05_wp, 6.13516237631e+00_wp, 3.93020969621960e+03_wp),&
         terms(1.324292e-05_wp, 7.42463563520e-01_wp, 1.15067697697936e+04_wp),&
         terms(1.273166e-05_wp, 2.03709655772e+00_wp, 5.29690965094600e+02_wp),&
         terms(1.199167e-05_wp, 1.10962944315e+00_wp, 1.57734354244780e+03_wp),&
         terms(9.902500e-06_wp, 5.23268129594e+00_wp, 5.88492684658320e+03_wp),&
         terms(9.018550e-06_wp, 2.04505443513e+00_wp, 2.62983197998000e+01_wp),&
         terms(8.572230e-06_wp, 3.50849156957e+00_wp, 3.98149003408200e+02_wp),&
         terms(7.797860e-06_wp, 1.17882652114e+00_wp, 5.22369391980220e+03_wp),&
         terms(7.531410e-06_wp, 2.53339053818e+00_wp, 5.50755323866740e+03_wp),&
         terms(5.052640e-06_wp, 4.58292563052e+00_wp, 1.88492275499742e+04_wp),&
         terms(4.923790e-06_wp, 4.20506639861e+00_wp, 7.75522611324000e+02_wp),&
         terms(3.566550e-06_wp, 2.91954116867e+00_wp, 6.73103028000000e-02_wp),&
         terms(3.170870e-06_wp, 5.84901952218e+00_wp, 1.17906290886588e+04_wp),&
         terms(2.841250e-06_wp, 1.89869034186E+00_wp, 7.96298006816400e+02_wp),&
         terms(2.710390e-06_wp, 3.14886076490e-01_wp, 1.09770788046990e+04_wp),&
         terms(2.428100e-06_wp, 3.44811409060e-01_wp, 5.48677784317500e+03_wp),&
         terms(2.061600e-06_wp, 4.80646606059e+00_wp, 2.54431441988340e+03_wp),&
         terms(2.053850e-06_wp, 1.86947813692e+00_wp, 5.57314280143310e+03_wp),&
         terms(2.022610e-06_wp, 2.45767795458e+00_wp, 6.06977675455340e+03_wp),&
         terms(1.555160e-06_wp, 8.33060738070e-01_wp, 2.13299095438000e+02_wp),&
         terms(1.322120e-06_wp, 3.41118275555e+00_wp, 2.94246342329160e+03_wp),&
         terms(1.261840e-06_wp, 1.08302630210e+00_wp, 2.07753954924000e+01_wp),&
         terms(1.151320e-06_wp, 6.45449116830e-01_wp, 9.80321068200000e-01_wp),&
         terms(1.028510e-06_wp, 6.35998467270e-01_wp, 4.69400295470760e+03_wp),&
         terms(1.018950e-06_wp, 9.75692218240e-01_wp, 1.57208387848784e+04_wp)/)
    TYPE (terms), PARAMETER :: L1B(32) = (/ &
         terms(1.017240e-06_wp, 4.26679821365e+00_wp, 7.11354700080000e+00_wp),&
         terms(9.920600e-07_wp, 6.20992940258e+00_wp, 2.14616541647520e+03_wp),&
         terms(9.760700e-07_wp, 6.81012722700e-01_wp, 1.55420399434200e+02_wp),&
         terms(8.580300e-07_wp, 5.98322631256E+00_wp,1.610006857376741E+05_wp),&
         terms(8.512800e-07_wp, 1.29870743025e+00_wp, 6.27596230299060e+03_wp),&
         terms(8.471100e-07_wp, 3.67080093025e+00_wp,7.143069561812909e+04_wp),&
         terms(7.963700e-07_wp, 1.80791330700e+00_wp, 1.72601546546904e+04_wp),&
         terms(7.875600e-07_wp, 3.03698313141e+00_wp, 1.20364607348882e+04_wp),&
         terms(7.465100e-07_wp, 1.75508916159e+00_wp, 5.08862883976680e+03_wp),&
         terms(7.387400e-07_wp, 3.50319443167e+00_wp, 3.15468708489560e+03_wp),&
         terms(7.354700e-07_wp, 4.67926565481e+00_wp,8.018209311238001e+02_wp),&
         terms(6.962700e-07_wp, 8.32975969660e-01_wp, 9.43776293488700e+03_wp),&
         terms(6.244900e-07_wp, 3.97763880587E+00_wp, 8.82739026987480e+03_wp),&
         terms(6.114800e-07_wp, 1.81839811024e+00_wp, 7.08489678111520e+03_wp),&
         terms(5.696300e-07_wp, 2.78430398043e+00_wp, 6.28659896834040e+03_wp),&
         terms(5.611600e-07_wp, 4.38694880779e+00_wp, 1.41434952424306e+04_wp),&
         terms(5.557700e-07_wp, 3.47006009062e+00_wp, 6.27955273164240e+03_wp),&
         terms(5.199200e-07_wp, 1.89149458340e-01_wp, 1.21395535091068e+04_wp),&
         terms(5.160500e-07_wp, 1.33282746983e+00_wp, 1.74801641306700e+03_wp),&
         terms(5.114500e-07_wp, 2.83068645010e-01_wp, 5.85647765911540e+03_wp),&
         terms(4.900000e-07_wp, 4.87350650330e-01_wp, 1.19444701022460e+03_wp),&
         terms(4.103600e-07_wp, 5.36817351402E+00_wp,8.429241266466601e+03_wp),&
         terms(4.093800e-07_wp, 2.39850881707e+00_wp, 1.96510484810980e+04_wp),&
         terms(3.920000e-07_wp, 6.16832995016e+00_wp, 1.04473878396044e+04_wp),&
         terms(3.677000e-07_wp, 6.04133859347e+00_wp, 1.02132855462110e+04_wp),&
         terms(3.659600e-07_wp, 2.56955238628e+00_wp, 1.05938193018920e+03_wp),&
         terms(3.595400e-07_wp, 1.70876111898e+00_wp, 2.35286615377180e+03_wp),&
         terms(3.556600e-07_wp, 1.77597314691e+00_wp, 6.81276681508600e+03_wp),&
         terms(3.329100e-07_wp, 5.93094994590E-01_wp, 1.77898456197850e+04_wp),&
         terms(3.041200e-07_wp, 4.42944641350e-01_wp,8.399684731811189e+04_wp),&
         terms(3.004700e-07_wp, 2.73975123935e+00_wp, 1.34986740965880e+03_wp),&
         terms(2.535200e-07_wp, 3.16470953405e+00_wp, 4.69047983635860e+03_wp)/)

    TYPE (terms), PARAMETER :: L1(64) = (/L1A,L1B/)

    TYPE (terms), PARAMETER :: L2(34) = (/ &
         terms(6.28331966747491e+03_wp, 0.0000000e+00_wp, 0.0000000000e+00_wp),&
         terms(2.06058863e-03_wp,2.67823455584e+00_wp, 6.2830758499914e+03_wp),&
         terms(4.30343e-05_wp, 2.63512650414e+00_wp, 1.256615169998280e+04_wp),&
         terms(4.25264e-06_wp, 1.59046980729e+00_wp, 3.523118349000000e+00_wp),&
         terms(1.19261e-06_wp, 5.79557487799e+00_wp, 2.629831979980000e+01_wp),&
         terms(1.08977e-06_wp, 2.96618001993e+00_wp, 1.577343542447800e+03_wp),&
         terms(9.34780e-07_wp, 2.59212835365e+00_wp, 1.884922754997420e+04_wp),&
         terms(7.21220e-07_wp, 1.13846158196e+00_wp, 5.296909650946000e+02_wp),&
         terms(6.77680e-07_wp, 1.87472304791e+00_wp, 3.981490034082000e+02_wp),&
         terms(6.73270e-07_wp, 4.40918235168e+00_wp, 5.507553238667400e+03_wp),&
         terms(5.90270e-07_wp, 2.88797038460e+00_wp, 5.223693919802200e+03_wp),&
         terms(5.59760e-07_wp, 2.17471680261e+00_wp, 1.554203994342000e+02_wp),&
         terms(4.54070e-07_wp, 3.98030798050e-01_wp, 7.962980068163999e+02_wp),&
         terms(3.63690e-07_wp, 4.66247398350e-01_wp, 7.755226113240000e+02_wp),&
         terms(2.89580e-07_wp, 2.64707383882e+00_wp, 7.113547000800000e+00_wp),&
         terms(2.08440e-07_wp, 5.34138275149e+00_wp, 9.803210682000000e-01_wp),&
         terms(1.90970e-07_wp, 1.84628332577e+00_wp, 5.486777843175000e+03_wp),&
         terms(1.85080e-07_wp, 4.96855124577e+00_wp, 2.132990954380000e+02_wp),&
         terms(1.72930e-07_wp, 2.99116864949e+00_wp, 6.275962302990600e+03_wp),&
         terms(1.62330e-07_wp, 3.21648304700e-02_wp, 2.544314419883400e+03_wp),&
         terms(1.58320e-07_wp, 1.43049285325e+00_wp, 2.146165416475200e+03_wp),&
         terms(1.46150e-07_wp, 1.20532366323e+00_wp, 1.097707880469900e+04_wp),&
         terms(1.24610e-07_wp, 2.83432285512e+00_wp, 1.748016413067000e+03_wp),&
         terms(1.18770e-07_wp, 3.25804815607e+00_wp, 5.088628839766800e+03_wp),&
         terms(1.18080e-07_wp, 5.27379790480e+00_wp, 1.194447010224600e+03_wp),&
         terms(1.15140e-07_wp, 2.07502418155e+00_wp, 4.694002954707600e+03_wp),&
         terms(1.06410e-07_wp, 7.66141992020e-01_wp, 5.535694028424000e+02_wp),&
         terms(9.96900e-08_wp, 1.30262991097e+00_wp, 6.286598968340400e+03_wp),&
         terms(9.72100e-08_wp, 4.23925472239E+00_wp, 1.349867409658800e+03_wp),&
         terms(9.45200e-08_wp, 2.69957062864e+00_wp, 2.427286039740000e+02_wp),&
         terms(8.57700e-08_wp, 5.64475868067e+00_wp, 9.517184062506000e+02_wp),&
         terms(7.57600e-08_wp, 5.30062664886e+00_wp, 2.352866153771800e+03_wp),&
         terms(6.38500e-08_wp, 2.65033984967E+00_wp, 9.437762934887000e+03_wp),&
         terms(6.10100e-08_wp, 4.66632584188e+00_wp, 4.690479836358600e+03_wp)/)

    TYPE (terms), PARAMETER :: L3(20) = (/ &
         terms(5.291887e-04_wp, 0.00000000000e+00_wp, 0.0000000000000e+00_wp), &
         terms(8.719837e-05_wp, 1.07209665242e+00_wp, 6.2830758499914e+03_wp), &
         terms(3.09125e-06_wp, 8.67288188320e-01_wp, 1.256615169998280e+04_wp),&
         terms(2.73390e-07_wp, 5.29787169100e-02_wp, 3.523118349000000e+00_wp),&
         terms(1.63340e-07_wp, 5.18826691036e+00_wp, 2.629831979980000e+01_wp),&
         terms(1.57520e-07_wp, 3.68457889430e+00_wp, 1.554203994342000e+02_wp),&
         terms(9.54100e-08_wp, 7.57422976750e-01_wp, 1.884922754997420e+04_wp),&
         terms(8.93700e-08_wp, 2.05705419118e+00_wp, 7.771377146812050e+04_wp),&
         terms(6.95200e-08_wp, 8.26733054100e-01_wp, 7.755226113240000e+02_wp),&
         terms(5.06400e-08_wp, 4.66284525271e+00_wp, 1.577343542447800e+03_wp),&
         terms(4.06100e-08_wp, 1.03057162962e+00_wp, 7.113547000800000e+00_wp),&
         terms(3.81000e-08_wp, 3.44050803490e+00_wp, 5.573142801433100e+03_wp),&
         terms(3.46300e-08_wp, 5.14074632811e+00_wp, 7.962980068163999e+02_wp),&
         terms(3.16900e-08_wp, 6.05291851171e+00_wp, 5.507553238667400e+03_wp),&
         terms(3.02000e-08_wp, 1.19246506441e+00_wp, 2.427286039740000e+02_wp),&
         terms(2.88600e-08_wp, 6.11652627155e+00_wp, 5.296909650946000e+02_wp),&
         terms(2.71400e-08_wp, 3.06378810250e-01_wp, 3.981490034082000e+02_wp),&
         terms(2.53800e-08_wp, 2.27992810679e+00_wp, 5.535694028424000e+02_wp),&
         terms(2.37100e-08_wp, 4.38118838167e+00_wp, 5.223693919802200e+03_wp),&
         terms(2.07900e-08_wp, 3.75435330484e+00_wp, 9.803210682000000e-01_wp)/)

    TYPE (terms), PARAMETER :: L4(7) = (/ &
         terms(2.89226e-06_wp, 5.84384198723e+00_wp, 6.283075849991400e+03_wp),&
         terms(3.49550e-07_wp, 0.00000000000e+00_wp, 0.000000000000000e+00_wp),&
         terms(1.68190e-07_wp, 5.48766912348e+00_wp, 1.256615169998280e+04_wp),&
         terms(2.96200e-08_wp, 5.19577265202e+00_wp, 1.554203994342000e+02_wp),&
         terms(1.28800e-08_wp, 4.72200252235e+00_wp, 3.523118349000000e+00_wp),&
         terms(7.14000e-09_wp, 5.30045809128e+00_wp, 1.884922754997420e+04_wp),&
         terms(6.35000e-09_wp, 5.96925937141e+00_wp, 2.427286039740000e+02_wp)/)

    TYPE (terms), PARAMETER :: L5(3) = (/ &
         terms(1.1408400e-06_wp, 3.14159265359e+00_wp, 0.0000000000000e+00_wp),&
         terms(7.7170000e-08_wp, 4.13446589358e+00_wp, 6.2830758499914e+03_wp),&
         terms(7.650000e-09_wp, 3.83803776214E+00_wp, 1.25661516999828e+04_wp)/)

    TYPE (terms), PARAMETER :: L6(1) = (/ &
         terms(8.78000e-09_wp, 3.14159265359e+00_wp, 0.00000000000000e+00_wp)/)

    TYPE (terms), PARAMETER :: B1(5) = (/ &
         terms(2.79620e-06_wp, 3.19870156017e+00_wp, 8.433466158130829e+04_wp),&
         terms(1.01643e-06_wp, 5.42248619256e+00_wp, 5.507553238667400e+03_wp),&
         terms(8.04450e-07_wp, 3.88013204458e+00_wp, 5.223693919802200e+03_wp),&
         terms(4.38060e-07_wp, 3.70444689758e+00_wp, 2.352866153771800e+03_wp),&
         terms(3.19330e-07_wp, 4.00026369781e+00_wp, 1.577343542447800e+03_wp)/)

    TYPE (terms), PARAMETER :: B2(2) = (/ &
         terms(9.0300000e-08_wp, 3.89729061890e+00_wp, 5.5075532386674e+03_wp),&
         terms(6.1770000e-08_wp, 1.73038850355e+00_wp, 5.2236939198022e+03_wp)/)

    TYPE (terms), PARAMETER :: R1A(20) = (/ &
         terms(1.00013988799e+00_wp, 0.000000000e+00_wp, 0.0000000000e+00_wp),&
         terms(1.670699626e-02_wp,3.09846350771e+00_wp,6.2830758499914e+03_wp),&
         terms(1.3956023e-04_wp,3.05524609620e+00_wp, 1.25661516999828e+04_wp),&
         terms(3.083720e-05_wp, 5.19846674381e+00_wp, 7.77137714681205e+04_wp),&
         terms(1.628461e-05_wp, 1.17387749012e+00_wp, 5.75338488489680e+03_wp),&
         terms(1.575568e-05_wp, 2.84685245825e+00_wp, 7.86041939243920e+03_wp),&
         terms(9.24799e-06_wp, 5.45292234084e+00_wp, 1.150676976979360e+04_wp),&
         terms(5.42444e-06_wp, 4.56409149777e+00_wp, 3.930209696219600e+03_wp),&
         terms(4.72110e-06_wp, 3.66100022149e+00_wp, 5.884926846583200e+03_wp),&
         terms(3.45983e-06_wp, 9.63686176870e-01_wp, 5.507553238667400e+03_wp),&
         terms(3.28780e-06_wp, 5.89983646482e+00_wp, 5.223693919802200e+03_wp),&
         terms(3.06784e-06_wp, 2.98671395120e-01_wp, 5.573142801433100e+03_wp),&
         terms(2.43189e-06_wp, 4.27349536153e+00_wp, 1.179062908865880e+04_wp),&
         terms(2.11829e-06_wp, 5.84714540314e+00_wp, 1.577343542447800e+03_wp),&
         terms(1.85752e-06_wp, 5.02194447178e+00_wp, 1.097707880469900e+04_wp),&
         terms(1.74844e-06_wp, 3.01193636534e+00_wp, 1.884922754997420e+04_wp),&
         terms(1.09835e-06_wp, 5.05510636285e+00_wp, 5.486777843175000e+03_wp),&
         terms(9.83160e-07_wp, 8.86813112770e-01_wp, 6.069776754553400e+03_wp),&
         terms(8.64990e-07_wp, 5.68959778254e+00_wp, 1.572083878487840e+04_wp),&
         terms(8.58250e-07_wp, 1.27083733351e+00_wp, 1.610006857376741e+05_wp)/)

    TYPE (terms), PARAMETER :: R1B(20) = (/ &
         terms(6.49030e-07_wp, 2.72506137870e-01_wp, 1.726015465469040e+04_wp),&
         terms(6.29160e-07_wp, 9.21771088320e-01_wp, 5.296909650946000e+02_wp),&
         terms(5.70560e-07_wp, 2.01374292014e+00_wp, 8.399684731811189e+04_wp),&
         terms(5.57360e-07_wp, 5.24159798933e+00_wp, 7.143069561812909e+04_wp),&
         terms(4.93840e-07_wp, 3.24501240359e+00_wp, 2.544314419883400e+03_wp),&
         terms(4.69630e-07_wp, 2.57805070386e+00_wp, 7.755226113240000e+02_wp),&
         terms(4.46610e-07_wp, 5.53715807302e+00_wp, 9.437762934887000e+03_wp),&
         terms(4.25150e-07_wp, 6.01110242003e+00_wp, 6.275962302990600e+03_wp),&
         terms(3.89680e-07_wp, 5.36071738169e+00_wp, 4.694002954707600e+03_wp),&
         terms(3.82450e-07_wp, 2.39255343974e+00_wp, 8.827390269874801e+03_wp),&
         terms(3.74900e-07_wp, 8.29529223320e-01_wp, 1.965104848109800e+04_wp),&
         terms(3.69570e-07_wp, 4.90107591914e+00_wp, 1.213955350910680e+04_wp),&
         terms(3.56600e-07_wp, 1.67468058995e+00_wp, 1.203646073488820e+04_wp),&
         terms(3.45370e-07_wp, 1.84270693282e+00_wp, 2.942463423291600e+03_wp),&
         terms(3.31930e-07_wp, 2.43703000980e-01_wp, 7.084896781115200e+03_wp),&
         terms(3.19210e-07_wp, 1.83682297810e-01_wp, 5.088628839766800e+03_wp),&
         terms(3.18460e-07_wp, 1.77775642085e+00_wp, 3.981490034082000e+02_wp),&
         terms(2.84640e-07_wp, 1.21344868176e+00_wp, 6.286598968340400e+03_wp),&
         terms(2.77930e-07_wp, 1.89934330904e+00_wp, 6.279552731642400e+03_wp),&
         terms(2.62750e-07_wp, 4.58896850401e+00_wp, 1.044738783960440e+04_wp)/)

    TYPE (terms), PARAMETER :: R1(40) = (/R1A,R1B/)

    TYPE (terms), PARAMETER :: R2(10) = (/ &
         terms(1.03018608e-03_wp,1.10748969588e+00_wp, 6.2830758499914e+03_wp),&
         terms(1.721238e-05_wp, 1.06442301418e+00_wp, 1.25661516999828e+04_wp),&
         terms(7.022150e-06_wp, 3.14159265359e+00_wp, 0.00000000000000e+00_wp),&
         terms(3.234600e-07_wp, 1.02169059149e+00_wp, 1.88492275499742e+04_wp),&
         terms(3.079900e-07_wp, 2.84353804832e+00_wp, 5.50755323866740e+03_wp),&
         terms(2.497100e-07_wp, 1.31906709482e+00_wp, 5.22369391980220e+03_wp),&
         terms(1.848500e-07_wp, 1.42429748614e+00_wp, 1.57734354244780e+03_wp),&
         terms(1.007800e-07_wp, 5.91378194648e+00_wp, 1.09770788046990e+04_wp),&
         terms(8.654000e-08_wp, 1.42046854427e+00_wp, 6.27596230299060e+03_wp),&
         terms(8.634000e-08_wp, 2.71461506020e-01_wp, 5.48677784317500e+03_wp)/)
    TYPE (terms), PARAMETER :: R3(6) = (/ &
         terms(4.359385e-05_wp,5.78455133738e+00_wp, 6.283075849991400e+03_wp),&
         terms(1.23633e-06_wp, 5.57934722157e+00_wp, 1.256615169998280e+04_wp),&
         terms(1.23410e-07_wp, 3.14159265359e+00_wp, 0.000000000000000e+00_wp),&
         terms(8.79200e-08_wp, 3.62777733395e+00_wp, 7.771377146812050e+04_wp),&
         terms(5.68900e-08_wp, 1.86958905084e+00_wp, 5.573142801433100e+03_wp),&
         terms(3.30100e-08_wp, 5.47027913302e+00_wp, 1.884922754997420e+04_wp)/)
    TYPE (terms), PARAMETER :: R4(2) = (/ &
         terms(1.44595e-06_wp, 4.27319435148e+00_wp, 6.283075849991400e+03_wp),&
         terms(6.72900e-08_wp, 3.91697608662e+00_wp, 1.256615169998280e+04_wp)/)
    TYPE (terms), PARAMETER :: R5(1) = (/ &
         terms(3.85800e-08_wp, 2.56384387339e+00_wp, 6.283075849991400e+03_wp)/)

    REAL(wp) :: l, b, r !< longitude [rad], latitude [rad] and radius [AU]
    REAL(wp) :: ld
    !
    ! --- Initialize and then calculate longitude, l, latitude, b and radius a
    !
    l = sum_vsop87 (t, L1, 1)
    l = l+sum_vsop87 (t, L2, 2)
    l = l+sum_vsop87 (t, L3, 3)
    l = l+sum_vsop87 (t, L4, 4)
    l = l+sum_vsop87 (t, L5, 5)
    l = l+sum_vsop87 (t, L6, 6)
    l = MODULO(l, twopi)
    !
    b = sum_vsop87 (t, B1, 1)
    b = b+sum_vsop87 (t, B2, 2)
    !
    r = sum_vsop87 (t, R1, 1)
    r = r+sum_vsop87 (t, R2, 2)
    r = r+sum_vsop87 (t, R3, 3)
    r = r+sum_vsop87 (t, R4, 4)
    r = r+sum_vsop87 (t, R5, 5)
    !
    ! --- Convert from Dynamic to FK5 equator & ecliptic (need extra factor 0.1
    !     for t in the equation for ld compared to Meeus original development). 
    !
    ld = l-0.1_wp*t*(1.397_wp+0.000031_wp*t)*deg2rad;
    l  = l+(-0.09033_wp+0.03916_wp*TAN(b)*(COS(ld)+SIN(ld)))*sec2rad;
    b  = b+0.03916_wp*(COS(ld)-SIN(ld))*sec2rad;
    !
    earth_position = (/l,b,r/)

  END FUNCTION earth_position
  !-----------------------------------------------------------------------------
  !>
  !! @brief Evaluates Polynomials for Coefficients given by VSOP87
  !
  PURE FUNCTION sum_vsop87 (t, term, order)

    REAL(wp) :: sum_vsop87

    REAL(wp), INTENT(in) :: t
    INTEGER , INTENT(in) :: order
    TYPE (terms), INTENT(in) :: term(:)

    REAL(wp) :: td10, tn
    INTEGER :: k

    td10       = t*0.1_wp ! convert time to decades since J2000.0
    sum_vsop87 = 0.0_wp
    IF (.NOT.(order == 1 .AND. td10 == 0.0_wp)) THEN
      tn = td10**(order-1)
    ELSE
      tn = 1.0_wp
    ENDIF

    DO k = 1, UBOUND(term,1)
      sum_vsop87 = sum_vsop87+term(k)%A*COS(term(k)%B+term(k)%C*td10)
    END DO
    sum_vsop87   = sum_vsop87*tn

  END FUNCTION sum_vsop87
  !-----------------------------------------------------------------------------
  !>
  !! @brief Returns declination calculated in last orbit call, or default value
  !
  SUBROUTINE inquire_declination(xdec)

    REAL(wp), INTENT(out) :: xdec !< declination of the sun

    IF (.NOT.initialized .AND. decl_warning) THEN
      CALL warning    ('mo_radiation_radiation/inquire_declination',                                  &
           &           'The declination of the sun has not been initialized --> use default value')
      CALL print_value('default declination',declination)
      decl_warning = .FALSE.
    END IF

    xdec = declination

  END SUBROUTINE inquire_declination


  !-----------------------------------------------------------------------------
  !>
  !! @brief Returns orbit time
  !
  SUBROUTINE get_orbit_times( l_orbvsop87, current_datetime,  &
                           & lyr_perp, yr_perp,  &
                           & time_of_day, orbit_date    )

    LOGICAL, INTENT(IN) :: l_orbvsop87
    TYPE(datetime), POINTER, INTENT(IN) :: current_datetime
    LOGICAL, INTENT(in) :: lyr_perp
    INTEGER, INTENT(in) :: yr_perp
    REAL (wp), INTENT (OUT) :: time_of_day, orbit_date

    TYPE(julianday), POINTER :: jd, jd_pal
    TYPE(datetime), POINTER  :: valid_datetime, pal_datetime
    REAL(wp)    :: zdy, zdy_mar0, zscr, zyearlen_days

       valid_datetime => newDateTime(current_datetime)
       IF (lyr_perp) THEN
         valid_datetime%date%year = yr_perp
       END IF

    IF (l_orbvsop87) THEN
      jd => newJulianday(0_i8, 0_i8)
      CALL getJulianDayFromDatetime(valid_datetime, jd) 
      orbit_date = REAL(jd%day,wp) + REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp)
      time_of_day = (REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp)-0.5_wp)*twopi
      CALL deallocateJulianday(jd)
      CALL deallocateDateTime(valid_datetime)
    ELSE
      ! Kepler orbit that needs fraction with respect to vernal equinox in rad
      ! mtime routines should be able to handle different calendars...
      jd => newJulianday(0_i8, 0_i8)
      CALL getJulianDayFromDatetime(valid_datetime, jd)
      time_of_day = (REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp)-0.5_wp)*twopi
      jd_pal => newJulianday(0_i8, 0_i8)
      pal_datetime => newDateTime('1900-01-01T00:00:00Z')
      CALL getJulianDayFromDatetime(pal_datetime, jd_pal)
      zyearlen_days=REAL(getNoOfDaysInYearDateTime(valid_datetime),wp)
      zdy = (REAL(jd%day,wp)+REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp))-       &
          &(REAL(jd_pal%day,wp)+REAL(jd_pal%ms,wp)/REAL(no_of_ms_in_a_day,wp))
      !
      ! Here is we convert to days since vernal equinox
      ! --------------------------------
      zdy_mar0   = 78.41_wp - 0.0078_wp*(1900-1987) + 0.25_wp*MOD(1900,4)
      zscr       = zdy + zyearlen_days - zdy_mar0
      orbit_date = MOD(zscr/zyearlen_days,1.0_wp)*twopi
      CALL deallocateJulianday(jd)
      CALL deallocateJulianday(jd_pal)
      CALL deallocateDateTime(valid_datetime)
    END IF
      
  END SUBROUTINE get_orbit_times

END MODULE mo_orbit
