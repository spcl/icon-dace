!#define G3X

!
! Computation of EUV heating
!
! ATTENTION
! This version includes the possibility to compute EUV heating for diffderent solar activity conditions.
! The F10.7 index is taken from the canadian website
! http://www.drao-ofr.hia-iha.nrc-cnrc.gc.ca/icarus/www/sol_home.shtml
! For solar high the value of 235.1 (Nov. 1989) and for solar low the
! value of 68.7 (Sep. 86) is taken. Judith Lean's solar variability data for UV
! and visible are taken at these dates.
! The third possibility is to use the mean value for the months of Jan to Jun 1990 (187.3)
! which were used in the simulations of the 27day cycle.
! If the Rottman data were used, the following values should be taken:
! solar high: 192.8 (29.3.1992); solar low: 77.3 (2.1.1995)
! H. Schmidt, June 2003
!
!=============================================================================
!
! - Description:
!
! This module serves to calculate the solar heating due to absorption
! in Extreme UV region (5-105nm). Data on fluxes and absorption coefficients
! are taken from Richards et al., 1994, JGR, vol 99, No A7, 13,283
!
! The module contains:
! A) status indicators for allocation and initialization
! B) constants and variables
! C) subroutine to compute solar flux in EUV
! D) subroutine to compute solar heating in EUV
!
!-----------------------------------------------------------------------------
!
! V. Fomichev   , November, 1997: original source
! M.A. Giorgetta, MPI, June 2001: rewrite for ECHAM5
!
! G. Zhou, MPI, June 2016: rewrite for ICON
! - make use of vertically-varying gravity
! Guidi Zhou, MPI-M, 2016-09-08:
! - disabled considering upper level grid points that can still be shined when surface
!   solar zenith angle is greater than 90deg: this is not crucial but nice to have, but
!   the original code produces large errors
! - cleaned up a bit
! Guidi Zhou, MPI-M, 2016-09-13:
! - fixed a bug not accumulating the heating rate from each spectral band
! Guidi Zhou, MPI-M, 2017-03-03:
! - added the ability to compute EUV heating only above a certain altitude for performance
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

MODULE mo_upatmo_phy_euv

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: ana => avo
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_upatmo_impl_const,  ONLY: isolvar, iorbit, icycle
  USE mo_math_constants,     ONLY: pi2
  USE mo_orbit,              ONLY: orbit_kepler, orbit_vsop87, get_orbit_times
  USE mtime,                 ONLY: datetime, timedelta, newDateTime, operator(-), getTotalSecondsTimeDelta

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  euv_heating       ! subroutine to compute solar heating in EUV

  ! =====================================================================
  ! constants and variables for heting efficiency, cross sections,
  ! solar EUV flux and solar sctivity
  ! =====================================================================

  ! heating efficiency in EUV
  ! fixed at 0.366 (Roble ?Ref.?)
  ! This value includes the loss of energy due to NO cooling. In this model version, however,
  ! NO cooling is calculated explicitly. According to Richards et al. (1982, Plant Sp. Sc.;
  ! cited by Wells et al. (1997, Ann. Geophys.), the consideration of NO cooling reduces
  ! the peak efficiency of EUV heating from 0.6 to 0.45. These values are used to scale the
  ! efficiency used in this model, resulting in a value of 0.49.
! REAL(wp), PARAMETER   :: euveff= 0.366_wp ! without explicit NO cooling in the model
  REAL(wp), PARAMETER   :: euveff= 0.49_wp

  INTEGER, PARAMETER    :: neuv=37

  ! centers of the EUV bands (nm):
  REAL(wp), PARAMETER   :: al(neuv)=                                    &
       & (/ 7.5,   12.5,  17.5,  22.5,   25.63, 28.415, 27.5, 30.331,   &
       &    30.378,32.5,  36.807,37.5,   42.5,  46.522, 47.5, 52.5,     &
       &    55.437,58.433,57.5,  60.976, 62.973,62.5,   67.5, 70.336,   &
       &    72.5,  76.515,77.041,78.936, 77.5,  82.5,   87.5, 92.5,     &
       &    97.702,97.5, 102.572,103.191,102.5/)

  ! cross section of O2 (in 1.e18 cm2)
  REAL(wp), PARAMETER   :: so2(neuv)= 1.e-18_wp*                        &
       & (/  1.316, 3.806, 7.509,10.900,13.370,15.790,14.387,16.800,    &
       &    16.810,17.438,18.320,18.118,20.310,21.910,23.101,24.606,    &
       &    26.040,22.720,26.610,28.070,32.060,26.017,21.919,27.440,    &
       &    28.535,20.800,18.910,26.668,22.145,16.631, 8.562,12.817,    &
       &    18.730,21.108, 1.630, 1.050, 1.346/)

  ! cross section of N2 (in 1.e18 cm2)
  REAL(wp), PARAMETER   :: sn2(neuv)= 1.e-18_wp*                        &
       & (/  0.720, 2.261, 4.958, 8.392,10.210,10.900,10.493,11.670,    &
       &    11.700,13.857,16.910,16.395,21.675,23.160,23.471,24.501,    &
       &    24.130,22.400,22.787,22.790,23.370,23.339,31.755,26.540,    &
       &    24.662,120.49,14.180,16.487,33.578,16.992,20.249, 9.680,    &
       &     2.240,50.988, 0.0,   0.0,  0.0/)

  ! cross section of O  (in 1.e18 cm2)
  REAL(wp), PARAMETER   :: so(neuv)= 1.e-18_wp*                         &
       & (/  0.730, 1.839, 3.732, 5.202, 6.050, 7.080, 6.461, 7.680,    &
       &     7.700, 8.693, 9.840, 9.687,11.496,11.930,12.127,12.059,    &
       &    12.590,13.090,13.024,13.400,13.400,13.365,17.245,11.460,    &
       &    10.736, 4.000, 3.890, 3.749, 5.091, 3.498, 4.554, 1.315,    &
       &     0.0,   0.0,   0.0,   0.0,   0.0/)

  ! parameters for the solar EUV flux model (flux in 1.e-9 photon/cm2/s)
  REAL(wp), PARAMETER   :: f74113(neuv)=                                &
       & (/ 1.200, 0.450, 4.800, 3.100, 0.460, 0.210, 1.679, 0.800,     &
       &    6.900, 0.965, 0.650, 0.314, 0.383, 0.290, 0.285, 0.452,     &
       &    0.720, 1.270, 0.357, 0.530, 1.590, 0.342, 0.230, 0.360,     &
       &    0.141, 0.170, 0.260, 0.702, 0.758, 1.625, 3.537, 3.000,     &
       &    4.400, 1.475, 3.500, 2.100, 2.467/)
  REAL(wp), PARAMETER   :: ai(neuv)=                                    &
       & (/ 1.0017e-02, 7.1250e-03, 1.3375e-02, 1.9450e-02, 2.7750e-03, &
       &    1.3768e-01, 2.6467e-02, 2.5000e-02, 3.3333e-03, 2.2450e-02, &
       &    6.5917e-03, 3.6542e-02, 7.4083e-03, 7.4917e-03, 2.0225e-02, &
       &    8.7583e-03, 3.2667e-03, 5.1583e-03, 3.6583e-03, 1.6175e-02, &
       &    3.3250e-03, 1.1800e-02, 4.2667e-03, 3.0417e-03, 4.7500e-03, &
       &    3.8500e-03, 1.2808e-02, 3.2750e-03, 4.7667e-03, 4.8167e-03, &
       &    5.6750e-03, 4.9833e-03, 3.9417e-03, 4.4167e-03, 5.1833e-03, &
       &    5.2833e-03, 4.3750e-03/)

  ! for error handling
  INTEGER, PARAMETER :: IERR_NO     = SUCCESS      ! = 0
  INTEGER, PARAMETER :: IERR_ORBIT  = SUCCESS + 1
  INTEGER, PARAMETER :: IERR_SOLVAR = SUCCESS + 2

CONTAINS

  !>
  !! Compute solar heating in Extreme UltraViolet
  !!
  !! Literature:
  !! - Richards, P. G., Fennelly, J. A., and Torr, D. G. (1994) 
  !!   EUVAC: A solar EUV flux model for aeronomic calculations. 
  !!   J. Geophys. Res.-Space, 99, 8981-8992.
  !!
  SUBROUTINE euv_heating(jcs, jce, kbdim, klev, prmu0, zo2, zn2, zo, tto2, ttn2, ttox, amu, cp, ptte, &
    &                    this_datetime, orbit_type, solvar_type, solcyc_type, cecc, cobld, clonp,     &
    &                    lyr_perp, yr_perp, opt_sunlit_idx, opt_nsunlit, opt_istartlev, opt_iendlev, opt_error)

    ! in/out variables
    INTEGER , INTENT(IN)  :: jcs,jce, kbdim, klev ! longitude and latitude dimensions
    REAL(wp), INTENT(IN)  :: prmu0(kbdim)     ! cos of solar zenith angle
    REAL(wp), INTENT(IN)  :: zo2(kbdim,klev)  ! o2 vmr
    REAL(wp), INTENT(IN)  :: zn2(kbdim,klev)  ! n2 vmr
    REAL(wp), INTENT(IN)  :: zo(kbdim,klev)   ! o vmr
    REAL(wp), INTENT(IN)  :: tto2(kbdim,klev) ! o2 column density (1/cm2)
    REAL(wp), INTENT(IN)  :: ttn2(kbdim,klev) ! n2 column density (1/cm2)
    REAL(wp), INTENT(IN)  :: ttox(kbdim,klev) ! o column density (1/cm2)
    REAL(wp), INTENT(IN)  :: amu(kbdim,klev)  ! molecular mass of air
    REAL(wp), INTENT(IN)  :: cp(kbdim,klev)   ! specific heat
    REAL(wp), INTENT(OUT) :: ptte(kbdim,klev) ! tendency dT/dt (K/s)
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER,  INTENT(IN)  :: orbit_type       ! orbit model 
    INTEGER,  INTENT(IN)  :: solvar_type      ! solar activity
    INTEGER,  INTENT(IN)  :: solcyc_type      ! solar cycle 
    REAL(wp), INTENT(IN)  :: cecc             ! eccentricity of orbit
    REAL(wp), INTENT(IN)  :: cobld            ! obliquity of Earth axis
    REAL(wp), INTENT(IN)  :: clonp            ! longitude of perihelion
    LOGICAL,  INTENT(IN)  :: lyr_perp         ! switch for perpetual Earth orbit of year yr_perp
    INTEGER,  INTENT(IN)  :: yr_perp          ! year for which orbit is perpetuated
    INTEGER,  OPTIONAL, TARGET, INTENT(IN) :: opt_sunlit_idx(:)   ! optional list with indices of sunlit grid columns
    INTEGER,  OPTIONAL, INTENT(IN)  :: opt_nsunlit                ! optional number of sunlit grid columns
    INTEGER,  OPTIONAL, INTENT(IN)  :: opt_istartlev, opt_iendlev ! optional vertical start and end indices
    INTEGER,  OPTIONAL, INTENT(OUT) :: opt_error                  ! for optional error handling

    ! local variables
    REAL(wp) :: s, s1, e, aux, cumtte
    INTEGER  :: jk, jb, jl
    INTEGER  :: istartlev, iendlev

    REAL(wp) :: flux(neuv)
    REAL(wp) :: sqrd_prmu0(kbdim)

    REAL(wp) :: rasc_sun, decl_sun, dist_sun
    REAL(wp) :: time_of_day, orbit_date

    INTEGER, ALLOCATABLE, TARGET :: sunlit_idx(:)
    INTEGER,             POINTER :: idxlist(:)
    INTEGER  :: nsunlit, jsunlit

    LOGICAL  :: l_present_error, l_present_nsunlit

    !---------------------------------------------------------

    ! please do not limit range of assignment 
    ! (e.g., ptte(jcs:jce,istartlev:iendlev) = 0._wp)), 
    ! since tendencies have attribute INTENT(OUT)
    ptte(:, :) = 0._wp

    ! we are within openMP-threading, 
    ! so only rudimentary error handling is possible
    IF (PRESENT(opt_error)) THEN 
      opt_error       = IERR_NO
      l_present_error = .TRUE.
    ELSE
      l_present_error = .FALSE.
    ENDIF

    ! determine start and end indices of vertical grid layers,
    ! for which tendencies should be computed
    IF (PRESENT(opt_istartlev)) THEN
      istartlev = MIN(MAX(1, opt_istartlev), klev)
    ELSE
      istartlev = 1
    ENDIF

    IF (PRESENT(opt_iendlev)) THEN
      iendlev = MIN(MAX(1, opt_iendlev), klev)
    ELSE
      iendlev = klev
    ENDIF

    IF (istartlev > iendlev) RETURN

    IF (PRESENT(opt_nsunlit)) THEN
      ! no further computations are necessary, 
      ! if all grid cell columns are dark
      IF (opt_nsunlit < 1) RETURN
      nsunlit           = opt_nsunlit
      l_present_nsunlit = .TRUE.
    ELSE
      l_present_nsunlit = .FALSE.
    ENDIF

    IF (PRESENT(opt_sunlit_idx)) THEN
      IF (.NOT. l_present_nsunlit) nsunlit = SIZE(opt_sunlit_idx)
      idxlist => opt_sunlit_idx
    ELSE
      ! we determine the index list ourselves
      nsunlit = 0
      ! for convenience, we allocate the index list with kbdim 
      ! and not with nsunlit
      ALLOCATE(sunlit_idx(kbdim))
      sunlit_idx(:) = 0
      DO jl = jcs, jce
        IF(prmu0(jl) > 0._wp) THEN
          nsunlit             = nsunlit + 1
          sunlit_idx(nsunlit) = jl
        ENDIF
      ENDDO  !jl
      idxlist => sunlit_idx
    ENDIF
    IF (nsunlit < 1) RETURN  ! all grid cell columns are dark

    ! solar flux
    flux = euv_flux(this_datetime, solvar_type, solcyc_type, opt_error)

    
    IF (orbit_type == iorbit%vsop87) THEN
      ! standard orbit model
      CALL get_orbit_times(.TRUE.,this_datetime, lyr_perp, yr_perp, time_of_day, orbit_date)
      CALL orbit_vsop87(orbit_date, rasc_sun, decl_sun, dist_sun)
    ELSEIF (orbit_type == iorbit%kepler) THEN
      ! idealized Kepler orbit
      CALL get_orbit_times(.FALSE.,this_datetime, lyr_perp, yr_perp, time_of_day, orbit_date)
      CALL orbit_kepler(cecc, cobld, clonp, orbit_date, rasc_sun, decl_sun, dist_sun)      
    ELSE
      ! invalid orbit type
      IF (l_present_error) opt_error = IERR_ORBIT
    ENDIF

    ! auxiliary factor
    aux = euveff * ana / ( 10000._wp * dist_sun**2 )

    ! precompute square of cosine of solar zenith angle
    DO jsunlit = 1, nsunlit
      jl = idxlist(jsunlit)
      sqrd_prmu0(jl) = prmu0(jl)**2
    ENDDO  !jsunlit

    DO jk = istartlev, iendlev
      
      DO jsunlit = 1, nsunlit
        jl = idxlist(jsunlit)
        
        cumtte = 0._wp

        DO jb = 1, neuv  ! loop over spectral bands
          
          s  = so2(jb) * tto2(jl, jk) + sn2(jb) * ttn2(jl, jk) + so(jb) * ttox(jl, jk)
          e  = -s * 35._wp / SQRT(1224.0_wp * sqrd_prmu0(jl) + 1._wp)
          s1 = flux(jb) * ( so2(jb) * zo2(jl,jk) + sn2(jb) * zn2(jl,jk) + so(jb) * zo(jl,jk) )
          cumtte = cumtte + s1 * EXP(e)
          
        END DO  !jb
        
        ! temperature tendency
        ptte(jl, jk) = aux * cumtte / ( amu(jl,jk) * cp(jl,jk) )
        
      END DO  !jsunlit
      
    END DO  !jk

    ! clean-up
    idxlist => NULL()
    IF (ALLOCATED(sunlit_idx)) DEALLOCATE(sunlit_idx)

  END SUBROUTINE euv_heating

  !>
  !! Compute solar flux in Extreme UltraViolet
  !!
  FUNCTION euv_flux(this_datetime, solvar_type, solcyc_type, opt_error) RESULT(flux)

    ! in/out variables
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    INTEGER, INTENT(IN) :: solvar_type
    INTEGER, INTENT(IN) :: solcyc_type
    INTEGER, OPTIONAL, INTENT(INOUT) :: opt_error
    REAL(wp) :: flux(neuv)

    ! local variables
    REAL(wp) :: p
    REAL(wp) :: f107, f107a, sin_27d
    INTEGER  :: jb
 
    ! The factor for the solar 27 day variability is the result of following assumption:
    ! Amplitude of the 27 day cycle is 41.8 sfu. This was evaluated by fitting a sinus
    ! for a period of 27 days to F10.7 data from the first 6 months of 1990.
    ! 11.16% of 187.3 sfu (mean value for these 6 months) gives half of this amplitude.
    REAL(wp), PARAMETER  :: fact_27d = .1116_wp
    REAL(wp), PARAMETER  :: fact_aux = 1._wp / 86400._wp / 27.0_wp * pi2

    TYPE(datetime),  TARGET :: date900101
    TYPE(timedelta), TARGET :: dt

    !---------------------------------------------------------

    IF (solvar_type == isolvar%low) THEN       ! solar low
      f107  = 68.7_wp
      f107a = 68.7_wp
    ELSEIF (solvar_type == isolvar%high) THEN  ! solar high
      f107  = 235.1_wp 
      f107a = 235.1_wp
    ELSEIF (solvar_type == isolvar%norm) THEN  ! normal conditions
      f107  = 150._wp  ! original 
      f107a = 150._wp
    ELSE
      ! no valid solar activity type
      IF (PRESENT(opt_error)) opt_error = IERR_SOLVAR
    ENDIF

    IF (solcyc_type == icycle%day27) THEN
       date900101 = newDateTime(1990, 1, 1, 0, 0, 0, 0)
       dt = this_datetime - date900101

       sin_27d = SIN(fact_aux * REAL(getTotalSecondsTimeDelta(dt, date900101), wp))

       f107  = f107 * (1._wp + sin_27d * fact_27d)
    ENDIF

    ! Computation in detail:
    !-----------------------

    ! ! to calculate the flux for given F107 and F107A:
    ! p = 0.5_wp * (f107 + f107a)
    ! ! fluxes in photon/cm2/s for given F107 and F107A:
    ! DO jb = 1, neuv
    !   photon(jb) = f74113(jb) * (1._wp + ai(jb)*(p-80._wp)) * 1.e9_wp
    ! END DO
    ! ! fluxes in erg/cm2/s, flux = h*c/al * photon (1nm = 1.e-7 cm):
    ! DO jb = 1, neuv
    !   flux(jb) = 1.98648e-09_wp / al(jb) * photon(jb)
    ! ENDDO

    ! Computation compressed:
    !------------------------

    p = 1.98648_wp * (0.5_wp * (f107 + f107a) - 80._wp)
    
    DO jb = 1, neuv
      flux(jb) = f74113(jb) * (1.98648_wp + p * ai(jb)) / al(jb)
    ENDDO

  END FUNCTION euv_flux

END MODULE mo_upatmo_phy_euv
