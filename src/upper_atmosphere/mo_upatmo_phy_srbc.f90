!
! Solar heating due to absorption in Schumann-Runge continuum (SRC) and
! Schumann-Runge bands (SRB)
!
! This module allows to compute solar heating due to absorption in
! schumann-runge continuum (src) and schumann-runge bands (srb) by the method
! of strobel (jrg, vol 83, p 6225, 1978) with taking into account the heating
! efficiency for src from mlynczack&solomon (jgr, vol 98, p 10517, 1993)
! (efficiency for srb is unit).
!
! Modified by Guidi Zhou, MPI-M (2016-02-09)
! - adapted for ICON
! - reorganized code structure for efficiency and readability
! Bug fix by Guidi Zhou, MPI-M, 2016-03-16
! - pressure depends on horizontal grid, so does the efficiency factor
! Modification by Guidi Zhou, MPI-M (2016-04-06)
! - enabled using height-dependent specific heat and molecular mass of air
! Modification by Guidi Zhou, MPI-M (2016-06-02)
! - make use of the new solvar_(low/high/norm) variables in mo_impl_constants
! Modification by Guidi Zhou, MPI-M (2017-03-03)
! - added the ability to compute SRBC heating only above a certain altitude for performance
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

MODULE mo_upatmo_phy_srbc

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: SUCCESS
  USE mo_upatmo_impl_const,    ONLY: isolvar, isolvardat

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: srbc_heating  ! subroutine to compute SRB and SRC heating

  ! for error handling
  INTEGER, PARAMETER :: IERR_NO     = SUCCESS      ! = 0
  INTEGER, PARAMETER :: IERR_SOLVAR = SUCCESS + 1

CONTAINS

  !>
  !! compute SRB and SRC heating
  !! subroutine o2strob in original source of V. Fomichev
  !!
  SUBROUTINE srbc_heating(jcs, jce, kbdim, klev, ppf, prmu0, am, cp, zo2, tto2, heato2, &
    &                     solvar_type, solvar_data, opt_sunlit_idx, opt_nsunlit,        &
    &                     opt_istartlev, opt_iendlev, opt_error)

    ! IN/OUT
    INTEGER , INTENT(IN)  :: jcs, jce, kbdim, klev ! dimensions
    REAL(wp), INTENT(IN)  :: ppf(kbdim, klev)    ! full level pressure       [Pa]
    REAL(wp), INTENT(IN)  :: prmu0(kbdim)        ! cos of solar zenith angle []
    REAL(wp), INTENT(IN)  :: am(kbdim, klev)     ! molecular mass of air     [g]
    REAL(wp), INTENT(IN)  :: cp(kbdim, klev)     ! specific heat of air      [J/K/kg]
    REAL(wp), INTENT(IN)  :: zo2(kbdim,klev)     ! o2 vmr                    [m3/m3]
    REAL(wp), INTENT(IN)  :: tto2(kbdim,klev)    ! o2 column density         [molecules/cm2]
    REAL(wp), INTENT(OUT) :: heato2(kbdim,klev)  ! tendency dT/dt (K/s)
    INTEGER,  INTENT(IN)  :: solvar_type         ! solar activity
    INTEGER,  INTENT(IN)  :: solvar_data         ! solar activity data type
    INTEGER,  OPTIONAL, TARGET, INTENT(IN) :: opt_sunlit_idx(:)   ! optional list with indices of sunlit grid columns
    INTEGER,  OPTIONAL, INTENT(IN)  :: opt_nsunlit                ! optional number of sunlit grid columns
    INTEGER,  OPTIONAL, INTENT(IN)  :: opt_istartlev, opt_iendlev ! optional vertical start and end indices
    INTEGER,  OPTIONAL, INTENT(OUT) :: opt_error                  ! for optional error handling

    ! LOCAL
    INTEGER  :: jl, jk, istartlev, iendlev
    REAL(wp) :: n2, src1, src2, src, srb
    REAL(wp) :: f_svar(4)          ! coefficients for solar activity conditions
    REAL(wp) :: effsrc             ! efficiency factor
    REAL(wp) :: x
    REAL(wp) :: inv_prmu0(kbdim)

    INTEGER, ALLOCATABLE, TARGET :: sunlit_idx(:)
    INTEGER,             POINTER :: idxlist(:)
    INTEGER  :: nsunlit, jsunlit

    LOGICAL  :: l_present_error, l_present_nsunlit

    ! coefficients for the Chebyshev polynomial fit for the efficiency factor
    REAL(wp), PARAMETER :: cho2(4) = [0.75349_wp, 0.0036_wp, 0.059468_wp, -0.022795_wp]
    REAL(wp), PARAMETER :: prmu0_min = 1.e-10_wp

    !---------------------------------------------------------

    ! please do not limit range of assignment 
    ! (e.g., heato2(jcs:jce,istartlev:iendlev) = 0._wp)), 
    ! since tendencies have attribute INTENT(OUT)
    heato2(:,:) = 0._wp

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
      IF (nsunlit < 1) RETURN  ! all grid cell columns are dark
      idxlist => sunlit_idx
    ENDIF
    IF (nsunlit < 1) RETURN  ! all grid cell columns are dark

    ! solar activity
    IF (solvar_type == isolvar%low .AND. solvar_data == isolvardat%rottman) THEN
      ! G. Rottman data for solar low
      f_svar = [1.378_wp, 1.108_wp, 0.951_wp, 0.97_wp]
    ELSEIF (solvar_type == isolvar%low .AND. solvar_data == isolvardat%lean) THEN
      ! J. Lean data for solar low
      f_svar = [1.384_wp, 1.145_wp, 0.956_wp, 0.949_wp]
    ELSEIF (solvar_type == isolvar%high .AND. solvar_data == isolvardat%rottman) THEN
      ! G. Rottman data for solar high
      f_svar = [1.666_wp, 1.277_wp, 1.045_wp, 1.03_wp]
    ELSEIF (solvar_type == isolvar%high .AND. solvar_data == isolvardat%lean) THEN
      ! J. Lean data for solar high
      f_svar = [1.717_wp, 1.318_wp, 1.068_wp, 1.021_wp]
    ELSEIF (solvar_type == isolvar%norm) THEN  ! normal conditions
      ! data for normal conditions
      f_svar = 1.0_wp
    ELSE
      ! no valid solar activity type
      IF (l_present_error) opt_error = IERR_SOLVAR
      RETURN
    ENDIF

    ! precompute inverse of cosine of solar zenith angle
    DO jsunlit = 1, nsunlit
      jl = idxlist(jsunlit)
      inv_prmu0(jl) = 1._wp / MAX(prmu0_min, prmu0(jl))
    ENDDO  !jsunlit

    DO jk = istartlev, iendlev
      DO jsunlit = 1, nsunlit
        jl = idxlist(jsunlit)
        ! efficiency factor for SRC
        IF(ppf(jl, jk) > 1.0_wp) THEN
          effsrc = 0.7938_wp
        ELSE IF(ppf(jl, jk) > 1.e-2_wp) THEN
          x = LOG10(ppf(jl, jk)) + 1._wp
          effsrc = cho2(1) + x * (cho2(2) + x * (cho2(3) + cho2(4) * x))
        ELSE
          effsrc = 0.8320_wp
        END IF
        
        ! O2 amount: downward looking paths
        n2 = tto2(jl, jk) * inv_prmu0(jl)
        
        !=================================================================
        ! Coefficients according to the original Strobel 1978 paper
        !=================================================================
        ! SRC:
        ! src1 =  1.1e-7_wp    * f_svar(1)                           * EXP(-1.e-17_wp   * n2)
        ! src2 = ( 3.267e-7_wp * f_svar(3)                           * EXP(-2.9e-19_wp  * n2) + &
        !   &     (1.433e-7_wp * f_svar(2) - 3.267e-7_wp * f_svar(3))* EXP(-1.7e-18_wp  * n2) - &
        !   &      1.433e-7_wp * f_svar(2)                           * EXP(-1.15e-17_wp * n2)   &
        !   &    ) / n2
        ! src = (src1 + src2) * effsrc
        
        ! ! SRB:
        ! IF(n2 > 1.e18_wp) THEN
        !    srb = 1._wp / (0.67e7_wp * n2 + 3.44e16_wp * SQRT(n2)) * f_svar(4)
        ! ELSE
        !    srb = 2.43e-26 * f_svar(4)
        ! ENDIF
        
        !=================================================================
        ! Coefficients used in HAMMONIA
        !=================================================================
        ! SRC:
        src1 =  2.716e6_wp  * f_svar(1)                * EXP(-1.e-17_wp   * n2)
        src2 = 5.902e23_wp * (f_svar(3)                * EXP(-2.9e-19_wp  * n2) + &
          &    (0.43883429_wp * f_svar(2) - f_svar(3)) * EXP(-1.7e-18_wp  * n2) - &
          &     0.43883429_wp * f_svar(2)              * EXP(-1.15e-17_wp * n2)   &
          &    ) / n2
        src = (src1 + src2) * effsrc
        
        ! SRB:
        IF(n2 > 1.e18_wp) THEN
          srb = 1._wp / (1.113e-24_wp * n2 + 5.712e-15_wp * SQRT(n2)) * f_svar(4)
        ELSE
          srb = 1.463e5 * f_svar(4)
        ENDIF

        ! O2 heating (K/s):
        heato2(jl, jk) = zo2(jl, jk) * (src + srb) / ( 10000._wp * am(jl, jk) * cp(jl, jk) )
      ENDDO  !jsunlit
    ENDDO  !jk

    ! clean-up
    idxlist => NULL()
    IF (ALLOCATED(sunlit_idx)) DEALLOCATE(sunlit_idx)

  END SUBROUTINE srbc_heating

END MODULE mo_upatmo_phy_srbc
