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

! Computes ion drag that is exerted on the horizontal wind fields u and v above 100 km of altitude
!
! This routine computes the physical tendencies of the horizontal
! wind fields and temperature due to ion drag above 100 km of altitude.
! The tendencies for winds are obtained from a semi-implicit time-stepping
! procedure that is precise to order (dtime)^2.
!
!      f(t+dtime)-f(t-dtime)    - M * ( f(t+dtime) + f(t-dtime) )
!      --------------------- =        ---------------------------
!             2*dtime                              2
!
! f is is the 2D wind vector (u,v) and M is a 2x2 matrix representing the
! drag coefficients.
! The values for the matrix M are calculated using the method of
! Hong and Lindzen, 1976: JAS, 33, 135-153.  for minimum solar forcing.
! The Lorenz terms (zldrag) are from m.charron and mimic those of Hong and
! Lindzen (fig. 4).
!
!           _                     _
!          |                       |
!          | zdrag        zldrag   |
!      M = |                       |
!          | -zldrag   zcoef*zdrag |
!          |_                     _|
!           _   _
!          |     |
!          |  u  |
!      f = |     |
!          |  v  |
!          |_   _|
!
! INPUT ARGUMENTS:
! ---------------
!
!  pum1     : zonal wind (t-dt)
!  pvm1     : meridional wind (t-dt)
!  pqm1     : humidity (t-dt)
!  pgeom1   : geopotential above surface (t-dt)
!
! INPUT/OUTPUT ARGUMENTS:
! ----------------------
!
!  pvol     : tendency of meridional wind
!  pvom     : tendency of zonal wind
!  ptte     : tendency of temperature
!
! Modification by G. Zhou - MPI - 20160608
! - rewrite for ICON
! Modification by G. Zhou - MPI - 20160620
! - make use of vertically-varying gravity

MODULE mo_upatmo_phy_iondrag

  USE mo_kind,                 ONLY: wp
  USE mo_physical_constants,   ONLY: vtmpc2
  USE mo_impl_constants,       ONLY: SUCCESS
  USE mo_upatmo_impl_const,    ONLY: isolvar
  USE mo_turb_vdiff_params,    ONLY: cvdifts

  IMPLICIT NONE

  PUBLIC :: iondrag

  ! for error handling
  INTEGER, PARAMETER :: IERR_NO     = SUCCESS      ! = 0
  INTEGER, PARAMETER :: IERR_SOLVAR = SUCCESS + 1

CONTAINS

  !>
  !! Compute ion drag and corresponding Joule heating
  !!
  !! Literature:
  !! - Hong, S.-S. and Lindzen, R. S. (1976) Solar semidiurnal tide in the thermosphere. 
  !!   J. Atmos. Sci., 33, 135-153.
  !!
  SUBROUTINE iondrag(jcs, jce, kbdim, klev, solvar_type, psteplen, lat, pum1, pvm1, pqm1, &
    &                grav, pgeom1, pcp, pvom, pvol, ptte, opt_istartlev, opt_iendlev, opt_error)

    !----- SUBROUTINE ARGUMENTS -----
    
    INTEGER,  INTENT(IN)  :: jcs, jce, kbdim, klev
    INTEGER,  INTENT(IN)  :: solvar_type             ! identifier for solar activity
    REAL(wp), INTENT(IN)  :: psteplen                ! time step length, usually 2*delta_time
    REAL(wp), INTENT(IN)  :: lat(kbdim)              ! latitude              [rad]
    REAL(wp), INTENT(IN)  :: pum1(kbdim,klev),    &  ! u-velocity            [u]
      &                      pvm1(kbdim,klev),    &  ! v-velocity            [v]
      &                      pqm1(kbdim,klev),    &  ! specific humidity     [q]
      &                      grav(kbdim,klev),    &  ! gravity
      &                      pgeom1(kbdim,klev),  &  ! geopotential
      &                      pcp(kbdim,klev)         ! specific heat of air  [cp]
    REAL(wp), INTENT(OUT) :: pvom(kbdim,klev),    &  ! u-velocity tendency   [du/dt]
      &                      pvol(kbdim,klev),    &  ! v-velocity tendency   [dv/dt]
      &                      ptte(kbdim,klev)        ! temperature tendency  [dT/dt]
    INTEGER, OPTIONAL, INTENT(IN)  :: opt_istartlev, opt_iendlev  ! optional vertical start and end indices
    INTEGER, OPTIONAL, INTENT(OUT) :: opt_error      ! for optional error handling
    
    !----- INTERNAL VARIABLES -----
    
    INTEGER                    :: jl, jk, istartlev, iendlev
    REAL(wp)                   :: zlat, zcons1, zcons2, zcons4(3), ztmp(3)
    REAL(wp)                   :: zalpha, zcons5, ztmst, zup1, zvp1
    REAL(wp), DIMENSION(kbdim) :: zcoef, zcoef1
    REAL(wp)                   :: zdrag, zldrag, zdenum
    REAL(wp), DIMENSION(3)     :: za,zb,zc,zd,ze
    REAL(wp)                   :: zal, zdgeom, zaux0, zaux1, zaux2, zaux3, zdum1, zdvm1

    LOGICAL  :: l_present_error
    
    !----- TABLE FOR ION DRAG COEFFICIENTS TAKEN FROM -----
    !_____   HONG AND LINDZEN, 1976: JAS, 33, p. 152  -----
    
    REAL(wp), PARAMETER, DIMENSION(3) :: zamin = (/  6.6E4_wp  , 1.56E5_wp ,  3.0E5_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zamax = (/ 1.15E5_wp  , 2.75E5_wp , 1.05E6_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zbmin = (/    1.4_wp  ,    1.0_wp ,   0.35_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zbmax = (/    1.4_wp  ,    1.0_wp ,   0.20_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zcmin = (/ 150.E3_wp  , 225.E3_wp , 275.E3_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zcmax = (/ 150.E3_wp  , 240.E3_wp , 300.E3_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zdmin = (/    0.2_wp  ,    0.0_wp ,    0.1_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zdmax = (/    0.2_wp  ,    0.0_wp ,    0.1_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zemin = (/   1.E3_wp  ,  42.E3_wp ,   1.E3_wp  /)
    REAL(wp), PARAMETER, DIMENSION(3) :: zemax = (/   1.E3_wp  ,  52.E3_wp ,   1.E3_wp  /)
    REAL(wp), PARAMETER               :: zcons3 = 5.0E-10_wp
    REAL(wp), PARAMETER               :: zalmin = 37.3_wp
    REAL(wp), PARAMETER               :: zalmax = 36.7_wp
    
    !---------------------------------------------------------
    
    !----- SET INITIAL VALUES -----

    ! please do not limit range of assignment 
    ! (e.g., ptte(jcs:jce,istartlev:iendlev) = 0._wp)), 
    ! since tendencies have attribute INTENT(OUT)
    pvom(:,:) = 0._wp
    pvol(:,:) = 0._wp
    ptte(:,:) = 0._wp

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
    
    IF (solvar_type == isolvar%low) THEN       ! solar low
      za  = zamin
      zb  = zbmin
      zc  = zcmin
      zd  = zdmin
      ze  = zemin
      zal = zalmin
    ELSEIF (solvar_type == isolvar%high) THEN  ! solar high
      za  = zamax
      zb  = zbmax
      zc  = zcmax
      zd  = zdmax
      ze  = zemax
      zal = zalmax
    ELSEIF (solvar_type == isolvar%norm) THEN  ! normal conditions
      za  = 0.5_wp * ( zamin + zamax )
      zb  = 0.5_wp * ( zbmin + zbmax )
      zc  = 0.5_wp * ( zcmin + zcmax )
      zd  = 0.5_wp * ( zdmin + zdmax )
      ze  = 0.5_wp * ( zemin + zemax )
      zal = 0.5_wp * ( zalmin + zalmax )
    ELSE
      ! no valid solar activity type
      IF (l_present_error) opt_error = IERR_SOLVAR
      RETURN
    ENDIF
    
    ztmst  = psteplen
    zalpha = cvdifts
    zcons5 = ztmst * zalpha
    zdrag  = 0._wp
    zldrag = 0._wp

    DO jl= jcs, jce
      zlat       = lat(jl)
      zcons1     = -2._wp * TAN(zlat)
      zcons2     = ATAN(zcons1)
      zcoef1(jl) = SIN(zcons2)
      zcoef(jl)  = zcoef1(jl) * zcoef1(jl)
    ENDDO  !jl

    DO jk = istartlev, iendlev
      DO jl = jcs, jce

        zaux0     = pgeom1(jl,jk) / grav(jl,jk)
        ztmp(:)   = ( zaux0 - zc(:) ) / ( zd(:) * zaux0 + ze(:) )
        zcons4(:) = 1._wp - ztmp(:) - EXP( -ztmp(:) )
        zdrag     = zcons3 * SUM( za(:) * EXP( zb(:) * zcons4(:) ) )
        zdgeom    = pgeom1(jl,jk) - 55.E3_wp * grav(jl,jk)
        zldrag    = zcoef1(jl) * EXP( 1.9E-17_wp * ( 3.5_wp * zdgeom )**3 / &
          &         ( EXP ( 3.5E-6_wp * zdgeom ) - 1._wp ) - zal )
        zaux1     = zcons5 * zdrag
        zaux2     = zcons5 * zldrag
        zaux3     = 1._wp + zaux1 * zcoef(jl)
        zdenum    = 1._wp / ( zaux3 * ( 1._wp + zaux1 ) + zaux2**2 )

        pvom(jl,jk) = zdenum * ( -zdrag * zaux3 * pum1(jl,jk) &
          &         - zldrag * ( pvm1(jl,jk) + zaux2 * pum1(jl,jk) ) )
        pvol(jl,jk) = zdenum * ( -zcoef(jl) * zdrag * ( 1._wp + zaux1 ) * pvm1(jl,jk) &
          &         + zldrag * ( pum1(jl,jk) - zaux2 * pvm1(jl,jk) ) )
        zdum1       = ztmst * pvom(jl,jk)
        zdvm1       = ztmst * pvol(jl,jk)
        zup1        = ( pum1(jl,jk) + 0.5_wp * zdum1 ) * zdum1
        zvp1        = ( pvm1(jl,jk) + 0.5_wp * zdvm1 ) * zdvm1
        ptte(jl,jk) = -( zup1 + zvp1 ) / ( ztmst * ( 1._wp + vtmpc2 * pqm1(jl,jk) ) * pcp(jl,jk) )

      ENDDO  !jl
    ENDDO  !jk
    
  END SUBROUTINE iondrag

END MODULE mo_upatmo_phy_iondrag
