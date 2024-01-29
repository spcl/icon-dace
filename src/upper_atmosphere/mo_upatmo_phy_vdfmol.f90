!
! Computation of molecular diffusion of momentum, heat and amount of water vapor
!
! vdf_mol calculates the temperature, humidity, u-wind, and v-wind
! tendencies due to molecular conductivity and viscosity.  It also
! calculates frictional heating in a physically meaningful way.
!
! Method:
!
! A semi-implicit scheme is employed to solve
!                              _             _
!     F(t+dt)-F(t-dt)      d  |          dF*  |
!     --------------- =  g -- | g mu rho --   |
!          2 dt            wp |          wp   |
!                              -             -
!
! where F* is alpha*F(t+dt) + (1-alpha)*F(t-dt) and mu is the
! dynamic (not kinematic!) viscosity.
!
! The boundary conditions are: 1) zero flux at the model top;
!                              2) zero flux at the model surface.
!
! This means solving a tridiagonal matrix system.
!
! M. Charron - MPI - September 6 2001.
!
! Modifications: H. Schmidt - MPI - 20020418
! - setting of tri-diag matrix cleaned
! H. Schmidt - MPI - 20020702
! - bug fix: msis variable index counts from bottom to top
! Th. Schoenemeyer/H. Schmidt - NEC/MPI - 20020708
! - optimized for vector architecture
!
! Rewrote for ICON by Guidi Zhou, MPI, 03.06.2016
!
! Modification by Guidi Zhou, MPI-M, 2016-06-20
! - make use of vertically-varying gravity
! Modification by Guidi Zhou, MPI-M, 2016-08-17
! - seperated frictional heating to mo_upatmo_fric
! Modification by Guidi Zhou, MPI-M (2017-03-07)
! - added the ability to compute molecular diffusion only above a certain altitude for performance
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

MODULE mo_upatmo_phy_vdfmol

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: SUCCESS
  USE mo_physical_constants,   ONLY: argas, tmelt
  USE mo_math_constants,       ONLY: dbl_eps
  USE mo_turb_vdiff_params,    ONLY: cvdifts

  IMPLICIT NONE

  PUBLIC :: vdf_mol
  ! for error handling
  INTEGER, PARAMETER :: IERR_NO     = SUCCESS      ! = 0
  INTEGER, PARAMETER :: IERR_SOLVER = SUCCESS + 1

CONTAINS

  !>
  !! Compute molecular diffusion
  !!
  !! Literature:
  !! - Huang, T., Walters, S., Brasseur, G., Hauglustaine, D., and Wu, W. (1998) 
  !!   Description of SOCRATES: A chemical dynamical radiative tow-dimensional model. 
  !!   Tech. Rep. NCAR/TN-440+EDD, NCAR, Boulder, CO.
  !! - Banks, P. M., and Kockarts, G. (1973) Aeronomy. Part B, Elsevier.
  !!
  SUBROUTINE vdf_mol(jcs, jce, kbdim, klev, ktracer, psteplen, ptvm1, ptm1, pqm1, pum1, pvm1, papm1, paphm1, grav, amu, & 
    &                ptte, pvom, pvol, pqte, opt_istartlev, opt_iendlev, opt_error)

    ! in/out variables
    
    INTEGER,  INTENT(IN)  :: jcs, jce, kbdim, klev, ktracer
    REAL(wp), INTENT(IN)  :: psteplen                ! time step length, usually 2*delta_time
    REAL(wp), INTENT(IN)  :: ptvm1(kbdim,klev)       ! virtual temperature   [TV]
    REAL(wp), INTENT(IN)  :: ptm1(kbdim,klev)        ! temperature           [T]
    REAL(wp), INTENT(IN)  :: pqm1(kbdim,klev,ktracer)! tracer concentration  [qtrc]
    REAL(wp), INTENT(IN)  :: pum1(kbdim,klev)        ! u-velocity            [u]
    REAL(wp), INTENT(IN)  :: pvm1(kbdim,klev)        ! v-velocity            [v]
    REAL(wp), INTENT(IN)  :: papm1(kbdim,klev)       ! full-level pressure   [pf]
    REAL(wp), INTENT(IN)  :: paphm1(kbdim,klev+1)    ! half-level pressure   [ph]
    REAL(wp), INTENT(IN)  :: grav(kbdim,klev)        ! gravity acceleration  [g]
    REAL(wp), INTENT(IN)  :: amu(kbdim,klev)         ! molar weight of air   [am]
    
    REAL(wp), INTENT(OUT) :: ptte(kbdim,klev)        ! temperature tendency  [dT/dt]
    REAL(wp), INTENT(OUT) :: pvom(kbdim,klev)        ! u-velocity tendency   [du/dt]
    REAL(wp), INTENT(OUT) :: pvol(kbdim,klev)        ! v-velocity tendency   [dv/dt]
    REAL(wp), INTENT(OUT) :: pqte(kbdim,klev,ktracer)! tracer concentration tendency [dqtrc/dt]

    INTEGER,  OPTIONAL, INTENT(IN)  :: opt_istartlev, opt_iendlev ! optional vertical start and end indices
    INTEGER,  OPTIONAL, INTENT(OUT) :: opt_error                  ! for optional error handling

    ! local variables
    
    INTEGER  :: jl, jk, jtr, istartlev, iendlev, iendlevm1, istartlevp1
    REAL(wp) :: zalpha, zmalpha, ztmst, zpr, zinvpr, zinvtmst
    REAL(wp) :: zmah, zrstar, zgvh, zrhoh, zth, ztvh
    REAL(wp) :: zinvmalpha, zmalpha_invpr, zinvmalpha_pr
    REAL(wp) :: zbet(kbdim)

    REAL(wp), ALLOCATABLE :: zgmurhoh(:,:), zgam(:,:), zhh(:,:), zrwp(:,:), zr(:,:), &
      &                      za(:,:), zb(:,:), zc(:,:), zla(:,:), zlb(:,:), zlc(:,:), zut(:,:)

    LOGICAL  :: lerror, l_present_error

    REAL(wp), PARAMETER :: eps = dbl_eps * 100._wp
    REAL(wp), PARAMETER :: inv_tmelt = 1._wp / tmelt ! (tmelt=273.15_wp)

    !---------------------------------------------------------
    
    ! please do not limit range of assignment 
    ! (e.g., ptte(jcs:jce,istartlev:iendlev) = 0._wp)), 
    ! since tendencies have attribute INTENT(OUT)
    ptte(:,:)   = 0._wp
    pvom(:,:)   = 0._wp
    pvol(:,:)   = 0._wp
    pqte(:,:,:) = 0._wp

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

    IF (istartlev >= iendlev) RETURN 

    istartlevp1 = istartlev + 1
    iendlevm1   = iendlev - 1

    ALLOCATE(zgmurhoh(kbdim,istartlev:iendlev), zgam(kbdim,istartlev:iendlev), &
      &      zhh(kbdim,istartlev:iendlev), zrwp(kbdim,istartlev:iendlev),      &
      &      zr(kbdim,istartlev:iendlev), za(kbdim,istartlev:iendlev),         &
      &      zb(kbdim,istartlev:iendlev), zc(kbdim,istartlev:iendlev),         &
      &      zla(kbdim,istartlev:iendlev), zlb(kbdim,istartlev:iendlev),       &
      &      zlc(kbdim,istartlev:iendlev), zut(kbdim,istartlev:iendlev))

    !** SET PARAMETERS
    
    zrstar        = argas * 1.e3_wp    ! universal gas constant 8314 J/kmol/K
    zpr           = 0.72_wp            ! Prandtl number
    zinvpr        = 1._wp / zpr
    zalpha        = cvdifts
    zmalpha       = 1._wp - zalpha
    zinvmalpha    = 1._wp / zmalpha
    zmalpha_invpr = zmalpha * zinvpr
    zinvmalpha_pr = zinvmalpha * zpr
    ztmst         = psteplen
    zinvtmst      = 1._wp / ztmst
   
    !** SET LAYERS' THICKNESS, (DYNAMIC VISCOSITY)*(pg)/(RT) 
    !   AND (DENSITY * TRACER DIFFUSION COEFF)
    
    DO jk = istartlevp1, iendlev
      DO jl = jcs, jce
        zhh(jl,jk)      = 1._wp / ( papm1(jl,jk-1) - papm1(jl,jk) )
        zrwp(jl,jk)     = ztmst * grav(jl,jk) / ( paphm1(jl,jk+1) - paphm1(jl,jk) )
        zgvh            = 0.5_wp * ( grav(jl,jk) + grav(jl,jk-1) )
        zmah            = 0.5_wp * ( amu(jl,jk) + amu(jl,jk-1) )
        zth             = 0.5_wp * ( ptm1(jl,jk) + ptm1(jl,jk-1) )
        ztvh            = 0.5_wp * ( ptvm1(jl,jk) + ptvm1(jl,jk-1) )
        zrhoh           = paphm1(jl,jk) * zmah / ( zrstar * ztvh )
        zgmurhoh(jl,jk) = 1.87E-5_wp * EXP( 0.69_wp * LOG( inv_tmelt * zth ) ) * zgvh * zrhoh
      ENDDO  !jl
    ENDDO  !jk

    zrwp(jcs:jce,istartlev) = ztmst * grav(jcs:jce,istartlev) / &
      &                       ( paphm1(jcs:jce,istartlevp1) - paphm1(jcs:jce,istartlev) )

    ! Please note: zla, zlb and zlc have to multiplied with the inverse Prandtl number, zinvpr,
    ! in case of the temperature equation!

    DO jk = istartlevp1, iendlevm1
      zla(jcs:jce,jk) = -zrwp(jcs:jce,jk) * ( zgmurhoh(jcs:jce,jk) * zhh(jcs:jce,jk) )
      zlc(jcs:jce,jk) = -zrwp(jcs:jce,jk) * ( zgmurhoh(jcs:jce,jk+1) * zhh(jcs:jce,jk+1) )
      zlb(jcs:jce,jk) = -( zla(jcs:jce,jk) + zlc(jcs:jce,jk) )
    ENDDO

    zlb(jcs:jce,istartlev) = zrwp(jcs:jce,istartlev) * zgmurhoh(jcs:jce,istartlevp1) * zhh(jcs:jce,istartlevp1)
    zlc(jcs:jce,istartlev) = -zlb(jcs:jce,istartlev)
    
    zla(jcs:jce,iendlev) = -zrwp(jcs:jce,iendlev) * zgmurhoh(jcs:jce,iendlev) * zhh(jcs:jce,iendlev)
    zlb(jcs:jce,iendlev) = -zla(jcs:jce,iendlev)

    !** SET TRI-DIAGONAL MATRIX: zc IS UPPER  DIAGONAL
    !                            zb IS MIDDLE DIAGONAL
    !                            za IS LOWER  DIAGONAL

    ! Please note: za and zc have to multiplied with the inverse Prandtl number, zinvpr,
    ! in case of the temperature equation!

    za(jcs:jce,istartlevp1:iendlev) = -zalpha * zla(jcs:jce,istartlevp1:iendlev)
    zc(jcs:jce,istartlev:iendlevm1) = -zalpha * zlc(jcs:jce,istartlev:iendlevm1)

    DEALLOCATE(zgmurhoh, zhh, zrwp)

    !   -----------------------
    !** 1: TEMPERATURE EQUATION
    !   -----------------------

    !** SET MIDDLE DIAGONAL OF TRI-DIAGONAL MATRIX: zb

    zb(jcs:jce,istartlev:iendlev) = 1._wp - zinvpr * zalpha * zlb(jcs:jce,istartlev:iendlev)
    
    !** SET RIGHT-HAND SIDE zr
    
    DO jk = istartlevp1, iendlevm1
      zr(jcs:jce,jk) = zmalpha_invpr * ( zla(jcs:jce,jk) * ptm1(jcs:jce,jk-1) &
        &            + ( zinvmalpha_pr + zlb(jcs:jce,jk) ) * ptm1(jcs:jce,jk) &
        &            + zlc(jcs:jce,jk) * ptm1(jcs:jce,jk+1) ) 
    ENDDO
    
    zr(jcs:jce,istartlev) = ( 1._wp + zmalpha_invpr * zlb(jcs:jce,istartlev) ) * ptm1(jcs:jce,istartlev) &
      &                   + zmalpha_invpr * zlc(jcs:jce,istartlev) * ptm1(jcs:jce,istartlevp1)
    zr(jcs:jce,iendlev)   = zmalpha_invpr * zla(jcs:jce,iendlev) * ptm1(jcs:jce,iendlevm1) &
      &                   + ( 1._wp + zmalpha_invpr * zlb(jcs:jce,iendlev) ) * ptm1(jcs:jce,iendlev)

    !** START OF MATRIX SOLVER (zut IS T(t+dt))
    
    lerror = .FALSE.

    DO jl = jcs, jce
      ! avoid division by zero
      IF (ABS(zb(jl,istartlev)) < eps) THEN
        lerror = .TRUE.
        EXIT 
      ENDIF
      zbet(jl)          = zb(jl,istartlev)
      zut(jl,istartlev) = zr(jl,istartlev) / zbet(jl)
    ENDDO  !jl

    IF (lerror) THEN
      IF (l_present_error) opt_error = IERR_SOLVER
      RETURN
    ENDIF

    !** DECOMPOSITION AND FORWARD SUBSTITUTION

    OUT_T: DO jk = istartlevp1, iendlev
      DO jl = jcs, jce
        zgam(jl,jk) = zinvpr * zc(jl,jk-1) / zbet(jl)
        zbet(jl)    = zb(jl,jk) - zinvpr * za(jl,jk) * zgam(jl,jk)
        ! avoid division by zero
        IF (ABS(zbet(jl)) < eps) THEN
          lerror = .TRUE.
          EXIT OUT_T
        ENDIF
        zut(jl,jk) = ( zr(jl,jk) - zinvpr * za(jl,jk) * zut(jl,jk-1) ) / zbet(jl)
      ENDDO  !jl
    ENDDO OUT_T

    IF (lerror) THEN
      IF (l_present_error) opt_error = IERR_SOLVER
      RETURN
    ENDIF

    !** BACKSUBSTITUTION

    DO jk = iendlevm1, istartlev, -1
      DO jl = jcs, jce
        zut(jl,jk) = zut(jl,jk) - zgam(jl,jk+1) * zut(jl,jk+1)
      ENDDO  !jl
    ENDDO  !jk

    ! temperature tendency
    ptte(jcs:jce,istartlev:iendlev) = zinvtmst * ( zut(jcs:jce,istartlev:iendlev) - ptm1(jcs:jce,istartlev:iendlev) )

    !   ------------------
    !** 2: U-WIND EQUATION
    !   ------------------

    !** RESET MIDDLE DIAGONAL OF TRI-DIAGONAL MATRIX: zb

    zb(jcs:jce,istartlev:iendlev) = 1._wp - zalpha * zlb(jcs:jce,istartlev:iendlev)

    !** SET RIGHT-HAND SIDE zr
    
    DO jk = istartlevp1, iendlevm1
      zr(jcs:jce,jk) = zmalpha * ( zla(jcs:jce,jk) * pum1(jcs:jce,jk-1)    &
        &            + ( zinvmalpha + zlb(jcs:jce,jk) ) * pum1(jcs:jce,jk) &
        &            + zlc(jcs:jce,jk) * pum1(jcs:jce,jk+1) )
    ENDDO
    
    zr(jcs:jce,istartlev) = ( 1._wp + zmalpha * zlb(jcs:jce,istartlev) ) * pum1(jcs:jce,istartlev) &
      &                   + zmalpha * zlc(jcs:jce,istartlev) * pum1(jcs:jce,istartlevp1)
    zr(jcs:jce,iendlev)   = zmalpha * zla(jcs:jce,iendlev) * pum1(jcs:jce,iendlevm1) &
      &                   + ( 1._wp + zmalpha * zlb(jcs:jce,iendlev) ) * pum1(jcs:jce,iendlev)

    !** START OF MATRIX SOLVER (zut IS u(t+dt))
    
    lerror = .FALSE.

    DO jl = jcs, jce
      ! avoid division by zero
      IF (ABS(zb(jl,istartlev)) < eps) THEN
        lerror = .TRUE.
        EXIT 
      ENDIF
      zbet(jl)          = zb(jl,istartlev)
      zut(jl,istartlev) = zr(jl,istartlev) / zbet(jl)
    ENDDO  !jl

    IF (lerror) THEN
      IF (l_present_error) opt_error = IERR_SOLVER
      RETURN
    ENDIF

    !** DECOMPOSITION AND FORWARD SUBSTITUTION

    OUT_U: DO jk = istartlevp1, iendlev
      DO jl = jcs, jce
        zgam(jl,jk) = zc(jl,jk-1) / zbet(jl)
        zbet(jl)    = zb(jl,jk) - za(jl,jk) * zgam(jl,jk)
        ! avoid division by zero
        IF (ABS(zbet(jl)) < eps) THEN
          lerror = .TRUE.
          EXIT OUT_U
        ENDIF
        zut(jl,jk) = ( zr(jl,jk) - za(jl,jk) * zut(jl,jk-1) ) / zbet(jl)
      ENDDO  !jl
    ENDDO OUT_U

    IF (lerror) THEN
      IF (l_present_error) opt_error = IERR_SOLVER
      RETURN
    ENDIF

    !** BACKSUBSTITUTION

    DO jk = iendlevm1, istartlev, -1
      DO jl = jcs, jce
        zut(jl,jk) = zut(jl,jk) - zgam(jl,jk+1) * zut(jl,jk+1)
      ENDDO  !jl
    ENDDO  !jk

    ! tendency of zonal wind component
    pvom(jcs:jce,istartlev:iendlev) = zinvtmst * ( zut(jcs:jce,istartlev:iendlev) - pum1(jcs:jce,istartlev:iendlev) )

    !   ------------------
    !** 3: V-WIND EQUATION
    !   ------------------
    
    !** SET RIGHT-HAND SIDE zr
    
    DO jk = istartlevp1, iendlevm1
      zr(jcs:jce,jk) = zmalpha * ( zla(jcs:jce,jk) * pvm1(jcs:jce,jk-1)    &
        &            + ( zinvmalpha + zlb(jcs:jce,jk) ) * pvm1(jcs:jce,jk) &
        &            + zlc(jcs:jce,jk) * pvm1(jcs:jce,jk+1) )
    ENDDO
    
    zr(jcs:jce,istartlev) = ( 1._wp + zmalpha * zlb(jcs:jce,istartlev) ) * pvm1(jcs:jce,istartlev) &
      &                   + zmalpha * zlc(jcs:jce,istartlev) * pvm1(jcs:jce,istartlevp1)
    zr(jcs:jce,iendlev)   = zmalpha * zla(jcs:jce,iendlev) * pvm1(jcs:jce,iendlevm1) &
      &                   + ( 1._wp + zmalpha * zlb(jcs:jce,iendlev) ) * pvm1(jcs:jce,iendlev)
    
    !** START OF MATRIX SOLVER (zut IS v(t+dt))
    
    lerror = .FALSE.

    DO jl = jcs, jce
      ! avoid division by zero
      IF (ABS(zb(jl,istartlev)) < eps) THEN
        lerror = .TRUE.
        EXIT 
      ENDIF
      zbet(jl)          = zb(jl,istartlev)
      zut(jl,istartlev) = zr(jl,istartlev) / zbet(jl)
    ENDDO  !jl

    IF (lerror) THEN
      IF (l_present_error) opt_error = IERR_SOLVER
      RETURN
    ENDIF

    !** DECOMPOSITION AND FORWARD SUBSTITUTION

    OUT_V: DO jk = istartlevp1, iendlev
      DO jl = jcs, jce
        zgam(jl,jk) = zc(jl,jk-1) / zbet(jl)
        zbet(jl)    = zb(jl,jk) - za(jl,jk) * zgam(jl,jk)
        ! avoid division by zero
        IF (ABS(zbet(jl)) < eps) THEN
          lerror = .TRUE.
          EXIT OUT_V
        ENDIF
        zut(jl,jk) = ( zr(jl,jk) - za(jl,jk) * zut(jl,jk-1) ) / zbet(jl)
      ENDDO  !jl
    ENDDO OUT_V

    IF (lerror) THEN
      IF (l_present_error) opt_error = IERR_SOLVER
      RETURN
    ENDIF

    !** BACKSUBSTITUTION

    DO jk = iendlevm1, istartlev, -1
      DO jl = jcs, jce
        zut(jl,jk) = zut(jl,jk) - zgam(jl,jk+1) * zut(jl,jk+1)
      ENDDO  !jl
    ENDDO  !jk

    ! tendency of meridional wind component
    pvol(jcs:jce,istartlev:iendlev) = zinvtmst * ( zut(jcs:jce,istartlev:iendlev) - pvm1(jcs:jce,istartlev:iendlev) )

    !   --------------------
    !** 4: TRACER EQUATION
    !   --------------------

    !** SET RIGHT-HAND SIDE zr
    
    DO jtr = 1, ktracer

      DO jk = istartlevp1, iendlevm1
        zr(jcs:jce,jk) = zmalpha * ( zla(jcs:jce,jk) * pqm1(jcs:jce,jk-1,jtr)    &
          &            + ( zinvmalpha + zlb(jcs:jce,jk) ) * pqm1(jcs:jce,jk,jtr) & 
          &            + zlc(jcs:jce,jk) * pqm1(jcs:jce,jk+1,jtr) )
      ENDDO
      
      zr(jcs:jce,istartlev) = ( 1._wp + zmalpha * zlb(jcs:jce,istartlev) ) * pqm1(jcs:jce,istartlev,jtr) &
        &                   + zmalpha * zlc(jcs:jce,istartlev) * pqm1(jcs:jce,istartlevp1,jtr)
      zr(jcs:jce,iendlev)   = zmalpha * zla(jcs:jce,iendlev) * pqm1(jcs:jce,iendlevm1,jtr) &
        &                   + ( 1._wp + zmalpha * zlb(jcs:jce,iendlev) ) * pqm1(jcs:jce,iendlev,jtr)

      !** START OF MATRIX SOLVER (zut IS q(t+dt))

      lerror = .FALSE.
      
      DO jl = jcs, jce
        ! avoid division by zero
        IF (ABS(zb(jl,istartlev)) < eps) THEN
          lerror = .TRUE.
          EXIT 
        ENDIF
        zbet(jl)          = zb(jl,istartlev)
        zut(jl,istartlev) = zr(jl,istartlev) / zbet(jl)
      ENDDO  !jl
      
      IF (lerror) THEN
        IF (l_present_error) opt_error = IERR_SOLVER
        RETURN
      ENDIF
      
      !** DECOMPOSITION AND FORWARD SUBSTITUTION
      
      OUT_Q: DO jk = istartlevp1, iendlev
        DO jl = jcs, jce
          zgam(jl,jk) = zc(jl,jk-1) / zbet(jl)
          zbet(jl)    = zb(jl,jk) - za(jl,jk) * zgam(jl,jk)
          ! avoid division by zero
          IF (ABS(zbet(jl)) < eps) THEN
            lerror = .TRUE.
            EXIT OUT_Q
          ENDIF
          zut(jl,jk) = ( zr(jl,jk) - za(jl,jk) * zut(jl,jk-1) ) / zbet(jl)
        ENDDO  !jl
      ENDDO OUT_Q
      
      IF (lerror) THEN
        IF (l_present_error) opt_error = IERR_SOLVER
        RETURN
      ENDIF
      
      !** BACKSUBSTITUTION
      
      DO jk = iendlevm1, istartlev, -1
        DO jl = jcs, jce
          zut(jl,jk) = zut(jl,jk) - zgam(jl,jk+1) * zut(jl,jk+1)
        ENDDO  !jl
      ENDDO  !jk

      pqte(jcs:jce,istartlev:iendlev,jtr) = zinvtmst * ( zut(jcs:jce,istartlev:iendlev) &
        &                                 - pqm1(jcs:jce,istartlev:iendlev,jtr) )

    END DO ! ktracer
    
    DEALLOCATE(zgam, zr, za, zb, zc, zla, zlb, zlc, zut)
    
  END SUBROUTINE vdf_mol
  
END MODULE mo_upatmo_phy_vdfmol
