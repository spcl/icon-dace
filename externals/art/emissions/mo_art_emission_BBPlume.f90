!NEC$ options "-finline-max-depth=3 -finline-max-function-size=4000 "
!
! mo_art_emission_BBPlume
! This module provides plume heights for a biomass burning aerosol model that considers
! plume rise. This is a reimplementation of the numerics while microphysics still mostly
! equals to the original version of mo_art_emission_plumerise .
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
! 33, L17808, doi:10.1029/2006GL026608, 2006.                                               !
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

MODULE mo_art_emission_BBPlume
!#ifdef __ICON_ART

  ! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
  USE mo_math_constants,                ONLY: pi
  USE mo_physical_constants,            ONLY: rd, cpd, cvd, rv, p0ref, grav, earth_radius, &
                                            & stbo, tmelt, alv, als, amw, amd

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine = 'mo_art_emission_BBPlume'

  INTEGER :: max_minutes = 40       ! steady state is usually reached after 50 minutes - Freitas

  REAL(wp), PARAMETER :: &
       & vtr = -4._wp,           &  ! initial rain velocity - small variation with qr
       & vti = -3._wp,           &  ! inital ice velocity - small variation with qi
       & gama = 0.5_wp,          &  ! mass virtual coeff.
       & eps = amw/amd,          &  ! fuel moisture fraction
       & fmoist = 0.1_wp,        &  ! fuel moisture fraction
       & alpha = 0.05_wp,        &  ! entrainment constant
       & HEAT = 19.3E6_wp,       &  ! J/kg - rain forest heat fuel
       & s_p_m = 60._wp             ! seconds per minute


  PUBLIC :: BB_plume, a_priori_BB_plume

CONTAINS

  SUBROUTINE  a_priori_BB_plume(heat_f,area,pe,l_cpr)
    ! this routine check whether a plume rise computation is necessary or not
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: &
         & pe,              & ! heatflux
         & heat_f,          & ! heatflux
         & area               ! burnt area
    LOGICAL, INTENT(OUT) :: &
         l_cpr                ! logical if true a Computation of Plume Rise is necessary
    INTEGER :: &
         & j                  ! horizontal loop index
    REAL(wp) :: &
         & R,               & ! plume radius
         & eflux,           & ! energy flux by the fire
         & f_b,             & ! buoyancy flux
         & w                  ! vertical velocity near ground

    l_cpr = .TRUE.

    R     = SQRT(area / pi)                         ! entrainment surface radius (m)
    eflux = heat_f * 0.55_wp                        ! W/m**2
    f_b   = grav * rd * eflux / (pe * cpd) * R * R  ! buoyancy flux
    w  = 5._wp/(6._wp*alpha)*((0.9_wp*alpha*f_b)**(1._wp/3._wp))/((5._wp / (6._wp * alpha)*R)**(1._wp/3._wp))
    ! check vertical velocity w near ground
    IF (w .LE. 1._wp) THEN
       l_cpr = .FALSE.
    END IF

  END SUBROUTINE a_priori_BB_plume

  SUBROUTINE  BB_plume(nlev,heat_f,area,Te_in,pe_in,ue, &
       &            ve,qve_in,z_mc_in,z_ifc_in,ztop)
    IMPLICIT NONE
    ! plume height calculation

    INTEGER, INTENT(IN)  :: &
         & nlev               ! nlev - index of the highest layer (from bottom to top)
    REAL(wp), INTENT(IN) :: &
         & heat_f,          & ! heatflux
         & area               ! burnt area
    REAL(wp), INTENT(IN) :: &
         & Te_in(nlev),        & ! temperatur
         & pe_in(nlev),        & ! pressure
         & ue(nlev),           & ! zonal wind (m/s)
         & ve(nlev),           & ! meridional wind (m/s)
         & z_mc_in(nlev),      & ! full level heights
         & z_ifc_in(nlev+1),     & ! half level heights
         & qve_in(nlev)          ! water vapor

    REAL(wp), INTENT(OUT) :: &
         & ztop               ! height for biomass emissions
    INTEGER  ::  &
         & k_max,           & ! index of ztop in vertical grid z_mc
         & k                  ! vertical index
    REAL(wp) ::  &
         & T(nlev),         & ! temperatur
         & qc(nlev),        & ! cloud droplet concentration
         & qr(nlev),        & ! rain drop concentration
         & qi(nlev),        & ! cloud ice concentration
         & qv(nlev),        & ! water vapor concentration
         & w(nlev+1),       & ! vertical velocity
         & KT(nlev),        & ! viscosity constant
         & R(nlev),         & ! entrainment plume radius
         & nu(nlev),        & ! horizontal wind
         & Te(nlev),        & ! environmental temperatur
         & pe(nlev),        & ! environmental pressure
         & z_mc(nlev),      & ! full level heights
         & z_ifc(nlev+1),   & ! half level heights
         & h_mc(nlev-1),    & ! full level thickness
         & h_ifc(nlev),     & ! half level thickness
         & qve(nlev),       & ! environmentalwater vapor
         & nue(nlev),       & ! environmental horizontal wind
         & maxtime,         & ! max time in seconds
         & dampH1,          & ! initial Height for gravity wave damping
         & damp_delta,      & ! gravity wave height difference
         & dt = 1._wp,      & ! time step
         & time = 0._wp,    & ! time coordinate
         & ztop_temp = 0._wp, & !
         & ztop_min(max_minutes+1), & !
         & zsurf              ! surface height
    INTEGER :: &
         & nmax               ! top level of current plume

    ! set ouput variable and check input
    ztop       = 0._wp
    IF((area < 1.e-6_wp) .OR. (heat_f < 1.e-6_wp)) THEN
       RETURN
    ENDIF

    ! initialize qv, vertical grid z_mc, z_ifc and effective horizontal wind nu
    ! index e marks the environmental conditions
    zsurf = z_ifc_in(1)                         ! surface height
    z_ifc(1) = 0._wp
    DO k=1,nlev
       Te(k)      = Te_in(k)
       T(k)       = Te(k)
       pe(k)      = pe_in(k)
       qve(k)     = MAX(qve_in(k),1.e-8_wp)     ! minimal watercontent
       qv(k)      = qve(k)                      ! qv equals qve initially
       nue(k)     = SQRT(ue(k)**2+ve(k)**2)     ! effectiv horizontal wind velocity
       z_mc(k)    = z_mc_in(k) - zsurf          ! thermo and water levels from the surface
       z_ifc(k+1) = z_ifc_in(k+1) - zsurf       ! dynamical levels from tthe surface
    ENDDO

    h_ifc(1) = z_ifc(2) - z_ifc(1)
    DO k=1,nlev-1
       h_mc(k) = z_mc(k+1)-z_mc(k)
       h_ifc(k+1) = z_ifc(k+2)-z_ifc(k+1)
    END DO
    ! initialize plume radius R, vertial velocity w and viscosity K
    R(1) = SQRT (area / pi)                  ! entrainment surface radius (m)
    nu = 0._wp
    w  = 0._wp
    qc = 0._wp
    qr = 0._wp
    qi = 0._wp
    KT(1) = 500._wp                          ! viscosity constant
    DO k=2,nlev
       R(k) = R(k-1)+(6._wp/5._wp)*alpha*h_mc(k-1)
       KT(k) = MAX(1.e-3_wp,KT(1) * (1._wp - z_mc(k-1)/2.e4_wp)) ! from the old code - no explanation so far
    ENDDO

    ! initialize time stepping and damping layer for Rayleigh friction
    maxtime = REAL(max_minutes,wp)*s_p_m
    ! if dt is chosen larger than 20._wp seconds
    ! the implementation of The Rayleigh friction layershould be rivised
    dt   = 10._wp  ! for crank-nicolson one would choose dt .le. dx with dx = 20m in ICON
    time = dt
    ztop_temp  = 0._wp
    dampH1     = 1000._wp  ! height in meter
    damp_delta = 2000._wp ! height in meter

    k = 1
    DO WHILE ((z_mc(k) < dampH1 + damp_delta) .AND. (k < nlev))
       k = k+1
    END DO
    k_max = 1
    nmax = MAX(1,k-1)

    ! time loop with time dependent early exit criterion
    DO WHILE (time < maxtime)
       ! define index for integration domain
       k = k_max
       DO WHILE ((z_mc(k) <  z_ifc(k_max+1) + damp_delta) .AND. (k < nlev))
          k = k+1
       END DO
       nmax = MAX(nmax,k -1)
       nmax = MIN(nmax,nlev)

       ! call time stepping - implicit or explicit method can be chosen inside
       CALL time_step(nlev =       nmax, &
            &   damp_delta = damp_delta, &
            &           dt =         dt, &
            &         time =       time, &
            &       heat_f =     heat_f, &
            &            w =          w, &
            &           te =         te, &
            &           pe =         pe, &
            &          qve =        qve, &
            &          nue =        nue, &
            &            R =          R, &
            &            T =          T, &
            &           qv =         qv, &
            &           qc =         qc, &
            &           qr =         qr, &
            &           qi =         qi, &
            &           nu =         nu, &
            &         z_mc =       z_mc, &
            &        z_ifc =      z_ifc, &
            &         h_mc =       h_mc, &
            &        h_ifc =      h_ifc, &
            &           KT =         KT)

       !-- try to find the plume top (above surface height)
       k = 1
       DO WHILE ((w (k) > 1._wp) .AND. (k < nmax-1))
          k = k + 1
       ENDDO

       ztop_temp =  z_mc(k)

       k_max = MAX( k, k_max)
       ztop  = MAX (ztop, ztop_temp)

       ztop_min(INT(time/s_p_m)+1:INT((time+dt)/s_p_m)+1) = ztop

       ! early exit criterion
       ! if the solution is going to a stationary phase, exit
       IF(time > 10._wp * s_p_m) THEN
          IF( ABS(ztop_min(INT(time/s_p_m))-ztop_min(INT(time/s_p_m)-5)) < 100._wp ) THEN
             EXIT
          ENDIF
       ENDIF

       ! update time
       time = time + dt
    END DO

  END SUBROUTINE BB_plume

  SUBROUTINE time_step(nlev,damp_delta,dt,time,heat_f,w,te,pe,qve,nue,R,T,qv,qc,qr,qi,nu,z_mc,z_ifc,h_mc,h_ifc,KT)
    ! we apply the Crank-Nicolson algorithm for all but the lowest and uppermost level
    ! (x_i+1 - x_i)/(t_i+1 -t_i) = .5*A*x_i+1 + .5*A*x_i + F where F is the PDE right hand side
    ! with a tridiagonal matrix A = [A(1,:), A(2,:), A(3,:)]  with the three diagonals A1, A2 and A3
    ! and i indexing time level, the vertical index is omitted here
    ! this results in a linear equation Mx=B with a tridiagonal matrix M = [M(1,:), M(2,:), M(3,:)]
    ! the Thomas algorithm is then applied to solve this system
    ! some simplifications apply:
    ! - we consider w_i instead of w_i+1 outside the vertical derivativs
    ! - the i+1 state is only considered for each respective variable (no implicit cross relations)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: &
       & nlev              ! number of levels
  REAL(wp),INTENT(IN) :: &
       & damp_delta,     & ! height difference for gravity wave damping
       & dt,             & ! time step
       & time,           & ! time coordinate
       & heat_f,         & ! heat flux
       & pe(:),          & ! air pressure
       & te(:),          & !
       & qve(:),         & !
       & nue(:),         & !
       & z_mc(:),        & ! full level heights
       & z_ifc(:),       & ! half level heights
       & h_mc(:),        & ! full level thickness
       & h_ifc(:),       & ! half level thickness
       & KT(:)             ! viscosity constant
  REAL(wp), INTENT(INOUT) :: &
       & nu(:),          & !
       & w(:),           & !
       & T(:),           & !
       & qc(:),          & !
       & qr(:),          & !
       & qi(:),          & !
       & qv(:),          & !
       & R(:)              !
  REAL(wp) :: &
       & cvi,            & !
       & dqsdz,          & !
       & entrain,        & !
       & entrain_R,      & !
       & entrain_w,      & !
       & h_p,            & ! h plus
       & h_m,            & ! h minus
       & h_1,            & !
       & h_2,            & !
       & h_3,            & !
       & wbar,           & !
       & nu_nue_bar,     & !
       & R_bar,          & !
       & F_w(nlev),      & !
       & F_T(nlev-1),    & !
       & F_qv(nlev-1),   & !
       & F_qc(nlev-1),   & !
       & F_qr(nlev-1),   & !
       & F_qi(nlev-1),   & !
       & F_nu(nlev-1),   & !
       & F_R(nlev-1),    & !
       & F_micro_t(nlev-1),             & ! temperature
       & F_micro_qv(nlev-1),            & ! water vapor fraction
       & F_micro_qc(nlev-1),            & ! cloud water fraction
       & F_micro_qr(nlev-1),            & ! rain  water fraction
       & F_micro_qi(nlev-1),            & ! ice         fraction
       & A_w(3,nlev),    & !
       & A_T(3,nlev-1),  & !
       & A_qv(3,nlev-1), & !
       & A_qc(3,nlev-1), & !
       & A_qr(3,nlev-1), & !
       & A_qi(3,nlev-1), & !
       & A_nu(3,nlev-1), & !
       & A_R(3,nlev-1),  & !
       & eflux,          & ! energy flux by the fire
       ! & rho(0:nlev),      & ! air density
       & rho(nlev),      & ! air density
       & ps(nlev),       & ! saturation vapor pressure
       & qsat(nlev),     & ! qv at water saturation
       & Buoy(nlev),     & ! buoyancy
       & f_b,            & ! temporary variable
       & denscor           ! temporary variable
  INTEGER :: &
       & damp_k,         & ! index for gravity wave damping
       & k                 ! loop index
  LOGICAL :: &
       & cond1,          & ! temporary variabl e to store conditions for Godunov coefficients
       & cond2             ! temporary variabl e to store conditions for Godunov coefficients

  ! heat ramp up over 5 minutes
  IF (time <= 60._wp) THEN
     eflux = 0.1_wp
  ELSEIF ((time <= 120._wp) .AND. (time > 60._wp)) THEN
     eflux = 0.25_wp * heat_f  * 0.55_wp     ! W/m**2 (0.55 conversion for convective energy)
  ELSEIF ((time <= 180._wp) .AND. (time > 120._wp)) THEN
     eflux = 0.5_wp  * heat_f  * 0.55_wp     ! W/m**2
  ELSEIF ((time <= 240._wp) .AND. (time > 180._wp)) THEN
     eflux = 0.75_wp * heat_f  * 0.55_wp     ! W/m**2
  ELSE
     eflux =           heat_f  * 0.55_wp     ! W/m**2
  END IF

  ! find gravity wave damping index
  k = 1
  DO WHILE ((z_mc(k) < z_mc(nlev) - damp_delta) .AND. (k < nlev))
     k = k+1
  END DO
  damp_k = MAX(1,k-1)

  !---------------------------
  ! lower boundary condition
  ! from old code - no explanation so far
  qr (1) = qr (2)   !soak up hydrometeors
  qi (1) = qi (2)

  f_b      = grav * rd * eflux / (pe(1) * cpd) * R(1) * R(1)  !buoyancy flux

  w (1)    = 5._wp/(6._wp*alpha)*((0.9_wp*alpha*f_b)**(1._wp/3._wp))/((5._wp / (6._wp * alpha)*R(1))**(1._wp/3._wp))

  denscor  = 5._wp / (6._wp * alpha) * f_b / grav / (0.9_wp * alpha * f_b) **(1._wp/3._wp) &
       &     / (5._wp / (6._wp * alpha)*R(1))**(5._wp/3._wp)   !density correction
  T(1)     = Te (1) / (1._wp - denscor)    !temperature of virtual plume at zsurf

  rho(1)   = pe(1) / rd / T(1)  ! air density
  qv(1)    = (eflux * (dt / heat) * (0.5_wp + fmoist) /0.55_wp) / (w(1) * dt * rho(1) ) + qve(1)
  ps(1)    = psat_water(T=T(1))
  qsat (1) = (eps * ps(1)) / pe (1)   !blob saturation lwc kg/kg dry air

  IF (qv (1) > qsat (1) ) THEN
     qc (1) = qv   (1) - qsat (1)  !remainder goes into cloud drops
     qv (1) = qsat (1)
  ENDIF

  !---------------------------
  ! update density and saturation as algebraic conditions
  DO k=2,nlev
     ps(k)   = psat_water(T=T(k))
     rho(k)  = pe(k) / rd / T(k)  ! air density
     qsat(k) = (eps * ps(k)) / pe (k)   !blob saturation lwc kg/kg dry air
  END DO

  !---------------------------
  ! update lower boundary condition with microphysics
  cvi      = 1.6_wp + 0.57e-3_wp * (ABS (w (1) + vti) ) **1.5_wp / 0.75_wp
  wbar     = ((z_mc(1)-z_ifc(1))*w(1)+(z_ifc(2)-z_mc(1))*w(2))/h_ifc(1)
  dqsdz    = (rho(2)*qsat(2)/rho(1) - qsat(1))/ h_mc(1)
  CALL microphysics(dt =            dt, &
       &          t_in =          T(1), &
       &           rho =        rho(1), &
       &          qsat =       qsat(1), &
       &          wbar =          wbar, &
       &         dqsdz =         dqsdz, &
       &           est =         ps(1), &
       &           cvi =           cvi, &
       &         qv_in =         qv(1), &
       &         qc_in =         qc(1), &
       &         qh_in =         qr(1), &
       &         qi_in =         qi(1), &
       &             t =  F_micro_t(1), &
       &            qv = F_micro_qv(1), &
       &            qc = F_micro_qc(1), &
       &            qh = F_micro_qr(1), &
       &            qi = F_micro_qi(1))

  t (1) = t (1) + dt *F_micro_t (1)
  qv(1) = qv(1) + dt *F_micro_qv(1)
  qc(1) = qc(1) + dt *F_micro_qc(1)
  qr(1) = qr(1) + dt *F_micro_qr(1)
  qi(1) = qi(1) + dt *F_micro_qi(1)

  !---------------------------
  ! prepare Micro physics
  DO k=2,nlev-1
     ! w (k) + vti relative velocity to surrounding cloud
     ! ice ventilation coefficient for sublimation
     cvi      = 1.6_wp + 0.57e-3_wp * (ABS (w (k) + vti) ) **1.5_wp / 0.75_wp
     wbar     = ((z_mc(k)-z_ifc(k))*w(k)+(z_ifc(k+1)-z_mc(k))*w(k+1))/h_ifc(k)
     dqsdz    = (rho(k+1)*qsat(k+1) - rho(k-1)*qsat(k-1))/ (rho(k) * (z_mc(k+1) - z_mc(k-1)))
     CALL microphysics(dt =            dt, &
          &          t_in =          T(k), &
          &           rho =        rho(k), &
          &          qsat =       qsat(k), &
          &          wbar =          wbar, &
          &         dqsdz =         dqsdz, &
          &           est =         ps(k), &
          &           cvi =           cvi, &
          &         qv_in =         qv(k), &
          &         qc_in =         qc(k), &
          &         qh_in =         qr(k), &
          &         qi_in =         qi(k), &
          &             t =  F_micro_t(k), &
          &            qv = F_micro_qv(k), &
          &            qc = F_micro_qc(k), &
          &            qh = F_micro_qr(k), &
          &            qi = F_micro_qi(k))
  END DO

  !---------------------------
  ! prepare buoyancy
  CALL buoyancy(nlev = nlev, &
       &           T =    T, &
       &          Te =   Te, &
       &          qv =   qv, &
       &         qve =  qve, &
       &          qr =   qr, &
       &          qi =   qi, &
       &          qc =   qc, &
       &        Buoy = Buoy)

  !---------------------------
  ! prepare coefficents for w
  DO k=2,nlev
     ! diffusion coefficients
     h_p        = KT(k)*2._wp/(h_ifc(k)  *(h_ifc(k)+h_ifc(k-1)))
     h_m        = KT(k)*2._wp/(h_ifc(k-1)*(h_ifc(k)+h_ifc(k-1)))
     ! advection coefficients
     cond1 = (w(k-1) >= w(k)).AND.(ABS(w(k-1)) >= ABS(w(k)))
     cond2 = (w(k-1) <  w(k)).AND.(ABS(w(k-1)) <  ABS(w(k)))
     IF (cond1.OR.cond2) THEN ! chose coefficients according to Godunov
        h_1        =  w(k-1)/h_mc(k-1)
        h_2        =  0._wp
        h_3        =  0._wp
     ELSE
        h_1        =  0._wp
        h_2        =  w(k)/h_mc(k-1)
        h_3        =  0._wp
     END IF
     cond1 = (w(k) >= w(k+1)).AND.(ABS(w(k)) >= ABS(w(k+1)))
     cond2 = (w(k) <  w(k+1)).AND.(ABS(w(k)) <  ABS(w(k+1)))
     IF (cond1.OR.cond2) THEN !
        h_1        =  h_1 - 0._wp
        h_2        =  h_2 - w(k)/h_mc(k-1)
        h_3        =  h_3 - 0._wp
     ELSE
        h_1        =  h_1 - 0._wp
        h_2        =  h_2 - 0._wp
        h_3        =  h_3 - w(k+1)/h_mc(k-1)
     END IF
     ! prep A and F
     nu_nue_bar = ((z_ifc(k)-z_mc(k-1))*(nu(k-1)-nue(k-1))+(z_mc(k)-z_ifc(k))*(nu(k)-nue(k)))/h_mc(k-1)
     R_bar      = ((z_ifc(k)-z_mc(k-1))*R(k-1)+(z_mc(k)-z_ifc(k))*R(k))/h_mc(k-1)
     entrain_w  = (2._wp*alpha*ABS(w(k))+ABS(nu_nue_bar)*2._wp/pi)/(0.5_wp*(R(k-1)+R(k)))
     F_w(k)     = ((z_ifc(k)-z_mc(k-1))*Buoy(k-1) + (z_mc(k)-z_ifc(k))*Buoy(k))/h_mc(k-1)
     A_w(2,k)   = -entrain_w + h_2 - h_p - h_m
     A_w(1,k)   =  h_1 + h_m
     A_w(3,k)   =  h_3 + h_p
     IF (z_mc(k) >= z_mc(damp_k)) THEN
        ! Rayleigh friction from the original implementation
        ! but without time step dt dependence and a 20 second relaxation time
        ! A_w(2,k) = A_w(2,k) - (z_mc(k) - z_mc(damp_k)) / (z_mc(nlev) - z_mc(damp_k)) / (MIN(3._wp*dt,60._wp))
        A_w(2,k) = A_w(2,k) - (z_mc(k) - z_mc(damp_k)) / (z_mc(nlev) - z_mc(damp_k)) / (3._wp*20._wp)
     END IF
  END DO

  !---------------------------
  ! prepare coefficients for remaining prognostic variables
  DO k=2,nlev-1
     ! diffusion coefficients
     h_p        = KT(k)*2._wp/(h_mc(k)  *(h_mc(k)+h_mc(k-1)))
     h_m        = KT(k)*2._wp/(h_mc(k-1)*(h_mc(k)+h_mc(k-1)))
     wbar       = ((z_mc(k)-z_ifc(k))*w(k)+(z_ifc(k+1)-z_mc(k))*w(k+1))/h_ifc(k)
     ! advection coefficients
     IF (wbar >= 0._wp) THEN ! Up wind depending on the wind
        h_1        =  wbar/h_mc(k-1)
        h_2        = -wbar/h_mc(k-1)
        h_3        =  0._wp
     ELSE
        h_1        =  0._wp
        h_2        =  wbar/h_mc(k)
        h_3        = -wbar/h_mc(k)
     END IF
     entrain    = (2._wp*alpha*ABS(wbar)+ABS(nu(k)-nue(k-1))*2._wp/pi)/R(k)
     ! coefficients for T
     F_T(k)     = -wbar*grav/cpd + entrain * Te(k) + F_micro_t(k)
     A_T(2,k)   = - entrain + h_2 - h_p - h_m
     A_T(1,k)   =  h_1 + h_m
     A_T(3,k)   =  h_3 + h_p
     ! gravity wave damping
     IF (z_mc(k) >= z_mc(damp_k)) THEN
        ! Rayleigh friction taken from the original implementation
        ! but without time step dt dependence and a 20 second relaxation time
        ! F_T(k)   = F_T(k)   + (z_mc(k) - z_mc(damp_k)) / (z_mc(nlev) - z_mc(damp_k)) / (MIN(3._wp*dt,60._wp)) *Te(k) !*Te(nlev)
        ! A_T(2,k) = A_T(2,k) - (z_mc(k) - z_mc(damp_k)) / (z_mc(nlev) - z_mc(damp_k)) / (MIN(3._wp*dt,60._wp))
        F_T(k)   = F_T(k)   + (z_mc(k) - z_mc(damp_k)) / (z_mc(nlev) - z_mc(damp_k)) / (3._wp*20._wp) *Te(k) !*Te(nlev)
        A_T(2,k) = A_T(2,k) - (z_mc(k) - z_mc(damp_k)) / (z_mc(nlev) - z_mc(damp_k)) / (3._wp*20._wp)
     END IF
     ! coefficients for qv
     F_qv(k)    =  entrain * qve(k) + F_micro_qv(k)
     A_qv(2,k)  = -entrain + h_2 - h_p - h_m
     A_qv(1,k)  =  h_1 + h_m
     A_qv(3,k)  =  h_3 + h_p
     ! coefficients for qc
     F_qc(k)    =  0._wp            + F_micro_qc(k)
     A_qc(2,k)  = -entrain + h_2 - h_p - h_m
     A_qc(1,k)  =  h_1 + h_m
     A_qc(3,k)  =  h_3 + h_p
     ! coefficients for qr with sedimentation
     F_qr(k)    =  0._wp            + F_micro_qr(k)
     A_qr(2,k)  = -entrain + h_2 - h_p - h_m + vtr / h_mc(k)
     A_qr(1,k)  =  h_1 + h_m
     A_qr(3,k)  =  h_3 + h_p - vtr * rho(k+1) / rho(k) /h_mc(k)
     ! coefficients for qi with sedimentation
     F_qi(k)    =  0._wp            + F_micro_qi(k)
     A_qi(2,k)  = -entrain + h_2 - h_p - h_m + vti * h_mc(k)
     A_qi(1,k)  =  h_1 + h_m
     A_qi(3,k)  =  h_3 + h_p - vti * rho(k+1) / rho(k) /h_mc(k)
     ! coefficients for nu
     F_nu(k)    =  entrain* nue(k)
     A_nu(2,k)  = -entrain + h_2 - h_p - h_m
     A_nu(1,k)  =  h_1 + h_m
     A_nu(3,k)  =  h_3 + h_p
     ! coefficients for R
     entrain_R  =  (6._wp/5._wp*alpha*ABS(wbar)+ABS(nu(k)-nue(k-1))/pi)/R(k)
     F_R(k)     =  0._wp
     A_R(2,k)   =  entrain_R + h_2 - h_p - h_m
     A_R(1,k)   =  h_1 + h_m
     A_R(3,k)   =  h_3 + h_p
  END DO

  ! solve the prognostic vairables implicitly
  CALL solve_impl_time_step( nlev = nlev+1, &
       &                       dt =     dt, &
       &                        x =      w, &
       &                      F_x =    F_w, &
       &                      A_x =    A_w)
  CALL solve_impl_time_step( nlev = nlev, &
       &                       dt =   dt, &
       &                        x =    T, &
       &                      F_x =  F_T, &
       &                      A_x =  A_T)
  CALL solve_impl_time_step( nlev = nlev, &
       &                       dt =   dt, &
       &                        x =   qv, &
       &                      F_x = F_qv, &
       &                      A_x = A_qv)
  CALL solve_impl_time_step( nlev = nlev, &
       &                       dt =   dt, &
       &                        x =   qc, &
       &                      F_x = F_qc, &
       &                      A_x = A_qc)
  CALL solve_impl_time_step( nlev = nlev, &
       &                       dt =   dt, &
       &                        x =   qr, &
       &                      F_x = F_qr, &
       &                      A_x = A_qr)
  CALL solve_impl_time_step( nlev = nlev, &
       &                       dt =   dt, &
       &                        x =   qi, &
       &                      F_x = F_qi, &
       &                      A_x = A_qi)
  CALL solve_impl_time_step( nlev = nlev, &
       &                       dt =   dt, &
       &                        x =   nu, &
       &                      F_x = F_nu, &
       &                      A_x = A_nu)
  CALL solve_impl_time_step( nlev = nlev, &
       &                       dt =   dt, &
       &                        x =    R, &
       &                      F_x =  F_R, &
       &                      A_x =  A_R)

END SUBROUTINE time_step


SUBROUTINE solve_impl_time_step(nlev,dt,x,F_x,A_x)
  ! implicit time stepping with crank-nicolson and thomas algosrithm
  ! solve Mx=B
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: &
       & nlev               ! number of levels
  REAL(wp), INTENT(IN) :: &
       & dt,              & ! time step
       & F_x(:),          & ! sinks and sources
       & A_x(:,:)           ! coefficients
  REAL(wp), INTENT(INOUT) :: &
       & x(:)               ! prognositc variable
  REAL(wp) :: &
       & M(3,nlev),       & ! M*x = B
       & B(nlev),         & ! M*x = B
       & C(nlev),         & ! Thomas algorithm coefficients
       & D(nlev)            ! Thomas algorithm coefficients
  INTEGER :: &
       & k                  ! loop index

  C(1)    = 0._wp
  D(1)    = x(1)

  DO k = 2,nlev-1
     ! Crank Nicolson solve
     M(1,k) =        - dt / 2._wp * A_x(1,k)
     M(2,k) =  1._wp - dt / 2._wp * A_x(2,k)
     M(3,k) =        - dt / 2._wp * A_x(3,k)
     B(k)   = (1._wp + dt / 2._wp * A_x(2,k))*x(k) + dt/2._wp*A_x(1,k)*x(k-1) &
            + dt/2._wp*A_x(3,k)*x(k+1) + dt * F_x(k)
     ! only implicit solve - higher order cannot be achieved anyway
     ! M(1,k) =        - dt * A_x(1,k)
     ! M(2,k) =  1._wp - dt * A_x(2,k)
     ! M(3,k) =        - dt * A_x(3,k)
     ! B(k)   = x(k) + dt * F_x(k)
  END DO
!NEC$ novector
  DO k = 2,nlev-1
     C(k)   = M(3,k) / (M(2,k)-C(k-1)*M(1,k))
     D(k)   = (B(k)-D(k-1)*M(1,k)) / (M(2,k)-C(k-1)*M(1,k))
  END DO
!NEC$ novector
  DO k=nlev-1,2,-1
     x(k) = D(k)-C(k)*x(k+1)
  END DO

END SUBROUTINE solve_impl_time_step


SUBROUTINE  buoyancy(nlev, T, Te, qv, qve, qr, qi, qc,Buoy)

  IMPLICIT NONE

  INTEGER :: k ! loop counter
  INTEGER, INTENT(IN) :: nlev ! index upper boundary plume
  REAL(wp), DIMENSION(nlev), INTENT(IN) :: &
       & T,     & !
       & Te,    & !
       & qv,    & !
       & qve,   & !
       & qr,    & !
       & qi,    & !
       & qc       !
  REAL(wp), DIMENSION(nlev), INTENT(OUT) :: &
       & Buoy     !
  REAL(wp) :: &
       & Tv,     & !
       & Tve,    & !
       & qwtotl, & !
       & umgamai   !

  !- orig
  umgamai = 1._wp/(1._wp+gama) ! compensa a falta do termo de aceleracao associado `as
  ! das pertubacoes nao-hidrostaticas no campo de pressao
  DO k = 1,nlev

     Tv =   T(k) * (1._wp + (qv(k) /eps))/(1._wp + qv(k) )  !blob virtual temp.
     Tve = Te(k) * (1._wp + (qve(k)/eps))/(1._wp + qve(k))  !and environment

     qwtotl = qr(k) + qi(k) + qc(k)                       ! qwtotl*g is drag

     Buoy(k)= grav*  umgamai*( (tv - tve) / tve   - qwtotl)
  ENDDO

END SUBROUTINE  buoyancy




SUBROUTINE microphysics(dt,t_in,rho,qsat,wbar,dqsdz,est,cvi,qv_in,qc_in,qh_in,qi_in, &
     &                  t,qv,qc,qh,qi)
  ! qh = qr
  ! this code part is mostly from the original implementation and only changed to give back tendencies
  ! in order to work consistently with implicit time stepping
  ! the varialble names t,qv,qc,qh,qi might be somewhat missleading since in the end they are tendencies
  ! while they are states in the beginning of this routine
  IMPLICIT NONE
  REAL(wp), INTENT(IN)       :: &
       & dt,                    & ! time step
       & rho,                   & ! air density
       & qsat,                  & ! water vapor fraction at saturation
       & wbar,                  & ! vertical wind interpolated at cell center
       & dqsdz,                 & ! central vertical difference quotient of qsat
       & est,                   & ! saturation vapor pressure, em kPa
       & cvi                      ! ice ventilation coefficient
  REAL(wp), INTENT(IN)    :: &
       & t_in,                  & ! temperature
       & qv_in,                 & ! water vapor fraction
       & qc_in,                 & ! cloud water fraction
       & qh_in,                 & ! rain  water fraction
       & qi_in                    ! ice         fraction
  REAL(wp), INTENT(OUT)    :: &
       & t,                     & ! temperature for computations and tendency as output
       & qv,                    & ! water vapor fraction for computations and tendency as output
       & qc,                    & ! cloud water fraction for computations and tendency as output
       & qh,                    & ! rain  water fraction for computations and tendency as output
       & qi                       ! ice         fraction for computations and tendency as output

  t  =  t_in
  qv = qv_in
  qc = qc_in
  qh = qh_in
  qi = qi_in

  IF (qc <=1.0e-10_wp) qc = 0._wp  !DEFEAT UNDERFLOW PROBLEM
  IF (qh <=1.0e-10_wp) qh = 0._wp
  IF (qi <=1.0e-10_wp) qi = 0._wp

  CALL evaporate( dt =    dt, &
       &           t =     t, &
       &         rho =   rho, &
       &        qsat =  qsat, &
       &        wbar =  wbar, &
       &       dqsdz = dqsdz, &
       &         est =   est, &
       &         cvi =   cvi, &
       &          qv =    qv, &
       &          qc =    qc, &
       &          qh =    qh, &
       &          qi =    qi)    !vapor to cloud,cloud to vapor

  CALL sublimate( dt =   dt, &
       &           t =    t, &
       &         rho =  rho, &
       &         est =  est, &
       &        qsat = qsat, &
       &          qv =   qv, &
       &          qi =   qi)     !vapor to ice

  CALL glaciate( dt =   dt, &
       &          t =    t, &
       &         qv =   qv, &
       &         qh =   qh, &
       &         qi =   qi, &
       &       qsat = qsat)      !rain to ice

  CALL melt( dt =  dt, &
       &    rho = rho, &
       &    cvi = cvi, &
       &      t =   t, &
       &     qi =  qi, &
       &     qh =  qh)           !ice to rain

  CALL convert( dt =  dt, &
       &         t =   t, &
       &       rho = rho, &
       &        qc =  qc, &
       &        qh =  qh)        !(auto)conversion and accretion

  t  = ( t -  t_in)/dt
  qv = (qv - qv_in)/dt
  qc = (qc - qc_in)/dt
  qh = (qh - qh_in)/dt
  qi = (qi - qi_in)/dt

END SUBROUTINE microphysics


SUBROUTINE evaporate(dt,t,rho,qsat,wbar,dqsdz,est,cvi,qv,qc,qh,qi)
  !- evaporates cloud,rain and ice to saturation

  IMPLICIT NONE
  REAL(wp), INTENT(IN)       :: &
       & dt,                    & ! time step
       & rho,                   & ! air density
       & qsat,                  & ! water vapor fraction at saturation
       & wbar,                  & ! vertical wind interpolated at cell center
       & dqsdz,                 & ! central vertical difference quotient of qsat
       & est,                   & ! saturation vapor pressure, em kPa
       & cvi                      ! ice ventilation coefficient
  REAL(wp), INTENT(INOUT)    :: &
       & t,                     & ! temperature
       & qv,                    & ! water vapor fraction
       & qc,                    & ! cloud water fraction
       & qh,                    & ! rain  water fraction
       & qi                       ! ice         fraction
  !     XNO=10.0E06
  !     HERC = 1.93*1.E-6*XN035        !evaporation constant
  REAL(wp), PARAMETER        :: &
       & herc = 5.44e-4_wp,     & !
       & heatsubl = 2834._wp,   & !
       & heatcond = 2.501e6_wp, & !
       & frc = heatcond / cpd,  & !
       & src = heatsubl / cpd     !
  REAL(wp)                   :: &
       & evhdt,                 & !
       & evidt,                 & !
       & evrate,                & !
       & evap,                  & !
       & sd,                    & !
       & quant,                 & !
       & dividend,              & !
       & divisor,               & !
       & devidt                   !

  sd = qsat - qv  !vapor deficit
  IF (sd == 0.0_wp)  RETURN

  evhdt = 0._wp
  evidt = 0._wp
  ! print *, TRIM(thisroutine), '::wbar,dqsdz    :', wbar, dqsdz
  evrate = ABS (wbar * dqsdz)   ! evaporation rate (Kessler 8.32)
  evap = evrate * dt   ! what we can get in dt

  IF (SD <= 0.0_wp) THEN  ! condense. SD is negative
     IF (EVAP >= ABS (sd) ) THEN    !we get it all
        qc = qc   - sd  !deficit,remember?
        qv = qsat       !set the vapor to saturation
        t  = t    - sd * frc  !heat gained through condensation per gram of dry air
        RETURN
     ELSE
        qc = qc + evap         !get what we can in dt
        qv = qv - evap         !remove it from the vapor
        t  = t  + evap * frc   !get some heat
        RETURN
     ENDIF
  ELSE                               !SD is positive, need some water
     ! not saturated. saturate if possible. use everything in order
     ! cloud, rain, ice. SD is positive
     IF (evap <= qc) THEN        !enough cloud to last DT
        IF (sd <= evap) THEN          !enough time to saturate
           qc = qc   - sd       !remove cloud
           qv = qsat          !saturate
           t  = t    - sd * frc   !cool the parcel
           RETURN  !done
        ELSE   !not enough time
           sd = sd - evap               !use what there is
           qv = qv + evap     !add vapor
           t  = t  - evap * frc !lose heat
           qc = qc - evap     !lose cloud
           !go on to rain.
        ENDIF
     ELSE                !not enough cloud to last dt
        IF (sd <= qc) THEN   !but there is enough to sat
           qv = qsat  !use it
           qc = qc   - sd
           t  = t    - sd * frc
           RETURN
        ELSE            !not enough to sat
           sd = sd - qc
           qv = qv + qc
           t  = t  - qc * frc
           qc = 0.0_wp  !all gone
        ENDIF       !on to rain
     ENDIF          !finished with cloud
     !  but still not saturated, so try to use some rain
     !  this is tricky, because we only have time DT to evaporate. if there
     !  is enough rain, we can evaporate it for dt. ice can also sublimate
     !  at the same time. there is a compromise here.....use rain first, then
     !  ice. saturation may not be possible in one DT time.
     !  rain evaporation rate (W12),(OT25),(K Table 4). evaporate rain first
     !  sd is still positive or we wouldn't be here.

     IF (qh > 1.e-10_wp) THEN
        quant = (qsat - qc - qv) * rho   !g/m**3
        evhdt = (dt * herc * quant * (qh * rho) **.65_wp) / rho
        !             rain evaporation in time DT
        IF (evhdt <= qh) THEN           !enough rain to last dt
           IF (sd <= evhdt) THEN         !enough time to saturate
              qh = qh - sd    !remove rain
              qv = qsat       !saturate
              t  = t  - sd * frc    !cool the parcel
              !if(mintime>40) print*,'1',L,T(L)-273.15,QV(L)*1000,QH(L)*1000
              RETURN              !done
           ELSE                               !not enough time
              sd = sd - evhdt        !use what there is
              qv = qv + evhdt      !add vapor
              t  = t  - evhdt * frc      !lose heat
              qh = qh - evhdt      !lose rain
           ENDIF                   !go on to ice.

        ELSE  !not enough rain to last dt
           IF (sd <= qh) THEN           !but there is enough to sat
              qv = qsat               !use it
              qh = qh - sd
              t  = t  - sd * frc
              RETURN
           ELSE                              !not enough to sat
              sd = sd - qh
              qv = qv + qh
              t  = t  - qh * frc
              qh = 0.0_wp                !all gone
           ENDIF                             !on to ice
        ENDIF                                !finished with rain

     ENDIF

     IF (qi <= 1.e-10_wp) RETURN            !no ice there

     dividend = ( (1.e6_wp / rho) **0.475_wp) * (sd / qsat &
          - 1) * (qi **0.525_wp) * 1.13_wp
     divisor = 7.e5_wp + 4.1e6_wp / (10._wp * est)

     devidt = -cvi * dividend / divisor   !rate of change

     evidt = devidt * dt                      !what we could get

     ! logic here is identical to rain. could get fancy and make subroutine
     ! but duplication of code is easier. God bless the screen editor.

     IF (evidt <= qi) THEN               !enough ice to last DT
        IF (sd <= evidt) THEN                  !enough time to saturate
           qi = qi - sd                 !remove ice
           qv = qsat                    !saturate
           t  = t  - sd * src           !cool the parcel
           RETURN                                 !done
        ELSE                                   !not enough time
           sd = sd - evidt                        !use what there is
           qv = qv + evidt              !add vapor
           t  = t  - evidt * src        !lose heat
           qi = qi - evidt              !lose ice
        ENDIF                                  !go on,unsatisfied
     ELSE                                     !not enough ice to last DT
        IF (sd <= qi) THEN                !but there is enough to sat
           qv = qsat                    !use it
           qi = qi - sd
           t  = t  - sd * src
           RETURN
        ELSE                                   !not enough to sat
           sd = sd - qi
           qv = qv + qi
           t  = t  - qi * src
           qi = 0.0_wp                         !all gone
        ENDIF                                  !on to better things
        !finished with ice
     ENDIF
  ENDIF                                      !finished with the SD decision

END SUBROUTINE evaporate

SUBROUTINE convert(dt,t,rho,qc,qh)

  IMPLICIT NONE

  !- ACCRETION AND AUTOCONVERSION
  REAL(wp), INTENT(IN)    :: &
       & dt,                 & ! time step
       & t,                  & ! temperature
       & rho                   ! air density
  REAL(wp), INTENT(INOUT) :: &
       & qc,                 & ! cloud water fraction
       & qh                    ! rain water fraction
  INTEGER, PARAMETER      :: &
       &iconv = 1              !Kessler conversion
  REAL(wp), PARAMETER     :: &
       & ak1 = 0.001_wp,     & !conversion rate constant
       & ak2 = 0.0052_wp,    & !collection (accretion) rate
       & th  = 0.5_wp,       & !Kessler threshold
       & anbase =100000._wp, & !Berry-number at cloud base #/m^3(continental)
       & bdisp = 0.146_wp,   & !Berry--size dispersion (continental)
       & tfreeze = 269.3_wp    !ice formation temperature
  REAL(wp)                :: &
       & accrete,            & !
       & con,                & !
       & q,                  & !
       & h,                  & !
       & total                 !

  !     selection rules
  IF (t  <= tfreeze) RETURN  !process not allowed above ice
  IF (qc <= 0._wp  ) RETURN

  accrete = 0._wp
  con = 0._wp
  q = rho * qc
  h = rho * qh

  IF (qh > 0._wp) accrete = ak2 * q * (h**.875_wp)  !accretion, Kessler

  IF (iconv /= 0) THEN   !select Berry or Kessler
     con = q*q*q*bdisp/(60._wp*(5._wp*q*bdisp+0.0366_wp*anbase))
  ELSE
     con = MAX(0._wp,ak1 * (q - th)) ! versao otimizada
  ENDIF

  total = (con + accrete) * dt / rho

  IF (total < qc) THEN
     qc = qc - total
     qh = qh + total    !no phase change involved
     RETURN
  ELSE
     qh = qh + qc    !uses all there is
     qc = 0.0_wp
  ENDIF

END SUBROUTINE convert

SUBROUTINE sublimate(dt, t, rho, est, qsat, qv, qi)

  ! ********************* VAPOR TO ICE (USE EQUATION OT22)***************
  REAL(wp), INTENT(IN)      :: &
       & dt,                   & ! time step
       & rho,                  & ! air density
       & est,                  & !
       & qsat                    ! water vapor fraction at saturation
  REAL(wp), INTENT(INOUT)   :: &
       & t,                    & ! temperature
       & qv,                   & ! water vapor fraction
       & qi                      ! ice fraction
  REAL(wp), PARAMETER       :: &
       & heatsubl = 2834._wp,  & ! LM: from global parameters
       & src = heatsubl / cpd, & !
       & tfreeze = 269.3_wp      !
  REAL(wp) :: &
       & dtsubh,   & !
       & dividend, & !
       & divisor,  & !
       & subl        !

  dtsubh = 0._wp
  !selection criteria for sublimation
  IF (t  > tfreeze) RETURN
  IF (qv <= qsat   ) RETURN

  !     from (OT); correction factors for units applied
  dividend = ( (1.e6_wp / rho) **0.475_wp) * (qv / qsat &
       - 1._wp) * (qi **0.525_wp) * 1.13_wp
  divisor = 7.e5_wp + 4.1e6_wp / (10._wp * est)
  dtsubh = ABS (dividend / divisor)   !sublimation rate
  subl = dtsubh * dt                  !and amount possible

  !     again check the possibilities

  IF (subl < qv) THEN
     qv = qv  - subl             !lose vapor
     qi = qi  + subl             !gain ice
     t  = t   + subl * src       !energy change, warms air
     RETURN
  ELSE
     qi = qv                     !use what there is
     t  = t  + qv * src          !warm the air
     qv = 0.0_wp
  ENDIF

END SUBROUTINE sublimate

SUBROUTINE glaciate(dt,t,qv,qh,qi, qsat)

  ! *********************** CONVERSION OF RAIN TO ICE *******************
  !     uses equation OT 16, simplest. correction from W not applied, but
  !     vapor pressure differences are supplied.
  IMPLICIT NONE
  REAL(wp), INTENT(IN)     :: &
       & dt,                  & ! time step
       & qv,                  & ! water vapor fraction
       & qsat                   ! water vapor fraction at saturation
  REAL(wp), INTENT(INOUT)  :: &
       & t,                   & ! temperature
       & qh,                  & ! rain water fraction
       & qi                     ! ice fraction
  REAL(wp), PARAMETER      :: &
       & heatfus = 334._wp,   & !
       & frc = heatfus / cpd, & !
       & tfreeze =  269.3_wp, & !
       & glconst = 0.025_wp     !glaciation time constant, 1/sec
  REAL(wp) :: dfrzh

  dfrzh = 0._wp    !rate of mass gain in ice
  !selection rules for glaciation
  IF (qh <= 0._wp  ) RETURN
  IF (qv <  qsat   ) RETURN
  IF (t  >  tfreeze) RETURN

  dfrzh = dt * glconst * qh    ! from OT(16)
  IF (dfrzh < qh) THEN
     qi = qi + dfrzh
     qh = qh - dfrzh
     t  = t  + frc * dfrzh  !warms air
     RETURN
  ELSE
     qi = qi + qh
     t  = t  + frc * qh
     qh = 0.0_wp
  ENDIF

END SUBROUTINE glaciate

SUBROUTINE melt(dt,rho,cvi,t,qi,qh)
  ! ******************* MAKES WATER OUT OF ICE **************************
  REAL(wp), INTENT(IN)   :: &
       & dt,                & ! time step
       & rho,               & ! air density
       & cvi                  ! ice ventilation coefficient
  REAL(wp), INTENT(INOUT):: &
       & t,                 & ! temperature
       & qi,                & ! ice fraction
       & qh                   ! rain water fraction
  REAL(wp), PARAMETER    :: &
       & frc = 332.27_wp,   & !
       & tmelt = 273._wp,   & !
       & f0 = 0.75_wp         ! ice velocity factor
  REAL(wp) :: dtmelt          ! conversion,ice to rain

  dtmelt = 0._wp   !conversion,ice to rain
  !selection rules
  IF (qi <= 0.0_wp  ) RETURN
  IF (t  <  tmelt) RETURN

  dtmelt = dt * (2.27_wp / rho) * cvi * (t - tmelt) * ( (rho  &
       * qi * 1.e-6_wp) **0.525_wp) * (f0** ( - 0.42_wp) )
  !after Mason,1956
  !     check the possibilities
  IF (dtmelt < qi) THEN
     qh = qh + dtmelt
     qi = qi - dtmelt
     t  = t  - frc * dtmelt     !cools air
     RETURN
  ELSE
     qh = qh + qi   !get all there is to get
     t  = t  - frc * qi
     qi = 0.0_wp
  ENDIF

END SUBROUTINE melt

PURE ELEMENTAL FUNCTION psat_water( T ) RESULT ( psatw )
  ! psat_water is taken from McSnow code to replace the old formutaltion
  REAL(wp), INTENT(in) :: T
  REAL(wp)             :: psatw
  REAL(wp), PARAMETER ::    &
       A_w  = 1.72693882e1_wp, &  !..constant in saturation pressure - liquid
       B_w  = 3.58600000e1_wp, &  !..constant in saturation pressure - liquid
       T_3  = 273.15_wp,       &  ! [K]     melting temperature of ice/snow
       e_3  = 6.10780000e2_wp     !..saturation pressure at T = T_3
  ! Beheng, lecture scrpt
  ! psatw = 610.7_wp * EXP( 17.15_wp * Tc / (Tc + 234.9_wp) )
  psatw = e_3 * EXP( A_w * (T - T_3) / (T - B_w) )
  ! alternative Goff-Gratch equation
  ! tts = 373.15_wp / atmo%T
  ! psatw = 101325_wp * 10**(-7.90298_wp*(tts-1._wp) + 5.02808_wp*LOG10(tts) - &
  !  1.3816e-7_wp*(10**(11.344*(1._wp-1._wp/tts))-1._wp) + &
  !  8.1328e-3_wp*(10**(3.49149_wp*(1._wp-tts)-1._wp)))
END FUNCTION psat_water

END MODULE mo_art_emission_BBPlume
