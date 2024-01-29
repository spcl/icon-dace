!
! Module containing thermodynamic functions used by the AES department in MPI-M
!
! Code slightly refactored from gscp_graupel (DWD)
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

MODULE mo_aes_graupel

USE mo_kind,               ONLY: wp         , &
                                 i4

USE mo_physical_constants, ONLY: alv   , & !! latent heat of vapourization
                                 als   , & !! latent heat of sublimation
                                 rv    , & !! vapor gas constant
                                 rd    , & !! vapor gas constant
                                 cpv   , & !! isobaric specific heat of vapor
                                 cvd   , & !! isometric specific heat of dry air
                                 cvv   , & !! isometric specific heat of vapor
                                 clw   , & !! specific heat of liquid water
                                 alv   , & !! invariant part of vaporization enthalpy
                                 als   , & !! invariant part of sublimation enthalpy
                                 tmelt     !! melting temperature of ice/snow

USE mo_aes_thermo, ONLY:        & !! 
        qsat_rho               ,&
        qsat_ice_rho           ,&
        internal_energy        ,&
        T_from_internal_energy

IMPLICIT NONE
PRIVATE
PUBLIC :: graupel, snow_number, snow_lambda, ice_number

LOGICAL, PARAMETER :: &
  lrain        = .TRUE.  , & ! switch for disabling rain
  lcold        = .TRUE.      ! switch for disabling freezing processes

INTEGER, PARAMETER :: &
  nx  = 6 , & ! number of water species
  np  = 4 , & ! number of precipitating water species
  lqr = 1 , & ! index for rain
  lqi = 2 , & ! index for ice
  lqs = 3 , & ! index for snow
  lqg = 4 , & ! index for graupel
  lqc = 5 , & ! index for cloud
  lqv = 6     ! index for vapor

INTEGER, PARAMETER ::                          &
  qx_ind(nx) = [lqv, lqc, lqr, lqs, lqi, lqg] ,& !
  qp_ind(np) = [lqr, lqi, lqs, lqg]

REAL(wp), PARAMETER, DIMENSION(3,np) :: &
   params = RESHAPE([ 14.58_wp,  0.111_wp     , 1.e-12_wp, &
                       1.25_wp,  0.160_wp     , 1.e-12_wp, &
                      57.80_wp,  0.5_wp/3.0_wp, 1.e-12_wp, &
                      12.24_wp,  0.217_wp     , 1.e-08_wp] ,[3,np])

#ifdef __INLINE_RESHAPE_WAR
! Workaround for buggy -Minline=reshape
REAL(wp), PARAMETER, DIMENSION(3) :: params_qr = params(:, lqr)
REAL(wp), PARAMETER, DIMENSION(3) :: params_qi = params(:, lqi)
REAL(wp), PARAMETER, DIMENSION(3) :: params_qs = params(:, lqs)
REAL(wp), PARAMETER, DIMENSION(3) :: params_qg = params(:, lqg)
#endif

REAL(wp), PARAMETER :: &
   rho_00 = 1.225_wp        , & ! reference air density
   q1     = 8.e-6_wp        , &               
   qmin   = 1.0E-15_wp      , & ! threshold for computation
   ams    = 0.069_wp        , & ! Formfactor in the mass-size relation of snow particles 
   bms    = 2.0_wp          , & ! Exponent in the mass-size relation of snow particles
   v0s    = 25.0_wp         , & ! prefactor in snow fall speed
   v1s    = 0.5_wp          , & ! Exponent in the terminal velocity for snow
   m0_ice = 1.0E-12_wp      , & ! initial crystal mass for cloud ice nucleation
   ci     = 2108._wp        , & ! specific heat of ice
   tx     = 3339.5_wp       , & !
   tfrz_het1 = tmelt- 6.0_wp, & ! temperature for het. freezing of cloud water with supersat
   tfrz_het2 = tmelt-25.0_wp, & ! temperature for het. freezing of cloud water
   tfrz_hom  = tmelt-37.0_wp, & ! temperature for hom. freezing of cloud water
   lvc = alv-(cpv-clw)*tmelt, & ! invariant part of vaporization enthalpy
   lsc = als-(cpv-ci )*tmelt    ! invariant part of vaporization enthalpy

TYPE t_qx_ptr                   ! type for pointer vector
  REAL(wp), POINTER    :: p(:), x(:,:)
END TYPE t_qx_ptr
  
CONTAINS

  SUBROUTINE graupel(nvec, ke, ivstart, ivend, kstart,        & !! start/end indicies
             dt, dz, t, p, rho, qv, qc, qi, qr, qs, qg, qnc,  & !! prognostic variables
             prr_gsp, pri_gsp, prs_gsp, prg_gsp, pflx        )  !  total precipitation flux

  INTEGER, INTENT(IN) ::  &
    nvec      , & !> number of horizontal points
    ke        , & !! number of grid points in vertical direction
    ivstart   , & !! start index for horizontal direction
    ivend     , & !! end index   for horizontal direction
    kstart        !! start index for the vertical 

  REAL(KIND=wp), INTENT(IN) :: &
    dt            !> time step for integration of microphysics   (  s  )

  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
    dz        , & !> layer thickness of full levels                (  m  )
    rho       , & !! density of moist air                          (kg/m3)
    p             !! pressure                                      ( Pa  )

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    t             !> temperature                                   (  K  )

  REAL(KIND=wp), TARGET, DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    qv        , & !! specific water vapor content                  (kg/kg)
    qc        , & !! specific cloud water content                  (kg/kg)
    qi        , & !! specific cloud ice   content                  (kg/kg)
    qr        , & !! specific rain content                         (kg/kg)
    qs        , & !! specific snow content                         (kg/kg)
    qg            !! specific graupel content                      (kg/kg)

  REAL(KIND=wp), DIMENSION(:), INTENT(IN) ::   &   ! dim (ie)
    qnc           !! cloud number concentration

  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT) ::   &   ! dim (ie,ke)
    pflx          ! total precipitation flux 

  REAL(KIND=wp), TARGET, DIMENSION(:), INTENT(OUT) ::   &   ! dim (ie)
    prr_gsp   , & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    pri_gsp   , & !> precipitation rate of ice, grid-scale         (kg/(m2*s))
    prs_gsp   , & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
    prg_gsp       !! precipitation rate of graupel, grid-scale     (kg/(m2*s))

  LOGICAL :: is_sig_present(nvec*ke) ! is snow, ice or graupel present? 

  INTEGER (KIND=i4)  :: iv, k, kp1, j, jmx, jmx_, ix, iqx, &   !> loop indices
                        ind_k(nvec*ke), & ! k index of gathred point
                        ind_i(nvec*ke), & ! iv index of gathered point
                        kmin(nvec,np)     ! first level with condensate

  REAL(KIND=wp)  :: cv, vc, eta, zeta, qvsi, qice, qliq,  qtot, dvsw, dvsw0, dvsi ,&
                    n_ice, m_ice, x_ice,n_snow,l_snow, ice_dep, e_int, stot, xrho

  REAL(KIND=wp) ::   &
    update(3)       ,& !> scratch array with output from precipitation step
    sink(nx)        ,& !! tendencies
    dqdt(nx)        ,& !! tendencies
    sx2x(nx,nx)     ,& !! conversion rates
    eflx(nvec)      ,& !! internal energy flux from precipitation      (W/m2 )
    vt(nvec,np)        !! terminal velocity for for different hydrometeor categories

  TYPE(t_qx_ptr) :: q(nx) ! vector of pointers to point to four hydrometeor inouts

  q(lqr)%x => qr(:,:); q(lqr)%p => prr_gsp(:)
  q(lqi)%x => qi(:,:); q(lqi)%p => pri_gsp(:)
  q(lqs)%x => qs(:,:); q(lqs)%p => prs_gsp(:)
  q(lqg)%x => qg(:,:); q(lqg)%p => prg_gsp(:)
  q(lqc)%x => qc(:,:)
  q(lqv)%x => qv(:,:)

  jmx = 0
  jmx_ = jmx

  !$ACC DATA &
  !$ACC   PRESENT(dz, t, p, rho, qv, qc, qi, qr, qs, qg, qnc) &
  !$ACC   PRESENT(prr_gsp, prs_gsp, pri_gsp, prg_gsp, pflx) &
  !$ACC   COPYIN(jmx_) &
  !$ACC   CREATE(is_sig_present, ind_k, ind_i, kmin, eflx, vt)

  !$ACC ENTER DATA COPYIN(q(1:nx))
  DO ix=1,nx
    !$ACC ENTER DATA COPYIN(q(ix)%x)
  END DO
  DO ix=1,np
    !$ACC ENTER DATA COPYIN(q(ix)%p)
  END DO

  !jmx=0
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP SEQ
  DO  k = ke,kstart,-1
    !$ACC LOOP GANG VECTOR PRIVATE(jmx, iqx)
    DO iv = ivstart, ivend
      IF ( (MAX(q(lqc)%x(iv,k),q(lqr)%x(iv,k),q(lqs)%x(iv,k),q(lqi)%x(iv,k),q(lqg)%x(iv,k)) > qmin)  &
              .OR. (t(iv,k)<tfrz_het2 .AND. q(lqv)%x(iv,k)> qsat_ice_rho(t(iv,k),rho(iv,k))) ) THEN
         !$ACC ATOMIC CAPTURE
         jmx_ = jmx_ + 1
         jmx  = jmx_
         !$ACC END ATOMIC
         ind_k(jmx) = k
         ind_i(jmx) = iv
         is_sig_present(jmx) = MAX(q(lqs)%x(iv,k),q(lqi)%x(iv,k),q(lqg)%x(iv,k)) > qmin
      ENDIF
      !$ACC LOOP SEQ
      DO ix=1,np
        iqx = qp_ind(ix)
        IF (k == ke) THEN
          kmin(iv,iqx) = ke+1
          q(iqx)%p(iv) = 0.0_wp
          vt(iv,ix)    = 0.0_wp
        ENDIF
        IF (q(iqx)%x(iv,k)>qmin) kmin(iv,iqx) = k
      END DO
    END DO
  END DO
  !$ACC END PARALLEL

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR &
  !$ACC   PRIVATE(sink, dqdt, sx2x) &
  !$ACC   PRIVATE(k, iv) &
  !$ACC   PRIVATE(dvsw, qvsi, dvsi, n_snow, l_snow) &
  !$ACC   PRIVATE(n_ice, m_ice, x_ice, eta, ice_dep) &
  !$ACC   PRIVATE(dvsw0, ix, iqx, stot, qice, qliq, qtot, cv)
  DO j=1,jmx_
    k  = ind_k(j)
    iv = ind_i(j)

    dvsw   = q(lqv)%x(iv,k)-qsat_rho(t(iv,k),rho(iv,k))
    qvsi   = qsat_ice_rho(t(iv,k),rho(iv,k))
    dvsi   = q(lqv)%x(iv,k)-qvsi
    n_snow = snow_number  (t(iv,k),rho(iv,k),q(lqs)%x(iv,k))
    l_snow = snow_lambda  (rho(iv,k),q(lqs)%x(iv,k),n_snow)

    sx2x(:,:)     = 0.0_wp
    sx2x(lqc,lqr) = cloud_to_rain   (t(iv,k),q(lqc)%x(iv,k),q(lqr)%x(iv,k),qnc(1))
    sx2x(lqr,lqv) = rain_to_vapor   (t(iv,k),rho(iv,k),q(lqc)%x(iv,k),q(lqr)%x(iv,k),dvsw,dt)
    sx2x(lqc,lqi) = cloud_x_ice (t(iv,k),q(lqc)%x(iv,k),q(lqi)%x(iv,k),dt)
    sx2x(lqi,lqc) = -MIN(sx2x(lqc,lqi),0.0_wp)
    sx2x(lqc,lqi) =  MAX(sx2x(lqc,lqi),0.0_wp)
    sx2x(lqc,lqs) = cloud_to_snow (t(iv,k),q(lqc)%x(iv,k),q(lqs)%x(iv,k),n_snow,l_snow )
    sx2x(lqc,lqg) = cloud_to_graupel (t(iv,k),rho(iv,k),q(lqc)%x(iv,k),q(lqg)%x(iv,k))

    IF (t(iv,k)<tmelt) THEN
      n_ice   = ice_number   (t(iv,k),rho(iv,k))
      m_ice   = ice_mass     (q(lqi)%x(iv,k),n_ice)
      x_ice   = ice_sticking (t(iv,k))
    
      IF (is_sig_present(j)) THEN
        eta           = deposition_factor(t(iv,k),qvsi) ! neglect cloud depth cor. from gcsp_graupel
        sx2x(lqv,lqi) = vapor_x_ice    (q(lqi)%x(iv,k),m_ice,eta,dvsi,dt)
        sx2x(lqi,lqv) =-MIN(sx2x(lqv,lqi),0.0_wp)
        sx2x(lqv,lqi) = MAX(sx2x(lqv,lqi),0.0_wp)
        ice_dep       = MIN(sx2x(lqv,lqi),dvsi/dt)

        sx2x(lqi,lqs) = deposition_auto_conversion(q(lqi)%x(iv,k),m_ice,ice_dep)
        sx2x(lqi,lqs) = sx2x(lqi,lqs) + ice_to_snow(q(lqi)%x(iv,k),n_snow,l_snow,x_ice)
        sx2x(lqi,lqg) = ice_to_graupel  (rho(iv,k),q(lqr)%x(iv,k),q(lqg)%x(iv,k),q(lqi)%x(iv,k),x_ice)
        sx2x(lqs,lqg) = snow_to_graupel (t(iv,k),rho(iv,k),q(lqc)%x(iv,k),q(lqs)%x(iv,k))
        sx2x(lqr,lqg) = rain_to_graupel(t(iv,k),rho(iv,k),q(lqc)%x(iv,k),q(lqr)%x(iv,k),q(lqi)%x(iv,k),q(lqs)%x(iv,k),m_ice,dvsw,dt)
      ENDIF
      sx2x(lqv,lqi)   = sx2x(lqv,lqi) + ice_deposition_nucleation (t(iv,k),q(lqc)%x(iv,k),q(lqi)%x(iv,k),n_ice,dvsi,dt)
    ELSE
      sx2x(lqc,lqr)   = sx2x(lqc,lqr) + sx2x(lqc,lqs) + sx2x(lqc,lqg)
      sx2x(lqc,lqs)   = 0.0_wp
      sx2x(lqc,lqg)   = 0.0_wp
      ice_dep         = 0.0_wp
      eta             = 0.0_wp
    ENDIF

    IF (is_sig_present(j)) THEN
      dvsw0         = q(lqv)%x(iv,k)-qsat_rho(tmelt,rho(iv,k))
      sx2x(lqv,lqs) = vapor_x_snow(t(iv,k),p(iv,k),rho(iv,k),q(lqs)%x(iv,k),n_snow,l_snow,eta,ice_dep,dvsw,dvsi,dvsw0,dt)
      sx2x(lqs,lqv) = -MIN(sx2x(lqv,lqs),0.0_wp)
      sx2x(lqv,lqs) =  MAX(sx2x(lqv,lqs),0.0_wp)
      sx2x(lqv,lqg) = vapor_x_graupel(t(iv,k),p(iv,k),rho(iv,k),q(lqg)%x(iv,k),dvsw,dvsi,dvsw0,dt)
      sx2x(lqg,lqv) = -MIN(sx2x(lqv,lqg),0.0_wp)
      sx2x(lqv,lqg) =  MAX(sx2x(lqv,lqg),0.0_wp)
      sx2x(lqs,lqr) = snow_to_rain    (t(iv,k),p(iv,k),rho(iv,k),dvsw0,q(lqs)%x(iv,k))
      sx2x(lqg,lqr) = graupel_to_rain (t(iv,k),p(iv,k),rho(iv,k),dvsw0,q(lqg)%x(iv,k))
    ENDIF

    !$ACC LOOP SEQ
    DO ix=1,nx
      iqx = qx_ind(ix)
      sink(iqx) = 0.0_wp
      IF (is_sig_present(j) .OR. iqx == lqc .OR. iqx==lqv .OR. iqx==lqr) THEN
        sink(iqx) = SUM(sx2x(iqx,:))
        stot = q(iqx)%x(iv,k) / dt
        IF (sink(iqx) > stot .AND. q(iqx)%x(iv,k)>qmin) THEN
          sx2x(iqx,:) = sx2x(iqx,:) * stot/sink(iqx)
          sink(iqx)   = SUM(sx2x(iqx,:))
        ENDIF
      ENDIF 
    END DO

    !$ACC LOOP SEQ
    DO ix=1,nx
      iqx = qx_ind(ix)
      dqdt(iqx)      = SUM(sx2x(:,iqx)) - sink(iqx)
      q(iqx)%x(iv,k) = MAX(0.0_wp, q(iqx)%x(iv,k) + dqdt(iqx)*dt)
    END DO
      
    qice     = q(lqs)%x(iv,k) + q(lqi)%x(iv,k) + q(lqg)%x(iv,k)
    qliq     = q(lqc)%x(iv,k) + q(lqr)%x(iv,k)
    qtot     = q(lqv)%x(iv,k)  + qice + qliq
    cv       = cvd + (cvv-cvd)*qtot + (clw-cvv)*qliq + (ci-cvv)*qice ! qtot? or qv?
    t(iv,k)  = t(iv,k) + dt *  ( (dqdt(lqc) + dqdt(lqr)) * (lvc-(clw-cvv)*t(iv,k))  &
             + (dqdt(lqi) + dqdt(lqs) + dqdt(lqg)) * (lsc-(ci -cvv)*t(iv,k)) )/cv
  END DO 
  !$ACC END PARALLEL

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP SEQ
  DO  k = kstart, MERGE(ke,kstart-1,lrain)  
    !$ACC LOOP GANG VECTOR &
    !$ACC   PRIVATE(kp1, qliq, qice, e_int, zeta, xrho, vc) &
    !$ACC   PRIVATE(ix, iqx, update)
    DO iv = ivstart, ivend  
      IF (k==kstart) eflx(iv) = 0.0_wp
      kp1 = min(ke    ,k+1)

      IF( k >= MINVAL(kmin(iv,:)) )THEN
        qliq  = q(lqc)%x(iv,k) + q(lqr)%x(iv,k)
        qice  = q(lqs)%x(iv,k) + q(lqi)%x(iv,k) + q(lqg)%x(iv,k)
        e_int = internal_energy(t(iv,k),qv(iv,k),qliq,qice,rho(iv,k),dz(iv,k) ) + eflx(iv)

        zeta  = dt/(2.0_wp*dz(iv,k))
        xrho  = SQRT(rho_00/rho(iv,k))

#ifdef __INLINE_RESHAPE_WAR
        IF (k >= kmin(iv,lqr)) THEN
          vc     = vel_scale_factor(lqr, xrho, rho(iv,k), t(iv,k), q(lqr)%x(iv,k))
          update = precip(params_qr,zeta,vc,q(lqr)%p(iv),vt(iv,1),q(lqr)%x(iv,k),q(lqr)%x(iv,kp1),rho(iv,k))
          q(lqr)%x(iv,k) = update(1); q(lqr)%p(iv) = update(2); vt(iv,1) = update(3)
        END IF
        IF (k >= kmin(iv,lqi)) THEN
          vc     = vel_scale_factor(lqi, xrho, rho(iv,k), t(iv,k), q(lqi)%x(iv,k))
          update = precip(params_qi,zeta,vc,q(lqi)%p(iv),vt(iv,2),q(lqi)%x(iv,k),q(lqi)%x(iv,kp1),rho(iv,k))
          q(lqi)%x(iv,k) = update(1); q(lqi)%p(iv) = update(2); vt(iv,2) = update(3)
        END IF
        IF (k >= kmin(iv,lqs)) THEN
          vc     = vel_scale_factor(lqs, xrho, rho(iv,k), t(iv,k), q(lqs)%x(iv,k))
          update = precip(params_qs,zeta,vc,q(lqs)%p(iv),vt(iv,3),q(lqs)%x(iv,k),q(lqs)%x(iv,kp1),rho(iv,k))
          q(lqs)%x(iv,k) = update(1); q(lqs)%p(iv) = update(2); vt(iv,3) = update(3)
        END IF
        IF (k >= kmin(iv,lqg)) THEN
          vc     = vel_scale_factor(lqg, xrho, rho(iv,k), t(iv,k), q(lqg)%x(iv,k))
          update = precip(params_qg,zeta,vc,q(lqg)%p(iv),vt(iv,4),q(lqg)%x(iv,k),q(lqg)%x(iv,kp1),rho(iv,k))
          q(lqg)%x(iv,k) = update(1); q(lqg)%p(iv) = update(2); vt(iv,4) = update(3)
        END IF
#else
        !$ACC LOOP SEQ
        DO ix=1,np
          iqx = qp_ind(ix)
          IF (k >= kmin(iv,iqx)) THEN
            vc     = vel_scale_factor(iqx, xrho, rho(iv,k), t(iv,k), q(iqx)%x(iv,k))
            update = precip(params(:,iqx),zeta,vc,q(iqx)%p(iv),vt(iv,ix),q(iqx)%x(iv,k),q(iqx)%x(iv,kp1),rho(iv,k))
            q(iqx)%x(iv,k) = update(1); q(iqx)%p(iv) = update(2); vt(iv,ix) = update(3)
          END IF
        END DO
#endif

        pflx(iv,k) = q(lqs)%p(iv) + q(lqi)%p(iv) + q(lqg)%p(iv)
        eflx(iv)   = dt*( q(lqr)%p(iv) * (clw*t(iv,k)-cvd*t(iv,kp1) - lvc)  &
                 &      +  pflx(iv,k)  * (ci *t(iv,k)-cvd*t(iv,kp1) - lsc))
        pflx(iv,k) = pflx(iv,k) + q(lqr)%p(iv)
        qliq       = q(lqc)%x(iv,k) + q(lqr)%x(iv,k)
        qice       = q(lqs)%x(iv,k) + q(lqi)%x(iv,k) + q(lqg)%x(iv,k)
        e_int      = e_int - eflx(iv)
        t(iv,k)    = T_from_internal_energy(e_int,qv(iv,k),qliq,qice,rho(iv,k),dz(iv,k))
      ENDIF
    END DO
  END DO
  !$ACC END PARALLEL

  !$ACC WAIT

  DO ix=1,nx
    !$ACC EXIT DATA COPYOUT(q(ix)%x)
  END DO
  DO ix=1,np
    !$ACC EXIT DATA COPYOUT(q(ix)%p)
  END DO
  !$ACC EXIT DATA COPYOUT(q(1:nx))

  !$ACC END DATA

END SUBROUTINE graupel

!!!=============================================================================================

PURE FUNCTION precip(params,zeta,vc,flx,vt,q,q_kp1,rho)
       
  REAL(KIND=wp) precip(3)       !> time step for integration of microphysics  (  s  )
  REAL(KIND=wp), INTENT(IN) :: &
    params(3) , &               !> fall speed parameters
    zeta      , &               !> dt/(2dz)
    vc        , &               !> state dependent fall speed correction
    flx       , &               !> flux into cell from above
    vt        , &               !> terminal velocity
    q         , &               !> specific mass of hydrometeor
    q_kp1     , &               !> specific mass in next lower cell
    rho                         !> density

  REAL(KIND=wp) :: rho_x, flx_eff, flx_partial

  !$ACC ROUTINE SEQ
  rho_x        = q*rho
  flx_eff      = rho_x/zeta + 2.0_wp*flx
  flx_partial  = rho_x * vc * fall_speed(rho_x, params) 
  flx_partial  = MIN( flx_partial, flx_eff )
  precip(1)    = zeta*(flx_eff-flx_partial) / ((1.0_wp + zeta*vt)*rho)  ! q update
  precip(2)    = (precip(1)*rho*vt + flx_partial)*0.5                   ! flx
  rho_x        = (precip(1)+q_kp1)*0.5_wp*rho
  precip(3)    = vc*fall_speed(rho_x, params)                           ! vt

END FUNCTION precip

!!!=============================================================================================

PURE FUNCTION vel_scale_factor(iqx,xrho,rho,t,qx)

  REAL(KIND=wp)                   :: vel_scale_factor
  INTEGER (KIND=i4), INTENT(IN)   :: iqx
  REAL(KIND=wp), INTENT(IN)       :: &
         xrho  , & ! sqrt(rho_00/rho)
         rho   , & ! density of condensate
         t     , & ! temperature
         qx        ! specific mass

  REAL (KIND=wp), PARAMETER ::  &
    b_i    =  2.0_wp/3.0_wp   , &
    b_s    = -1.0_wp/6.0_wp
   
   !$ACC ROUTINE SEQ
   SELECT CASE(iqx)
   CASE (lqi)
      vel_scale_factor = xrho**b_i
   CASE (lqs)
      vel_scale_factor = xrho * snow_number(t, rho, qx)**b_s
   CASE DEFAULT
      vel_scale_factor = xrho
   END SELECT

END FUNCTION vel_scale_factor

!!!=============================================================================================

PURE FUNCTION fall_speed(density, params)
  REAL(KIND=wp)             :: fall_speed
  REAL(KIND=wp), INTENT(IN) :: density   , & ! density of condensate
          &                    params(3)     ! fall speed parameters    

  !$ACC ROUTINE SEQ
  fall_speed  =  params(1) * ((density+params(3)) ** params(2))

END FUNCTION fall_speed

!!!=============================================================================================

PURE FUNCTION snow_number(t,rho,qs)
  REAL(KIND=wp)             :: snow_number
  REAL(KIND=wp), INTENT(IN) :: t     , & ! temperature
          &                    rho   , & ! ambient air density
          &                    qs        ! snow  specific mass

  REAL (KIND=wp),     PARAMETER ::   &
    tmin    =  tmelt-40._wp        , &
    tmax    =  tmelt               , &
    qsmin   =  2.0e-6_wp           , &
    xa1     = -1.65e+0_wp          , &
    xa2     =  5.45e-2_wp          , &
    xa3     =  3.27e-4_wp          , &
    xb1     =  1.42e+0_wp          , &
    xb2     =  1.19e-2_wp          , &
    xb3     =  9.60e-5_wp          , &
    n0s0    =  8.00e+5_wp          , &
    n0s1    =  13.5_wp* 5.65e+05_wp, &
    n0s2    = -0.107_wp            , &
    n0s3    =  13.5_wp             , &
    n0s4    =  0.5_wp*n0s1         , &
    n0s5    =  1.e6_wp             , &
    n0s6    =  1.e2_wp*n0s1        , &
    n0s7    =  1.e9_wp       

  REAL(KIND=wp) ::                   &
    tc                             , &
    alf                            , &
    bet                            , &
    n0s                            , &
    y                              , &
    n0smn                          , &
    n0smx                          

    !$ACC ROUTINE SEQ
    IF (qs > qmin) THEN
      tc    = MAX(MIN(t,tmax),tmin) - tmelt
      alf   = 10.**( xa1 + tc*(xa2 + tc*xa3) )
      bet   = xb1 + tc*(xb2 + tc*xb3)
      n0s   = n0s3 * ((qs+qsmin) * rho / ams )**(4-3*bet) / (alf*alf*alf)

      y     = EXP(n0s2*tc)
      n0smn = MAX(n0s4*y,n0s5)
      n0smx = MIN(n0s6*y,n0s7)
      snow_number = MIN(n0smx,MAX(n0smn,n0s))
    ELSE
      snow_number = n0s0
    ENDIF 

END FUNCTION snow_number

!!!=============================================================================================

PURE FUNCTION snow_lambda(rho,qs,ns)
  REAL(KIND=wp)             :: snow_lambda     ! returns riming snow rate
  REAL(KIND=wp), INTENT(IN) :: rho         , & ! ambient density
          &                    qs          , & ! snow specific mass
          &                    ns              ! snow number

  REAL(KIND=wp), PARAMETER  ::        &
      !a1       = ams/bms            , & ! -- used constants in expression 
      a2       = ams*2.0_wp         , & ! '' (with ams*gam(bms+1.0_wp) where gam(3) = 2) 
      lmd_0    = 1.0e+10_wp         , & ! no snow value of lambda
      bx       = 1.0_wp/(bms+1.0_wp), & ! ''
      qsmin    = 0.0e-6_wp ! previous had 2.0e-6         

  !$ACC ROUTINE SEQ
  IF (qs > qmin) THEN
    snow_lambda =  (a2*ns/((qs+qsmin)*rho)) ** bx
  ELSE
    snow_lambda =  lmd_0
  ENDIF
      
END FUNCTION snow_lambda

!!!=============================================================================================

PURE FUNCTION ice_number(t,rho)
  REAL(KIND=wp)             :: ice_number    ! ice number following cooper
  REAL(KIND=wp), INTENT(IN) :: t          , &! ambient temperature, in kelvin
                               rho           ! ambient density

  REAL (KIND=wp), PARAMETER ::    &
    a     = 5.000_wp            , & ! parameter in cooper fit
    b     = 0.304_wp            , & ! parameter in cooper fit
    nimax = 250.E+3_wp              ! maximal number of ice crystals 

  !$ACC ROUTINE SEQ
  ice_number =  MIN(nimax,a * EXP(b * (tmelt - t))) / rho

END FUNCTION ice_number

!!!=============================================================================================

PURE FUNCTION ice_mass(qi,ni)
  REAL(KIND=wp)             :: ice_mass   ! conversion rate of ice to snow
  REAL(KIND=wp), INTENT(IN) :: qi             , & ! ice specific mass 
          &                    ni                 ! ice crystal number

  REAL(KIND=wp), PARAMETER  :: &
      mi_max = 1.0E-09_wp        ! maximum mass of cloud ice crystals

  !$ACC ROUTINE SEQ
  ice_mass = MAX( m0_ice, MIN( qi/ni, mi_max ))
      
END FUNCTION ice_mass

!!!=============================================================================================

PURE FUNCTION ice_sticking(t)
  REAL(KIND=wp)             :: ice_sticking ! returns sticking efficiency of ic3
  REAL(KIND=wp), INTENT(IN) :: t            ! temperature

  REAL(KIND=wp), PARAMETER  ::      &
      a       = 0.09_wp           , & ! scale factor for freezing depression
      b       = 1.00_wp           , & ! maximum for exponential temperature factor
      eff_min = 0.075_wp          , & ! minimum sticking efficiency
      eff_fac = 3.5E-3_wp         , & ! Scaling factor [1/K] for cloud ice sticking efficiency
      tcrit   = tmelt-85._wp          ! Temperature at which cloud ice autoconversion starts

  ! per original code seems like aggregation is allowed even with no snow present
  !
  !$ACC ROUTINE SEQ
  ice_sticking = MAX(MIN(EXP(a*(t-tmelt)),b), eff_min, eff_fac*(t-tcrit)) 
      
END FUNCTION ice_sticking

!!!=============================================================================================

PURE FUNCTION deposition_factor(t,qvsi)
  REAL(KIND=wp)             :: deposition_factor  ! returns deposition factor 
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    qvsi            ! saturation (ice) specific vapor mass

  REAL(KIND=wp), PARAMETER  ::      &
      kappa  = 2.40E-2_wp         , & ! thermal conductivity of dry air
      b      = 1.94_wp            , & 
      a      = als*als/(kappa*rv) , &
      cx     = 2.22E-5_wp * tmelt**(-b) * 101325.0_wp

  REAL(KIND=wp) :: x
      
  !$ACC ROUTINE SEQ
  x =  cx/rd * t**(b-1.0_wp)
  deposition_factor    = x / (1.0_wp + a * x * qvsi/(t*t))
      
END FUNCTION deposition_factor

!!!=============================================================================================

PURE FUNCTION cloud_to_rain(t,qc,qr,nc)
  REAL(KIND=wp)             :: cloud_to_rain      ! mass from qc to qr 
  REAL(KIND=wp), INTENT(IN) :: t         , & ! temperature
          &                    qc        , & ! cloud water specific mass
          &                    qr        , & ! rain water specific mass
          &                    nc            ! cloud water number concentration 

  REAL(KIND=wp), PARAMETER  ::  &
          qmin_ac    = 1.00e-06_wp,  & ! threshold for auto conversion
          tau_max    = 0.90e+00_wp,  & ! maximum allowed value of tau
          tau_min    = 1.00e-30_wp,  & ! maximum allowed value of tau
          a          = 6.00e+02_wp,  & ! constant in phi-function for autoconversion
          b          = 0.68e+00_wp,  & ! exponent in phi-function for autoconversion
          c          = 5.00e-05_wp,  & ! exponent in phi-function for accretion
          ac_kernel  = 5.25e+00_wp,  & ! kernel coeff for SB2001 accretion
          x3         = 2.00e+00_wp,  & ! gamma exponent for cloud distribution
          x2         = 2.60e-10_wp,  & ! separating mass between cloud and rain
          x1         = 9.44e+09_wp,  & ! kernel coeff for SB2001 autoconversion
          au_kernel  = x1 / (20.0_wp*x2) * (x3+2.0_wp)*(x3+4.0_wp)/(x3+1.0_wp)**2

  REAL(KIND=wp) :: tau, & ! time-scale
                   phi, & ! similarity function for autoconversion
                   xau, & ! autoconversion rate
                   xac    ! accretion rate
    !
    ! Kessler (1969) autoconversion rate
    !    scau = zccau * MAX( qc_ik - qc0, 0.0_wp )
    !    scac = zcac  * qc_ik * zeln7o8qrk
    ! 
    ! Seifert and Beheng (2001) autoconversion rate
    ! with constant cloud droplet number concentration qnc
    !
    !$ACC ROUTINE SEQ
    cloud_to_rain = 0.0_wp
    IF (qc > qmin_ac .AND. t > tfrz_hom) THEN
      tau  = MAX(tau_min,MIN(1.0_wp-qc/(qc+qr),tau_max))
      phi  = tau**b
      phi  = a * phi * (1.0_wp - phi)**3
      xau  = au_kernel * (qc*qc/nc)**2  * (1.0_wp + phi/(1.0_wp - tau)**2)
      xac  = ac_kernel * qc * qr * (tau/(tau+c))**4
      cloud_to_rain = xau + xac
    ENDIF 

END FUNCTION cloud_to_rain

!!!=============================================================================================

PURE FUNCTION cloud_x_ice(t,qc,qi,dt)
  REAL(KIND=wp)             :: cloud_x_ice      ! returns homogeneous freezing rate
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    qc          , & ! cloud specific mass
          &                    qi          , & ! ice specific mass
          &                    dt              ! time step

  !$ACC ROUTINE SEQ
  cloud_x_ice = 0.0_wp
  IF (qc > qmin .AND. t < tfrz_hom) cloud_x_ice  =   qc / dt
  IF (qi > qmin .AND. t > tmelt)    cloud_x_ice  = - qi / dt
      
END FUNCTION cloud_x_ice
!!!=============================================================================================

PURE FUNCTION cloud_to_snow(t,qc,qs,ns,lambda)
  REAL(KIND=wp)             :: cloud_to_snow    ! returns riming snow rate
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    qc          , & ! cloud specific mass
          &                    qs          , & ! snow specific mass
          &                    ns          , & ! snow number
          &                    lambda          ! snow slope parameter (lambda)

  REAL(KIND=wp), PARAMETER  ::      &
      ecs    = 0.9_wp             , & ! Collection efficiency for snow collecting cloud water
      b_rim  =-(v1s+3.0_wp)       , & ! ''
      c_rim  = 2.61_wp*ecs*v0s        ! '' (with pi*gam(v1s+3)/4 = 2.610)

  !$ACC ROUTINE SEQ
  cloud_to_snow = 0.0_wp
  IF (min(qc,qs) > qmin .AND. t > tfrz_hom) THEN
    cloud_to_snow  = (c_rim * ns)  * qc * lambda**b_rim  
  END IF

END FUNCTION cloud_to_snow

!!!=============================================================================================

PURE FUNCTION cloud_to_graupel(t,rho,qc,qg)
  REAL(KIND=wp)             :: cloud_to_graupel  ! returns graupel riming rate
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    rho         , & ! ambient density
          &                    qc          , & ! snow specific mass
          &                    qg              ! graupel specific mass

  REAL(KIND=wp), PARAMETER  :: &
      a_rim  = 4.43_wp,        & ! Constants in riming formula
      b_rim  = 0.94878_wp        ! ''
  
  !$ACC ROUTINE SEQ
  cloud_to_graupel = 0.0_wp
  IF (min(qc,qg) > qmin .AND. t > tfrz_hom) THEN
    cloud_to_graupel = a_rim * qc * (qg*rho)**b_rim
  END IF
      
END FUNCTION cloud_to_graupel

!!!=============================================================================================

PURE FUNCTION rain_to_vapor(t,rho,qc,qr,dvsw,dt)
  REAL(KIND=wp)             :: rain_to_vapor   ! mass from qc to qr 
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    rho         , & ! ambient density
          &                    qc          , & ! specific humidity of cloud
          &                    qr          , & ! specific humidity of rain
          &                    dvsw        , & ! qv - qsat_water(T)
          &                    dt              ! time-step

  REAL(KIND=wp), PARAMETER  ::     &
          b1    =  0.16667_wp    , & ! exponent in power-law relation for mass density
          b2    =  0.55555_wp    , & ! ''
          c1    =  0.61_wp       , & ! ''
          c2    = -0.0163_wp     , & ! ''
          c3    =  1.111e-4_wp   , & ! ''
          a1    =  1.536e-3_wp   , & ! ''
          a2    =  1.0E+0_wp     , & ! constant in rain evap formula
          a3    =  19.0621E+0_wp     ! prefactor (from gamma dist. and properties of air/water)

  REAL(KIND=wp) :: tc, evap_max

  !$ACC ROUTINE SEQ
  rain_to_vapor = 0.0_wp
  IF( qr>qmin .AND. (dvsw+qc <= 0.0_wp)) THEN
    tc            = t - tmelt
    evap_max      = (c1 + tc*(c2 + c3*tc)) * (-dvsw)/dt
    rain_to_vapor = MIN(a1*( a2 + a3 * (qr*rho)**b1 ) * (-dvsw) * (qr*rho)**b2, evap_max)
  ENDIF

END FUNCTION rain_to_vapor

!!!=============================================================================================

PURE FUNCTION rain_to_graupel(t,rho,qc,qr,qi,qs,mi,dvsw,dt)
  REAL(KIND=wp)             :: rain_to_graupel ! freezing rain
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    rho         , & ! ambient density
          &                    qr          , & ! specific humidity of rain
          &                    qc          , & ! cloud liquid specific mass
          &                    qi          , & ! cloud ice specific mass
          &                    qs          , & ! snow specific mass
          &                    mi          , & ! ice crystal mass
          &                    dvsw        , & ! qv-qsat_water(T) 
          &                    dt              ! time-step

  REAL(KIND=wp), PARAMETER  ::     &
      tfrz_rain = tmelt-2.0_wp, &
      a1 = 9.95e-5_wp         , & !FR: 1. coefficient for immersion raindrop freezing: alpha_if
      b1 = 7.0_wp/4.0_wp      , & !FR: 2. coefficient for immersion raindrop freezing: a_if
      c1 = 1.68_wp            , & ! coefficient for raindrop freezing
      c2 = 0.66_wp            , & !FR: 2. coefficient for immersion raindrop freezing: a_if
      c3 = 1.0_wp             , & !FR: 2. coefficient for immersion raindrop freezing: a_if
      c4 = 0.1_wp             , & !FR: 2. coefficient for immersion raindrop freezing: a_if
      a2 = 1.24E-3_wp         , & ! (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)
      b2 = 13.0_wp/8.0_wp     , & ! ''
      qs_crit = 1.e-7_wp

  !$ACC ROUTINE SEQ
  rain_to_graupel  = 0.0_wp
  IF (qr  > qmin .AND. t< tfrz_rain) THEN
    IF( t > tfrz_hom ) THEN
      IF ( dvsw+qc <= 0.0_wp .OR. qr > c4*qc ) THEN
        rain_to_graupel  = (EXP(c2*(tfrz_rain-t))-c3) * (a1*(qr*rho)**b1)
      ENDIF
    ELSE 
      rain_to_graupel  = qr/dt
    ENDIF 
  ENDIF 
  IF (min(qi,qr) > qmin .AND. qs > qs_crit) THEN ! rain + ice creating graupel
    rain_to_graupel = rain_to_graupel +  a2 * (qi/mi) * qr**b2 
  END IF

END FUNCTION rain_to_graupel

!!!=============================================================================================

PURE FUNCTION deposition_auto_conversion(qi,m_ice,ice_dep)

  REAL(KIND=wp)             :: deposition_auto_conversion  
  REAL(KIND=wp), INTENT(IN) :: qi             , & ! ice specific mass
          &                    m_ice          , & ! ice crystal mass
          &                    ice_dep            ! rate of ice deposition (some to snow)

  REAL(KIND=wp), PARAMETER  :: &
      m0_s  = 3.0E-9_wp      , & ! initial mass of snow crystals      
      b     = 2.0_wp/3.0_wp  , & ! 2/3
      xcrit = 1.0_wp             ! threshold parameter

  REAL(KIND=wp) tau_inv

  !$ACC ROUTINE SEQ
  deposition_auto_conversion = 0.0_wp
  IF (qi > qmin) THEN
    tau_inv     = b /((m0_s/m_ice)**b - xcrit)
    deposition_auto_conversion = MAX(0.0_wp,ice_dep)*tau_inv
  END IF
      
END FUNCTION deposition_auto_conversion

!!!=============================================================================================

PURE FUNCTION ice_to_snow(qi,ns,lambda,sticking_eff)
  REAL(KIND=wp)             :: ice_to_snow        ! conversion rate of ice to snow
  REAL(KIND=wp), INTENT(IN) :: qi             , & ! ice specific mass
          &                    ns             , & ! snow number
          &                    lambda         , & ! snow intercept parameter, lambda
          &                    sticking_eff       ! ice sticking effiency

  REAL(KIND=wp), PARAMETER  :: &
      qi0   = 0.0_wp         , & ! critical ice required for autoconversion
      c_iau = 1.0E-3_wp      , & ! coefficient of auto conversion
      c_agg = 2.61_wp*v0s    , & ! coeff of aggregation (2.610 = pi*gam(v1s+3)/4)
      b_agg = -(v1s+3.0_wp)      ! aggregation exponent

  !$ACC ROUTINE SEQ
  ice_to_snow = 0.0_wp
  IF (qi > qmin) THEN
    ice_to_snow = sticking_eff  * (c_iau*MAX(0.0_wp,(qi-qi0)) + qi*(c_agg*ns)*(lambda)**b_agg)
  END IF
      
END FUNCTION ice_to_snow

!!!=============================================================================================

PURE FUNCTION ice_to_graupel(rho,qr,qg,qi,sticking_eff)
  REAL(KIND=wp)             :: ice_to_graupel      ! returns aggregation of ice by graupel
  REAL(KIND=wp), INTENT(IN) :: rho             , & ! density
          &                    qr              , & ! rain specific mass
          &                    qg              , & ! graupel specific mass
          &                    qi              , & ! ice specific mass
          &                    sticking_eff        ! sticking efficiency

  REAL(KIND=wp), PARAMETER  :: &
      a = 1.72_wp            , & ! (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
      b = 7.0_wp/8.0_wp      , & ! ''
      c_agg  = 2.46_wp       , & !
      b_agg  = 0.94878_wp

  !$ACC ROUTINE SEQ
  ice_to_graupel  = 0.0_wp
  IF (qi > qmin) THEN
    IF (qg > qmin) ice_to_graupel = sticking_eff*qi*c_agg*((rho*qg)**b_agg)
    IF (qr > qmin) ice_to_graupel = ice_to_graupel + a*qi*((rho*qr)**b )
  END IF
      
END FUNCTION ice_to_graupel

!!!=============================================================================================

PURE FUNCTION snow_to_rain(t,p,rho,dvsw0,qs)
  REAL(KIND=wp)             :: snow_to_rain    ! melting of snow to form rain
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    p           , & ! ambient  pressure
          &                    rho         , & ! ambient  density
          &                    dvsw0       , & ! qv-qsat_water(T0)
          &                    qs              ! snow specific mass

  REAL(KIND=wp), PARAMETER  :: &
      c1= 79.6863_wp         , & ! Constants in melting formula
      c2= 0.612654E-3_wp     , & ! Constants in melting formula
      a = tx - 389.5_wp      , & ! melting prefactor
      b = 4.0_wp/5.0_wp          ! melting exponent

  !$ACC ROUTINE SEQ
  snow_to_rain = 0.0_wp
  IF (t>MAX(tmelt,tmelt - tx*dvsw0) .AND. qs>qmin) THEN
    snow_to_rain =  (c1/p + c2) *(t - tmelt + a*dvsw0) * (qs*rho)**b
  END IF 

END FUNCTION snow_to_rain

!!!=============================================================================================

PURE FUNCTION snow_to_graupel(t,rho,qc,qs)
  REAL(KIND=wp)             :: snow_to_graupel ! returns conversion rate
  REAL(KIND=wp), INTENT(IN) :: t           , & ! ambient temperature
          &                    rho         , & ! ambient density
          &                    qc          , & ! cloud specific mass
          &                    qs              ! snow specific mass

  REAL(KIND=wp), PARAMETER  :: &
      a = 0.5_wp,              & ! Constants in riming formula
      b = 3.0_wp/4.0_wp          ! ''

  !$ACC ROUTINE SEQ
  snow_to_graupel =  0.0_wp
  IF (min(qc,qs) > qmin .AND. t > tfrz_hom ) THEN
     snow_to_graupel =  a * qc * (qs*rho)**b
  END IF

END FUNCTION snow_to_graupel

!!!=============================================================================================

PURE FUNCTION graupel_to_rain(t,p,rho,dvsw0,qg)
  REAL(KIND=wp)             :: graupel_to_rain ! melting of graupel to form rain
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    p           , & ! ambient pressure
          &                    rho         , & ! ambient density
          &                    dvsw0       , & ! qv-qsat_water(T0)
          &                    qg              ! snow specific mass

  REAL(KIND=wp), PARAMETER  :: &
      c1= 12.31698_wp        , & ! Constants in melting formula
      c2= 7.39441e-05_wp     , & ! Constants in melting formula
      a = tx - 389.5_wp      , & ! melting prefactor
      b = 3.0_wp/5.0_wp          ! melting exponent

  !$ACC ROUTINE SEQ
  graupel_to_rain = 0.0_wp
  IF (t > MAX(tmelt,tmelt - tx*dvsw0) .AND. qg>qmin) THEN
    graupel_to_rain =  (c1/p + c2) * (t - tmelt + a * dvsw0) * (qg*rho)**b
  END IF 

END FUNCTION graupel_to_rain

!!!=============================================================================================

PURE FUNCTION ice_deposition_nucleation(t,qc,qi,ni,dvsi,dt)
  REAL(KIND=wp)             :: ice_deposition_nucleation  ! rate of vapor deposition for new ice
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    qc          , & ! specific humidity of ice
          &                    qi          , & ! specific humidity of ice
          &                    ni          , & ! ice crystal number
          &                    dvsi        , & ! vapor excess with respect to ice sat
          &                    dt              ! time-step

  !$ACC ROUTINE SEQ
  ice_deposition_nucleation =  0.0_wp
  IF (qi <= qmin .AND. (( t<tfrz_het2 .AND. dvsi>0.0_wp) .OR. (t<=tfrz_het1 .AND. qc>qmin))) THEN
    ice_deposition_nucleation =  MIN(m0_ice * ni, MAX(0.0_wp,dvsi)) / dt
  ENDIF

END FUNCTION ice_deposition_nucleation

!!!=============================================================================================

PURE FUNCTION vapor_x_ice(qi,mi,eta,dvsi,dt)
  REAL(KIND=wp)             :: vapor_x_ice     ! returns rate of vapor deposition to ice
  REAL(KIND=wp), INTENT(IN) :: qi          , & ! specific humidity of ice
          &                    mi          , & ! ice crystal mass
          &                    eta         , & ! deposition factor
          &                    dvsi        , & ! vapor excess with respect to ice sat
          &                    dt              ! time-step

  REAL(KIND=wp), PARAMETER  ::                 &
      ami   = 130.0_wp                       , & ! Formfactor for mass-size relation of cld ice
      a     = 4.0_wp * ami**(-1.0_wp/3.0_wp) , & !
      b     = -0.67_wp                           ! exp. for conv. (-1 + 0.33) of ice mass to sfc area 

  !$ACC ROUTINE SEQ
  vapor_x_ice = 0.0_wp
  IF (qi>qmin) THEN
    vapor_x_ice = (a * eta) * qi* (mi**b) * dvsi
    IF (vapor_x_ice > 0._wp) THEN
      vapor_x_ice = MIN(vapor_x_ice, dvsi/dt)
    ELSE
      vapor_x_ice = MAX(vapor_x_ice, dvsi/dt)
      vapor_x_ice = MAX(vapor_x_ice,  -qi/dt)
    END IF
  ENDIF

END FUNCTION vapor_x_ice

!!!=============================================================================================

PURE FUNCTION vapor_x_snow(t,p,rho,qs,ns,lambda,eta,ice_dep,dvsw,dvsi,dvsw0,dt)
  REAL(KIND=wp)             :: vapor_x_snow    ! returns rate of vapor deposition to graupel
  REAL(KIND=wp), INTENT(IN) :: t           , & ! temperature
          &                    p           , & ! ambient pressure
          &                    rho         , & ! ambient density
          &                    qs          , & ! snow specific mass
          &                    ns          , & ! snow number
          &                    lambda      , & ! slope parameter (lambda) snow
          &                    eta         , & ! deposition factor
          &                    ice_dep     , & ! limiter for vapor dep on snow
          &                    dvsw        , & ! qv-qsat_water(T)
          &                    dvsi        , & ! qv-qsat_ice(T)
          &                    dvsw0       , & ! qv-qsat_water(T0)
          &                    dt              ! time-step

  REAL(KIND=wp), PARAMETER  :: &
      nu     = 1.75e-5_wp                , &! kinematic viscosity of air
      a0     = 1.0_wp                    , & !
      a1     = 0.4182_wp * SQRT(v0s/nu)  , & ! 0.26*gam((v1s+5)/2)=0.4182
      a2     = -(v1s+1.0_wp)/2.0_wp      , &
      eps    = 1.e-15_wp                 , & ! ''
      qs_lim = 1.e-7_wp                  , & !
      cnx    = 4.0_wp                    , & !
      b      = 0.8_wp                    , & !
      c1     = 31282.3_wp                , & !
      c2     = 0.241897_wp               , & !
      c3     = 0.28003_wp                , & !
      c4     =-0.146293E-6_wp            

  !$ACC ROUTINE SEQ
  vapor_x_snow = 0.0_wp
  IF (qs>qmin) THEN
    IF ( t < tmelt ) THEN
      vapor_x_snow = (cnx*ns*eta/rho) * (a0 + a1 * lambda**a2) * dvsi  / (lambda*lambda+eps)
      ! 
      ! GZ: This limitation, which was missing in the original graupel scheme,
      ! is crucial for numerical stability in the tropics!
      ! a meaningful distiction between cloud ice and snow
      !
      IF (vapor_x_snow > 0.0_wp) vapor_x_snow = MIN(vapor_x_snow, dvsi/dt-ice_dep)
      IF (qs <= qs_lim) vapor_x_snow = MIN(vapor_x_snow, 0.0_wp)
    ELSE 
      IF ( t > (tmelt-tx*dvsw0) ) THEN
        vapor_x_snow = (c1/p+c2) * MIN(0.0_wp,dvsw0) * (qs*rho)**b
      ELSE
        vapor_x_snow = (c3+c4*p) * dvsw  * (qs*rho)**b
      ENDIF
    ENDIF
    vapor_x_snow = MAX(vapor_x_snow,-qs/dt)
  ENDIF

END FUNCTION vapor_x_snow

!!!=============================================================================================

PURE FUNCTION vapor_x_graupel(t,p,rho,qg,dvsw,dvsi,dvsw0,dt)
  REAL(KIND=wp)             :: vapor_x_graupel   ! graupel-vapor exchange rate
  REAL(KIND=wp), INTENT(IN) :: t             , & ! temperature
          &                    p             , & ! ambient pressure
          &                    rho           , & ! ambient density
          &                    qg            , & ! graupel specific mass
          &                    dvsw          , & ! qv-qsat_water(T)
          &                    dvsi          , & ! qv-qsat_ice(T)
          &                    dvsw0         , & ! qv-qsat_water(T0)
          &                    dt                ! timestep

  REAL(KIND=wp), PARAMETER  :: &
      a1 = 0.398561_wp       , & ! Constants in vapor deposition formula
      a2 =-0.00152398_wp     , & !
      a3 = 2554.99_wp        , & !
      a4 = 2.6531E-7_wp      , & !
      a5 = 0.153907_wp       , & !
      a6 =-7.86703e-07_wp    , & !
      a7 = 0.0418521_wp      , & !
      a8 =-4.7524E-8_wp      , & !
      b  = 0.6_wp            

  !$ACC ROUTINE SEQ
  vapor_x_graupel = 0.0_wp
  IF (qg>qmin) THEN
    IF ( t < tmelt ) THEN
      vapor_x_graupel = (a1 +a2*t+ a3/p + a4*p) * dvsi * (qg*rho)**b
    ELSE 
      IF ( t > (tmelt-tx*dvsw0 ) ) THEN
        vapor_x_graupel = (a5+a6*p) * MIN(0.0_wp,dvsw0) * (qg*rho)**b
      ELSE
        vapor_x_graupel = (a7+a8*p) * dvsw * (qg*rho)**b
      ENDIF
    ENDIF
    vapor_x_graupel = MAX(vapor_x_graupel,-qg/dt)
  ENDIF

END FUNCTION vapor_x_graupel

END MODULE mo_aes_graupel
