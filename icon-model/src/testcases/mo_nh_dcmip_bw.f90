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

! Initialization routines for the (moist) baroclinic instability as used in DCMIP2016
!
! Initialization routines for the Ullrich, Melvin, Staniforth and Jablonowski
! baroclinic instability, including moisture and terminator toy chemistry.
! This baroclinic instability test is part of the DCMIP2016 test case suite.
!
! Literature
! Dynamical Core Model Intercomparison Project (DCMIP2016) Test Case Document

MODULE mo_nh_dcmip_bw

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: omega=>earth_angular_velocity !, &
!    &                               a => earth_radius,             &
!    &                               rd,                            &
!    &                               g => grav,                     &
!    &                               cp=> cpd,                      &
!    &                               Mvap => vtmpc1
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp,                ONLY: cells2edges_scalar
  USE mo_grid_config,         ONLY: grid_rescale_factor
  USE mo_nh_testcases_nml,    ONLY: dcmip_bw
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: finish
  USE mo_run_config,          ONLY: iqv, ntracer
  USE mo_sync,                ONLY: sync_patch_array, SYNC_E

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_dcmip_bw'

  !=======================================================================
  !    Physical constants
  !=======================================================================

  REAL(wp), PARAMETER ::              &
       a     = 6371220.0_wp,          & ! Reference Earth's Radius (m)
       Rd    = 287.0_wp,              & ! Ideal gas const dry air (J kg^-1 K^1)
       g     = 9.80616_wp,            & ! Gravity (m s^2)
       cp    = 1004.5_wp,             & ! Specific heat capacity (J kg^-1 K^1)
       Mvap  = 0.608_wp,              & ! Ratio of molar mass of dry air/water
       p0    = 100000.0_wp              ! surface pressure (Pa)

  !=======================================================================
  !    Test case parameters
  !=======================================================================
  REAL(wp), PARAMETER ::              &
       T0E        = 310._wp    ,      & ! temperature at equatorial surface (K)
       T0P        = 240._wp    ,      & ! temperature at polar surface (K)
       B          = 2._wp      ,      & ! jet half-width parameter
       K          = 3._wp      ,      & ! jet width parameter
       lapse      = 0.005_wp            ! lapse rate parameter

  REAL(wp), PARAMETER ::              &
       pertu0     = 0.5_wp     ,      & ! SF Perturbation wind velocity (m/s)
       pertr      = 1._wp/6._wp,      & ! SF Perturbation radius (Earth radii)
       pertup     = 1.0_wp     ,      & ! Exp. perturbation wind velocity (m/s)
       pertexpr   = 0.1_wp     ,      & ! Exp. perturbation radius (Earth radii)
       pertlon    = pi/9._wp   ,      & ! Perturbation longitude
       pertlat    = 2._wp*pi/9._wp,   & ! Perturbation latitude
       pertz      = 15000._wp  ,      & ! Perturbation height cap
       dxepsilon  = 1.e-5_wp            ! Small value for numerical derivatives

  REAL(wp), PARAMETER ::              &
       moistqlat  = 2._wp*pi/9._wp,   & ! Humidity latitudinal width
       moistqp    = 34000._wp,        & ! Humidity vertical pressure width
       moisttr    = 0.1_wp,           & ! Vertical cut-off pressure for humidity
       moistqs    = 1.e-12_wp,        & ! Humidity above cut-off
       moistq0    = 0.018_wp,         & ! Maximum specific humidity
       moistqr    = 0.9_wp,           & ! Maximum saturation ratio
       moisteps   = 0.622_wp,         & ! Ratio of gas constants
       moistT0    = 273.16_wp,        & ! Reference temperature (K)
       moistE0Ast = 610.78_wp           ! Saturation vapor pressure at T0 (Pa) 


  PUBLIC :: init_nh_dcmip_bw


CONTAINS


  !>
  !! Setup idealized conditions for DCMIP2016 baroclinic instability.
  !!
  !! Setup idealized conditions for DCMIP2016 baroclinic instability 
  !! test case.
  !!
  SUBROUTINE init_nh_dcmip_bw (p_patch, p_nh_prog, p_nh_diag, p_int, p_metrics)

    TYPE(t_patch),        INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_int_state),    INTENT(IN)    :: &  !< interpolation state vector
      &  p_int

    TYPE(t_nh_metrics),   INTENT(INOUT) :: &  !< NH metrics state
      &  p_metrics

    ! local
    CHARACTER(LEN = *), PARAMETER :: routine = modname//':init_nh_dcmip_bw'
    INTEGER :: jc, je, jk, jb              ! loop indices for cell, edge, level, block
    INTEGER :: i_rlstart, i_rlend
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER :: nlev                        ! number of vertical (full) levels

    REAL(wp) :: zlon, zlat                 ! lat/lon of cell circumcenter
    REAL(wp) :: zu, zv                     ! zonal and meridional velocity component
    REAL(wp) :: zphis(nproma,p_patch%nlev) ! surface geopotential (m^2 s^-2)

    REAL(wp) :: zqv(nproma,p_patch%nlev)                      ! specific humidity (kg/kg)
    REAL(wp) :: ztempv_e(nproma,p_patch%nlev,p_patch%nblks_e) ! virtual temperature at cell edge midpoints
    REAL(wp) :: z_me(nproma,p_patch%nlev,p_patch%nblks_e)     ! geometric height at cell edge points
    !
!!$    INTEGER, PARAMETER :: deep=0       ! deep atmosphere (1 = yes or 0 = no)
!!$    INTEGER, PARAMETER :: moist=1      ! include moisture (1 = yes or 0 = no)
!!$    INTEGER, PARAMETER :: pertt=0      ! type of perturbation (0 = exponential, 1 = stream function)
    INTEGER, PARAMETER :: zcoords=1    ! 1 if z is specified, 0 if p is specified


    ! Sanity check
    IF (dcmip_bw%moist == 1) THEN
      ! make sure that sufficient tracer fields are allocated
      IF (.NOT. ASSOCIATED(p_nh_prog%tracer)) THEN
        CALL finish (routine, 'Tracer field not allocated')
      ENDIF
      IF (SIZE(p_nh_prog%tracer,4) < 1) THEN
        CALL finish (routine, 'Testcase necessitates allocation of at least 1 tracer field')
      ENDIF
    ENDIF


    nlev = p_patch%nlev

    i_rlstart = 1
    i_rlend   = min_rlcell

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    !
    ! Init prognostic variables rho and exner (theta_v)
    ! The initial velocities provided for the cell center are ignored 
    ! at this point. They will be recomputed for the cell edges lateron.
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zlon,zlat,zu,zv,zphis,zqv)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev

        DO jc = i_startidx, i_endidx

          ! get geographical coordinates of cell circumcenter
          !
          zlon = p_patch%cells%center(jc,jb)%lon
          zlat = p_patch%cells%center(jc,jb)%lat
          CALL baroclinic_wave_test(deep    = dcmip_bw%deep               ,& !in
            &                       moist   = dcmip_bw%moist              ,& !in
            &                       pertt   = dcmip_bw%pertt              ,& !in
            &                       X       = 1._wp/grid_rescale_factor   ,& !in
            &                       lon     = zlon                        ,& !in
            &                       lat     = zlat                        ,& !in
            &                       p       = p_nh_diag%pres(jc,jk,jb)    ,& !out
            &                       z       = p_metrics%z_mc(jc,jk,jb)    ,& !in
            &                       zcoords = zcoords                     ,& !in
            &                       u       = zu                          ,& !out (not used)
            &                       v       = zv                          ,& !out (not used)
            &                       t       = p_nh_diag%temp (jc,jk,jb)   ,& !out
            &                       thetav  = p_nh_prog%theta_v(jc,jk,jb) ,& !out
            &                       phis    = zphis(jc,jk)                ,& !out (not used)
            &                       ps      = p_nh_diag%pres_sfc(jc,jb)   ,& !out
            &                       rho     = p_nh_prog%rho(jc,jk,jb)     ,& !out
            &                       q       = zqv(jc,jk)                   & !out
            &                                                              )

          ! compute virtual temperature from temperature and qv
          p_nh_diag%tempv(jc,jk,jb) = p_nh_diag%temp(jc,jk,jb) * (1._wp + Mvap*zqv(jc,jk))
          !
          ! compute exner pressure
          p_nh_prog%exner(jc,jk,jb) = (p_nh_diag%pres(jc,jk,jb) / p0)**(Rd / cp)
        ENDDO  ! jc
      ENDDO  ! jk
      IF ( (dcmip_bw%moist==1) .AND. (ntracer>=1)) THEN
        ! initialize qv-Tracer
        p_nh_prog%tracer(i_startidx:i_endidx,1:nlev,jb,iqv) = zqv(i_startidx:i_endidx,1:nlev) 
      ENDIF
    ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL


    ! Interpolate virtual temperature and geometric height to edge midpoints 
    ! (required for velocity initialization) 
    CALL cells2edges_scalar(p_nh_diag%tempv, p_patch, p_int%c_lin_e, ztempv_e)
    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, p_int%c_lin_e, z_me)
    CALL sync_patch_array(SYNC_E,p_patch,ztempv_e)
    CALL sync_patch_array(SYNC_E,p_patch,z_me)


    !
    ! Init prognostic normal velocity field
    !
    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,zlon,zlat,zu,zv)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev

        DO je = i_startidx, i_endidx
          ! get geographical coordinates of edge midpoint
          !
          zlon = p_patch%edges%center(je,jb)%lon
          zlat = p_patch%edges%center(je,jb)%lat
          CALL init_velocity_field(deep = dcmip_bw%deep,             & !in
            &                      pertt= dcmip_bw%pertt,            & !in
            &                      X    = 1._wp/grid_rescale_factor, & !in
            &                      lon  = zlon,                      & !in
            &                      lat  = zlat,                      & !in
            &                      z    = z_me(je,jk,jb),            & !in
            &                      tempv= ztempv_e(je,jk,jb),        & !in
            &                      u    = zu,                        & !out
            &                      v    = zv                         ) !out

          p_nh_prog%vn(je,jk,jb) = zu * p_patch%edges%primal_normal(je,jb)%v1   &
            &                    + zv * p_patch%edges%primal_normal(je,jb)%v2
        ENDDO  ! je

      ENDDO  ! jk
    ENDDO  !jb
!$OMP ENDDO
!$OMP END PARALLEL

  END SUBROUTINE init_nh_dcmip_bw



  ! The following routines were provided by the DCMIP2016 organizers.
  !=======================================================================
  !
  !  Date:  July 29, 2015
  !
  !  Functions for setting up idealized initial conditions for the
  !  Ullrich, Melvin, Staniforth and Jablonowski baroclinic instability.
  !
  !  SUBROUTINE baroclinic_wave_sample(
  !    deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q)
  !
  !  Options:
  !     deep    deep atmosphere (1 = yes or 0 = no)
  !    moist    include moisture (1 = yes or 0 = no)
  !    pertt    type of perturbation (0 = exponential, 1 = stream function)
  !        X    Earth scaling factor
  !
  !  Given a point specified by: 
  !      lon    longitude (radians) 
  !      lat    latitude (radians) 
  !      p/z    pressure (Pa) / height (m)
  !  zcoords    1 if z is specified, 0 if p is specified
  !
  !  the functions will return:
  !        p    pressure if z is specified and zcoords = 1 (Pa)
  !        u    zonal wind (m s^-1)
  !        v    meridional wind (m s^-1)
  !        t    temperature (K)
  !   thetav    virtual potential temperature (K)
  !     phis    surface geopotential (m^2 s^-2)
  !       ps    surface pressure (Pa)
  !      rho    density (kj m^-3)
  !        q    water vapor mixing ratio (kg/kg)
  !
  !
  !  Author: Paul Ullrich
  !          University of California, Davis
  !          Email: paullrich@ucdavis.edu
  !
  !=======================================================================

  !=======================================================================
  !    Generate the baroclinic instability initial conditions
  !=======================================================================
  SUBROUTINE baroclinic_wave_test(deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho,q)

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !     input/output params parameters at given location
    !-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: &
         deep,       & ! Deep (1) or Shallow (0) test case
         moist,      & ! Moist (1) or Dry (0) test case
         pertt         ! Perturbation type

    REAL(wp), INTENT(IN)  :: &
         lon,        & ! Longitude (radians)
         lat,        & ! Latitude (radians)
         X             ! Earth scaling parameter

    REAL(wp), INTENT(INOUT) :: &
         p,            & ! Pressure (Pa)
         z               ! Altitude (m)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                       ! 0 if p coordinates are specified

    REAL(wp), INTENT(OUT) :: &
         u,          & ! Zonal wind (m s^-1)
         v,          & ! Meridional wind (m s^-1)
         t,          & ! Temperature (K)
         thetav,     & ! Virtual potential temperature (K)
         phis,       & ! Surface Geopotential (m^2 s^-2)
         ps,         & ! Surface Pressure (Pa)
         rho,        & ! density (kg m^-3)
         q             ! water vapor mixing ratio (kg/kg)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    REAL(wp) :: eta

    !------------------------------------------------
    !   Pressure and temperature
    !------------------------------------------------
    if (zcoords .eq. 1) then
       CALL evaluate_pressure_temperature(deep, X, lat, z, p, t)
    else
       CALL evaluate_z_temperature(deep, X, lat, p, z, t)
    end if

    !-----------------------------------------------------
    !   Initialize surface pressure
    !-----------------------------------------------------
    ps = p0

    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    CALL init_velocity_field(deep, pertt, X, lon, lat, z, t, u, v)


    !-----------------------------------------------------
    !   Initialize surface geopotential
    !-----------------------------------------------------
    phis = 0.d0

    !-----------------------------------------------------
    !   Initialize density
    !-----------------------------------------------------
    rho = p / (Rd * t)

    !-----------------------------------------------------
    !   Initialize specific humidity
    !-----------------------------------------------------
    if (moist .eq. 1) then
       eta = p/p0

       if (eta .gt. moisttr) then
          q = moistq0 * exp(- (lat/moistqlat)**4)          &
               * exp(- ((eta-1.d0)*p0/moistqp)**2)
       else
          q = moistqs
       end if

       ! Convert virtual temperature to temperature
       t = t / (1.d0 + Mvap * q)

    else
       q = 0.d0
    end if

    !-----------------------------------------------------
    !   Initialize virtual potential temperature
    !-----------------------------------------------------
    thetav = t * (1.d0 + 0.61d0 * q) * (p0 / p)**(Rd / cp)
!DR In order to be consistent with the above computation of t, 
!DR 0.61 should be replaced by Mvap 
!DR    thetav = t * (1.d0 + Mvap * q) * (p0 / p)**(Rd / cp)

  END SUBROUTINE baroclinic_wave_test


  !-----------------------------------------------------------------------
  !  Init velocity field
  !-----------------------------------------------------------------------
  SUBROUTINE init_velocity_field(deep, pertt, X, lon, lat, z, tempv, u, v)
    INTEGER , INTENT(IN)  :: deep      ! Deep (1) or Shallow (0) test case
    INTEGER , INTENT(IN)  :: pertt     ! Perturbation type
    REAL(wp), INTENT(IN)  :: X         ! Earth scaling parameter
    REAL(wp), INTENT(IN)  :: lon       ! Longitude (radians)
    REAL(wp), INTENT(IN)  :: lat       ! Latitude (radians)
    REAL(wp), INTENT(IN)  :: z         ! Altitude (m)
    REAL(wp), INTENT(IN)  :: tempv     ! virtual temperature (K)
    REAL(wp), INTENT(OUT) :: u         ! zonal velocity
    REAL(wp), INTENT(OUT) :: v         ! meridional velocity

    ! local
    REAL(wp) :: aref, omegaref
    REAL(wp) :: inttermU
    REAL(wp) :: bigU
    REAL(wp) :: rcoslat
    REAL(wp) :: omegarcoslat
    REAL(wp) :: rratio
    REAL(wp) :: T0, constH, constC, scaledZ, inttau2

    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5d0 * (T0E + T0P)
    constH = Rd * T0 / g
    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)
    scaledZ = z / (B * constH)
    inttau2 = constC * z * exp(- scaledZ**2)

    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
       rratio = 1.d0
    else
       rratio = (z + aref) / aref;
    end if

    inttermU = (rratio * cos(lat))**(K - 1.d0) - (rratio * cos(lat))**(K + 1.d0)
    bigU = g / aref * K * inttau2 * inttermU * tempv  !t
    if (deep .eq. 0) then
       rcoslat = aref * cos(lat)
    else
       rcoslat = (z + aref) * cos(lat)
    end if

    omegarcoslat = omegaref * rcoslat

    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.d0

    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------

    ! Exponential type
    if (pertt .eq. 0) then
       u = u + evaluate_exponential(lon, lat, z)

       ! Stream function type
    elseif (pertt .eq. 1) then
       u = u - 1.d0 / (2.d0 * dxepsilon) *                       &
            ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
            - evaluate_streamfunction(lon, lat - dxepsilon, z))

       v = v + 1.d0 / (2.d0 * dxepsilon * cos(lat)) *            &
            ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
            - evaluate_streamfunction(lon - dxepsilon, lat, z))
    end if

  END SUBROUTINE init_velocity_field

  !-----------------------------------------------------------------------
  !    Calculate pointwise pressure and temperature
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_pressure_temperature(deep, X, lat, z, p, t)

    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(wp), INTENT(IN)  :: &
         X,          & ! Earth scaling ratio
         lat,        & ! Latitude (radians)
         z             ! Altitude (m)

    REAL(wp), INTENT(OUT) :: &
         p,          & ! Pressure (Pa)
         t             ! Temperature (K)

    REAL(wp) :: aref
    REAL(wp) :: T0, constA, constB, constC, constH, scaledZ
    REAL(wp) :: tau1, tau2, inttau1, inttau2
    REAL(wp) :: rratio, inttermT

    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = a / X

    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rd * T0 / g

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)

    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
         + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)

    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
       rratio = 1.d0
    else
       rratio = (z + aref) / aref;
    end if

    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**K &
         - K / (K + 2.d0) * (rratio * cos(lat))**(K + 2.d0)

    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    t = 1.d0 / (rratio**2 * (tau1 - tau2 * inttermT))

    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    p = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))

  END SUBROUTINE evaluate_pressure_temperature

  !-----------------------------------------------------------------------
  !    Calculate pointwise z and temperature given pressure
  !-----------------------------------------------------------------------
  SUBROUTINE evaluate_z_temperature(deep, X, lat, p, z, t)

    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(wp), INTENT(IN)  :: &
         X,          & ! Earth scaling ratio
         lat,        & ! Latitude (radians)
         p             ! Pressure (Pa)

    REAL(wp), INTENT(OUT) :: &
         z,          & ! Altitude (m)
         t             ! Temperature (K)

    INTEGER :: ix

    REAL(wp) :: z0, z1, z2
    REAL(wp) :: p0, p1, p2

    z0 = 0.d0
    z1 = 10000.d0

    CALL evaluate_pressure_temperature(deep, X, lat, z0, p0, t)
    CALL evaluate_pressure_temperature(deep, X, lat, z1, p1, t)

    DO ix = 1, 100
       z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)

       CALL evaluate_pressure_temperature(deep, X, lat, z2, p2, t)

       IF (ABS((p2 - p)/p) .lt. 1.0d-13) THEN
          EXIT
       END IF

       z0 = z1
       p0 = p1

       z1 = z2
       p1 = p2
    END DO

    z = z2

    CALL evaluate_pressure_temperature(deep, X, lat, z, p0, t)

  END SUBROUTINE evaluate_z_temperature

  !-----------------------------------------------------------------------
  !    Exponential perturbation function
  !-----------------------------------------------------------------------
  REAL(wp) FUNCTION evaluate_exponential(lon, lat, z)

    REAL(wp), INTENT(IN)  :: &
         lon,        & ! Longitude (radians)
         lat,        & ! Latitude (radians)
         z             ! Altitude (meters)

    REAL(wp) :: greatcircler, perttaper

    ! Great circle distance
    greatcircler = 1.d0 / pertexpr &
         * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
       perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
       perttaper = 0.d0
    end if

    ! Zonal velocity perturbation
    if (greatcircler < 1.d0) then
       evaluate_exponential = pertup * perttaper * exp(- greatcircler**2)
    else
       evaluate_exponential = 0.d0
    end if

  END FUNCTION evaluate_exponential

  !-----------------------------------------------------------------------
  !    Stream function perturbation function
  !-----------------------------------------------------------------------
  REAL(wp) FUNCTION evaluate_streamfunction(lon, lat, z)

    REAL(wp), INTENT(IN)  :: &
         lon,        & ! Longitude (radians)
         lat,        & ! Latitude (radians)
         z             ! Altitude (meters)

    REAL(wp) :: greatcircler, perttaper, cospert

    ! Great circle distance
    greatcircler = 1.d0 / pertr &
         * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
       perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
       perttaper = 0.d0
    end if

    ! Horizontal tapering of stream function
    if (greatcircler .lt. 1.d0) then
       cospert = cos(0.5d0 * pi * greatcircler)
    else
       cospert = 0.d0
    end if

    evaluate_streamfunction = &
         (- pertu0 * pertr * perttaper * cospert**4)

  END FUNCTION evaluate_streamfunction


END MODULE mo_nh_dcmip_bw

