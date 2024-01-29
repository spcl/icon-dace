! Function for setting up idealized tropical cyclone initial conditions
!
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
!
!
! Given a point specified by:
!
!          longitude (radians)
!         latitude (radians)
!         height
!
! the functions will return:
!        u        zonal wind (m s^-1)
!        v        meridional wind (m s^-1)
!        w        vertical velocity (m s^-1)
!        t        temperature (K)
!        phis        surface geopotential (m^2 s^-2)
!        ps        surface pressure (Pa)
!        rho        density (kj m^-3)
!        q        specific humidity (kg/kg)
!        qi        tracers (kg/kg)
!      p       pressure (Pa)
!
! Initial data are currently identical to:
!
!                 Reed, K. A., and C. Jablonowski, 2011: An analytic
!                 vortex initialization technique for idealized tropical
!                 cyclone studies in AGCMs. Mon. Wea. Rev., 139, 689-710.
!
! Author: Kevin A. Reed (University of Michigan, kareed@umich.edu)
!          version 1
!          5/21/2012

MODULE mo_nh_dcmip_tc

  ! items of the infrastructure
  ! ---------------------------
  !
  ! precision
  USE mo_kind,                ONLY: wp
  !
  ! type definitions
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,      ONLY: t_int_state
  !
  ! dimensions and indices
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv
  !
  ! interpolation
  USE mo_intp,                ONLY: cells2edges_scalar


  ! items for the computations
  ! --------------------------
  !
  USE mo_grid_config,         ONLY: grid_rescale_factor ,& ! (m/m)     Small planet scaling factor (=1/x)
    &                               grid_sphere_radius  ,& ! (m)       Earth radius/x
    &                               grid_angular_velocity  ! (1/s)     Earth rotation*x
  ! 
  USE mo_math_constants,      ONLY: deg2rad                ! (rad/deg) Convert degree to radian
  !
  USE mo_physical_constants,  ONLY: rd                  ,& ! (J/K/kg)  Ideal gas const dry air (J kg^-1 K^1)
    &                               rd_o_cpd            ,& ! ()        rd/cpd
    &                               p0ref               ,& ! (Pa)      Reference pressure for Exner function
    &                               grav                   ! (m s^2)   Gravity (m s^2)

  USE mo_nh_diagnose_pres_temp,ONLY:diagnose_pres_temp

  USE mo_sync,                ONLY: SYNC_E, sync_patch_array

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: init_nh_dcmip_tc

  ! parameters defining the initial atmospheric state of
  ! the background and the tropical cyclone
  !-----------------------------------------------------

  REAL(wp), PARAMETER :: rp1        = 282000.0_wp       ,& ! (m)   Radius for calculation of PS
    &                    dp         = 1115.0_wp         ,& ! (Pa)  Delta P for calculation of PS
    &                    zp         = 7000.0_wp         ,& ! (m)   Height for calculation of P
    &                    q0         = 0.0210_wp         ,& !       q at surface from Jordan
    &                    gamma      = 0.0070_wp         ,& ! (K/m) Lapse rate
    &                    ts0        = 302.15_wp         ,& ! (K)   Surface temperature (SST)
    &                    p00        = 101500.0_wp       ,& ! (Pa)  Global mean surface pressure
    &                    p0         = 100000.0_wp       ,& ! (Pa)  p for model level calculation
    &                    cen_lat    =  10.0_wp*deg2rad  ,& ! (rad) Center latitude of initial vortex
    &                    cen_lon    = 180.0_wp*deg2rad  ,& ! (rad) Center longitufe of initial vortex
    &                    zq1        = 3000.0_wp         ,& ! (m)   Height 1 for q calculation
    &                    zq2        = 8000.0_wp         ,& ! (m)   Height 2 for q calculation
    &                    exppr      = 1.5_wp            ,& !       Exponent for r dependence of p
    &                    exppz      = 2.0_wp            ,& !       Exponent for z dependence of p
    &                    ztrop      = 15000.0_wp        ,& ! (m)   Tropopause Height
    &                    qtrop      = 1.e-11_wp         ,& !       Tropopause specific humidity
    &                    consttv    = 0.6080_wp         ,& !       Constant for Virtual Temp Conversion
    &                    epsilon    = 1.e-25_wp            !       Small number to aviod dividing by zero in wind calc

CONTAINS

  SUBROUTINE init_nh_dcmip_tc( p_patch, p_nh_prog, p_nh_diag, p_metrics, p_int )

    TYPE(t_patch),      INTENT(inout), TARGET :: p_patch   !< patch on which computation is performed
    TYPE(t_nh_prog),    INTENT(inout)         :: p_nh_prog !< prognostic state vector
    TYPE(t_nh_diag),    INTENT(inout)         :: p_nh_diag !< diagnostic state vector
    TYPE(t_nh_metrics), INTENT(in)            :: p_metrics !< NH metrics state
    TYPE(t_int_state),  INTENT(in)            :: p_int     !< interpolation state

    ! dimensions, loops and indices
    !
    ! grids
    INTEGER  :: i_rlstart, i_rlend   !< grid levels
    INTEGER  :: i_nchdom             !< children domains per grid level
    !
    ! rows per block
    INTEGER  :: nblks_c, nblks_e     !< number of rows of cells and edges
    INTEGER  :: i_startblk, i_endblk !< loop index range
    INTEGER  :: jb                   !< loop index
    !
    ! loop in row
    INTEGER  :: i_startidx, i_endidx !< loop index range
    INTEGER  :: jc, je               !< cell and edge index
    !
    ! levels
    INTEGER  :: nlev                 !< number of full levels
    INTEGER  :: jk                   !< index of level


    ! variables
    !
    ! coordinates
    ! -----------
    REAL(wp), POINTER     :: lon(:,:)          !< (rad)   longitude
    REAL(wp), POINTER     :: lat(:,:)          !< (rad)   latitude
    REAL(wp), POINTER     :: z(:,:,:)          !< (m)     height
    !
    ! prognostic state variables
    ! --------------------------
    REAL(wp), POINTER     :: exner(:,:,:)      !< ()      Exner pressure
    REAL(wp), POINTER     :: theta_v(:,:,:)    !< (K)     virtual potential temperature
    REAL(wp), POINTER     :: rho(:,:,:)        !< (kg/m3) density
    REAL(wp), POINTER     :: x(:,:,:,:)        !< (kg/kg) tracer
    REAL(wp), POINTER     :: q(:,:,:)          !< (kg/kg) specific humidity
    !
    REAL(wp), POINTER     :: vn(:,:,:)         !< (m/s)   normal wind
    REAL(wp), POINTER     :: w(:,:,:)          !< (m/s)   vertical velocity
    !
    ! diagnostic state variables
    ! --------------------------
    REAL(wp), POINTER     :: p(:,:,:)          !< (Pa)    pressure
    REAL(wp), POINTER     :: ps(:,:)           !< (Pa)    pressure at surface
    REAL(wp), POINTER     :: tv(:,:,:)         !< (K)     virtual temperature


    ! auxiliary variables
    ! -------------------
    REAL(wp)              :: rp                !< (m)     scaled radius for ps, rp1/x
    REAL(wp)              :: f                 !< (1/s)   scaled Coriolis parameter at cen_lat, f*x
    !
    REAL(wp), ALLOCATABLE :: gr(:,:)           !< (m)     great circle distance to center of vortex
    !
    REAL(wp), ALLOCATABLE :: d1(:,:)           !< ()      for computing d
    REAL(wp), ALLOCATABLE :: d2(:,:)           !< ()      for computing d
    REAL(wp), ALLOCATABLE :: d(:,:)            !< ()      for computing ufac and vfac
    REAL(wp), ALLOCATABLE :: ufac(:,:)         !< ()      for computing u anv v from tangential wind of cyclone
    REAL(wp), ALLOCATABLE :: vfac(:,:)         !< ()      for computing u anv v from tangential wind of cyclone
    !
    REAL(wp), ALLOCATABLE :: vtc_c(:,:,:)      !< (m/s)   tangential wind of cyclone at cell centers
    REAL(wp), ALLOCATABLE :: vtc_e(:,:,:)      !< (m/s)   tangential wind of cyclone at edge centers

    REAL(wp)              :: exponent          !<         exponent
    REAL(wp)              :: t0                !< (K)     virtual temperature at surface
    REAL(wp)              :: ttrop             !< (K)     virtual temperature at tropopause
    REAL(wp)              :: ptrop             !< (Pa)    pressure at tropopause

    ! auxiliary constants
    rp       = rp1 * grid_rescale_factor
    f        = 2._wp*grid_angular_velocity*SIN(cen_lat)
    exponent = rd*gamma/grav
    t0       = ts0*(1.0_wp+consttv*q0)
    ttrop    = t0 - gamma*ztrop
    ptrop    = p00*(ttrop/t0)**(1.0_wp/exponent)

    ! number of rows in a block
    nblks_c  = p_patch%nblks_c ! for cells
    nblks_e  = p_patch%nblks_e ! for edges

    ! number of full levels
    nlev = p_patch%nlev

    ! 1. initialize variables defined at cell centers
    ! ===============================================

    ! coordinates
    ! -----------
    lon        => p_patch%cells%center%lon
    lat        => p_patch%cells%center%lat
    z          => p_metrics%z_mc
    !
    ! prognostic state variables
    ! --------------------------
    exner      => p_nh_prog%exner
    theta_v    => p_nh_prog%theta_v
    rho        => p_nh_prog%rho
    x          => p_nh_prog%tracer(:,:,:,:)
    q          => p_nh_prog%tracer(:,:,:,iqv)
    !
    w          => p_nh_prog%w
    !
    ! diagnostic state variables
    ! --------------------------
    p          => p_nh_diag%pres
    ps         => p_nh_diag%pres_sfc
    tv         => p_nh_diag%tempv
    !
    ! auxiliary variables
    ! -------------------
    ALLOCATE( gr   (nproma,nblks_c) )
    !
    ALLOCATE( d1   (nproma,nblks_c) )
    ALLOCATE( d2   (nproma,nblks_c) )
    ALLOCATE( d    (nproma,nblks_c) )
    ALLOCATE( ufac (nproma,nblks_c) )
    ALLOCATE( vfac (nproma,nblks_c) )
    !
    ALLOCATE( vtc_c(nproma,nlev,nblks_c) )
    !
    gr(:,:)      = 0.0_wp
    !
    d1  (:,:)    = 0.0_wp
    d2  (:,:)    = 0.0_wp
    d   (:,:)    = 0.0_wp
    ufac(:,:)    = 0.0_wp
    vfac(:,:)    = 0.0_wp
    !
    vtc_c(:,:,:) = 0.0_wp

    ! loop index ranges
    ! -----------------
    i_rlstart  = 1
    i_rlend    = min_rlcell
    !
    i_nchdom   = MAX(1,p_patch%n_childdom)
    !
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx, i_endidx)
    rows_c: DO jb = i_startblk,  i_endblk

      CALL get_indices_c(p_patch   , jb      , &
        &                i_startblk, i_endblk, &
        &                i_startidx, i_endidx, &
        &                i_rlstart , i_rlend  )

      cells_2d: DO jc = i_startidx, i_endidx

        ! great circle distance (gr) to vortex center                    (eq.97)
        !-----------------------------------------------------------------------
        gr(jc,jb) = grid_sphere_radius                                           &
          &        *ACOS( SIN(cen_lat)*SIN(lat(jc,jb))                           &
          &              +COS(cen_lat)*COS(lat(jc,jb))*COS(lon(jc,jb)-cen_lon) )

        ! pressure at surface                                           (eqs.99)
        !-----------------------------------------------------------------------
        ps(jc,jb) = p00-dp*EXP(-(gr(jc,jb)/rp)**exppr)

        ! factors for u and v components of the tangential wind    (eqs.109-113)
        !-----------------------------------------------------------------------
        d1(jc,jb) =  SIN(cen_lat)*COS(lat(jc,jb))                         &
          &         -COS(cen_lat)*SIN(lat(jc,jb))*COS(lon(jc,jb)-cen_lon)
        d2(jc,jb) =  COS(cen_lat)*SIN(lon(jc,jb)-cen_lon)
        d (jc,jb) =  MAX(epsilon, SQRT(d1(jc,jb)**2+d2(jc,jb)**2))

        ufac(jc,jb) = d1(jc,jb)/d(jc,jb)
        vfac(jc,jb) = d2(jc,jb)/d(jc,jb)

      END DO cells_2d

      ! initialize atmosphere ...
      !-----------------------------------------------------------------------

      ! tracers
      x(:,:,jb,:) = 0.0_wp

      levels_c: DO jk = 1, nlev

        cells_3d: DO jc = i_startidx, i_endidx

          IF (z(jc,jk,jb) > ztrop) THEN ! ... above tropopause ...

            ! pressure                                                   (eqs.94+98)
            !-----------------------------------------------------------------------
            p(jc,jk,jb) = ptrop*EXP(-(grav*(z(jc,jk,jb)-ztrop))/(rd*ttrop))

            ! virtual temperature                                       (eqs.92,102)
            !-----------------------------------------------------------------------
            tv(jc,jk,jb) = ttrop

            ! specific humidity                                              (eq.91)
            !-----------------------------------------------------------------------
            q(jc,jk,jb)   = qtrop

            ! tangential wind of cyclone                                    (eq.108)
            !-----------------------------------------------------------------------
            vtc_c(jc,jk,jb) = 0.0_wp

          ELSE ! ... below tropopause ...

            ! pressure                                                   (eqs.94+98)
            !-----------------------------------------------------------------------
            p(jc,jk,jb) = (p00-dp*EXP(-(gr(jc,jb)/rp)**exppr-(z(jc,jk,jb)/zp)**exppz)) &
              &          *((t0-gamma*z(jc,jk,jb))/t0)**(1.0_wp/exponent)

            ! virtual temperature                                       (eqs.92,102)
            !-----------------------------------------------------------------------
            tv(jc,jk,jb) =  (t0-gamma*z(jc,jk,jb))                                                     &
              &            /(1.0_wp + (exppz*rd*(t0-gamma*z(jc,jk,jb))*z(jc,jk,jb))                    &
              &                      /(grav*zp**exppz*(1.0_wp-p00/dp*EXP( (gr(jc,jb)/rp)**exppr        &
              &                                                          +(z(jc,jk,jb)/zp)**exppz) )))

            ! specific humidity                                              (eq.91)
            !-----------------------------------------------------------------------
            q(jc,jk,jb)   = q0*EXP(-(z(jc,jk,jb)/zq1)-(z(jc,jk,jb)/zq2)**exppz)

            ! tangential wind of cyclone                                    (eq.108)
            !-----------------------------------------------------------------------
            vtc_c(jc,jk,jb) = -f*gr(jc,jb)/2.0_wp                                                           &
              &               +SQRT( (f*gr(jc,jb)/2.0_wp)**2                                                &
              &                     - ( exppr*(gr(jc,jb)/rp)**exppr*rd*(t0-gamma*z(jc,jk,jb)))              &
              &                      /( 1.0_wp+exppz*rd*(t0-gamma*z(jc,jk,jb))*z(jc,jk,jb)/(grav*zp**exppz) &
              &                        -p00/dp*EXP((gr(jc,jb)/rp)**exppr+(z(jc,jk,jb)/zp)**exppz)))

          END IF

          ! ... and at all levels

          ! Exner pressure
          !-----------------------------------------------------------------------
          exner(jc,jk,jb) = (p(jc,jk,jb)/p0ref)**rd_o_cpd

          ! virtual potential temparature
          !-----------------------------------------------------------------------
          theta_v(jc,jk,jb) = tv(jc,jk,jb)/exner(jc,jk,jb)

          ! density of moist air                                          (eq.106)
          !-----------------------------------------------------------------------
          rho(jc,jk,jb) = p(jc,jk,jb)/(rd*tv(jc,jk,jb) )

          ! vertical velocity
          !-----------------------------------------------------------------------
          w(jc,jk,jb) = 0.0_wp

        END DO cells_3d
      END DO levels_c
    END DO rows_c
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL diagnose_pres_temp ( p_metrics, p_nh_prog,   &
      &                       p_nh_prog, p_nh_diag,   &
      &                       p_patch               )


    DEALLOCATE(gr)
    DEALLOCATE(d1, d2, d)
    DEALLOCATE(ufac, vfac)

    ! 2. interpolate from cell centers to edge centers
    ! ================================================

    ALLOCATE(vtc_e(nproma,nlev,nblks_e))

    CALL cells2edges_scalar(vtc_c,                 &
      &                     p_patch,p_int%c_lin_e, &
      &                     vtc_e                  )

    CALL sync_patch_array(SYNC_E,p_patch,vtc_e)

    DEALLOCATE(vtc_c)

    ! 3. initialize variables defined at edge centers
    ! ===============================================

    ! coordinates
    ! -----------
    lon        => p_patch%edges%center%lon
    lat        => p_patch%edges%center%lat
    !
    ! prognostic state variables
    ! --------------------------
    vn         => p_nh_prog%vn

    ! auxiliary variables
    ! -------------------
    ALLOCATE( d1   (nproma,nblks_e) )
    ALLOCATE( d2   (nproma,nblks_e) )
    ALLOCATE( d    (nproma,nblks_e) )
    ALLOCATE( ufac (nproma,nblks_e) )
    ALLOCATE( vfac (nproma,nblks_e) )
    !
    d1(:,:)      = 0.0_wp
    d2(:,:)      = 0.0_wp
    d(:,:)       = 0.0_wp
    !
    ufac(:,:)    = 0.0_wp
    vfac(:,:)    = 0.0_wp

    ! loop index ranges
    ! -----------------
    i_rlstart = 2
    i_rlend   = min_rledge
    !
    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)
    !
    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,i_startidx, i_endidx)
    rows_e: DO jb = i_startblk,  i_endblk

      CALL get_indices_e(p_patch   , jb      , &
        &                i_startblk, i_endblk, &
        &                i_startidx, i_endidx, &
        &                i_rlstart , i_rlend  )

      edges_2d: DO je = i_startidx, i_endidx

        ! factors for u and v components of the tangential wind    (eqs.109-111)
        !-----------------------------------------------------------------------

        d1(je,jb) =  SIN(cen_lat)*COS(lat(je,jb))                             &
          &         -COS(cen_lat)*SIN(lat(je,jb))*COS(lon(je,jb)-cen_lon)
        d2(je,jb) =  COS(cen_lat)*SIN(lon(je,jb)-cen_lon)
        d (je,jb) =  MAX(epsilon, SQRT(d1(je,jb)**2+d2(je,jb)**2))

        ufac(je,jb) = d1(je,jb)/d(je,jb)
        vfac(je,jb) = d2(je,jb)/d(je,jb)

      END DO edges_2d

      levels_e: DO jk = 1, nlev

        edges_3d: DO je = i_startidx, i_endidx

          ! compute wind normal to edges from
          ! tangential wind of cyclone at edges                  (eqs.108,112,113)
          !-----------------------------------------------------------------------
          vn(je,jk,jb) = vtc_e(je,jk,jb)*ufac(je,jb) * p_patch%edges%primal_normal(je,jb)%v1   &
            &           +vtc_e(je,jk,jb)*vfac(je,jb) * p_patch%edges%primal_normal(je,jb)%v2

        END DO edges_3d
      END DO levels_e
    END DO rows_e
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DEALLOCATE(d1, d2, d)
    DEALLOCATE(ufac, vfac)
    DEALLOCATE(vtc_e)

  END SUBROUTINE init_nh_dcmip_tc

END MODULE mo_nh_dcmip_tc

