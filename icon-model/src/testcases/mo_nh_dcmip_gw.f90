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

! Subroutine to initialize testcase 31 (gravity wave on a small planet) proposed
! for the DCMIP summer school
!
! The non-hydrostatic gravity wave test examines the response of models to
! short time-scale wavemotion triggered by a localized perturbation. The
! formulation presented in this document is new, but is based on previous
! approaches by Skamarock et al. (JAS 1994), Tomita and Satoh (FDR 2004), and
! Jablonowski et al. (NCAR Tech Report 2008)
!
! Literature
! - Dynamical Core Model Intercomparison Project (DCMIP)
!   Test Case Document (P. Ullrich et al, 2012)
! - Baldauf, M. et al. (2013): in preparation

MODULE mo_nh_dcmip_gw

   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd, grav, p0ref, cpd, cvd_o_rd
   USE mo_hydro_adjust,         ONLY: hydro_adjust, hydro_adjust_iterative
   USE mo_math_constants,       ONLY: pi, deg2rad
   USE mo_impl_constants,       ONLY: min_rlcell, min_rledge, min_rlvert, MAX_CHAR_LENGTH
   USE mo_parallel_config,      ONLY: nproma
   USE mo_loopindices,          ONLY: get_indices_c, get_indices_e, get_indices_v
   USE mo_model_domain,         ONLY: t_patch
   USE mo_intp_data_strc,       ONLY: t_int_state
   USE mo_io_config,            ONLY: inextra_3d
   USE mo_grid_config,          ONLY: grid_sphere_radius, grid_angular_velocity
   USE mo_dynamics_config,      ONLY: lcoriolis
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_math_utilities,       ONLY: rotate_latlon
   USE mo_math_divrot,          ONLY: div
   USE mo_exception,            ONLY: finish


   IMPLICIT NONE


   PRIVATE

   PUBLIC :: init_nh_dcmip_gw
   PUBLIC :: init_nh_gw_analyt
   PUBLIC :: gw_clat           ! namelist variable
   PUBLIC :: gw_u0             ! namelist variable
   PUBLIC :: gw_delta_temp     ! namelist variable

   REAL(wp) :: gw_clat                    ! Lat of perturbation center [deg]
   REAL(wp) :: gw_u0                      ! maximum amplitude
                                          ! of the zonal wind          [m s^-1]
   REAL(wp) :: gw_delta_temp              ! Max amplitude of perturbation [K]


!--------------------------------------------------------------------

CONTAINS
!-------------------------------------------------------------------------
!

!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state vector for the DCMIP nh gravity 
  !! wave test 
  !!
  !! The non-hydrostatic gravity wave test examines the response of models to 
  !! short time-scale wavemotion triggered by a localized perturbation. The 
  !! formulation presented in this document is new, but is based on previous 
  !! approaches by Skamarock et al. (JAS 1994), Tomita and Satoh (FDR 2004), and
  !! Jablonowski et al. (NCAR Tech Report 2008) 
  !!
  SUBROUTINE init_nh_dcmip_gw( p_patch, p_nh_prog, p_nh_diag, p_metrics, l_hydro_adjust)

    TYPE(t_patch),        INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics

    LOGICAL,              INTENT(IN)    :: &  !< TRUE: compute hydrostatically balanced state 
      & l_hydro_adjust                        !  by integrating the discretized third equation 
                                              !  of motion.

    REAL(wp) :: z_lon, z_lat          !< geographical coordinates

    REAL(wp) ::           &           !< surface temperature              [K]
      & z_tsfc(nproma,p_patch%nblks_c) 

    REAL(wp) ::           &           !< zonal and meridional velocity
      & zu(nproma), zv(nproma)


    REAL(wp) :: dist                  !< great circle distance            [m]
    REAL(wp) :: shape_func            !< value of horiz. shape function
    REAL(wp) :: theta_pert            !< theta_perturbation               [K]
    REAL(wp) :: zsin, zcos
    INTEGER  :: jc, je, jk, jb        !< loop indices
    INTEGER  :: ist                   !< error status
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend  
    INTEGER  :: nlev, nlevp1          !< number of full and half levels

    REAL(wp) :: brunt2                !< Brunt-Vaisala frequency squared  [s^-1]
    REAL(wp) :: big_g                 !< auxiliary constant
    REAL(wp) :: phic                  !< Lat of perturbation center       [rad]

    REAL(wp), ALLOCATABLE :: zrelhum(:,:,:) !< dummy relative humidity field
    REAL(wp), ALLOCATABLE :: zqv    (:,:,:) !< dummy specific moisture field 
                                            !  both needed in case of hydrostatic adjustment. 
    ! test case parameters
    !
    REAL(wp), PARAMETER :: peq = 100000._wp !< reference surface pressure
                                            !< at the equator             [Pa]
    REAL(wp), PARAMETER :: teq = 300._wp    !< surface temperature
                                            !< at the equator             [K]
    REAL(wp), PARAMETER :: brunt= 0.01_wp   !< Brunt-Vaisala frequency    [s^-1]

    REAL(wp), PARAMETER :: lambdac = 2.0_wp*pi/3.0_wp !< Lon of perturbation center

    REAL(wp), PARAMETER :: width= 5000._wp  !< width for perturbation     [m]

    REAL(wp), PARAMETER :: delta_theta = 1.0_wp  !< Max amplitude of perturbation [K]

    REAL(wp), PARAMETER :: lz = 20000._wp   !< vert. wavelength of perturbation [m]

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_dcmip_gw:init_nh_dcmip_gw' 
!--------------------------------------------------------------------
!

    ! center of temperature/density perturbation in radians
    phic = deg2rad * gw_clat


    ! initialize some constants:
    !
    brunt2  = brunt * brunt    ! Brunt-Vaisala frequency squared 
    big_g   = (grav*grav)/(brunt2*cpd)


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!    I: Init background state     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !
    ! Init prognostic variables vn, w
    !
!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,z_lat,zu,zv)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jk = 1, nlev
        DO je = i_startidx, i_endidx


          ! get geographical coordinates of edge midpoint
          !
          z_lat = p_patch%edges%center(je,jb)%lat

          ! init velocity field
          !
          ! zonal velocity
          zu(je) = gw_u0 * COS(z_lat)

          ! meridional velocity
          zv(je) = 0._wp

          ! compute normal wind component
          p_nh_prog%vn(je,jk,jb) = zu(je) * p_patch%edges%primal_normal(je,jb)%v1  &
            &                    + zv(je) * p_patch%edges%primal_normal(je,jb)%v2

        ENDDO  ! je
      ENDDO  ! jk
    ENDDO ! jb
!$OMP ENDDO NOWAIT


    i_rlstart = 1
    i_rlend   = min_rlcell

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    !
    ! Init prognostic variables exner, theta_v and rho
    !
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,z_lat)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jc = i_startidx, i_endidx

        ! get geographical coordinates of cell circumcenter
        !
        z_lat = p_patch%cells%center(jc,jb)%lat


        ! init surface temperature
        !
        z_tsfc(jc,jb) = big_g + (teq-big_g) * exp( -(gw_u0*brunt2/(4.0_wp*grav*grav)) &
          &           *(gw_u0+2.0_wp*grid_angular_velocity*grid_sphere_radius)  &
          &           * (COS(2.0_wp*z_lat)-1.0_wp) )


        ! init surface pressure field
        !
        p_nh_diag%pres_sfc(jc,jb) = peq*exp( (gw_u0/(4.0_wp*big_g*rd))                 &
          &                       * (gw_u0+2.0_wp*grid_angular_velocity*grid_sphere_radius)&
          &                       * (COS(2.0_wp*z_lat)-1.0_wp))&
          &                       * (z_tsfc(jc,jb)/teq)**(cpd/rd)

      ENDDO  ! jc



      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx


          ! init temperature
          !
          p_nh_diag%temp(jc,jk,jb) = big_g*(1.0_wp - exp(brunt2*p_metrics%z_mc(jc,jk,jb)/grav)) &
            &            + z_tsfc(jc,jb)*exp(brunt2*p_metrics%z_mc(jc,jk,jb)/grav)


          ! init pressure field
          !
          p_nh_diag%pres(jc,jk,jb) = p_nh_diag%pres_sfc(jc,jb)                  &
            &                      * ( (big_g/z_tsfc(jc,jb))                    &
            &                      * exp(-brunt2*p_metrics%z_mc(jc,jk,jb)/grav) &
            &                      + 1.0_wp - (big_g/z_tsfc(jc,jb)) )**(cpd/Rd)


          ! init exner pressure
          !
          p_nh_prog%exner(jc,jk,jb) = (p_nh_diag%pres(jc,jk,jb)/p0ref)**(rd/cpd)


          ! init virtual potential temperature
          !
          p_nh_prog%theta_v(jc,jk,jb) = p_nh_diag%temp(jc,jk,jb)/p_nh_prog%exner(jc,jk,jb)


          ! store theta_v background state for visualization purposes
          !
          IF (inextra_3d > 0) THEN
            p_nh_diag%extra_3d(jc,jk,jb,1) = p_nh_prog%theta_v(jc,jk,jb)
          ENDIF


          ! init density of moist air
          !
          p_nh_prog%rho(jc,jk,jb) = p_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref  &
            &                     /rd/p_nh_prog%theta_v(jc,jk,jb)

        ENDDO !jc
      ENDDO !jk


      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx

          ! vertical velocity
          !
          p_nh_prog%w(jc,jk,jb) = 0._wp

        ENDDO
      ENDDO

    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL


    IF (l_hydro_adjust) THEN
      !
      ALLOCATE(zrelhum(SIZE(p_nh_diag%temp,1),SIZE(p_nh_diag%temp,2),SIZE(p_nh_diag%temp,3)), &
        &      zqv    (SIZE(p_nh_diag%temp,1),SIZE(p_nh_diag%temp,2),SIZE(p_nh_diag%temp,3)), &
        &      STAT=ist)
      IF (ist /= 0) CALL finish(TRIM(routine),'Allocation of zrelhum and zqv failed')
      !
      ! as this testcase is dry, we set relhum to zero
      zrelhum(:,:,:) = 0._wp
      zqv    (:,:,:) = 0._wp

      !
      ! remark: calling hydro_adjust would be equally fine here
      ! 
      CALL hydro_adjust_iterative ( p_patch       = p_patch,                      & !in
        &                           p_nh_metrics  = p_metrics,                    & !in
        &                           temp_ini      = p_nh_diag%temp,               & !in
        &                           rh_ini        = zrelhum,                      & !in
        &                           exner         = p_nh_prog%exner,              & !inout
        &                           theta_v       = p_nh_prog%theta_v,            & !inout
        &                           rho           = p_nh_prog%rho,                & !inout
        &                           qv            = zqv,                          & !inout
        &                           luse_exner_fg = .TRUE.,                       & !in
        &                           opt_exner_lbc = p_nh_prog%exner(:,nlev,:)     ) !in

      ! cleanup
      DEALLOCATE(zrelhum, zqv, STAT=ist)
      IF (ist /= 0) CALL finish(TRIM(routine),'Deallocation of zrelhum and zqv failed')
    END IF   



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!    II. Add potential temperature perturbation  !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    i_rlstart = 1
    i_rlend   = min_rlcell

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,z_lat,z_lon,zsin,zcos,dist, &
!$OMP            shape_func,theta_pert)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! get geographical coordinates of cell circumcenter
          !
          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lon = p_patch%cells%center(jc,jb)%lon

          zsin = SIN(z_lat) * SIN(phic)
          zcos = COS(z_lat) * COS(phic)

          ! great circle distance with 'a/X' 
          dist  = grid_sphere_radius * ACOS (zsin + zcos*COS(z_lon-lambdac)) 
          shape_func = (width**2)/(width**2 + dist**2)

          theta_pert = delta_theta*shape_func*SIN(2.0_wp*pi*p_metrics%z_mc(jc,jk,jb)/lz)

          ! add perturbation to virtual potential temperature field
          !
          p_nh_prog%theta_v(jc,jk,jb) = p_nh_prog%theta_v(jc,jk,jb) + theta_pert


          ! add perturbation to temperature field 
          ! (not strictly necessary, only for plotting purposes) 
          p_nh_diag%temp(jc,jk,jb) = p_nh_diag%temp(jc,jk,jb) + theta_pert    &
            &                      * (p_nh_diag%pres(jc,jk,jb)/p0ref)**(rd/cpd)


! In order to match the equation of state, a perturbation should be added to the density 
! field as well. However, the DCMIP tets case uses an initially unperturbed density field.
!          ! add perturbation to density field
!          p_nh_prog%rho(jc,jk,jb) = p_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref  &
!            &                     /rd/p_nh_prog%theta_v(jc,jk,jb)
! end Test


        ENDDO !jc
      ENDDO !jk
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE init_nh_dcmip_gw





  !>
  !! Initialization of prognostic state vector for the nh gravity 
  !! wave test. 
  !!
  !! The non-hydrostatic gravity wave test examines the response of models to 
  !! short time-scale wavemotion triggered by a localized perturbation. The 
  !! formulation presented here is based on work by Baldauf et al. (2012). 
  !! For this particular setup an analytical reference solution is available.  
  !!
  !! Available options:
  !! - gravity wave without coriolis force and without background (soild body) flow
  !! - gravity wave with coriolis force (f-plane) approximation and without 
  !!   background (solid body) flow
  !! - gravity wave with coriolis force and with background (soild body) flow.
  !!   This third option is somewhat special in the sense, that the sum of the 
  !!   coriolis force and centrifugal force exactly cancel the additional 
  !!   centrifugal force due to the background flow. I.e. in an inertial 
  !!   frame, the atmosphere is at rest. However, in the rotating frame, 
  !!   a wind speed of gw_u0*cos(\phi) is observed.  
  !!
  !! Literature
  !! - Baldauf, M. et al. (2013): in preparation
  !!
  SUBROUTINE init_nh_gw_analyt( p_patch, p_nh_prog, p_nh_diag, p_metrics, p_int)

    TYPE(t_patch),        INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics

    TYPE(t_int_state),    INTENT(IN)    :: &  !< interpolation state
      &  p_int

    REAL(wp) :: z_lat, z_lon          !< geographical latitude and longitude

    REAL(wp) :: delta
    REAL(wp) :: rhos                  !< surface density                  [kg/m**3]
    REAL(wp) :: shape_func            !< shape function
    REAL(wp) :: temp_pert             !< temperature perturbation         [K]
    REAL(wp) :: rho_pert              !< density perturbation             [kg/m**3]
    REAL(wp) :: temp_b, rho_b
    REAL(wp) :: zpsi(nproma,p_patch%nblks_v) ! horizontal stream function
    INTEGER  :: jc, je, jv, jk, jb    !< loop indices
    INTEGER  :: ilv1, ibv1, ilv2, ibv2 !< vertex line and block indices
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom  
    INTEGER  :: nlev, nlevp1          !< number of full and half levels
    INTEGER  :: iorient
    REAL(wp) :: phic                  !< Lat of perturbation center       [rad]
    REAL(wp) :: zomega                !< earth angular velocity in the case gw_u0/=0

    ! test case parameters
    !
    REAL(wp), PARAMETER :: t0    = 250._wp    !< surface temperature        [K]

    REAL(wp), PARAMETER :: kappa = 100._wp    !< perturbation width

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_dcmip_gw:init_nh_gw_analyt'

    ! Note:
    ! p0ref = 100000.0_wp   !> [Pa]  (mo_physical_constants.f90)   
!--------------------------------------------------------------------
!
    ! center of temperature/density perturbation in radians
    phic = deg2rad * gw_clat

    ! initialize some constants:
    !
    delta = grav/(rd * t0)

    ! surface density
    rhos  = p0ref/(rd * t0) 


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!    I: Init background state     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !
    ! Init prognostic variables vn, w
    !
!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)


    i_rlstart = 1
    i_rlend   = min_rlvert

    i_startblk = p_patch%verts%start_blk(i_rlstart,1)
    i_endblk   = p_patch%verts%end_blk(i_rlend,i_nchdom)


    !
    ! Init prognostic variables vn, w
    !
    ! Use stream function initialization, in order to get a discretely 
    ! non-divergent horizontal wind field
    ! 
    ! compute velocity stream function at vertices
    !
!$OMP DO PRIVATE(jv,jb,i_startidx,i_endidx,z_lat)
    DO jb = i_startblk, i_endblk

      CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


        DO jv = i_startidx, i_endidx

          ! get geographical coordinates of vertices
          !
          z_lat = p_patch%verts%vertex(jv,jb)%lat

          ! get streamfunction
          zpsi(jv,jb) = -grid_sphere_radius * gw_u0 * SIN(z_lat)

        ENDDO  ! jv
    ENDDO ! jb
!$OMP ENDDO



    ! compute normal velocities at edge center
    !
    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,ilv1,ibv1,ilv2,ibv2,iorient)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          ! get edge vertices
          ilv1 = p_patch%edges%vertex_idx(je,jb,1)
          ibv1 = p_patch%edges%vertex_blk(je,jb,1)
          ilv2 = p_patch%edges%vertex_idx(je,jb,2)
          ibv2 = p_patch%edges%vertex_blk(je,jb,2)

          iorient = p_patch%edges%tangent_orientation(je,jb)

          ! compute normal wind component
          p_nh_prog%vn(je,jk,jb) = iorient                                &
            &                  * (zpsi(ilv2,ibv2) - zpsi(ilv1,ibv1))      &
            &                  / p_patch%edges%primal_edge_length(je,jb)

        ENDDO  ! je
      ENDDO  ! jk
    ENDDO  ! jb
!$OMP ENDDO





    i_rlstart = 1
    i_rlend   = min_rlcell

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    !
    ! Init prognostic variables exner, theta_v and rho
    !
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jc = i_startidx, i_endidx

        ! init surface pressure field (for visualization purposes only)
        !
        p_nh_diag%pres_sfc(jc,jb) = p0ref

      ENDDO  ! jc



      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx


          ! init temperature (background)
          !
          p_nh_diag%temp(jc,jk,jb) = t0


          ! init pressure field (background)
          !
          p_nh_diag%pres(jc,jk,jb) = p0ref * exp(-delta * p_metrics%z_mc(jc,jk,jb))


          ! init exner pressure (background)
          !
          p_nh_prog%exner(jc,jk,jb) = (p_nh_diag%pres(jc,jk,jb)/p0ref)**(rd/cpd)


          ! init virtual potential temperature (background)
          !
          p_nh_prog%theta_v(jc,jk,jb) = p_nh_diag%temp(jc,jk,jb)/p_nh_prog%exner(jc,jk,jb)


          ! store theta_v background state for visualization purposes
          !
          IF (inextra_3d > 0) THEN
            p_nh_diag%extra_3d(jc,jk,jb,1) = p_nh_prog%theta_v(jc,jk,jb)
          ENDIF


          ! init density of moist air (background)
          !
          p_nh_prog%rho(jc,jk,jb) = p_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref  &
            &                     /rd/p_nh_prog%theta_v(jc,jk,jb)

        ENDDO !jc
      ENDDO !jk


      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx

          ! vertical velocity
          !
          p_nh_prog%w(jc,jk,jb) = 0._wp

        ENDDO
      ENDDO

    ENDDO !jb
!$OMP END DO


   CALL hydro_adjust ( p_patch, p_metrics, p_nh_prog%rho,  &
                     & p_nh_prog%exner, p_nh_prog%theta_v  )



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!    II. Add temperature and density perturbation  !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,z_lat,z_lon,shape_func,temp_b, &
!$OMP            temp_pert,rho_b,rho_pert)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! get geographical latitude of cell circumcenter
          !
          z_lat  = p_patch%cells%center(jc,jb)%lat
          z_lon  = p_patch%cells%center(jc,jb)%lon

          CALL rotate_latlon(z_lat, z_lon, phic, 0._wp)

          ! note that from now on, z_lat and z_lon are given with respect to 
          ! the rotated north pole (at (lat,lon)=(phic,0.0))

          shape_func = exp(kappa*(SIN(z_lat)-1._wp))                             &
            &        * SIN(pi * p_metrics%z_mc(jc,jk,jb)/p_metrics%z_ifc(jc,1,jb))


          temp_b   = gw_delta_temp * shape_func
          temp_pert= temp_b * exp(0.5_wp * delta * p_metrics%z_mc(jc,jk,jb))


          rho_b    = - rhos * temp_b/t0
          rho_pert = rho_b * exp(-0.5_wp * delta * p_metrics%z_mc(jc,jk,jb))


          ! add perturbation to temperature field 
          ! (not strictly necessary, only done for plotting) 
          p_nh_diag%temp(jc,jk,jb) = p_nh_diag%temp(jc,jk,jb) + temp_pert


          ! add perturbation to virtual potential temperature field
          !
          p_nh_prog%theta_v(jc,jk,jb) = p_nh_prog%theta_v(jc,jk,jb)                       &
            &                         + temp_pert * (p0ref/p_nh_diag%pres(jc,jk,jb))**(rd/cpd)



          ! In order to match the equation of state, a perturbation should be 
          ! added to the density field as well. 
          p_nh_prog%rho(jc,jk,jb) = p_nh_prog%rho(jc,jk,jb) + rho_pert


        ENDDO !jc
      ENDDO !jk
    ENDDO !jb
!$OMP END DO



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!    III. Initialize Coriolis parameter, if Coriolis effect is considered !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! f-plane approximation, with f_0 = 2 \Omega SIN(\phi_0)
    ! \phi_0 = 45 \deg
    !
    IF ( lcoriolis ) THEN

      ! test case version without background flow
      IF ( gw_u0==0._wp ) THEN 

        ! center of f-plane
        z_lat = 0.25_wp * pi


        i_rlstart = 1
        i_rlend   = min_rlcell

        i_startblk = p_patch%cells%start_blk(i_rlstart,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jc = i_startidx, i_endidx
! Note that a modified scaling factor (10 instead of 50) is used for the Coriolis parameter
            p_patch%cells%f_c(jc,jb) = 0.2_wp * 2._wp * grid_angular_velocity * SIN(z_lat)
          ENDDO  ! jc
        ENDDO  ! jb
!$OMP ENDDO NOWAIT


        i_rlstart = 1
        i_rlend   = min_rledge

        i_startblk = p_patch%edges%start_blk(i_rlstart,1)
        i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx)
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, i_rlstart, i_rlend)

          DO je = i_startidx, i_endidx
! Note that a modified scaling factor (10 instead of 50) is used for the Coriolis parameter
            p_patch%edges%f_e(je,jb) = 0.2_wp * 2._wp * grid_angular_velocity * SIN(z_lat)
          ENDDO  ! je
        ENDDO  ! jb
!$OMP ENDDO NOWAIT


        i_rlstart = 1
        i_rlend   = min_rlvert

        i_startblk = p_patch%verts%start_blk(i_rlstart,1)
        i_endblk   = p_patch%verts%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jv,i_startidx,i_endidx)
        DO jb = i_startblk, i_endblk

          CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jv = i_startidx, i_endidx
! Note that a modified scaling factor (10 instead of 50) is used for the Coriolis parameter
            p_patch%verts%f_v(jv,jb) = 0.2_wp * 2._wp * grid_angular_velocity * SIN(z_lat)
          ENDDO  ! jv
        ENDDO  ! jb
!$OMP ENDDO NOWAIT

      ELSE ! gw_u0 .NE. 0._wp


        ! earth angular velocity is chosen such, that the metric term (u**/r*tan(\phi)) 
        ! is balanced by the sum of coriolis and centrifugal force.
        ! \OMEGA = -gw_u0/r
        !
        ! no f-plane approximation !!
        !
        zomega = -gw_u0/grid_sphere_radius ! note that the scaled radius must be used


        i_rlstart = 1
        i_rlend   = min_rlcell

        i_startblk = p_patch%cells%start_blk(i_rlstart,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,z_lat)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jc = i_startidx, i_endidx
            z_lat = p_patch%cells%center(jc,jb)%lat
            p_patch%cells%f_c(jc,jb) = 2._wp * zomega * SIN(z_lat)
          ENDDO  ! jc
        ENDDO  ! jb
!$OMP ENDDO NOWAIT


        i_rlstart = 1
        i_rlend   = min_rledge

        i_startblk = p_patch%edges%start_blk(i_rlstart,1)
        i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,z_lat)
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, i_rlstart, i_rlend)


          DO je = i_startidx, i_endidx
            z_lat = p_patch%edges%center(je,jb)%lat
            p_patch%edges%f_e(je,jb) = 2._wp * zomega * SIN(z_lat)
          ENDDO  ! je
        ENDDO  ! jb
!$OMP ENDDO NOWAIT


        i_rlstart = 1
        i_rlend   = min_rlvert

        i_startblk = p_patch%verts%start_blk(i_rlstart,1)
        i_endblk   = p_patch%verts%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jv,i_startidx,i_endidx,z_lat)
        DO jb = i_startblk, i_endblk

          CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jv = i_startidx, i_endidx
            z_lat = p_patch%verts%vertex(jv,jb)%lat
            p_patch%verts%f_v(jv,jb) = 2._wp * zomega * SIN(z_lat)
          ENDDO  ! jv
        ENDDO  ! jb
!$OMP ENDDO NOWAIT


        !
        ! compute centrifugal force at edge (normal component): stored on ddt_vn_phy in order to be
        ! added as an external term to the right-hand-side of the horizontal momentum equation
        !
        i_rlstart = 1
        i_rlend   = min_rledge

        i_startblk = p_patch%edges%start_blk(i_rlstart,1)
        i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,z_lat)
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              z_lat = p_patch%edges%center(je,jb)%lat
              p_nh_diag%ddt_vn_phy(je,jk,jb) = -(zomega**2)*grid_sphere_radius*COS(z_lat)*SIN(z_lat) &
                &            * p_patch%edges%primal_normal(je,jb)%v2
            ENDDO
          ENDDO
        ENDDO

      ENDIF  ! gw_u0

    ENDIF  ! lcoriolis
!$OMP END PARALLEL


    ! cross check: analytic solutions are only available for the following setups_
    ! 1) gw_u0  = 0 and lcorio = .FALSE.
    ! 2) gw_u0  = 0 and lcorio = .TRUE.
    ! 3) gw_u0 \= 0 and lcorio = .TRUE.
    !
    IF ( (gw_u0 /= 0._wp) .AND. .NOT. lcoriolis ) THEN
      CALL finish(TRIM(routine),'lcoriolis=True required, if gw_u0/=0')
    ENDIF


    ! diag for Output
    CALL div(p_nh_prog%vn, p_patch, p_int, p_nh_diag%div)


  END SUBROUTINE init_nh_gw_analyt

!--------------------------------------------------------------------
END MODULE mo_nh_dcmip_gw
