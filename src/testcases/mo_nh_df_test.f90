! Module for class of deformational flow tests
!
! Test case for horizontal advection problems.
!
! Literature:
! - Nair and Lauritzen (2010): A Class of Deformational Flow Test-Cases
!   for the Advection Problems on the Sphere. Submitted to JCP
! Modification by Daniel Reinert, DWD (2011-03-11)
! - adapted to the NH-core
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

MODULE mo_nh_df_test

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH
  USE mo_physical_constants,  ONLY: rd, cpd, p0ref
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_math_constants,      ONLY: pi
  USE mo_math_types,          ONLY: t_geographical_coordinates, t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: gnomonic_proj, gc2cc!, az_eqdist_proj
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ntracer
  USE mo_grid_config,         ONLY: grid_sphere_radius

  IMPLICIT NONE

  PRIVATE

  !DEFINED PARAMETERS:
  REAL(wp), PARAMETER :: zp0     = 100000._wp !< surface pressure
  REAL(wp), PARAMETER :: zt0     = 300._wp    !< temperature (isothermal)

  REAL(wp), PARAMETER :: deg2rad = pi/180._wp  !< convert degrees to radians
  REAL(wp), PARAMETER :: tottime = 1036800._wp !< total time 12 days

                                         !< case 4
  REAL(wp), PARAMETER ::   &             !< flow field amplitude [m/s]
    &   u1 = 73.741_wp, u2 = 61.451_wp
!    &   u1 = 2.4_wp, u2 = 2.0_wp, u3 = 1._wp !< original values
                                              !< for unit sphere
  ! position of north pole Nair et al. (2010) (no rotation)
  REAL(wp), PARAMETER :: npole_lon = 0._wp
  REAL(wp), PARAMETER :: npole_lat = 0.5_wp * pi

  ! peak and background tracer values as in Nair et al. (2010)
  REAL(wp), PARAMETER :: c_bell      = 0.9_wp,  & !< peak value bells
    &                    c_slotted   = 1.0_wp,  & !< peak value slotted cylinder
    &                    c_gauss     = 0.95_wp, & !< peak value gaussian bells
    &                    c_lin_bell  = -2._wp,  & !< peak value lin. correl. bells
    &                    c_nlin_bell = -0.8_wp, & !< peak value nlin bells
    &                    b_bell      = 0.1_wp,  & !< background value bells
    &                    b_slotted   = 0.1_wp,  & !< background value slotted cylinder
    &                    b_gauss     = 5._wp,   & !< exp. decay gaussian bells
    &                    b_lin_bell  = 2.1_wp,  & !< background lin. correlated bells
    &                    b_nlin_bell = 0.9_wp     !< background nlin. correlated bells

  REAL(wp), ALLOCATABLE ::   &      !< distance vector cell center --> barycenter
    &  df_distv_barycenter(:,:,:,:) !< of departure region
  INTEGER, ALLOCATABLE ::   &       !< line and block indices of cell centers in
    &  df_cell_indices(:,:,:,:)     !< which the calculated barycenters are located


  PUBLIC :: init_nh_state_prog_dftest
  PUBLIC :: get_nh_df_velocity
  PUBLIC :: df_distv_barycenter
  PUBLIC :: df_cell_indices


CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE init_nh_state_prog_dftest
  !!
  !! Short description:
  !! Initialization of prognostic state vector for deformational flow
  !! test case.
  !!
  SUBROUTINE init_nh_state_prog_dftest( ptr_patch, ptr_nh_prog, ptr_nh_diag,  &
    &                                   ptr_int, ptr_ext_data,                &
    &                                   p_rotate_axis_deg, p_ctest_name,      &
    &                                   linit_tracer_fv, tracer_inidist_list )

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,INTENT(IN) :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_int_state),  INTENT(IN)  :: &
      &  ptr_int
    TYPE(t_nh_prog), INTENT(INOUT)  :: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag), INTENT(INOUT)  :: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_external_data), INTENT(INOUT) :: & !< external data
      &  ptr_ext_data

    CHARACTER(len=MAX_CHAR_LENGTH), INTENT(IN) :: p_ctest_name
    REAL(wp), INTENT(IN)  :: p_rotate_axis_deg !< Earths rotation axis pitch
                                               !< angle in deg
    LOGICAL, INTENT(IN)  :: linit_tracer_fv    !< fv init. for tracer fields
    INTEGER, INTENT(IN)  :: tracer_inidist_list(:)  !< selected initial tracer distributions

    REAL(wp) :: rovcp

    INTEGER  :: jk,jb,jt   !< loop indices
    INTEGER  :: nblks_c,npromz_c
    INTEGER  :: nlev, nlevp1                   !< number of full and half levels
    INTEGER  :: nlen
  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    rovcp  = rd/cpd   ! kappa

    ! topography
    ptr_ext_data%atm%topography_c(:,:) = 0.0_wp

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zp0

    ! init contravariant correction (no topography)
    ptr_nh_diag%w_concorr_c(:,:,:) = 0.0_wp


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      ! bottom level vertical velocity
      ptr_nh_prog%w(1:nlen,nlevp1,jb) = 0.0_wp

      ! bottom level contravariant correction (no topography)
      ptr_nh_diag%w_concorr_c(1:nlen,nlevp1,jb) = 0.0_wp

      ! bottom interface density field rho_ic
      ptr_nh_diag%rho_ic(1:nlen,nlevp1,jb) = 1._wp

      DO jk = 1, nlev

        ! init vertical velocity
        ptr_nh_prog%w(1:nlen,jk,jb) = 0.0_wp

        ! init contravariant correction (no topography)
        ptr_nh_diag%w_concorr_c(1:nlen,jk,jb) = 0.0_wp

        ! init full level pressure
        ptr_nh_diag%pres(1:nlen,jk,jb) = zp0 

        ! init density field rho
        ptr_nh_prog%rho(1:nlen,jk,jb) = 1._wp

        ! init density field rho_ic (half level)
        ptr_nh_diag%rho_ic(1:nlen,jk,jb) = 1._wp

        ! init virtual potential temperature
        ptr_nh_prog%theta_v(1:nlen,jk,jb) = zt0

        ! init exner pressure
        ptr_nh_prog%exner(1:nlen,jk,jb)= (ptr_nh_diag%pres(1:nlen,jk,jb)/p0ref)**rovcp

      ENDDO !jk
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    ! get initial velocity field
    CALL get_nh_df_velocity( ptr_patch, ptr_nh_prog, p_ctest_name, p_rotate_axis_deg, &
      &                     0._wp )

    ! get initial tracer fields for time t=0
    DO jt=1,ntracer
      CALL init_df_tracer( ptr_patch, ptr_int, p_ctest_name,      &! in
        &                  p_rotate_axis_deg, linit_tracer_fv,    &! in
        &                  tracer_inidist_list(jt),               &! in
                           ptr_nh_prog%tracer(:,:,:,jt)           )! inout
    ENDDO

  END SUBROUTINE init_nh_state_prog_dftest



  !---------------------------------------------------------------------------
  !>
  !! SUBROUTINE init_df_velocity
  !!
  !! Short description:
  !! get new horizontal velocity field for deformational flow test
  !! case.
  !!
  SUBROUTINE get_nh_df_velocity( ptr_patch, ptr_prog, p_ctest_name, p_rotate_axis_deg, &
    &                      p_sim_time )

    !INPUT PARAMETERS:
    TYPE(t_patch), INTENT(IN)      :: ptr_patch
    TYPE(t_nh_prog), INTENT(INOUT) :: ptr_prog

    CHARACTER(len=MAX_CHAR_LENGTH), INTENT(IN) :: p_ctest_name

    REAL(wp), INTENT(IN) :: p_rotate_axis_deg !< Earths rotation axis pitch
                                              !< angle in deg
    REAL(wp), INTENT(IN) :: p_sim_time  !< simulation time in seconds
                                        !< since start

    INTEGER  :: nblks_e, nblks_v
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startblk, i_startidx, i_endidx, i_rcstartlev
    INTEGER  :: jb,je,jv,jk         !< loop indices
    REAL(wp) :: z_alpha             !< Earths rotation axis pitch angle in rad
    REAL(wp) :: z_timing_func       !< timing function (responsible for flow
                                    !< reversal)
    REAL(wp) :: u_wind, v_wind      !< zonal and meridional velocity component
    REAL(wp) :: zlon, zlat          !< lon/lat in unrotated system
    REAL(wp) :: zlon_rot, zlat_rot  !< lon/lat in rotated system

    REAL(wp) :: &                   !< streamfunction at vertices
      &  zpsi(nproma,ptr_patch%nblks_v)
    INTEGER  :: ilv1, ilv2, ibv1, ibv2
    REAL(wp) :: iorient             !< system orientation
    
    REAL(wp) ::  u0, u3, u4       !< flow field amplitude [m/s]
    !---------------------------------------------------------------------------
    ! constants
    u0  = (2._wp*pi*grid_sphere_radius)/(tottime)  !< circumference / 12 days [m/s]
    u3  = (5._wp*grid_sphere_radius)/(tottime)     !< for case 3 (divergent flow)
    u4  = (10._wp*grid_sphere_radius)/(tottime)    !< flow field amplitude [m/s]
    ! values for the blocking
    nblks_e  = ptr_patch%nblks_e
    nblks_v  = ptr_patch%nblks_v

   ! number of vertical levels
    nlev = ptr_patch%nlev

    i_rcstartlev = 2
    i_startblk   = ptr_patch%edges%start_blk(i_rcstartlev,1)

    z_alpha = p_rotate_axis_deg*deg2rad


    ! Since this is a 2D test case, the velocity field
    ! only needs to be calculated for 1 full level. Then the
    ! solution is copied to each vertical level.

    ! compute timing function
    z_timing_func = COS(pi*p_sim_time/tottime)


    SELECT CASE (TRIM(p_ctest_name))

      !
      ! Case 1 of Nair (2010)
      !
      CASE ('DF1')
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,zlon,zlat,zlon_rot,zlat_rot, &
!$OMP            u_wind,v_wind)
      DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO je = i_startidx, i_endidx

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  zlon_rot,zlat_rot             ) !<inout
          zlon_rot = zlon
          zlat_rot = zlat

          ! velocity in local east direction
          u_wind = u1 * SIN(0.5_wp*zlon_rot)**2 * SIN(2._wp*zlat_rot) &
            &    * z_timing_func

          ! velocity in local north direction
          v_wind = 0.5_wp * u1 * SIN(zlon_rot) * COS(zlat_rot) * z_timing_func

          ! compute normal wind component
          ptr_prog%vn(je,1,jb) =                                   &
            &     u_wind * ptr_patch%edges%primal_normal(je,jb)%v1 &
            &   + v_wind * ptr_patch%edges%primal_normal(je,jb)%v2


        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      !
      ! Case 2 of Nair (2010)
      !
      CASE ('DF2')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,zlon,zlat,zlon_rot,zlat_rot, &
!$OMP            u_wind,v_wind)
      DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, i_rcstartlev)

        DO je = i_startidx, i_endidx

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  zlon_rot,zlat_rot             ) !<inout

          ! velocity in local east direction
          u_wind = u2 * SIN(zlon_rot)**2 * SIN(2._wp*zlat_rot)     &
            &      * z_timing_func

          ! velocity in local north direction
          v_wind = u2 * SIN(2._wp*zlon_rot) * COS(zlat_rot) * z_timing_func

          ! compute normal wind component
          ptr_prog%vn(je,1,jb) =                                   &
            &     u_wind * ptr_patch%edges%primal_normal(je,jb)%v1 &
            &   + v_wind * ptr_patch%edges%primal_normal(je,jb)%v2

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      !
      ! Case 3 of Nair (2010) (divergent flow)
      !
      CASE ('DF3')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,zlon,zlat,zlon_rot,zlat_rot, &
!$OMP            u_wind,v_wind)
      DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, i_rcstartlev)

        DO je = i_startidx, i_endidx

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  zlon_rot,zlat_rot             ) !<inout

          zlon_rot = zlon_rot - 2._wp*pi*p_sim_time/tottime

          ! velocity in local east direction
          u_wind = -u3 * SIN(0.5_wp*zlon_rot)**2 * SIN(2._wp*zlat_rot) &
            &      * COS(zlat_rot)**2 * z_timing_func + u0*cos(zlat_rot)

          ! velocity in local north direction
          v_wind = 0.5_wp * u3 * SIN(zlon_rot) * COS(zlat_rot)**3  &
            &      * z_timing_func

          ! compute normal wind component
          ptr_prog%vn(je,1,jb) =                                   &
            &     u_wind * ptr_patch%edges%primal_normal(je,jb)%v1 &
            &   + v_wind * ptr_patch%edges%primal_normal(je,jb)%v2

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      !
      ! Case 4 of Nair (2010)
      !
      CASE ('DF4')

!$OMP PARALLEL PRIVATE(i_startblk)

    i_startblk   = ptr_patch%verts%start_blk(i_rcstartlev,1)

!$OMP DO PRIVATE(jb,jv,i_startidx,i_endidx,zlon,zlat,zlon_rot,zlat_rot)
      DO jb = i_startblk, nblks_v

        CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                           i_startidx, i_endidx, i_rcstartlev)

        ! Evaluate stream function at triangle vertices
        DO jv = i_startidx, i_endidx

          ! location of vertices
          zlon = ptr_patch%verts%vertex(jv,jb)%lon
          zlat = ptr_patch%verts%vertex(jv,jb)%lat

          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  zlon_rot,zlat_rot             ) !<inout

          zlon_rot = zlon_rot - 2._wp*pi*p_sim_time/tottime

          zpsi(jv,jb) = u4 * SIN(zlon_rot)**2 * COS(zlat_rot)**2  &
            &           * z_timing_func - u0*sin(zlat_rot)
        ENDDO

      ENDDO
!$OMP END DO


    i_startblk   = ptr_patch%edges%start_blk(i_rcstartlev,1)

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilv1,ilv2,ibv1,ibv2,iorient)
      DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, i_rcstartlev)

        ! deduce normal velocity components at edge midpoints from 
        ! stream function at vertices
        DO je = i_startidx, i_endidx

          ! get edge vertices
          ilv1 = ptr_patch%edges%vertex_idx(je,jb,1)
          ibv1 = ptr_patch%edges%vertex_blk(je,jb,1)
          ilv2 = ptr_patch%edges%vertex_idx(je,jb,2)
          ibv2 = ptr_patch%edges%vertex_blk(je,jb,2)

          ! compute normal wind component
          ! this is done by computing the tangential derivative of the 
          ! streamfunction.
          ! Note that the tangential direction is defined by
          ! iorient*(vertex2 - vertex1)
          !
          ! DR: I am not quite sure, why I have to multiply with -1 once 
          ! again, in order to get the correct result. This may have 
          ! something to do with the fact, that for the system orientation 
          ! the vector product between the dual normal and the primal normal 
          ! (dn x pn) is computed and not (pn x dn). 
          iorient = ptr_patch%edges%tangent_orientation(je,jb)

          ptr_prog%vn(je,1,jb) = grid_sphere_radius * iorient                             &
            &                  * (zpsi(ilv2,ibv2) - zpsi(ilv1,ibv1))      &
            &                  / ptr_patch%edges%primal_edge_length(je,jb)
        ENDDO

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SELECT


    ! copy to all vertical levels
    DO jk = 2, nlev
      ptr_prog%vn(:,jk,:) =  ptr_prog%vn(:,jk-1,:)
    END DO

  END SUBROUTINE get_nh_df_velocity


  !-----------------------------------------------------------------------
  !!>
  !! SUBROUTINE init_df_tracer
  !!
  !! Short description:
  !! Initialization of tracer distribution for deformational flow
  !! test. Independent of the chosen flow field, the user can choose 
  !! between a C1 cosine bell, a slotted cylinder and a C_infty gaussian 
  !! hill.
  !!
  SUBROUTINE init_df_tracer( ptr_patch, ptr_int, p_ctest_name,                   &
    &                        p_rotate_axis_deg, linit_tracer_fv, tracer_inidist, &
    &                        p_cc )

    !INPUT PARAMETERS:
    TYPE(t_patch), INTENT(IN)     :: &
      &  ptr_patch
    TYPE(t_int_state), INTENT(IN) :: &
      &  ptr_int
    REAL(wp), INTENT(IN)          :: & !< Earths rotation axis pitch
      &  p_rotate_axis_deg             !< angle in deg
                                             
    CHARACTER(len=MAX_CHAR_LENGTH), INTENT(IN) :: &
      &  p_ctest_name
    LOGICAL, INTENT(IN)           :: & !< fv init. for tracer fields
      &  linit_tracer_fv    
    INTEGER, INTENT(IN)           :: & !< selected initial tracer distribution
      &  tracer_inidist   
    REAL(wp), INTENT(INOUT)       :: & !< tracer array
      &  p_cc(:,:,:)    

    INTEGER  :: jb,jc,jk            !> loop indices
    INTEGER  :: nblks_c
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startblk, i_startidx, i_endidx, i_rcstartlev
    REAL(wp) :: z_alpha             !< Earths rotation axis pitch angle in rad
    REAL(wp) :: zlon, zlat          !< lon/lat of cell centers (unrotated)
    REAL(wp) :: z_dist_1, z_dist_2  !< distance from ic-center
    REAL(wp) :: d1,d2               !< normalized distance from ic-center
    REAL(wp) :: zcos_bell           !< cosine bell dummy field
    TYPE(t_geographical_coordinates) :: & !< gg. coords of ic-centers
      &  ic_c1, ic_c2, gc_rot             !< rotated gg coords. of current point 
    TYPE(t_cartesian_coordinates) :: &  !< cart. coords of ic-centers
      &  ic_cc_c1, ic_cc_c2,         &
      &  cc_rot                         !< rotated point in cart. coords
    REAL(wp) :: bell_radius  !< tracer distribution
    !-----------------------------------------------------------------------
    ! constants
    bell_radius = 0.5_wp * grid_sphere_radius  !< tracer distribution
    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c

    ! number of vertical levels
    nlev = ptr_patch%nlev

    i_rcstartlev = 2
    i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)


    z_alpha = p_rotate_axis_deg*deg2rad

    ! Set ic-centers depending on the chosen testcase
    SELECT CASE( TRIM(p_ctest_name) )

      ! Case 1 of Nair (2010)
      CASE('DF1')
        ic_c1%lon = pi
        ic_c1%lat = pi/3._wp

        ic_c2%lon = pi
        ic_c2%lat = -pi/3._wp

      ! Case 2 of Nair (2010)
      CASE('DF2')
        ic_c1%lon = (5._wp/6._wp)*pi
        ic_c1%lat = 0._wp

        ic_c2%lon = (7._wp/6._wp)*pi
        ic_c2%lat = 0._wp

      ! Case 3 of Nair (2010)
      CASE('DF3')
!!$        ic_c1%lon = (3._wp/4._wp)*pi
!!$        ic_c1%lat = 0._wp
!!$
!!$        ic_c2%lon = (5._wp/4._wp)*pi
!!$        ic_c2%lat = 0._wp

        ! Following Lauritzen (2011), JCP
        ic_c1%lon = (5._wp/6._wp)*pi
        ic_c1%lat = 0._wp

        ic_c2%lon = (7._wp/6._wp)*pi
        ic_c2%lat = 0._wp

      ! Case 4 of Nair (2010)
      CASE('DF4')
        ic_c1%lon = (5._wp/6._wp)*pi
        ic_c1%lat = 0._wp

        ic_c2%lon = (7._wp/6._wp)*pi
        ic_c2%lat = 0._wp

    END SELECT


    ! Since this is a 2D test case, the analytic tracer distribution
    ! only needs to be computed for 1 full level. Then the solution
    ! is copied to each vertical level.
    !
    SELECT CASE(tracer_inidist)
    !
    ! C1 cosine bells
    !
    CASE(5)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,gc_rot, &
!$OMP            z_dist_1,z_dist_2,d1,d2)
      DO jb = i_startblk, nblks_c

        CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO jc = i_startidx, i_endidx

          ! location of quadrature points
          IF (linit_tracer_fv) THEN
            zlon = ptr_int%gquad%qpts_tri_l(jc,jb)%lon
            zlat = ptr_int%gquad%qpts_tri_l(jc,jb)%lat
          ELSE
            zlon = ptr_patch%cells%center(jc,jb)%lon
            zlat = ptr_patch%cells%center(jc,jb)%lat
          ENDIF


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  gc_rot%lon,gc_rot%lat         ) !<inout

          ! initial distribution
          ! distance (arclength) to ic-center1
          z_dist_1 = grid_sphere_radius * &
            & ACOS( SIN(ic_c1%lat) * SIN(gc_rot%lat) + COS(ic_c1%lat)   &
            & * COS(gc_rot%lat) * COS(gc_rot%lon - ic_c1%lon) )

          ! distance (arclength) to ic-center2
          z_dist_2 = grid_sphere_radius * &
            & ACOS( SIN(ic_c2%lat) * SIN(gc_rot%lat) + COS(ic_c2%lat)   &
            & * COS(gc_rot%lat) * COS(gc_rot%lon - ic_c2%lon) )

          d1 = z_dist_1/bell_radius
          d2 = z_dist_2/bell_radius

          IF ( d1 < 1._wp ) THEN
            p_cc(jc,1,jb) = b_bell + c_bell * (0.5_wp * ( 1._wp + COS(pi*d1) ) )
          ELSE IF ( d2 < 1._wp ) THEN
            p_cc(jc,1,jb) = b_bell + c_bell * (0.5_wp * ( 1._wp + COS(pi*d2) ) )
          ELSE
            p_cc(jc,1,jb) = b_bell
          ENDIF

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    !
    ! slotted cylinders
    ! 
    CASE(6)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,gc_rot, &
!$OMP            z_dist_1,z_dist_2)
      DO jb = i_startblk, nblks_c

        CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO jc = i_startidx, i_endidx

          ! location of quadrature points
          IF (linit_tracer_fv) THEN
            zlon = ptr_int%gquad%qpts_tri_l(jc,jb)%lon
            zlat = ptr_int%gquad%qpts_tri_l(jc,jb)%lat
          ELSE
            zlon = ptr_patch%cells%center(jc,jb)%lon
            zlat = ptr_patch%cells%center(jc,jb)%lat
          ENDIF


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  gc_rot%lon,gc_rot%lat         ) !<inout


          ! initial distribution
          ! distance (arclength) to ic-center1
          z_dist_1 = grid_sphere_radius * &
            &  ACOS( SIN(ic_c1%lat) * SIN(gc_rot%lat) + COS(ic_c1%lat)   &
            &  * COS(gc_rot%lat) * COS(gc_rot%lon - ic_c1%lon) )

          ! distance (arclength) to ic-center2
          z_dist_2 = grid_sphere_radius * &
            & ACOS( SIN(ic_c2%lat) * SIN(gc_rot%lat) + COS(ic_c2%lat)   &
            & * COS(gc_rot%lat) * COS(gc_rot%lon - ic_c2%lon) )


          IF (z_dist_1 <= bell_radius &
            & .AND. ABS(gc_rot%lon - ic_c1%lon)>= bell_radius/(6._wp*grid_sphere_radius)) &
            & THEN
              p_cc(jc,1,jb) = c_slotted
          ELSE IF (z_dist_2 <= bell_radius                           & 
            & .AND. ABS(gc_rot%lon - ic_c2%lon)>= bell_radius/(6._wp*grid_sphere_radius)) &
            & THEN
              p_cc(jc,1,jb) = c_slotted
          ELSE IF (z_dist_1 <= bell_radius                           &
            & .AND. ABS(gc_rot%lon - ic_c1%lon)< bell_radius/(6._wp*grid_sphere_radius)   &
            & .AND. (gc_rot%lat-ic_c1%lat) < -5._wp*bell_radius/(12._wp*grid_sphere_radius)) &
            & THEN
              p_cc(jc,1,jb) = c_slotted
          ELSE IF(z_dist_2 <= bell_radius                            &
            & .AND. ABS(gc_rot%lon - ic_c2%lon) < bell_radius/(6._wp*grid_sphere_radius)  &
            & .AND. (gc_rot%lat-ic_c2%lat) > 5._wp*bell_radius/(12._wp*grid_sphere_radius)) &
            &  THEN
              p_cc(jc,1,jb) = c_slotted
          ELSE
              p_cc(jc,1,jb) = b_slotted
          ENDIF  

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    !
    ! C_infty gaussian hills (see Levy et. al., 2007)
    ! 
    CASE(7)

    ! Transform ic_center to cartesian coordinates
    ic_cc_c1 = gc2cc(ic_c1)
    ic_cc_c2 = gc2cc(ic_c2)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,gc_rot,cc_rot, &
!$OMP            d1,d2)
      DO jb = i_startblk, nblks_c

        CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,  &
                           i_startidx, i_endidx, i_rcstartlev)


        DO jc = i_startidx, i_endidx

          ! location of quadrature points
          IF (linit_tracer_fv) THEN
            zlon = ptr_int%gquad%qpts_tri_l(jc,jb)%lon
            zlat = ptr_int%gquad%qpts_tri_l(jc,jb)%lat
          ELSE
            zlon = ptr_patch%cells%center(jc,jb)%lon
            zlat = ptr_patch%cells%center(jc,jb)%lat
          ENDIF


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  gc_rot%lon,gc_rot%lat         ) !<inout

          ! Transform geographical coords to cartesian coords
          cc_rot = gc2cc(gc_rot)

          !
          ! get distance between ic center and current point in space
          !
          ! distance to ic-center1
          d1 = (cc_rot%x(1)-ic_cc_c1%x(1))**2 + (cc_rot%x(2)-ic_cc_c1%x(2))**2  &
            & + (cc_rot%x(3)-ic_cc_c1%x(3))**2

          ! distance to ic-center2
          d2 = (cc_rot%x(1)-ic_cc_c2%x(1))**2 + (cc_rot%x(2)-ic_cc_c2%x(2))**2  &
            & + (cc_rot%x(3)-ic_cc_c2%x(3))**2

          !
          ! initial distribution
          !
          p_cc(jc,1,jb) = c_gauss * ( exp(-b_gauss*d1) + exp(-b_gauss*d2) )

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL


    !
    ! linearly correlated C1 cosine bells
    !
    CASE(8)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,gc_rot, &
!$OMP            z_dist_1,z_dist_2,d1,d2,zcos_bell)
      DO jb = i_startblk, nblks_c

        CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO jc = i_startidx, i_endidx

          ! location of quadrature points
          IF (linit_tracer_fv) THEN
            zlon = ptr_int%gquad%qpts_tri_l(jc,jb)%lon
            zlat = ptr_int%gquad%qpts_tri_l(jc,jb)%lat
          ELSE
            zlon = ptr_patch%cells%center(jc,jb)%lon
            zlat = ptr_patch%cells%center(jc,jb)%lat
          ENDIF


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  gc_rot%lon,gc_rot%lat         ) !<inout

          ! initial distribution
          ! distance (arclength) to ic-center1
          z_dist_1 = grid_sphere_radius * &
            & ACOS( SIN(ic_c1%lat) * SIN(gc_rot%lat) + COS(ic_c1%lat)   &
            &  * COS(gc_rot%lat) * COS(gc_rot%lon - ic_c1%lon) )

          ! distance (arclength) to ic-center2
          z_dist_2 = grid_sphere_radius * &
            & ACOS( SIN(ic_c2%lat) * SIN(gc_rot%lat) + COS(ic_c2%lat)   &
            & * COS(gc_rot%lat) * COS(gc_rot%lon - ic_c2%lon) )

          d1 = z_dist_1/bell_radius
          d2 = z_dist_2/bell_radius

          IF ( d1 < 1._wp ) THEN
            zcos_bell = b_bell + c_bell * (0.5_wp * ( 1._wp + COS(pi*d1) ) )
          ELSE IF ( d2 < 1._wp ) THEN
            zcos_bell = b_bell + c_bell * (0.5_wp * ( 1._wp + COS(pi*d2) ) )
          ELSE
            zcos_bell = b_bell
          ENDIF

          p_cc(jc,1,jb) = b_lin_bell + c_lin_bell*zcos_bell

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL


    !
    ! nonlinearly correlated cosine bells
    !
    CASE(9)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,gc_rot, &
!$OMP            z_dist_1,z_dist_2,d1,d2,zcos_bell)
      DO jb = i_startblk, nblks_c

        CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO jc = i_startidx, i_endidx

          ! location of quadrature points
          IF (linit_tracer_fv) THEN
            zlon = ptr_int%gquad%qpts_tri_l(jc,jb)%lon
            zlat = ptr_int%gquad%qpts_tri_l(jc,jb)%lat
          ELSE
            zlon = ptr_patch%cells%center(jc,jb)%lon
            zlat = ptr_patch%cells%center(jc,jb)%lat
          ENDIF


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  gc_rot%lon,gc_rot%lat         ) !<inout

          ! initial distribution
          ! distance (arclength) to ic-center1
          z_dist_1 = grid_sphere_radius * &
            & ACOS( SIN(ic_c1%lat) * SIN(gc_rot%lat) + COS(ic_c1%lat)   &
            &    * COS(gc_rot%lat) * COS(gc_rot%lon - ic_c1%lon) )

          ! distance (arclength) to ic-center2
          z_dist_2 = grid_sphere_radius * &
            & ACOS( SIN(ic_c2%lat) * SIN(gc_rot%lat) + COS(ic_c2%lat)   &
            &    * COS(gc_rot%lat) * COS(gc_rot%lon - ic_c2%lon) )

          d1 = z_dist_1/bell_radius
          d2 = z_dist_2/bell_radius

          IF ( d1 < 1._wp ) THEN
            zcos_bell = b_bell + c_bell * (0.5_wp * ( 1._wp + COS(pi*d1) ) )
          ELSE IF ( d2 < 1._wp ) THEN
            zcos_bell = b_bell + c_bell * (0.5_wp * ( 1._wp + COS(pi*d2) ) )
          ELSE
            zcos_bell = b_bell
          ENDIF

          p_cc(jc,1,jb) = b_nlin_bell + c_nlin_bell*zcos_bell**2

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SELECT 

    ! copy to all vertical levels
    DO jk = 2, nlev
        p_cc(:,jk,:) =  p_cc(:,jk-1,:)
    ENDDO

  END SUBROUTINE init_df_tracer



  !---------------------------------------------------------------------------
  !>
  !! SUBROUTINE get_departure_points
  !!
  !! Short description:
  !! Compute backward trajectories and corresponding departure points based
  !! on a third order Taylor series approximation. All backward trajectories
  !! start at edge midpoints.
  !!
  !! To be more precise, instead of the departure points, the barycenter of
  !! the departure region is calculated. This particular point is needed for the
  !! MIURA scheme, when combined with a second order reconstruction. Once
  !! the barycenters are known in geographical coordinates, they are projected
  !! onto a tangent plane (tangent to the circumcenter of the upwind cell) using
  !! the gnomonic projection. Then the distance vector between the circumcenter
  !! of the upwind cell and the barycenter of the departure region is computed
  !! and stored. These vectors may directly be applied to our Miura-scheme in
  !! order to assess the importance of accurate departure points. Moreover these
  !! departure points may directly be compared to the departure points derived
  !! with our operational methods.
  !!
  SUBROUTINE get_departure_points( ptr_patch, ptr_prog, p_ctest_name, p_rotate_axis_deg, &
    &                      p_dtime, p_sim_time )

    !INPUT PARAMETERS:
    TYPE(t_patch), INTENT(IN)      :: ptr_patch
    TYPE(t_nh_prog), INTENT(INOUT) :: ptr_prog

    CHARACTER(len=MAX_CHAR_LENGTH), INTENT(IN) :: p_ctest_name

    REAL(wp), INTENT(IN) :: p_rotate_axis_deg !< Earths rotation axis pitch
                                              !< angle in deg
    REAL(wp), INTENT(IN) :: p_dtime     !< time step

    REAL(wp), INTENT(IN) :: p_sim_time  !< simulation time in seconds
                                        !< since start

    INTEGER  :: nblks_e
    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startblk, i_startidx, i_endidx, i_rcstartlev
    INTEGER  :: jb,je,jk            !< loop indices
    INTEGER  :: ilc0, ibc0          !< line and block index of upwind cell
    REAL(wp) :: z_alpha             !< Earths rotation axis pitch angle in rad
    REAL(wp) :: z_timing_func       !< timing function (responsible for flow
                                    !< reversal)
    REAL(wp) :: z_dthalf            !< 0.5 pdtime
    REAL(wp) :: u_wind, v_wind      !< zonal and meridional velocity component
    REAL(wp) :: zlon, zlat          !< lon/lat in unrotated system
    REAL(wp) :: zlon_rot, zlat_rot  !< lon/lat in rotated system
    REAL(wp) :: c_1, c_2, c_3, c_4  !< terms of Taylor series approx.
    REAL(wp) :: xloc, yloc
    REAL(wp) :: z_dist_g(2)
    REAL(wp) ::   &                 !< geographical coordinates of departure region
     &  z_barycenter(nproma,ptr_patch%nlev,ptr_patch%nblks_e,2) !< barycenter

    REAL(wp) :: u3, u4      !< flow field amplitude [m/s]
    !---------------------------------------------------------------------------
    ! constants
    u3  = (5._wp*grid_sphere_radius)/(tottime)     !< for case 3 (divergent flow)
    u4  = (10._wp*grid_sphere_radius)/(tottime)    !< flow field amplitude [m/s]

    ! values for the blocking
    nblks_e  = ptr_patch%nblks_e

    ! number of vertical levels
    nlev = ptr_patch%nlev

    i_rcstartlev = 2
    i_startblk   = ptr_patch%edges%start_blk(i_rcstartlev,1)

    z_alpha = p_rotate_axis_deg*deg2rad

    ! compute timing function
    z_timing_func = COS(pi*p_sim_time/tottime)

    ! since we are interested in the barycenter and not the
    ! departure point, we use:
    z_dthalf = 0.5_wp * p_dtime



    ! Since this is a 2D test case, the departure points only need
    ! to be calculated for 1 full level. Then the solution is copied
    ! to each vertical level.
    !
    SELECT CASE (TRIM(p_ctest_name))

      !
      ! Case 1 of Nair (2010)
      !
      CASE ('DF1')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,zlon,zlat,zlon_rot,zlat_rot, &
!$OMP            u_wind,v_wind,c_1,c_2,c_3,c_4)
      DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO je = i_startidx, i_endidx

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  zlon_rot,zlat_rot             ) !<inout


          ! zonal velocity component at edge midpoint
          u_wind = u1 * SIN(0.5_wp*zlon_rot) * SIN(0.5_wp*zlon_rot)  &
            &      * SIN(2._wp*zlat_rot) * z_timing_func

          ! meridional velocity component at edge midpoint
          v_wind = 0.5_wp * u1 * SIN(zlon_rot) * COS(zlat_rot) * z_timing_func

          c_1 = u_wind/(grid_sphere_radius * COS(zlat_rot))

          c_3 = v_wind/grid_sphere_radius

          c_2 = (u1/grid_sphere_radius) * SIN(0.5_wp*zlon_rot) * ( -SIN(0.5_wp*zlon_rot)   &
            & * SIN(zlat_rot) * SIN(pi*p_sim_time/tottime) * (pi/tottime)  &
            & + z_timing_func * (c_3 * COS(zlat_rot)*SIN(0.5_wp*zlon_rot)  &
            & + COS(0.5_wp*zlon_rot) * SIN(zlat_rot) * c_1) )

          c_4 = (0.25_wp*u1/grid_sphere_radius) * (-SIN(zlon_rot) * COS(zlat_rot)          &
            & * SIN(pi*p_sim_time/tottime) * (pi/tottime)                  &
            & + z_timing_func * (COS(zlat_rot) * COS(zlon_rot)             &
            & * c_1 - c_3 * SIN(zlon_rot) * SIN(zlat_rot)) )

          ! longitudinal position of barycenter
          z_barycenter(je,1,jb,1) = zlon_rot - z_dthalf * c_1    &
            &                     + z_dthalf * z_dthalf * c_2

          ! latitudinal position of barycenter
          z_barycenter(je,1,jb,2) = zlat_rot - z_dthalf * c_3    &
            &                     + z_dthalf * z_dthalf * c_4

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      !
      ! Case 2 of Nair (2010)
      !
      CASE ('DF2')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,zlon,zlat,zlon_rot,zlat_rot, &
!$OMP            u_wind,v_wind,c_1,c_2,c_3,c_4)
      DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO je = i_startidx, i_endidx

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  zlon_rot,zlat_rot             ) !<inout


          ! zonal velocity component at edge midpoint
          u_wind = u2 * SIN(zlon_rot) * SIN(zlon_rot) * SIN(2._wp*zlat_rot) &
            &      * z_timing_func

          ! velocity in local north direction
          v_wind = u2 * SIN(2._wp*zlon_rot) * COS(zlat_rot) * z_timing_func

          c_1 = u_wind/(grid_sphere_radius * COS(zlat_rot))

          c_3 = v_wind/grid_sphere_radius

          c_2 = (u2/grid_sphere_radius) * SIN(zlon_rot) * ( -SIN(zlon_rot)  &
            & * SIN(zlat_rot) * SIN(pi*p_sim_time/tottime) * (pi/tottime)  &
            & + z_timing_func * (c_3 * COS(zlat_rot)*SIN(zlon_rot)         &
            & + 2._wp * COS(zlon_rot) * SIN(zlat_rot) * c_1) )

          c_4 = (0.5_wp*u2/grid_sphere_radius) * (-SIN(2._wp*zlon_rot) * COS(zlat_rot) &
            & * SIN(pi*p_sim_time/tottime) * (pi/tottime)                        &
            & + z_timing_func * (2._wp*COS(zlat_rot) * COS(2._wp*zlon_rot) * c_1 &
            & - c_3 * SIN(2._wp*zlon_rot) * SIN(zlat_rot)) )

          ! longitudinal position of barycenter
          z_barycenter(je,1,jb,1) = zlon_rot - z_dthalf * c_1      &
            &                     + z_dthalf * z_dthalf * c_2

          ! latitudinal position of barycenter
          z_barycenter(je,1,jb,2) = zlat_rot - z_dthalf * c_3      &
            &                     + z_dthalf * z_dthalf * c_4

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL


      !
      ! Case 3 of Nair (2010)
      !
      CASE ('DF3')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,zlon,zlat,zlon_rot,zlat_rot, &
!$OMP            u_wind,v_wind,c_1,c_2,c_3,c_4)
      DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO je = i_startidx, i_endidx

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  zlon_rot,zlat_rot             ) !<inout


          ! zonal velocity component at edge midpoint
          u_wind = -u3 * SIN(0.5_wp*zlon_rot) * SIN(0.5_wp*zlon_rot)  &
            &      * SIN(2._wp*zlat_rot) * z_timing_func

          ! meridional velocity component at edge midpoint
          v_wind = 0.5_wp * u3 * SIN(zlon_rot) * COS(zlat_rot) * z_timing_func

          c_1 = u_wind/(grid_sphere_radius * COS(zlat_rot))

          c_3 = v_wind/grid_sphere_radius

          c_2 = (u3/grid_sphere_radius) * SIN(0.5_wp*zlon_rot) * ( SIN(0.5_wp * zlon_rot) &
            & * SIN(zlat_rot) * SIN(pi*p_sim_time/tottime) * (pi/tottime)  &
            & - z_timing_func * (c_3 * COS(zlat_rot)*SIN(0.5_wp*zlon_rot)  &
            & + COS(0.5_wp*zlon_rot) * SIN(zlat_rot) * c_1) )

          c_4 = (0.25_wp*u3/grid_sphere_radius) * (-SIN(zlon_rot) * COS(zlat_rot)       &
            & * SIN(pi*p_sim_time/tottime) * (pi/tottime)               &
            & + z_timing_func * (COS(zlat_rot) * COS(zlon_rot) * c_1    &
            & - c_3 * SIN(zlon_rot) * SIN(zlat_rot)) )

          ! longitudinal position of barycenter
          z_barycenter(je,1,jb,1) = zlon_rot - z_dthalf * c_1      &
            &                     + z_dthalf * z_dthalf * c_2

          ! latitudinal position of barycenter
          z_barycenter(je,1,jb,2) = zlat_rot - z_dthalf * c_3      &
            &                     + z_dthalf * z_dthalf * c_4

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      !
      ! Case 4 of Nair (2010)
      !
      CASE ('DF4')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,zlon,zlat,zlon_rot,zlat_rot, &
!$OMP            u_wind,v_wind,c_1,c_2,c_3,c_4)
      DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO je = i_startidx, i_endidx

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat


          ! get zlat, zlon in rotated system with npole at
          ! (npole_lon,npole_lat)
          CALL rotated_sphere( npole_lon,npole_lat-z_alpha,  & !<in
            &                  zlon,zlat,                    & !<in
            &                  zlon_rot,zlat_rot             ) !<inout


          ! zonal velocity component at edge midpoint
          u_wind = - u4 * SIN(0.5_wp*zlon_rot) * SIN(0.5_wp*zlon_rot)  &
            &      * SIN(2._wp*zlat_rot) * z_timing_func

          ! meridional velocity component at edge midpoint
          v_wind = - u4 * SIN(zlon_rot) * COS(zlat_rot) * COS(zlat_rot) &
            &      * z_timing_func

          c_1 = u_wind/(grid_sphere_radius * COS(zlat_rot))

          c_3 = v_wind/grid_sphere_radius

          c_2 = (u4/grid_sphere_radius) * SIN(0.5_wp*zlon_rot) * ( SIN(0.5_wp * zlon_rot)&
            & * SIN(zlat_rot) * SIN(pi*p_sim_time/tottime) * (pi/tottime)       &
            & - z_timing_func * (c_3 * COS(zlat_rot)*SIN(0.5_wp*zlon_rot)       &
            & + COS(0.5_wp*zlon_rot) * SIN(zlat_rot) * c_1) )

          c_4 = (0.5_wp*u4/grid_sphere_radius) * (SIN(zlon_rot) * COS(zlat_rot) &
            & * COS(zlat_rot) * SIN(pi*p_sim_time/tottime) * (pi/tottime)            &
            & - z_timing_func * (COS(zlat_rot) * COS(zlat_rot) * COS(zlon_rot)  &
            & * c_1 - SIN(2._wp*zlat_rot) * SIN(zlon_rot) * c_3 ) )

          ! longitudinal position of barycenter
          z_barycenter(je,1,jb,1) = zlon_rot - z_dthalf * c_1      &
            &                     + z_dthalf * z_dthalf * c_2

          ! latitudinal position of barycenter
          z_barycenter(je,1,jb,2) = zlat_rot - z_dthalf * c_3      &
            &                     + z_dthalf * z_dthalf * c_4

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SELECT


    ! project barycenter onto tangent plane and compute distance vector
    ! circumcenter of cell --> barycenter of departure region.
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc0,ibc0,xloc,yloc,z_dist_g)
    DO jb = i_startblk, nblks_e

       CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, i_rcstartlev)

       DO je = i_startidx, i_endidx
         ! Determine the cell in which the barycenter is located (upwind cell)
         ! the cell indices are chosen such that the direction from cell 1 to cell 2
         ! is the positive direction of the normal vector N
         IF (ptr_prog%vn(je,1,jb) >= 0._wp) THEN
           ! line and block indices of neighboring cell with barycenter
           ilc0 = ptr_patch%edges%cell_idx(je,jb,1)
           ibc0 = ptr_patch%edges%cell_blk(je,jb,1)
         ELSE
           ! line and block indices of neighboring cell with barycenter
           ilc0 = ptr_patch%edges%cell_idx(je,jb,2)
           ibc0 = ptr_patch%edges%cell_blk(je,jb,2)
         ENDIF
         ! save line and block indices
         df_cell_indices(je,1,jb,1) = ilc0
         df_cell_indices(je,1,jb,2) = ibc0

         ! calculate distance vector between barycenter and cell center
         ! in geographical coordinates
         ! Apply gnomonic projection. (Projection onto tangential plane touching
         ! the sphere at the cell center)
         xloc = ptr_patch%cells%center(ilc0,ibc0)%lon
         yloc = ptr_patch%cells%center(ilc0,ibc0)%lat

!!$         CALL az_eqdist_proj( xloc, yloc, z_barycenter(je,1,jb,1), & !in
!!$           &                 z_barycenter(je,1,jb,2),             & !in
!!$           &                 z_dist_g(1), z_dist_g(2) )             !out
         CALL gnomonic_proj( xloc, yloc, z_barycenter(je,1,jb,1), & !in
           &                 z_barycenter(je,1,jb,2),             & !in
           &                 z_dist_g(1), z_dist_g(2) )             !out

         ! Save distance vector cell center -> barycenter of advected area
         ! in an edge based data structure. This is possible, since one
         ! unique distance vector corresponds to one edge.
         df_distv_barycenter(je,1,jb,1:2) = grid_sphere_radius * z_dist_g(1:2)

      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! copy results to all vertical levels
    DO jk = 2, nlev
      df_cell_indices(:,jk,:,:)     = df_cell_indices(:,jk-1,:,:)
      df_distv_barycenter(:,jk,:,:) = df_distv_barycenter(:,jk-1,:,:)
    END DO

  END SUBROUTINE get_departure_points



  !-----------------------------------------------------------------------
  !!>
  !! SUBROUTINE prep_departure_points_err
  !!
  !! Short description:
  !! prepares output of error for departure points.
  !!
  !! For each cell the difference between the analytic and approximate
  !! distance vector is computed by summing over the edges.
  !!
    SUBROUTINE prep_departure_points_err( ptr_patch, p_cc)

    !INPUT PARAMETERS:
    TYPE(t_patch), TARGET, INTENT(IN) :: ptr_patch

    REAL(wp), INTENT(INOUT) :: p_cc(:,:,:)   !< tracer array

    INTEGER  :: jb,jc,jk                     !< loop indices
    INTEGER  :: nblks_c
    INTEGER  :: nlev                         !< number of full levels
    INTEGER  :: i_startblk, i_startidx, i_endidx, i_rcstartlev

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block
      &  iile, iibe                          !< indices (array)

    !-----------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c

    ! number of vertical levels
    nlev = ptr_patch%nlev

    i_rcstartlev = 2
    i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)

    ! pointer to line and block indices of cell edges
    iile => ptr_patch%cells%edge_idx
    iibe => ptr_patch%cells%edge_blk

    DO jb = i_startblk, nblks_c

      CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,  &
                         i_startidx, i_endidx, i_rcstartlev)

      DO jc = i_startidx, i_endidx


        p_cc(jc,1,jb) = df_distv_barycenter(iile(jc,jb,1),1,iibe(jc,jb,1),1)  &
          &           + df_distv_barycenter(iile(jc,jb,2),1,iibe(jc,jb,2),1)  &
          &           + df_distv_barycenter(iile(jc,jb,3),1,iibe(jc,jb,3),1)

      ENDDO
    ENDDO


    ! copy to all vertical levels
    DO jk = 2, nlev
        p_cc(:,jk,:) =  p_cc(:,jk-1,:)
    ENDDO


  END SUBROUTINE prep_departure_points_err



  !---------------------------------------------------------------------------
  !!>
  !! SUBROUTINE rotated_sphere
  !!
  !! Short description:
  !! Rotate a sphere so that its NP is at (lc,tc). The rotated coordinates
  !!(rol,rot) corresponding to (la,th) in the regular sph. coordinates
  !!
  SUBROUTINE rotated_sphere(lc,tc,la,th,rol,rot)

    IMPLICIT NONE
    REAL(wp), INTENT(IN)  :: lc,tc,la,th
    REAL(wp), INTENT(OUT) :: rol,rot

    REAL(wp) :: cost,sint,sinc,cosc
    REAL(wp) :: trm,trm1,trm2,trm3

    !---------------------------------------------------------------------------

    sinc = sin(tc)
    cosc = cos(tc)
    cost = cos(th)
    sint = sin(th)

    trm  = cost * cos(la - lc)
    trm1 = cost * sin(la - lc)
    trm2 = sinc * trm  - cosc * sint
    trm3 = sinc * sint + cosc * trm

    rol = atan2(trm1,trm2)
    IF ( rol > 2._wp*pi  ) rol = rol - 2._wp*pi
    IF ( rol < 0._wp     ) rol = rol + 2._wp*pi
    rot = asin(trm3)

  END SUBROUTINE rotated_sphere


END MODULE mo_nh_df_test

