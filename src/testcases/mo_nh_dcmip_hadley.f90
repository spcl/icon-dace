! Subroutine to initialize testcase 12 (Hadley-like meridional circulation) proposed
! for the DCMIP summer school
!
! Literature
! - Dynamical Core Model Intercomparison Project (DCMIP)
!   Test Case Document (P. Ullrich et al, 2012)
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

MODULE mo_nh_dcmip_hadley

   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd, grav, p0ref, cpd, cvd_o_rd
   USE mo_math_constants,       ONLY: pi
   USE mo_impl_constants,       ONLY: min_rlcell, min_rledge
   USE mo_parallel_config,      ONLY: nproma
   USE mo_loopindices,          ONLY: get_indices_c, get_indices_e
   USE mo_model_domain,         ONLY: t_patch
   USE mo_grid_config,          ONLY: grid_sphere_radius
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_intp_data_strc,       ONLY: t_int_state
   USE mo_intp,                 ONLY: cells2edges_scalar
   USE mo_sync,                 ONLY: sync_patch_array, SYNC_E
   USE mo_fortran_tools,        ONLY: set_acc_host_or_device

   IMPLICIT NONE


   PRIVATE

   PUBLIC :: init_nh_dcmip_hadley
   PUBLIC :: set_nh_velocity_hadley


   ! test case parameters
   !
   REAL(wp), PARAMETER :: t0    = 300.0_wp    !< temperature           [K]
   REAL(wp), PARAMETER :: scale_hgt = rd*t0/grav !< scale height       [m] 
   REAL(wp), PARAMETER :: u0    = 40.0_wp     !< maximum amplitude 
                                              !< of the zonal wind     [m s^-1]
   REAL(wp), PARAMETER :: w0    = 0.15_wp     !< maximum amplitude 
                                              !< of the vertical wind  [m s^-1]
   REAL(wp), PARAMETER :: z1    = 2000.0_wp   !< lower tracer bound    [m]
   REAL(wp), PARAMETER :: z2    = 5000.0_wp   !< upper tracer bound    [m]
   REAL(wp), PARAMETER :: z_mid = 0.5_wp*(z1+z2) !< midpoint              [m]
   REAL(wp), PARAMETER :: ztop  = 12000._wp   !< model top             [m]
   REAL(wp), PARAMETER :: tau   = 1.0_wp * 86400.0_wp ! period of motion 1 day [s]
   INTEGER , PARAMETER :: nhadley = 5         !< number of overturning hadley cells

!--------------------------------------------------------------------

CONTAINS
!-------------------------------------------------------------------------
!

!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state vector for the DCMIP Hadley-like 
  !! meridional circulation test 
  !!
  SUBROUTINE init_nh_dcmip_hadley( p_patch, p_nh_prog, p_nh_diag, &
    &                         p_int, p_metrics )

    TYPE(t_patch),        INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_int_state),    INTENT(IN)    :: &  !< interpolation state vector
      &  p_int

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics


    INTEGER  :: jc, jk, jb                    !< loop indices
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom  
    INTEGER  :: nlev, nlevp1                  !< number of full/half levels
    INTEGER  :: ntracer_alloc                 !< number of allocated tracer fields

 !--------------------------------------------------------------------


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1


    i_rlstart = 1
    i_rlend   = min_rlcell

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


    ntracer_alloc=SIZE(p_nh_prog%tracer,4)

    !
    ! Init prognostic variables
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! init tracer field
          !
          IF (p_metrics%z_mc(jc,jk,jb) < z2 .AND. p_metrics%z_mc(jc,jk,jb) > z1) THEN
            p_nh_prog%tracer(jc,jk,jb,1) = 0.5_wp * (1.0_wp + cos( 2.0_wp*pi      &
              &                          *(p_metrics%z_mc(jc,jk,jb)-z_mid)/(z2-z1) ) )
          ELSE
            p_nh_prog%tracer(jc,jk,jb,1) = 0.0_wp
          ENDIF



          ! temperature is constant
          !
          p_nh_diag%temp(jc,jk,jb) = t0


          ! init pressure field
          !
          p_nh_diag%pres(jc,jk,jb) = p0ref * exp(-p_metrics%z_mc(jc,jk,jb)/scale_hgt)


          ! init exner pressure
          !
          p_nh_prog%exner(jc,jk,jb) = (p_nh_diag%pres(jc,jk,jb)/p0ref)**(rd/cpd)


          ! init virtual potential temperature
          !
          p_nh_prog%theta_v(jc,jk,jb) = p_nh_diag%temp(jc,jk,jb)/p_nh_prog%exner(jc,jk,jb)


          ! init density of moist air
          !
          p_nh_prog%rho(jc,jk,jb) = p_nh_diag%pres(jc,jk,jb) / (rd * t0)

        ENDDO  ! jc
      ENDDO  ! jk



      ! half level initialization
      !
      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx

          ! init pressure field (half levels)
          !
          p_nh_diag%pres_ifc(jc,jk,jb) = p0ref * exp(-p_metrics%z_ifc(jc,jk,jb)/scale_hgt)

          ! init density field (half levels)
          !
          p_nh_diag%rho_ic(jc,jk,jb)   = p_nh_diag%pres_ifc(jc,jk,jb) / (rd * t0)

        ENDDO  ! jc
      ENDDO  ! jk


      IF (ntracer_alloc > 1) THEN
        ! constant tracer field in order to check whether the given velocity field is non-divergent
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
           p_nh_prog%tracer(jc,jk,jb,2) = 1._wp
          ENDDO
        ENDDO
      ENDIF

    ENDDO ! jb
!$OMP ENDDO
!$OMP END PARALLEL



   ! set initial velocity field for t=0s
   CALL set_nh_velocity_hadley( p_patch, p_nh_prog, p_nh_diag, p_int, &
     &                          p_metrics, 0.0_wp )


  END SUBROUTINE init_nh_dcmip_hadley

!--------------------------------------------------------------------


!--------------------------------------------------------------------

  !>
  !! Initialization of horizontal and vertical velocity field for 
  !! the DCMIP Hadley-like meridional circulation test 
  !!
  SUBROUTINE set_nh_velocity_hadley( p_patch, p_nh_prog, p_nh_diag, p_int,  &
    &                                p_metrics, time, lacc )

    TYPE(t_patch),        INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_int_state),    INTENT(IN)    :: &  !< interpolation state vector
      &  p_int

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics

    REAL(wp), INTENT(IN) :: time              !< simulation time   [s]

    LOGICAL, INTENT(IN), OPTIONAL :: lacc     ! if true use OpenACC

    REAL(wp) ::  zu, zv                       !< zonal and meridional velocity

    REAL(wp) ::        &                      !< geometric height at edge points
      &  z_me(nproma,p_patch%nlev,p_patch%nblks_e), &
      &  z_ife(nproma,p_patch%nlevp1,p_patch%nblks_e)

    REAL(wp) ::        &                      !< density at edge points
      &  z_rho_e(nproma,p_patch%nlev,p_patch%nblks_e)

    !REAL(wp), DIMENSION(:,:), POINTER :: edges_center
    !REAL(wp), DIMENSION(:,:), POINTER :: edges_primal_normal

    REAL(wp) :: rho0                          !< surface density
    REAL(wp) :: z_lat                         !< geographical latitude
    REAL(wp) :: ztop                          !< model top
    INTEGER  :: jc, je, jk, jb                !< loop indices
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom  
    INTEGER  :: nlev, nlevp1                  !< number of full and half levels
    LOGICAL  :: lzacc

    !-----------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    !edges_center        => p_patch%edges%center
    !edges_primal_normal => p_patch%edges%primal_normal 


    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    !$ACC DATA CREATE(zu, zv, z_me, z_ife, z_rho_e) &
    !$ACC   PRESENT(p_patch, p_nh_prog, p_nh_diag, p_int, p_metrics) IF(lzacc)
    !!$ACC     PRESENT(edges_center, edges_primal_normal) IF(lzacc)

    ! Compute geometric height of full levels at edge midpoints
    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, p_int%c_lin_e, z_me, lacc=lzacc)

    ! Compute geometric height of half levels at edge midpoints
    CALL cells2edges_scalar(p_metrics%z_ifc, p_patch, p_int%c_lin_e, z_ife, lacc=lzacc)

    ! Compute rho at full level edge midpoints
    CALL cells2edges_scalar(p_nh_prog%rho, p_patch, p_int%c_lin_e, z_rho_e, lacc=lzacc)

    !$ACC WAIT

    ! syncs
    CALL sync_patch_array(SYNC_E, p_patch, z_me)
    CALL sync_patch_array(SYNC_E, p_patch, z_ife)
    CALL sync_patch_array(SYNC_E, p_patch, z_rho_e)


    ! density at surface
    rho0 = p0ref/(rd * t0)

    !
    ! set normal velocity field
    !

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,z_lat,zu,zv,ztop)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_lat, ztop, zu, zv)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          z_lat  = p_patch%edges%center(je,jb)%lat

          zu = u0*cos(z_lat)

          ztop   = z_ife(je,1,jb)

          zv = -(rho0/z_rho_e(je,jk,jb))                                 &
            &    * (grid_sphere_radius*w0*pi) /(nhadley*ztop)                &
            &    * cos(z_lat)*sin(nhadley*z_lat)*cos(pi*z_me(je,jk,jb)/ztop) &
            &    * cos(pi*time/tau)
          p_nh_prog%vn(je,jk,jb) = zu * p_patch%edges%primal_normal(je,jb)%v1  &
            &                    + zv * p_patch%edges%primal_normal(je,jb)%v2

        ENDDO  !je
      ENDDO  !jk
      !$ACC END PARALLEL

    ENDDO  !jb
!$OMP END DO


    !
    ! set vertical velocity field
    !

    i_rlstart = 1
    i_rlend   = min_rlcell

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,z_lat)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_lat)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          z_lat = p_patch%cells%center(jc,jb)%lat

          p_nh_prog%w(jc,jk,jb) = (w0/nhadley)                                    &
            &       * (rho0/p_nh_diag%rho_ic(jc,jk,jb))                           & 
            &       * (-2.0_wp*SIN(nhadley*z_lat)*SIN(z_lat)                      &
            &       + nhadley*COS(z_lat)*COS(nhadley*z_lat))                      &
            &       * SIN(pi*p_metrics%z_ifc(jc,jk,jb)/p_metrics%z_ifc(jc,1,jb))  &
            &       * COS(pi*time/tau)

        ENDDO  !jc
      ENDDO  !jk
      !$ACC END PARALLEL

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      p_nh_prog%w(i_startidx:i_endidx,nlevp1,jb) = 0._wp
      !$ACC END KERNELS

    ENDDO  !jb
    !$ACC WAIT(1)

    !$ACC END DATA

!$OMP ENDDO
!$OMP END PARALLEL


  END SUBROUTINE set_nh_velocity_hadley


END MODULE mo_nh_dcmip_hadley
