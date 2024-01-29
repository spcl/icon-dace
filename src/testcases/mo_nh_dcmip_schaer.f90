! Subroutine to initialize testcases 21 22 (Mountain waves over a Schaer
! type mountain on a small planet without and with wind shear) proposed
! for the DCMIP summer school
!
! The tests in this section examine the response of atmospheric models to flow over a
! mountain profile, with and without wind shear. In order to ensure the simulated response
! contains both hydrostatic and non-hydrostatic features, the radius of the Earth is scaled so that the
! simulation is in the non-hydrostatic domain. We chose a non-rotating Earth with angular velocity
! = 0 s^-1 and select a reduced-size Earth with radius as = a/X. The reduction factor is set to
! X = 500. This choice leads to an Earth with a circumference at the equator of about 2pi a/X, which is
! approximately 80 km. These underlying ideas behind the tests are based on the work of Schaer et al. (MWR 2002),
! Wedi and Smolarkiewicz (QJR 2009), and Wedi et al. (ECMWF Tech Report 2009)
! Note however that in the presence of vertical wind shear we do not recommend the isothermal conditions
! suggested in the literature. They lead to imbalances of the initial conditions in spherical geometry
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

MODULE mo_nh_dcmip_schaer

   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd, grav
   USE mo_math_constants,       ONLY: pi
   USE mo_impl_constants,       ONLY: min_rlcell, min_rledge
   USE mo_parallel_config,      ONLY: nproma
   USE mo_loopindices,          ONLY: get_indices_c, get_indices_e
   USE mo_model_domain,         ONLY: t_patch
   USE mo_grid_config,          ONLY: grid_sphere_radius
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_ref, t_nh_metrics
   USE mo_intp_data_strc,       ONLY: t_int_state
   USE mo_intp,                 ONLY: cells2edges_scalar
   USE mo_sync,                 ONLY: sync_patch_array, sync_patch_array_mult, &
     &                                SYNC_C, SYNC_E
   USE mo_nh_init_utils,        ONLY: convert_thdvars  !, init_w 
   USE mo_hydro_adjust,         ONLY: hydro_adjust

   IMPLICIT NONE


   PRIVATE

   PUBLIC :: init_nh_prog_dcmip_schaer
   PUBLIC :: init_nh_topo_dcmip_schaer

   LOGICAL, PUBLIC :: lshear_dcmip

!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of topograpphy for the nh schaer-type DCMIP test cases 
  !!
  SUBROUTINE init_nh_topo_dcmip_schaer( p_patch, topo_c, fis)

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

  
   REAL(wp), INTENT(INOUT)  :: topo_c    (:,:)
   REAL(wp), INTENT(INOUT)  :: fis       (:,:)

   ! local variables
   INTEGER     :: jc, jb
   REAL(wp)    :: r, z_lat, z_lon
   REAL(wp)    :: sin_tmp, cos_tmp
   INTEGER     :: i_startidx, i_endidx, i_startblk, i_endblk
   INTEGER     :: i_rlstart, i_rlend, i_nchdom 

!  !DEFINED PARAMETERS for the Schaer-type testcase (DCMIP):
   REAL(wp), PARAMETER :: h0 = 250._wp  ! maximum schaer-type mountain height(m)
   REAL(wp), PARAMETER :: d  = 5000._wp ! Schaer-type mountain half-width (m)
   REAL(wp), PARAMETER :: xi = 4000._wp ! Schaer-type mountain wavelenght (m)
   REAL(wp), PARAMETER :: lambdac = pi/4._wp ! longitude of mountain centerpoint
   REAL(wp), PARAMETER :: phic = 0._wp ! Latitude of mountain centerpoint

!--------------------------------------------------------------------

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk, i_rlstart, i_rlend)

    i_rlstart = 1
    i_rlend   = min_rlcell

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,z_lat,z_lon,sin_tmp,cos_tmp,r)
   DO jb = i_startblk, i_endblk

     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

     DO jc = i_startidx, i_endidx

       z_lat =  p_patch%cells%center(jc,jb)%lat
       z_lon =  p_patch%cells%center(jc,jb)%lon
       sin_tmp = SIN(z_lat) * SIN(phic)
       cos_tmp = COS(z_lat) * COS(phic)
       !great circle distance
       r = grid_sphere_radius * ACOS (sin_tmp + cos_tmp *COS(z_lon-lambdac))
       topo_c(jc,jb) = h0 * EXP(- (r**2)/(d**2)) * (COS(pi*r/xi)**2)
       fis(jc,jb)    = grav * topo_c(jc,jb)

     ENDDO  !jc
   ENDDO  !jb
!$OMP END DO 
!$OMP END PARALLEL


   END SUBROUTINE init_nh_topo_dcmip_schaer
 
!-------------------------------------------------------------------------

  !>
  !! Initialization of prognostic state vector for the DCMIP nh  schaer-type 
  !! test case.
  !!
  !!
  !!  The tests in this section examine the response of atmospheric models to flow over a 
  !! mountain profile, with and without wind shear. In order to ensure the simulated response
  !! contains both hydrostatic and non-hydrostatic features, the radius of the Earth is scaled so that the
  !! simulation is in the non-hydrostatic domain. We chose a non-rotating Earth with angular velocity
  !! = 0 s^-1 and select a reduced-size Earth with radius as = a/X. The reduction factor is set to
  !! X = 500. This choice leads to an Earth with a circumference at the equator of about 2pi a/X, which is 
  !! approximately 80 km. These underlying ideas behind the tests are based on the work of Schaer et al. (MWR 2002),
  !! Wedi and Smolarkiewicz (QJR 2009), and Wedi et al. (ECMWF Tech Report 2009)
  !! Note however that in the presence of vertical wind shear we do not recommend the isothermal conditions 
  !! suggested in the literature. They lead to imbalances of the initial conditions in spherical geometry

  SUBROUTINE init_nh_prog_dcmip_schaer(p_patch, p_nh_prog, p_nh_diag, p_nh_ref, &
    &                                  p_metrics, p_int, l_hydro_adjust )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_nh_ref),       INTENT(INOUT) :: &  !< reference state vector
      &  p_nh_ref

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics
    TYPE(t_int_state),    INTENT(IN)    :: &  !< interpolation state
      &  p_int
    LOGICAL,              INTENT(IN)    :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                          ! initial condition

    !local variables
    INTEGER  :: jc, je, jk, jb        !< loop indices
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom  
    INTEGER  :: nlev, nlevp1          !< number of full/half levels
    REAL(wp) :: z_lat                 !< geographical coordinates
    REAL(wp) :: temp_e, zu, zv, c
    REAL(wp) :: &                     !< height at model level edges
      &  z_me(nproma,p_patch%nlev,p_patch%nblks_e)

!   !DEFINED PARAMETERS for the Schaer-type testcase (DCMIP):
    REAL(wp), PARAMETER :: peq = 100000._wp  ! Reference surface pressure at the equator (Pa)
    REAL(wp), PARAMETER :: teq = 300._wp     ! Reference surface temperature at the equator (K)
    REAL(wp), PARAMETER :: ueq = 20._wp      ! Reference zonal wind velocity (m/s)
    REAL(wp), PARAMETER :: cs  = 2.5e-4_wp   ! equatorial surface wind shear 
                                             ! (lshear_dcmip .TRUE.) (m-1)

!--------------------------------------------------------------------


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    ! Compute geometric height at edge points
    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, &
             p_int%c_lin_e, z_me)

    CALL sync_patch_array(SYNC_E,p_patch,z_me)


    ! Mountain waves with or without vertical shear
    IF (lshear_dcmip) THEN
      c= cs
    ELSE
      c= 0._wp
    END IF


    !
    ! Init prognostic variables vn, w
    !
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,i_rlstart,i_rlend)

    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,z_lat,zu,zv,temp_e)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          ! get geographical coordinates of edge midpoint
          !
          z_lat = p_patch%edges%center(je,jb)%lat
          ! the temperature at the edge is needed
          temp_e = teq *(1.0_wp - (c*ueq*ueq/(grav))*(SIN(z_lat)**2) )

          zu = ueq * COS(z_lat) * SQRT( (2.0_wp*teq/temp_e)*c*z_me(je,jk,jb) + temp_e/teq) 
          ! meridional velocity
          zv = 0._wp

          ! compute normal wind component
          p_nh_prog%vn(je,jk,jb) = zu * p_patch%edges%primal_normal(je,jb)%v1  &
            &                    + zv * p_patch%edges%primal_normal(je,jb)%v2


          ! copy vn to reference state vector (needed by Rayleigh damping mechanism)
          p_nh_ref%vn_ref(je,jk,jb) = p_nh_prog%vn(je,jk,jb)

        ENDDO  ! je
      ENDDO  ! jk
    ENDDO ! jb
!$OMP ENDDO




! initialized vertical velocity

  ! CALL init_w(p_patch, p_int, p_nh_prog%vn, p_metrics%z_ifc, p_nh_prog%w)
  ! CALL sync_patch_array(SYNC_C, p_patch, p_nh_prog%w)

   i_rlstart = 1
   i_rlend   = min_rlcell

   i_startblk = p_patch%cells%start_blk(i_rlstart,1)
   i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,z_lat)
   DO jb = i_startblk, i_endblk

     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

     DO jk = 1, nlev
       DO jc = i_startidx, i_endidx

         z_lat = p_patch%cells%center(jc,jb)%lat

         ! init temperature, does not depend on height
         !
         p_nh_diag%temp(jc,jk,jb) = teq*(1.0_wp - (c*ueq*ueq/(grav))*(SIN(z_lat)**2) )

         p_nh_diag%pres(jc,jk,jb) = peq*EXP( -(ueq*ueq/(2.0_wp*rd*teq))*(SIN(z_lat)**2) & 
                      & - grav*p_metrics%z_mc(jc,jk,jb)/(rd*p_nh_diag%temp(jc,jk,jb))    )

         ! initialized vertical velocity
         p_nh_prog%w(jc,jk,jb) = 0._wp

         ! copy w to reference state vector (needed by Rayleigh damping mechanism)
         p_nh_ref%w_ref(jc,jk,jb) = p_nh_prog%w(jc,jk,jb)

       ENDDO !jc
     ENDDO !jk

     ! values at lower boundary
     p_nh_prog%w(i_startidx:i_endidx,nlevp1,jb)    = 0._wp
     p_nh_ref%w_ref(i_startidx:i_endidx,nlevp1,jb) = p_nh_prog%w(i_startidx:i_endidx,nlevp1,jb)


     DO jc = i_startidx, i_endidx

       z_lat = p_patch%cells%center(jc,jb)%lat

       p_nh_diag%pres_sfc(jc,jb) = peq*EXP( -(ueq*ueq/(2.0_wp*rd*teq))*(SIN(z_lat)**2) & 
                &  - grav* p_metrics%z_ifc(jc,nlevp1,jb)/(rd*p_nh_diag%temp(jc,nlev,jb)) )
     ENDDO !jc

   ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL




! As long as we do not have water vapour, ptr_nh_diag%temp is also the virtual temperature

   CALL convert_thdvars(p_patch, p_nh_diag%pres, p_nh_diag%temp, &
                      & p_nh_prog%rho, p_nh_prog%exner, p_nh_prog%theta_v  )


   CALL sync_patch_array_mult(SYNC_C, p_patch, 2, p_nh_prog%w, p_nh_ref%w_ref)
   CALL sync_patch_array_mult(SYNC_C, p_patch, 3,                               &
     &                        p_nh_prog%rho, p_nh_prog%exner, p_nh_prog%theta_v )


   IF (l_hydro_adjust) THEN

     CALL hydro_adjust ( p_patch, p_metrics, p_nh_prog%rho, &
                       & p_nh_prog%exner, p_nh_prog%theta_v )

     CALL sync_patch_array_mult(SYNC_C, p_patch, 3, p_nh_prog%rho, p_nh_prog%exner, p_nh_prog%theta_v)
   END IF


  END SUBROUTINE init_nh_prog_dcmip_schaer
!--------------------------------------------------------------------
  END MODULE mo_nh_dcmip_schaer
