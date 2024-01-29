! Subroutine to initialize testcase 20  (examine the impact of 3D Schaer-like
! circular mountain profiles on an  atmosphere at rest ) proposed
! for the DCMIP summer school
!
! A non-rotating planet is used. Test 2-0 is conducted on an unscaled
! regular-size planet and primarily examines the accuracy of the pressure
! gradient calculation in a steady-state hydrostatically-balanced
! atmosphere at rest. This test is especially appealing for models with
! orography-following vertical coordinates. It increases the complexity
! of test 1-3, that investigated the impact of the same Schaer-type
! orographic profile on the accuracy of purely-horizontal passive
! tracer advection.
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

MODULE mo_nh_dcmip_rest_atm

   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd, grav
   USE mo_math_constants,       ONLY: pi
   USE mo_impl_constants,       ONLY: min_rlcell, min_rledge
   USE mo_loopindices,          ONLY: get_indices_c, get_indices_e
   USE mo_model_domain,         ONLY: t_patch
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_sync,                 ONLY: sync_patch_array, sync_patch_array_mult, &
     &                                SYNC_C
   USE mo_nh_init_utils,        ONLY: convert_thdvars  !, init_w
   USE mo_hydro_adjust,         ONLY: hydro_adjust 

   IMPLICIT NONE


   PRIVATE

   PUBLIC :: init_nh_prog_dcmip_rest_atm
   PUBLIC :: init_nh_topo_dcmip_rest_atm


!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of topograpphy for the nh schaer-type DCMIP test cases 
  !!
  SUBROUTINE init_nh_topo_dcmip_rest_atm( p_patch, topo_c, fis)

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

!  !DEFINED PARAMETERS for the 2-0 testcase (DCMIP):
   REAL(wp), PARAMETER :: h0 = 2000._wp ! maximum  mountain height(m)
   REAL(wp), PARAMETER :: lambdam = 3._wp*pi/2._wp ! Mountain longitud centerpoint
                                                   ! in radians
   REAL(wp), PARAMETER :: phim = 0._wp             ! Mountain latitud centerpoint
                                                   ! in radians
   REAL(wp), PARAMETER :: rm = 3._wp*pi/4._wp ! Mountain radius (radians)
   REAL(wp), PARAMETER :: zetam = pi/16._wp   ! Mountain oscillation half-width (radians)


!--------------------------------------------------------------------

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,i_rlstart,i_rlend)

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
       sin_tmp = SIN(z_lat) * SIN(phim)
       cos_tmp = COS(z_lat) * COS(phim)
       !great circle distance in radians
       r = ACOS (sin_tmp + cos_tmp *COS(z_lon-lambdam))
       IF (r .LT. rm) THEN
        topo_c(jc,jb) = (h0/2.0_wp) * (1._wp+COS(pi*r/rm)) * (COS(pi*r/zetam)**2)
       ELSE
        topo_c(jc,jb) = 0._wp
       END IF
       fis(jc,jb)    = grav * topo_c(jc,jb)

     ENDDO  !jc
   ENDDO  !jb
!$OMP END DO 

!$OMP END PARALLEL


   END SUBROUTINE init_nh_topo_dcmip_rest_atm
 
!-------------------------------------------------------------------------

  !>
  !! Initialization of prognostic state vector for the DCMIP nh 2-0 test case :
  !! steady state atmosphere at rest in the presence of orography.
  !!
  !!

  SUBROUTINE init_nh_prog_dcmip_rest_atm(p_patch, p_nh_prog, p_nh_diag, &
    &                                  p_metrics, l_hydro_adjust )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics
    LOGICAL,              INTENT(IN)    :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                          ! initial condition

    !local variables
    INTEGER  :: jc, je, jk, jb        !< loop indices
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom  
    INTEGER  :: nlev, nlevp1          !< number of full/half levels

!   !DEFINED PARAMETERS for the Schaer-type testcase (DCMIP):
    REAL(wp), PARAMETER :: p0 = 100000._wp     ! Reference surface pressure (Pa)
    REAL(wp), PARAMETER :: t0 = 300._wp        ! Temperature (K)
    REAL(wp), PARAMETER :: gamma   = 0.0065_wp ! temperature lapse rate (K/m)  

!--------------------------------------------------------------------


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    !
    ! Init prognostic variables vn, w
    !
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,i_rlstart,i_rlend)

    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          ! normal wind component
          p_nh_prog%vn(je,jk,jb) = 0._wp

        ENDDO  ! je
      ENDDO  ! jk
    ENDDO ! jb
!$OMP ENDDO




! initialized vertical velocity

  ! CALL sync_patch_array(SYNC_C, p_patch, p_nh_prog%w)

   i_rlstart = 1
   i_rlend   = min_rlcell

   i_startblk = p_patch%cells%start_blk(i_rlstart,1)
   i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,fac)
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
   DO jb = i_startblk, i_endblk

     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

     DO jk = 1, nlev
       DO jc = i_startidx, i_endidx

         ! init temp, pres
         !
         p_nh_diag%temp(jc,jk,jb) = t0 - gamma * p_metrics%z_mc(jc,jk,jb)
         p_nh_diag%pres(jc,jk,jb) = p0 *                                      &
           &      (1._wp -gamma/t0*p_metrics%z_mc(jc,jk,jb))**(grav/rd/gamma)

!         ! ALTERNATIVE implementation, using the vertical-mean temperature and 
!         ! pressure between two half-levels.
!         p_nh_diag%temp(jc,jk,jb) = t0 - 0.5_wp * gamma                           &
!           &    * (p_metrics%z_ifc(jc,jk,jb)**2 - p_metrics%z_ifc(jc,jk+1,jb)**2) &
!           &    / (p_metrics%z_ifc(jc,jk,jb)    - p_metrics%z_ifc(jc,jk+1,jb))

!         fac = 1._wp + grav/(rd*gamma)

!         p_nh_diag%pres(jc,jk,jb) = (p0*t0*rd)                                               &
!           &  /((rd*gamma + grav)*(p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb))) &
!           &  *((1._wp - (gamma*p_metrics%z_ifc(jc,jk+1,jb)/t0))**fac                        &
!           &  - (1._wp - (gamma*p_metrics%z_ifc(jc,jk,jb)  /t0))**fac )

         ! initialized vertical velocity
         p_nh_prog%w(jc,jk,jb) = 0._wp


       ENDDO !jc
     ENDDO !jk

     ! values at lower boundary
     p_nh_prog%w(i_startidx:i_endidx,nlevp1,jb)    = 0._wp

     DO jc = i_startidx, i_endidx
      p_nh_diag%pres_sfc(jc,jb) = p0 *                                      &
                 &  (1._wp -gamma/t0*p_metrics%z_ifc(jc,nlevp1,jb))**(grav/rd/gamma)
     ENDDO !jc

   ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

  CALL sync_patch_array_mult(SYNC_C, p_patch, 2,                                &
                           & p_nh_diag%temp, p_nh_diag%pres)



   ! As long as we do not have water vapour, ptr_nh_diag%temp is also the 
   ! virtual temperature

   CALL convert_thdvars(p_patch, p_nh_diag%pres, p_nh_diag%temp, &
                      & p_nh_prog%rho, p_nh_prog%exner, p_nh_prog%theta_v  )

   CALL sync_patch_array(SYNC_C, p_patch, p_nh_prog%w)
   CALL sync_patch_array_mult(SYNC_C, p_patch, 3,                               &
     &                        p_nh_prog%rho, p_nh_prog%exner, p_nh_prog%theta_v )


   IF (l_hydro_adjust) THEN

     CALL hydro_adjust ( p_patch, p_metrics, p_nh_prog%rho, &
                     & p_nh_prog%exner, p_nh_prog%theta_v   )

     CALL sync_patch_array_mult(SYNC_C, p_patch, 3, p_nh_prog%rho, p_nh_prog%exner, p_nh_prog%theta_v )
   END IF


  END SUBROUTINE init_nh_prog_dcmip_rest_atm
!--------------------------------------------------------------------
  END MODULE mo_nh_dcmip_rest_atm
