! Subroutine to initialized the Mountain Rossby Wave related test cases
! for the NH-Core (mrw_nh, mrw2_nh and mwbr_const)
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

MODULE mo_nh_mrw_exp

!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

   USE mo_kind,                ONLY: wp
   USE mo_model_domain,        ONLY: t_patch
   USE mo_parallel_config,     ONLY: nproma
   USE mo_math_constants,      ONLY: pi
   USE mo_physical_constants,  ONLY: rd, cpd, grav
   USE mo_extpar_config,       ONLY: itopo
   USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics

   USE mo_nh_init_utils,       ONLY: convert_thdvars, init_w
   USE mo_hydro_adjust,        ONLY: hydro_adjust
   USE mo_util_phys,           ONLY: virtual_temp
   USE mo_nh_jabw_exp,         ONLY: init_nh_inwp_tracers
   USE mo_loopindices,         ONLY: get_indices_e
   USE mo_run_config,          ONLY: iqv
   USE mo_intp_data_strc,      ONLY: t_int_state
   USE mo_exception,           ONLY: finish
   USE mo_impl_constants,      ONLY: inwp
   USE mo_sync,                ONLY: sync_patch_array, SYNC_C
   USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
   USE mo_grid_config,         ONLY: grid_sphere_radius, grid_angular_velocity

   IMPLICIT NONE


   PRIVATE


   REAL(wp), PUBLIC :: u0_mrw                 ! (m/s) wind speed for mrw and mwbr_const cases 
   REAL(wp), PUBLIC :: mount_height_mrw       ! (m) maximum mount height in mrw and mwbr
   REAL(wp), PUBLIC :: mount_half_width       ! (m) half width of mountain in mrw, mwbr and bell
   REAL(wp), PUBLIC :: mount_lonctr_mrw_deg   ! (deg) lon of mountain center in mrw and mwbr
   REAL(wp), PUBLIC :: mount_latctr_mrw_deg   ! (deg) lat of mountain center in mrw and mwbr
   REAL(wp), PUBLIC :: p_int_mwbr_const       ! pressure at interface in mwbr_const test case
   REAL(wp), PUBLIC :: temp_i_mwbr_const      ! temp in isothermal lower layer in mwbr_const
   REAL(wp), PUBLIC :: bruntvais_u_mwbr_const ! brunt vaisala freq in upper layer in mwbr_const

   PUBLIC :: init_nh_topo_mrw, init_nh_state_prog_mrw, init_nh_prog_mwbr_const

!  !DEFINED PARAMETERS for mountain induced Rossby wave train:
   REAL(wp), PARAMETER :: pres_sp  = 93000.0_wp  !pressure surface at the south pole
   REAL(wp), PARAMETER :: temp_mrw = 288._wp     !temperature of isothermal atmosphere


!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------

  !>
  !! Initialization of topography for the nh mrw test cases 
  !!
  SUBROUTINE init_nh_topo_mrw( ptr_patch, topo_c, nblks_c, npromz_c, l_modified )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

  
   INTEGER, INTENT (IN)     :: nblks_c, npromz_c
   REAL(wp), INTENT(INOUT)  :: topo_c    (nproma,nblks_c)
   LOGICAL , INTENT (IN)    :: l_modified

    ! local variables

   INTEGER        :: jc, jb, nlen
   REAL(wp)       :: z_lon, z_lat
   REAL(wp)       :: z_lon_ctr, z_lat_ctr
   REAL(wp)       :: zexp, zr
 

!--------------------------------------------------------------------


    z_lon_ctr = mount_lonctr_mrw_deg*pi/180.0_wp
    z_lat_ctr = mount_latctr_mrw_deg*pi/180.0_wp

 !$OMP PARALLEL
 !$OMP DO PRIVATE(jb,nlen,jc,z_lat,z_lon,zr,zexp)
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  npromz_c
        ENDIF
        DO jc = 1, nlen
          z_lat   = ptr_patch%cells%center(jc,jb)%lat
          z_lon   = ptr_patch%cells%center(jc,jb)%lon

          zr = SIN(z_lat_ctr)*SIN(z_lat)+COS(z_lat_ctr)*COS(z_lat)*COS(z_lon-z_lon_ctr) 
          zexp = grid_sphere_radius*ACOS(zr)/mount_half_width

          IF ( itopo==0 ) THEN
            IF (.NOT. l_modified ) THEN
                topo_c(jc,jb) =  mount_height_mrw*EXP( - zexp*zexp )
            ELSE
                topo_c(jc,jb) = &
                         mount_height_mrw*EXP( - zexp*zexp )*&
                         0.5_wp*(1._wp+COS(pi*zexp*2._wp))
            ENDIF
          ENDIF

        ENDDO
      ENDDO
!$OMP END DO 
!$OMP END PARALLEL

  END SUBROUTINE init_nh_topo_mrw
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state vector for the nh mrw test case 
  !! Tracers can also be initialized in case of inwp
  !!
  SUBROUTINE init_nh_state_prog_mrw( ptr_patch, ptr_nh_prog, ptr_nh_diag,      &
    &                                topo_c, p_metrics, p_int, l_hydro_adjust, &
    &                                iforcing, l_moist,  opt_rh_at_1000hpa,    &
    &                                opt_qv_max, opt_global_moist              ) 

   TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
     &  ptr_patch

   TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
     &  ptr_nh_prog

   TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
     &  ptr_nh_diag

   REAL(wp),  INTENT(IN)               :: topo_c(:,:)
   TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
   TYPE(t_int_state), INTENT(IN)       :: p_int
   LOGICAL, INTENT(IN)                 :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                         ! initial condition
   INTEGER, INTENT(IN)                 :: iforcing
   LOGICAL, INTENT(IN)                 :: l_moist !if .TRUE. tracers are initialized
   REAL(wp),INTENT(IN), OPTIONAL       :: opt_rh_at_1000hpa, opt_qv_max
   REAL(wp),INTENT(IN), OPTIONAL       :: opt_global_moist



   INTEGER                            :: je, jc, jk, jb, nlen, &
                                       & nblks_e, nblks_c, npromz_c
   INTEGER                            :: i_startidx, i_endidx, i_startblk
   INTEGER                            :: nlev        !< number of full levels
   REAL(wp)                           :: bruntvaissq, kappa, zhelp1, zhelp2, &
                                       & z_sfc, z_pres, z_klev, z_u, zcoslat, &
                                       & zlat, rh_at_1000hpa, qv_max
   REAL(wp), ALLOCATABLE              :: z_qv(:,:,:)
   LOGICAL        :: l_rediag

!-------------------------------------------------------------------------
   ! number of vertical levels
   nlev   = ptr_patch%nlev

   IF (l_moist) THEN
      ALLOCATE ( z_qv(nproma,nlev,ptr_patch%nblks_c) )
      IF (PRESENT(opt_rh_at_1000hpa)) THEN
        rh_at_1000hpa = opt_rh_at_1000hpa
      ELSE
        rh_at_1000hpa = 0.7_wp
      END IF
      IF (PRESENT(opt_qv_max)) THEN
       qv_max  = opt_qv_max
      ELSE
        qv_max = 20.e-3_wp
      END IF     
   END IF


   bruntvaissq = grav*grav/cpd/temp_mrw
   kappa       = rd/cpd
   zhelp1      = bruntvaissq/grav/grav/kappa
   zhelp2      = grav/rd/temp_mrw

   nblks_c   = ptr_patch%nblks_c
   npromz_c  = ptr_patch%npromz_c
   nblks_e   = ptr_patch%nblks_e


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,zlat,z_pres,z_sfc,z_klev,zcoslat)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
!     set the surface pressure. This is a diagnostic field
      DO jc = 1, nlen
        zlat= ptr_patch%cells%center(jc,jb)%lat
        zcoslat=COS(zlat)
        ptr_nh_diag%pres_sfc(jc,jb) = pres_sp * EXP( zhelp1 * ( u0_mrw *&
          & ( 0.5_wp*u0_mrw + grid_sphere_radius*grid_angular_velocity) * &
          & zcoslat*zcoslat - topo_c(jc,jb)*grav))
      ENDDO !jc
      DO jk = nlev, 1, -1
          ! Use analytic expressions at lowest model level
            DO jc = 1, nlen
              z_sfc  = topo_c(jc,jb)
              z_klev = 0.5_wp * ( p_metrics%z_ifc(jc,jk,jb) + &
                       p_metrics%z_ifc(jc,jk+1,jb) )
              z_pres = ptr_nh_diag%pres_sfc(jc,jb) * &
                       EXP(- zhelp2 * ( z_klev  - z_sfc )  )   !isothermal atm.

              ! initialized diagnostic fields
              ptr_nh_diag%temp(jc,jk,jb) = temp_mrw    
              ptr_nh_diag%pres(jc,jk,jb) = z_pres
            ENDDO !jc
      ENDDO !jk     
     ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

! As long as we do not have water vapour, ptr_nh_diag%temp is also the virtual temperature

  CALL convert_thdvars(ptr_patch, ptr_nh_diag%pres, ptr_nh_diag%temp, &
                     & ptr_nh_prog%rho, ptr_nh_prog%exner, ptr_nh_prog%theta_v  )


! initialized horizontal velocities
    i_startblk = ptr_patch%edges%start_blk(2,1)
    ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,zlat,z_u)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            zlat = ptr_patch%edges%center(je,jb)%lat
            z_u = u0_mrw * COS(zlat)  !v component is zero
            ptr_nh_prog%vn(je,jk,jb) = &
             z_u * ptr_patch%edges%primal_normal(je,jb)%v1
          ENDDO !je
        ENDDO !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

! initialized vertical velocity

   CALL init_w(ptr_patch, p_int, ptr_nh_prog%vn, p_metrics%z_ifc, ptr_nh_prog%w)
   CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)

! if physics, some fields like the pressure at the interface levels have to be initialized
  IF ( iforcing == inwp ) THEN

    CALL diagnose_pres_temp ( p_metrics, ptr_nh_prog, ptr_nh_prog, ptr_nh_diag,     &
                              ptr_patch, opt_calc_temp=.TRUE., opt_calc_pres=.TRUE. )

  END IF


! IF l_moist is .TRUE. the tracers are initialized similar as in jabw test case with moisture 
!  In this case the temp and pres fields should be kept and the virtual temperature and 
!  the NH prognostic variables have to be recalculated

  IF (l_moist) THEN

   l_rediag = .FALSE. 

   IF (PRESENT(opt_global_moist)) THEN
     CALL init_nh_inwp_tracers (ptr_patch, ptr_nh_prog, ptr_nh_diag, &
                              & p_metrics, rh_at_1000hpa, qv_max,    &
                              & l_rediag, opt_global_moist           )
   ELSE
     CALL init_nh_inwp_tracers (ptr_patch, ptr_nh_prog, ptr_nh_diag, &
                              & p_metrics, rh_at_1000hpa, qv_max,    &
                              & l_rediag                             )
 
   END IF
   !Calculate virtual temperature, pres field has not changed
   z_qv(:,:,:) = ptr_nh_prog%tracer(:,:,:,iqv)

   CALL virtual_temp ( ptr_patch, ptr_nh_diag%temp, z_qv,             &
                     & temp_v= ptr_nh_diag%tempv)
   
   !Calculate again the nh prognostic variables with the new tempv
   CALL convert_thdvars(ptr_patch, ptr_nh_diag%pres, ptr_nh_diag%tempv, &
               & ptr_nh_prog%rho, ptr_nh_prog%exner, ptr_nh_prog%theta_v  )


   DEALLOCATE(z_qv)

  END IF

  IF (l_hydro_adjust) THEN

   CALL hydro_adjust ( ptr_patch, p_metrics, ptr_nh_prog%rho,  &
                     & ptr_nh_prog%exner, ptr_nh_prog%theta_v  )

  END IF

   
  END SUBROUTINE init_nh_state_prog_mrw

!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state vector for the nh mrw test case 
  !! Tracers can also be initialized in case of inwp
  !!
  SUBROUTINE init_nh_prog_mwbr_const( ptr_patch, ptr_nh_prog, ptr_nh_diag,     &
    &                                topo_c, p_metrics, p_int, l_hydro_adjust, &
    &                                iforcing, l_moist,  opt_rh_at_1000hpa,    &
    &                                opt_qv_max, opt_global_moist              ) 

   TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
     &  ptr_patch

   TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
     &  ptr_nh_prog

   TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
     &  ptr_nh_diag

   REAL(wp),  INTENT(IN)               :: topo_c(:,:)
   TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
   TYPE(t_int_state), INTENT(IN)       :: p_int
   LOGICAL, INTENT(IN)                 :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                         ! initial condition
   INTEGER, INTENT(IN)                 :: iforcing
   LOGICAL, INTENT(IN)                 :: l_moist !if .TRUE. tracers are initialized
   REAL(wp),INTENT(IN), OPTIONAL       :: opt_rh_at_1000hpa, opt_qv_max
   REAL(wp),INTENT(IN), OPTIONAL       :: opt_global_moist



   INTEGER                            :: je, jc, jk, jb, nlen, &
                                       & nblks_e, nblks_c, npromz_c
   INTEGER                            :: i_startidx, i_endidx, i_startblk, icount
   INTEGER                            :: nlev        !< number of full levels
   REAL(wp)                           :: bruntvaissq_i,bruntvaissq_u, kappa, &
                                       & rkappa, zhelp3, zhelp4, &
                                       & z_sfc, z_pres, z_temp, z_klev, z_u, & 
                                       & zcoslat, zlat, zhelp1_i, zhelp2_i,  &
                                       & zhelp1_u, rh_at_1000hpa, qv_max
   REAL(wp), ALLOCATABLE              :: z_qv(:,:,:)
   LOGICAL        :: l_rediag


   REAL(wp), ALLOCATABLE              :: z_int_c(:,:)    ! z at interface in mwbr_const
   CHARACTER(len=*), PARAMETER :: &
    &  routine = '(mo_nh_mrw_exp) init_nh_prog_mwbr_const:'

!-------------------------------------------------------------------------

  bruntvaissq_i = grav*grav/cpd/temp_i_mwbr_const
  bruntvaissq_u = bruntvais_u_mwbr_const * bruntvais_u_mwbr_const
  kappa       = rd/cpd
  zhelp1_i     = bruntvaissq_i/grav/grav/kappa
  zhelp2_i     = grav/rd/temp_i_mwbr_const
  zhelp1_u     = bruntvaissq_u/grav/grav/kappa
  zhelp3       = u0_mrw/grid_sphere_radius + 2.0_wp*grid_angular_velocity
 ! theta_v_int  = temp_i_mwbr_const * (p0ref/p_int_mwbr_const)**kappa
  rkappa       = 1.0_wp/kappa

  nblks_c   = ptr_patch%nblks_c
  npromz_c  = ptr_patch%npromz_c
  nblks_e   = ptr_patch%nblks_e

  ! number of vertical levels
  nlev   = ptr_patch%nlev

  IF (l_moist) THEN
      ALLOCATE ( z_qv(nproma,nlev,ptr_patch%nblks_c) )
      IF (PRESENT(opt_rh_at_1000hpa)) THEN
        rh_at_1000hpa = opt_rh_at_1000hpa
      ELSE
        rh_at_1000hpa = 0.7_wp
      END IF
      IF (PRESENT(opt_qv_max)) THEN
       qv_max  = opt_qv_max
      ELSE
        qv_max = 20.e-3_wp
      END IF 
  END IF

  ALLOCATE (z_int_c(nproma,ptr_patch%nblks_c))

!$OMP PARALLEL PRIVATE(icount)
!$OMP DO PRIVATE(jb,nlen,jk,jc,zlat,z_pres,z_temp,z_sfc,z_klev,zcoslat,zhelp4)
  DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
!     set z_int_c, z at the interface
      icount = 0
      DO jc = 1, nlen
        zlat= ptr_patch%cells%center(jc,jb)%lat
        zcoslat=COS(zlat)

        z_sfc  = topo_c(jc,jb)
        z_int_c(jc,jb) = LOG(pres_sp/p_int_mwbr_const)/grav/zhelp1_i + &
                         zcoslat*zcoslat*zhelp3*grid_sphere_radius*u0_mrw/2.0_wp/grav
        IF (z_int_c(jc,jb) <  0._wp ) icount = icount + 1
        IF (z_int_c(jc,jb) >= z_sfc ) THEN

          ptr_nh_diag%pres_sfc(jc,jb) = pres_sp * EXP( zhelp1_i * ( u0_mrw * &
                              zhelp3*grid_sphere_radius * zcoslat*zcoslat/2.0_wp - &
                              z_sfc*grav))
        ELSE
          zhelp4 = z_sfc - z_int_c(jc,jb)
          ptr_nh_diag%pres_sfc(jc,jb) = p_int_mwbr_const * ( 1.0_wp +       &
                              (EXP(-bruntvaissq_u*zhelp4/grav) - 1.0_wp) *   &
                               bruntvaissq_i/bruntvaissq_u )**rkappa
        ENDIF
      ENDDO !jc
      DO jk = nlev, 1, -1
            DO jc = 1, nlen
              z_sfc  = topo_c(jc,jb)
              z_klev = 0.5_wp * ( p_metrics%z_ifc(jc,jk,jb) + &
                       p_metrics%z_ifc(jc,jk+1,jb) )

              IF (z_klev < z_int_c(jc,jb) ) THEN

                z_pres = ptr_nh_diag%pres_sfc(jc,jb) * &
                         EXP(- zhelp2_i * ( z_klev  - z_sfc )  )   !isothermal atm.

                ptr_nh_diag%temp(jc,jk,jb) = temp_i_mwbr_const    
                ptr_nh_diag%pres(jc,jk,jb) = z_pres
               ELSEIF (z_klev >= z_int_c(jc,jb) ) THEN
                zhelp4 = z_klev - z_int_c(jc,jb)
                z_pres = p_int_mwbr_const * ( 1.0_wp +       &
                              (EXP(-bruntvaissq_u*zhelp4/grav) - 1.0_wp) *   &
                               bruntvaissq_i/bruntvaissq_u )**rkappa
                z_temp = temp_i_mwbr_const * EXP (bruntvaissq_u*zhelp4/grav) * &
                             (z_pres/p_int_mwbr_const)**kappa
                ptr_nh_diag%temp(jc,jk,jb) =  z_temp  
                ptr_nh_diag%pres(jc,jk,jb) = z_pres 
              ENDIF
            ENDDO !jc
      ENDDO !jk     
  ENDDO !jb
!$OMP END DO
  IF (icount > 0) CALL finish(routine, &
    & 'z at interface is negative, needed p_int_mwbr_const < pres_sp')
!$OMP END PARALLEL

! As long as we do not have water vapour, ptr_nh_diag%temp is also the virtual temperature

  CALL convert_thdvars(ptr_patch, ptr_nh_diag%pres, ptr_nh_diag%temp, &
                     & ptr_nh_prog%rho, ptr_nh_prog%exner, ptr_nh_prog%theta_v  )

    i_startblk = ptr_patch%edges%start_blk(2,1)
    ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,zlat,z_u)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            zlat = ptr_patch%edges%center(je,jb)%lat
            z_u = u0_mrw * COS(zlat)  !v component is zero
            ptr_nh_prog%vn(je,jk,jb) = &
             z_u * ptr_patch%edges%primal_normal(je,jb)%v1
          ENDDO
        ENDDO
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
! initialized vertical velocity

   CALL init_w(ptr_patch, p_int, ptr_nh_prog%vn, p_metrics%z_ifc, ptr_nh_prog%w)
   CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)

! if physics, some fields like the pressure at the interface levels and others have to be initialized
  IF ( iforcing == inwp ) THEN

    CALL diagnose_pres_temp ( p_metrics, ptr_nh_prog, ptr_nh_prog, ptr_nh_diag,     &
                              ptr_patch, opt_calc_temp=.TRUE., opt_calc_pres=.TRUE. )

  END IF

! IF l_moist is .TRUE. the tracers are initialized similar as in jabw test case with moisture 
!  In this case the temp and pres fields should be kept and the virtual temperature and 
!  the NH prognostic variables have to be recalculated

  IF (l_moist) THEN

   l_rediag = .FALSE. 

   IF (PRESENT(opt_global_moist)) THEN
     CALL init_nh_inwp_tracers (ptr_patch, ptr_nh_prog, ptr_nh_diag, &
                              & p_metrics, rh_at_1000hpa, qv_max,    &
                              & l_rediag, opt_global_moist           )
   ELSE
     CALL init_nh_inwp_tracers (ptr_patch, ptr_nh_prog, ptr_nh_diag, &
                              & p_metrics, rh_at_1000hpa, qv_max,    &
                              & l_rediag                             )
 
   END IF
   !Calculate virtual temperature, pres field has not changed
   z_qv(:,:,:) = ptr_nh_prog%tracer(:,:,:,iqv)

   CALL virtual_temp ( ptr_patch, ptr_nh_diag%temp, z_qv,             &
                     & temp_v= ptr_nh_diag%tempv)
   
   !Calculate again the nh prognostic variables with the new tempv
   CALL convert_thdvars(ptr_patch, ptr_nh_diag%pres, ptr_nh_diag%tempv, &
               & ptr_nh_prog%rho, ptr_nh_prog%exner, ptr_nh_prog%theta_v  )
   DEALLOCATE(z_qv)

  END IF


  IF (l_hydro_adjust) THEN

   CALL hydro_adjust ( ptr_patch, p_metrics, ptr_nh_prog%rho, &
                     & ptr_nh_prog%exner, ptr_nh_prog%theta_v )

  END IF

  DEALLOCATE (z_int_c)

  END SUBROUTINE init_nh_prog_mwbr_const
!--------------------------------------------------------------------
  END MODULE mo_nh_mrw_exp
