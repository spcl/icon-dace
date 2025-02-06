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

! Pure 3D-Advection test case adapted to the NH-Core
!
! This module contains parameters, initialization subroutines and
! functions to be used in the Pure 3D-Advection test case of the
! non-hydrostatic dynamical core.
!
! Literature
! - Jablonowski et al. (2008): Idealized test cases for the dynamical cores
!   of Atmospheric General Circulation Models: A proposal for the NCAR ASP
!   2008 summer colloquium

MODULE mo_nh_pa_test

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
USE mo_impl_constants,      ONLY: SUCCESS
USE mo_physical_constants,  ONLY: rd, cpd, p0ref
USE mo_math_constants,      ONLY: pi_2, pi
USE mo_model_domain,        ONLY: t_patch
USE mo_ext_data_types,      ONLY: t_external_data
USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: nsteps, ntracer
USE mo_ncar_testcases,      ONLY: init_pure_adv_wind, init_pure_adv_tracers, &
  & init_ncar_testcases_domain
USE mo_io_units,            ONLY: find_next_free_unit
USE mo_exception,           ONLY: finish
USE mo_mpi,                 ONLY: my_process_is_stdio
USE mo_time_base,           ONLY: rdaylen
USE mo_grid_config,         ONLY: grid_angular_velocity

IMPLICIT NONE

PRIVATE

REAL(wp), PARAMETER :: zp0      = 100000._wp          !< surface pressure
REAL(wp), PARAMETER :: zt0      = 300._wp             !< atmospheric temperature
REAL(wp), PARAMETER :: zscale_h = 8781.419026_wp      !< scale height
REAL(wp), PARAMETER :: ztau     = 4._wp * rdaylen     !< 4 days period
REAL(wp), PARAMETER :: zomega0  = pi*40000._wp / ztau !< maximum vertical pressure velocity

PUBLIC :: init_nh_state_prog_patest, set_nh_w_rho

!--------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!
  !>
  !!               Initialization of prognostic state vector.
  !!
  !!
  SUBROUTINE init_nh_state_prog_patest( ptr_patch, ptr_int, ptr_nh_prog,      &
    &                                   ptr_nh_diag, ptr_ext_data, p_metrics, &
    &                                   p_rotate_axis_deg, linit_tracer_fv,   &
    &                                   tracer_inidist_list )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_int_state), INTENT(IN)       :: ptr_int

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag

    TYPE(t_external_data), INTENT(INOUT):: &  !< external data
      &  ptr_ext_data

    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state

    REAL(wp), INTENT(IN)              :: &  !< rotation angle
      &  p_rotate_axis_deg

    LOGICAL, INTENT(IN)  :: linit_tracer_fv  !< tracer finite volume initialization

    INTEGER, INTENT(IN) :: tracer_inidist_list(:) !< selected initial tracer distributions

    INTEGER  :: nblks_e, nblks_c, nblks_v, npromz_e, npromz_c, npromz_v, &
                nlen, it4, it5, it6, it7, it8
    INTEGER  :: nlev                    !< number of full levels
    INTEGER  :: je, jk, jc, jt, jb, jv  ! loop variables

    REAL(wp) :: zlon, zlat, zheight     !< location
    REAL(wp) :: zu, zv, zq4, zq5, zq6, zq7, zq8   !< initialized variables
    REAL(wp) :: z_aleph, rovcp


    INTEGER :: ilc1, ibc1  !< line and block indices of cell1 adjacent 
                           !< to the current edge
!--------------------------------------------------------------------
!
    CALL init_ncar_testcases_domain()

    rovcp   = rd/cpd                            !< kappa
    z_aleph = p_rotate_axis_deg * pi/180.0_wp   !< deg2rad rotation angle


    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e
    npromz_e = ptr_patch%npromz_e
    nblks_v  = ptr_patch%nblks_v
    npromz_v = ptr_patch%npromz_v

    nlev = ptr_patch%nlev

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

      DO jk = 1, nlev

        ! compute full level pressure
        ptr_nh_diag%pres(1:nlen,jk,jb) = zp0                                      &
          &                  * exp(-p_metrics%z_mc(1:nlen,jk,jb)/zscale_h)

        ! init virtual potential temperature
        ptr_nh_prog%theta_v(1:nlen,jk,jb) = zt0                                   &
          & * (ptr_nh_diag%pres_sfc(1:nlen,jb)/ptr_nh_diag%pres(1:nlen,jk,jb))**rovcp

        ! init density field rho
        ptr_nh_prog%rho(1:nlen,jk,jb) = ptr_nh_diag%pres(1:nlen,jk,jb)            &
          &                           / (rd * zt0)

        ! init exner pressure
        ptr_nh_prog%exner(1:nlen,jk,jb) = (ptr_nh_diag%pres(1:nlen,jk,jb)/p0ref)**rovcp

      ENDDO !jk
    ENDDO !jb
!$OMP END DO


    !
    ! initialize horizontal velocity field (time independent)
    !
!$OMP DO PRIVATE(jb,jk,je,nlen,zlon,zlat,zu,zv,ilc1,ibc1)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
         nlen = nproma
      ELSE
         nlen = npromz_e
      ENDIF
      DO jk = 1, nlev

        DO je = 1, nlen

          ! location of edge midpoint
          zlon = ptr_patch%edges%center(je,jb)%lon
          zlat = ptr_patch%edges%center(je,jb)%lat


          CALL init_pure_adv_wind ( zlon, zlat,                      &
                                    p_rotate_axis_deg, zu, zv )


          ! calculate normal wind component
          ptr_nh_prog%vn(je,jk,jb) =                                 &
                      zu * ptr_patch%edges%primal_normal(je,jb)%v1   &
                    + zv * ptr_patch%edges%primal_normal(je,jb)%v2

          ! Coriolis parameter
          ptr_patch%edges%f_e(je,jb) = 2.0_wp*grid_angular_velocity &
            & *(SIN(zlat)*COS(z_aleph) - COS(zlon)*COS(zlat)*SIN(z_aleph))

          ! get line and block indices of 
          ilc1 = ptr_patch%edges%cell_idx(je,jb,1)
          ibc1 = ptr_patch%edges%cell_blk(je,jb,1)

          ! Since rho is constant in horizontal direction, there is no need 
          ! for interpolation.
          ptr_nh_diag%mass_fl_e(je,jk,jb) = ptr_nh_prog%rho(ilc1,jk,ibc1)  &
            * ptr_nh_prog%vn(je,jk,jb) * p_metrics%ddqz_z_full_e(je,jk,jb)

        ENDDO  ! edge loop
      ENDDO  ! level loop
    ENDDO  ! block loop
!$OMP END DO
!$OMP END PARALLEL


    !
    ! initialize vertical velocity field (time level n)
    !
    CALL set_nh_w_rho( ptr_patch, p_metrics, 0, 0._wp, 0._wp, ptr_nh_prog%w, &
      &                ptr_nh_diag%pres, ptr_nh_diag%rho_ic )



    !
    ! initialize tracer fields
    !
    it4 = 0
    it5 = 0
    it6 = 0
    it7 = 0
    it8 = 0
    DO jt = 1, ntracer
      SELECT CASE(tracer_inidist_list(jt))
      CASE(4)
        it4 = jt
      CASE(5)
        it5 = jt
      CASE(6)
        it6 = jt
      CASE(7)
        it7 = jt
      CASE(8)
        it8 = jt
      END SELECT
    ENDDO
    IF ( it4==0 .AND. it5==0 .AND. it6==0 .AND. it7==0 .AND. it8==0 ) THEN
      CALL finish('mo_nh_pa_test:init_nh_state_prog_patest',  &
         &         'incorrect tracer list')
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen,zlon,zlat,zheight,zq4,zq5,zq6,zq7,zq8)
    DO jb = 1, nblks_c
      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)
      DO jk = 1, nlev


        DO jc = 1, nlen

          ! location of quadrature points
          IF (linit_tracer_fv) THEN
            zlon = ptr_int%gquad%qpts_tri_l(jc,jb)%lon
            zlat = ptr_int%gquad%qpts_tri_l(jc,jb)%lat
          ELSE
            zlon = ptr_patch%cells%center(jc,jb)%lon
            zlat = ptr_patch%cells%center(jc,jb)%lat
          ENDIF

          zheight = p_metrics%z_mc(jc,jk,jb)


          CALL init_pure_adv_tracers ( tracer_inidist_list(1:5), zlon, zlat, zheight, &
            &                          zq4, zq5, zq6, zq7, zq8 )

          IF (it4 /= 0) ptr_nh_prog%tracer(jc,jk,jb,it4) = zq4
          IF (it5 /= 0) ptr_nh_prog%tracer(jc,jk,jb,it5) = zq5
          IF (it6 /= 0) ptr_nh_prog%tracer(jc,jk,jb,it6) = zq6
          IF (it7 /= 0) ptr_nh_prog%tracer(jc,jk,jb,it7) = zq7
          IF (it8 /= 0) ptr_nh_prog%tracer(jc,jk,jb,it8) = zq8


          zlon = ptr_patch%cells%center(jc,jb)%lon
          zlat = ptr_patch%cells%center(jc,jb)%lat

          ! Coriolis parameter
          ptr_patch%cells%f_c(jc,jb) = 2.0_wp*grid_angular_velocity &
            & *(SIN(zlat)*COS(z_aleph)-COS(zlon)*COS(zlat)*SIN(z_aleph))
        ENDDO ! cell loop

      ENDDO ! vertical level loop

    ENDDO ! block loop
!$OMP ENDDO


!$OMP DO PRIVATE(jb,jv,nlen,zlat,zlon)
    DO jb = 1, nblks_v
      IF (jb /= nblks_v) THEN
         nlen = nproma
      ELSE
         nlen = npromz_v
      ENDIF
      DO jv = 1, nlen
        zlat   = ptr_patch%verts%vertex(jv,jb)%lat
        zlon   = ptr_patch%verts%vertex(jv,jb)%lon
        ! Coriolis parameter
        ptr_patch%verts%f_v(jv,jb) = 2.0_wp*grid_angular_velocity &
          & *(SIN(zlat)*COS(z_aleph)-COS(zlon)*COS(zlat)*SIN(z_aleph))
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_nh_state_prog_patest

!-------------------------------------------------------------------------
!
  !>
  !!  Set the time-variant vertical velocity for pure 3D advection test
  !!
  !!  Set the time-variant vertical velocity for pure 3D advection test
  !!  (non-hydrostatic version). Optionally, an updated density
  !!  field may be computed (commented out, so far)
  !!
  SUBROUTINE set_nh_w_rho( ptr_patch, p_metrics, k_step, p_dtime, p_sim_time, &
    &                      p_w_prog, p_diag_pres, p_diag_rho_ic )


    TYPE(t_patch), INTENT(IN)       :: ptr_patch
    TYPE(t_nh_metrics), INTENT(IN)  :: p_metrics !< NH metrics state

    INTEGER,  INTENT(IN) :: k_step      !< actual time step
    REAL(wp), INTENT(IN) :: p_sim_time  !< simulation time in seconds since start
    REAL(wp), INTENT(IN) :: p_dtime     !< time step
    REAL(wp), INTENT(INOUT) ::   &      !< prognostic vertical velocity in height
      &  p_w_prog(:,:,:)                !< coordinates
                                        !< dim:(nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(INOUT) ::   &      !< diagnostic full level pressure
      &  p_diag_pres(:,:,:)

    REAL(wp), INTENT(INOUT) ::   &      !< diagnostic half level density
      &  p_diag_rho_ic(:,:,:)

    INTEGER  :: jb, jk, jc              !< loop indices
    INTEGER  :: nlen, nblks_c, npromz_c, ist, ikcenter
    INTEGER  :: nlev, nlevp1            !< number of full and half levels

    REAL(wp) :: &                       !< pressure at interface levels
      &  z_pres_ic(nproma,ptr_patch%nlevp1,ptr_patch%nblks_c)

!!$    REAL(wp) :: &                       !< updated pressure
!!$      &  zpres_help(nproma,ptr_patch%nlev,ptr_patch%nblks_c)

    REAL(wp) ::  &                      !< vertical velocity in pressure coordinates
     &  zomega(nproma,ptr_patch%nlevp1,ptr_patch%nblks_c)
    REAL(wp) :: zomega_t                !< time dependent part of the vertical
                                        !< velocity in pressure coordinates
    REAL(wp) :: zshape                  !< shape function
    REAL(wp) :: zomega_mid              !< vertical velocity in pressure coordinates
                                        !< at vertical mid-point of model domain
    REAL(wp) :: zw                      !< vertical velocity in height coordinates

    ! for writing data
    INTEGER, SAVE :: n_file_ti          ! file identifier
    CHARACTER(LEN=*), PARAMETER :: file_ti = 'vertical_velocity.dat'
    REAL(wp) :: zsim_days

!--------------------------------------------------------------------
!

!     write(0,*) "set_nh_w_rho: k_step=", k_step
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev     = ptr_patch%nlev
    nlevp1   = ptr_patch%nlevp1

    zomega   = 0.0_wp
    ikcenter = MAX(2,nlev/2)

    ! time dependent part of vertical velocity field
    zomega_t = zomega0 * COS( 2._wp*pi/ztau * (p_sim_time+p_dtime) )


    ! Open the datafile in first time step
    IF ( my_process_is_stdio() .AND. k_step == 1 ) THEN
      n_file_ti = find_next_free_unit(10,20)
      OPEN( UNIT=n_file_ti, FILE=TRIM(file_ti), FORM='FORMATTED', IOSTAT=ist )
      IF (ist/=SUCCESS) THEN
         CALL finish('mo_nh_pa_test:set_nh_w_rho',&
                     'could not open datafile')
      ENDIF
      WRITE (n_file_ti,'(4A22,A8,I3,A4,F12.5)') &
        &  ' TIMESTEP            ,'    ,       &
        &  ' ELAPSED TIME  (days),'    ,       &
        &  ' grid_angular_velocity         (Pa/s),'    ,       &
        &  ' W             (m/s) '     ,       &
        &  ' at Z(', ikcenter, ') = '  ,       &
        &   p_metrics%z_mc(1,ikcenter,1)
    ENDIF


    ! re-diagnose the pressure at full levels (equals initialization)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen)
    DO jb = 1, nblks_c
      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)

      DO jk = 1, nlev

        DO jc=1, nlen

          p_diag_pres(jc,jk,jb) = zp0                                          &
            &                * exp(-p_metrics%z_mc(jc,jk,jb)/zscale_h)

        ENDDO

      ENDDO
    ENDDO
!$OMP END DO


    ! compute pressure and density at interface levels (time level n)
!$OMP DO PRIVATE(jb,jk,jc,nlen)
    DO jb = 1, nblks_c
      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)

      ! lower boundary
      z_pres_ic(1:nlen,nlevp1,jb) = zp0
      ! upper boundary
      z_pres_ic(1:nlen,1,jb) =  zp0                                             &
          &               * exp(-p_metrics%z_ifc(1:nlen,1,jb)/zscale_h)

      ! lower boundary
      p_diag_rho_ic(1:nlen,nlevp1,jb) = z_pres_ic(1:nlen,nlevp1,jb)/(rd*zt0)
      ! upper boundary
      p_diag_rho_ic(1:nlen,1,jb) =  z_pres_ic(1:nlen,1,jb)/(rd*zt0)

      DO jk = 2, nlev

        DO jc= 1, nlen

         z_pres_ic(jc,jk,jb) = p_metrics%wgtfac_c(jc,jk,jb)         &
           &                 * p_diag_pres(jc,jk,jb)                &
           &                 + (1._wp-p_metrics%wgtfac_c(jc,jk,jb)) &
           &                 * p_diag_pres(jc,jk-1,jb)

         ! density at half levels
         p_diag_rho_ic(jc,jk,jb) = z_pres_ic(jc,jk,jb)/(rd*zt0)

        ENDDO

      ENDDO
    ENDDO
!$OMP END DO


    !
    ! compute vertical velocity field (time level n+1)
    !

    ! add time independent part to time dependent part zomega_t
!$OMP DO PRIVATE(jb,jk,jc,nlen,zshape) LASTPRIVATE(zomega_mid,zw)
    DO jb = 1, nblks_c
      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)

      ! top and bottom boundary values
      p_w_prog(:,1,jb)        = 0._wp
      p_w_prog(:,nlevp1,jb)   = 0._wp

      DO jk = 2, nlev

        DO jc=1, nlen

          ! shape function
          zshape = MIN( 1._wp,2._wp*SQRT( SIN( (z_pres_ic(jc,jk,jb)-z_pres_ic(jc,1,jb))&
            &    / (zp0-z_pres_ic(jc,1,jb))*pi) ) )

          ! compute vertical velocity
          zomega(jc,jk,jb)  = zomega_t * SIN( zshape * pi_2 )

          ! get vertical velocity (height based)
          p_w_prog(jc,jk,jb) = - zscale_h/z_pres_ic(jc,jk,jb) * zomega(jc,jk,jb)

        ENDDO ! cell loop

        jc = MAX(1,nlen/2)

        IF ( jk == ikcenter ) THEN
          ! vertical velocity in pressure coordinates
          zomega_mid     = p_metrics%wgtfac_c(jc,jk,jb) * zomega(jc,jk,jb) &
            &            + (1._wp-p_metrics%wgtfac_c(jc,jk,jb))            &
            &            * zomega(jc,jk-1,jb)

          ! vertical velocity in height coordinates
          zw             = p_metrics%wgtfac_c(jc,jk,jb) * p_w_prog(jc,jk,jb) &
            &            + (1._wp-p_metrics%wgtfac_c(jc,jk,jb))              &
            &            * p_w_prog(jc,jk-1,jb)
        ENDIF

      ENDDO ! vertical level loop

    ENDDO ! block loop
!$OMP END DO
!$OMP END PARALLEL

!DR Test in order to check correctness of CFL-independent PPM-scheme
!DR    p_w_prog(:,2,:) = 0.6 * p_w_prog(:,2,:)


! It is suggested by Jablonowski to update the density, which is needed in order
! to convert between \rho*q and q. So far this has been neglected.

!!! ATTENTION: To avoid unused dummy arguments, p_rho_prog_new and 
!!! p_rho_prog_now have been removed from the argument list. These must be 
!!! re-included in case the following code is used.
!    !
!    ! compute updated pressure and density at full levels
!    ! (time level n+1)
!    !
!    DO jb = 1, nblks_c
!      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)

!      DO jk = 1, nlev

!        DO jc=1, nlen

!          zpres_help(jc,jk,jb) = p_diag_pres(jc,jk,jb) + p_dtime               &
!            &   * zomega0 * COS( 2._wp*pi/ztau * (p_sim_time+0.5_wp*p_dtime) ) &
!            &   * SIN(zshape * pi_2)

!          p_rho_prog_new(jc,jk,jb) = zpres_help(jc,jk,jb)/(rd * zt0)

!          ! reset density at time n
!          p_rho_prog_now(jc,jk,jb) = p_diag_pres(jc,jk,jb)/(rd * zt0)
!        ENDDO ! cell loop

!      ENDDO ! vertical level loop

!    ENDDO

    ! write vertical velocity to data file
    IF(my_process_is_stdio() .AND. k_step >= 1) THEN
      zsim_days = p_sim_time / rdaylen
      WRITE(n_file_ti,'(i21,f22.8,e22.8,e22.8)') &
            k_step, zsim_days, zomega_mid, zw

      IF ( k_step == nsteps ) THEN
        CLOSE(n_file_ti)
      ENDIF
    ENDIF


  END SUBROUTINE set_nh_w_rho

!--------------------------------------------------------------------
  END MODULE mo_nh_pa_test
