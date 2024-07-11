! Offline integration of density for idealized test cases
!
! Integrates density for idealized test cases, i.e. when the dynamical core
! is switched off.
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

MODULE mo_integrate_density_pa

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge_int, min_rlcell_int, &
    &                               MIURA, MIURA3, ippm_v
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_advection_config,    ONLY: advection_config 
  USE mo_advection_traj,      ONLY: t_back_traj, btraj_compute_o1
  USE mo_advection_hflux,     ONLY: upwind_hflux_miura, upwind_hflux_miura3
  USE mo_advection_vflux,     ONLY: upwind_vflux_ppm
  USE mo_math_divrot,         ONLY: div
  USE mo_intp,                ONLY: cells2edges_scalar
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_edge
  USE mo_sync,                ONLY: sync_patch_array, SYNC_E, SYNC_C
  USE mo_exception,           ONLY: finish


  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_integrate_density_pa'



  PUBLIC :: integrate_density_pa


CONTAINS

  !---------------------------------------------------------------------------
  !>
  !! SUBROUTINE integrate_density_pa
  !!
  !! Short description:
  !! Integrates mass equation when running in pure advection mode
  !!
  SUBROUTINE integrate_density_pa( p_patch, p_int, p_prog_now, p_prog_new, &
    &                             p_metrics, p_diag, p_dtime, k_step,      &
    &                             lcoupled_rho )

    TYPE(t_patch),      INTENT(INOUT):: p_patch
    TYPE(t_int_state),  INTENT(IN)   :: p_int
    TYPE(t_nh_prog),    INTENT(IN)   :: p_prog_now
    TYPE(t_nh_prog),    INTENT(INOUT):: p_prog_new
    TYPE(t_nh_metrics), INTENT(IN)   :: p_metrics
    TYPE(t_nh_diag),    INTENT(INOUT):: p_diag

    REAL(wp),           INTENT(IN)   :: p_dtime

    INTEGER,            INTENT(IN)   :: k_step       !< time step counter [1]
    LOGICAL,            INTENT(IN)   :: lcoupled_rho !< integrate mass equation (TRUE/FALSE)
 
    REAL(wp) ::   &                   !< density edge value
      &  z_rho_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) ::   &                   !< time averaged normal velocities
      &  z_vn_traj(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) ::   &                   !< time averaged tangential velocities
      &  z_vt_traj(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) ::   &                   !< time averaged vertical velocities
      &  z_w_traj(nproma,p_patch%nlevp1,p_patch%nblks_c)
    REAL(wp) ::  &                    !< flux divergence at cell center
      &  z_fluxdiv_rho(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) ::   &                   !< vertical mass flux
      &  z_mflx_contra_v(nproma,p_patch%nlevp1)

    INTEGER  :: nlev, nlevp1          !< number of full/half levels
    INTEGER  :: i_startblk, i_startidx, i_endblk, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

    INTEGER  :: jb, je, jc, jk        !< loop indices
    INTEGER  :: ikp1                  !< jk +1

    TYPE(t_back_traj) :: btraj
    LOGICAL  :: lcompute, lcleanup

    INTEGER  :: pid                   !< patch ID

    REAL(wp), POINTER ::  &
      & ptr_current_rho(:,:,:) => NULL()  !< pointer to density field

    CHARACTER(len=*), PARAMETER :: routine = modname//':integrate_density_pa'

    !---------------------------------------------------------------------------

!#ifdef _OPENACC
!    CALL finish (routine,': OpenACC version currently not implemented')
!#endif

    ! get patch ID
    pid = p_patch%id


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1


    i_nchdom = MAX(1,p_patch%n_childdom)

    !$ACC DATA CREATE(z_rho_e, z_vn_traj, z_vt_traj, z_w_traj, z_fluxdiv_rho, z_mflx_contra_v, btraj, ptr_current_rho) &
    !$ACC   PRESENT(p_patch, p_int, p_prog_now, p_prog_new, p_metrics, p_diag)

    IF (lcoupled_rho) THEN  ! integrate mass equation

#ifdef _OPENACC
    CALL finish (routine,': OpenACC version in case of lcoupled_rho currently not implemented')
#endif
            
      lcompute =.TRUE.
      lcleanup =.TRUE.


      ! point to correct starting time level for upcoming integration
      !
      ptr_current_rho => p_prog_now%rho

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

      i_rlstart = 1
      i_rlend   = min_rledge_int

      i_startblk = p_patch%edges%start_blk(i_rlstart,1)
      i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

      ! preparations for horizontal and vertical integration
      !
!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        ! time-averaged velocity field for mass flux and trajectory 
        ! computation.
        DO jk = 1,nlev
          DO je = i_startidx, i_endidx

            z_vn_traj(je,jk,jb) = 0.5_wp * (p_prog_now%vn(je,jk,jb) &
              &                           + p_prog_new%vn(je,jk,jb)) 

          ENDDO  ! je
        ENDDO  ! jk
      ENDDO  ! jb
!$OMP END DO NOWAIT


      IF ( advection_config(pid)%lvadv_tracer ) THEN

        i_rlstart = 1
        i_rlend   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(i_rlstart,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jk = 1,nlevp1
            DO jc = i_startidx, i_endidx

              z_w_traj(jc,jk,jb) = 0.5_wp * (p_prog_now%w(jc,jk,jb) &
                &                         +  p_prog_new%w(jc,jk,jb)) 

            ENDDO  ! jc
          ENDDO  ! jk
        ENDDO  ! jb
!$OMP END DO

      ENDIF
!$OMP END PARALLEL



      !***************************************************
      ! Vertical integration of mass continuity equation
      !***************************************************
      IF ( advection_config(pid)%lvadv_tracer .AND. MOD( k_step, 2 ) == 0 ) THEN

        i_rlstart = grf_bdywidth_c-1
        i_rlend   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(i_rlstart,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


!DR replaced 
!DR \rho^{n+1/2}*w^{n+1/2} by w^{n+1/2}
!DR \rho^{n} \Delta z by \Delta z

        ! CALL third order PPM (unrestricted timestep-version) (i.e. CFL>1)
        CALL upwind_vflux_ppm( p_patch, ptr_current_rho,                   &! in
          &                  advection_config(pid)%iubc_adv,               &! in
          &                  z_w_traj, p_dtime,                            &! in
          &                  lcompute, lcleanup,                           &! in
          &                  advection_config(pid)%itype_vlimit(1),        &! in
          &                  advection_config(pid)%ivlimit_selective(1),   &! in
          &                  p_metrics%ddqz_z_full,                        &! in
          &                  p_metrics%ddqz_z_full,                        &! in
          &                  .TRUE.,                                       &! in
          &                  ippm_v,                                       &! in           
          &                  p_diag%rho_ic,                                &! out
          &                  opt_lout_edge = .TRUE.,                       &! in
          &                  opt_rlstart=i_rlstart,                        &! in
          &                  opt_rlend=i_rlend                             )! in




!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_mflx_contra_v,ikp1)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend )


          !
          ! compute vertical mass flux
          !
          DO jk = 1, nlevp1

            DO jc = i_startidx, i_endidx

              z_mflx_contra_v(jc,jk) = z_w_traj(jc,jk,jb) * p_diag%rho_ic(jc,jk,jb)

            ENDDO  ! jc

          ENDDO  ! jk


          !
          ! compute vertical flux divergence
          !
          DO jk = 1, nlev

            ! index of top half level
            ikp1 = jk + 1

            DO jc = i_startidx, i_endidx

              p_prog_new%rho(jc,jk,jb) = ptr_current_rho(jc,jk,jb) - p_dtime  &
                &                      * ( z_mflx_contra_v(jc,  jk)           &
                &                      -   z_mflx_contra_v(jc,ikp1)  )        &
                &                      * p_metrics%inv_ddqz_z_full(jc,jk,jb)
            ENDDO  ! jc
          ENDDO  ! jk
        ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

        ! point to correct starting time level for upcoming integration
        !
        ptr_current_rho => p_prog_new%rho

      ENDIF  ! lvadv_tracer


      !***************************************************
      ! Horizontal integration of mass continuity equation
      !***************************************************

      i_rlstart = grf_bdywidth_e-1
      i_rlend   = min_rledge_int-1

      i_startblk = p_patch%edges%start_blk(i_rlstart,1)
      i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)


      ! reconstruct tangential velocity component at edge midpoints
      CALL rbf_vec_interpol_edge( z_vn_traj, p_patch, p_int,    &! in
        &                         z_vt_traj, opt_rlend=i_rlend  )! inout
      !
      ! get new Density 'edge-value'
      !

      ! Select desired flux calculation method
      !
      ! Note: It is implicitly assumed, that only one tracer is advected.
      !       The namelist settings of the first tracer are applied to 
      !       the density as well.
      SELECT CASE( advection_config(pid)%ihadv_tracer(1) )
      CASE( MIURA )

        ! 1st order backward trajectory
        CALL btraj_compute_o1 ( btraj    = btraj,            & !inout
          &                  ptr_p       = p_patch,          & !in
          &                  ptr_int     = p_int,            & !in
          &                  p_vn        = z_vn_traj,        & !in
          &                  p_vt        = z_vt_traj,        & !in
          &                  p_dthalf    = 0.5_wp * p_dtime, & !in
          &                  opt_rlstart = 5,                & !in
          &                  opt_rlend   = i_rlend,          & !in
          &                  opt_slev    = 1,                & !in
          &                  opt_elev    = p_patch%nlev      ) !in

        ! CALL MIURA with second order accurate reconstruction
        CALL upwind_hflux_miura( p_patch, ptr_current_rho, p_prog_new%vn,  &! in
          &                 p_dtime, p_int, btraj,                         &! in
          &                 advection_config(pid)%igrad_c_miura,           &! in
          &                 advection_config(pid)%itype_hlimit(1),         &! in
          &                 z_rho_e,                                       &! inout
          &                 opt_lout_edge=.TRUE.,                          &! in
          &                 opt_rlstart_e= 5,                              &! in
          &                 opt_rlend_e  = i_rlend,                        &! in
          &                 opt_slev     = 1,                              &! in
          &                 opt_elev     = p_patch%nlev                    )! in

      CASE( MIURA3 )

        CALL upwind_hflux_miura3(p_patch, ptr_current_rho, p_prog_new%vn,&
          &                      z_vn_traj, z_vt_traj, p_dtime, p_int,   &
          &                      lcompute, lcleanup,                     &
          &                      advection_config(pid)%itype_hlimit(1),  &
          &                      z_rho_e, opt_lout_edge=.TRUE.,          &
          &                      opt_rlend=i_rlend  )
      END SELECT




      !
      ! massflux at edges
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jk = 1,nlev
          DO je = i_startidx, i_endidx
            p_diag%mass_fl_e(je,jk,jb) = z_vn_traj(je,jk,jb)                & 
              &                         * p_metrics%ddqz_z_full_e(je,jk,jb) &
              &                         * z_rho_e(je,jk,jb)   
          ENDDO  ! je
        ENDDO  ! jk
      ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


      !
      ! Compute updated density field
      !
      CALL div( p_diag%mass_fl_e(:,:,:), p_patch, &! in
        &       p_int, z_fluxdiv_rho(:,:,:)       )! in,inout


      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell

      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        ! New density field
        DO jk = 1,nlev
          DO jc = i_startidx, i_endidx
            p_prog_new%rho(jc,jk,jb)= ptr_current_rho(jc,jk,jb)         &
              &                    - p_dtime * z_fluxdiv_rho(jc,jk,jb)  &
              &                    * p_metrics%inv_ddqz_z_full(jc,jk,jb)

          ENDDO  ! jc
        ENDDO  ! jk
      ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


      ! point to correct starting time level for upcoming integration
      !
      ptr_current_rho => p_prog_new%rho


      !**************************************************
      ! Vertical integration of mass continuity equation
      !**************************************************

      IF ( advection_config(pid)%lvadv_tracer .AND. MOD( k_step, 2 ) == 1 ) THEN

      ! Integrate mass equation also in the vertical

        i_rlstart = grf_bdywidth_c+1
        i_rlend   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(i_rlstart,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)



!DR replaced 
!DR \rho^{n+1/2}*w^{n+1/2} by w^{n+1/2}
!DR \rho^{n} \Delta z by \Delta z

        ! CALL third order PPM (unrestricted timestep-version) (i.e. CFL>1)
        CALL upwind_vflux_ppm( p_patch, ptr_current_rho,                   &! in
          &                  advection_config(pid)%iubc_adv,               &! in
          &                  z_w_traj, p_dtime,                            &! in
          &                  lcompute, lcleanup,                           &! in
          &                  advection_config(pid)%itype_vlimit(1),        &! in
          &                  advection_config(pid)%ivlimit_selective(1),   &! in
          &                  p_metrics%ddqz_z_full,                        &! in
          &                  p_metrics%ddqz_z_full,                        &! in
          &                  .TRUE.,                                       &! in
          &                  ippm_v,                                       &! in           
          &                  p_diag%rho_ic,                                &! out
          &                  opt_lout_edge = .TRUE.,                       &! in
          &                  opt_rlstart=i_rlstart,                        &! in
          &                  opt_rlend=i_rlend                             )! in



        !
        ! compute vertical mass flux
        !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_mflx_contra_v,ikp1)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend )

          DO jk = 1, nlevp1

            DO jc = i_startidx, i_endidx

              z_mflx_contra_v(jc,jk) = z_w_traj(jc,jk,jb) * p_diag%rho_ic(jc,jk,jb)

            ENDDO  ! jc

          ENDDO  ! jk



          DO jk = 1, nlev

            ! index of top half level
            ikp1 = jk + 1

            DO jc = i_startidx, i_endidx

              p_prog_new%rho(jc,jk,jb) =  ptr_current_rho(jc,jk,jb) - p_dtime  &
                &                      * ( z_mflx_contra_v(jc,  jk)            &
                &                      -   z_mflx_contra_v(jc,ikp1)     )      &
                &                      * p_metrics%inv_ddqz_z_full(jc,jk,jb)
            ENDDO  ! jc
          ENDDO  ! jk
        ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

      ENDIF  ! lvadv_tracer


    ELSE   ! IF not lcoupled_rho

      ! mass continuity equation will not be re-integrated.
      ! Only the updated mass fluxes are provided.
      ! This version is only recommended, if:
      ! - the density field is constant
      ! - the velocity field is divergence free


      i_rlstart = grf_bdywidth_c-1
      i_rlend   = min_rledge_int

      i_startblk = p_patch%edges%start_blk(i_rlstart,1)
      i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)


      ! Compute density at cell edges
      CALL cells2edges_scalar(p_prog_now%rho, p_patch, p_int%c_lin_e, z_rho_e, lacc=.TRUE.)

      !$ACC WAIT

      CALL sync_patch_array(SYNC_E, p_patch, z_rho_e)


      !
      ! massflux at edges
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1,nlev
          DO je = i_startidx, i_endidx

            p_diag%mass_fl_e(je,jk,jb) = z_rho_e(je,jk,jb)                 &
                                       * 0.5_wp * (p_prog_now%vn(je,jk,jb) &
              &                        + p_prog_new%vn(je,jk,jb))          &
              &                        * p_metrics%ddqz_z_full_e(je,jk,jb)

          ENDDO  ! je
        ENDDO  ! jk
        !$ACC END PARALLEL

      ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    ENDIF  ! lcoupled_rho


    CALL sync_patch_array(SYNC_E, p_patch, p_diag%mass_fl_e)
    CALL sync_patch_array(SYNC_C, p_patch, p_diag%rho_ic )
!$ACC WAIT(1)
  !$ACC END DATA

  END SUBROUTINE integrate_density_pa



END MODULE mo_integrate_density_pa

