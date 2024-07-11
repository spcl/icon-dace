!
! Calculates 3D subgrid-scale (with metric terms) for viscosity and diffusivity
! in the nonhydrostatic model
!
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
#include "omp_definitions.inc"
!----------------------------


MODULE mo_sgs_turbmetric

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_intp,                ONLY: cells2verts_scalar, cells2edges_scalar
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: iqv, iqc, iqtke, msg_level
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge_int, min_rlcell_int, min_rlvert_int,   &
                                    iprog
  USE mo_math_utilities,      ONLY: tdma_solver, tdma_solver_vec
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
                                    sync_patch_array_mult
  USE mo_physical_constants,  ONLY: cpd, rcvd, rcpd, alv, grav
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag
  USE mo_surface_les,         ONLY: surface_conditions
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_les_config,          ONLY: les_config
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_turbulent_diagnostic,ONLY: is_sampling_time, idx_sgs_th_flx, &
                                    idx_sgs_qv_flx, idx_sgs_qc_flx,   &
                                    idx_sgs_u_flx, idx_sgs_v_flx
  USE mo_statistics,          ONLY: levels_horizontal_mean
  USE mo_nh_vert_interp_les,  ONLY: brunt_vaisala_freq, vert_intp_full2half_cell_3d
  USE mo_fortran_tools,       ONLY: copy, init, assert_acc_device_only
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config 


  IMPLICIT NONE

  PRIVATE
  REAL(wp),         PARAMETER :: z_1by3  = 1._wp/3._wp, z_2by3  = 2._wp/3._wp, &
                                 z_4by3  = 4._wp/3._wp

  !Parameter for vertical scheme type
  INTEGER, PARAMETER :: iexplicit = 1
  INTEGER, PARAMETER :: iimplicit = 2

  ! parameter to distinguish different tracers
  INTEGER, PARAMETER :: tracer_theta = 1
  INTEGER, PARAMETER :: tracer_qv = 2
  INTEGER, PARAMETER :: tracer_qc = 3

  PUBLIC :: drive_subgrid_diffusion_m

  !Variables for the module
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: km_iv, km_ie, km_c, kh_ie, theta_v,    &
                                             rho_ic, div_c, u_vert, v_vert,         &
                                             w_ie, w_vert, u_iv, v_iv, vn_ie, vt_ie

  !Store variables for metric-terms implementation

  CHARACTER(len=*), PARAMETER :: inmodule = 'mo_sgs_turbmetric:'

  CONTAINS

  !-------------------------------------------------------------------------------------
  !>
  !! drive_subgrid_diffusion_m
  !!------------------------------------------------------------------------
  !! Driver for computing the sgs viscosity and diffusivity using Smagorisnky model
  !! and evaluating diffusion term on triangles
  !!------------------------------------------------------------------------
  SUBROUTINE drive_subgrid_diffusion_m(p_sim_time, p_nh_prog, p_nh_prog_now_rcf, p_nh_prog_rcf,   &
                                       p_nh_diag, p_nh_metrics, p_patch, p_int, p_prog_lnd_now,   &
                                       p_prog_lnd_new, p_diag_lnd, prm_diag, prm_nwp_tend, dt, lacc)

    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog     !< single nh prognostic state
    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog_now_rcf !< old state for tke
    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog_rcf     !< rcf nh prognostic state    
    TYPE(t_nh_diag),   INTENT(inout)     :: p_nh_diag     !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    TYPE(t_patch),  INTENT(inout),TARGET :: p_patch       !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int         !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(in)        :: p_prog_lnd_now!<land prog state
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_new!<land prog state
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars
    REAL(wp),          INTENT(in)        :: dt
    REAL(wp),          INTENT(in)        :: p_sim_time    !current sim time

    REAL(wp), ALLOCATABLE :: theta(:,:,:)
    REAL(wp), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_e) :: &
         D_11_ie, D_12_ie, D_13_ie

    INTEGER :: nlev, nlevp1
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jb, jc, jk, jg

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc


    CALL assert_acc_device_only("drive_subgrid_diffusion_m", lacc)

    IF (msg_level >= 15) &
         CALL message(inmodule, 'drive_subgrid_diffusion_m')

    jg = p_patch%id

    nlev   = p_patch%nlev
    nlevp1 = nlev+1

    ALLOCATE( u_vert(nproma,nlev,p_patch%nblks_v),   &
              v_vert(nproma,nlev,p_patch%nblks_v),   &
              u_iv(nproma,nlevp1,p_patch%nblks_v),   &
              v_iv(nproma,nlevp1,p_patch%nblks_v),   &
              w_vert(nproma,nlevp1,p_patch%nblks_v), &
              vn_ie(nproma,nlevp1,p_patch%nblks_e),  &
              vt_ie(nproma,nlevp1,p_patch%nblks_e),  &
              w_ie(nproma,nlevp1,p_patch%nblks_e),   &
              km_iv(nproma,nlevp1,p_patch%nblks_v),  &
              km_c(nproma,nlev,p_patch%nblks_c),     &
              km_ie(nproma,nlevp1,p_patch%nblks_e),  &
              theta(nproma,nlev,p_patch%nblks_c),    &
              theta_v(nproma,nlev,p_patch%nblks_c),  &
              div_c(nproma,nlev,p_patch%nblks_c),    &
              rho_ic(nproma,nlevp1,p_patch%nblks_c)  &
            )

    !$ACC DATA CREATE(u_vert, v_vert, w_vert, u_iv, v_iv, vn_ie, vt_ie, w_ie) &
    !$ACC   CREATE(km_iv, km_c, km_ie, theta_v, div_c, rho_ic) &
    !$ACC   CREATE(theta, D_11_ie, D_12_ie, D_13_ie)

    !Initialize

!$OMP PARALLEL
    CALL init(km_iv(:,:,:), lacc=.TRUE.)
    CALL init(km_c(:,:,:), lacc=.TRUE.)
    CALL init(km_ie(:,:,:), lacc=.TRUE.)

    IF(p_test_run)THEN
     CALL init(u_vert(:,:,:), lacc=.TRUE.)
     CALL init(v_vert(:,:,:), lacc=.TRUE.)
     CALL init(w_vert(:,:,:), lacc=.TRUE.)
     CALL init(u_iv(:,:,:), lacc=.TRUE.)
     CALL init(v_iv(:,:,:), lacc=.TRUE.)
     CALL init(vn_ie(:,:,:), lacc=.TRUE.)
     CALL init(vt_ie(:,:,:), lacc=.TRUE.)
    END IF
!$OMP END PARALLEL

    !Convert temperature to potential temperature: all routines within
    !use theta.

    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,nlev
        DO jc = i_startidx, i_endidx
          theta(jc,jk,jb)   = p_nh_diag%temp(jc,jk,jb)  / p_nh_prog%exner(jc,jk,jb)

          theta_v(jc,jk,jb) = p_nh_diag%tempv(jc,jk,jb) / p_nh_prog%exner(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !Get rho at interfaces to be used later
    CALL vert_intp_full2half_cell_3d(p_patch, p_nh_metrics, p_nh_prog%rho, rho_ic,    &
                                     2, min_rlcell_int-2, lacc=.TRUE.)

    CALL surface_conditions(p_nh_metrics, p_patch, p_nh_diag, p_prog_lnd_now,         &
                            p_prog_lnd_new, p_diag_lnd, prm_diag, theta,              &
                            p_nh_prog%tracer(:,:,:,iqv), p_sim_time, lacc=.TRUE.)

    IF ( atm_phy_nwp_config(jg)%inwp_turb == iprog) THEN 
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'inwp_turb=iprog: OpenACC version currently not implemented')
#endif
      CALL prognostic_tke(p_nh_prog, p_nh_prog_now_rcf, p_nh_prog_rcf, p_nh_diag,                 &
                          p_nh_metrics, p_patch, p_int, prm_diag, dt, D_11_ie, D_12_ie, D_13_ie)
    ELSE
      CALL brunt_vaisala_freq(p_patch, p_nh_metrics, nproma, theta_v, prm_diag%bruvais, lacc=.TRUE.)

      CALL smagorinsky_model(p_nh_prog, p_nh_metrics, p_patch, p_int, prm_diag%tkvh,              &
                             prm_diag%mech_prod, prm_diag%tkvm, prm_diag%bruvais,                 &
                             D_11_ie, D_12_ie, D_13_ie)
    END IF

    CALL diffuse_hori_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, prm_diag,      &
                               prm_nwp_tend%ddt_u_turb, prm_nwp_tend%ddt_v_turb,                  &
                               D_11_ie, D_12_ie, D_13_ie, dt)

    !Vertical velocity is updated here
    CALL diffuse_vert_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int,                &
                               prm_diag%tkvm, prm_nwp_tend%ddt_w_turb, dt)

    CALL diffuse_scalar(theta, p_nh_metrics, p_patch, p_int, prm_nwp_tend%ddt_temp_turb,          &
                        p_nh_prog%exner, prm_diag, p_nh_prog%rho, dt, tracer_theta)
                        

    !For qv and qc: implement for qr as well
    IF(.NOT.les_config(jg)%is_dry_cbl)THEN
      CALL diffuse_scalar(p_nh_prog%tracer(:,:,:,iqv), p_nh_metrics, p_patch, p_int,  &
                          prm_nwp_tend%ddt_tracer_turb(:,:,:,iqv), p_nh_prog%exner,   &
                          prm_diag, p_nh_prog%rho, dt, tracer_qv)

      CALL diffuse_scalar(p_nh_prog%tracer(:,:,:,iqc), p_nh_metrics, p_patch, p_int,  &
                          prm_nwp_tend%ddt_tracer_turb(:,:,:,iqc), p_nh_prog%exner,   &
                          prm_diag, p_nh_prog%rho, dt, tracer_qc)
    ELSE
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'only moist cbl is suppported on GPU')
#endif

!$OMP PARALLEL
      CALL init(prm_nwp_tend%ddt_tracer_turb(:,:,:,iqv), lacc=.TRUE.)
      CALL init(prm_nwp_tend%ddt_tracer_turb(:,:,:,iqc), lacc=.TRUE.)
!$OMP END PARALLEL
    END IF

    !$ACC WAIT
    !$ACC END DATA ! CREATE module & subroutine
    DEALLOCATE(u_vert, v_vert, w_vert, vn_ie, vt_ie, w_ie, km_iv, km_ie,  &
               theta, theta_v, km_c, rho_ic, div_c, u_iv, v_iv)

  END SUBROUTINE drive_subgrid_diffusion_m
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! smagorinsky_model
  !!------------------------------------------------------------------------
  !! Computes the sgs viscosity and diffusivity using Smagorisnky model
  !! \tau_ij = KD_ij where D_ij = du_i/dx_j+du_j/dx_i
  !!
  !! and  K = cs*\Delta**2*D/sqrt(2) where D = sqrt(D_ijD_ij)
  !!
  !! and D**2 = D_11**2 + D_22**2 + D_33**2 + 2D_12**2 + 2D_13**2 + 2D_23**2
  !!
  !! where, D_11 = 2du_1/dx_1
  !!        D_22 = 2d_u2/dx_2
  !!        D_33 = 2d_u3/dx_3
  !!        D_12 = du_1/dx_2 + du_2/dx_1
  !!        D_13 = du_1/dx_3 + du_3/dx_1
  !!        D_23 = du_2/dx_3 + du_3/dx_2
  !! For triangles: 1=normal, 2=tangential, and 3 = z directions
  !!------------------------------------------------------------------------
  SUBROUTINE smagorinsky_model(p_nh_prog, p_nh_metrics, p_patch, p_int, &
                               kh_ic, mech_prod, km_ic, bruvais, &
                               D_11_ie, D_12_ie, D_13_ie)

    TYPE(t_patch),  INTENT(inout),TARGET :: p_patch    !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    REAL(wp), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_c), INTENT(inout) :: &
         kh_ic, mech_prod, km_ic, bruvais      !< atm phys vars
    REAL(wp), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_e), INTENT(out) :: &
         D_11_ie, D_12_ie, D_13_ie

    ! local variables
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end, jg
    INTEGER :: jk, jb, jc, je
    INTEGER  :: nlev, nlevp1             !< number of full levels
    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk

    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4 
    REAL(wp) :: vt_vert1, vt_vert2, vt_vert3, vt_vert4
    REAL(wp) :: w_full_c1, w_full_c2, w_full_v1, w_full_v2
    REAL(wp) :: D_11, D_12, D_13, D_22, D_23, D_33

    REAL(wp) :: shear(nproma,p_patch%nlev,p_patch%nblks_e)           
    REAL(wp) :: div_of_stress(nproma,p_patch%nlev,p_patch%nblks_e)    

    !--------------------------------------------------------------------------

    !$ACC DATA CREATE(shear, div_of_stress)

    IF (msg_level >= 18) &
         CALL message(inmodule, 'smagorinsky_model')

    jg = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = nlev+1

    !Initialize
    IF(p_test_run)THEN
!$OMP PARALLEL
      CALL init(kh_ic, lacc=.TRUE.)
      CALL init(shear, lacc=.TRUE.)
      CALL init(div_of_stress, lacc=.TRUE.)
      CALL init(div_c, lacc=.TRUE.)
!$OMP END PARALLEL
    END IF

    !--------------------------------------------------------------------------
    !1) Interpolate velocities at desired locations- mostly around the quadrilateral
    !
    !<It assumes that prog values are all synced at this stage while diag values might not>
    !--------------------------------------------------------------------------

    CALL cells2verts_scalar(p_nh_prog%w, p_patch, p_int%cells_aw_verts, w_vert,   &
                            opt_rlend=min_rlvert_int, opt_acc_async=.TRUE.)
    CALL cells2edges_scalar(p_nh_prog%w, p_patch, p_int%c_lin_e, w_ie,            &
                            opt_rlend=min_rledge_int-2, lacc=.TRUE.)

    ! RBF reconstruction of velocity at vertices: include halos
    CALL rbf_vec_interpol_vertex(p_nh_prog%vn, p_patch, p_int, u_vert, v_vert,    &
                            opt_rlend=min_rlvert_int, opt_acc_async=.TRUE.)

    !$ACC WAIT
    !sync them
    CALL sync_patch_array_mult(SYNC_V, p_patch, 3, u_vert, v_vert, w_vert)

    !Get vn at interfaces and then get vt at interfaces
    !Boundary values are extrapolated like dynamics although
    !they are not required in current implementation

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    rl_start   = 2
    rl_end     = min_rledge_int-3
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
#endif
          vn_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * p_nh_prog%vn(je,jk,jb) +            &
                            ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * p_nh_prog%vn(je,jk-1,jb)
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO je = i_startidx, i_endidx
          vn_ie(je,1,jb)      = p_nh_metrics%wgtfacq1_e(je,1,jb) * p_nh_prog%vn(je,1,jb) +        &
                                p_nh_metrics%wgtfacq1_e(je,2,jb) * p_nh_prog%vn(je,2,jb) +        &
                                p_nh_metrics%wgtfacq1_e(je,3,jb) * p_nh_prog%vn(je,3,jb)

          vn_ie(je,nlevp1,jb) = p_nh_metrics%wgtfacq_e(je,1,jb) * p_nh_prog%vn(je,nlev,jb)   +    &
                                p_nh_metrics%wgtfacq_e(je,2,jb) * p_nh_prog%vn(je,nlev-1,jb) +    &
                                p_nh_metrics%wgtfacq_e(je,3,jb) * p_nh_prog%vn(je,nlev-2,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL rbf_vec_interpol_edge(vn_ie, p_patch, p_int, vt_ie, opt_rlstart=3, &
                               opt_rlend=min_rledge_int-2, opt_acc_async=.TRUE.)

    ! RBF reconstruction of velocity on half levels at vertices: include halos
    CALL rbf_vec_interpol_vertex( vn_ie, p_patch, p_int, &
                                  u_iv, v_iv, opt_rlend=min_rlvert_int, opt_acc_async=.TRUE.)

    !$ACC WAIT
    CALL sync_patch_array_mult(SYNC_V, p_patch, 2, u_iv, v_iv)

    !--------------------------------------------------------------------------
    !2) Compute horizontal strain rate tensor at full levels
    !--------------------------------------------------------------------------
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk


!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    rl_start   = 4
    rl_end     = min_rledge_int-2
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,  &
!$OMP            vt_vert1,vt_vert2,vt_vert3,vt_vert4,w_full_c1,w_full_c2,w_full_v1,  &
!$OMP            w_full_v2,D_11, D_12, D_13, D_22, D_23, D_33)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO je = i_startidx, i_endidx
        D_11_ie(je,1,jb)       = 0.0_wp
        D_11_ie(je,nlev+1,jb)  = 0.0_wp
        D_12_ie(je,1,jb)       = 0.0_wp
        D_12_ie(je,nlev+1,jb)  = 0.0_wp
        D_13_ie(je,1,jb)       = 0.0_wp
        D_13_ie(je,nlev+1,jb)  = 0.0_wp
      END DO
      !$ACC END PARALLEL

      ! top layer approximation
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO je = i_startidx, i_endidx
        !$ACC LOOP SEQ
        DO jk = 1, 3
#else
      !$ACC LOOP SEQ
      DO jk = 1, 3
        !$ACC LOOP GANG VECTOR
        DO je = i_startidx, i_endidx
#endif
#         include "smagorinsky_strain_rate.inc"
          D_11_ie(je,1,jb)    = D_11_ie(je,1,jb) + p_nh_metrics%wgtfacq1_e(je,jk,jb) * D_11
          D_12_ie(je,1,jb)    = D_12_ie(je,1,jb) + p_nh_metrics%wgtfacq1_e(je,jk,jb) * D_12
          D_13_ie(je,1,jb)    = D_13_ie(je,1,jb) + p_nh_metrics%wgtfacq1_e(je,jk,jb) * D_13
          D_11_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_11
          D_12_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_12
          D_13_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_13
          IF (jk >= 2) THEN
            D_11_ie(je,jk,jb) = D_11_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_11
            D_12_ie(je,jk,jb) = D_12_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_12
            D_13_ie(je,jk,jb) = D_13_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_13
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO je = i_startidx, i_endidx
        !$ACC LOOP SEQ
        DO jk = 4, nlev-3
#else
      !$ACC LOOP SEQ
      DO jk = 4, nlev-3
        !$ACC LOOP GANG VECTOR
        DO je = i_startidx, i_endidx
#endif
#         include "smagorinsky_strain_rate.inc"
          D_11_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_11
          D_12_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_12
          D_13_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_13
          D_11_ie(je,jk,jb)   = D_11_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_11
          D_12_ie(je,jk,jb)   = D_12_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_12
          D_13_ie(je,jk,jb)   = D_13_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_13
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO je = i_startidx, i_endidx
        !$ACC LOOP SEQ
        DO jk = nlev-2, nlev
#else
      !$ACC LOOP SEQ
      DO jk = nlev-2, nlev
        !$ACC LOOP GANG VECTOR
        DO je = i_startidx, i_endidx
#endif
#         include "smagorinsky_strain_rate.inc"
          IF (jk <= nlev-1) THEN
            D_11_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_11
            D_12_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_12
            D_13_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_13
          END IF
          D_11_ie(je,jk,jb)     = D_11_ie(je,jk,jb) +                                             &
                                  ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_11
          D_12_ie(je,jk,jb)     = D_12_ie(je,jk,jb) +                                             &
                                  ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_12
          D_13_ie(je,jk,jb)     = D_13_ie(je,jk,jb) +                                             &
                                  ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_13
          D_11_ie(je,nlev+1,jb) = D_11_ie(je,nlev+1,jb) +                                         &
                                  p_nh_metrics%wgtfacq_e(je,nlev-jk+1,jb) * D_11
          D_12_ie(je,nlev+1,jb) = D_12_ie(je,nlev+1,jb) +                                         &
                                  p_nh_metrics%wgtfacq_e(je,nlev-jk+1,jb) * D_12
          D_13_ie(je,nlev+1,jb) = D_13_ie(je,nlev+1,jb) +                                         &
                                  p_nh_metrics%wgtfacq_e(je,nlev-jk+1,jb) * D_13
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO

    !Interpolate mech production term from mid level edge to interface level cell
    !except top and bottom boundaries
    rl_start   = 3
    rl_end     = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          mech_prod(jc,jk,jb) = p_nh_metrics%wgtfac_c(jc,jk,jb) * (                               &
                          shear(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb)   +  &
                          shear(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb)   +  &
                          shear(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) ) +  &
                          ( 1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb) ) * (                         &
                          shear(ieidx(jc,jb,1),jk-1,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +  &
                          shear(ieidx(jc,jb,2),jk-1,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +  &
                          shear(ieidx(jc,jb,3),jk-1,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT

    !div_of_stress from edge to cell-scalar interpolation
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int-1 !-1 for its use in hor diffusion
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                        i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
        div_c(jc,jk,jb) =                                                                         &
              ( div_of_stress(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +      &
                div_of_stress(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +      &
                div_of_stress(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

    !--------------------------------------------------------------------------
    !3) Classical Smagorinsky model with stability correction due to Lilly 1962
    !   at interface cell centers. At this point buoyancy_prod contains Brunt Vaisala
    !   Frequency and mech_prod is twice the actual mechanical production term. These
    !   are only converted to proper terms during output.
    !--------------------------------------------------------------------------
    !MP = Mechanical prod term calculated above
    !visc = mixing_length_sq * SQRT(MP/2) * SQRT(1-Ri/Pr) where
    !Ri = (g/theta)*d_theta_dz/(MP/2), where Brunt_vaisala_freq (byncy prod term/kh)
    !   = (g/theta)*d_theta_dz.
    !After simplification: visc = mixing_length_sq * SQRT[ MP/2 - (Brunt_vaisala_frq/Pr) ]

    rl_start   = 3
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
        kh_ic(jc,jk,jb) = rho_ic(jc,jk,jb) * les_config(jg)%rturb_prandtl *   &
                          p_nh_metrics%mixing_length_sq(jc,jk,jb)         *   &
                          SQRT( MAX( 0._wp, mech_prod(jc,jk,jb) * 0.5_wp  -   &
                          les_config(jg)%rturb_prandtl*bruvais(jc,jk,jb) ) )

        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
        kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    CALL sync_patch_array(SYNC_C, p_patch, kh_ic)

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    !--------------------------------------------------------------------------
    !4) Interpolate difusivity (viscosity) to different locations: calculate them for
    !   halos also because they will be used later in diffusion
    !--------------------------------------------------------------------------

    !4a) visc at cell center
    rl_start   = grf_bdywidth_c
    rl_end     = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          km_c(jc,jk,jb) = MAX( les_config(jg)%km_min, &
                                ( kh_ic(jc,jk,jb) + kh_ic(jc,jk+1,jb) ) * &
                                0.5_wp * les_config(jg)%turb_prandtl )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! 4b) visc at vertices
    CALL cells2verts_scalar(kh_ic, p_patch, p_int%cells_aw_verts, km_iv, &
                            opt_rlstart=5, opt_rlend=min_rlvert_int-1, opt_acc_async=.TRUE.)

    DO jb = 1, p_patch%nblks_v
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlevp1
        DO jc = 1, nproma
          km_iv(jc,jk,jb) = MAX( les_config(jg)%km_min, km_iv(jc,jk,jb) * les_config(jg)%turb_prandtl )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO

    ! 4c) Now calculate km at half levels at edge
    CALL cells2edges_scalar(kh_ic, p_patch, p_int%c_lin_e, km_ie, &
                            opt_rlstart=grf_bdywidth_e, opt_rlend=min_rledge_int-1, lacc=.TRUE.)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = 1, p_patch%nblks_e
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = 1, nproma
        DO jk = 1, nlev+1
#else
      DO jk = 1, nlev+1
        DO jc = 1, nproma
#endif
          km_ie(jc,jk,jb) = MAX( les_config(jg)%km_min,                        &
                                 km_ie(jc,jk,jb) * les_config(jg)%turb_prandtl )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! 4d) Get km_ic
    DO jb = 1, p_patch%nblks_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlevp1
        DO jc = 1, nproma
          km_ic(jc,jk,jb) = MAX( les_config(jg)%km_min, kh_ic(jc,jk,jb) * les_config(jg)%turb_prandtl )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO

    !$ACC WAIT
    !$ACC END DATA ! shear, div_of_stress

  END SUBROUTINE smagorinsky_model
  !-------------------------------------------------------------------------------------

  SUBROUTINE prognostic_tke(p_nh_prog, p_prog_now_rcf, p_prog_rcf, p_nh_diag, p_nh_metrics,       &
                            p_patch, p_int, p_diag, dtime, D_11_ie, D_12_ie, D_13_ie)
                            
    TYPE(t_nh_prog),      INTENT(inout)     :: p_nh_prog            !< single nh prognostic state
    TYPE(t_nh_prog),      INTENT(in)        :: p_prog_now_rcf       !< old state for tke
    TYPE(t_nh_prog),      INTENT(inout)     :: p_prog_rcf           !< progs w. red. call frequency
    TYPE(t_nh_diag),      INTENT(in)        :: p_nh_diag            !< single nh diagnostic state
    TYPE(t_nh_metrics),   INTENT(in),TARGET :: p_nh_metrics         !< single nh metric state
    TYPE(t_patch),        INTENT(inout),TARGET :: p_patch              !< single patch
    TYPE(t_int_state),    INTENT(in),TARGET :: p_int                !< single interpolation state
    TYPE(t_nwp_phy_diag), INTENT(inout)     :: p_diag               !< atm phys vars
    REAL(wp),             INTENT(in)        :: dtime                !< time-step  
    REAL(wp), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_e), INTENT(out) :: &
                                               D_11_ie, D_12_ie, D_13_ie             

    INTEGER ::  i_startblk, i_endblk, i_startidx, i_endidx          !< loop variables
    INTEGER ::  rl_start, rl_end, nlev, nlevp1                      !< loop & level variables 
    INTEGER ::  jb, jc, jg, je, jk                                  !< indices 

    INTEGER, DIMENSION(:,:,:), POINTER ::                   & 
                ividx, ivblk, iecidx, iecblk, ieidx, ieblk          !< index conversion arrays

    REAL(wp) :: D_11, D_12, D_13, D_22, D_23, D_33
    REAL(wp) :: ddt_tke, dthvdz, l_stable                           !< tke-tend, thv grad, dummy
    REAL(wp) :: l_grid, mixlen_by_l_grid, mixlen                    !< mix. length & avg spacing 
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4              !< nor. velocities at vortex
    REAL(wp) :: vt_vert1, vt_vert2, vt_vert3, vt_vert4              !< tang. velociies at vortex
    REAL(wp) :: w_full_c1, w_full_c2, w_full_v1, w_full_v2          !< vert. velocities
    REAL(wp) :: ddt_tke_adv(nproma,p_patch%nlev+1,p_patch%nblks_c)  !< advection tendency 
    REAL(wp) :: diff_tke(nproma,p_patch%nlev+1,p_patch%nblks_c)     !< advection tendency 
    REAL(wp) :: div_of_stress(nproma,p_patch%nlev,p_patch%nblks_e)  !< divergence of stress
    REAL(wp) :: shear(nproma,p_patch%nlev+1,p_patch%nblks_e)        !< shear in mech. prod.          
    REAL(wp) :: thetav_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)    !< thv at half level
    REAL(wp) :: vn_ie(nproma,p_patch%nlev+1,p_patch%nblks_e)        !< vn at half lev. edges
    REAL(wp) :: vt_ie(nproma,p_patch%nlev+1,p_patch%nblks_e)        !< vt at half lev. edges  

    REAL(wp), POINTER :: km_ic(:,:,:), kh_ic(:,:,:)

    !--------------------------------------------------------------------------

    IF (msg_level >= 18) CALL message(TRIM(inmodule), 'prognostic_tke')

    jg     = p_patch%id
    nlev   = p_patch%nlev
    nlevp1 = nlev+1

    km_ic => p_diag%tkvm
    kh_ic => p_diag%tkvh

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! setting boundary conditions; should be moved to set up of LEM
    km_ic(:,1,:) = les_config(jg)%km_min
    kh_ic(:,1,:) = les_config(jg)%km_min

    CALL vert_intp_full2half_cell_3d(p_patch, p_nh_metrics, theta_v, thetav_ic,                   &
                                     rl_start, rl_end)

    !--------------------------------------------------------------------------
    !  Advection of TKE
    !
    !  Utilizes the advection mechanism used for tracers; defined at cell center.
    !  Tendency (mo_step_advection) is interpolated to half-level. 
    !--------------------------------------------------------------------------
    CALL vert_intp_full2half_cell_3d(p_patch, p_nh_metrics,                                       &
                                     p_nh_diag%ddt_tracer_adv(:,:,:,iqtke), ddt_tke_adv,          &
                                     rl_start, rl_end)

    !--------------------------------------------------------------------------
    !  Diffusion term
    !
    !  The diffusion term is parameterized (Deardorff, 1980), but with variable
    !  density since non-hydrostatic (as in eq. 37 in Bryan's CM1(2017)).
    !  Adapted from diffusion of scalars on full level, has included interpolation
    !  to half-level.
    !--------------------------------------------------------------------------
    CALL diffuse_tke(p_prog_now_rcf%tke, p_prog_rcf%tracer(:,:,:,iqtke), p_nh_metrics,            &
                     p_patch,p_int, diff_tke, p_diag, dtime)

    !--------------------------------------------------------------------------
    !  Shear term
    !
    !  The method to compute the mechanical production is adopted from the   
    !  Smagorinsky scheme: 
    !  mech_prod   = K_m * D**2     with:
    !  D**2        = 2 D_ij * D_ij =  (D_11**2 + D_22**2 + D_33**2 + 
    !                                 2 * (D_12**2 + D_13**2 + D_23**2) 
    !  where, D_11 = 2 * du_1/dx_1
    !         D_22 = 2 * d_u2/dx_2
    !         D_33 = 2 * d_u3/dx_3
    !         D_12 = du_1/dx_2 + du_2/dx_1
    !         D_13 = du_1/dx_3 + du_3/dx_1
    !         D_23 = du_2/dx_3 + du_3/dx_2
    !  For triangles: 1=normal, 2=tangential, and 3 = z directions 
    !--------------------------------------------------------------------------
    CALL cells2verts_scalar(p_nh_prog%w, p_patch, p_int%cells_aw_verts, w_vert,                   &
                            opt_rlend=min_rlvert_int)
    CALL cells2edges_scalar(p_nh_prog%w, p_patch, p_int%c_lin_e, w_ie,                            &
                            opt_rlend=min_rledge_int-2)
    CALL rbf_vec_interpol_vertex(p_nh_prog%vn, p_patch, p_int, u_vert, v_vert,                    &
                                  opt_rlend=min_rlvert_int)
    CALL sync_patch_array_mult(SYNC_V, p_patch, 3, w_vert, u_vert, v_vert)
 
    rl_start   = grf_bdywidth_e+1
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    rl_start   = 2
    rl_end     = min_rledge_int-3
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
#endif
          vn_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * p_nh_prog%vn(je,jk,jb) +            &
                            ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * p_nh_prog%vn(je,jk-1,jb)
        END DO
      END DO
      DO je = i_startidx, i_endidx
          vn_ie(je,1,jb)      = p_nh_metrics%wgtfacq1_e(je,1,jb) * p_nh_prog%vn(je,1,jb) +        &
                                p_nh_metrics%wgtfacq1_e(je,2,jb) * p_nh_prog%vn(je,2,jb) +        &
                                p_nh_metrics%wgtfacq1_e(je,3,jb) * p_nh_prog%vn(je,3,jb)

          vn_ie(je,nlevp1,jb) = p_nh_metrics%wgtfacq_e(je,1,jb) * p_nh_prog%vn(je,nlev,jb)   +    &
                                p_nh_metrics%wgtfacq_e(je,2,jb) * p_nh_prog%vn(je,nlev-1,jb) +    &
                                p_nh_metrics%wgtfacq_e(je,3,jb) * p_nh_prog%vn(je,nlev-2,jb)
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E, p_patch, vn_ie)

    CALL rbf_vec_interpol_edge(vn_ie, p_patch, p_int, vt_ie, opt_rlstart=3, &
                               opt_rlend=min_rledge_int-2)

    ! RBF reconstruction of velocity on half levels at vertices: include halos
    CALL rbf_vec_interpol_vertex( vn_ie, p_patch, p_int, &
                                  u_iv, v_iv, opt_rlend=min_rlvert_int )

    CALL sync_patch_array_mult(SYNC_V, p_patch, 2, u_iv, v_iv)

    !--------------------------------------------------------------------------
    !2) Compute horizontal strain rate tensor at full levels
    !--------------------------------------------------------------------------
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk


!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    rl_start   = 4
    rl_end     = min_rledge_int-2
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,  &
!$OMP            vt_vert1,vt_vert2,vt_vert3,vt_vert4,w_full_c1,w_full_c2,w_full_v1,  &
!$OMP            w_full_v2)
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)

      DO je = i_startidx, i_endidx
        D_11_ie(je,1,jb)       = 0.0_wp
        D_11_ie(je,nlev+1,jb)  = 0.0_wp
        D_12_ie(je,1,jb)       = 0.0_wp
        D_12_ie(je,nlev+1,jb)  = 0.0_wp
        D_13_ie(je,1,jb)       = 0.0_wp
        D_13_ie(je,nlev+1,jb)  = 0.0_wp
      END DO
      ! top layer approximation
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, 3
#else
      DO jk = 1, 3
        DO je = i_startidx, i_endidx
#endif
#         include "smagorinsky_strain_rate.inc"
          D_11_ie(je,1,jb)    = D_11_ie(je,1,jb) + p_nh_metrics%wgtfacq1_e(je,jk,jb) * D_11
          D_12_ie(je,1,jb)    = D_12_ie(je,1,jb) + p_nh_metrics%wgtfacq1_e(je,jk,jb) * D_12
          D_13_ie(je,1,jb)    = D_13_ie(je,1,jb) + p_nh_metrics%wgtfacq1_e(je,jk,jb) * D_13
          D_11_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_11
          D_12_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_12
          D_13_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_13
          IF (jk >= 2) THEN
            D_11_ie(je,jk,jb) = D_11_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_11
            D_12_ie(je,jk,jb) = D_12_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_12
            D_13_ie(je,jk,jb) = D_13_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_13
          END IF
        END DO
      END DO

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 4, nlev-3
#else
      DO jk = 4, nlev-3
        DO je = i_startidx, i_endidx
#endif
#         include "smagorinsky_strain_rate.inc"
          D_11_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_11
          D_12_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_12
          D_13_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_13
          D_11_ie(je,jk,jb)   = D_11_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_11
          D_12_ie(je,jk,jb)   = D_12_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_12
          D_13_ie(je,jk,jb)   = D_13_ie(je,jk,jb) +                                               &
                                ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_13
        ENDDO
      ENDDO
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = nlev-2, nlev
#else
      DO jk = nlev-2, nlev
        DO je = i_startidx, i_endidx
#endif
#         include "smagorinsky_strain_rate.inc"
          IF (jk <= nlev-1) THEN
            D_11_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_11
            D_12_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_12
            D_13_ie(je,jk+1,jb) = p_nh_metrics%wgtfac_e(je,jk+1,jb) * D_13
          END IF
          D_11_ie(je,jk,jb)     = D_11_ie(je,jk,jb) +                                             &
                                  ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_11
          D_12_ie(je,jk,jb)     = D_12_ie(je,jk,jb) +                                             &
                                  ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_12
          D_13_ie(je,jk,jb)     = D_13_ie(je,jk,jb) +                                             &
                                  ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * D_13
          D_11_ie(je,nlev+1,jb) = D_11_ie(je,nlev+1,jb) +                                         &
                                  p_nh_metrics%wgtfacq_e(je,nlev-jk+1,jb) * D_11
          D_12_ie(je,nlev+1,jb) = D_12_ie(je,nlev+1,jb) +                                         &
                                  p_nh_metrics%wgtfacq_e(je,nlev-jk+1,jb) * D_12
          D_13_ie(je,nlev+1,jb) = D_13_ie(je,nlev+1,jb) +                                         &
                                  p_nh_metrics%wgtfacq_e(je,nlev-jk+1,jb) * D_13
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL  

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int !-1 for its use in hor. diffusion
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! Interpolation of div_of_stress from edge to cell-center 
!$OMP PARALLEL    
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 1, nlev
#else
       DO jk = 1, nlev
         DO jc = i_startidx, i_endidx
#endif
           div_c(jc,jk,jb) =                                                                      &
                  ( div_of_stress(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +  &
                    div_of_stress(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +  &
                    div_of_stress(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
         END DO
       END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,l_grid,l_stable,mixlen,    &
!$OMP            mixlen_by_l_grid,ddt_tke,dthvdz )
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          ddt_tke = km_ic(jc,jk,jb) * ( p_nh_metrics%wgtfac_c(jc,jk,jb) * (                   &
                    shear(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb)   +    &
                    shear(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb)   +    &
                    shear(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) ) +    &
                    ( 1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb) ) * (                           &  
                    shear(ieidx(jc,jb,1),jk-1,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +    &
                    shear(ieidx(jc,jb,2),jk-1,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +    &
                    shear(ieidx(jc,jb,3),jk-1,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) ) )

          !--------------------------------------------------------------------------
          !  Buoyancy term
          !
          !  TKE production by buoyancy: buoy = - K_h * g * dtheta_v/dz / theta_v
          !-------------------------------------------------------------------------- 
          ! Temperature gradient over two levels due to better numerical stability
          dthvdz  = ( thetav_ic(jc,jk-1,jb) - thetav_ic(jc,jk+1,jb) ) /                         &
                    ( p_nh_metrics%ddqz_z_half(jc,jk,jb) + p_nh_metrics%ddqz_z_half(jc,jk+1,jb) )   

          ddt_tke = ddt_tke - ( kh_ic(jc,jk,jb) * grav * dthvdz / thetav_ic(jc,jk,jb) )      

          !--------------------------------------------------------------------------
          !  Dissipation term 
          !
          !  Takes into account stratification and grid spacing. 
          !--------------------------------------------------------------------------
          l_grid = ( p_nh_metrics%ddqz_z_half(jc,jk,jb) * p_patch%cells%area(jc,jb) )**0.333333_wp       

          ! mixing length either avg. grid spacing or stratification dependent
          IF ( dthvdz > 0.0_wp ) THEN  
            l_stable   = 0.76_wp * SQRT( p_prog_now_rcf%tke(jc,jk,jb) /                           &
                         ( grav * dthvdz / thetav_ic(jc,jk,jb) ) ) + 1E-5_wp
          ELSE
            l_stable = l_grid
          END IF

          IF (jk == nlev) THEN
            mixlen   = MIN( l_grid, l_stable, 1.8_wp * p_nh_metrics%ddqz_z_half(jc,nlev,jb) )
            l_grid   = MIN( l_grid, 1.8_wp * p_nh_metrics%ddqz_z_half(jc,nlev,jb))         
          ELSE
            mixlen           = MIN( l_grid, l_stable)
          END IF

          mixlen_by_l_grid = mixlen / l_grid
          ddt_tke          = ddt_tke - ( 0.19_wp + 0.74_wp * mixlen_by_l_grid ) *                 &
                             p_prog_now_rcf%tke(jc,jk,jb) *                                       &
                             SQRT( p_prog_now_rcf%tke(jc,jk,jb) ) / mixlen

          !--------------------------------------------------------------------------
          !  Time-stepping of the TKE
          !--------------------------------------------------------------------------
          p_prog_rcf%tke(jc,jk,jb)   = p_prog_now_rcf%tke(jc,jk,jb) + dtime * (                   &
                                       ddt_tke_adv(jc,jk,jb) + ddt_tke + diff_tke(jc,jk,jb) )

          IF (p_prog_rcf%tke(jc,jk,jb) < 0._wp) THEN
            p_prog_rcf%tke(jc,jk,jb) = 0.1_wp * p_prog_now_rcf%tke(jc,jk,jb)
          END IF

          km_ic(jc,jk,jb) = 0.1_wp * mixlen * SQRT( p_prog_now_rcf%tke(jc,jk,jb) ) 
          kh_ic(jc,jk,jb) = ( 1._wp + 2._wp * mixlen_by_l_grid ) * km_ic(jc,jk,jb) 
        END DO
      END DO 

      DO jc = i_startidx, i_endidx
        p_prog_rcf%tke(jc,nlevp1,jb) = p_prog_rcf%tke(jc,nlev,jb)                                
        km_ic(jc,nlevp1,jb)          = km_ic(jc,nlev,jb)                    
        kh_ic(jc,nlevp1,jb)          = kh_ic(jc,nlev,jb)
        p_prog_rcf%tke(jc,1,jb)      = p_prog_rcf%tke(jc,2,jb)
      ENDDO 

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          p_prog_rcf%tracer(jc,jk,jb,iqtke) = 0.5_wp * ( p_prog_rcf%tke(jc,jk,jb) +               &
                                                         p_prog_rcf%tke(jc,jk+1,jb) )

          km_c(jc,jk,jb) = MAX( les_config(jg)%km_min,                                            &
                                0.5_wp * ( km_ic(jc,jk,jb) + km_ic(jc,jk+1,jb) ) )                    
        ENDDO  
      ENDDO     
    END DO !jb
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array_mult(SYNC_C, p_patch, 3, km_ic, kh_ic, p_prog_rcf%tke)

    ! Interpolate diffusivity (viscosity) to different locations: calculate  
    ! them for halos also because they will be used in diffusion
    CALL cells2verts_scalar(km_ic, p_patch, p_int%cells_aw_verts, km_iv, opt_rlstart=5,           &
                            opt_rlend=min_rlvert_int-1)
    CALL cells2edges_scalar(km_ic, p_patch, p_int%c_lin_e, km_ie, opt_rlstart=grf_bdywidth_e,     &
                            opt_rlend=min_rledge_int-1)
    CALL cells2edges_scalar(kh_ic, p_patch, p_int%c_lin_e, kh_ie, opt_rlstart=grf_bdywidth_e,     &
                            opt_rlend=min_rledge_int-1)

    ! Assure that diffusivity and viscosity have minimum values
    ! May be obsolete if km_ic & kh_ic are initialized with les_config(jg)%km_min
    km_ic = MAX( les_config(jg)%km_min, km_ic )
    km_iv = MAX( les_config(jg)%km_min, km_iv )
    km_ie = MAX( les_config(jg)%km_min, km_ie )
    kh_ie = MAX( les_config(jg)%km_min, kh_ie )

  END SUBROUTINE prognostic_tke
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_hori_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for normal velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !!
  !! d_vn/d_t =  d_tau_11/d_x1 + d_tau_12/d_x2 + d_tau_13/d_x3
  !!------------------------------------------------------------------------
  SUBROUTINE diffuse_hori_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
                                   prm_diag, ddt_u, ddt_v, &
                                   D_11_ie, D_12_ie, D_13_ie, dt)

    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog    !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in)        :: p_nh_metrics !< single nh metric state
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    TYPE(t_nwp_phy_diag),INTENT(inout)   :: prm_diag     !< atm phys vars
    REAL(wp),   TARGET, INTENT(inout)    :: ddt_u(:,:,:) !< u tendency
    REAL(wp),   TARGET, INTENT(inout)    :: ddt_v(:,:,:) !< v tendency
    REAL(wp), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_e), INTENT(in) :: &
                                            D_11_ie, D_12_ie, D_13_ie
    REAL(wp),           INTENT(in)       :: dt           !< dt turb

    REAL(wp) :: flux_up_e, flux_dn_e, flux_up_v, flux_dn_v, flux_up_c,      &
                flux_dn_c, vflux_up, vflux_dn
    REAL(wp) :: norm_metr, tang_metr, vert_metr, div_of_stress, div_of_stress_p1
    REAL(wp) :: stress_c1n, stress_c2n
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt, dwdn, inv_dt,  &
                vn_vert1i, vn_vert1i_p1, vn_vert2i, vn_vert2i_p1,           &
                vn_vert3i, vn_vert4i, vn_vert3i_p1, vn_vert4i_p1,           &
                vt_vert1i, vt_vert1i_p1, vt_vert2i, vt_vert2i_p1,           &
                vt_vert3i, vt_vert3i_p1, vt_vert4i, vt_vert4i_p1
    REAL(wp) :: inv_rhoe(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: vn_new(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: unew(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: vnew(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: tot_tend(nproma,p_patch%nlev,p_patch%nblks_e)

    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs, var_new
    REAL(wp), DIMENSION(p_patch%nlev) :: outvar

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jcn, jbn, jvn, jc
    INTEGER :: nlev, jg, nlevp1

    IF (msg_level >= 18) &
         CALL message(inmodule, 'diffuse_hori_velocity')

    !patch id
    jg = p_patch%id

    ! number of vertical levels
    nlev     = p_patch%nlev
    nlevp1   = nlev+1

    inv_dt  = 1._wp / dt

    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    !$ACC DATA CREATE(inv_rhoe, vn_new, unew, vnew, tot_tend, a, b, c, rhs, var_new, outvar)

!$OMP PARALLEL
    CALL init(a, lacc=.TRUE.)
    CALL init(c, lacc=.TRUE.)
    CALL init(tot_tend, lacc=.TRUE.)
    CALL copy(p_nh_prog%vn, vn_new, lacc=.TRUE.)
    IF (p_test_run) THEN
      CALL init(inv_rhoe, lacc=.TRUE.)
    END IF
!$OMP END PARALLEL

    rl_start   = grf_bdywidth_e+1
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

    !density at edge
    CALL cells2edges_scalar(p_nh_prog%rho, p_patch, p_int%c_lin_e, inv_rhoe, &
                            opt_rlstart=rl_start, opt_rlend=rl_end, lacc=.TRUE.)

    ! compute inv_rho
    ! compute D_11, D_12 on interface edges
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE   
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          inv_rhoe(je,jk,jb) = 1._wp / inv_rhoe(je,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! 1) First get the horizontal tendencies

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,&
!$OMP   vn_vert1i,vn_vert2i,vn_vert3i,vn_vert4i,vt_vert1i,vt_vert2i,&
!$OMP   vn_vert1i_p1,vn_vert2i_p1,vn_vert3i_p1,vn_vert4i_p1,vt_vert1i_p1,&
!$OMP   vt_vert2i_p1,vt_vert3i, vt_vert3i_p1, vt_vert4i, vt_vert4i_p1,&
!$OMP   dvt,jcn,jbn,flux_up_c,flux_dn_c,jvn,flux_up_v,flux_dn_v,&
!$OMP   norm_metr,tang_metr,div_of_stress,div_of_stress_p1)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(vn_vert1, vn_vert1i, vn_vert1i_p1, vn_vert2, vn_vert2i, vn_vert2i_p1) &
      !$ACC   PRIVATE(vn_vert3, vn_vert3i, vn_vert3i_p1, vn_vert4, vn_vert4i, vn_vert4i_p1) &
      !$ACC   PRIVATE(vt_vert1i, vt_vert1i_p1, vt_vert2i, vt_vert2i_p1) &
      !$ACC   PRIVATE(vt_vert3i, vt_vert3i_p1, vt_vert4i, vt_vert4i_p1) &
      !$ACC   PRIVATE(dvt, jcn, jbn, jvn, flux_up_c, flux_dn_c, flux_up_v, flux_dn_v) &
      !$ACC   PRIVATE(div_of_stress, div_of_stress_p1, norm_metr, tang_metr)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
        vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))         * &
                   p_patch%edges%primal_normal_vert(je,jb,1)%v1     + &
                   v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))         * &
                   p_patch%edges%primal_normal_vert(je,jb,1)%v2

        vn_vert1i = u_iv(ividx(je,jb,1),jk,ivblk(je,jb,1))          * &
                    p_patch%edges%primal_normal_vert(je,jb,1)%v1    + &
                    v_iv(ividx(je,jb,1),jk,ivblk(je,jb,1))          * &
                    p_patch%edges%primal_normal_vert(je,jb,1)%v2

        vn_vert1i_p1 = u_iv(ividx(je,jb,1),jk+1,ivblk(je,jb,1))     * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                       v_iv(ividx(je,jb,1),jk+1,ivblk(je,jb,1))     * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v2

        vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))         * &
                   p_patch%edges%primal_normal_vert(je,jb,2)%v1     + &
                   v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))         * &
                   p_patch%edges%primal_normal_vert(je,jb,2)%v2

        vn_vert2i = u_iv(ividx(je,jb,2),jk,ivblk(je,jb,2))          * &
                    p_patch%edges%primal_normal_vert(je,jb,2)%v1    + &
                    v_iv(ividx(je,jb,2),jk,ivblk(je,jb,2))          * &
                    p_patch%edges%primal_normal_vert(je,jb,2)%v2

        vn_vert2i_p1 = u_iv(ividx(je,jb,2),jk+1,ivblk(je,jb,2))     * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                       v_iv(ividx(je,jb,2),jk+1,ivblk(je,jb,2))     * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v2

        vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))         * &
                   p_patch%edges%primal_normal_vert(je,jb,3)%v1     + &
                   v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))         * &
                   p_patch%edges%primal_normal_vert(je,jb,3)%v2

        vn_vert3i = u_iv(ividx(je,jb,3),jk,ivblk(je,jb,3))          * &
                    p_patch%edges%primal_normal_vert(je,jb,3)%v1    + &
                    v_iv(ividx(je,jb,3),jk,ivblk(je,jb,3))          * &
                    p_patch%edges%primal_normal_vert(je,jb,3)%v2

        vn_vert3i_p1 = u_iv(ividx(je,jb,3),jk+1,ivblk(je,jb,3))     * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                       v_iv(ividx(je,jb,3),jk+1,ivblk(je,jb,3))     * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v2

        vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))         * &
                   p_patch%edges%primal_normal_vert(je,jb,4)%v1     + &
                   v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))         * &
                   p_patch%edges%primal_normal_vert(je,jb,4)%v2

        vn_vert4i = u_iv(ividx(je,jb,4),jk,ivblk(je,jb,4))          * &
                    p_patch%edges%primal_normal_vert(je,jb,4)%v1    + &
                    v_iv(ividx(je,jb,4),jk,ivblk(je,jb,4))          * &
                    p_patch%edges%primal_normal_vert(je,jb,4)%v2

        vn_vert4i_p1 = u_iv(ividx(je,jb,4),jk+1,ivblk(je,jb,4))     * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                       v_iv(ividx(je,jb,4),jk+1,ivblk(je,jb,4))     * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v2

        vt_vert1i = u_iv(ividx(je,jb,1),jk,ivblk(je,jb,1))          * &
                    p_patch%edges%dual_normal_vert(je,jb,1)%v1      + &
                    v_iv(ividx(je,jb,1),jk,ivblk(je,jb,1))          * &
                    p_patch%edges%dual_normal_vert(je,jb,1)%v2

        vt_vert1i_p1 = u_iv(ividx(je,jb,1),jk+1,ivblk(je,jb,1))     * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v1   + &
                       v_iv(ividx(je,jb,1),jk+1,ivblk(je,jb,1))     * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v2

        vt_vert2i = u_iv(ividx(je,jb,2),jk,ivblk(je,jb,2))          * &
                    p_patch%edges%dual_normal_vert(je,jb,2)%v1      + &
                    v_iv(ividx(je,jb,2),jk,ivblk(je,jb,2))          * &
                    p_patch%edges%dual_normal_vert(je,jb,2)%v2

        vt_vert2i_p1 = u_iv(ividx(je,jb,2),jk+1,ivblk(je,jb,2))     * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v1   + &
                       v_iv(ividx(je,jb,2),jk+1,ivblk(je,jb,2))     * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v2

        vt_vert3i = u_iv(ividx(je,jb,3),jk,ivblk(je,jb,3))          * &
                    p_patch%edges%dual_normal_vert(je,jb,3)%v1      + &
                    v_iv(ividx(je,jb,3),jk,ivblk(je,jb,3))          * &
                    p_patch%edges%dual_normal_vert(je,jb,3)%v2

        vt_vert3i_p1 = u_iv(ividx(je,jb,3),jk+1,ivblk(je,jb,3))     * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v1   + &
                       v_iv(ividx(je,jb,3),jk+1,ivblk(je,jb,3))     * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v2

        vt_vert4i = u_iv(ividx(je,jb,4),jk,ivblk(je,jb,4))          * &
                    p_patch%edges%dual_normal_vert(je,jb,4)%v1      + &
                    v_iv(ividx(je,jb,4),jk,ivblk(je,jb,4))          * &
                    p_patch%edges%dual_normal_vert(je,jb,4)%v2

        vt_vert4i_p1 = u_iv(ividx(je,jb,4),jk+1,ivblk(je,jb,4))     * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v1   + &
                       v_iv(ividx(je,jb,4),jk+1,ivblk(je,jb,4))     * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v2

        dvt      = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))         * &
                   p_patch%edges%dual_normal_vert(je,jb,4)%v1       + &
                   v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))         * &
                   p_patch%edges%dual_normal_vert(je,jb,4)%v2       - &
                   ( u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))       * &
                   p_patch%edges%dual_normal_vert(je,jb,3)%v1       + &
                   v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))         * &
                   p_patch%edges%dual_normal_vert(je,jb,3)%v2 )


        !tendency in normal direction:
        !flux = visc*(D_11-2/3DIV) = visc*(2*delta_v/(vert_vert_len/2)-2/3*div_of_stress)

        jcn     = iecidx(je,jb,2)
        jbn     = iecblk(je,jb,2)
        flux_up_c = km_c(jcn,jk,jbn) * &
                   ( 4._wp * ( vn_vert4 - p_nh_prog%vn(je,jk,jb) )  * &
                     p_patch%edges%inv_vert_vert_length(je,jb)      - &
                     z_2by3 * ( vn_vert4i + 2._wp * vn_ie(je,jk,jb) - &
                     vn_vert4i_p1 - 2._wp * vn_ie(je,jk+1,jb) )     * &
                     p_nh_metrics%inv_ddqz_z_full(jcn,jk,jbn)       * &
                     p_nh_metrics%ddxn_z_full_c(jcn,jk,jbn)         - &
                     z_1by3*div_c(jcn,jk,jbn) )

         jcn     = iecidx(je,jb,1)
         jbn     = iecblk(je,jb,1)
         flux_dn_c = km_c(jcn,jk,jbn) * &
                   ( 4._wp*( p_nh_prog%vn(je,jk,jb) - vn_vert3 )    * &
                     p_patch%edges%inv_vert_vert_length(je,jb)      - &
                     z_2by3 * ( vn_vert3i + 2._wp * vn_ie(je,jk,jb) - &
                     vn_vert3i_p1 - 2._wp * vn_ie(je,jk+1,jb) )     * &
                     p_nh_metrics%inv_ddqz_z_full(jcn,jk,jbn)       * &
                     p_nh_metrics%ddxn_z_full_c(jcn,jk,jbn)         - &
                     z_1by3*div_c(jcn,jk,jbn) )

         !tendency in tangential direction

         !D_12 between edge center and the vertex: delta_v/(primal_edge_len/2) +
         ! ((vt4+vt2)/2-(vt3+vt2)/2)/(distance_opp_edges)
         !flux = D_12*visc

         !Note that the tangential velocities at vertices are used in D_12 is an
         !approximation for speed. Better way is to use vt reconsttructed from vn at
         !each edges. Also, visc at somewhere between edge mid point and the vertex
         !should be used but this is a good approximation
         jvn     = ividx(je,jb,2)
         jbn     = ivblk(je,jb,2)
         flux_up_v = 0.5_wp * ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) ) *                       &
            ( p_patch%edges%tangent_orientation(je,jb) * ( vn_vert2 - p_nh_prog%vn(je,jk,jb) )*   &
              p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp - &
              ( vn_vert2i - vn_vert2i_p1 )                * &
              p_nh_metrics%inv_ddqz_z_full_v(jvn,jk,jbn)  * &
              p_nh_metrics%ddxt_z_full_v(jvn,jk,jbn)      - &
              ( vt_vert2i - vt_vert2i_p1 )                * &
              p_nh_metrics%inv_ddqz_z_full_v(jvn,jk,jbn)  * &
              p_nh_metrics%ddxn_z_full_v(jvn,jk,jbn)      + &
              dvt * p_patch%edges%inv_vert_vert_length(je,jb) )

         jvn     = ividx(je,jb,1)
         jbn     = ivblk(je,jb,1)
         flux_dn_v = 0.5_wp * (km_iv(jvn,jk,jbn)+km_iv(jvn,jk+1,jbn)) *                           &
            ( p_patch%edges%tangent_orientation(je,jb) * ( p_nh_prog%vn(je,jk,jb) - vn_vert1 ) *  &
              p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp - &
              ( vn_vert1i-vn_vert1i_p1 ) * &
              p_nh_metrics%inv_ddqz_z_full_v(jvn,jk,jbn)  * &
              p_nh_metrics%ddxt_z_full_v(jvn,jk,jbn)      - &
              ( vt_vert1i - vt_vert1i_p1 )                * &
              p_nh_metrics%inv_ddqz_z_full_v(jvn,jk,jbn)  * &
              p_nh_metrics%ddxn_z_full_v(jvn,jk,jbn)      + &
              dvt*p_patch%edges%inv_vert_vert_length(je,jb) )

         div_of_stress = z_1by3 * ( D_11_ie(je,jk,jb) + D_12_ie(je,jk,jb) + D_13_ie(je,jk,jb) )
         div_of_stress_p1 = z_1by3 * ( D_11_ie(je,jk+1,jb) + D_12_ie(je,jk+1,jb) +                &
                            D_13_ie(je,jk+1,jb) )

         norm_metr = -p_nh_metrics%ddxn_z_full(je,jk,jb)              * &
                      p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)        * &
                      ( ( km_ie(je,jk,jb) * div_of_stress             - &
                          km_ie(je,jk+1,jb) * div_of_stress_p1 )      + &
                      p_patch%edges%inv_dual_edge_length(je,jb)       * &
                      ( km_ie(je,jk,jb) * ( vn_vert4i - vn_vert3i )   - &
                        km_ie(je,jk+1,jb) * ( vn_vert4i_p1 - vn_vert3i_p1 ) ) )

         tang_metr = -p_nh_metrics%ddxt_z_full(je,jk,jb)              * &          
                      p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)        * &
                      ( p_patch%edges%inv_primal_edge_length(je,jb)   * &
                      ( ( km_ie(je,jk,jb) * ( vn_vert2i - vn_vert1i ) - &
                          km_ie(je,jk+1,jb) * ( vn_vert2i_p1 - vn_vert1i_p1 ) ) ) + &
                      p_patch%edges%inv_vert_vert_length(je,jb)       * &
                      ( ( km_ie(je,jk,jb) * ( vt_vert4i - vt_vert3i ) - &
                          km_ie(je,jk+1,jb) * ( vt_vert4i_p1 - vt_vert3i_p1) ) ) )

         tot_tend(je,jk,jb) = ( ( flux_up_c - flux_dn_c ) *                 &
                                p_patch%edges%inv_dual_edge_length(je,jb) + &
                                p_patch%edges%tangent_orientation(je,jb)  * &
                                ( flux_up_v - flux_dn_v )                 * &
                                p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp - &
                                norm_metr - tang_metr ) * inv_rhoe(je,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !
    ! 2) Vertical tendency
    !
    SELECT CASE(les_config(jg)%vert_scheme_type)

    CASE(iexplicit)
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'vert_scheme_type=iexplicit: OpenACC version currently not implemented')
#endif

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,flux_up_e,flux_dn_e,stress_c1n,stress_c2n,&
!$OMP vert_metr,vflux_up,vflux_dn)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev-1
#else
      DO jk = 2, nlev-1
        DO je = i_startidx, i_endidx
#endif
         !tendency in vertical direction
          flux_up_e = km_ie(je,jk,jb) *                                        &
                      ( ( p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2))  -   &
                          p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1))) *   &
                          p_patch%edges%inv_dual_edge_length(je,jb)        -   &
                         0.5_wp * ( w_ie(je,jk-1,jb) - w_ie(je,jk+1,jb) )  *   &
                         p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb)          *   &
                         p_nh_metrics%ddxn_z_half_e(je,jk,jb) )

          flux_dn_e = km_ie(je,jk+1,jb) *                                      &
                      ( ( p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2) ) - &
                          p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1))) * &
                          p_patch%edges%inv_dual_edge_length(je,jb)          - &
                          0.5_wp * ( w_ie(je,jk,jb) - w_ie(je,jk+2,jb) )     * &
                          p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb)           * &
                          p_nh_metrics%ddxn_z_half_e(je,jk+1,jb) )

          vflux_up = km_ie(je,jk,jb) *                                         &
                     ( ( p_nh_prog%vn(je,jk-1,jb) - p_nh_prog%vn(je,jk,jb) ) * &
                        p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) )

          vflux_dn = km_ie(je,jk+1,jb) *                                       &
                     ( ( p_nh_prog%vn(je,jk,jb) - p_nh_prog%vn(je,jk+1,jb) ) * &
                        p_nh_metrics%inv_ddqz_z_half_e(je,jk+1,jb) )

          vert_metr = 0.5_wp * p_nh_metrics%ddxt_z_full(je,jk,jb) *              &
                      p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)    *              &
                      ( p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb)  *              &
                        km_ie(je,jk,jb) * p_nh_metrics%ddxn_z_half_e(je,jk,jb) * &
                        ( vt_ie(je,jk-1,jb) - vt_ie(je,jk+1,jb) ) -              &
                        p_nh_metrics%inv_ddqz_z_half_e(je,jk+1,jb) *             &
                        km_ie(je,jk+1,jb) * p_nh_metrics%ddxn_z_half_e(je,jk+1,jb) * &
                        ( vt_ie(je,jk,jb) - vt_ie(je,MIN(nlev+1,jk+2),jb) ) )

          tot_tend(je,jk,jb) = tot_tend(je,jk,jb) + p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) *    &
                               ( flux_up_e - flux_dn_e + vflux_up - vflux_dn +                    &
                                 2._wp * ( vflux_up * p_nh_metrics%ddxn_z_half_e(je,jk,jb) -      &
                                 vflux_dn * p_nh_metrics%ddxn_z_half_e(je,jk+1,jb) ) *            &
                                 p_nh_metrics%ddxn_z_full(je,jk,jb) +                             &
                                 ( vflux_up*p_nh_metrics%ddxt_z_half_e(je,jk,jb) -                &
                                   vflux_dn*p_nh_metrics%ddxt_z_half_e(je,jk+1,jb) ) *            &
                                 p_nh_metrics%ddxt_z_full(je,jk,jb) + vert_metr ) *               &
                               inv_rhoe(je,jk,jb)
        END DO
      END DO

      !-----------------------------------------------------------------
      !3) Boundary treatment in vertical: use surface fluxes from
      !   surface_conditions
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !jk = 1
      !-----------------------------------------------------------------
      DO je = i_startidx, i_endidx
        flux_dn_e = km_ie(je,2,jb) *                                      &
                    ( ( p_nh_prog%vn(je,1,jb) - p_nh_prog%vn(je,2,jb)) *  &
                      p_nh_metrics%inv_ddqz_z_half_e(je,2,jb)  +          &
                      ( p_nh_prog%w(iecidx(je,jb,2),2,iecblk(je,jb,2)) -  &
                      p_nh_prog%w(iecidx(je,jb,1),2,iecblk(je,jb,1)) ) *  &
                      p_patch%edges%inv_dual_edge_length(je,jb) -         &
                      0.5_wp * ( w_ie(je,1,jb) - w_ie(je,3,jb) ) *        &
                      p_nh_metrics%inv_ddqz_z_half_e(je,2,jb) *           &
                      p_nh_metrics%ddxn_z_half_e(je,2,jb) )

        vflux_dn = km_ie(je,2,jb) * ( ( p_nh_prog%vn(je,1,jb) - p_nh_prog%vn(je,2,jb) ) *         &
                                       p_nh_metrics%inv_ddqz_z_half_e(je,2,jb) )

        vert_metr = 0.5_wp * p_nh_metrics%ddxt_z_full(je,1,jb) *              &
                     p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) *                &
                     ( p_nh_metrics%inv_ddqz_z_half_e(je,1,jb) *              &
                     km_ie(je,1,jb) * p_nh_metrics%ddxn_z_half_e(je,1,jb) *   &
                     ( vt_ie(je,1,jb)-vt_ie(je,2,jb) ) -                      &
                     p_nh_metrics%inv_ddqz_z_half_e(je,2,jb) *                &
                     km_ie(je,2,jb) * p_nh_metrics%ddxn_z_half_e(je,2,jb) *   &
                     ( vt_ie(je,1,jb) - vt_ie(je,3,jb) ) )

        tot_tend(je,1,jb) = tot_tend(je,1,jb) + p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) *         &
                            ( - flux_dn_e - vflux_dn -                                            &
                              2._wp * vflux_dn*p_nh_metrics%ddxn_z_half_e(je,2,jb) *              &
                              p_nh_metrics%ddxn_z_full(je,1,jb) -                                 &
                              vflux_dn * p_nh_metrics%ddxt_z_half_e(je,2,jb) *                    &
                              p_nh_metrics%ddxt_z_full(je,1,jb) +                                 &
                              vert_metr ) * inv_rhoe(je,1,jb)
        END DO

        !-----------------------------------------------------------------
        ! jk = nlev
        !-----------------------------------------------------------------
        DO je = i_startidx, i_endidx
          flux_up_e = km_ie(je,nlev,jb) *                                           &
                      ( ( p_nh_prog%vn(je,nlev-1,jb) - p_nh_prog%vn(je,nlev,jb) ) * &
                        p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) +                &
                        ( p_nh_prog%w(iecidx(je,jb,2),nlev,iecblk(je,jb,2) ) -      &
                        p_nh_prog%w(iecidx(je,jb,1),nlev,iecblk(je,jb,1))) *        &
                        p_patch%edges%inv_dual_edge_length(je,jb) -                 &
                        0.5_wp * ( w_ie(je,nlev-1,jb) - w_ie(je,nlevp1,jb) *        &
                        p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) *                &
                        p_nh_metrics%ddxn_z_half_e(je,nlev,jb) ) )

          !Get net shear stress in the direction of vn at surface

          !shear stress in normal direction from cell 1
          stress_c1n = prm_diag%umfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_cell(je,jb,1)%v1     + &
                       prm_diag%vmfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_cell(je,jb,1)%v2

          !shear stress in normal direction from cell 2
          stress_c2n = prm_diag%umfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_cell(je,jb,2)%v1     + &
                       prm_diag%vmfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_cell(je,jb,2)%v2

          !Net stress at the edge
          flux_dn_e    = stress_c1n * p_int%c_lin_e(je,1,jb) + &
                         stress_c2n * p_int%c_lin_e(je,2,jb)

          vflux_up = km_ie(je,1,jb) * &
                     ( ( p_nh_prog%vn(je,nlev-1,jb) - p_nh_prog%vn(je,nlev,jb) ) * &
                       p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) )

          vert_metr = 0.5_wp * p_nh_metrics%ddxt_z_full(je,nlev,jb) *                    &
                      p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) *                       &
                      ( p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) *                     &
                        km_ie(je,nlev,jb) * p_nh_metrics%ddxn_z_half_e(je,nlev,jb) *     &
                        ( vt_ie(je,nlev-1,jb)-vt_ie(je,nlev+1,jb) ) -                    &
                        p_nh_metrics%inv_ddqz_z_half_e(je,nlev+1,jb) *                   &
                        km_ie(je,nlev+1,jb) * p_nh_metrics%ddxn_z_half_e(je,nlev+1,jb) * &
                        ( vt_ie(je,nlev,jb)-vt_ie(je,nlev+1,jb) ) )

          tot_tend(je,nlev,jb) = tot_tend(je,nlev,jb) +                                           &
                                 p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) *                     &
                                 ( flux_up_e - flux_dn_e + vflux_up +                             &
                                   2._wp * vflux_up * p_nh_metrics%ddxn_z_half_e(je,jk,jb) *      &
                                   p_nh_metrics%ddxn_z_full(je,jk,jb) +                           &
                                   vflux_up * p_nh_metrics%ddxt_z_half_e(je,jk,jb) *              &
                                   p_nh_metrics%ddxt_z_full(je,jk,jb) +                           &
                                   vert_metr ) * inv_rhoe(je,jk,jb)
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CASE(iimplicit)

      !vertical tendency- implicit solver
      !a*x(k-1)+b*x(k)+c*x(k+1)=rhs(k)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,a,b,c,rhs,var_new,&
!$OMP    flux_dn_e,stress_c1n,stress_c2n,dwdn)
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(dwdn)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 2 , nlev-1
#else
        DO jk = 2, nlev-1
          DO je = i_startidx, i_endidx
#endif
            a(je,jk)    = -km_ie(je,jk,jb) * ( z_4by3 *                                           &
                    p_nh_metrics%ddxn_z_half_e(je,jk,jb) * p_nh_metrics%ddxn_z_full(je,jk,jb) +   &
                    p_nh_metrics%ddxt_z_half_e(je,jk,jb) * p_nh_metrics%ddxt_z_full(je,jk,jb) +   &
                    1._wp ) * p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) *                          &
                    p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)

            c(je,jk)    = -km_ie(je,jk+1,jb) * (z_4by3 *                                          &
                    p_nh_metrics%ddxn_z_half_e(je,jk+1,jb) * p_nh_metrics%ddxn_z_full(je,jk,jb) + &
                    p_nh_metrics%ddxt_z_half_e(je,jk+1,jb) * p_nh_metrics%ddxt_z_full(je,jk,jb) + &
                    1._wp ) * p_nh_metrics%inv_ddqz_z_half_e(je,jk+1,jb) *                        &
                    p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)

             b(je,jk)   =  inv_dt - a(je,jk) - c(je,jk)

             !term due to dwdn- goes to RHS
             dwdn       = ( km_ie(je,jk,jb) * p_patch%edges%inv_dual_edge_length(je,jb) *         &
                              ( p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2))  -                &
                                p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1)) ) -               &
                              km_ie(je,jk+1,jb) * p_patch%edges%inv_dual_edge_length(je,jb) *     &
                              ( p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2))  -              &
                                p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) ) ) *           &
                          p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)

             rhs(je,jk) = p_nh_prog%vn(je,jk,jb) * inv_dt + dwdn

          END DO
        END DO
        !$ACC END PARALLEL

        !TOP Boundary treatment
        ! jk = 1
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(dwdn)
        DO je = i_startidx, i_endidx
          c(je,1)   = -km_ie(je,2,jb) * ( z_4by3 *                                                &
                   p_nh_metrics%ddxn_z_half_e(je,2,jb) * p_nh_metrics%ddxn_z_full(je,1,jb) +      &
                   p_nh_metrics%ddxt_z_half_e(je,2,jb)*p_nh_metrics%ddxt_z_full(je,1,jb)   +      &
                   1._wp ) * p_nh_metrics%inv_ddqz_z_half_e(je,2,jb) *                            &
                   p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) * inv_rhoe(je,1,jb)

          b(je,1)   =  inv_dt - c(je,1)

          !term due to dwdn- goes to RHS
          dwdn      =  -km_ie(je,2,jb) * p_patch%edges%inv_dual_edge_length(je,jb) *              &
                        ( p_nh_prog%w(iecidx(je,jb,2),2,iecblk(je,jb,2))  -                       &
                          p_nh_prog%w(iecidx(je,jb,1),2,iecblk(je,jb,1)) ) *                      &
                        p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) * inv_rhoe(je,1,jb)

          rhs(je,1) =  p_nh_prog%vn(je,1,jb) * inv_dt + dwdn
        END DO
        !$ACC END PARALLEL

        !SFC Boundary treatment
        !jk = nlev
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(dwdn, stress_c1n, stress_c2n, flux_dn_e)
        DO je = i_startidx, i_endidx
          a(je,nlev) = -km_ie(je,nlev,jb) * ( z_4by3 *                                            &
                  p_nh_metrics%ddxn_z_half_e(je,nlev,jb) * p_nh_metrics%ddxn_z_full(je,nlev,jb) + &
                  p_nh_metrics%ddxt_z_half_e(je,nlev,jb) * p_nh_metrics%ddxt_z_full(je,nlev,jb) + &
                  1._wp ) * p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) *                          &
                  p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)

          b(je,nlev) = inv_dt - a(je,nlev)

          !term due to dwdn- goes to RHS
          dwdn       =  km_ie(je,nlev,jb) * p_patch%edges%inv_dual_edge_length(je,jb) *           &
                            ( p_nh_prog%w(iecidx(je,jb,2),nlev,iecblk(je,jb,2))  -                &
                             p_nh_prog%w(iecidx(je,jb,1),nlev,iecblk(je,jb,1)) ) *                &
                            p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)

          !Get net shear stress in the direction of vn at surface

          !shear stress in normal direction from cell 1
          stress_c1n = prm_diag%umfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) *   &
                       p_patch%edges%primal_normal_cell(je,jb,1)%v1     +   &
                       prm_diag%vmfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) *   &
                       p_patch%edges%primal_normal_cell(je,jb,1)%v2

          !shear stress in normal direction from cell 2
          stress_c2n = prm_diag%umfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) *   &
                       p_patch%edges%primal_normal_cell(je,jb,2)%v1     +   &
                       prm_diag%vmfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) *   &
                       p_patch%edges%primal_normal_cell(je,jb,2)%v2


          !Net stress at the edge
          flux_dn_e    = stress_c1n * p_int%c_lin_e(je,1,jb) + stress_c2n * p_int%c_lin_e(je,2,jb)

          rhs(je,nlev) = p_nh_prog%vn(je,nlev,jb) * inv_dt + dwdn - flux_dn_e *                   &
                         p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)
        END DO
        !$ACC END PARALLEL

        !$ACC WAIT
        CALL tdma_solver_vec(a,b,c,rhs,1,nlev,i_startidx,i_endidx,var_new)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            tot_tend(je,jk,jb) = tot_tend(je,jk,jb) +                                  &
                                     ( var_new(je,jk) - p_nh_prog%vn(je,jk,jb) ) * inv_dt
          END DO
        END DO
        !$ACC END PARALLEL

      END DO!block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END SELECT !vert_scheme

   ! 4) Update vn: it makes more sense to first apply diffusion on vn
   ! and then get ddt_u/v than to interpolating tot_tend directly to
   ! ddt_u/v. Proof: during the test phase it was found that the latter slowed
   ! down the computation by 15%, although the results look same.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
       DO je = i_startidx, i_endidx
#endif
        vn_new(je,jk,jb) = p_nh_prog%vn(je,jk,jb) + dt * tot_tend(je,jk,jb)
       END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !5) Get turbulent tendency at cell center
    !$ACC WAIT
    CALL sync_patch_array(SYNC_E, p_patch, vn_new)
    CALL rbf_vec_interpol_cell(vn_new, p_patch, p_int, unew, vnew, opt_rlend=min_rlcell_int, opt_acc_async=.TRUE.)

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          ddt_u(jc,jk,jb) = ( unew(jc,jk,jb) - p_nh_diag%u(jc,jk,jb) ) * inv_dt
          ddt_v(jc,jk,jb) = ( vnew(jc,jk,jb) - p_nh_diag%v(jc,jk,jb) ) * inv_dt
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !subgrid fluxes: using vn_new, unew, vnew to store sgs_fluxes

    IF(is_sampling_time)THEN
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'is_sampling_time: OpenACC version currently not implemented')
#endif

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
      CALL init(unew(:,:,:), lacc=.TRUE.)
      CALL init(vnew(:,:,:), lacc=.TRUE.)
      CALL init(vn_new(:,:,:), lacc=.TRUE.)
!$OMP BARRIER

      rl_start   = grf_bdywidth_e+1
      rl_end     = min_rledge_int
      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,stress_c1n,stress_c2n)
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO je = i_startidx, i_endidx
#endif
            vn_new(je,jk-1,jb) = -km_ie(je,jk,jb) *                                               &
                                  ( ( p_nh_prog%vn(je,jk-1,jb) - p_nh_prog%vn(je,jk,jb) ) *       &
                                    p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) +                    &
                                    ( p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -           &
                                      p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1)) ) *         &
                                    p_patch%edges%inv_dual_edge_length(je,jb) )
          END DO
        END DO

        !-----------------------------------------------------------------
        ! jk = nlev+1
        !-----------------------------------------------------------------
        DO je = i_startidx, i_endidx

           !Get net shear stress in the direction of vn at surface

           !shear stress in normal direction from cell 1
           stress_c1n = prm_diag%umfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_cell(je,jb,1)%v1     + &
                        prm_diag%vmfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_cell(je,jb,1)%v2

           !shear stress in normal direction from cell 2
           stress_c2n = prm_diag%umfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_cell(je,jb,2)%v1     + &
                        prm_diag%vmfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_cell(je,jb,2)%v2

           !Net stress at the edge
           vn_new(je,nlev,jb)  = -stress_c1n * p_int%c_lin_e(je,1,jb) - &
                                  stress_c2n * p_int%c_lin_e(je,2,jb)
         END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      !Get sgs flux at cell center
      CALL sync_patch_array(SYNC_E, p_patch, vn_new)
      CALL rbf_vec_interpol_cell(vn_new, p_patch, p_int, unew, vnew, &
                                 opt_rlend=min_rlcell_int)

      !u sgs flux
      CALL levels_horizontal_mean(unew, p_patch%cells%area, p_patch%cells%owned, outvar)
      prm_diag%turb_diag_1dvar(1,idx_sgs_u_flx) = 0._wp
      DO jk = 2, nlevp1
        prm_diag%turb_diag_1dvar(jk,idx_sgs_u_flx) =  &
              prm_diag%turb_diag_1dvar(jk,idx_sgs_u_flx)+outvar(jk-1)
      END DO

      !v sgs flux
      CALL levels_horizontal_mean(vnew, p_patch%cells%area, p_patch%cells%owned, outvar)
      prm_diag%turb_diag_1dvar(1,idx_sgs_v_flx) = 0._wp
      DO jk = 2, nlevp1
        prm_diag%turb_diag_1dvar(jk,idx_sgs_v_flx) =  &
              prm_diag%turb_diag_1dvar(jk,idx_sgs_v_flx)+outvar(jk-1)
      END DO

    END IF !is_sampling_time

    !$ACC WAIT
    !$ACC END DATA ! local CREATE

  END SUBROUTINE diffuse_hori_velocity
  !-------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_vert_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for vertical velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !! - only solves for jk=2 to nlev. The bottom and top boundaries are left untouched
  !!------------------------------------------------------------------------
  SUBROUTINE diffuse_vert_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, &
                                   p_patch, p_int, km_ic, ddt_w, dt)

    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog    !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    REAL(wp),          INTENT(in)        :: km_ic(:,:,:)
    REAL(wp),  TARGET, INTENT(inout)     :: ddt_w(:,:,:) !< w tendency
    REAL(wp),          INTENT(in)        :: dt

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end, jg
    INTEGER :: jk, jb, je, jc, jcn, jbn, jvn, jcb
    INTEGER :: nlev

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk

    REAL(wp) :: flux_up_c, flux_dn_c, dvn1, dvn2, dvt1, dvt2, flux_up_v, flux_dn_v
    REAL(wp) :: norm_metr, tang_metr, w2mw1, w2mw1_p1

    REAL(wp) :: vt_e(nproma,p_patch%nlev,p_patch%nblks_e), inv_dt
    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs, var_new
    !interface level variables but only nlev quantities are needed
    REAL(wp) :: hor_tend(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp), POINTER :: tot_tend(:,:,:)
    REAL(wp) :: inv_rho_ic(nproma,2:p_patch%nlev,p_patch%nblks_c)!not necessary to allocate for nlev+1

    IF (msg_level >= 18) &
         CALL message(inmodule, 'diffuse_vert_velocity')

    !patch id
    jg = p_patch%id

    ! number of vertical levels
    nlev     = p_patch%nlev

    inv_dt  = 1._wp / dt

    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    tot_tend => ddt_w

    !$ACC DATA CREATE(vt_e, a, b, c, rhs, var_new, hor_tend, inv_rho_ic)

    !Some initializations
!$OMP PARALLEL
    CALL init(a(:,:), lacc=.TRUE.)
    CALL init(c(:,:), lacc=.TRUE.)
    IF(p_test_run)THEN
      CALL init(tot_tend, lacc=.TRUE.)
      CALL init(hor_tend, lacc=.TRUE.)
    END IF
!$OMP END PARALLEL


    CALL rbf_vec_interpol_edge(p_nh_prog%vn, p_patch, p_int, vt_e, opt_rlend=min_rledge_int-1, &
                               opt_acc_async=.TRUE.)

    !Calculate rho at interface for vertical diffusion
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          inv_rho_ic(jc,jk,jb) = 1._wp / rho_ic(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! 1) First get the horizontal tendencies at the half level edges

    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,jcn,jcb,dvn1,dvn2,flux_up_c,flux_dn_c,&
!$OMP jvn,jbn,dvt1,dvt2,flux_up_v,flux_dn_v,norm_metr,tang_metr,w2mw1,w2mw1_p1)
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4) PRIVATE(jcn, jcb, jvn, jbn, dvn1, dvn2, dvt1, dvt2) &
      !$ACC   PRIVATE(flux_up_c, flux_dn_c, flux_up_v, flux_dn_v, w2mw1, w2mw1_p1, norm_metr, tang_metr)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
#endif
          !tendency in normal direction
          !flux = visc_c * D_31_c where D_31(=D_13) is calculated at half level
          !cell center

          jcn     = iecidx(je,jb,2)
          jcb     = iecblk(je,jb,2)

          dvn2  = p_nh_diag%u(jcn,jk-1,jcb)*p_patch%edges%primal_normal_cell(je,jb,2)%v1 +  &
                  p_nh_diag%v(jcn,jk-1,jcb)*p_patch%edges%primal_normal_cell(je,jb,2)%v2 -  &
                  p_nh_diag%u(jcn,jk,jcb)*p_patch%edges%primal_normal_cell(je,jb,2)%v1   -  &
                  p_nh_diag%v(jcn,jk,jcb)*p_patch%edges%primal_normal_cell(je,jb,2)%v2

          flux_up_c = km_ic(jcn,jk,jcb) * (                                                 &
                      dvn2 * p_nh_metrics%inv_ddqz_z_half(jcn,jk,jcb) +                     &
                      ( w_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) - w_ie(je,jk,jb) ) *       &
                      p_patch%edges%inv_dual_edge_length(je,jb) -                           &
                      0.5_wp * ( p_nh_prog%w(jcn,jk-1,jcb) - p_nh_prog%w(jcn,jk+1,jcb) ) *  &
                      p_nh_metrics%inv_ddqz_z_half(jcn,jk,jcb) *                            &
                      p_nh_metrics%ddxn_z_half_c(jcn,jk,jcb) )


          jcn     = iecidx(je,jb,1)
          jcb     = iecblk(je,jb,1)

          dvn1  = p_nh_diag%u(jcn,jk-1,jcb)*p_patch%edges%primal_normal_cell(je,jb,1)%v1 +  &
                  p_nh_diag%v(jcn,jk-1,jcb)*p_patch%edges%primal_normal_cell(je,jb,1)%v2 -  &
                  p_nh_diag%u(jcn,jk,jcb)*p_patch%edges%primal_normal_cell(je,jb,1)%v1   -  &
                  p_nh_diag%v(jcn,jk,jcb)*p_patch%edges%primal_normal_cell(je,jb,1)%v2


          flux_dn_c = km_ic(jcn,jk,jcb) * (                                                 &
                      dvn1 * p_nh_metrics%inv_ddqz_z_half(jcn,jk,jcb) +                     &
                      ( w_ie(je,jk,jb) - w_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) ) *       &
                      p_patch%edges%inv_vert_vert_length(je,jb) * 2.0_wp -                  &
                      0.5_wp * ( p_nh_prog%w(jcn,jk-1,jcb) - p_nh_prog%w(jcn,jk+1,jcb) ) *  &
                      p_nh_metrics%inv_ddqz_z_half(jcn,jk,jcb) * &
                      p_nh_metrics%ddxn_z_half_c(jcn,jk,jcb) )

          ! km_c should be replaced with km_e
          jcn = iecidx(je,jb,2)
          jcb = iecblk(je,jb,2)
          w2mw1    = p_patch%edges%inv_vert_vert_length(je,jb) * km_c(jcn,jk-1,jcb) *       &
                     ( p_nh_prog%w(jcn,jk-1,jcb) + p_nh_prog%w(jcn,jk,jcb) -                &
                       p_nh_prog%w(iecidx(je,jb,1),jk-1,iecblk(je,jb,1))   -                &
                       p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1)) )
          w2mw1_p1 = p_patch%edges%inv_vert_vert_length(je,jb) * km_c(jcn,jk,jcb) *         &
                     ( p_nh_prog%w(jcn,jk,jcb) + p_nh_prog%w(jcn,jk+1,jcb) -                &
                       p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1)) -                    &
                       p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) )

          ! km_c should be replaced with km_e
          norm_metr = -p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) *                           &
                       p_nh_metrics%ddxn_z_half_e(je,jk,jb) * ( km_c(jcn,jk-1,jcb) *        &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk-1,jb) *                         &
                       ( vn_ie(je,jk-1,jb)-vn_ie(je,jk,jb) ) -                              &
                       km_c(jcn,jk,jcb) * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) *        &
                       ( vn_ie(je,jk,jb)-vn_ie(je,jk+1,jb) ) +                              &
                      p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) * ( w2mw1 - w2mw1_p1 ) )

          !tendency in tangential direction
          !flux = visc_v * D_32_v where D_32(=D_23) is calculated at half level
          !between vertex and edge center

          jvn     = ividx(je,jb,2)
          jbn     = ivblk(je,jb,2)

          dvt2    = ( u_vert(jvn,jk-1,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                      v_vert(jvn,jk-1,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v2 + &
                      vt_e(je,jk-1,jb) ) * 0.5_wp           - &
                    ( u_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v1 +   &
                      v_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v2 +   &
                      vt_e(je,jk,jb) ) * 0.5_wp

          flux_up_v = km_iv(jvn,jk,jbn) * (                                             &
                      dvt2 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) +               &
                      p_patch%edges%tangent_orientation(je,jb) *                        &
                      (w_vert(jvn,jk,jbn) - w_ie(je,jk,jb) ) /                          &
                      p_patch%edges%edge_vert_length(je,jb,2) -                         & !check orig. code
                      0.5_wp * ( w_vert(jvn,jk-1,jbn) - w_vert(jvn,jk+1,jbn) ) *        &
                      p_nh_metrics%inv_ddqz_z_half_v(jvn,jk-1,jbn) *                    &
                      p_nh_metrics%ddxt_z_half_v(jvn,jk,jbn) )

          jvn     = ividx(je,jb,1)
          jbn     = ivblk(je,jb,1)

          dvt1    = ( u_vert(jvn,jk-1,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                      v_vert(jvn,jk-1,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v2 + &
                      vt_e(je,jk-1,jb) ) * 0.5_wp           - &
                    ( u_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v1 +   &
                      v_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v2 +   &
                      vt_e(je,jk,jb) ) * 0.5_wp

          flux_dn_v = km_iv(jvn,jk,jbn) * (                                             &
                      dvt1 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) +               &
                      p_patch%edges%tangent_orientation(je,jb) *                        &
                      ( w_ie(je,jk,jb) - w_vert(jvn,jk,jbn) ) /                         &
                      p_patch%edges%edge_vert_length(je,jb,1) -                         & !check orig. code
                      0.5_wp * ( w_vert(jvn,jk-1,jbn) - w_vert(jvn,jk+1,jbn) ) *        &
                      p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) *                      &
                      p_nh_metrics%ddxt_z_half_v(jvn,jk,jbn) )

          ! km_c should be replaced with km_iv
          jvn = ividx(je,jb,2)
          jbn = ivblk(je,jb,2)
          w2mw1 = 0.5_wp * p_patch%edges%inv_primal_edge_length(je,jb) *                &
                  p_patch%edges%tangent_orientation(je,jb) *                            &
                  km_iv(jvn,jk-1,jbn) * ( w_vert(jvn,jk-1,jbn) +                        &
                  w_vert(jvn,jk,jbn) - w_vert(ividx(je,jb,1),jk-1,ivblk(je,jb,1)) -     &
                  w_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))  )
          w2mw1_p1 = 0.5_wp * p_patch%edges%inv_primal_edge_length(je,jb) *             & 
                     p_patch%edges%tangent_orientation(je,jb) *                         &
                     km_iv(jvn,jk,jbn) * ( w_vert(jvn,jk,jbn) +                         &
                     w_vert(jvn,jk+1,jbn) - w_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) -  &
                     w_vert(ividx(je,jb,1),jk+1,ivblk(je,jb,1)) )

          ! km_c should be replaced with km_e
          tang_metr = -p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) *                           &
                       p_nh_metrics%ddxt_z_half_e(je,jk,jb) *                               &
                       ( km_c(jcn,jk-1,jcb) * p_nh_metrics%inv_ddqz_z_full_e(je,jk-1,jb) *  &
                         ( vt_ie(je,jk-1,jb) - vt_ie(je,jk,jb) ) -                          &
                         km_c(jcn,jk,jcb) * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) *      &
                         ( vt_ie(je,jk,jb) -vt_ie(je,jk+1,jb) ) -                           &
                         p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) * (w2mw1 - w2mw1_p1) )
                     

          hor_tend(je,jk,jb) = ( flux_up_c - flux_dn_c ) *                                  &
                               p_patch%edges%inv_dual_edge_length(je,jb) +                  &
                               p_patch%edges%tangent_orientation(je,jb) *                   &
                               ( flux_up_v - flux_dn_v ) *                                  &
                               p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp -        &
                               norm_metr - tang_metr
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    CALL sync_patch_array(SYNC_E, p_patch, hor_tend)

    !Interpolate horizontal tendencies to w point: except top and bottom boundaries
    !w==0 at these boundaries

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          tot_tend(jc,jk,jb) =  inv_rho_ic(jc,jk,jb) *                                            &
                         ( hor_tend(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
                           hor_tend(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
                           hor_tend(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb) )

        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !
    ! 2) Vertical tendency: evaluated at w point
    !
    SELECT CASE(les_config(jg)%vert_scheme_type)

    CASE(iexplicit)
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'vert_scheme_type=iexplicit: OpenACC version currently not implemented')
#endif

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,flux_up_c,flux_dn_c)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
#endif
            flux_up_c = km_c(jc,jk-1,jb) * ( ( p_nh_prog%w(jc,jk-1,jb) - p_nh_prog%w(jc,jk,jb) ) * &
                        p_nh_metrics%inv_ddqz_z_full(jc,jk-1,jb) - z_1by3 * div_c(jc,jk-1,jb) )

            flux_dn_c = km_c(jc,jk,jb) * ( ( p_nh_prog%w(jc,jk,jb) - p_nh_prog%w(jc,jk+1,jb) ) *   &
                        p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) - z_1by3 * div_c(jc,jk,jb) )

            tot_tend(jc,jk,jb) = tot_tend(jc,jk,jb) + 2._wp *                                        &
                                  ( flux_up_c - flux_dn_c + ( flux_up_c *                            &
                                    p_nh_metrics%ddxn_z_full_c(jc,jk-1,jb) - flux_dn_c *             &
                                    p_nh_metrics%ddxn_z_full_c(jc,jk,jb) ) *                         &
                                    p_nh_metrics%ddxn_z_half_c(jc,jk,jb) +                           &
                                  ( flux_up_c * p_nh_metrics%ddxt_z_full_c(jc,jk-1,jb) - flux_dn_c * &
                                    p_nh_metrics%ddxt_z_full_c(jc,jk,jb) ) *                         &
                                    p_nh_metrics%ddxt_z_half_c(jc,jk,jb) +                           &
                                    km_c(jc,jk-1,jb) * z_1by3 * div_c(jc,jk-1,jb) -                  &
                                    km_c(jc,jk,jb) * z_1by3 * div_c(jc,jk,jb) ) *                    &
                                    p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CASE(iimplicit)
    !vertical tendency- implicit solver
    !a*x(k-1)+b*x(k)+c*x(k+1)=rhs(k)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,a,b,c,rhs,var_new)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 3, nlev-1
#else
        DO jk = 3, nlev-1
          DO jc = i_startidx, i_endidx
#endif

            a(jc,jk)   = - km_c(jc,jk-1,jb) * &
                ( p_nh_metrics%ddxn_z_full_c(jc,jk-1,jb) * p_nh_metrics%ddxn_z_half_c(jc,jk,jb) + &
                  p_nh_metrics%ddxt_z_full_c(jc,jk-1,jb) * p_nh_metrics%ddxt_z_half_c(jc,jk,jb) + &
                  2._wp ) * p_nh_metrics%inv_ddqz_z_full(jc,jk-1,jb) *                            &
                p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

            c(jc,jk)   = - km_c(jc,jk,jb) *   &
                ( p_nh_metrics%ddxn_z_full_c(jc,jk,jb) * p_nh_metrics%ddxn_z_half_c(jc,jk,jb) +   &
                  p_nh_metrics%ddxt_z_full_c(jc,jk,jb) * p_nh_metrics%ddxt_z_half_c(jc,jk,jb) +   &
                  2._wp ) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) *                              &
                p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

             b(jc,jk)   =  inv_dt - a(jc,jk) - c(jc,jk)

             rhs(jc,jk) =  p_nh_prog%w(jc,jk,jb) * inv_dt +                                       &
                           2._wp * ( km_c(jc,jk,jb) * z_1by3 * div_c(jc,jk,jb) -                  &
                                     km_c(jc,jk-1,jb) * z_1by3 * div_c(jc,jk-1,jb) ) *            &
                           p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)
          END DO
        END DO
        !$ACC END PARALLEL

        !TOP Boundary treatment:w==0
        !jk = 2
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          c(jc,2)   = - km_c(jc,2,jb) *                                                           & 
                    ( p_nh_metrics%ddxn_z_full_c(jc,2,jb) * p_nh_metrics%ddxn_z_half_c(jc,2,jb) + &
                      p_nh_metrics%ddxt_z_full_c(jc,2,jb) * p_nh_metrics%ddxt_z_half_c(jc,2,jb) + &
                      2._wp ) * p_nh_metrics%inv_ddqz_z_full(jc,2,jb) *                           &
                    p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

          b(jc,2)   = inv_dt - c(jc,2) + 2._wp * km_c(jc,1,jb) *                                  &
                       p_nh_metrics%inv_ddqz_z_full(jc,1,jb) *                                    &
                       p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

          rhs(jc,2) = p_nh_prog%w(jc,2,jb) * inv_dt +  &
                       2._wp * ( km_c(jc,2,jb) * z_1by3 * div_c(jc,2,jb) -                        &
                                 km_c(jc,1,jb) * z_1by3 * div_c(jc,1,jb) ) *                      &
                               p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)
        END DO
        !$ACC END PARALLEL

         !SFC Boundary treatment:w==0
         !jk = nlev
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          a(jc,nlev)   = - km_c(jc,nlev-1,jb) * &
            ( p_nh_metrics%ddxn_z_full_c(jc,nlev-1,jb) * p_nh_metrics%ddxn_z_half_c(jc,nlev,jb) + &
              p_nh_metrics%ddxt_z_full_c(jc,nlev-1,jb) * p_nh_metrics%ddxt_z_half_c(jc,nlev,jb) + &
              2._wp ) * p_nh_metrics%inv_ddqz_z_full(jc,nlev-1,jb) *                              &
              p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)

          b(jc,nlev)  = inv_dt - a(jc,nlev) + 2._wp * km_c(jc,nlev,jb) *                          &
                        p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) *                                &
                        p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)

          rhs(jc,nlev)= p_nh_prog%w(jc,nlev,jb) * inv_dt +                                        &
                        2._wp * ( km_c(jc,nlev,jb) * z_1by3 * div_c(jc,nlev,jb) -                 &
                                  km_c(jc,nlev-1,jb) * z_1by3 * div_c(jc,nlev-1,jb) ) *           &
                                p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)
        END DO
        !$ACC END PARALLEL

        !$ACC WAIT
        CALL tdma_solver_vec(a,b,c,rhs,2,nlev,i_startidx,i_endidx,var_new)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
#endif
            tot_tend(jc,jk,jb) = tot_tend(jc,jk,jb) +                                  &
                                     ( var_new(jc,jk) - p_nh_prog%w(jc,jk,jb) ) * inv_dt
          END DO
        END DO
        !$ACC END PARALLEL
      END DO !block
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END SELECT !vert_scheme

    !Now update w

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          p_nh_prog%w(jc,jk,jb) = p_nh_prog%w(jc,jk,jb) + dt * tot_tend(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    CALL sync_patch_array(SYNC_C, p_patch, p_nh_prog%w)

    !$ACC WAIT
    !$ACC END DATA ! local CREATE

  END SUBROUTINE diffuse_vert_velocity
  !-------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_scalar
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for cell based scalars
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !!------------------------------------------------------------------------
  SUBROUTINE diffuse_scalar(var, p_nh_metrics, p_patch, p_int, tot_tend,  &
                            exner, prm_diag, rho, dt, scalar_name)

    REAL(wp),          INTENT(in)        :: var(:,:,:)   !input scalar
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch), INTENT(inout), TARGET :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    REAL(wp),        INTENT(inout)       :: tot_tend(:,:,:)!<total tendency
    REAL(wp),          INTENT(in)        :: exner(:,:,:)   !
    REAL(wp),          INTENT(in)        :: rho(:,:,:)     !density at cell center
    TYPE(t_nwp_phy_diag),  INTENT(inout) :: prm_diag       !< atm phys vars
    REAL(wp),          INTENT(in)        :: dt
    INTEGER, INTENT(in)                  :: scalar_name

    !Local variables
    REAL(wp) :: flux_up, flux_dn, inv_dt, norm_metr, tang_metr
    REAL(wp) :: nabla2_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: fac(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: exner_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: sgs_flux(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: exner_me(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: exner_ie(nproma,p_patch%nlev+1,p_patch%nblks_e)
    REAL(wp) :: sflux(nproma,p_patch%nblks_c), outvar(p_patch%nlev+1)

    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs
    REAL(wp)                                 :: var_new(nproma, p_patch%nlev)
    REAL(wp), POINTER :: kh_ic(:,:,:)

    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, &
      &                                    ieblk, ievidx, ievblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jc, jg, jbn, jvn, jcn
    INTEGER  :: nlev, nlevp1

    ! helper variables for metric terms
    REAL(wp) :: &
      var_ic        (nproma, p_patch%nlevp1, p_patch%nblks_c),  &
      var_ie        (nproma, p_patch%nlevp1, p_patch%nblks_e),  &
      var_iv        (nproma, p_patch%nlevp1, p_patch%nblks_v),  &
      metric_tend_e (nproma, p_patch%nlev,   p_patch%nblks_e)

    IF (msg_level >= 18) &
         CALL message(inmodule, 'diffuse_scalar')

    !patch id
    jg = p_patch%id

    ! number of vertical levels
    nlev = p_patch%nlev
    nlevp1 = nlev+1
    inv_dt  = 1._wp / dt

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ievidx => p_patch%edges%vertex_idx
    ievblk => p_patch%edges%vertex_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk
    kh_ic => prm_diag%tkvh

    !$ACC DATA CREATE(nabla2_e, fac, exner_ic, sgs_flux, exner_me, exner_ie, sflux, a, b, c, rhs) &
    !$ACC   CREATE(var_new, var_ic, var_ie, var_iv, metric_tend_e)

    !multiply by exner to convert from theta tend to temp tend
    !assuming that exner perturbation are small compared to temp

    !1) First set exner local vars to 1 for other scalars
    !   Soon get different routines for different scalars
!$OMP PARALLEL
    CALL init(exner_me(:,:,:),1._wp, lacc=.TRUE.)
    CALL init(exner_ie(:,:,:),1._wp, lacc=.TRUE.)
    CALL init(exner_ic(:,:,:),1._wp, lacc=.TRUE.)
    CALL init(a(:,:), lacc=.TRUE.)
    CALL init(c(:,:), lacc=.TRUE.)
    CALL init(var_ic(:,:,:), lacc=.TRUE.)
!$OMP END PARALLEL

    !2) Calculate exner at edge for horizontal diffusion

    IF (scalar_name == tracer_theta) &
      CALL cells2edges_scalar(exner, p_patch, p_int%c_lin_e, exner_me, opt_rlend=min_rledge_int-2, lacc=.TRUE.)

    !3) Calculate exner at interface for vertical diffusion

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
    IF(scalar_name == tracer_theta) THEN

      rl_start   = grf_bdywidth_e
      rl_end     = min_rledge_int-1
      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(je,jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO je = i_startidx, i_endidx
#endif
            exner_ie(je,jk,jb) = ( p_nh_metrics%wgtfac_e(je,jk,jb) * exner_me(je,jk,jb) + &
                                 ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) *            &
                                 exner_me(je,jk-1,jb) )
          END DO
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO je = i_startidx, i_endidx
          exner_ie(je,nlev+1,jb) = p_nh_metrics%wgtfacq_e(je,1,jb) * exner_me(je,nlev,jb)   +     &
                                   p_nh_metrics%wgtfacq_e(je,2,jb) * exner_me(je,nlev-1,jb) +     &
                                   p_nh_metrics%wgtfacq_e(je,3,jb) * exner_me(je,nlev-2,jb)

          exner_ie(je,1,jb)      = p_nh_metrics%wgtfacq1_e(je,1,jb) * exner_me(je,1,jb) +         &
                                   p_nh_metrics%wgtfacq1_e(je,2,jb) * exner_me(je,2,jb) +         &
                                   p_nh_metrics%wgtfacq1_e(je,3,jb) * exner_me(je,3,jb)
        END DO
        !$ACC END PARALLEL
      END DO
!$OMP END DO NOWAIT


      rl_start   = grf_bdywidth_c+1
      rl_end     = min_rlcell_int
      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
#endif
            exner_ic(jc,jk,jb) = ( p_nh_metrics%wgtfac_c(jc,jk,jb) * exner(jc,jk,jb) +            &
                                 ( 1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb)) * exner(jc,jk-1,jb) )
          END DO
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
            exner_ic(jc,nlev+1,jb) = p_nh_metrics%wgtfacq_c(jc,1,jb) * exner(jc,nlev,jb)   +      &
                                     p_nh_metrics%wgtfacq_c(jc,2,jb) * exner(jc,nlev-1,jb) +      &
                                     p_nh_metrics%wgtfacq_c(jc,3,jb) * exner(jc,nlev-2,jb)
        ENDDO
        !$ACC END PARALLEL
      END DO
!$OMP END DO

    END IF ! var == theta

    !4) Calculate var at interface for vertical diffusion

    rl_start   = grf_bdywidth_c
    rl_end     = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          var_ic(jc,jk,jb) = ( p_nh_metrics%wgtfac_c(jc,jk,jb) * var(jc,jk,jb) +                  &
                             ( 1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb) ) * var(jc,jk-1,jb) )
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        var_ic(jc,nlev+1,jb) = p_nh_metrics%wgtfacq_c(jc,1,jb) * var(jc,nlev,jb)   +              &
                               p_nh_metrics%wgtfacq_c(jc,2,jb) * var(jc,nlev-1,jb) +              &
                               p_nh_metrics%wgtfacq_c(jc,3,jb) * var(jc,nlev-2,jb)
        var_ic(jc,1,jb)      = p_nh_metrics%wgtfacq1_c(jc,1,jb) * var(jc,1,jb) +                  &
                               p_nh_metrics%wgtfacq1_c(jc,2,jb) * var(jc,2,jb) +                  &
                               p_nh_metrics%wgtfacq1_c(jc,3,jb) * var(jc,3,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !$ACC WAIT
    CALL sync_patch_array(SYNC_C, p_patch, var_ic)

    CALL cells2edges_scalar(var_ic, p_patch, p_int%c_lin_e, var_ie, lacc=.TRUE.)
    CALL cells2verts_scalar(var_ic, p_patch, p_int%cells_aw_verts, var_iv, opt_acc_async=.TRUE.)


!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !Special boundary treatment for different scalars

    IF (scalar_name ==tracer_theta) THEN
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
#endif
            fac(jc,jk,jb) = cpd * rcvd / rho(jc,jk,jb)
          END DO
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          sflux(jc,jb) = prm_diag%shfl_s(jc,jb) * rcpd
        END DO
        !$ACC END PARALLEL
      END DO
!$OMP END DO
    ELSEIF (scalar_name == tracer_qv) THEN
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
#endif
            fac(jc,jk,jb) = 1._wp / rho(jc,jk,jb)
          END DO
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          sflux(jc,jb) = prm_diag%lhfl_s(jc,jb) / alv
        END DO
        !$ACC END PARALLEL
      END DO
!$OMP END DO
    ELSEIF (scalar_name == tracer_qc) THEN
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
#endif
            fac(jc,jk,jb) = 1._wp / rho(jc,jk,jb)
          END DO
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
            sflux(jc,jb) = 0._wp
        END DO
        !$ACC END PARALLEL
      END DO
!$OMP END DO
    END IF
!$OMP END PARALLEL


    ! compute metric contribution on main level edges
    !
!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
        rl_start   = grf_bdywidth_e
        rl_end     = min_rledge_int-1
        i_startblk = p_patch%edges%start_block(rl_start)
        i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(je,jb,jk,i_startidx,i_endidx,norm_metr,tang_metr,jcn,jvn,jbn)
        DO jb = i_startblk,i_endblk
          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(jcn, jbn, jvn, norm_metr, tang_metr)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = 1, nlev
#else
          DO jk = 1, nlev-1
            DO je = i_startidx, i_endidx
#endif
              jcn     = iecidx(je,jb,2)
              jbn     = iecblk(je,jb,2)

              norm_metr = - les_config(jg)%rturb_prandtl *                                        &
                          p_patch%edges%inv_dual_edge_length(je,jb) * ( km_ie(je,jk,jb) *         &
                          ( var_ic(jcn,jk,jbn) - var_ic(iecidx(je,jb,1),jk,iecblk(je,jb,1)) ) *   &
                          exner_ie(je,jk,jb) - km_ie(je,jk+1,jb) *                                &
                          ( var_ic(jcn,jk+1,jbn) -                                                &
                            var_ic(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) ) * exner_ie(je,jk+1,jb) )

              jvn     = ievidx(je,jb,2)
              jbn     = ievblk(je,jb,2)

              tang_metr = - les_config(jg)%rturb_prandtl *                                        &
                          p_patch%edges%inv_primal_edge_length(je,jb) *                           &
                          p_patch%edges%tangent_orientation(je,jb) * ( km_ie(je,jk,jb) *          & 
                          ( var_iv(jvn,jk,jbn) - var_iv(ievidx(je,jb,1),jk,ievblk(je,jb,1)) ) *   &
                            exner_ie(je,jk,jb) - km_ie(je,jk+1,jb) *                              &
                          ( var_iv(jvn,jk+1,jbn) -                                                &
                            var_iv(ievidx(je,jb,1),jk+1,ievblk(je,jb,1)) ) * exner_ie(je,jk+1,jb) )

              metric_tend_e(je,jk,jb) = ( norm_metr * p_nh_metrics%ddxn_z_full(je,jk,jb) +        &
                                          tang_metr * p_nh_metrics%ddxt_z_full(je,jk,jb) ) *      &
                                          p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

            ENDDO
          ENDDO
          !$ACC END PARALLEL

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO je = i_startidx, i_endidx
            metric_tend_e(je,nlev,jb) = 0._wp
          END DO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO

    ! use conservative discretization div(k*grad(var))-horizontal part
    ! following mo_nh_diffusion

    !include halo points and boundary points because these values will be
    !used in next loop
    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

    IF ( atm_phy_nwp_config(jg)%inwp_turb == iprog) THEN 
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'inwp_turb=iprog: OpenACC version currently not implemented')
#endif
!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! compute kh_smag_e * grad(var)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = 1, nlev
#else
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
#endif
                nabla2_e(je,jk,jb) = 0.5_wp * ( kh_ie(je,jk,jb) + kh_ie(je,jk+1,jb) ) *           &
                                     exner_me(je,jk,jb) * (                                       &
                                     p_patch%edges%inv_dual_edge_length(je,jb) *                  &
                                     ( var(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -                  &
                                       var(iecidx(je,jb,1),jk,iecblk(je,jb,1)) ) -                &
                                     ( var_ie(je,jk,jb) - var_ie(je,jk+1,jb) ) *                  &
                                     p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) *                   &
                                     p_nh_metrics%ddxn_z_full(je,jk,jb)  )
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
    ELSE
!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! compute kh_smag_e * grad(var)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = 1, nlev
#else
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
#endif
                nabla2_e(je,jk,jb) = 0.5_wp * ( km_ie(je,jk,jb) + km_ie(je,jk+1,jb) ) *           &
                                     les_config(jg)%rturb_prandtl * exner_me(je,jk,jb) * (        &
                                     p_patch%edges%inv_dual_edge_length(je,jb) *                  &
                                     ( var(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -                  &
                                       var(iecidx(je,jb,1),jk,iecblk(je,jb,1)) ) -                &
                                     ( var_ie(je,jk,jb) - var_ie(je,jk+1,jb) ) *                  &
                                     p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) *                   &
                                     p_nh_metrics%ddxn_z_full(je,jk,jb)  )
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO
    ENDIF

        ! now compute the divergence of the quantity above
        ! and add metric-terms contribution
        rl_start   = grf_bdywidth_c+1
        rl_end     = min_rlcell_int
        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1, nlev
#else
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
#endif
              !horizontal tendency
              tot_tend(jc,jk,jb) = fac(jc,jk,jb) * (                                              &
                    nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) +      &
                    nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) +      &
                    nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) +      &
                    metric_tend_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +  &
                    metric_tend_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +  &
                    metric_tend_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO
!$OMP END PARALLEL


       !---------------------------------------------------------------
       !Vertical diffusion
       !---------------------------------------------------------------

       rl_start   = grf_bdywidth_c+1
       rl_end     = min_rlcell_int
       i_startblk = p_patch%cells%start_block(rl_start)
       i_endblk   = p_patch%cells%end_block(rl_end)

       SELECT CASE(les_config(jg)%vert_scheme_type)

       CASE(iexplicit)
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'vert_scheme_type=iexplicit: OpenACC version currently not implemented')
#endif

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,flux_up,flux_dn)
       DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 2, nlev-1
#else
          DO jk = 2, nlev-1
            DO jc = i_startidx, i_endidx
#endif
              flux_up = kh_ic(jc,jk,jb) * ( var(jc,jk-1,jb) - var(jc,jk,jb) ) *           &
                        p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * exner_ic(jc,jk,jb)

              flux_dn = kh_ic(jc,jk+1,jb) * ( var(jc,jk,jb) - var(jc,jk+1,jb) ) *         &
                        p_nh_metrics%inv_ddqz_z_half(jc,jk+1,jb) * exner_ic(jc,jk+1,jb)

              tot_tend(jc,jk,jb) = tot_tend(jc,jk,jb) + &
                             p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * fac(jc,jk,jb) *             &
                             ( flux_up - flux_dn + ( flux_up *                                    &
                                p_nh_metrics%ddxn_z_half_c(jc,jk,jb) -                            &
                               flux_dn * p_nh_metrics%ddxn_z_half_c(jc,jk+1,jb) ) *               &
                               p_nh_metrics%ddxn_z_full_c(jc,jk,jb) +                             &
                               ( flux_up * p_nh_metrics%ddxt_z_half_c(jc,jk,jb) -                 &
                                 flux_dn  *p_nh_metrics%ddxt_z_half_c(jc,jk+1,jb) ) *             &
                               p_nh_metrics%ddxt_z_full_c(jc,jk,jb) )
            ENDDO
          ENDDO

          !Boundary treatment

          !--------------------------------------------------------
          !jk = 1
          !--------------------------------------------------------
          DO jc = i_startidx, i_endidx
            flux_up = 0._wp

            flux_dn = kh_ic(jc,2,jb) * (var(jc,1,jb) - var(jc,2,jb)) *                &
                      p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * exner_ic(jc,2,jb)

            tot_tend(jc,1,jb) = tot_tend(jc,1,jb) + &
                                p_nh_metrics%inv_ddqz_z_full(jc,1,jb) * fac(jc,1,jb) *            &
                                ( flux_up - flux_dn + ( flux_up *                                 &
                                  p_nh_metrics%ddxn_z_half_c(jc,1,jb) -                           &
                                  flux_dn*p_nh_metrics%ddxn_z_half_c(jc,2,jb) ) *                 &
                                p_nh_metrics%ddxn_z_full_c(jc,1,jb) )
          ENDDO
          !--------------------------------------------------------
          !jk = nlev
          !--------------------------------------------------------
          DO jc = i_startidx, i_endidx
            flux_up = kh_ic(jc,nlev,jb) * ( var(jc,nlev-1,jb) - var(jc,nlev,jb) ) *   &
                      p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * exner_ic(jc,nlev,jb)

            flux_dn  = - sflux(jc,jb) * exner_ic(jc,nlevp1,jb)

            tot_tend(jc,nlev,jb) = tot_tend(jc,nlev,jb) + &
                                    p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) * fac(jc,nlev,jb) *  &
                                    ( flux_up - flux_dn + ( flux_up *                             &
                                      p_nh_metrics%ddxn_z_half_c(jc,nlev,jb) -                    &
                                      flux_dn * p_nh_metrics%ddxn_z_half_c(jc,nlev+1,jb) ) *      &
                                      p_nh_metrics%ddxn_z_full_c(jc,nlev,jb) )
          ENDDO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CASE(iimplicit)

        !vertical tendency- implicit solver
        !a*x(k-1)+b*x(k)+c*x(k+1)=rhs(k)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,a,b,c,rhs,var_new)
          DO jb = i_startblk,i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                               i_startidx, i_endidx, rl_start, rl_end)
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
            DO jc = i_startidx, i_endidx
              DO jk = 2, nlev-1
#else
            DO jk = 2, nlev-1
              DO jc = i_startidx, i_endidx
#endif
                a(jc,jk)   = - kh_ic(jc,jk,jb) * p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) *         &
                  p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) *                                        &
                  ( p_nh_metrics%ddxn_z_half_c(jc,jk,jb) * p_nh_metrics%ddxn_z_full_c(jc,jk,jb) + &
                    p_nh_metrics%ddxt_z_half_c(jc,jk,jb) * p_nh_metrics%ddxt_z_full_c(jc,jk,jb) + & 
                    1._wp ) * exner_ic(jc,jk,jb) * fac(jc,jk,jb)

                c(jc,jk)   = -kh_ic(jc,jk+1,jb) * p_nh_metrics%inv_ddqz_z_half(jc,jk+1,jb) *      &
                  p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) *                                        &
                  ( p_nh_metrics%ddxn_z_half_c(jc,jk+1,jb) *                                      &
                    p_nh_metrics%ddxn_z_full_c(jc,jk,jb) +                                        &
                    p_nh_metrics%ddxt_z_half_c(jc,jk+1,jb) *                                      &
                    p_nh_metrics%ddxt_z_full_c(jc,jk,jb) + 1._wp ) * exner_ic(jc,jk+1,jb) *       &
                  fac(jc,jk,jb)

                b(jc,jk)   =  inv_dt - a(jc,jk) - c(jc,jk)

                rhs(jc,jk) =  var(jc,jk,jb) * inv_dt

              END DO
            END DO
            !$ACC END PARALLEL

            !TOP Boundary treatment
            ! jk = 1
            !
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx, i_endidx
              c(jc,1)   = - kh_ic(jc,2,jb) * p_nh_metrics%inv_ddqz_z_half(jc,2,jb) *              &
               p_nh_metrics%inv_ddqz_z_full(jc,1,jb) *                                            &
               ( p_nh_metrics%ddxn_z_half_c(jc,2,jb) * p_nh_metrics%ddxn_z_full_c(jc,1,jb) +      &
                 p_nh_metrics%ddxt_z_half_c(jc,2,jb) * p_nh_metrics%ddxt_z_full_c(jc,1,jb) +      &
                 1._wp ) * exner_ic(jc,2,jb) * fac(jc,1,jb)

              b(jc,1)   = inv_dt - c(jc,1)

              rhs(jc,1) = var(jc,1,jb) * inv_dt
            END DO
            !$ACC END PARALLEL

            !SFC Boundary treatment
            !jk = nlev
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx, i_endidx
              a(jc,nlev)   = - kh_ic(jc,nlev,jb) * p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) *     &
                              p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) *                          &
                              ( p_nh_metrics%ddxn_z_half_c(jc,nlev,jb) *                          &
                                p_nh_metrics%ddxn_z_full_c(jc,nlev,jb) +                          &
                                p_nh_metrics%ddxt_z_half_c(jc,nlev,jb) *                          &
                                p_nh_metrics%ddxt_z_full_c(jc,nlev,jb) + 1._wp ) *                &
                                exner_ic(jc,nlev,jb) * fac(jc,nlev,jb)

              b(jc,nlev)   = inv_dt - a(jc,nlev)

              rhs(jc,nlev) = var(jc,nlev,jb) * inv_dt + sflux(jc,jb) * exner_ic(jc,nlevp1,jb) *   &
                             p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) * fac(jc,nlev,jb)
            END DO
            !$ACC END PARALLEL

            !$ACC WAIT
            !CALL TDMA
            CALL tdma_solver_vec(a,b,c,rhs,1,nlev,i_startidx,i_endidx,var_new)
            !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                tot_tend(jc,jk,jb) = tot_tend(jc,jk,jb) + ( var_new(jc,jk) - var(jc,jk,jb) ) * inv_dt
              END DO
            END DO
            !$ACC END PARALLEL LOOP
          END DO!block
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END SELECT !vert_scheme


    !subgrid fluxes for different scalars

    IF(is_sampling_time)THEN
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'is_sampling_time: OpenACC version currently not implemented')
#endif

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
#endif
            sgs_flux(jc,jk,jb) = -kh_ic(jc,jk,jb) * (var(jc,jk-1,jb) - var(jc,jk,jb)) * &
                                   p_nh_metrics%inv_ddqz_z_half(jc,jk,jb)
          END DO
        END DO
        DO jc = i_startidx, i_endidx
          sgs_flux(jc,1,jb) = 0._wp
          !extrapolate for surface density
          sgs_flux(jc,nlevp1,jb) = sflux(jc,jb)
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (scalar_name == tracer_theta) THEN
        CALL levels_horizontal_mean(sgs_flux, p_patch%cells%area, p_patch%cells%owned, &
                                    outvar)

        outvar = outvar * cpd
        prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_th_flx) =  &
               prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_th_flx)+outvar(1:nlevp1)

      ELSEIF (scalar_name == tracer_qv) THEN
        CALL levels_horizontal_mean(sgs_flux, p_patch%cells%area, p_patch%cells%owned, &
                                     outvar)

        outvar = outvar * alv
        prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_qv_flx) =  &
               prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_qv_flx)+outvar(1:nlevp1)

      ELSEIF (scalar_name == tracer_qc) THEN
        CALL levels_horizontal_mean(sgs_flux, p_patch%cells%area, p_patch%cells%owned, &
                                    outvar)

        outvar = outvar * alv
        prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_qc_flx) =  &
               prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_qc_flx)+outvar(1:nlevp1)
      END IF
    END IF !is_sampling_time

    !$ACC WAIT
    !$ACC END DATA ! CREATE local

  END SUBROUTINE diffuse_scalar

!---------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_tke
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for cell based TKE with interpolation afterwards
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !! - Should be mergeable with diffuse_scalar
  !!------------------------------------------------------------------------
  SUBROUTINE diffuse_tke(var_ic, var, p_nh_metrics, p_patch, p_int, tot_tend_ic, p_diag, dt)
                        
    REAL(wp),             INTENT(in)        :: var_ic(:,:,:)          !< SGS-TKE half level
    REAL(wp),             INTENT(in)        :: var(:,:,:)             !< SGS-TKE full level
    TYPE(t_nh_metrics),   INTENT(in),TARGET :: p_nh_metrics           !< single nh metric state
    TYPE(t_patch),        INTENT(in),TARGET :: p_patch                !< single patch
    TYPE(t_int_state),    INTENT(in),TARGET :: p_int                  !< single interpolation state
    REAL(wp),             INTENT(inout)     :: tot_tend_ic(:,:,:)     !< total tendency
    TYPE(t_nwp_phy_diag), INTENT(inout)     :: p_diag                 !< atm phys vars
    REAL(wp),             INTENT(in)        :: dt

    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: rl_start, rl_end, jk, jb, je, jc, jg, nlev
    INTEGER  :: jcn, jvn, jbn

    INTEGER, DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk
    INTEGER, DIMENSION(:,:,:), POINTER :: ievidx, ievblk

    REAL(wp) :: flux_up, flux_dn, inv_dt, norm_metr, tang_metr
    REAL(wp) :: inv_rho_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)     
    REAL(wp) :: var_e(nproma,p_patch%nlevp1,p_patch%nblks_e)
    REAL(wp) :: var_v(nproma,p_patch%nlevp1,p_patch%nblks_v)
    REAL(wp) :: km_e(nproma,p_patch%nlevp1,p_patch%nblks_e)
    REAL(wp) :: km_v(nproma,p_patch%nlevp1,p_patch%nblks_v)   
    REAL(wp) :: metric_tend_ie(nproma,p_patch%nlev,p_patch%nblks_e)     

    REAL(wp), DIMENSION(p_patch%nlev+1)      :: var_new
    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs
    REAL(wp), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_e) :: nabla2_ie, rho_ie 
    REAL(wp), POINTER :: km_ic(:,:,:)

    !--------------------------------------------------------------------------    

    IF (msg_level >= 18) &
         CALL message(TRIM(inmodule), 'diffuse_tke')

    jg     = p_patch%id
    nlev   = p_patch%nlev
    inv_dt = 1._wp / dt

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk
    ieidx  => p_patch%cells%edge_idx
    ieblk  => p_patch%cells%edge_blk

    ievidx => p_patch%edges%vertex_idx
    ievblk => p_patch%edges%vertex_blk

    km_ic  => p_diag%tkvm

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
    CALL init(a(:,:), lacc=.FALSE.)
    CALL init(c(:,:), lacc=.FALSE.)

!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
#endif
            inv_rho_ic(jc,jk,jb) = 1._wp / rho_ic(jc,jk,jb)
            km_c(jc,jk,jb)       = MAX( les_config(jg)%km_min,  &
                                        0.5_wp * ( km_ic(jc,jk,jb) + km_ic(jc,jk+1,jb) ) ) 
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL cells2edges_scalar(rho_ic, p_patch, p_int%c_lin_e, rho_ie, opt_rlstart=grf_bdywidth_e,   &
                            opt_rlend=min_rledge_int-1)
    CALL cells2edges_scalar(km_c, p_patch, p_int%c_lin_e, km_e, opt_rlstart=grf_bdywidth_e,       &
                            opt_rlend=min_rledge_int-1)
    CALL cells2edges_scalar(var, p_patch, p_int%c_lin_e, var_e)
    CALL cells2verts_scalar(var, p_patch, p_int%cells_aw_verts, var_v)

    var_e(:,nlev+1,:) = var_e(:,nlev,:)
    var_v(:,nlev+1,:) = var_v(:,nlev,:)
    km_e  = MAX( les_config(jg)%km_min, km_e )

    !---------------------------------------------------------------
    ! Horizontal diffusion (conservative; following mo_nh_diffusion)
    !---------------------------------------------------------------
    ! compute metric contribution on half level edges
!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,norm_metr,tang_metr,jcn,jvn,jbn)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          jcn       =   iecidx(je,jb,2)
          jbn       =   iecblk(je,jb,2)

          norm_metr = - p_patch%edges%inv_dual_edge_length(je,jb) * ( km_e(je,jk,jb) *            &
                        ( var(jcn,jk,jbn) - var(iecidx(je,jb,1),jk,iecblk(je,jb,1)) ) -           &
                        km_e(je,jk+1,jb) * ( var(jcn,jk+1,jbn) -                                  &
                        var(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) ) )

          jvn       =   ievidx(je,jb,2)
          jbn       =   ievblk(je,jb,2)

          tang_metr = - p_patch%edges%inv_primal_edge_length(je,jb) *                             &
                        p_patch%edges%tangent_orientation(je,jb) * ( km_e(je,jk,jb) *             &
                        ( var_v(jvn,jk,jbn) - var_v(ievidx(je,jb,1),jk,ievblk(je,jb,1)) ) -       &
                        km_e(je,jk+1,jb) * ( var_v(jvn,jk+1,jbn) -                                &
                        var_v(ievidx(je,jb,1),jk+1,ievblk(je,jb,1)) ) )

          metric_tend_ie(je,jk,jb) = ( norm_metr * p_nh_metrics%ddxn_z_half_e(je,jk,jb)   +       &
                                       tang_metr * p_nh_metrics%ddxt_z_half_e(je,jk,jb) ) *       &
                                       p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! include halo points and boundary points because these values will be
    ! used in next loop
    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev-1
#else
      DO jk = 2, nlev-1
        DO je = i_startidx, i_endidx
#endif
          ! compute 2 * rho_ie * km_ie * grad_horiz(tke)
          nabla2_ie(je,jk,jb) = 2._wp * rho_ie(je,jk,jb) * km_ie(je,jk,jb)     *   &
                                ( p_patch%edges%inv_dual_edge_length(je,jb)    *   & 
                                ( var_ic(iecidx(je,jb,2),jk,iecblk(je,jb,2))   -   &
                                  var_ic(iecidx(je,jb,1),jk,iecblk(je,jb,1)) ) -   &
                                ( var_e(je,jk,jb) - var_e(je,jk+1,jb) )        *   &
                                p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)       *   &
                                p_nh_metrics%ddxn_z_full(je,jk,jb) )
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! divergence of km_ie * grad_horiz(e) at interface center
!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          tot_tend_ic(jc,jk,jb) = inv_rho_ic(jc,jk,jb) * (                                        &
                    nabla2_ie(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) +     &
                    nabla2_ie(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) +     &
                    nabla2_ie(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) +     &
                    metric_tend_ie(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) + &
                    metric_tend_ie(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) + &
                    metric_tend_ie(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    
    !---------------------------------------------------------------
    ! Vertical diffusion
    !---------------------------------------------------------------
    SELECT CASE(les_config(jg)%vert_scheme_type)

    ! Vertical tendency - explicit solver
    CASE(iexplicit)
#ifdef _OPENACC
      CALL finish ('mo_sgs_turbmetric:', 'vert_scheme_type=iexplicit: OpenACC version currently not implemented')
#endif

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,flux_up,flux_dn)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
#endif
            flux_up = km_c(jc,jk-1,jb) * ( var(jc,jk-1,jb) - var(jc,jk,jb) ) *                    &
                      p_nh_metrics%inv_ddqz_z_full(jc,jk-1,jb)

            flux_dn = km_c(jc,jk,jb) * ( var(jc,jk,jb) - var(jc,jk+1,jb) ) *                      &
                      p_nh_metrics%inv_ddqz_z_full(jc,jk,jb)

            tot_tend_ic(jc,jk,jb) = tot_tend_ic(jc,jk,jb) + 2._wp * inv_rho_ic(jc,jk,jb) *        &
                                    p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * ( flux_up - flux_dn +&
                                    ( flux_up * p_nh_metrics%ddxn_z_full_c(jc,jk-1,jb) -          &
                                      flux_dn * p_nh_metrics%ddxn_z_full_c(jc,jk,jb)) *           &
                                      p_nh_metrics%ddxn_z_half_c(jc,jk,jb) +                      &
                                      ( flux_up * p_nh_metrics%ddxt_z_full_c(jc,jk-1,jb) -        &
                                        flux_dn * p_nh_metrics%ddxt_z_full_c(jc,jk,jb) ) *        &
                                      p_nh_metrics%ddxt_z_half_c(jc,jk,jb) )
          ENDDO
        ENDDO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! vertical tendency - implicit solver
    ! a*x(k-1) + b*x(k) + c*x(k+1) = rhs(k)
    CASE(iimplicit)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,a,b,c,rhs,var_new)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 3, nlev-1
#else
        DO jk = 3, nlev-1
          DO jc = i_startidx, i_endidx
#endif
            a(jc,jk)   = - 2._wp * km_c(jc,jk-1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk-1,jb) *  &
                           p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb) *        &
                           ( p_nh_metrics%ddxn_z_full_c(jc,jk-1,jb) *                             & 
                             p_nh_metrics%ddxn_z_half_c(jc,jk,jb)   +                             &
                             p_nh_metrics%ddxt_z_full_c(jc,jk-1,jb) *                             &
                             p_nh_metrics%ddxt_z_half_c(jc,jk,jb) + 2._wp ) 

            c(jc,jk)   = - 2._wp * km_c(jc,jk,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) *      &
                           p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb) *        &
                           ( p_nh_metrics%ddxn_z_full_c(jc,jk,jb) *                               &
                             p_nh_metrics%ddxn_z_half_c(jc,jk,jb) +                               &
                             p_nh_metrics%ddxt_z_full_c(jc,jk,jb) *                               &
                             p_nh_metrics%ddxt_z_half_c(jc,jk,jb) + 2._wp )

            b(jc,jk)   =   inv_dt - a(jc,jk) - c(jc,jk)

            rhs(jc,jk) =   var_ic(jc,jk,jb) * inv_dt
          END DO
        END DO

        ! Boundary treatment
        !--------------------------------------------------------
        ! jk = 2
        !--------------------------------------------------------
        DO jc = i_startidx, i_endidx
          c(jc,2)   = - 2._wp * km_c(jc,2,jb) * p_nh_metrics%inv_ddqz_z_full(jc,2,jb) *           &
                        p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb) *             &
                        ( p_nh_metrics%ddxn_z_full_c(jc,2,jb) *                                   &
                          p_nh_metrics%ddxn_z_half_c(jc,2,jb) +                                   &
                          p_nh_metrics%ddxt_z_full_c(jc,2,jb) *                                   &
                          p_nh_metrics%ddxt_z_half_c(jc,2,jb) + 2._wp )

          b(jc,2)   =   inv_dt - c(jc,1) + 2._wp * km_c(jc,1,jb) *                                &
                        p_nh_metrics%inv_ddqz_z_full(jc,1,jb) *                                   &
                        p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

          rhs(jc,2) =   var_ic(jc,2,jb) * inv_dt
        END DO
        !--------------------------------------------------------
        ! jk = nlev
        !--------------------------------------------------------
        DO jc = i_startidx, i_endidx
          a(jc,nlev)   = - 2._wp * km_c(jc,nlev-1,jb) *                                           &
                           p_nh_metrics%inv_ddqz_z_full(jc,nlev-1,jb) *                           &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb) *    &
                           ( p_nh_metrics%ddxn_z_full_c(jc,nlev-1,jb) *                           &  
                             p_nh_metrics%ddxn_z_half_c(jc,nlev,jb)   +                           &
                             p_nh_metrics%ddxt_z_full_c(jc,nlev-1,jb) *                           &
                             p_nh_metrics%ddxt_z_half_c(jc,nlev,jb) + 2._wp )

          b(jc,nlev)   =   inv_dt - a(jc,nlev) + 2._wp * km_c(jc,nlev,jb) *                       &
                           p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) *                             &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)

          rhs(jc,nlev) =   var_ic(jc,nlev,jb) * inv_dt 
        END DO
        ! CALL TDMA
        DO jc = i_startidx, i_endidx
            CALL tdma_solver( a(jc,2:nlev), b(jc,2:nlev), c(jc,2:nlev), rhs(jc,2:nlev),           &
                              nlev-1,var_new(2:nlev) )

            tot_tend_ic(jc,2:nlev,jb) = tot_tend_ic(jc,2:nlev,jb) + inv_dt *                      &
                                        ( var_new(2:nlev) - var_ic(jc,2:nlev,jb) )
        END DO
      END DO !jb
!$OMP END DO
!$OMP END PARALLEL

    END SELECT !vert_scheme

  END SUBROUTINE diffuse_tke
  
!---------------------------------------------------------------------------------

END MODULE mo_sgs_turbmetric
