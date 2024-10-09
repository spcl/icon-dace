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

! Type definition for the GPU implementation of the dynamical core of ICONAM.

MODULE mo_nonhydro_gpu_types

#if defined( _OPENACC )

  USE mo_kind,                 ONLY: wp, vp
  USE mo_impl_constants,       ONLY: iaes
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d, assert_acc_device_only
  USE mo_math_types,           ONLY: t_geographical_coordinates
  USE mo_model_domain,         ONLY: t_patch, t_tangent_vectors
  USE mo_nonhydro_types,       ONLY: t_nh_state, t_nh_diag, t_nh_prog
  USE mo_prepadv_types,        ONLY: t_prepare_adv
  USE mo_advection_config,     ONLY: t_advection_config
  USE mo_les_config,           ONLY: t_les_config
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_grf_intp_data_strc,   ONLY: t_gridref_single_state, t_gridref_state
  USE mo_var_list_gpu,         ONLY: gpu_update_var_list
  USE mo_run_config,           ONLY: ltestcase, num_lev

  USE mo_intp_lonlat_types,   ONLY: lonlat_grids
  USE mo_interpol_config,     ONLY: support_baryctr_intp
  USE mo_grid_config,         ONLY: n_dom
  IMPLICIT NONE
  PRIVATE 

  PUBLIC :: h2d_icon, d2h_icon, devcpy_grf_state

CONTAINS

  SUBROUTINE h2d_icon( p_int_state, p_int_state_local_parent, p_patch, p_patch_local_parent, &
                       p_nh_state, prep_adv, advection_config, les_config, iforcing, lacc )

    TYPE ( t_int_state ),       INTENT(INOUT) :: p_int_state(:)
    TYPE ( t_int_state ),       INTENT(INOUT) :: p_int_state_local_parent(:)
    TYPE ( t_patch ),           INTENT(INOUT) :: p_patch(:)
    TYPE ( t_patch ),           INTENT(INOUT) :: p_patch_local_parent(:)
    TYPE ( t_nh_state ),        INTENT(INOUT) :: p_nh_state(:)
    TYPE ( t_prepare_adv),      INTENT(INOUT) :: prep_adv(:)
    TYPE ( t_advection_config), INTENT(INOUT) :: advection_config(:)
    TYPE ( t_les_config),       INTENT(INOUT) :: les_config(:)
    INTEGER, INTENT(IN)                       :: iforcing 
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
    INTEGER :: jg
!
! Copy all data need on GPU from host to device
!

    CALL assert_acc_device_only("h2d_icon", lacc)

    !$ACC ENTER DATA COPYIN(p_int_state, p_int_state_local_parent, p_patch, p_patch_local_parent) &
    !$ACC   COPYIN(p_nh_state, prep_adv, advection_config, les_config, num_lev)

    CALL transfer_int_state( p_int_state, .TRUE. )
    CALL transfer_int_state( p_int_state_local_parent, .TRUE. )

    CALL transfer_lonlat_grids(.TRUE.)

    CALL transfer_patch( p_patch, .TRUE. )
    CALL transfer_patch( p_patch_local_parent, .TRUE. )

    CALL transfer_prep_adv( prep_adv, .TRUE. )

    CALL transfer_nh_state( p_nh_state, .TRUE. )

    CALL transfer_advection_config( advection_config, .TRUE. )

    CALL transfer_les_config( les_config, .TRUE. )

    IF( iforcing == iaes ) THEN
      CALL transfer_aes( p_patch, .TRUE. )
    END IF

  END SUBROUTINE h2d_icon

  SUBROUTINE d2h_icon( p_int_state, p_int_state_local_parent, p_patch, p_patch_local_parent, &
                       p_nh_state, prep_adv, advection_config, les_config, iforcing, lacc )

    TYPE ( t_int_state ),  INTENT(INOUT)      :: p_int_state(:)
    TYPE ( t_int_state ),  INTENT(INOUT)      :: p_int_state_local_parent(:)
    TYPE ( t_patch ),      INTENT(INOUT)      :: p_patch(:)
    TYPE ( t_patch ),      INTENT(INOUT)      :: p_patch_local_parent(:)
    TYPE ( t_nh_state ),   INTENT(INOUT)      :: p_nh_state(:)
    TYPE ( t_prepare_adv), INTENT(INOUT)      :: prep_adv(:)
    TYPE ( t_advection_config), INTENT(INOUT) :: advection_config(:)
    TYPE ( t_les_config),       INTENT(INOUT) :: les_config(:)
    INTEGER, INTENT(IN)                       :: iforcing 
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    !
    ! Delete all data on GPU
    !
    CALL assert_acc_device_only("d2h_icon", lacc)
    !$ACC WAIT

    CALL transfer_nh_state( p_nh_state, .FALSE. )
    CALL transfer_prep_adv( prep_adv, .FALSE. )
    CALL transfer_patch( p_patch, .FALSE. )
    CALL transfer_patch( p_patch_local_parent, .FALSE. )
    CALL transfer_lonlat_grids(.FALSE.)    
    CALL transfer_int_state( p_int_state, .FALSE. )
    CALL transfer_int_state( p_int_state_local_parent, .FALSE. )
    CALL transfer_advection_config( advection_config, .FALSE. )
    CALL transfer_les_config( les_config, .FALSE. )

    IF( iforcing == iaes ) THEN
      CALL transfer_aes( p_patch, .FALSE. )
    END IF

    !$ACC WAIT(1)
#if defined(_CRAYFTN) && _RELEASE_MAJOR <= 16
    ! ACCWA (Cray Fortran 16.0.1) : p_int_state and p_int_state_local_parent are interpreted
    !  as already deleted; I think this is bug so treat as a workaround for time being
    !$ACC EXIT DATA DELETE(p_patch, p_patch_local_parent) &
    !$ACC   DELETE(p_nh_state, prep_adv, advection_config, les_config, num_lev)
#else
    !$ACC EXIT DATA DELETE(p_int_state, p_int_state_local_parent, p_patch, p_patch_local_parent) &
    !$ACC   DELETE(p_nh_state, prep_adv, advection_config, les_config, num_lev)
#endif

  END SUBROUTINE d2h_icon

  SUBROUTINE transfer_int_state( p_int, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_int_state ), TARGET,  INTENT(INOUT) :: p_int(:)

    INTEGER  :: j

    DO j= LBOUND(p_int,1), UBOUND(p_int,1) ! some p_int start at 2.

      IF ( host_to_device ) THEN

        !$ACC ENTER DATA CREATE(p_int(j))
        !$ACC ENTER DATA &
        !$ACC   COPYIN(p_int(j)%cell_environ, p_int(j)%lsq_high) &
        !$ACC   COPYIN(p_int(j)%lsq_lin)

        !$ACC ENTER DATA &
        !$ACC   COPYIN(p_int(j)%c_bln_avg, p_int(j)%c_lin_e, p_int(j)%cells_aw_verts) &
        !$ACC   COPYIN(p_int(j)%e_bln_c_s, p_int(j)%e_flx_avg, p_int(j)%geofac_div) &
        !$ACC   COPYIN(p_int(j)%geofac_grdiv, p_int(j)%geofac_grg, p_int(j)%geofac_n2s) &
        !$ACC   COPYIN(p_int(j)%geofac_rot, p_int(j)%lsq_high%lsq_blk_c) &
        !$ACC   COPYIN(p_int(j)%lsq_high%lsq_dim_stencil, p_int(j)%lsq_high%lsq_idx_c) &
        !$ACC   COPYIN(p_int(j)%lsq_high%lsq_moments, p_int(j)%lsq_high%lsq_moments_hat) &
        !$ACC   COPYIN(p_int(j)%lsq_high%lsq_pseudoinv, p_int(j)%lsq_high%lsq_qtmat_c) &
        !$ACC   COPYIN(p_int(j)%lsq_high%lsq_rmat_utri_c, p_int(j)%lsq_high%lsq_weights_c) &
        !$ACC   COPYIN(p_int(j)%lsq_lin%lsq_blk_c) &
        !$ACC   COPYIN(p_int(j)%lsq_lin%lsq_dim_stencil, p_int(j)%lsq_lin%lsq_idx_c) &
        !$ACC   COPYIN(p_int(j)%lsq_lin%lsq_moments, p_int(j)%lsq_lin%lsq_moments_hat) &
        !$ACC   COPYIN(p_int(j)%lsq_lin%lsq_pseudoinv, p_int(j)%lsq_lin%lsq_qtmat_c) &
        !$ACC   COPYIN(p_int(j)%lsq_lin%lsq_rmat_utri_c, p_int(j)%lsq_lin%lsq_weights_c) &
        !$ACC   COPYIN(p_int(j)%nudgecoeff_c, p_int(j)%nudgecoeff_e, p_int(j)%pos_on_tplane_e) &
        !$ACC   COPYIN(p_int(j)%rbf_c2grad_blk, p_int(j)%rbf_c2grad_idx, p_int(j)%rbf_c2grad_coeff) &
        !$ACC   COPYIN(p_int(j)%rbf_vec_blk_c, p_int(j)%rbf_vec_idx_c, p_int(j)%rbf_vec_coeff_c) &
        !$ACC   COPYIN(p_int(j)%rbf_vec_blk_e, p_int(j)%rbf_vec_idx_e, p_int(j)%rbf_vec_coeff_e) &
        !$ACC   COPYIN(p_int(j)%rbf_vec_blk_v, p_int(j)%rbf_vec_idx_v, p_int(j)%rbf_vec_coeff_v) &
        !$ACC   COPYIN(p_int(j)%verts_aw_cells, p_int(j)%cell_environ%idx) &
        !$ACC   COPYIN(p_int(j)%cell_environ%blk, p_int(j)%cell_environ%area_norm) &
        !$ACC   COPYIN(p_int(j)%pos_on_tplane_c_edge)

      ELSE

        !$ACC WAIT(1)
        !$ACC EXIT DATA &
        !$ACC   DELETE(p_int(j)%c_bln_avg, p_int(j)%c_lin_e, p_int(j)%cells_aw_verts) &
        !$ACC   DELETE(p_int(j)%e_bln_c_s, p_int(j)%e_flx_avg, p_int(j)%geofac_div) &
        !$ACC   DELETE(p_int(j)%geofac_grdiv, p_int(j)%geofac_grg, p_int(j)%geofac_n2s) &
        !$ACC   DELETE(p_int(j)%geofac_rot, p_int(j)%lsq_high%lsq_blk_c) &
        !$ACC   DELETE(p_int(j)%lsq_high%lsq_dim_stencil, p_int(j)%lsq_high%lsq_idx_c) &
        !$ACC   DELETE(p_int(j)%lsq_high%lsq_moments, p_int(j)%lsq_high%lsq_moments_hat) &
        !$ACC   DELETE(p_int(j)%lsq_high%lsq_pseudoinv, p_int(j)%lsq_high%lsq_qtmat_c) &
        !$ACC   DELETE(p_int(j)%lsq_high%lsq_rmat_utri_c, p_int(j)%lsq_high%lsq_weights_c) &
        !$ACC   DELETE(p_int(j)%lsq_lin%lsq_blk_c) &
        !$ACC   DELETE(p_int(j)%lsq_lin%lsq_dim_stencil, p_int(j)%lsq_lin%lsq_idx_c) &
        !$ACC   DELETE(p_int(j)%lsq_lin%lsq_moments, p_int(j)%lsq_lin%lsq_moments_hat) &
        !$ACC   DELETE(p_int(j)%lsq_lin%lsq_pseudoinv, p_int(j)%lsq_lin%lsq_qtmat_c) &
        !$ACC   DELETE(p_int(j)%lsq_lin%lsq_rmat_utri_c, p_int(j)%lsq_lin%lsq_weights_c) &
        !$ACC   DELETE(p_int(j)%nudgecoeff_c, p_int(j)%nudgecoeff_e, p_int(j)%pos_on_tplane_e) &
        !$ACC   DELETE(p_int(j)%rbf_c2grad_blk, p_int(j)%rbf_c2grad_idx, p_int(j)%rbf_c2grad_coeff) &
        !$ACC   DELETE(p_int(j)%rbf_vec_blk_c, p_int(j)%rbf_vec_idx_c, p_int(j)%rbf_vec_coeff_c) &
        !$ACC   DELETE(p_int(j)%rbf_vec_blk_e, p_int(j)%rbf_vec_idx_e, p_int(j)%rbf_vec_coeff_e) &
        !$ACC   DELETE(p_int(j)%rbf_vec_blk_v, p_int(j)%rbf_vec_idx_v, p_int(j)%rbf_vec_coeff_v) &
        !$ACC   DELETE(p_int(j)%verts_aw_cells) &
        !$ACC   DELETE(p_int(j)%cell_environ%idx, p_int(j)%cell_environ%blk) &
        !$ACC   DELETE(p_int(j)%cell_environ%area_norm, p_int(j)%pos_on_tplane_c_edge) &
        !$ACC   DELETE(p_int(j))

        !$ACC EXIT DATA &
        !$ACC   DELETE(p_int(j)%lsq_high, p_int(j)%lsq_lin, p_int(j)%cell_environ)

        !$ACC EXIT DATA DELETE(p_int(j))

      ENDIF

    ENDDO

  END SUBROUTINE transfer_int_state


  SUBROUTINE transfer_patch( p_patch, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_patch ), TARGET, INTENT(INOUT)    :: p_patch(:)

    INTEGER :: j

!
! Copy the static data structures in p_patch to the device -- this is a small subset of all the components
! The communication patterns are copied over in mo_communication_orig.
!
      DO j=1,SIZE(p_patch)

        IF ( host_to_device ) THEN
        
        !$ACC ENTER DATA &
        !$ACC   COPYIN(p_patch(j)%cells, p_patch(j)%edges, p_patch(j)%verts)

        !$ACC ENTER DATA &
        !$ACC   COPYIN(p_patch(j)%cells%decomp_info)
                
        !$ACC ENTER DATA &
        !$ACC   COPYIN(p_patch(j)%cells%decomp_info%owner_mask) &
        !$ACC   COPYIN(p_patch(j)%cells%ddqz_z_full, p_patch(j)%cells%area) &
        !$ACC   COPYIN(p_patch(j)%cells%edge_idx, p_patch(j)%cells%edge_blk) &
        !$ACC   COPYIN(p_patch(j)%cells%neighbor_idx, p_patch(j)%cells%neighbor_blk) &
        !$ACC   COPYIN(p_patch(j)%cells%center, p_patch(j)%cells%refin_ctrl, p_patch(j)%cells%f_c) &
        !$ACC   COPYIN(p_patch(j)%cells%vertex_blk, p_patch(j)%cells%vertex_idx) &
        !$ACC   COPYIN(p_patch(j)%cells%child_blk, p_patch(j)%cells%child_idx) &
        !$ACC   COPYIN(p_patch(j)%cells%start_index, p_patch(j)%cells%end_index) &
        !$ACC   COPYIN(p_patch(j)%edges%area_edge, p_patch(j)%edges%cell_idx) &
        !$ACC   COPYIN(p_patch(j)%edges%cell_blk, p_patch(j)%edges%edge_cell_length, p_patch(j)%edges%f_e) &
        !$ACC   COPYIN(p_patch(j)%edges%quad_idx, p_patch(j)%edges%quad_blk, p_patch(j)%edges%vertex_idx) &
        !$ACC   COPYIN(p_patch(j)%edges%vertex_blk, p_patch(j)%edges%primal_normal_cell) &
        !$ACC   COPYIN(p_patch(j)%edges%start_index, p_patch(j)%edges%end_index) &
        !$ACC   COPYIN(p_patch(j)%edges%dual_normal_cell, p_patch(j)%edges%primal_normal_vert) &
        !$ACC   COPYIN(p_patch(j)%edges%dual_normal_vert, p_patch(j)%edges%inv_vert_vert_length) &
        !$ACC   COPYIN(p_patch(j)%edges%inv_dual_edge_length, p_patch(j)%edges%inv_primal_edge_length) &
        !$ACC   COPYIN(p_patch(j)%edges%primal_edge_length, p_patch(j)%edges%edge_vert_length) &
        !$ACC   COPYIN(p_patch(j)%edges%child_idx, p_patch(j)%edges%child_blk) &
        !$ACC   COPYIN(p_patch(j)%edges%tangent_orientation, p_patch(j)%edges%refin_ctrl) &
        !$ACC   COPYIN(p_patch(j)%edges%parent_loc_idx, p_patch(j)%edges%parent_loc_blk) &
        !$ACC   COPYIN(p_patch(j)%edges%butterfly_idx, p_patch(j)%edges%butterfly_blk) &
        !$ACC   COPYIN(p_patch(j)%verts%cell_idx, p_patch(j)%verts%cell_blk) &
        !$ACC   COPYIN(p_patch(j)%verts%start_index, p_patch(j)%verts%end_index, p_patch(j)%edges%pc_idx) &
        !$ACC   COPYIN(p_patch(j)%verts%edge_idx, p_patch(j)%verts%edge_blk, p_patch(j)%verts%refin_ctrl) &
        !$ACC   COPYIN(p_patch(j)%edges%center, p_patch(j)%edges%primal_normal)

      ELSE

        !$ACC WAIT(1)
        !$ACC EXIT DATA &
        !$ACC   DELETE(p_patch(j)%cells%decomp_info%owner_mask) &
        !$ACC   DELETE(p_patch(j)%cells%ddqz_z_full, p_patch(j)%cells%area) &
        !$ACC   DELETE(p_patch(j)%cells%edge_idx, p_patch(j)%cells%edge_blk) &
        !$ACC   DELETE(p_patch(j)%cells%neighbor_idx, p_patch(j)%cells%neighbor_blk) &
        !$ACC   DELETE(p_patch(j)%cells%center, p_patch(j)%cells%refin_ctrl, p_patch(j)%cells%f_c) &
        !$ACC   DELETE(p_patch(j)%cells%vertex_blk, p_patch(j)%cells%vertex_idx) &
        !$ACC   DELETE(p_patch(j)%cells%child_blk, p_patch(j)%cells%child_idx) &
        !$ACC   DELETE(p_patch(j)%cells%start_index, p_patch(j)%cells%end_index) &
        !$ACC   DELETE(p_patch(j)%edges%area_edge, p_patch(j)%edges%cell_idx) &
        !$ACC   DELETE(p_patch(j)%edges%cell_blk, p_patch(j)%edges%edge_cell_length, p_patch(j)%edges%f_e) &
        !$ACC   DELETE(p_patch(j)%edges%quad_idx, p_patch(j)%edges%quad_blk, p_patch(j)%edges%vertex_idx) &
        !$ACC   DELETE(p_patch(j)%edges%vertex_blk, p_patch(j)%edges%primal_normal_cell) &
        !$ACC   DELETE(p_patch(j)%edges%start_index, p_patch(j)%edges%end_index) &
        !$ACC   DELETE(p_patch(j)%edges%dual_normal_cell, p_patch(j)%edges%primal_normal_vert) &
        !$ACC   DELETE(p_patch(j)%edges%dual_normal_vert, p_patch(j)%edges%inv_vert_vert_length) &
        !$ACC   DELETE(p_patch(j)%edges%inv_dual_edge_length, p_patch(j)%edges%inv_primal_edge_length) &
        !$ACC   DELETE(p_patch(j)%edges%primal_edge_length, p_patch(j)%edges%edge_vert_length) &
        !$ACC   DELETE(p_patch(j)%edges%child_idx, p_patch(j)%edges%child_blk) &
        !$ACC   DELETE(p_patch(j)%edges%tangent_orientation, p_patch(j)%edges%refin_ctrl) &
        !$ACC   DELETE(p_patch(j)%edges%parent_loc_idx, p_patch(j)%edges%parent_loc_blk) &
        !$ACC   DELETE(p_patch(j)%edges%butterfly_idx, p_patch(j)%edges%butterfly_blk) &
        !$ACC   DELETE(p_patch(j)%verts%cell_idx, p_patch(j)%verts%cell_blk) &
        !$ACC   DELETE(p_patch(j)%verts%start_index, p_patch(j)%verts%end_index, p_patch(j)%edges%pc_idx) &
        !$ACC   DELETE(p_patch(j)%verts%edge_idx, p_patch(j)%verts%edge_blk, p_patch(j)%verts%refin_ctrl) &
        !$ACC   DELETE(p_patch(j)%edges%center, p_patch(j)%edges%primal_normal)

        !$ACC EXIT DATA &
        !$ACC   DELETE(p_patch(j)%cells%decomp_info)

        !$ACC EXIT DATA &
        !$ACC   DELETE(p_patch(j)%cells, p_patch(j)%edges, p_patch(j)%verts)

      ENDIF

    ENDDO

  END SUBROUTINE transfer_patch


  SUBROUTINE transfer_prep_adv( prep_adv, host_to_device )

    LOGICAL, INTENT(IN)                           :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_prepare_adv ), TARGET, INTENT(INOUT) :: prep_adv(:)

    INTEGER :: jg

    DO jg=1, SIZE(prep_adv)
      CALL gpu_update_var_list('prepadv_of_domain_', host_to_device, domain=jg, lacc=.TRUE. )
    ENDDO    

  END SUBROUTINE transfer_prep_adv


  SUBROUTINE transfer_advection_config( advection_config, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_advection_config ), TARGET, INTENT(INOUT) :: advection_config(:)

    INTEGER :: j

    !$ACC ENTER DATA COPYIN(advection_config) IF(host_to_device)
    DO j=1, SIZE(advection_config)

      IF ( host_to_device ) THEN      
        !$ACC ENTER DATA COPYIN(advection_config(j)%trHydroMass%list, advection_config(j)%trAdvect%list)
      ELSE
        !$ACC WAIT(1)
        !$ACC EXIT DATA DELETE(advection_config(j)%trHydroMass%list, advection_config(j)%trAdvect%list)
      ENDIF

    ENDDO
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(advection_config) IF(.NOT. host_to_device)

  END SUBROUTINE transfer_advection_config

  SUBROUTINE transfer_les_config( les_config, host_to_device )

    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_les_config ), TARGET, INTENT(INOUT) :: les_config(:)

    INTEGER :: j

    IF ( host_to_device ) THEN
        !$ACC ENTER DATA COPYIN(les_config)
    ELSE
        !$ACC WAIT(1)
        !$ACC EXIT DATA DELETE(les_config)
    ENDIF

  END SUBROUTINE transfer_les_config

  SUBROUTINE transfer_nh_state( p_nh, host_to_device )
    LOGICAL, INTENT(IN)                        :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    TYPE ( t_nh_state ), TARGET, INTENT(INOUT) :: p_nh(:)
    INTEGER :: istep, jg
    CHARACTER(*), PARAMETER :: &
      & metrics = 'nh_state_metrics_of_domain_', diag = 'nh_state_diag_of_domain_', &
      & ref = 'nh_state_ref_of_domain_', prog = 'nh_state_prog_of_domain_'

! At this point, p_nh and all its underlying subtypes have been created on the device
! HB: merged interfaces of gpu_XXX_var_list... therefore the IF condition
! WS:  currently it appears to be unnecessary to update any of these values back to the host 
!      after the end of the time loop.  Dycore variables are updated in ACC_VALIDATE mode individually
!      BUT: there should be a way to delete all the variables with DEL_VAR
    IF (.NOT.host_to_device) RETURN
    DO jg = 1, SIZE(p_nh)
      CALL gpu_update_var_list(metrics, host_to_device, domain=jg, lacc=.TRUE. )
      IF (ltestcase) CALL gpu_update_var_list(ref, host_to_device, domain=jg, lacc=.TRUE. )
      CALL gpu_update_var_list(diag, host_to_device, domain=jg, lacc=.TRUE. )
      DO istep = 1, SIZE(p_nh(jg)%prog)
        CALL gpu_update_var_list(prog, host_to_device, domain=jg, substr='_and_timelev_', timelev=istep, lacc=.TRUE. )
      ENDDO
    ENDDO
  END SUBROUTINE transfer_nh_state

  SUBROUTINE transfer_aes( p_patch, host_to_device )
    TYPE ( t_patch ),      INTENT(INOUT) :: p_patch(:)
    LOGICAL, INTENT(IN)                  :: host_to_device     !   .TRUE. : h2d   .FALSE. : d2h
    INTEGER :: jg

    DO jg = 1, SIZE(p_patch)
      CALL gpu_update_var_list('prm_field_D', host_to_device, domain=jg, lacc=.TRUE.)
      CALL gpu_update_var_list('prm_tend_D', host_to_device, domain=jg, lacc=.TRUE.)
    END DO
  END SUBROUTINE transfer_aes

  SUBROUTINE devcpy_grf_state( p_grf, l_h2d, lacc )

      TYPE ( t_gridref_state ), TARGET,  INTENT(INOUT) :: p_grf(:)
      LOGICAL, INTENT(IN) :: l_h2d    ! true host-to-device, false device-to-host
      LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc


      INTEGER  :: j,k

      CALL assert_acc_device_only("devcpy_grf_state", lacc)

      IF (l_h2d) THEN
        !$ACC ENTER DATA COPYIN(p_grf)

        DO j=1, SIZE(p_grf)

          !$ACC ENTER DATA &
          !$ACC   COPYIN(p_grf(j)%fbk_wgt_aw, p_grf(j)%fbk_wgt_bln, p_grf(j)%fbk_wgt_e) &
          !$ACC   COPYIN(p_grf(j)%mask_ovlp_c, p_grf(j)%mask_ovlp_ch, p_grf(j)%mask_ovlp_e, p_grf(j)%mask_ovlp_v) &
          !$ACC   COPYIN(p_grf(j)%idxlist_bdyintp_src_c, p_grf(j)%idxlist_bdyintp_src_e, p_grf(j)%blklist_bdyintp_src_c) &
          !$ACC   COPYIN(p_grf(j)%blklist_bdyintp_src_e, p_grf(j)%p_dom)

          CALL devcpy_grf_single_state( p_grf(j)%p_dom, l_h2d )
        ENDDO

      ELSE

        !$ACC WAIT(1)
        DO j=1, SIZE(p_grf)
          CALL devcpy_grf_single_state( p_grf(j)%p_dom, l_h2d )

          !$ACC EXIT DATA &
          !$ACC   DELETE(p_grf(j)%fbk_wgt_aw, p_grf(j)%fbk_wgt_bln, p_grf(j)%fbk_wgt_e) &
          !$ACC   DELETE(p_grf(j)%mask_ovlp_c, p_grf(j)%mask_ovlp_ch, p_grf(j)%mask_ovlp_e, p_grf(j)%mask_ovlp_v) &
          !$ACC   DELETE(p_grf(j)%idxlist_bdyintp_src_c, p_grf(j)%idxlist_bdyintp_src_e, p_grf(j)%blklist_bdyintp_src_c) &
          !$ACC   DELETE(p_grf(j)%blklist_bdyintp_src_e, p_grf(j)%p_dom)

        ENDDO

        !$ACC EXIT DATA DELETE(p_grf)
      ENDIF


    END SUBROUTINE devcpy_grf_state

    SUBROUTINE devcpy_grf_single_state( p_grf, l_h2d )

      TYPE ( t_gridref_single_state ), TARGET,  INTENT(INOUT) :: p_grf(:)
      LOGICAL, INTENT(IN) :: l_h2d    ! true host-to-device, false device-to-host

      INTEGER  :: j


      DO j=1, SIZE(p_grf)

        IF (l_h2d) THEN

          !$ACC ENTER DATA &
          !$ACC   COPYIN(p_grf(j)%grf_dist_pc2cc, p_grf(j)%grf_dist_pe2ce, p_grf(j)%idxlist_bdyintp_c) &
          !$ACC   COPYIN(p_grf(j)%idxlist_bdyintp_e, p_grf(j)%idxlist_ubcintp_c, p_grf(j)%idxlist_ubcintp_e, p_grf(j)%blklist_bdyintp_c) &
          !$ACC   COPYIN(p_grf(j)%blklist_bdyintp_e, p_grf(j)%blklist_ubcintp_c, p_grf(j)%blklist_ubcintp_e, p_grf(j)%idxlist_rbfintp_v) &
          !$ACC   COPYIN(p_grf(j)%blklist_rbfintp_v, p_grf(j)%edge_vert_idx, p_grf(j)%coeff_bdyintp_c, p_grf(j)%coeff_ubcintp_c) &
          !$ACC   COPYIN(p_grf(j)%dist_pc2cc_bdy, p_grf(j)%dist_pc2cc_ubc, p_grf(j)%prim_norm, p_grf(j)%coeff_bdyintp_e12) &
          !$ACC   COPYIN(p_grf(j)%coeff_bdyintp_e34, p_grf(j)%dist_pe2ce, p_grf(j)%coeff_ubcintp_e12, p_grf(j)%coeff_ubcintp_e34) &
          !$ACC   COPYIN(p_grf(j)%grf_vec_ind_2a, p_grf(j)%grf_vec_blk_2a, p_grf(j)%grf_vec_ind_2b, p_grf(j)%grf_vec_blk_2b) &
          !$ACC   COPYIN(p_grf(j)%grf_vec_coeff_2a, p_grf(j)%grf_vec_coeff_2b) &
          !$ACC   COPYIN(p_grf(j)%coeff_rbf_v)

        ELSE

          !$ACC WAIT(1)
          !$ACC EXIT DATA &
          !$ACC   DELETE(p_grf(j)%grf_dist_pc2cc, p_grf(j)%grf_dist_pe2ce, p_grf(j)%idxlist_bdyintp_c) &
          !$ACC   DELETE(p_grf(j)%idxlist_bdyintp_e, p_grf(j)%idxlist_ubcintp_c, p_grf(j)%idxlist_ubcintp_e, p_grf(j)%blklist_bdyintp_c) &
          !$ACC   DELETE(p_grf(j)%blklist_bdyintp_e, p_grf(j)%blklist_ubcintp_c, p_grf(j)%blklist_ubcintp_e, p_grf(j)%idxlist_rbfintp_v) &
          !$ACC   DELETE(p_grf(j)%blklist_rbfintp_v, p_grf(j)%edge_vert_idx, p_grf(j)%coeff_bdyintp_c, p_grf(j)%coeff_ubcintp_c) &
          !$ACC   DELETE(p_grf(j)%dist_pc2cc_bdy, p_grf(j)%dist_pc2cc_ubc, p_grf(j)%prim_norm, p_grf(j)%coeff_bdyintp_e12) &
          !$ACC   DELETE(p_grf(j)%coeff_bdyintp_e34, p_grf(j)%dist_pe2ce, p_grf(j)%coeff_ubcintp_e12, p_grf(j)%coeff_ubcintp_e34) &
          !$ACC   DELETE(p_grf(j)%grf_vec_ind_2a, p_grf(j)%grf_vec_blk_2a, p_grf(j)%grf_vec_ind_2b, p_grf(j)%grf_vec_blk_2b) &
          !$ACC   DELETE(p_grf(j)%grf_vec_coeff_2a, p_grf(j)%grf_vec_coeff_2b) &
          !$ACC   DELETE(p_grf(j)%coeff_rbf_v)

        ENDIF

      ENDDO

    END SUBROUTINE devcpy_grf_single_state


    SUBROUTINE transfer_lonlat_grids(l_h2d)
      LOGICAL, INTENT(IN) :: l_h2d    ! true host-to-device, false device-to-host
      INTEGER :: jg, i

      ! The derived type components are (de)allocated and created/deleted on
      ! the device in their *init / *finalize routines

      IF (l_h2d) THEN

        DO jg=1,n_dom
          DO i=1, lonlat_grids%ngrids
            IF (lonlat_grids%list(i)%l_dom(jg) .AND. lonlat_grids%list(i)%intp(jg)%l_initialized) THEN
              !$ACC ENTER DATA CREATE(lonlat_grids%list(i)%intp(jg:jg))
              !$ACC UPDATE &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_vec%stencil) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_vec%idx) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_vec%blk) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_vec%coeff) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_c2l%stencil) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_c2l%idx) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_c2l%blk) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_c2l%coeff) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%rbf_c2l%l_cutoff) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%nnb%stencil) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%nnb%idx) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%nnb%blk) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%nnb%coeff) &
              !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%nnb%l_cutoff) &
              !$ACC   ASYNC(1)

              IF (support_baryctr_intp) THEN
                !$ACC UPDATE &
                !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%baryctr%stencil) &
                !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%baryctr%idx) &
                !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%baryctr%blk) &
                !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%baryctr%coeff) &
                !$ACC   DEVICE(lonlat_grids%list(i)%intp(jg)%baryctr%l_cutoff) &
                !$ACC   ASYNC(1)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ELSE
        DO jg=1,n_dom
          DO i=1, lonlat_grids%ngrids
            IF (lonlat_grids%list(i)%l_dom(jg) .AND. lonlat_grids%list(i)%intp(jg)%l_initialized) THEN
              !ptr_int_lonlat => lonlat_grids%list(i)%intp(jg)%rbf_c2l
              !$ACC WAIT(1)
              !$ACC EXIT DATA DELETE(lonlat_grids%list(i)%intp(jg:jg))
            ENDIF
          ENDDO
        ENDDO
      ENDIF

    END SUBROUTINE transfer_lonlat_grids
#endif



END MODULE mo_nonhydro_gpu_types
