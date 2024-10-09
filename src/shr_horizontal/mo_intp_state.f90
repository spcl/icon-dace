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

! Contains the implementation of interpolation and reconstruction
! routines used by the shallow water model, including the RBF
! reconstruction routines.

! #ifdef __xlC__
! @PROCESS HOT
! #endif
#ifdef __PGI
!pgi$g opt=0
#endif

MODULE mo_intp_state
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message, finish
USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
USE mo_model_domain,        ONLY: t_patch
USE mo_grid_config,         ONLY: n_dom, n_dom_start, lplane, l_limited_area
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: ldynamics
USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_c2grad_dim,                &
  &                               rbf_vec_dim_v, rbf_vec_dim_e, lsq_lin_set,    &
  &                               lsq_high_set
USE mo_initicon_config,     ONLY: icpl_da_seaice, icpl_da_snowalb
USE mo_lnd_nwp_config,      ONLY: lterra_urb
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_intp_rbf_coeffs,     ONLY: rbf_vec_index_cell, rbf_vec_index_edge,                &
  &                               rbf_vec_index_vertex, rbf_vec_compute_coeff_cell,      &
  &                               rbf_vec_compute_coeff_edge,                            &
  &                               rbf_vec_compute_coeff_vertex, rbf_c2grad_index,        &
  &                               rbf_compute_coeff_c2grad, gen_index_list_radius
USE mo_intp_coeffs,         ONLY: init_cellavg_wgt,                                    &
  &                               init_geo_factors, complete_patchinfo, init_tplane_e, &
  &                               init_tplane_c, tri_quadrature_pts,                   &
  &                               init_nudgecoeffs
USE mo_intp_coeffs_lsq_bln, ONLY: lsq_stencil_create, lsq_compute_coeff_cell,          &
  &                               scalar_int_coeff, bln_int_coeff_e2c
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V
USE mo_communication,       ONLY: t_comm_pattern, blk_no, idx_no, idx_1d, &
  &                               delete_comm_pattern, exchange_data
  USE mo_communication_factory, ONLY: setup_comm_pattern
USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info, get_valid_local_index
USE mo_dist_dir,            ONLY: dist_dir_get_owners
USE mo_update_dyn_scm ,     ONLY: rbf_coeff_scm
USE mo_grid_config,         ONLY: l_scm_mode
USE mo_name_list_output_config, ONLY: is_variable_in_output
USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
USE mo_master_control,      ONLY: my_process_is_waves
#ifdef SERIALIZE
USE mo_ser_rbf_coefficients, ONLY: ser_rbf_coefficients
USE mo_ser_nml,             ONLY: ser_rbf
#endif

IMPLICIT NONE

PRIVATE

!!CJPUBLIC :: setup_interpol,
PUBLIC :: construct_2d_interpol_state, destruct_2d_interpol_state
PUBLIC :: transfer_interpol_state
PUBLIC :: allocate_int_state, deallocate_int_state

INTERFACE xfer_var
  MODULE PROCEDURE xfer_var_r2
  MODULE PROCEDURE xfer_var_r3
  MODULE PROCEDURE xfer_var_r4
  MODULE PROCEDURE xfer_var_i2
END INTERFACE

INTERFACE xfer_idx
  MODULE PROCEDURE xfer_idx_2
  MODULE PROCEDURE xfer_idx_3
END INTERFACE

CLASS(t_comm_pattern), POINTER :: comm_pat_glb_to_loc_c, &
  &                               comm_pat_glb_to_loc_e, &
  &                               comm_pat_glb_to_loc_v


CONTAINS

!-------------------------------------------------------------------------
!
!
!> Allocation of components of interpolation state.
SUBROUTINE allocate_int_state( ptr_patch, ptr_int)
!
  TYPE(t_patch), INTENT(IN) :: ptr_patch

  TYPE(t_int_state), INTENT(inout) :: ptr_int

  INTEGER :: nblks_c, nblks_e, nblks_v
  INTEGER :: ist
  INTEGER :: idummy
  LOGICAL :: lsdi   ,&  ! deleted the former expl. init with .FALSE.
             llpi   ,&  ! to avoid the implicit SAVE attribute, which
             llpim  ,&  ! would be dangerous for subsequent calls
             llsc   ,&
             llsd   ,&
             llde   ,&
             lmconv

!-----------------------------------------------------------------------

  !
  ! determine size of arrays, i.e.
  ! values for the blocking
  !
  nblks_c  = ptr_patch%nblks_c
  nblks_e  = ptr_patch%nblks_e
  nblks_v  = ptr_patch%nblks_v
  !
  !
  ! allocate interpolation state
  !
  ! c_lin_e
  !
  ALLOCATE (ptr_int%c_lin_e(nproma,2,nblks_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',  &
    &            'allocation for c_lin_e failed')
  ENDIF

  IF (ptr_patch%geometry_info%cell_type == 3) THEN
    !
    ! e_bln_c_s
    !
    ALLOCATE (ptr_int%e_bln_c_s(nproma,3,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_s failed')
    ENDIF
    !
    ! e_bln_c_u
    !
    ALLOCATE (ptr_int%e_bln_c_u(nproma,6,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_u failed')
    ENDIF
    !
    ! e_bln_c_v
    !
    ALLOCATE (ptr_int%e_bln_c_v(nproma,6,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_bln_c_v failed')
    ENDIF
    !
    ! c_bln_avg
    !
    ALLOCATE (ptr_int%c_bln_avg(nproma,4,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for c_bln_avg failed')
    ENDIF
    !
    ! gradc_bmat
    !
    ALLOCATE (ptr_int%gradc_bmat(nproma,2,3,nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for gradc_bmat failed')
    ENDIF
    !
    ! e_flx_avg
    !
    ALLOCATE (ptr_int%e_flx_avg(nproma,5,nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state', &
      &            'allocation for e_flx_avg failed')
    ENDIF

  ENDIF
  !
  ! v_1o2_e
  !
  ALLOCATE (ptr_int%v_1o2_e(nproma,2,nblks_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &            'allocation for v_1o2_e failed')
  ENDIF
  !
  ! e_inn_c
  !
  ALLOCATE (ptr_int%e_inn_c(nproma,ptr_patch%geometry_info%cell_type,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &            'allocation for e_inn_c failed')
  ENDIF
  !
  ! verts_aw_cells
  !
  ALLOCATE (ptr_int%verts_aw_cells(nproma,ptr_patch%geometry_info%cell_type,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for verts_aw_cells failed')
  ENDIF
  !
  ! cells_aw_verts
  !
  ALLOCATE (ptr_int%cells_aw_verts(nproma,9-ptr_patch%geometry_info%cell_type,nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                     &
      &             'allocation for cells_aw_verts failed')
  ENDIF


  IF (ptr_patch%geometry_info%cell_type == 3) THEN
    !
    ! rbf_vec_idx_c, rbf_vec_blk_c
    !
    ALLOCATE (ptr_int%rbf_vec_idx_c(rbf_vec_dim_c, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_c failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_c(rbf_vec_dim_c, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_blk_c failed')
    ENDIF
    !
    ! rbf_vec_stencil_c
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_c(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_c failed')
    ENDIF
    !
    ! rbf_vec_coeff_c
    !
    ALLOCATE (ptr_int%rbf_vec_coeff_c(rbf_vec_dim_c, 2, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_c failed')
    ENDIF
    !
    ! rbf_c2grad_idx, rbf_c2grad_blk
    !
    ALLOCATE (ptr_int%rbf_c2grad_idx(rbf_c2grad_dim, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_idx failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_c2grad_blk(rbf_c2grad_dim, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_blk failed')
    ENDIF
    !
    ! rbf_c2grad_coeff
    !
    ALLOCATE (ptr_int%rbf_c2grad_coeff(rbf_c2grad_dim, 2, nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_c2grad_coeff failed')
    ENDIF
    !
    ! rbf_vec_idx_v, rbf_vec_blk_v
    !
    ALLOCATE (ptr_int%rbf_vec_idx_v(rbf_vec_dim_v, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_v failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_v(rbf_vec_dim_v, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
        &          'allocation for rbf_vec_blk_v failed')
    ENDIF
    !
    ! rbf_vec_stencil_v
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_v(nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_v failed')
    ENDIF
    !
    ! rbf_vec_coeff_v
    !
    ALLOCATE (ptr_int%rbf_vec_coeff_v(rbf_vec_dim_v, 2, nproma, nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_v failed')
    ENDIF
    !
    ! rbf_vec_idx_e, rbf_vec_blk_e
    !
    ALLOCATE (ptr_int%rbf_vec_idx_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_idx_e failed')
    ENDIF
    ALLOCATE (ptr_int%rbf_vec_blk_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_blk_e failed')
    ENDIF
    !
    ! rbf_vec_stencil_e
    !
    ALLOCATE (ptr_int%rbf_vec_stencil_e(nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_stencil_e failed')
    ENDIF
    !
    ! rbf_vec_coeff_e
    !
    ALLOCATE (ptr_int%rbf_vec_coeff_e(rbf_vec_dim_e, nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for rbf_vec_coeff_e failed')
    ENDIF

    lsdi         = .FALSE.
    llpi         = .FALSE.
    llpim        = .FALSE.
    lmconv       = .FALSE.

    ! GZ: offloading 'is_variable_in_output' to vector hosts requires separate calls in order to
    !     avoid an MPI deadlock in p_bcast
                    lsdi         = is_variable_in_output(var_name="sdi2")
    IF (.NOT.lsdi)  llpi         = is_variable_in_output(var_name="lpi")
    IF (.NOT.llpi)  llpim        = is_variable_in_output(var_name="lpi_max")
    IF (.NOT.llpim) lmconv       = is_variable_in_output(var_name="mconv")
    llsc         = atm_phy_nwp_config(MAX(1,ptr_patch%id))%lstoch_expl 
    llsd         = atm_phy_nwp_config(MAX(1,ptr_patch%id))%lstoch_sde
    llde         = atm_phy_nwp_config(MAX(1,ptr_patch%id))%lstoch_deep
    
    ptr_int%cell_environ%is_used = lsdi .OR. llpi .OR. llpim .OR. lmconv .OR. llsc .OR. llsd .OR. llde .OR. &
                                   icpl_da_seaice >= 2 .OR. icpl_da_snowalb >= 2 .OR. lterra_urb

    IF ( ptr_int%cell_environ%is_used ) THEN
      !
      ! cell_environ
      !
      ALLOCATE( ptr_int%cell_environ%nmbr_nghbr_cells( nproma, nblks_c), STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ('mo_interpolation:construct_int_state',&
        &            'allocation for cell_environ%nmbr_nghbr_cells failed')
      ENDIF

      ALLOCATE( ptr_int%cell_environ%idx( nproma, nblks_c, ptr_int%cell_environ%nmbr_nghbr_cells_alloc), STAT=ist )
      IF (ist /= SUCCESS)  THEN
        CALL finish ('mo_interpolation:construct_int_state',&
        &            'allocation for cell_environ%idx failed')
      ENDIF
      ALLOCATE( ptr_int%cell_environ%blk( nproma, nblks_c, ptr_int%cell_environ%nmbr_nghbr_cells_alloc), STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ('mo_interpolation:construct_int_state',&
        &            'allocation for cell_environ%blk failed')
      ENDIF
      ALLOCATE( ptr_int%cell_environ%area_norm( nproma, nblks_c, ptr_int%cell_environ%nmbr_nghbr_cells_alloc), STAT=ist )
      IF (ist /= SUCCESS)  THEN
        CALL finish ('mo_interpolation:construct_int_state',&
        &            'allocation for cell_environ%area_norm failed')
      ENDIF

      ptr_int%cell_environ%nmbr_nghbr_cells(:,:) = 0 
      ptr_int%cell_environ%idx(:,:,:) = 0
      ptr_int%cell_environ%blk(:,:,:) = 0
      ptr_int%cell_environ%area_norm(:,:,:) = 0.0_wp

    END IF

  ENDIF


  !
  ! pos_on_tplane_e
  !
  ALLOCATE (ptr_int%pos_on_tplane_e(nproma, 4, 2, nblks_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
    &            'allocation for pos_on_tplane_e failed')
  ENDIF

  IF (ptr_patch%geometry_info%cell_type == 3) THEN
    !
    ! pos_on_tplane_c_edge
    !
    ALLOCATE (ptr_int%pos_on_tplane_c_edge(nproma, nblks_e, 2, 5), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',&
      &            'allocation for pos_on_tplane_c_edge failed')
    ENDIF
  ENDIF

  !
  ! Least squares reconstruction
  !
  ! *** linear ***
  !
  !
  ! lsq_dim_stencil
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_dim_stencil(nproma, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_dim_stencil failed')
  ENDIF
  !
  ! lsq_idx_c
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_idx_c(nproma, nblks_c, lsq_lin_set%dim_c),          &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_idx_c failed')
  ENDIF
  !
  ! lsq_blk_c
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_blk_c(nproma, nblks_c, lsq_lin_set%dim_c),          &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_blk_c failed')
  ENDIF
  !
  ! lsq_weights_c
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_weights_c(nproma, lsq_lin_set%dim_c, nblks_c),      &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_weights_c failed')
  ENDIF
  !
  ! lsq_qtmat_c
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_qtmat_c(nproma, lsq_lin_set%dim_unk, lsq_lin_set%dim_c, &
    &       nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_qtmat_c failed')
  ENDIF
  !
  ! lsq_rmat_rdiag_c
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_rmat_rdiag_c(nproma, lsq_lin_set%dim_unk, nblks_c), &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_rmat_rdiag_c failed')
  ENDIF
  !
  ! lsq_rmat_utri_c
  !
  idummy=(lsq_lin_set%dim_unk*lsq_lin_set%dim_unk - lsq_lin_set%dim_unk)/2
  ALLOCATE (ptr_int%lsq_lin%lsq_rmat_utri_c(nproma, idummy, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_rmat_utri_c failed')
  ENDIF
  !
  ! lsq_pseudoinv
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_pseudoinv(nproma, lsq_lin_set%dim_unk,              &
    &       lsq_lin_set%dim_c, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                            &
      &             'allocation for lsq_pseudoinv failed')
  ENDIF
  !
  ! lsq_moments
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_moments(nproma, nblks_c, lsq_lin_set%dim_unk),      &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                            &
      &             'allocation for lsq_moments failed')
  ENDIF
  !
  ! lsq_moments_hat
  !
  ALLOCATE (ptr_int%lsq_lin%lsq_moments_hat(nproma, nblks_c, lsq_lin_set%dim_c,     &
    &       lsq_lin_set%dim_unk), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                            &
      &             'allocation for lsq_moments_hat failed')
  ENDIF

  ! *** higher order ***
  !
  !
  ! lsq_dim_stencil
  !
  ALLOCATE (ptr_int%lsq_high%lsq_dim_stencil(nproma, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_dim_stencil failed')
  ENDIF
  !
  ! lsq_idx_c
  !
  ALLOCATE (ptr_int%lsq_high%lsq_idx_c(nproma, nblks_c, lsq_high_set%dim_c),        &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_idx_c failed')
  ENDIF
  !
  ! lsq_blk_c
  !
  ALLOCATE (ptr_int%lsq_high%lsq_blk_c(nproma, nblks_c, lsq_high_set%dim_c),        &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_blk_c failed')
  ENDIF
  !
  ! lsq_weights_c
  !
  ALLOCATE (ptr_int%lsq_high%lsq_weights_c(nproma, lsq_high_set%dim_c, nblks_c),    &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_weights_c failed')
  ENDIF
  !
  ! lsq_qtmat_c
  !
  ALLOCATE (ptr_int%lsq_high%lsq_qtmat_c(nproma, lsq_high_set%dim_unk, lsq_high_set%dim_c, &
    &       nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_qtmat_c failed')
  ENDIF
  !
  ! lsq_rmat_rdiag_c
  !
  ALLOCATE (ptr_int%lsq_high%lsq_rmat_rdiag_c(nproma, lsq_high_set%dim_unk, nblks_c), &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_rmat_rdiag_c failed')
  ENDIF
  !
  ! lsq_rmat_utri_c
  !
  idummy=(lsq_high_set%dim_unk*lsq_high_set%dim_unk - lsq_high_set%dim_unk)/2
  ALLOCATE (ptr_int%lsq_high%lsq_rmat_utri_c(nproma, idummy, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',&
      &          'allocation for lsq_rmat_utri_c failed')
  ENDIF
  !
  ! lsq_pseudoinv
  !
  ALLOCATE (ptr_int%lsq_high%lsq_pseudoinv(nproma, lsq_high_set%dim_unk,            &
    &       lsq_high_set%dim_c, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                            &
      &             'allocation for lsq_pseudoinv failed')
  ENDIF
  !
  ! lsq_moments
  !
  ALLOCATE (ptr_int%lsq_high%lsq_moments(nproma, nblks_c, lsq_high_set%dim_unk),    &
    &       STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                            &
      &             'allocation for lsq_moments failed')
  ENDIF
  !
  ! lsq_moments_hat
  !
  ALLOCATE (ptr_int%lsq_high%lsq_moments_hat(nproma, nblks_c, lsq_high_set%dim_c,   &
    &       lsq_high_set%dim_unk), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                            &
      &             'allocation for lsq_moments_hat failed')
  ENDIF


  IF (ptr_patch%geometry_info%cell_type == 3) THEN
    ALLOCATE (ptr_int%geofac_grdiv(nproma, 5, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for geofac_grdiv failed')
    ENDIF
    ALLOCATE (ptr_int%nudgecoeff_c(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for nudgecoeff_c failed')
    ENDIF
    ALLOCATE (ptr_int%nudgecoeff_e(nproma, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for nudgecoeff_e failed')
    ENDIF

    !
    ! Quadrature points and weights for integration over triangular element
    !
    ALLOCATE (ptr_int%gquad%qpts_tri_l(nproma, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_l failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%qpts_tri_q(nproma, nblks_c,3), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_q failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%qpts_tri_c(nproma, nblks_c,4), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for qpts_tri_c failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%weights_tri_q(3), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for weights_tri_q failed')
    ENDIF
    ALLOCATE (ptr_int%gquad%weights_tri_c(4), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:construct_int_state',                     &
        &             'allocation for weights_tri_c failed')
    ENDIF
  ENDIF

  ALLOCATE (ptr_int%geofac_div(nproma, ptr_patch%geometry_info%cell_type, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_div failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_rot(nproma, 9-ptr_patch%geometry_info%cell_type, nblks_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state', &
    &             'allocation for geofac_rot failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_n2s(nproma, ptr_patch%geometry_info%cell_type+1, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_n2s failed')
  ENDIF

  ALLOCATE (ptr_int%geofac_grg(nproma, ptr_patch%geometry_info%cell_type+1, nblks_c, 2), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for geofac_grg failed')
  ENDIF

  ALLOCATE (ptr_int%primal_normal_ec(nproma, nblks_c,ptr_patch%geometry_info%cell_type, 2), STAT=ist)
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for primal_normal_ec failed')
  ENDIF

  ALLOCATE (ptr_int%edge_cell_length(nproma, nblks_c, ptr_patch%geometry_info%cell_type), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for edge_cell_length failed')
  ENDIF

  ALLOCATE (ptr_int%cell_vert_dist(nproma, ptr_patch%geometry_info%cell_type, 2, nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                       &
      &             'allocation for cell_vert_dist failed')
  ENDIF


  !
  ! initialize all components
  !
  IF (ptr_patch%geometry_info%cell_type == 3 ) THEN
    ptr_int%e_bln_c_s     = 0._wp
    ptr_int%e_bln_c_u     = 0._wp
    ptr_int%e_bln_c_v     = 0._wp
    ptr_int%c_bln_avg     = 0._wp
    ptr_int%gradc_bmat    = 0._wp
    ptr_int%e_flx_avg     = 0._wp
  ENDIF

  ptr_int%v_1o2_e          = 0.5_wp
  ptr_int%c_lin_e          = 0._wp
  ptr_int%e_inn_c          = 0._wp
  ptr_int%verts_aw_cells   = 0._wp
  ptr_int%cells_aw_verts   = 0._wp

  IF( ptr_patch%geometry_info%cell_type == 3) THEN
    ptr_int%rbf_vec_idx_c     = 0
    ptr_int%rbf_vec_blk_c     = 0
    ptr_int%rbf_vec_stencil_c = 0
    ptr_int%rbf_vec_coeff_c   = 0._wp

    ptr_int%rbf_c2grad_idx    = 0
    ptr_int%rbf_c2grad_blk    = 0
    ptr_int%rbf_c2grad_coeff  = 0._wp

    ptr_int%rbf_vec_idx_v     = 0
    ptr_int%rbf_vec_blk_v     = 0
    ptr_int%rbf_vec_stencil_v = 0
    ptr_int%rbf_vec_coeff_v   = 0._wp

    ptr_int%rbf_vec_idx_e     = 0
    ptr_int%rbf_vec_blk_e     = 0
    ptr_int%rbf_vec_stencil_e = 0
    ptr_int%rbf_vec_coeff_e   = 0._wp
  ENDIF

  ptr_int%pos_on_tplane_e     = 0._wp

  IF (ptr_patch%geometry_info%cell_type == 3) THEN
    ptr_int%pos_on_tplane_c_edge(:,:,:,:)%lon = 0._wp
    ptr_int%pos_on_tplane_c_edge(:,:,:,:)%lat = 0._wp
  ENDIF

  ptr_int%lsq_lin%lsq_dim_stencil   = 0
  ptr_int%lsq_lin%lsq_idx_c         = 0
  ptr_int%lsq_lin%lsq_blk_c         = 0
  ptr_int%lsq_lin%lsq_weights_c     = 0._wp
  ptr_int%lsq_lin%lsq_qtmat_c       = 0._wp
  ptr_int%lsq_lin%lsq_rmat_rdiag_c  = 0._wp
  ptr_int%lsq_lin%lsq_rmat_utri_c   = 0._wp
  ptr_int%lsq_lin%lsq_pseudoinv     = 0._wp
  ptr_int%lsq_lin%lsq_moments       = 0._wp
  ptr_int%lsq_lin%lsq_moments_hat   = 0._wp

  ptr_int%lsq_high%lsq_dim_stencil  = 0
  ptr_int%lsq_high%lsq_idx_c        = 0
  ptr_int%lsq_high%lsq_blk_c        = 0
  ptr_int%lsq_high%lsq_weights_c    = 0._wp
  ptr_int%lsq_high%lsq_qtmat_c      = 0._wp
  ptr_int%lsq_high%lsq_rmat_rdiag_c = 0._wp
  ptr_int%lsq_high%lsq_rmat_utri_c  = 0._wp
  ptr_int%lsq_high%lsq_pseudoinv    = 0._wp
  ptr_int%lsq_high%lsq_moments      = 0._wp
  ptr_int%lsq_high%lsq_moments_hat  = 0._wp

  IF (ptr_patch%geometry_info%cell_type ==3) THEN
    ptr_int%geofac_grdiv = 0._wp
    ptr_int%nudgecoeff_c = 0._wp
    ptr_int%nudgecoeff_e = 0._wp

    ptr_int%gquad%qpts_tri_l(:,:)%lat    = 0._wp
    ptr_int%gquad%qpts_tri_l(:,:)%lon    = 0._wp
    ptr_int%gquad%qpts_tri_q(:,:,:)%lat  = 0._wp
    ptr_int%gquad%qpts_tri_q(:,:,:)%lon  = 0._wp
    ptr_int%gquad%qpts_tri_c(:,:,:)%lat  = 0._wp
    ptr_int%gquad%qpts_tri_c(:,:,:)%lon  = 0._wp
    ptr_int%gquad%weights_tri_q(:)       = 0._wp
    ptr_int%gquad%weights_tri_c(:)       = 0._wp
  ENDIF
  ptr_int%geofac_div = 0._wp
  ptr_int%geofac_rot = 0._wp
  ptr_int%geofac_n2s = 0._wp
  ptr_int%geofac_grg = 0._wp
  ptr_int%primal_normal_ec = 0._wp
  ptr_int%edge_cell_length = 0._wp
  ptr_int%cell_vert_dist = 0._wp

  CALL message ('mo_intp_state:allocate_int_state','memory allocation finished')

END SUBROUTINE allocate_int_state


!-------------------------------------------------------------------------
!
!
!! Allocation of components of interpolation state.
!!
!! Initialization of components.
!!
SUBROUTINE construct_2d_interpol_state(ptr_patch, ptr_int_state)
!

TYPE(t_patch), INTENT(INOUT) :: ptr_patch(n_dom_start:)

TYPE(t_int_state),     INTENT(INOUT) :: ptr_int_state(n_dom_start:)

INTEGER  :: jg
REAL(wp) :: search_radius

CHARACTER(len=MAX_CHAR_LENGTH) :: text

  !-----------------------------------------------------------------------

CALL message('mo_intp_state:construct_2d_interpol_state','start to construct int_state')

DO jg = n_dom_start, n_dom

  WRITE(text,'(a,i0)') 'constructing int_state for patch ',jg
  CALL message('mo_intp_state:construct_2d_interpol_state',text)

  CALL allocate_int_state( ptr_patch(jg), ptr_int_state(jg))

  !
  ! initializion of coefficients for averaging of scalars and kinetic energy
  !
  CALL scalar_int_coeff(ptr_patch(jg), ptr_int_state(jg))

  ! Initialization of coefficients for bilinear cell averaging, divergence, rotation
  ! and nabla_2_scalar; transformation of edge orientation vectors to the locations
  ! of cells and vertices; computation of coefficients for bilinear edge-to-cell
  ! interpolation

  CALL complete_patchinfo( ptr_patch(jg), ptr_int_state(jg))
  CALL init_geo_factors(ptr_patch(jg), ptr_int_state(jg))
  IF (ptr_patch(jg)%geometry_info%cell_type==3)THEN

    ! CALL of init_cellavg_wgt is skipped for the wave model
    !
    IF (.NOT.my_process_is_waves()) THEN
      !
      ! not needed for the wave model.
      ! Running this routine for the wave model with the NAG compiler
      ! results in a floating point exception in calculate_bilinear_cellavg_wgt.
      CALL init_cellavg_wgt(ptr_patch(jg), ptr_int_state(jg))
    ENDIF

    CALL bln_int_coeff_e2c( ptr_patch(jg), ptr_int_state(jg) )
    IF(jg>0) THEN
      IF (l_limited_area .AND. jg == 1 .OR. jg > 1) THEN
        CALL init_nudgecoeffs( ptr_patch(jg), ptr_int_state(jg) )
      ENDIF
    ENDIF
  ENDIF

  !
  ! initialization of indices and coefficients for vector rbf interpolation
  !
  IF (ptr_patch(jg)%geometry_info%cell_type == 3) THEN

    ! ... at cell centers
    CALL rbf_vec_index_cell (ptr_patch(jg), ptr_int_state(jg))
    !
    CALL rbf_vec_compute_coeff_cell (ptr_patch(jg), ptr_int_state(jg))
    !
    ! ... at triangle vertices
    CALL rbf_vec_index_vertex (ptr_patch(jg), ptr_int_state(jg))
    !
    CALL rbf_vec_compute_coeff_vertex (ptr_patch(jg), ptr_int_state(jg))
    !
    ! ... at edge midpoints
    CALL rbf_vec_index_edge (ptr_patch(jg), ptr_int_state(jg))
    !
    CALL rbf_vec_compute_coeff_edge (ptr_patch(jg), ptr_int_state(jg))
    !
    ! Compute coefficients needed for gradient reconstruction at cell midpoints
    CALL rbf_c2grad_index (ptr_patch(jg), ptr_int_state(jg))
    CALL rbf_compute_coeff_c2grad (ptr_patch(jg), ptr_int_state(jg))


    IF ( ptr_int_state(jg)%cell_environ%is_used ) THEN
      ! Determine ptr_int_state(..)%cell_environ
      search_radius = 1.0e20_wp  ! here, the search radius shall not be used -> choose a large enough value >> dx
      CALL gen_index_list_radius( ptr_int_state(jg), ptr_patch(jg), jg, search_radius, 1 )
    END IF

    ! initialization of quadrature points and weights
    !
    CALL tri_quadrature_pts (ptr_patch(jg), ptr_int_state(jg))
  ENDIF

  !
  ! - Initialization of a tangential plane at edge midpoints for the calculation
  !   of backward trajectories.
  ! - Initialization of a tangential plane at cell centers
  ! - stencil generation
  ! - initialization of coefficients for least squares polynomial
  !   reconstruction at cell centers
  !
  IF (.NOT. lplane) THEN

    CALL init_tplane_e(ptr_patch(jg), ptr_int_state(jg))

    CALL init_tplane_c(ptr_patch(jg), ptr_int_state(jg))

    ! 3-point stencil
    CALL lsq_stencil_create( ptr_patch(jg), ptr_int_state(jg)%lsq_lin,      &
      &                      lsq_lin_set%dim_c )
    CALL lsq_compute_coeff_cell( ptr_patch(jg), ptr_int_state(jg)%lsq_lin,  &
      &                      lsq_lin_set%l_consv, lsq_lin_set%dim_c,        &
      &                      lsq_lin_set%dim_unk, lsq_lin_set%wgt_exp )

    ! 9 or 12-point stencil for higher order reconstruction
    CALL lsq_stencil_create( ptr_patch(jg), ptr_int_state(jg)%lsq_high,     &
      &                   lsq_high_set%dim_c )
    CALL lsq_compute_coeff_cell( ptr_patch(jg), ptr_int_state(jg)%lsq_high, &
      &                       lsq_high_set%l_consv, lsq_high_set%dim_c,     &
      &                       lsq_high_set%dim_unk, lsq_high_set%wgt_exp )
  ENDIF


  ! SCM initialization of RBF coefficients
  !
  IF ( l_scm_mode .and. (.not.ldynamics) ) THEN
    CALL rbf_coeff_scm( ptr_patch(jg), ptr_int_state(jg) )
  ENDIF

ENDDO

#ifdef SERIALIZE
IF (ser_rbf) THEN
  CALL ser_rbf_coefficients(ptr_int_state)
ENDIF
#endif

CALL message('mo_intp_state:construct_2d_interpol_state', &
  & 'construction of interpolation state finished')

END SUBROUTINE construct_2d_interpol_state

!-------------------------------------------------------------------------
!> xfer_var family: transfer variables from parent to local parent
!-------------------------------------------------------------------------

SUBROUTINE xfer_var_r2(typ, pos_nproma, pos_nblks, p_p, p_lp, arri, arro)

  INTEGER, INTENT(IN) :: typ, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  REAL(wp), INTENT(IN)    :: arri(:,:)
  REAL(wp), INTENT(INOUT) :: arro(:,:)
  ! local variables

  IF(typ == SYNC_C) THEN
    CALL exchange_data(p_pat=comm_pat_glb_to_loc_c, lacc=.false., RECV=arro, SEND=arri)
  ELSEIF(typ == SYNC_E) THEN
    CALL exchange_data(p_pat=comm_pat_glb_to_loc_e, lacc=.false., RECV=arro, SEND=arri)
  ELSEIF(typ == SYNC_V) THEN
    CALL exchange_data(p_pat=comm_pat_glb_to_loc_v, lacc=.false., RECV=arro, SEND=arri)
  ELSE
    CALL finish ('mo_interpolation:xfer_var','Illegal type for sync')
  ENDIF

END SUBROUTINE xfer_var_r2

!-------------------------------------------------------------------------

SUBROUTINE xfer_var_r3(typ, pos_nproma, pos_nblks, p_p, p_lp, arri, arro)

  INTEGER, INTENT(IN) :: typ, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  REAL(wp), INTENT(IN)    :: arri(:,:,:)
  REAL(wp), INTENT(INOUT) :: arro(:,:,:)

  INTEGER :: j

  IF(pos_nproma==1 .AND. pos_nblks==2) THEN
    DO j = 1, UBOUND(arri,3)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,:,j),arro(:,:,j))
    ENDDO
  ELSEIF(pos_nproma==1 .AND. pos_nblks==3) THEN
    DO j = 1, UBOUND(arri,2)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,j,:),arro(:,j,:))
    ENDDO
  ELSEIF(pos_nproma==2 .AND. pos_nblks==3) THEN
    DO j = 1, UBOUND(arri,1)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(j,:,:),arro(j,:,:))
    ENDDO
  ELSE
    CALL finish ('mo_interpolation:xfer_var','Illegal value for pos_nproma/pos_nblks')
  ENDIF

END SUBROUTINE xfer_var_r3

!-------------------------------------------------------------------------

SUBROUTINE xfer_var_r4(typ, pos_nproma, pos_nblks, p_p, p_lp, arri, arro)

  INTEGER, INTENT(IN) :: typ, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  REAL(wp), INTENT(IN)    :: arri(:,:,:,:)
  REAL(wp), INTENT(INOUT) :: arro(:,:,:,:)

  INTEGER :: i,j

  ! Variable has 4 dimensions
  IF(pos_nproma == 1 .AND. pos_nblks == 2)  THEN
    DO j = 1, UBOUND(arri,4)
    DO i = 1, UBOUND(arri,3)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,:,i,j),arro(:,:,i,j))
    ENDDO
    ENDDO
  ELSEIF(pos_nproma == 1 .AND. pos_nblks == 3)  THEN
    DO j = 1, UBOUND(arri,4)
    DO i = 1, UBOUND(arri,2)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,i,:,j),arro(:,i,:,j))
    ENDDO
    ENDDO
  ELSEIF(pos_nproma == 1 .AND. pos_nblks == 4)  THEN
    DO j = 1, UBOUND(arri,3)
    DO i = 1, UBOUND(arri,2)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(:,i,j,:),arro(:,i,j,:))
    ENDDO
    ENDDO
  ELSEIF(pos_nproma == 3 .AND. pos_nblks == 4)  THEN
    DO j = 1, UBOUND(arri,2)
    DO i = 1, UBOUND(arri,1)
      CALL xfer_var_r2(typ,1,2,p_p,p_lp,arri(i,j,:,:),arro(i,j,:,:))
    ENDDO
    ENDDO
  ELSE
    ! Other pos_nproma/pos_nblks combinations are possible but currently not existing!
    CALL finish ('mo_interpolation:xfer_var','unsupported value for pos_nproma/pos_nblks')
  ENDIF


END SUBROUTINE xfer_var_r4
!-------------------------------------------------------------------------

SUBROUTINE xfer_var_i2(typ, pos_nproma, pos_nblks, p_p, p_lp, arri, arro)

  INTEGER, INTENT(IN) :: typ, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  INTEGER,INTENT(IN) :: arri(:,:)
  INTEGER,INTENT(INOUT) :: arro(:,:)

  REAL(wp) :: r_arri(UBOUND(arri,1),UBOUND(arri,2))
  REAL(wp) :: r_arro(UBOUND(arro,1),UBOUND(arro,2))

  r_arri(:,:) = REAL(arri,wp)
  r_arro(:,:) = 0._wp ! Safety only
  CALL xfer_var_r2(typ,pos_nproma,pos_nblks,p_p,p_lp,r_arri,r_arro)
  arro(:,:) = INT(r_arro)

END SUBROUTINE xfer_var_i2

!-------------------------------------------------------------------------
!> xfer_idx family: transfer index variables from parent to local parent
!-------------------------------------------------------------------------

SUBROUTINE xfer_idx_2(type_arr, type_idx, pos_nproma, pos_nblks, p_p, p_lp, idxi, blki, idxo, blko)

  INTEGER, INTENT(IN) :: type_arr, type_idx, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN), TARGET :: p_p, p_lp
  INTEGER, INTENT(IN)    :: idxi(:,:), blki(:,:)
  INTEGER, INTENT(INOUT) :: idxo(:,:), blko(:,:)

  INTEGER :: jb, jl, i_l, i_g, n_idx_l, n_idx_g
  REAL(wp) :: z_idxi(UBOUND(idxi,1),UBOUND(idxi,2))
  REAL(wp) :: z_idxo(UBOUND(idxo,1),UBOUND(idxo,2))
  INTEGER, POINTER :: glb_index(:)
  TYPE(t_grid_domain_decomp_info), POINTER :: local_decomp_info

  IF(type_idx == SYNC_C) THEN
    glb_index => p_p%cells%decomp_info%glb_index
    local_decomp_info => p_lp%cells%decomp_info
    n_idx_l = p_p%n_patch_cells
    n_idx_g = p_p%n_patch_cells_g
  ELSEIF(type_idx == SYNC_E) THEN
    glb_index => p_p%edges%decomp_info%glb_index
    local_decomp_info => p_lp%edges%decomp_info
    n_idx_l = p_p%n_patch_edges
    n_idx_g = p_p%n_patch_edges_g
  ELSEIF(type_idx == SYNC_V) THEN
    glb_index => p_p%verts%decomp_info%glb_index
    local_decomp_info => p_lp%verts%decomp_info
    n_idx_l = p_p%n_patch_verts
    n_idx_g = p_p%n_patch_verts_g
  ELSE
    CALL finish('xfer_idx','Unsupported type_idx')
  ENDIF

  z_idxi(:,:) = 0._wp
  z_idxo(:,:) = 0._wp

  DO jb = 1, UBOUND(idxi,2)
    DO jl = 1, nproma

      i_l = idx_1d(idxi(jl,jb),blki(jl,jb))

      IF(i_l <= 0 .OR. i_l > n_idx_l) THEN
        z_idxi(jl,jb) = 0._wp
      ELSE
        z_idxi(jl,jb) = glb_index(i_l)
      ENDIF

    END DO
  END DO

  CALL xfer_var_r2(type_arr,pos_nproma,pos_nblks,p_p,p_lp,z_idxi,z_idxo)

  DO jb = 1, UBOUND(idxo,2)
    DO jl = 1, nproma

      i_g = INT(z_idxo(jl,jb))

      IF(i_g <= 0 .OR. i_g > n_idx_g) THEN
        idxo(jl,jb) = 0
        blko(jl,jb) = 0
      ELSE
        i_l = get_valid_local_index(local_decomp_info%glb2loc_index, i_g)
        idxo(jl,jb) = idx_no(i_l)
        blko(jl,jb) = blk_no(i_l)
      ENDIF

    END DO
  END DO

END SUBROUTINE xfer_idx_2

!-------------------------------------------------------------------------

SUBROUTINE xfer_idx_3(type_arr, type_idx, pos_nproma, pos_nblks, p_p, p_lp, idxi, blki, idxo, blko)

  INTEGER, INTENT(IN) :: type_arr, type_idx, pos_nproma, pos_nblks
  TYPE(t_patch), INTENT(IN) :: p_p, p_lp
  INTEGER, INTENT(IN)    :: idxi(:,:,:), blki(:,:,:)
  INTEGER, INTENT(INOUT) :: idxo(:,:,:), blko(:,:,:)

  INTEGER :: j

  IF(pos_nproma==1 .AND. pos_nblks==2) THEN
    DO j = 1, UBOUND(idxi,3)
      CALL xfer_idx_2(type_arr,type_idx,1,2,p_p,p_lp,idxi(:,:,j),blki(:,:,j), &
                                                   & idxo(:,:,j),blko(:,:,j))
    ENDDO
  ELSEIF(pos_nproma==1 .AND. pos_nblks==3) THEN
    DO j = 1, UBOUND(idxi,2)
      CALL xfer_idx_2(type_arr,type_idx,1,2,p_p,p_lp,idxi(:,j,:),blki(:,j,:), &
                                                   & idxo(:,j,:),blko(:,j,:))
    ENDDO
  ELSEIF(pos_nproma==2 .AND. pos_nblks==3) THEN
    DO j = 1, UBOUND(idxi,1)
      CALL xfer_idx_2(type_arr,type_idx,1,2,p_p,p_lp,idxi(j,:,:),blki(j,:,:), &
                                                   & idxo(j,:,:),blko(j,:,:))
    ENDDO
  ELSE
    CALL finish ('mo_interpolation:xfer_idx','Illegal value for pos_nproma/pos_nblks')
  ENDIF

END SUBROUTINE xfer_idx_3

!-------------------------------------------------------------------------
!! Transfers interpolation state from parent to local parent
!-------------------------------------------------------------------------
!
SUBROUTINE transfer_interpol_state(p_p, p_lp, pi, po)
!
  TYPE(t_patch), INTENT(IN)    :: p_p   ! parent
  TYPE(t_patch), INTENT(INOUT) :: p_lp  ! local parent

  TYPE(t_int_state), INTENT(IN)    :: pi ! Interpolation state on parent
  TYPE(t_int_state), INTENT(INOUT) :: po ! Interpolation state on local parent

  INTEGER, ALLOCATABLE :: owner(:)

  ! Allocate interpolation state for local parent

  CALL allocate_int_state(p_lp, po)

  ! Set up communication patterns for transferring the data to local parents.
  ! Since these communication patterns are not used elsewhere, they are
  ! stored locally and deleted at the end of the routine

  ALLOCATE(owner(MAX(p_lp%n_patch_cells, p_lp%n_patch_verts, &
    &                p_lp%n_patch_edges)))

  owner(1:p_lp%n_patch_cells) = &
    dist_dir_get_owners(p_p%cells%decomp_info%owner_dist_dir, &
      &                 p_lp%cells%decomp_info%glb_index(1:p_lp%n_patch_cells))
  CALL setup_comm_pattern(p_lp%n_patch_cells, owner(1:p_lp%n_patch_cells), &
    &                     p_lp%cells%decomp_info%glb_index,  &
    &                     p_p%cells%decomp_info%glb2loc_index, &
    &                     p_p%n_patch_cells, &
    &                     p_p%cells%decomp_info%owner_local, &
    &                     p_p%cells%decomp_info%glb_index, &
    &                     comm_pat_glb_to_loc_c)

  owner(1:p_lp%n_patch_edges) = &
    dist_dir_get_owners(p_p%edges%decomp_info%owner_dist_dir, &
      &                 p_lp%edges%decomp_info%glb_index(1:p_lp%n_patch_edges))
  CALL setup_comm_pattern(p_lp%n_patch_edges, owner(1:p_lp%n_patch_edges), &
    &                     p_lp%edges%decomp_info%glb_index,  &
    &                     p_p%edges%decomp_info%glb2loc_index, &
    &                     p_p%n_patch_edges, &
    &                     p_p%edges%decomp_info%owner_local, &
    &                     p_p%edges%decomp_info%glb_index, &
    &                     comm_pat_glb_to_loc_e)

  owner(1:p_lp%n_patch_verts) = &
    dist_dir_get_owners(p_p%verts%decomp_info%owner_dist_dir, &
      &                 p_lp%verts%decomp_info%glb_index(1:p_lp%n_patch_verts))
  CALL setup_comm_pattern(p_lp%n_patch_verts, owner(1:p_lp%n_patch_verts), &
    &                     p_lp%verts%decomp_info%glb_index,  &
    &                     p_p%verts%decomp_info%glb2loc_index, &
    &                     p_p%n_patch_verts, &
    &                     p_p%verts%decomp_info%owner_local, &
    &                     p_p%verts%decomp_info%glb_index, &
    &                     comm_pat_glb_to_loc_v)

  DEALLOCATE(owner)

  ! Some edge related values of the patch are only set in
  ! construct_2d_interpol_state (complete_patchinfo) and have to
  ! be set here also

  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p%edges%area_edge,p_lp%edges%area_edge)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%inv_primal_edge_length, &
                                  & p_lp%edges%inv_primal_edge_length)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%inv_dual_edge_length, &
                                  & p_lp%edges%inv_dual_edge_length)
  IF (p_p%geometry_info%cell_type == 3) THEN
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%inv_vert_vert_length, &
                                  & p_lp%edges%inv_vert_vert_length)
  ENDIF
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%primal_normal_cell(:,:,:)%v1, &
                                  & p_lp%edges%primal_normal_cell(:,:,:)%v1)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%primal_normal_cell(:,:,:)%v2, &
                                  & p_lp%edges%primal_normal_cell(:,:,:)%v2)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%primal_normal_vert(:,:,:)%v1, &
                                  & p_lp%edges%primal_normal_vert(:,:,:)%v1)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%primal_normal_vert(:,:,:)%v2, &
                                  & p_lp%edges%primal_normal_vert(:,:,:)%v2)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%dual_normal_cell(:,:,:)%v1, &
                                  & p_lp%edges%dual_normal_cell(:,:,:)%v1)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%dual_normal_cell(:,:,:)%v2, &
                                  & p_lp%edges%dual_normal_cell(:,:,:)%v2)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%dual_normal_vert(:,:,:)%v1, &
                                  & p_lp%edges%dual_normal_vert(:,:,:)%v1)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%dual_normal_vert(:,:,:)%v2, &
                                  & p_lp%edges%dual_normal_vert(:,:,:)%v2)

  CALL xfer_idx(SYNC_E,SYNC_V,1,2,p_p,p_lp,p_p %edges%vertex_idx,p_p %edges%vertex_blk, &
                                         & p_lp%edges%vertex_idx,p_lp%edges%vertex_blk)

  ! The same for quad_area, quad_orientation, quad_idx which is calculated
  ! after the complete patch has been read

  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%quad_area, &
                                  & p_lp%edges%quad_area)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,p_p %edges%quad_orientation, &
                                  & p_lp%edges%quad_orientation)

  CALL xfer_idx(SYNC_E,SYNC_E,1,2,p_p,p_lp,p_p %edges%quad_idx,p_p %edges%quad_blk, &
                                         & p_lp%edges%quad_idx,p_lp%edges%quad_blk)


  ! Transfer interpolation state

  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%c_lin_e,po%c_lin_e)
  IF (p_p%geometry_info%cell_type == 3) THEN
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_bln_c_s,po%e_bln_c_s)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_bln_c_u,po%e_bln_c_u)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_bln_c_v,po%e_bln_c_v)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%c_bln_avg,po%c_bln_avg)
  CALL xfer_var(SYNC_C,1,4,p_p,p_lp,pi%gradc_bmat,po%gradc_bmat)
  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%e_flx_avg,po%e_flx_avg)
  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%v_1o2_e,po%v_1o2_e)
  ENDIF
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%e_inn_c,po%e_inn_c)

  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%verts_aw_cells,po%verts_aw_cells)
  CALL xfer_var(SYNC_V,1,3,p_p,p_lp,pi%cells_aw_verts,po%cells_aw_verts)

  IF (p_p%geometry_info%cell_type == 3) THEN
  CALL xfer_idx(SYNC_C,SYNC_E,2,3,p_p,p_lp,pi%rbf_vec_idx_c,pi%rbf_vec_blk_c, &
                                         & po%rbf_vec_idx_c,po%rbf_vec_blk_c)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%rbf_vec_stencil_c,po%rbf_vec_stencil_c)
  CALL xfer_var(SYNC_C,3,4,p_p,p_lp,pi%rbf_vec_coeff_c,po%rbf_vec_coeff_c)
  CALL xfer_idx(SYNC_C,SYNC_C,2,3,p_p,p_lp,pi%rbf_c2grad_idx,pi%rbf_c2grad_blk, &
                                         & po%rbf_c2grad_idx,po%rbf_c2grad_blk)

  IF ( pi%cell_environ%is_used ) THEN
    CALL xfer_var(SYNC_C,1,2,p_p,p_lp, pi%cell_environ%nmbr_nghbr_cells, &
      &                                po%cell_environ%nmbr_nghbr_cells)
    CALL xfer_idx(SYNC_C,SYNC_C,1,2,p_p,p_lp, pi%cell_environ%idx, pi%cell_environ%blk, &
      &                                       po%cell_environ%idx, po%cell_environ%blk)
    CALL xfer_var(SYNC_C,1,2,p_p,p_lp, pi%cell_environ%area_norm, &
      &                                po%cell_environ%area_norm)
    po%cell_environ%is_used              = pi%cell_environ%is_used
    po%cell_environ%max_nmbr_nghbr_cells = pi%cell_environ%max_nmbr_nghbr_cells
    po%cell_environ%radius               = pi%cell_environ%radius
    po%cell_environ%max_nmbr_iter        = pi%cell_environ%max_nmbr_iter
  END IF

  CALL xfer_var(SYNC_C,3,4,p_p,p_lp,pi%rbf_c2grad_coeff,po%rbf_c2grad_coeff)
  CALL xfer_idx(SYNC_V,SYNC_E,2,3,p_p,p_lp,pi%rbf_vec_idx_v,pi%rbf_vec_blk_v, &
                                         & po%rbf_vec_idx_v,po%rbf_vec_blk_v)
  CALL xfer_var(SYNC_V,1,2,p_p,p_lp,pi%rbf_vec_stencil_v,po%rbf_vec_stencil_v)
  CALL xfer_var(SYNC_V,3,4,p_p,p_lp,pi%rbf_vec_coeff_v,po%rbf_vec_coeff_v)
  CALL xfer_idx(SYNC_E,SYNC_E,2,3,p_p,p_lp,pi%rbf_vec_idx_e,pi%rbf_vec_blk_e, &
                                         & po%rbf_vec_idx_e,po%rbf_vec_blk_e)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,pi%rbf_vec_stencil_e,po%rbf_vec_stencil_e)
  CALL xfer_var(SYNC_E,2,3,p_p,p_lp,pi%rbf_vec_coeff_e,po%rbf_vec_coeff_e)
  ENDIF

  IF (p_p%geometry_info%cell_type == 3) THEN
  CALL xfer_var(SYNC_E,1,3,p_p,p_lp,pi%geofac_grdiv,po%geofac_grdiv)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%nudgecoeff_c,po%nudgecoeff_c)
  CALL xfer_var(SYNC_E,1,2,p_p,p_lp,pi%nudgecoeff_e,po%nudgecoeff_e)

  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_l(:,:)%lon,po%gquad%qpts_tri_l(:,:)%lon)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_l(:,:)%lat,po%gquad%qpts_tri_l(:,:)%lat)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_q(:,:,:)%lon,po%gquad%qpts_tri_q(:,:,:)%lon)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_q(:,:,:)%lat,po%gquad%qpts_tri_q(:,:,:)%lat)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_c(:,:,:)%lon,po%gquad%qpts_tri_c(:,:,:)%lon)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%gquad%qpts_tri_c(:,:,:)%lat,po%gquad%qpts_tri_c(:,:,:)%lat)

  po%gquad%weights_tri_q(:) = pi%gquad%weights_tri_q(:)
  po%gquad%weights_tri_c(:) = pi%gquad%weights_tri_c(:)

  ENDIF
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%geofac_div,po%geofac_div)
  CALL xfer_var(SYNC_V,1,3,p_p,p_lp,pi%geofac_rot,po%geofac_rot)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%geofac_n2s,po%geofac_n2s)
  CALL xfer_var(SYNC_C,1,3,p_p,p_lp,pi%geofac_grg,po%geofac_grg)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%primal_normal_ec,po%primal_normal_ec)
  CALL xfer_var(SYNC_C,1,2,p_p,p_lp,pi%edge_cell_length,po%edge_cell_length)
  CALL xfer_var(SYNC_C,1,4,p_p,p_lp,pi%cell_vert_dist,po%cell_vert_dist)

  ! clean up

!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_glb_to_loc_c)
!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_glb_to_loc_e)
!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_glb_to_loc_v)

END SUBROUTINE transfer_interpol_state
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!! Deallocation of components of a single 2d interpolation state.
!!
SUBROUTINE deallocate_int_state( ptr_int)
!
TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER :: ist

!-----------------------------------------------------------------------

  ! deallocate interpolation state
  !
  ! c_lin_e
  !
  DEALLOCATE (ptr_int%c_lin_e, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for c_lin_e failed')
  ENDIF
  !
  !
  ! e_bln_c_s
  !
  DEALLOCATE (ptr_int%e_bln_c_s, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_bln_c_s failed')
  ENDIF
  !
  ! e_bln_c_u
  !
  DEALLOCATE (ptr_int%e_bln_c_u, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_bln_c_u failed')
  ENDIF
  !
  ! e_bln_c_v
  !
  DEALLOCATE (ptr_int%e_bln_c_v, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_bln_c_v failed')
  ENDIF
  !
  ! c_bln_avg
  !
  DEALLOCATE (ptr_int%c_bln_avg, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for c_bln_avg failed')
  ENDIF
  !
  ! gradc_bmat
  !
  DEALLOCATE (ptr_int%gradc_bmat, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for gradc_bmat failed')
  ENDIF
  !
  ! e_flx_avg
  !
  DEALLOCATE (ptr_int%e_flx_avg, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_flx_avg failed')
  ENDIF
  !
  ! v_1o2_e
  !
  DEALLOCATE (ptr_int%v_1o2_e, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for v_1o2_e failed')
  ENDIF
  !
  ! e_inn_c
  !
  DEALLOCATE (ptr_int%e_inn_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for e_inn_c failed')
  ENDIF
  !
  !
  ! verts_aw_cells
  !
  DEALLOCATE (ptr_int%verts_aw_cells, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for verts_aw_cells failed')
  ENDIF
  !
  ! cells_aw_verts
  !
  DEALLOCATE (ptr_int%cells_aw_verts, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for cells_aw_verts failed')
  ENDIF

    !
    ! rbf_vec_idx_c, rbf_vec_blk_c
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_idx_c failed')
    ENDIF
    DEALLOCATE (ptr_int%rbf_vec_blk_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for rbf_vec_blk_c failed')
    ENDIF
    !
    ! rbf_vec_stencil_c
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_stencil_c failed')
    ENDIF
    !
    ! rbf_vec_coeff_c
    !
    DEALLOCATE (ptr_int%rbf_vec_coeff_c, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_coeff_c failed')
    ENDIF
    !
    ! rbf_c2grad_idx, rbf_c2grad_blk
    !
    DEALLOCATE (ptr_int%rbf_c2grad_idx, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_c2grad_idx failed')
    ENDIF
    DEALLOCATE (ptr_int%rbf_c2grad_blk, STAT=ist )
    IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for rbf_c2grad_blk failed')
    ENDIF
    !
    ! rbf_c2grad_coeff_c
    !
    DEALLOCATE (ptr_int%rbf_c2grad_coeff, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_c2grad_coeff failed')
    ENDIF
    !
    ! rbf_vec_idx_v, rbf_vec_blk_v
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_idx_v failed')
    ENDIF
    DEALLOCATE (ptr_int%rbf_vec_blk_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_idx_v failed')
    ENDIF
    !
    ! rbf_vec_stencil_v
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_stencil_v failed')
    ENDIF
    !
    ! rbf_vec_coeff_v
    !
    DEALLOCATE (ptr_int%rbf_vec_coeff_v, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_coeff_v failed')
    ENDIF
    !
    ! rbf_vec_idx_e, rbf_vec_blk_e
    !
    DEALLOCATE (ptr_int%rbf_vec_idx_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_idx_e failed')
    ENDIF
    DEALLOCATE (ptr_int%rbf_vec_blk_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_blk_e failed')
    ENDIF
    !
    ! rbf_vec_stencil_e
    !
    DEALLOCATE (ptr_int%rbf_vec_stencil_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
        &             'deallocation for rbf_vec_stencil_e failed')
    ENDIF
    !
    ! rbf_vec_coeff_e
    !
    DEALLOCATE (ptr_int%rbf_vec_coeff_e, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for rbf_vec_coeff_e failed')
    ENDIF

  IF ( ptr_int%cell_environ%is_used ) THEN
    !
    ! cell_environ
    !
    DEALLOCATE( ptr_int%cell_environ%nmbr_nghbr_cells, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',&
      &            'deallocation for cell_environ%nmbr_nghbr_cells failed')
    ENDIF
    DEALLOCATE( ptr_int%cell_environ%idx, STAT=ist )
    IF (ist /= SUCCESS)  THEN
      CALL finish ('mo_interpolation:destruct_int_state',&
      &            'deallocation for cell_environ%idx failed')
    ENDIF
    DEALLOCATE( ptr_int%cell_environ%blk, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:destruct_int_state',&
      &            'deallocation for cell_environ%blk failed')
    ENDIF
  END IF


  !
  ! pos_on_tplane_e
  !
  DEALLOCATE (ptr_int%pos_on_tplane_e, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for pos_on_tplane_e failed')
  ENDIF

  !
  ! pos_on_tplane_c_edge
  !
  DEALLOCATE (ptr_int%pos_on_tplane_c_edge, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for pos_on_tplane_c_edge failed')
  ENDIF

  !
  ! Least squares reconstruction
  !

  !
  ! *** linear ***
  !
  ! lsq_dim_stencil
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_dim_stencil, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                   &
      &          'deallocation for lsq_dim_stencil failed')
  ENDIF
  !
  ! lsq_idx_c
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_idx_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                   &
      &          'deallocation for lsq_idx_c failed')
  ENDIF
  !
  ! lsq_blk_c
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_blk_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                   &
      &          'deallocation for lsq_blk_c failed')
  ENDIF
  !
  ! lsq_weights_c
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_weights_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_weights_c failed')
  ENDIF
  !
  ! lsq_qtmat_c
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_qtmat_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_qtmat_c failed')
  ENDIF
  !
  ! lsq_rmat_rdiag_c
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_rmat_rdiag_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_rmat_rdiag_c failed')
  ENDIF
  !
  ! lsq_rmat_utri_c
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_rmat_utri_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_rmat_utri_c failed')
  ENDIF
  !
  ! lsq_pseudoinv
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_pseudoinv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_pseudoinv failed')
  ENDIF
  !
  ! lsq_moments
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_moments, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_moments failed')
  ENDIF
  !
  ! lsq_moments_hat
  !
  DEALLOCATE (ptr_int%lsq_lin%lsq_moments_hat, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_moments_hat failed')
  ENDIF

  !
  ! *** higher order ***
  !
  !
  ! lsq_dim_stencil
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_dim_stencil, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                   &
      &          'deallocation for lsq_dim_stencil failed')
  ENDIF
  !
  ! lsq_idx_c
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_idx_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                   &
      &          'deallocation for lsq_idx_c failed')
  ENDIF
  !
  ! lsq_blk_c
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_blk_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:construct_int_state',                   &
      &          'deallocation for lsq_blk_c failed')
  ENDIF
  !
  ! lsq_weights_c
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_weights_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_weights_c failed')
  ENDIF
  !
  ! lsq_qtmat_c
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_qtmat_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_qtmat_c failed')
  ENDIF
  !
  ! lsq_rmat_rdiag_c
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_rmat_rdiag_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_rmat_rdiag_c failed')
  ENDIF
  !
  ! lsq_rmat_utri_c
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_rmat_utri_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_rmat_utri_c failed')
  ENDIF
  !
  ! lsq_pseudoinv
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_pseudoinv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_pseudoinv failed')
  ENDIF
  !
  ! lsq_moments
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_moments, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_moments failed')
  ENDIF
  !
  ! lsq_moments_hat
  !
  DEALLOCATE (ptr_int%lsq_high%lsq_moments_hat, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                    &
      &          'deallocation for lsq_moments_hat failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_grdiv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_grdiv failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%qpts_tri_l, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for qpts_tri_l failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%qpts_tri_q, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for qpts_tri_q failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%qpts_tri_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for qpts_tri_q failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%weights_tri_q, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for weights_tri_q failed')
  ENDIF
  DEALLOCATE (ptr_int%gquad%weights_tri_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for weights_tri_c failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_div, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_div failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_rot, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_rot failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_n2s, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_n2s failed')
  ENDIF

  DEALLOCATE (ptr_int%geofac_grg, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                      &
      &             'deallocation for geofac_grg failed')
  ENDIF

  DEALLOCATE (ptr_int%primal_normal_ec, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                       &
      &             'deallocation for primal_normal_ec failed')
  ENDIF

  DEALLOCATE (ptr_int%edge_cell_length, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                       &
      &             'deallocation for edge_cell_length failed')
  ENDIF

  DEALLOCATE (ptr_int%cell_vert_dist, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:destruct_int_state',                       &
      &             'deallocation for cell_vert_dist failed')
  ENDIF

END SUBROUTINE deallocate_int_state


!-------------------------------------------------------------------------
!
!! Deallocation of components of 2d interpolation state.
!!
SUBROUTINE destruct_2d_interpol_state( ptr_int_state)
  !
  TYPE(t_int_state), INTENT(inout) :: ptr_int_state(n_dom_start:)
  ! local variables:
  CHARACTER(*), PARAMETER :: routine = TRIM("mo_interpolation:destruct_int_state")
  INTEGER                 :: jg

  !-----------------------------------------------------------------------

  CALL message(routine, 'start to destruct int state')

  DO jg = n_dom_start, n_dom
    CALL deallocate_int_state(ptr_int_state(jg))
  ENDDO

  CALL message (routine, 'destruction of interpolation state finished')

END SUBROUTINE destruct_2d_interpol_state


END MODULE mo_intp_state
