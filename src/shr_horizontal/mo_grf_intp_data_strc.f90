! Contains the interpolation routines needed for grid refinement.
!
! These had originally been included in mo_grf_interpolation but then were
! packed into a separate module to clean up the code
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

MODULE mo_grf_intp_data_strc

USE mo_kind, ONLY: wp

PUBLIC

!>
!----------------------------------------------------------------------------
!                                    Types
!----------------------------------------------------------------------------
TYPE t_gridref_single_state

  ! Arrays needed during the computation phase of the index relationships and the coefficients
  !
  INTEGER, ALLOCATABLE  :: grf_vec_ind_1a (:,:,:), & ! index arrays defining the stencil
                           grf_vec_ind_1b (:,:,:), & ! of surrounding edges for RBF
                           grf_vec_ind_2a (:,:,:), & ! interpolation to lateral boundaries of
                           grf_vec_ind_2b (:,:,:)    ! refined grids

  INTEGER, ALLOCATABLE  :: grf_vec_blk_1a (:,:,:), & ! index arrays defining the stencil
                           grf_vec_blk_1b (:,:,:), & ! of surrounding edges for RBF
                           grf_vec_blk_2a (:,:,:), & ! interpolation to lateral boundaries of
                           grf_vec_blk_2b (:,:,:)    ! refined grids

  INTEGER, ALLOCATABLE  :: grf_vec_stencil_1a (:,:), & ! arrays defining the number of surrounding
                           grf_vec_stencil_1b (:,:), & ! edges/cells in the stencil for RBF
                           grf_vec_stencil_2a (:,:), & ! interpolation to lateral boundaries of
                           grf_vec_stencil_2b (:,:)   ! refined grids

  REAL(wp), ALLOCATABLE :: grf_vec_coeff_1a(:,:,:), & ! arrays containing the coefficients
                           grf_vec_coeff_1b(:,:,:), & ! used for RBF vector interpolation
                           grf_vec_coeff_2a(:,:,:), & ! to lateral boundaries of refined grids
                           grf_vec_coeff_2b(:,:,:)

  ! distances between parent cells and child cells needed for gradient-based interpolation
  REAL(wp), ALLOCATABLE :: grf_dist_pc2cc(:,:,:,:), grf_dist_pe2ce(:,:,:)


  ! Arrays containing the index/block and coefficient lists allowing a flow control that does not
  ! rely on sorting the nest overlap points at parent grid level
  !
  INTEGER :: npoints_bdyintp_c, npoints_bdyintp_e, npoints_ubcintp_c, npoints_ubcintp_e, npoints_bdyintp_v

  INTEGER, ALLOCATABLE :: idxlist_bdyintp_c(:,:), idxlist_bdyintp_e(:,:),                  &
                          idxlist_ubcintp_c(:,:), idxlist_ubcintp_e(:,:),                  &
                          blklist_bdyintp_c(:,:), blklist_bdyintp_e(:,:),                  &
                          blklist_ubcintp_c(:,:), blklist_ubcintp_e(:,:),                  &
                          idxlist_rbfintp_v(:,:), blklist_rbfintp_v(:,:), edge_vert_idx(:,:)

  REAL(wp), ALLOCATABLE :: coeff_bdyintp_c(:,:,:), coeff_ubcintp_c(:,:,:),                   &
                           dist_pc2cc_bdy(:,:,:),  dist_pc2cc_ubc(:,:,:),  prim_norm(:,:,:), &
                           coeff_bdyintp_e12(:,:), coeff_bdyintp_e34(:,:), dist_pe2ce(:,:),  &
                           coeff_ubcintp_e12(:,:), coeff_ubcintp_e34(:,:), coeff_rbf_v(:,:,:)

END TYPE t_gridref_single_state

TYPE t_gridref_state

  ! p_dom holds fields that are allocated for the boundary interpolation zone only
  ! They therefore need to be defined separately for each child domain
  TYPE(t_gridref_single_state), ALLOCATABLE :: p_dom(:)

  ! These fields are allocated for the full parent domain and thus do not need
  ! to be held separately for each child domain
  REAL(wp), ALLOCATABLE :: fbk_wgt_aw(:,:,:)    ! Area-weighted feedback weights for cell-based variables
                                                !  dim1=1,nproma, dim2=nblks_c, dim3=4
  REAL(wp), ALLOCATABLE :: fbk_wgt_bln(:,:,:)   ! Bilinear feedback weights for cell-based variables
                                                !  dim1=1,nproma, dim2=nblks_c, dim3=4
  REAL(wp), ALLOCATABLE :: fbk_wgt_e(:,:,:)     ! Feedback weights for edge-based variables

  ! Mask fields defined at parent level for feedback overlap regions with a nested domain:
  ! dim1=nproma, dim2=nblks, dim3=n_childdom 
  LOGICAL, ALLOCATABLE :: mask_ovlp_c(:,:,:), mask_ovlp_ch(:,:,:), mask_ovlp_e(:,:,:), mask_ovlp_v(:,:,:)

  ! Number of grid points for which tendencies for lateral boundary interpolation need to be computed
  ! at the parent grid level
  INTEGER :: npoints_bdyintp_src_c, npoints_bdyintp_src_e

  ! Index lists of grid points for which tendencies for lateral boundary interpolation need to be computed
  ! at the parent grid level
  INTEGER, ALLOCATABLE :: idxlist_bdyintp_src_c(:), idxlist_bdyintp_src_e(:), &
                          blklist_bdyintp_src_c(:), blklist_bdyintp_src_e(:)

END TYPE t_gridref_state

  !----------------------------------------------------------------------------
  !                                Variables
  !----------------------------------------------------------------------------
  TYPE(t_gridref_state), TARGET, ALLOCATABLE :: p_grf_state(:), p_grf_state_local_parent(:)


END MODULE mo_grf_intp_data_strc
