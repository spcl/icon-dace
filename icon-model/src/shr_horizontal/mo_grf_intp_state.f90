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

! Contains the interpolation routines needed for grid refinement.
!
! These had originally been included in mo_grf_interpolation but then were
! packed into a separate module to clean up the code

! #ifdef __xlC__
! @PROCESS HOT
! #endif
#ifdef __PGI
!pgi$g opt=1
#endif

MODULE mo_grf_intp_state
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
USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, min_rledge_int, min_rlvert_int, &
                                  min_rlcell, min_rledge, min_rlvert
USE mo_model_domain,        ONLY: t_patch
USE mo_grid_config,         ONLY: n_dom, n_dom_start
USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e, grf_nudgintp_start_c, &
                                  grf_nudgintp_start_e, grf_bdyintp_end_c, grf_bdyintp_end_e,     &
                                  grf_fbk_start_c, grf_fbk_start_e, grf_bdywidth_c,               &
                                  grf_bdywidth_e
USE mo_parallel_config,     ONLY: nproma
USE mo_communication,       ONLY: t_p_comm_pattern, blk_no, idx_no, idx_1d, &
  &                               delete_comm_pattern, &
  &                               exchange_data, t_comm_pattern_collection, &
  &                               t_comm_pattern, &
  &                               delete_comm_pattern_collection
USE mo_communication_factory, ONLY: setup_comm_pattern, &
  &                                 setup_comm_pattern_collection
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup, get_valid_local_index, &
                                  init_glb2loc_index_lookup, set_inner_glb_index, &
                                  deallocate_glb2loc_index_lookup

USE mo_grf_intp_data_strc,  ONLY: t_gridref_single_state, t_gridref_state
USE mo_grf_intp_coeffs,     ONLY: compute_pc2cc_distances, compute_pe2ce_distances, &
  &                               init_fbk_wgt, grf_index, rbf_compute_coeff_grf_e
USE mo_dist_dir,            ONLY: dist_dir_get_owners, t_dist_dir



IMPLICIT NONE

PRIVATE

PUBLIC :: construct_2d_gridref_state, destruct_2d_gridref_state
PUBLIC :: allocate_grf_state, deallocate_grf_state
PUBLIC :: transfer_grf_state
PUBLIC :: create_grf_index_lists, destruct_interpol_patterns

CLASS(t_comm_pattern), POINTER :: comm_pat_loc_to_glb_c, comm_pat_loc_to_glb_e

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_grf_intp_state'

CONTAINS

!-------------------------------------------------------------------------
!
!! Allocation of fields needed for grid refinement.
!!
SUBROUTINE allocate_grf_state( ptr_patch, ptr_grf )
!
TYPE(t_patch), INTENT(IN) :: ptr_patch
TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: ptr_grf

CHARACTER(LEN=*), PARAMETER :: routine = modname//':allocate_grf_state'
TYPE(t_gridref_single_state), POINTER :: ptr_int    => NULL()

INTEGER :: jcd
INTEGER :: nblks_c, nblks_e, &
           nproma_grf, isb_e, ieb_e, isb_c, ieb_c, i_nchdom
INTEGER :: ist
INTEGER :: grf_vec_dim_1, grf_vec_dim_2

!-----------------------------------------------------------------------

  grf_vec_dim_1 = 6
  grf_vec_dim_2 = 5

  ! determine size of arrays, i.e.
  ! values for the blocking
  !
  nblks_c  = ptr_patch%nblks_c
  nblks_e  = ptr_patch%nblks_e

  ! determine number of child domains
  i_nchdom = MAX(1,ptr_patch%n_childdom)

  ! build data structure for grid refinement for each child domain
  ALLOCATE(ptr_grf%p_dom(i_nchdom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine,  &
      &             'allocation for p_dom failed')
  ENDIF

  ! Now loop over all child domains of the present domain
  CD_LOOP: DO jcd = 1, i_nchdom

    ptr_int => ptr_grf%p_dom(jcd)

    ! Determine dimensions needed for grid refinement fields
    IF (ptr_patch%n_childdom > 0) THEN  ! n_childdom = 0 in the innermost nest(s)

      isb_e      = 1
      ieb_e      = nblks_e
      isb_c      = 1
      ieb_c      = nblks_c
      nproma_grf = nproma

    ELSE

      isb_e      = 1
      ieb_e      = 1
      isb_c      = 1
      ieb_c      = 1
      nproma_grf = 1

    ENDIF

    ! index arrays (external points; global indices) for grid refinement interpolation

    ALLOCATE (ptr_int%grf_vec_ind_1a(nproma_grf, grf_vec_dim_1, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_ind_1a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_ind_1b(nproma_grf, grf_vec_dim_1, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_ind_1b failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_ind_2a(nproma_grf, grf_vec_dim_2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_ind_2a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_ind_2b(nproma_grf, grf_vec_dim_2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_ind_2b failed')
    ENDIF

    ! index arrays (external points; global indices) for grid refinement interpolation

    ALLOCATE (ptr_int%grf_vec_blk_1a(nproma_grf, grf_vec_dim_1, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_blk_1a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_blk_1b(nproma_grf, grf_vec_dim_1, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_blk_1b failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_blk_2a(nproma_grf, grf_vec_dim_2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_blk_2a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_blk_2b(nproma_grf, grf_vec_dim_2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_blk_2b failed')
    ENDIF

    ! arrays carrying the actual number of edges available for reconstruction

    ALLOCATE (ptr_int%grf_vec_stencil_1a(nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_stencil_1a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_stencil_1b(nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_stencil_1b failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_stencil_2a(nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_stencil_2a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_stencil_2b(nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_stencil_2b failed')
    ENDIF

    ! arrays carrying the interpolation coefficients for reconstruction related to grid refinement

    ALLOCATE (ptr_int%grf_vec_coeff_1a(grf_vec_dim_1, nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_coeff_1a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_coeff_1b(grf_vec_dim_1, nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_coeff_1b failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_coeff_2a(grf_vec_dim_2, nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_coeff_2a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_coeff_2b(grf_vec_dim_2, nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_vec_coeff_2b failed')
    ENDIF

    ! Distances from parent cell to child cells
    ALLOCATE (ptr_int%grf_dist_pc2cc(nproma_grf, 4, 2, isb_c:ieb_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_dist_pc2cc failed')
    ENDIF

    ! Normalized distances from parent edge to child edges 1 and 2
    ALLOCATE (ptr_int%grf_dist_pe2ce(nproma_grf, 2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'allocation for grf_dist_pe2ce failed')
    ENDIF

    ptr_int%grf_vec_ind_1a   = 0
    ptr_int%grf_vec_ind_1b   = 0
    ptr_int%grf_vec_ind_2a   = 0
    ptr_int%grf_vec_ind_2b   = 0

    ptr_int%grf_vec_blk_1a   = 0
    ptr_int%grf_vec_blk_1b   = 0
    ptr_int%grf_vec_blk_2a   = 0
    ptr_int%grf_vec_blk_2b   = 0

    ptr_int%grf_vec_stencil_1a   = 0
    ptr_int%grf_vec_stencil_1b   = 0
    ptr_int%grf_vec_stencil_2a   = 0
    ptr_int%grf_vec_stencil_2b   = 0

    ptr_int%grf_vec_coeff_1a   = 0._wp
    ptr_int%grf_vec_coeff_1b   = 0._wp
    ptr_int%grf_vec_coeff_2a   = 0._wp
    ptr_int%grf_vec_coeff_2b   = 0._wp

    ptr_int%grf_dist_pc2cc   = 0._wp
    ptr_int%grf_dist_pe2ce   = 0._wp

  ENDDO CD_LOOP

  ! Feedback weights for cell-based variables
  ALLOCATE (ptr_grf%fbk_wgt_aw(nproma,nblks_c,4), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine,                       &
      &             'allocation for fbk_wgt_aw failed')
  ENDIF
  ALLOCATE (ptr_grf%fbk_wgt_bln(nproma,nblks_c,4), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine,                       &
      &             'allocation for fbk_wgt_bln failed')
  ENDIF

  ! Feedback weights for edge-based variables
  ALLOCATE (ptr_grf%fbk_wgt_e(nproma,nblks_e,6), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine,                       &
      &             'allocation for fbk_wgt_e failed')
  ENDIF

  ptr_grf%fbk_wgt_aw    = 0._wp
  ptr_grf%fbk_wgt_bln   = 0._wp
  ptr_grf%fbk_wgt_e     = 0._wp

END SUBROUTINE allocate_grf_state

!-------------------------------------------------------------------------
!
!! Allocation of fields needed for grid refinement.
!!
!! Initialization of components.
!!
SUBROUTINE construct_2d_gridref_state( p_patch, p_patch_local_parent, &
  &                                    p_grf_state, p_grf_state_local_parent )
  !
  TYPE(t_patch),         INTENT(IN)    :: p_patch                  (n_dom_start:)
  TYPE(t_patch),         INTENT(IN)    :: p_patch_local_parent     (n_dom_start+1:)
  TYPE(t_gridref_state), INTENT(INOUT) :: p_grf_state              (n_dom_start:)
  TYPE(t_gridref_state), INTENT(INOUT) :: p_grf_state_local_parent (n_dom_start+1:)

  CHARACTER(LEN=*), PARAMETER :: routine = modname//':construct_2d_gridref_state'
  INTEGER :: jg

!-----------------------------------------------------------------------

  CALL message(routine, 'start to construct grf_state')

  DO jg = n_dom_start, n_dom
    CALL allocate_grf_state(p_patch(jg), p_grf_state(jg))
  ENDDO
  DO jg = n_dom_start+1, n_dom
    CALL allocate_grf_state(p_patch_local_parent(jg), p_grf_state_local_parent(jg))
  ENDDO

  CALL message (routine, 'memory allocation finished')


  CALL init_fbk_wgt (p_patch(n_dom_start:), p_patch_local_parent(n_dom_start+1:), &
    &                p_grf_state_local_parent(n_dom_start+1:))

  CALL compute_pc2cc_distances (p_patch(n_dom_start:), p_patch_local_parent(n_dom_start+1:), &
    &                           p_grf_state_local_parent(n_dom_start+1:))
  CALL compute_pe2ce_distances (p_patch(n_dom_start:), p_patch_local_parent(n_dom_start+1:), &
    &                           p_grf_state_local_parent(n_dom_start+1:))

  CALL grf_index (p_patch(n_dom_start:), p_patch_local_parent(n_dom_start+1:), &
    &             p_grf_state_local_parent(n_dom_start+1:))

  CALL rbf_compute_coeff_grf_e (p_patch(n_dom_start:), p_patch_local_parent(n_dom_start+1:), &
    &                           p_grf_state_local_parent(n_dom_start+1:))


  CALL message (routine, 'construction of interpolation state finished')

END SUBROUTINE construct_2d_gridref_state

!-------------------------------------------------------------------------
!!
!>
!! Get local idx_no and blk_no
!!
ELEMENTAL SUBROUTINE loc_idx_blk_no(glb2loc_index, glb_index, loc_idx, loc_blk)

  TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index
  INTEGER, INTENT(IN) :: glb_index
  INTEGER, INTENT(OUT) :: loc_idx, loc_blk

  INTEGER :: valid_local_index

  valid_local_index = get_valid_local_index(glb2loc_index, glb_index, .TRUE.)

  loc_idx = idx_no(valid_local_index)
  loc_blk = blk_no(valid_local_index)
END SUBROUTINE loc_idx_blk_no

!-------------------------------------------------------------------------
!!
!>
!! Get global idx_1d for edges
!!
ELEMENTAL INTEGER FUNCTION glb_idx_1d_e(p_p, idx, blk)
  TYPE(t_patch), INTENT(IN) :: p_p
  INTEGER, INTENT(IN) :: idx, blk

  IF(idx<=0 .or. blk<=0 .or. idx_1d(idx, blk)>p_p%n_patch_edges) THEN
    glb_idx_1d_e = 0
  ELSE
    glb_idx_1d_e = p_p%edges%decomp_info%glb_index(idx_1d(idx, blk))
  ENDIF
END FUNCTION glb_idx_1d_e
!-------------------------------------------------------------------------
!!
!! Transfers interpolation state from local parent to global parent
!! Some variables have to be transferred the other way round - see below
!!
SUBROUTINE transfer_grf_state(p_p, p_lp, p_grf, p_lgrf, jcd)
!
  TYPE(t_patch), INTENT(IN) :: p_p
  TYPE(t_patch), INTENT(IN) :: p_lp

  TYPE(t_gridref_state), INTENT(INOUT) :: p_grf
  TYPE(t_gridref_state), INTENT(INOUT) :: p_lgrf

  INTEGER, INTENT(IN) :: jcd

  INTEGER, ALLOCATABLE :: owner(:)
  LOGICAL, ALLOCATABLE :: mask(:)
  INTEGER :: j, k, l, m, n, i_nchdom, jc, je, jb, icid
  INTEGER :: isb_e, ieb_e, isb_le, ieb_le, is_e, ie_e
  INTEGER :: isb_c, ieb_c, isb_lc, ieb_lc, is_c, ie_c

  REAL(wp), ALLOCATABLE :: z_tmp_s(:,:,:), z_tmp_r(:,:,:)

  INTEGER :: grf_vec_dim_1, grf_vec_dim_2

  grf_vec_dim_1 = 6
  grf_vec_dim_2 = 5

  i_nchdom = MAX(1,p_p%n_childdom)
  icid     = p_p%child_id(jcd)

  isb_e      = p_p%edges%start_blk(1,1)
  ieb_e      = p_p%edges%end_blk(min_rledge_int-1,i_nchdom)
  isb_c      = p_p%cells%start_blk(1,1)
  ieb_c      = p_p%cells%end_blk(min_rlcell_int,i_nchdom)

  isb_le      = p_lp%edges%start_blk(grf_bdyintp_start_e,jcd)
  ieb_le      = p_lp%edges%end_blk(min_rledge_int,jcd)
  isb_lc      = p_lp%cells%start_blk(grf_bdyintp_start_c,jcd)
  ieb_lc      = p_lp%cells%end_blk(min_rlcell_int,jcd)

  is_e = idx_1d(p_p%edges%start_idx(1,1), &
                p_p%edges%start_blk(1,1))
  ie_e = idx_1d(p_p%edges%end_idx(min_rledge_int-1,i_nchdom), &
                p_p%edges%end_blk(min_rledge_int-1,i_nchdom))

  is_c = idx_1d(p_p%cells%start_idx(1,1), &
                p_p%cells%start_blk(1,1))
  ie_c = idx_1d(p_p%cells%end_idx(min_rlcell_int,i_nchdom), &
                p_p%cells%end_blk(min_rlcell_int,i_nchdom))

  ! Set up communication patterns for transferring the data to local parents.
  ! Since these communication patterns are not used elsewhere, they are
  ! stored locally and deleted at the end of the routine


  ALLOCATE(owner(MAX(p_p%n_patch_edges, p_p%n_patch_cells)), &
    &      mask(MAX(p_p%n_patch_edges, p_p%n_patch_cells)))

  mask(1:p_p%n_patch_edges) = .FALSE.
  DO j = is_e, ie_e
    jb = blk_no(j)
    je = idx_no(j)
    mask(j) = p_p%edges%refin_ctrl(je,jb) <= grf_bdyintp_start_e .AND. &
      &       p_p%edges%child_id(je,jb) == icid
  ENDDO
  owner(1:p_p%n_patch_edges) = &
       dist_dir_get_owners(p_lp%edges%decomp_info%owner_dist_dir, &
       &                   p_p%edges%decomp_info%glb_index(:), &
       &                   mask(1:p_p%n_patch_edges))
  CALL setup_comm_pattern(p_p%n_patch_edges, owner(1:p_p%n_patch_edges), &
    &                     p_p%edges%decomp_info%glb_index, &
    &                     p_lp%edges%decomp_info%glb2loc_index, &
    &                     p_lp%n_patch_edges, &
    &                     p_lp%edges%decomp_info%owner_local, &
    &                     p_lp%edges%decomp_info%glb_index, &
    &                     comm_pat_loc_to_glb_e)

  mask(1:p_p%n_patch_cells) = .FALSE.
  DO j = is_c, ie_c
    jb = blk_no(j)
    jc = idx_no(j)
    mask(j) = p_p%cells%refin_ctrl(jc,jb) <= grf_bdyintp_start_c .AND. &
      &       p_p%cells%child_id(jc,jb) == icid
  ENDDO
  owner(1:p_p%n_patch_cells) = &
       dist_dir_get_owners(p_lp%cells%decomp_info%owner_dist_dir, &
       &                   p_p%cells%decomp_info%glb_index(:), &
       &                   mask(1:p_p%n_patch_cells))
  CALL setup_comm_pattern(p_p%n_patch_cells, owner(1:p_p%n_patch_cells), &
    &                     p_p%cells%decomp_info%glb_index, &
    &                     p_lp%cells%decomp_info%glb2loc_index, &
    &                     p_lp%n_patch_cells, &
    &                     p_lp%cells%decomp_info%owner_local, &
    &                     p_lp%cells%decomp_info%glb_index, &
    &                     comm_pat_loc_to_glb_c)
  DEALLOCATE(owner, mask)

  ALLOCATE(z_tmp_s(nproma,4,p_lp%nblks_e))
  ALLOCATE(z_tmp_r(nproma,4,p_p%nblks_e))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp
  z_tmp_s(:,1,isb_le:ieb_le) = REAL(p_lgrf%p_dom(jcd)%grf_vec_stencil_1a(:,isb_le:ieb_le),wp)
  z_tmp_s(:,2,isb_le:ieb_le) = REAL(p_lgrf%p_dom(jcd)%grf_vec_stencil_1b(:,isb_le:ieb_le),wp)
  z_tmp_s(:,3,isb_le:ieb_le) = REAL(p_lgrf%p_dom(jcd)%grf_vec_stencil_2a(:,isb_le:ieb_le),wp)
  z_tmp_s(:,4,isb_le:ieb_le) = REAL(p_lgrf%p_dom(jcd)%grf_vec_stencil_2b(:,isb_le:ieb_le),wp)

  CALL exchange_data(p_pat=comm_pat_loc_to_glb_e, lacc=.FALSE., RECV=z_tmp_r, SEND=z_tmp_s)

  IF (ieb_e >= isb_e) THEN
    p_grf%p_dom(jcd)%grf_vec_stencil_1a(:,isb_e:ieb_e) = INT(z_tmp_r(:,1,isb_e:ieb_e))
    p_grf%p_dom(jcd)%grf_vec_stencil_1b(:,isb_e:ieb_e) = INT(z_tmp_r(:,2,isb_e:ieb_e))
    p_grf%p_dom(jcd)%grf_vec_stencil_2a(:,isb_e:ieb_e) = INT(z_tmp_r(:,3,isb_e:ieb_e))
    p_grf%p_dom(jcd)%grf_vec_stencil_2b(:,isb_e:ieb_e) = INT(z_tmp_r(:,4,isb_e:ieb_e))
  ENDIF
  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

  ALLOCATE(z_tmp_s(nproma,2*grf_vec_dim_1+2*grf_vec_dim_2,p_lp%nblks_e))
  ALLOCATE(z_tmp_r(nproma,2*grf_vec_dim_1+2*grf_vec_dim_2,p_p%nblks_e))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  n = 0
  DO k = 1, grf_vec_dim_1
    n = n+1
!CDIR IEXPAND
    z_tmp_s(:,n,isb_le:ieb_le) =REAL(glb_idx_1d_e(p_lp,p_lgrf%p_dom(jcd)%grf_vec_ind_1a(:,k,isb_le:ieb_le), &
                                                       p_lgrf%p_dom(jcd)%grf_vec_blk_1a(:,k,isb_le:ieb_le)),wp)
  ENDDO
  DO k = 1, grf_vec_dim_1
    n = n+1
!CDIR IEXPAND
    z_tmp_s(:,n,isb_le:ieb_le) =REAL(glb_idx_1d_e(p_lp,p_lgrf%p_dom(jcd)%grf_vec_ind_1b(:,k,isb_le:ieb_le), &
                                                       p_lgrf%p_dom(jcd)%grf_vec_blk_1b(:,k,isb_le:ieb_le)),wp)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
!CDIR IEXPAND
    z_tmp_s(:,n,isb_le:ieb_le) =REAL(glb_idx_1d_e(p_lp,p_lgrf%p_dom(jcd)%grf_vec_ind_2a(:,k,isb_le:ieb_le), &
                                                       p_lgrf%p_dom(jcd)%grf_vec_blk_2a(:,k,isb_le:ieb_le)),wp)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
!CDIR IEXPAND
    z_tmp_s(:,n,isb_le:ieb_le) =REAL(glb_idx_1d_e(p_lp,p_lgrf%p_dom(jcd)%grf_vec_ind_2b(:,k,isb_le:ieb_le), &
                                                       p_lgrf%p_dom(jcd)%grf_vec_blk_2b(:,k,isb_le:ieb_le)),wp)
  ENDDO

  CALL exchange_data(p_pat=comm_pat_loc_to_glb_e, lacc=.FALSE., RECV=z_tmp_r, SEND=z_tmp_s)

  IF (ieb_e >= isb_e) THEN
    n = 0
    DO k = 1, grf_vec_dim_1
      n = n+1
      DO l = isb_e, ieb_e
        DO m = 1, nproma
          CALL loc_idx_blk_no(p_p%edges%decomp_info%glb2loc_index, &
            &                 int(z_tmp_r(m,n,l)), &
            &                 p_grf%p_dom(jcd)%grf_vec_ind_1a(m,k,l), &
            &                 p_grf%p_dom(jcd)%grf_vec_blk_1a(m,k,l))
        ENDDO
      ENDDO
    ENDDO
    DO k = 1, grf_vec_dim_1
      n = n+1
      DO l = isb_e, ieb_e
        DO m = 1, nproma
          CALL loc_idx_blk_no(p_p%edges%decomp_info%glb2loc_index, &
            &                 int(z_tmp_r(m,n,l)), &
            &                 p_grf%p_dom(jcd)%grf_vec_ind_1b(m,k,l), &
            &                 p_grf%p_dom(jcd)%grf_vec_blk_1b(m,k,l))
        ENDDO
      ENDDO
    ENDDO
    DO k = 1, grf_vec_dim_2
      n = n+1
      DO l = isb_e, ieb_e
        DO m = 1, nproma
          CALL loc_idx_blk_no(p_p%edges%decomp_info%glb2loc_index, &
            &                 int(z_tmp_r(m,n,l)), &
            &                 p_grf%p_dom(jcd)%grf_vec_ind_2a(m,k,l), &
            &                 p_grf%p_dom(jcd)%grf_vec_blk_2a(m,k,l))
        ENDDO
      ENDDO
    ENDDO
    DO k = 1, grf_vec_dim_2
      n = n+1
      DO l = isb_e, ieb_e
        DO m = 1, nproma
          CALL loc_idx_blk_no(p_p%edges%decomp_info%glb2loc_index, &
            &                 int(z_tmp_r(m,n,l)), &
            &                 p_grf%p_dom(jcd)%grf_vec_ind_2b(m,k,l), &
            &                 p_grf%p_dom(jcd)%grf_vec_blk_2b(m,k,l))
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  n = 0
  DO k = 1, grf_vec_dim_1
    n = n+1
    z_tmp_s(:,n,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_vec_coeff_1a(k,:,isb_le:ieb_le)
  ENDDO
  DO k = 1, grf_vec_dim_1
    n = n+1
    z_tmp_s(:,n,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_vec_coeff_1b(k,:,isb_le:ieb_le)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
    z_tmp_s(:,n,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_vec_coeff_2a(k,:,isb_le:ieb_le)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
    z_tmp_s(:,n,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_vec_coeff_2b(k,:,isb_le:ieb_le)
  ENDDO

  CALL exchange_data(p_pat=comm_pat_loc_to_glb_e, lacc=.FALSE., RECV=z_tmp_r, SEND=z_tmp_s)

  IF (ieb_e >= isb_e) THEN
    n = 0
    DO k = 1, grf_vec_dim_1
      n = n+1
      p_grf%p_dom(jcd)%grf_vec_coeff_1a(k,:,isb_e:ieb_e) = z_tmp_r(:,n,isb_e:ieb_e)
    ENDDO
    DO k = 1, grf_vec_dim_1
      n = n+1
      p_grf%p_dom(jcd)%grf_vec_coeff_1b(k,:,isb_e:ieb_e) = z_tmp_r(:,n,isb_e:ieb_e)
    ENDDO
    DO k = 1, grf_vec_dim_2
      n = n+1
      p_grf%p_dom(jcd)%grf_vec_coeff_2a(k,:,isb_e:ieb_e) = z_tmp_r(:,n,isb_e:ieb_e)
    ENDDO
    DO k = 1, grf_vec_dim_2
      n = n+1
      p_grf%p_dom(jcd)%grf_vec_coeff_2b(k,:,isb_e:ieb_e) = z_tmp_r(:,n,isb_e:ieb_e)
    ENDDO
  ENDIF

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

  ALLOCATE(z_tmp_s(nproma,2,p_lp%nblks_e))
  ALLOCATE(z_tmp_r(nproma,2,p_p%nblks_e))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  z_tmp_s(:,:,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_dist_pe2ce(:,:,isb_le:ieb_le)

  CALL exchange_data(p_pat=comm_pat_loc_to_glb_e, lacc=.FALSE., RECV=z_tmp_r, SEND=z_tmp_s)

  IF (ieb_e >= isb_e) p_grf%p_dom(jcd)%grf_dist_pe2ce(:,:,isb_e:ieb_e) = z_tmp_r(:,:,isb_e:ieb_e)

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

  ALLOCATE(z_tmp_s(nproma,8,p_lp%nblks_c))
  ALLOCATE(z_tmp_r(nproma,8,p_p%nblks_c))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  z_tmp_s(:,1,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,1,1,isb_lc:ieb_lc)
  z_tmp_s(:,2,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,2,1,isb_lc:ieb_lc)
  z_tmp_s(:,3,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,3,1,isb_lc:ieb_lc)
  z_tmp_s(:,4,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,4,1,isb_lc:ieb_lc)
  z_tmp_s(:,5,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,1,2,isb_lc:ieb_lc)
  z_tmp_s(:,6,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,2,2,isb_lc:ieb_lc)
  z_tmp_s(:,7,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,3,2,isb_lc:ieb_lc)
  z_tmp_s(:,8,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,4,2,isb_lc:ieb_lc)

  CALL exchange_data(p_pat=comm_pat_loc_to_glb_c, lacc=.FALSE., RECV=z_tmp_r, SEND=z_tmp_s)

  IF (ieb_c >= isb_c) THEN
    p_grf%p_dom(jcd)%grf_dist_pc2cc(:,1,1,isb_c:ieb_c) = z_tmp_r(:,1,isb_c:ieb_c)
    p_grf%p_dom(jcd)%grf_dist_pc2cc(:,2,1,isb_c:ieb_c) = z_tmp_r(:,2,isb_c:ieb_c)
    p_grf%p_dom(jcd)%grf_dist_pc2cc(:,3,1,isb_c:ieb_c) = z_tmp_r(:,3,isb_c:ieb_c)
    p_grf%p_dom(jcd)%grf_dist_pc2cc(:,4,1,isb_c:ieb_c) = z_tmp_r(:,4,isb_c:ieb_c)
    p_grf%p_dom(jcd)%grf_dist_pc2cc(:,1,2,isb_c:ieb_c) = z_tmp_r(:,5,isb_c:ieb_c)
    p_grf%p_dom(jcd)%grf_dist_pc2cc(:,2,2,isb_c:ieb_c) = z_tmp_r(:,6,isb_c:ieb_c)
    p_grf%p_dom(jcd)%grf_dist_pc2cc(:,3,2,isb_c:ieb_c) = z_tmp_r(:,7,isb_c:ieb_c)
    p_grf%p_dom(jcd)%grf_dist_pc2cc(:,4,2,isb_c:ieb_c) = z_tmp_r(:,8,isb_c:ieb_c)
  ENDIF

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_loc_to_glb_c)
!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_loc_to_glb_e)

  ! Now create communication patterns for the complete patch

  ALLOCATE(owner(MAX(p_p%n_patch_edges, p_p%n_patch_cells)))
  owner(1:p_p%n_patch_edges) = &
    dist_dir_get_owners(p_lp%edges%decomp_info%owner_dist_dir, &
      &                 p_p%edges%decomp_info%glb_index(:))
  CALL setup_comm_pattern(p_p%n_patch_edges, owner(1:p_p%n_patch_edges), &
    &                     p_p%edges%decomp_info%glb_index, &
    &                     p_lp%edges%decomp_info%glb2loc_index, &
    &                     p_lp%n_patch_edges, &
    &                     p_lp%edges%decomp_info%owner_local, &
    &                     p_lp%edges%decomp_info%glb_index, &
    &                     comm_pat_loc_to_glb_e)

  owner(1:p_p%n_patch_cells) = &
    dist_dir_get_owners(p_lp%cells%decomp_info%owner_dist_dir, &
      &                 p_p%cells%decomp_info%glb_index(:))
  CALL setup_comm_pattern(p_p%n_patch_cells, owner(1:p_p%n_patch_cells), &
    &                     p_p%cells%decomp_info%glb_index, &
    &                     p_lp%cells%decomp_info%glb2loc_index, &
    &                     p_lp%n_patch_cells, &
    &                     p_lp%cells%decomp_info%owner_local, &
    &                     p_lp%cells%decomp_info%glb_index, &
    &                     comm_pat_loc_to_glb_c)
  DEALLOCATE(owner)

  ALLOCATE(z_tmp_s(nproma,8,p_lp%nblks_c))
  ALLOCATE(z_tmp_r(nproma,8,p_p%nblks_c))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  n = 0
  DO k = 1, 4
    n = n+1
    z_tmp_s(:,n,:) = p_lgrf%fbk_wgt_aw(:,:,k)
  ENDDO
  DO k = 1, 4
    n = n+1
    z_tmp_s(:,n,:) = p_lgrf%fbk_wgt_bln(:,:,k)
  ENDDO

  CALL exchange_data(p_pat=comm_pat_loc_to_glb_c, lacc=.FALSE., RECV=z_tmp_r, SEND=z_tmp_s)

  n = 0
  DO k = 1, 4
    n = n+1
    p_grf%fbk_wgt_aw(:,:,k) = z_tmp_r(:,n,:)
  ENDDO
  DO k = 1, 4
    n = n+1
    p_grf%fbk_wgt_bln(:,:,k) = z_tmp_r(:,n,:)
  ENDDO

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

  ALLOCATE(z_tmp_s(nproma,8,p_lp%nblks_e))
  ALLOCATE(z_tmp_r(nproma,8,p_p%nblks_e))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  n = 0
  DO k = 1, 6
    n = n+1
    z_tmp_s(:,n,:) = p_lgrf%fbk_wgt_e(:,:,k)
  ENDDO

  CALL exchange_data(p_pat=comm_pat_loc_to_glb_e, lacc=.FALSE., RECV=z_tmp_r, SEND=z_tmp_s)

  n = 0
  DO k = 1, 6
    n = n+1
    p_grf%fbk_wgt_e(:,:,k) = z_tmp_r(:,n,:)
  ENDDO

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_loc_to_glb_c)
!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_loc_to_glb_e)

END SUBROUTINE transfer_grf_state


!-------------------------------------------------------------------------

FUNCTION scal_grf_select_func(idx, blk, n, p_patch) RESULT(p)
  INTEGER, INTENT(IN) :: idx, blk, n
  TYPE(t_patch), INTENT(IN) :: p_patch
  LOGICAL :: p

  p = p_patch%cells%refin_ctrl(idx,blk) > 0 .AND. &
    & p_patch%cells%refin_ctrl(idx,blk) <= grf_bdywidth_c .AND. &
    & p_patch%cells%pc_idx(idx,blk) == n

END FUNCTION scal_grf_select_func

FUNCTION vec_grf_select_func(idx, blk, n, p_patch) RESULT(p)
  INTEGER, INTENT(IN) :: idx, blk, n
  TYPE(t_patch), INTENT(IN) :: p_patch
  LOGICAL :: p

  p = p_patch%edges%refin_ctrl(idx,blk) > 0 .AND. &
    & p_patch%edges%refin_ctrl(idx,blk) <= grf_bdywidth_e .AND. &
    & p_patch%edges%pc_idx(idx,blk) == n

END FUNCTION vec_grf_select_func

FUNCTION scal_ubc_select_func(idx, blk, n, p_patch) RESULT(p)
  INTEGER, INTENT(IN) :: idx, blk, n
  TYPE(t_patch), INTENT(IN) :: p_patch
  LOGICAL :: p

  p = (p_patch%cells%refin_ctrl(idx,blk) >= grf_bdywidth_c - 1 .OR. &
    &  p_patch%cells%refin_ctrl(idx,blk) <= 0) .AND. &
    &  p_patch%cells%pc_idx(idx,blk) == n

END FUNCTION scal_ubc_select_func

FUNCTION vec_ubc_select_func(idx, blk, n, p_patch) RESULT(p)
  INTEGER, INTENT(IN) :: idx, blk, n
  TYPE(t_patch), INTENT(IN) :: p_patch
  LOGICAL :: p

  p = (p_patch%edges%refin_ctrl(idx,blk) >= grf_bdywidth_e + 1 .OR. &
    &  p_patch%edges%refin_ctrl(idx,blk) <= 0) .AND. &
    &  p_patch%edges%pc_idx(idx,blk) == n

END FUNCTION vec_ubc_select_func

!-------------------------------------------------------------------------
!
!! Creates index lists for lateral/upper boundary interpolation needed to
!! get rid of the grid point reordering in the parent domain
!!
SUBROUTINE create_grf_index_lists( p_patch_all, p_grf_state, p_int_state )
  !
  TYPE(t_patch), TARGET, INTENT(INOUT)         :: p_patch_all(n_dom_start:)
  TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: p_grf_state(n_dom_start:)
  TYPE(t_int_state), TARGET, INTENT(IN)        :: p_int_state(n_dom_start:)

  CHARACTER(LEN=*), PARAMETER :: routine = modname//':create_grf_index_lists'
  TYPE(t_gridref_single_state), POINTER :: p_grf_s
  TYPE(t_gridref_state),        POINTER :: p_grf
  TYPE(t_patch),                POINTER :: p_patch
  TYPE(t_int_state),            POINTER :: p_int

  INTEGER :: jcd, icid, i_nchdom, ist, n
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: jg, jb, jc, je, jv, i, ic, ib, iv, iv1, iv2, ib1, ib2
  INTEGER :: npoints_lbc, npoints_ubc, icount_lbc, icount_ubc, npoints_lbc_src, icount_lbc_src
  LOGICAL :: lprocess, lfound(2,6)
  INTEGER, ALLOCATABLE :: inv_ind_v(:,:)
  TYPE(t_glb2loc_index_lookup) :: inv_glb2loc
  INTEGER, ALLOCATABLE :: inv_glb2loc_loc_index_c_grf(:), &
    &                     inv_glb2loc_loc_index_c_ubc(:), &
    &                     inv_glb2loc_loc_index_e_lbc(:), &
    &                     inv_glb2loc_loc_index_e_ubc(:)
  INTEGER, ALLOCATABLE :: owner_local(:)
  INTEGER, ALLOCATABLE :: glb_index(:)

!-----------------------------------------------------------------------

  DO jg = 1, n_dom

   p_patch => p_patch_all(jg)
   p_int   => p_int_state(jg)
   p_grf   => p_grf_state(jg)

   ! number of child domains
   i_nchdom = p_patch%n_childdom
   IF (i_nchdom == 0) CYCLE

   ! Loop over child domains
   DO jcd = 1, i_nchdom

    p_grf_s => p_grf%p_dom(jcd)
    icid    =  p_patch%child_id(jcd)

    ! Part 1: Determine the number of grid points entering into the index lists

    ! 1a) cell points

    npoints_lbc = 0
    npoints_ubc = 0

    i_startblk = p_patch%cells%start_blk(1,1)
    i_endblk   = p_patch%cells%end_blk(min_rlcell_int,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell_int)

      DO jc = i_startidx, i_endidx

        IF (p_patch%cells%refin_ctrl(jc,jb) <= grf_bdyintp_start_c .AND. &
            p_patch%cells%refin_ctrl(jc,jb) >= grf_bdyintp_end_c .AND. &
            p_patch%cells%child_id(jc,jb) == icid) THEN
          npoints_lbc = npoints_lbc + 1
        ENDIF
        IF (p_patch%cells%refin_ctrl(jc,jb) <= grf_nudgintp_start_c+1 .AND. &
            p_patch%cells%child_id(jc,jb) == icid) THEN
          npoints_ubc = npoints_ubc + 1
        ENDIF

      ENDDO
    ENDDO

    p_grf_s%npoints_bdyintp_c = npoints_lbc
    p_grf_s%npoints_ubcintp_c = npoints_ubc

    ! 1b) edge points

    npoints_lbc = 0
    npoints_ubc = 0

    i_startblk = p_patch%edges%start_blk(1,1)
    i_endblk   = p_patch%edges%end_blk(min_rledge_int-1,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rledge_int-1)

      DO je = i_startidx, i_endidx

        IF (p_patch%edges%refin_ctrl(je,jb) <= grf_bdyintp_start_e .AND. &
            p_patch%edges%refin_ctrl(je,jb) >= grf_bdyintp_end_e .AND. &
            p_patch%edges%child_id(je,jb) == icid) THEN
          npoints_lbc = npoints_lbc + 1
        ENDIF
        IF (p_patch%edges%refin_ctrl(je,jb) <= grf_nudgintp_start_e .AND. &
            p_patch%edges%child_id(je,jb) == icid) THEN
          npoints_ubc = npoints_ubc + 1
        ENDIF

      ENDDO
    ENDDO

    p_grf_s%npoints_bdyintp_e = npoints_lbc
    p_grf_s%npoints_ubcintp_e = npoints_ubc

    ! 1c) vertex points needed for lateral-boundary interpolation of edge-based variables
    !     when gradient-based interpolation is chosen for child edges lying on a parent edge

    npoints_lbc = 0

    i_startblk = p_patch%verts%start_blk(1,1)
    i_endblk   = p_patch%verts%end_blk(min_rlvert_int-1,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlvert_int-1)

      DO jv = i_startidx, i_endidx

        ! Note: the use of grf_bdyintp_start|end_c is intended here because the respective
        ! parameters do not exist for vertices. One row more than for cells is needed here
        !
        IF (p_patch%verts%refin_ctrl(jv,jb) <= grf_bdyintp_start_c .AND. &
            p_patch%verts%refin_ctrl(jv,jb) >= grf_bdyintp_end_c-1)  THEN
          DO i = 1, p_patch%verts%num_edges(jv,jb)
            ic = p_patch%verts%cell_idx(jv,jb,i)
            ib = p_patch%verts%cell_blk(jv,jb,i)
            IF (p_patch%cells%refin_ctrl(ic,ib) ==  p_patch%verts%refin_ctrl(jv,jb) .AND. &
                p_patch%cells%child_id(ic,ib) == icid) THEN
              npoints_lbc = npoints_lbc + 1
              EXIT
            ENDIF
          ENDDO
        ENDIF

      ENDDO
    ENDDO

    p_grf_s%npoints_bdyintp_v = npoints_lbc

    ! Allocate fields containing index lists and coefficients
    npoints_lbc = MAX(1,p_grf_s%npoints_bdyintp_c)
    npoints_ubc = MAX(1,p_grf_s%npoints_ubcintp_c)
    ALLOCATE (p_grf_s%idxlist_bdyintp_c(10,npoints_lbc),p_grf_s%idxlist_ubcintp_c(10,npoints_ubc), &
              p_grf_s%blklist_bdyintp_c(10,npoints_lbc),p_grf_s%blklist_ubcintp_c(10,npoints_ubc), &
              p_grf_s%coeff_bdyintp_c(10,2,npoints_lbc),p_grf_s%coeff_ubcintp_c(10,2,npoints_ubc), &
              p_grf_s%dist_pc2cc_bdy(4,2,npoints_lbc),  p_grf_s%dist_pc2cc_ubc(4,2,npoints_ubc),   &
              inv_glb2loc_loc_index_c_grf(p_patch%n_patch_cells),   &
              inv_glb2loc_loc_index_c_ubc(p_patch%n_patch_cells), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocation of cell index lists failed')
    ENDIF

    inv_glb2loc_loc_index_c_grf = -1
    inv_glb2loc_loc_index_c_ubc = -1
    p_grf_s%idxlist_bdyintp_c = -1
    p_grf_s%blklist_bdyintp_c = -1
    p_grf_s%coeff_bdyintp_c = 0._wp
    p_grf_s%dist_pc2cc_bdy = 0._wp
    p_grf_s%idxlist_ubcintp_c = -1
    p_grf_s%blklist_ubcintp_c = -1
    p_grf_s%coeff_ubcintp_c = 0._wp
    p_grf_s%dist_pc2cc_ubc = 0._wp

    npoints_lbc = MAX(1,p_grf_s%npoints_bdyintp_e)
    npoints_ubc = MAX(1,p_grf_s%npoints_ubcintp_e)
    ALLOCATE (p_grf_s%idxlist_bdyintp_e(15,npoints_lbc),p_grf_s%idxlist_ubcintp_e(15,npoints_ubc), &
              p_grf_s%blklist_bdyintp_e(15,npoints_lbc),p_grf_s%blklist_ubcintp_e(15,npoints_ubc), &
              p_grf_s%coeff_bdyintp_e12(12,npoints_lbc),p_grf_s%coeff_ubcintp_e12(12,npoints_ubc), &
              p_grf_s%coeff_bdyintp_e34(10,npoints_lbc),p_grf_s%coeff_ubcintp_e34(10,npoints_ubc), &
              p_grf_s%edge_vert_idx(2,npoints_lbc),     p_grf_s%dist_pe2ce(2,npoints_lbc),         &
              p_grf_s%prim_norm(2,2,npoints_lbc),                                                  &
              inv_glb2loc_loc_index_e_lbc(p_patch%n_patch_edges),                                  &
              inv_glb2loc_loc_index_e_ubc(p_patch%n_patch_edges), STAT=ist)

    inv_glb2loc_loc_index_e_lbc = -1
    inv_glb2loc_loc_index_e_ubc = -1
    p_grf_s%blklist_bdyintp_e = -1
    p_grf_s%coeff_bdyintp_e12 = 0._wp
    p_grf_s%coeff_bdyintp_e34 = 0._wp
    p_grf_s%edge_vert_idx = -1
    p_grf_s%prim_norm = 0._wp
    p_grf_s%blklist_ubcintp_e = -1
    p_grf_s%coeff_ubcintp_e12 = 0._wp
    p_grf_s%coeff_ubcintp_e34 = 0._wp
    p_grf_s%dist_pe2ce = 0._wp

    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocation of edge index lists failed')
    ENDIF

    p_grf_s%idxlist_bdyintp_e(:,:) = 0
    p_grf_s%idxlist_ubcintp_e(:,:) = 0

    npoints_lbc = MAX(1,p_grf_s%npoints_bdyintp_v)
    ALLOCATE (p_grf_s%idxlist_rbfintp_v(6,npoints_lbc),p_grf_s%blklist_rbfintp_v(6,npoints_lbc), &
              p_grf_s%coeff_rbf_v(6,2,npoints_lbc), inv_ind_v(nproma,p_patch%nblks_v), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocation of vertex index lists failed')
    ENDIF

    p_grf_s%idxlist_rbfintp_v = -1
    p_grf_s%blklist_rbfintp_v = -1
    p_grf_s%coeff_rbf_v = 0._wp

    ! Part 2: Compute index and coefficient lists

    ! 2a) lateral/upper boundary interpolation for cells
    !
    i_startblk = p_patch%cells%start_blk(1,1)
    i_endblk   = p_patch%cells%end_blk(min_rlcell_int,i_nchdom)

    icount_lbc = 0
    icount_ubc = 0

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell_int)

      DO jc = i_startidx, i_endidx

        IF (p_patch%cells%refin_ctrl(jc,jb) <= grf_bdyintp_start_c .AND. &
            p_patch%cells%refin_ctrl(jc,jb) >= grf_bdyintp_end_c .AND. &
            p_patch%cells%child_id(jc,jb) == icid) THEN
          icount_lbc = icount_lbc + 1
          p_grf_s%idxlist_bdyintp_c(1,icount_lbc) = jc
          p_grf_s%blklist_bdyintp_c(1,icount_lbc) = jb
          inv_glb2loc_loc_index_c_grf(idx_1d(jc, jb)) = icount_lbc
          p_grf_s%idxlist_bdyintp_c(2:10,icount_lbc)   = p_int%rbf_c2grad_idx(2:10,jc,jb)
          p_grf_s%blklist_bdyintp_c(2:10,icount_lbc)   = p_int%rbf_c2grad_blk(2:10,jc,jb)
          p_grf_s%coeff_bdyintp_c(1:10,1:2,icount_lbc) = p_int%rbf_c2grad_coeff(1:10,1:2,jc,jb)
          p_grf_s%dist_pc2cc_bdy(1:4,1:2,icount_lbc)   = p_grf_s%grf_dist_pc2cc(jc,1:4,1:2,jb)
        ENDIF
        IF (p_patch%cells%refin_ctrl(jc,jb) <= grf_nudgintp_start_c+1 .AND. &
            p_patch%cells%child_id(jc,jb) == icid) THEN
          icount_ubc = icount_ubc + 1
          p_grf_s%idxlist_ubcintp_c(1,icount_ubc) = jc
          p_grf_s%blklist_ubcintp_c(1,icount_ubc) = jb
          inv_glb2loc_loc_index_c_ubc(idx_1d(jc, jb)) = icount_ubc
          p_grf_s%idxlist_ubcintp_c(2:10,icount_ubc)   = p_int%rbf_c2grad_idx(2:10,jc,jb)
          p_grf_s%blklist_ubcintp_c(2:10,icount_ubc)   = p_int%rbf_c2grad_blk(2:10,jc,jb)
          p_grf_s%coeff_ubcintp_c(1:10,1:2,icount_ubc) = p_int%rbf_c2grad_coeff(1:10,1:2,jc,jb)
          p_grf_s%dist_pc2cc_ubc(1:4,1:2,icount_ubc)   = p_grf_s%grf_dist_pc2cc(jc,1:4,1:2,jb)
        ENDIF

      ENDDO
    ENDDO

    ! 2b) lateral/upper boundary interpolation for vertices
    !     this is done before the edges because of the remapping needed for the edge-vertex connectivity
    !
    i_startblk = p_patch%verts%start_blk(1,1)
    i_endblk   = p_patch%verts%end_blk(min_rlvert_int-1,i_nchdom)

    icount_lbc = 0
    inv_ind_v(:,:) = -1

    DO jb = i_startblk, i_endblk

      CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlvert_int-1)

      DO jv = i_startidx, i_endidx

        ! Note: the use of grf_bdyintp_start|end_c is intended here because the respective
        ! parameters do not exist for vertices. One row more than for cells is needed here
        !
        IF (p_patch%verts%refin_ctrl(jv,jb) <= grf_bdyintp_start_c .AND. &
            p_patch%verts%refin_ctrl(jv,jb) >= grf_bdyintp_end_c-1)  THEN
          !
          ! The following lines are quite awkward, but verts%child_id is not available
          ! from the grid generator because it has never been needed before
          !
          lprocess = .FALSE.
          DO i = 1, p_patch%verts%num_edges(jv,jb)
            ic = p_patch%verts%cell_idx(jv,jb,i)
            ib = p_patch%verts%cell_blk(jv,jb,i)
            IF (p_patch%cells%refin_ctrl(ic,ib) ==  p_patch%verts%refin_ctrl(jv,jb) .AND. &
                p_patch%cells%child_id(ic,ib) == icid) THEN
              lprocess = .TRUE.
              EXIT
            ENDIF
          ENDDO
          IF (lprocess) THEN
            icount_lbc = icount_lbc + 1
            p_grf_s%idxlist_rbfintp_v(1:6,icount_lbc)      = p_int%rbf_vec_idx_v(1:6,jv,jb)
            p_grf_s%blklist_rbfintp_v(1:6,icount_lbc)      = p_int%rbf_vec_blk_v(1:6,jv,jb)
            p_grf_s%coeff_rbf_v      (1:6,1:2,icount_lbc)  = p_int%rbf_vec_coeff_v(1:6,1:2,jv,jb)
            inv_ind_v(jv,jb) = icount_lbc
          ENDIF
        ENDIF
        
      ENDDO
    ENDDO


    ! 2c) lateral/upper boundary interpolation for edges
    !
    i_startblk = p_patch%edges%start_blk(1,1)
    i_endblk   = p_patch%edges%end_blk(min_rledge_int-1,i_nchdom)

    icount_lbc = 0
    icount_ubc = 0

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rledge_int-1)

      DO je = i_startidx, i_endidx

        IF (p_patch%edges%refin_ctrl(je,jb) <= grf_bdyintp_start_e .AND. &
            p_patch%edges%refin_ctrl(je,jb) >= grf_bdyintp_end_e .AND. &
            p_patch%edges%child_id(je,jb) == icid) THEN
          icount_lbc = icount_lbc + 1
          lfound(:,:) = .FALSE.
          p_grf_s%idxlist_bdyintp_e(1,icount_lbc) = je
          p_grf_s%blklist_bdyintp_e(1,icount_lbc) = jb
          inv_glb2loc_loc_index_e_lbc(idx_1d(je,jb)) = icount_lbc
          DO i = 1,p_grf_s%grf_vec_stencil_1a(je,jb)
            IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == je .AND. &
                p_grf_s%grf_vec_blk_1a(je,i,jb) == jb) THEN
              p_grf_s%coeff_bdyintp_e12(1,icount_lbc)  = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ENDIF
          ENDDO
          DO i = 1,p_grf_s%grf_vec_stencil_1b(je,jb)
            IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == je .AND. &
                p_grf_s%grf_vec_blk_1b(je,i,jb) == jb) THEN
              p_grf_s%coeff_bdyintp_e12(7,icount_lbc)  = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ENDIF
          ENDDO
          p_grf_s%coeff_bdyintp_e34(1,icount_lbc)  = p_grf_s%grf_vec_coeff_2a(1,je,jb)
          p_grf_s%coeff_bdyintp_e34(6,icount_lbc)  = p_grf_s%grf_vec_coeff_2b(1,je,jb)
          DO i = 1,p_grf_s%grf_vec_stencil_1a(je,jb)
            IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == p_grf_s%grf_vec_ind_2a(je,2,jb) .AND. &
                p_grf_s%grf_vec_blk_1a(je,i,jb) == p_grf_s%grf_vec_blk_2a(je,2,jb)) THEN
              p_grf_s%idxlist_bdyintp_e(2,icount_lbc) = p_grf_s%grf_vec_ind_2a(je,2,jb)
              p_grf_s%blklist_bdyintp_e(2,icount_lbc) = p_grf_s%grf_vec_blk_2a(je,2,jb)
              p_grf_s%coeff_bdyintp_e34(2,icount_lbc) = p_grf_s%grf_vec_coeff_2a(2,je,jb)
              p_grf_s%coeff_bdyintp_e12(2,icount_lbc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ELSE IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == p_grf_s%grf_vec_ind_2a(je,3,jb) .AND. &
                     p_grf_s%grf_vec_blk_1a(je,i,jb) == p_grf_s%grf_vec_blk_2a(je,3,jb)) THEN
              p_grf_s%idxlist_bdyintp_e(2,icount_lbc) = p_grf_s%grf_vec_ind_2a(je,3,jb)
              p_grf_s%blklist_bdyintp_e(2,icount_lbc) = p_grf_s%grf_vec_blk_2a(je,3,jb)
              p_grf_s%coeff_bdyintp_e34(2,icount_lbc) = p_grf_s%grf_vec_coeff_2a(3,je,jb)
              p_grf_s%coeff_bdyintp_e12(2,icount_lbc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ENDIF
            IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == p_grf_s%grf_vec_ind_2b(je,2,jb) .AND. &
                p_grf_s%grf_vec_blk_1a(je,i,jb) == p_grf_s%grf_vec_blk_2b(je,2,jb)) THEN
              p_grf_s%idxlist_bdyintp_e(6,icount_lbc) = p_grf_s%grf_vec_ind_2b(je,2,jb)
              p_grf_s%blklist_bdyintp_e(6,icount_lbc) = p_grf_s%grf_vec_blk_2b(je,2,jb)
              p_grf_s%coeff_bdyintp_e34(7,icount_lbc) = p_grf_s%grf_vec_coeff_2b(2,je,jb)
              p_grf_s%coeff_bdyintp_e12(3,icount_lbc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ELSE IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == p_grf_s%grf_vec_ind_2b(je,3,jb) .AND. &
                     p_grf_s%grf_vec_blk_1a(je,i,jb) == p_grf_s%grf_vec_blk_2b(je,3,jb)) THEN
              p_grf_s%idxlist_bdyintp_e(6,icount_lbc) = p_grf_s%grf_vec_ind_2b(je,3,jb)
              p_grf_s%blklist_bdyintp_e(6,icount_lbc) = p_grf_s%grf_vec_blk_2b(je,3,jb)
              p_grf_s%coeff_bdyintp_e34(7,icount_lbc) = p_grf_s%grf_vec_coeff_2b(3,je,jb)
              p_grf_s%coeff_bdyintp_e12(3,icount_lbc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ENDIF
          ENDDO
          DO i = 1,p_grf_s%grf_vec_stencil_1b(je,jb)
            IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == p_grf_s%grf_vec_ind_2a(je,2,jb) .AND. &
                p_grf_s%grf_vec_blk_1b(je,i,jb) == p_grf_s%grf_vec_blk_2a(je,2,jb)) THEN
              p_grf_s%idxlist_bdyintp_e(3,icount_lbc) = p_grf_s%grf_vec_ind_2a(je,2,jb)
              p_grf_s%blklist_bdyintp_e(3,icount_lbc) = p_grf_s%grf_vec_blk_2a(je,2,jb)
              p_grf_s%coeff_bdyintp_e34(3,icount_lbc) = p_grf_s%grf_vec_coeff_2a(2,je,jb)
              p_grf_s%coeff_bdyintp_e12(8,icount_lbc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ELSE IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == p_grf_s%grf_vec_ind_2a(je,3,jb) .AND. &
                     p_grf_s%grf_vec_blk_1b(je,i,jb) == p_grf_s%grf_vec_blk_2a(je,3,jb)) THEN
              p_grf_s%idxlist_bdyintp_e(3,icount_lbc) = p_grf_s%grf_vec_ind_2a(je,3,jb)
              p_grf_s%blklist_bdyintp_e(3,icount_lbc) = p_grf_s%grf_vec_blk_2a(je,3,jb)
              p_grf_s%coeff_bdyintp_e34(3,icount_lbc) = p_grf_s%grf_vec_coeff_2a(3,je,jb)
              p_grf_s%coeff_bdyintp_e12(8,icount_lbc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ENDIF
            IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == p_grf_s%grf_vec_ind_2b(je,2,jb) .AND. &
                p_grf_s%grf_vec_blk_1b(je,i,jb) == p_grf_s%grf_vec_blk_2b(je,2,jb)) THEN
              p_grf_s%idxlist_bdyintp_e(7,icount_lbc) = p_grf_s%grf_vec_ind_2b(je,2,jb)
              p_grf_s%blklist_bdyintp_e(7,icount_lbc) = p_grf_s%grf_vec_blk_2b(je,2,jb)
              p_grf_s%coeff_bdyintp_e34(8,icount_lbc) = p_grf_s%grf_vec_coeff_2b(2,je,jb)
              p_grf_s%coeff_bdyintp_e12(9,icount_lbc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ELSE IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == p_grf_s%grf_vec_ind_2b(je,3,jb) .AND. &
                     p_grf_s%grf_vec_blk_1b(je,i,jb) == p_grf_s%grf_vec_blk_2b(je,3,jb)) THEN
              p_grf_s%idxlist_bdyintp_e(7,icount_lbc) = p_grf_s%grf_vec_ind_2b(je,3,jb)
              p_grf_s%blklist_bdyintp_e(7,icount_lbc) = p_grf_s%grf_vec_blk_2b(je,3,jb)
              p_grf_s%coeff_bdyintp_e34(8,icount_lbc) = p_grf_s%grf_vec_coeff_2b(3,je,jb)
              p_grf_s%coeff_bdyintp_e12(9,icount_lbc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ENDIF
          ENDDO
          p_grf_s%idxlist_bdyintp_e(4:5,icount_lbc)  = p_grf_s%grf_vec_ind_2a(je,4:5,jb)
          p_grf_s%blklist_bdyintp_e(4:5,icount_lbc)  = p_grf_s%grf_vec_blk_2a(je,4:5,jb)
          p_grf_s%idxlist_bdyintp_e(8:9,icount_lbc)  = p_grf_s%grf_vec_ind_2b(je,4:5,jb)
          p_grf_s%blklist_bdyintp_e(8:9,icount_lbc)  = p_grf_s%grf_vec_blk_2b(je,4:5,jb)
          p_grf_s%coeff_bdyintp_e34(4:5,icount_lbc)  = p_grf_s%grf_vec_coeff_2a(4:5,je,jb)
          p_grf_s%coeff_bdyintp_e34(9:10,icount_lbc) = p_grf_s%grf_vec_coeff_2b(4:5,je,jb)
          ic = 9
          DO i = 1,6
            IF (.NOT. lfound(1,i)) THEN
              IF (p_grf_s%idxlist_bdyintp_e(6,icount_lbc) == 0) THEN ! this is the case for refin_ctrl = -1
                p_grf_s%idxlist_bdyintp_e(6,icount_lbc) = p_grf_s%grf_vec_ind_1a(je,i,jb)
                p_grf_s%blklist_bdyintp_e(6,icount_lbc) = p_grf_s%grf_vec_blk_1a(je,i,jb)
                p_grf_s%coeff_bdyintp_e12(3,icount_lbc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              ELSE
                ic = ic+1
                p_grf_s%idxlist_bdyintp_e(ic,icount_lbc) = p_grf_s%grf_vec_ind_1a(je,i,jb)
                p_grf_s%blklist_bdyintp_e(ic,icount_lbc) = p_grf_s%grf_vec_blk_1a(je,i,jb)
                p_grf_s%coeff_bdyintp_e12(ic-6,icount_lbc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              ENDIF
            ENDIF
          ENDDO
          DO i = 1,6
            IF (.NOT. lfound(2,i)) THEN
              IF (p_grf_s%idxlist_bdyintp_e(7,icount_lbc) == 0) THEN ! this is the case for refin_ctrl = -1
                p_grf_s%idxlist_bdyintp_e(7,icount_lbc) = p_grf_s%grf_vec_ind_1b(je,i,jb)
                p_grf_s%blklist_bdyintp_e(7,icount_lbc) = p_grf_s%grf_vec_blk_1b(je,i,jb)
                p_grf_s%coeff_bdyintp_e12(9,icount_lbc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              ELSE
                ic = ic+1
                p_grf_s%idxlist_bdyintp_e(ic,icount_lbc) = p_grf_s%grf_vec_ind_1b(je,i,jb)
                p_grf_s%blklist_bdyintp_e(ic,icount_lbc) = p_grf_s%grf_vec_blk_1b(je,i,jb)
                p_grf_s%coeff_bdyintp_e12(ic-3,icount_lbc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              ENDIF
            ENDIF
          ENDDO
          p_grf_s%dist_pe2ce(1:2,icount_lbc)         = p_grf_s%grf_dist_pe2ce(je,1:2,jb)
          p_grf_s%prim_norm(1:2,1,icount_lbc)        = p_patch%edges%primal_normal_vert(je,jb,1:2)%v1
          p_grf_s%prim_norm(1:2,2,icount_lbc)        = p_patch%edges%primal_normal_vert(je,jb,1:2)%v2

          iv = p_patch%edges%vertex_idx(je,jb,1)
          ib = p_patch%edges%vertex_blk(je,jb,1)
          p_grf_s%edge_vert_idx(1,icount_lbc)        = inv_ind_v(iv,ib)

          iv = p_patch%edges%vertex_idx(je,jb,2)
          ib = p_patch%edges%vertex_blk(je,jb,2)
          p_grf_s%edge_vert_idx(2,icount_lbc)        = inv_ind_v(iv,ib)
        ENDIF

        IF (p_patch%edges%refin_ctrl(je,jb) <= grf_nudgintp_start_e .AND. &
            p_patch%edges%child_id(je,jb) == icid) THEN
          icount_ubc = icount_ubc + 1
          lfound(:,:) = .FALSE.
          p_grf_s%idxlist_ubcintp_e(1,icount_ubc) = je
          p_grf_s%blklist_ubcintp_e(1,icount_ubc) = jb
          inv_glb2loc_loc_index_e_ubc(idx_1d(je,jb)) = icount_ubc
          DO i = 1,p_grf_s%grf_vec_stencil_1a(je,jb)
            IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == je .AND. &
                p_grf_s%grf_vec_blk_1a(je,i,jb) == jb) THEN
              p_grf_s%coeff_ubcintp_e12(1,icount_ubc)  = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ENDIF
          ENDDO
          DO i = 1,p_grf_s%grf_vec_stencil_1b(je,jb)
            IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == je .AND. &
                p_grf_s%grf_vec_blk_1b(je,i,jb) == jb) THEN
              p_grf_s%coeff_ubcintp_e12(7,icount_ubc)  = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ENDIF
          ENDDO
          p_grf_s%coeff_ubcintp_e34(1,icount_ubc)  = p_grf_s%grf_vec_coeff_2a(1,je,jb)
          p_grf_s%coeff_ubcintp_e34(6,icount_ubc)  = p_grf_s%grf_vec_coeff_2b(1,je,jb)
          DO i = 1,p_grf_s%grf_vec_stencil_1a(je,jb)
            IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == p_grf_s%grf_vec_ind_2a(je,2,jb) .AND. &
                p_grf_s%grf_vec_blk_1a(je,i,jb) == p_grf_s%grf_vec_blk_2a(je,2,jb)) THEN
              p_grf_s%idxlist_ubcintp_e(2,icount_ubc) = p_grf_s%grf_vec_ind_2a(je,2,jb)
              p_grf_s%blklist_ubcintp_e(2,icount_ubc) = p_grf_s%grf_vec_blk_2a(je,2,jb)
              p_grf_s%coeff_ubcintp_e34(2,icount_ubc) = p_grf_s%grf_vec_coeff_2a(2,je,jb)
              p_grf_s%coeff_ubcintp_e12(2,icount_ubc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ELSE IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == p_grf_s%grf_vec_ind_2a(je,3,jb) .AND. &
                     p_grf_s%grf_vec_blk_1a(je,i,jb) == p_grf_s%grf_vec_blk_2a(je,3,jb)) THEN
              p_grf_s%idxlist_ubcintp_e(2,icount_ubc) = p_grf_s%grf_vec_ind_2a(je,3,jb)
              p_grf_s%blklist_ubcintp_e(2,icount_ubc) = p_grf_s%grf_vec_blk_2a(je,3,jb)
              p_grf_s%coeff_ubcintp_e34(2,icount_ubc) = p_grf_s%grf_vec_coeff_2a(3,je,jb)
              p_grf_s%coeff_ubcintp_e12(2,icount_ubc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ENDIF
            IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == p_grf_s%grf_vec_ind_2b(je,2,jb) .AND. &
                p_grf_s%grf_vec_blk_1a(je,i,jb) == p_grf_s%grf_vec_blk_2b(je,2,jb)) THEN
              p_grf_s%idxlist_ubcintp_e(6,icount_ubc) = p_grf_s%grf_vec_ind_2b(je,2,jb)
              p_grf_s%blklist_ubcintp_e(6,icount_ubc) = p_grf_s%grf_vec_blk_2b(je,2,jb)
              p_grf_s%coeff_ubcintp_e34(7,icount_ubc) = p_grf_s%grf_vec_coeff_2b(2,je,jb)
              p_grf_s%coeff_ubcintp_e12(3,icount_ubc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ELSE IF (p_grf_s%grf_vec_ind_1a(je,i,jb) == p_grf_s%grf_vec_ind_2b(je,3,jb) .AND. &
                     p_grf_s%grf_vec_blk_1a(je,i,jb) == p_grf_s%grf_vec_blk_2b(je,3,jb)) THEN
              p_grf_s%idxlist_ubcintp_e(6,icount_ubc) = p_grf_s%grf_vec_ind_2b(je,3,jb)
              p_grf_s%blklist_ubcintp_e(6,icount_ubc) = p_grf_s%grf_vec_blk_2b(je,3,jb)
              p_grf_s%coeff_ubcintp_e34(7,icount_ubc) = p_grf_s%grf_vec_coeff_2b(3,je,jb)
              p_grf_s%coeff_ubcintp_e12(3,icount_ubc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
              lfound(1,i) = .TRUE.
            ENDIF
          ENDDO
          DO i = 1,p_grf_s%grf_vec_stencil_1b(je,jb)
            IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == p_grf_s%grf_vec_ind_2a(je,2,jb) .AND. &
                p_grf_s%grf_vec_blk_1b(je,i,jb) == p_grf_s%grf_vec_blk_2a(je,2,jb)) THEN
              p_grf_s%idxlist_ubcintp_e(3,icount_ubc) = p_grf_s%grf_vec_ind_2a(je,2,jb)
              p_grf_s%blklist_ubcintp_e(3,icount_ubc) = p_grf_s%grf_vec_blk_2a(je,2,jb)
              p_grf_s%coeff_ubcintp_e34(3,icount_ubc) = p_grf_s%grf_vec_coeff_2a(2,je,jb)
              p_grf_s%coeff_ubcintp_e12(8,icount_ubc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ELSE IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == p_grf_s%grf_vec_ind_2a(je,3,jb) .AND. &
                     p_grf_s%grf_vec_blk_1b(je,i,jb) == p_grf_s%grf_vec_blk_2a(je,3,jb)) THEN
              p_grf_s%idxlist_ubcintp_e(3,icount_ubc) = p_grf_s%grf_vec_ind_2a(je,3,jb)
              p_grf_s%blklist_ubcintp_e(3,icount_ubc) = p_grf_s%grf_vec_blk_2a(je,3,jb)
              p_grf_s%coeff_ubcintp_e34(3,icount_ubc) = p_grf_s%grf_vec_coeff_2a(3,je,jb)
              p_grf_s%coeff_ubcintp_e12(8,icount_ubc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ENDIF
            IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == p_grf_s%grf_vec_ind_2b(je,2,jb) .AND. &
                p_grf_s%grf_vec_blk_1b(je,i,jb) == p_grf_s%grf_vec_blk_2b(je,2,jb)) THEN
              p_grf_s%idxlist_ubcintp_e(7,icount_ubc) = p_grf_s%grf_vec_ind_2b(je,2,jb)
              p_grf_s%blklist_ubcintp_e(7,icount_ubc) = p_grf_s%grf_vec_blk_2b(je,2,jb)
              p_grf_s%coeff_ubcintp_e34(8,icount_ubc) = p_grf_s%grf_vec_coeff_2b(2,je,jb)
              p_grf_s%coeff_ubcintp_e12(9,icount_ubc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ELSE IF (p_grf_s%grf_vec_ind_1b(je,i,jb) == p_grf_s%grf_vec_ind_2b(je,3,jb) .AND. &
                     p_grf_s%grf_vec_blk_1b(je,i,jb) == p_grf_s%grf_vec_blk_2b(je,3,jb)) THEN
              p_grf_s%idxlist_ubcintp_e(7,icount_ubc) = p_grf_s%grf_vec_ind_2b(je,3,jb)
              p_grf_s%blklist_ubcintp_e(7,icount_ubc) = p_grf_s%grf_vec_blk_2b(je,3,jb)
              p_grf_s%coeff_ubcintp_e34(8,icount_ubc) = p_grf_s%grf_vec_coeff_2b(3,je,jb)
              p_grf_s%coeff_ubcintp_e12(9,icount_ubc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
              lfound(2,i) = .TRUE.
            ENDIF
          ENDDO
          p_grf_s%idxlist_ubcintp_e(4:5,icount_ubc)  = p_grf_s%grf_vec_ind_2a(je,4:5,jb)
          p_grf_s%blklist_ubcintp_e(4:5,icount_ubc)  = p_grf_s%grf_vec_blk_2a(je,4:5,jb)
          p_grf_s%idxlist_ubcintp_e(8:9,icount_ubc)  = p_grf_s%grf_vec_ind_2b(je,4:5,jb)
          p_grf_s%blklist_ubcintp_e(8:9,icount_ubc)  = p_grf_s%grf_vec_blk_2b(je,4:5,jb)
          p_grf_s%coeff_ubcintp_e34(4:5,icount_ubc)  = p_grf_s%grf_vec_coeff_2a(4:5,je,jb)
          p_grf_s%coeff_ubcintp_e34(9:10,icount_ubc) = p_grf_s%grf_vec_coeff_2b(4:5,je,jb)
          ic = 9
          DO i = 1,6
            IF (.NOT. lfound(1,i)) THEN
              ic = ic+1
              p_grf_s%idxlist_ubcintp_e(ic,icount_ubc) = p_grf_s%grf_vec_ind_1a(je,i,jb)
              p_grf_s%blklist_ubcintp_e(ic,icount_ubc) = p_grf_s%grf_vec_blk_1a(je,i,jb)
              p_grf_s%coeff_ubcintp_e12(ic-6,icount_ubc) = p_grf_s%grf_vec_coeff_1a(i,je,jb)
            ENDIF
          ENDDO
          DO i = 1,6
            IF (.NOT. lfound(2,i)) THEN
              ic = ic+1
              p_grf_s%idxlist_ubcintp_e(ic,icount_ubc) = p_grf_s%grf_vec_ind_1b(je,i,jb)
              p_grf_s%blklist_ubcintp_e(ic,icount_ubc) = p_grf_s%grf_vec_blk_1b(je,i,jb)
              p_grf_s%coeff_ubcintp_e12(ic-3,icount_ubc) = p_grf_s%grf_vec_coeff_1b(i,je,jb)
            ENDIF
          ENDDO

        ENDIF

      ENDDO
    ENDDO

    DEALLOCATE (inv_ind_v)

    ! Generate the communication patterns for lateral and upper boundary
    ! interpolation
    n = COUNT(inv_glb2loc_loc_index_c_grf > 0)
    CALL init_glb2loc_index_lookup(inv_glb2loc, p_patch%n_patch_cells_g)
    CALL set_inner_glb_index(inv_glb2loc, p_patch%cells%decomp_info%glb_index, &
                             inv_glb2loc_loc_index_c_grf)
    ALLOCATE(owner_local(n), glb_index(n))
    DO i = 1, p_patch%n_patch_cells
      IF (inv_glb2loc_loc_index_c_grf(i) > 0) THEN
        owner_local(inv_glb2loc_loc_index_c_grf(i)) = &
          p_patch%cells%decomp_info%owner_local(i)
        glb_index(inv_glb2loc_loc_index_c_grf(i)) = &
          p_patch%cells%decomp_info%glb_index(i)
      END IF
    END DO
    CALL generate_interpol_pattern(p_patch_all(icid)%n_patch_cells, &
      &                            p_patch_all(icid)%cells%parent_glb_idx, &
      &                            p_patch_all(icid)%cells%parent_glb_blk, &
      &                            p_patch_all(jg)%cells%decomp_info%owner_dist_dir, &
      &                            inv_glb2loc, n, &
      &                            owner_local, glb_index, scal_grf_select_func, &
      &                            p_patch_all(icid), &
      &                            p_patch_all(icid)%comm_pat_coll_interpol_scal_grf)
    CALL deallocate_glb2loc_index_lookup(inv_glb2loc)
    DEALLOCATE(owner_local, glb_index)

    n = COUNT(inv_glb2loc_loc_index_c_ubc > 0)
    CALL init_glb2loc_index_lookup(inv_glb2loc, p_patch%n_patch_cells_g)
    CALL set_inner_glb_index(inv_glb2loc, p_patch%cells%decomp_info%glb_index, &
                             inv_glb2loc_loc_index_c_ubc)
    ALLOCATE(owner_local(n), glb_index(n))
    DO i = 1, p_patch%n_patch_cells
      IF (inv_glb2loc_loc_index_c_ubc(i) > 0) THEN
        owner_local(inv_glb2loc_loc_index_c_ubc(i)) = &
          p_patch%cells%decomp_info%owner_local(i)
        glb_index(inv_glb2loc_loc_index_c_ubc(i)) = &
          p_patch%cells%decomp_info%glb_index(i)
      END IF
    END DO
    CALL generate_interpol_pattern(p_patch_all(icid)%n_patch_cells, &
      &                            p_patch_all(icid)%cells%parent_glb_idx, &
      &                            p_patch_all(icid)%cells%parent_glb_blk, &
      &                            p_patch_all(jg)%cells%decomp_info%owner_dist_dir, &
      &                            inv_glb2loc, n, &
      &                            owner_local, glb_index, scal_ubc_select_func, &
      &                            p_patch_all(icid), &
      &                            p_patch_all(icid)%comm_pat_coll_interpol_scal_ubc)
    CALL deallocate_glb2loc_index_lookup(inv_glb2loc)
    DEALLOCATE(owner_local, glb_index)

    n = COUNT(inv_glb2loc_loc_index_e_lbc > 0)
    CALL init_glb2loc_index_lookup(inv_glb2loc, p_patch%n_patch_edges_g)
    CALL set_inner_glb_index(inv_glb2loc, p_patch%edges%decomp_info%glb_index, &
                             inv_glb2loc_loc_index_e_lbc)
    ALLOCATE(owner_local(n), glb_index(n))
    DO i = 1, p_patch%n_patch_edges
      IF (inv_glb2loc_loc_index_e_lbc(i) > 0) THEN
        owner_local(inv_glb2loc_loc_index_e_lbc(i)) = &
          p_patch%edges%decomp_info%owner_local(i)
        glb_index(inv_glb2loc_loc_index_e_lbc(i)) = &
          p_patch%edges%decomp_info%glb_index(i)
      END IF
    END DO
    CALL generate_interpol_pattern(p_patch_all(icid)%n_patch_edges, &
      &                            p_patch_all(icid)%edges%parent_glb_idx, &
      &                            p_patch_all(icid)%edges%parent_glb_blk, &
      &                            p_patch_all(jg)%edges%decomp_info%owner_dist_dir, &
      &                            inv_glb2loc, n, &
      &                            owner_local, glb_index, vec_grf_select_func, &
      &                            p_patch_all(icid), &
      &                            p_patch_all(icid)%comm_pat_coll_interpol_vec_grf)
    CALL deallocate_glb2loc_index_lookup(inv_glb2loc)
    DEALLOCATE(owner_local, glb_index)

    n = COUNT(inv_glb2loc_loc_index_e_ubc > 0)
    CALL init_glb2loc_index_lookup(inv_glb2loc, p_patch%n_patch_edges_g)
    CALL set_inner_glb_index(inv_glb2loc, p_patch%edges%decomp_info%glb_index, &
                             inv_glb2loc_loc_index_e_ubc)
    ALLOCATE(owner_local(n), glb_index(n))
    DO i = 1, p_patch%n_patch_edges
      IF (inv_glb2loc_loc_index_e_ubc(i) > 0) THEN
        owner_local(inv_glb2loc_loc_index_e_ubc(i)) = &
          p_patch%edges%decomp_info%owner_local(i)
        glb_index(inv_glb2loc_loc_index_e_ubc(i)) = &
          p_patch%edges%decomp_info%glb_index(i)
      END IF
    END DO
    CALL generate_interpol_pattern(p_patch_all(icid)%n_patch_edges, &
      &                            p_patch_all(icid)%edges%parent_glb_idx, &
      &                            p_patch_all(icid)%edges%parent_glb_blk, &
      &                            p_patch_all(jg)%edges%decomp_info%owner_dist_dir, &
      &                            inv_glb2loc, n, &
      &                            owner_local, glb_index, vec_ubc_select_func, &
      &                            p_patch_all(icid), &
      &                            p_patch_all(icid)%comm_pat_coll_interpol_vec_ubc)
    CALL deallocate_glb2loc_index_lookup(inv_glb2loc)
    DEALLOCATE(owner_local, glb_index)

    DEALLOCATE(inv_glb2loc_loc_index_c_grf, inv_glb2loc_loc_index_c_ubc, &
      &        inv_glb2loc_loc_index_e_lbc, inv_glb2loc_loc_index_e_ubc)

   ENDDO ! child domains

   ! Index lists for grid points on which tendency fields for lateral boundary interpolation
   ! need to be computed at parent level (including neighbor plus halo points for
   ! gradient computation / interpolation stencil)

   ! cell points

    npoints_lbc_src = 0

    i_startblk = p_patch%cells%start_blk(1,1)
    i_endblk   = p_patch%cells%end_blk(min_rlcell,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell)

      DO jc = i_startidx, i_endidx

        IF (p_patch%cells%refin_ctrl(jc,jb) <= grf_bdyintp_start_c .AND. &
            p_patch%cells%refin_ctrl(jc,jb) >= grf_bdyintp_end_c-1 ) THEN
          npoints_lbc_src = npoints_lbc_src + 1
        ENDIF
        IF (p_patch%cells%refin_ctrl(jc,jb) == 0) THEN
          neighbor_points: DO i = 2, 10
            ic = p_int%rbf_c2grad_idx(i,jc,jb)
            ib = p_int%rbf_c2grad_blk(i,jc,jb)
            IF (p_patch%cells%refin_ctrl(ic,ib) < 0) THEN
              npoints_lbc_src = npoints_lbc_src + 1
              EXIT neighbor_points ! to avoid multiple counts of the same grid point
            ENDIF
          ENDDO neighbor_points
        ENDIF

      ENDDO
    ENDDO

    p_grf%npoints_bdyintp_src_c = npoints_lbc_src

    ! edge points

    npoints_lbc_src = 0

    i_startblk = p_patch%edges%start_blk(1,1)
    i_endblk   = p_patch%edges%end_blk(min_rledge,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rledge)

      DO je = i_startidx, i_endidx

        IF (p_patch%edges%refin_ctrl(je,jb) <= grf_bdyintp_start_e .AND. &
            p_patch%edges%refin_ctrl(je,jb) >= grf_bdyintp_end_e-2 ) THEN
          npoints_lbc_src = npoints_lbc_src + 1
        ENDIF
        IF (p_patch%edges%refin_ctrl(je,jb) == 0) THEN
          iv1 = p_patch%edges%vertex_idx(je,jb,1)
          ib1 = p_patch%edges%vertex_blk(je,jb,1)
          iv2 = p_patch%edges%vertex_idx(je,jb,2)
          ib2 = p_patch%edges%vertex_blk(je,jb,2)
          IF (p_patch%verts%refin_ctrl(iv1,ib1) < 0 .OR. p_patch%verts%refin_ctrl(iv2,ib2) < 0) THEN
            npoints_lbc_src = npoints_lbc_src + 1
          ENDIF
        ENDIF

      ENDDO
    ENDDO

    p_grf%npoints_bdyintp_src_e = npoints_lbc_src

    ALLOCATE (p_grf%idxlist_bdyintp_src_c(MAX(1,p_grf%npoints_bdyintp_src_c)),         &
              p_grf%blklist_bdyintp_src_c(MAX(1,p_grf%npoints_bdyintp_src_c)),         &
              p_grf%idxlist_bdyintp_src_e(MAX(1,p_grf%npoints_bdyintp_src_e)),         &
              p_grf%blklist_bdyintp_src_e(MAX(1,p_grf%npoints_bdyintp_src_e)), STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocation of parent-level index lists failed')
    ENDIF

    ! Now compute index lists

    ! cell points

    i_startblk = p_patch%cells%start_blk(1,1)
    i_endblk   = p_patch%cells%end_blk(min_rlcell,i_nchdom)

    icount_lbc_src = 0

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell)

      DO jc = i_startidx, i_endidx

        IF (p_patch%cells%refin_ctrl(jc,jb) <= grf_bdyintp_start_c .AND. &
            p_patch%cells%refin_ctrl(jc,jb) >= grf_bdyintp_end_c-1 ) THEN
          icount_lbc_src = icount_lbc_src + 1
          p_grf%idxlist_bdyintp_src_c(icount_lbc_src) = jc
          p_grf%blklist_bdyintp_src_c(icount_lbc_src) = jb
        ENDIF
        IF (p_patch%cells%refin_ctrl(jc,jb) == 0) THEN
          neighbor_points_l2: DO i = 2, 10
            ic = p_int%rbf_c2grad_idx(i,jc,jb)
            ib = p_int%rbf_c2grad_blk(i,jc,jb)
            IF (p_patch%cells%refin_ctrl(ic,ib) < 0) THEN
              icount_lbc_src = icount_lbc_src + 1
              p_grf%idxlist_bdyintp_src_c(icount_lbc_src) = jc
              p_grf%blklist_bdyintp_src_c(icount_lbc_src) = jb
              EXIT neighbor_points_l2 ! to avoid multiple counts of the same grid point
            ENDIF
          ENDDO neighbor_points_l2
        ENDIF

      ENDDO
    ENDDO

    ! edge points

    i_startblk = p_patch%edges%start_blk(1,1)
    i_endblk   = p_patch%edges%end_blk(min_rledge,i_nchdom)

    icount_lbc_src = 0

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rledge)

      DO je = i_startidx, i_endidx

        IF (p_patch%edges%refin_ctrl(je,jb) <= grf_bdyintp_start_e .AND. &
            p_patch%edges%refin_ctrl(je,jb) >= grf_bdyintp_end_e-2 ) THEN
            icount_lbc_src = icount_lbc_src + 1
            p_grf%idxlist_bdyintp_src_e(icount_lbc_src) = je
            p_grf%blklist_bdyintp_src_e(icount_lbc_src) = jb
        ENDIF
        IF (p_patch%edges%refin_ctrl(je,jb) == 0) THEN
          iv1 = p_patch%edges%vertex_idx(je,jb,1)
          ib1 = p_patch%edges%vertex_blk(je,jb,1)
          iv2 = p_patch%edges%vertex_idx(je,jb,2)
          ib2 = p_patch%edges%vertex_blk(je,jb,2)
          IF (p_patch%verts%refin_ctrl(iv1,ib1) < 0 .OR. p_patch%verts%refin_ctrl(iv2,ib2) < 0) THEN
            icount_lbc_src = icount_lbc_src + 1
            p_grf%idxlist_bdyintp_src_e(icount_lbc_src) = je
            p_grf%blklist_bdyintp_src_e(icount_lbc_src) = jb
          ENDIF
        ENDIF

      ENDDO
    ENDDO


   ! Compute mask fields needed for feedback computations executed at the (normal) parent grid level
   !
   ALLOCATE (p_grf%mask_ovlp_c(nproma,p_patch%nblks_c,i_nchdom),p_grf%mask_ovlp_ch(nproma,p_patch%nblks_c,i_nchdom),&
             p_grf%mask_ovlp_e(nproma,p_patch%nblks_e,i_nchdom),p_grf%mask_ovlp_v (nproma,p_patch%nblks_v,i_nchdom),&
             STAT=ist)
   IF (ist /= SUCCESS) THEN
     CALL finish (routine,'allocation of mask fields failed')
   ENDIF

   ! Initialization of mask fields with .FALSE.
   p_grf%mask_ovlp_c (:,:,:) = .FALSE.
   p_grf%mask_ovlp_ch(:,:,:) = .FALSE.
   p_grf%mask_ovlp_e (:,:,:) = .FALSE.
   p_grf%mask_ovlp_v (:,:,:) = .FALSE.

   ! Loop over child domains
   DO jcd = 1, i_nchdom

     icid    =  p_patch%child_id(jcd)

     ! masks for cells; 'ch' is for intermediate calculations that need to include a halo row
     !
     i_startblk = p_patch%cells%start_blk(1,1)
     i_endblk   = p_patch%cells%end_blk(min_rlcell,i_nchdom)

     DO jb = i_startblk, i_endblk

       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, min_rlcell)

       DO jc = i_startidx,i_endidx

         IF (p_patch%cells%refin_ctrl(jc,jb) <= grf_fbk_start_c .AND.                   &
             p_patch%cells%child_id(jc,jb) == icid) p_grf%mask_ovlp_c(jc,jb,jcd) = .TRUE.

         IF (p_patch%cells%refin_ctrl(jc,jb) <= grf_fbk_start_c + 1 .AND.                &
             p_patch%cells%child_id(jc,jb) == icid) p_grf%mask_ovlp_ch(jc,jb,jcd) = .TRUE.

       ENDDO
     ENDDO

     ! mask for edges
     !
     i_startblk = p_patch%edges%start_blk(1,1)
     i_endblk   = p_patch%edges%end_blk(min_rledge,i_nchdom)

     DO jb = i_startblk, i_endblk

       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, min_rledge)

       DO je = i_startidx,i_endidx

         IF (p_patch%edges%refin_ctrl(je,jb) <= grf_fbk_start_e .AND.                   &
             p_patch%edges%child_id(je,jb) == icid) p_grf%mask_ovlp_e(je,jb,jcd) = .TRUE.

       ENDDO
     ENDDO

     ! mask for vertices; note that a workaround is needed for determining the overlap with
     ! the right child domain because the child_id field is not available for vertices
     !
     i_startblk = p_patch%verts%start_blk(1,1)
     i_endblk   = p_patch%verts%end_blk(min_rlvert,i_nchdom)

     DO jb = i_startblk, i_endblk

       CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, min_rlvert)

       DO jv = i_startidx,i_endidx

         ! As the lateral boundary of the nest overlap region is excluded anyway, we can just
         ! choose any adjacent cell to determine the child id
         ic = p_patch%verts%cell_idx(jv,jb,1)
         ib = p_patch%verts%cell_blk(jv,jb,1)

         IF (ic <= 0 .OR. ib <= 0) CYCLE ! to avoid bound-checking errors
         IF (p_patch%verts%refin_ctrl(jv,jb) <= grf_fbk_start_c .AND.                   &
             p_patch%cells%child_id(ic,ib) == icid) p_grf%mask_ovlp_v(jv,jb,jcd) = .TRUE.

       ENDDO
     ENDDO

   ENDDO ! loop over child domains

  ENDDO ! domain ID

  ! When everything works, the 'old' index and coefficient lists can be deallocated here;
  ! in the deallocation routine, they need to be replaced with the new ones in this case

CONTAINS

  SUBROUTINE generate_interpol_pattern(n_patch_cve, parent_glb_idx, &
    &                                  parent_glb_blk, parent_owner_dist_dir, &
    &                                  glb2loc, n_src, owner_local_src, &
    &                                  glb_index_src, &
    &                                  select_func, p_patch, &
    &                                  comm_pat_coll_interpol)

    INTEGER, INTENT(IN) :: n_patch_cve, parent_glb_idx(:,:), parent_glb_blk(:,:)
    TYPE(t_dist_dir), INTENT(IN) :: parent_owner_dist_dir
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc
    INTEGER, INTENT(IN) :: n_src, owner_local_src(:), glb_index_src(:)
    TYPE(t_patch), INTENT(IN) :: p_patch
    CLASS(t_comm_pattern_collection), POINTER, INTENT(OUT) :: &
      comm_pat_coll_interpol
    INTERFACE
      FUNCTION select_func(idx, blk, n, p_patch) RESULT(p)
        USE mo_model_domain, ONLY : t_patch
        INTEGER, INTENT(IN) :: idx, blk, n
        TYPE(t_patch), INTENT(IN) :: p_patch
        LOGICAL :: p
      END FUNCTION select_func
    END INTERFACE

    TYPE(t_p_comm_pattern) :: comm_pat_interpol(4)
    INTEGER :: i, n, idx, blk
    INTEGER :: owner_local_dst(n_patch_cve), glb_index_dst(n_patch_cve)

    !--------------------------------------------------------------------
    ! Cells

    ! For our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    DO n = 1, 4

      glb_index_dst(:) = -1

      ! Communication to nest boundary points includes halo points in order to save subsequent synchronization
      DO i = 1, n_patch_cve
        idx = idx_no(i)
        blk = blk_no(i)
        IF (select_func(idx, blk, n, p_patch)) &
          glb_index_dst(i) = idx_1d(parent_glb_idx(idx, blk), &
            &                       parent_glb_blk(idx, blk))
      ENDDO

      owner_local_dst(:) = dist_dir_get_owners(parent_owner_dist_dir, &
        &                                      glb_index_dst(:), &
        &                                      glb_index_dst(:) /= -1)

      ! Set up communication pattern
      CALL setup_comm_pattern(n_patch_cve, owner_local_dst, glb_index_dst,  &
        &                     glb2loc, n_src, owner_local_src, glb_index_src, &
        &                     comm_pat_interpol(n)%p)

    ENDDO

    CALL setup_comm_pattern_collection(comm_pat_interpol, &
      &                                comm_pat_coll_interpol)

  END SUBROUTINE generate_interpol_pattern

END SUBROUTINE create_grf_index_lists

SUBROUTINE destruct_interpol_patterns(p_patch_all)
  !
  TYPE(t_patch), TARGET, INTENT(INOUT)         :: p_patch_all(n_dom_start:)

  TYPE(t_patch), POINTER :: p_patch
  INTEGER :: jg, jcd, i_nchdom, icid

!-----------------------------------------------------------------------

  DO jg = 1, n_dom

    p_patch => p_patch_all(jg)
    ! number of child domains
    i_nchdom = p_patch%n_childdom
    IF (i_nchdom == 0) CYCLE
    ! Loop over child domains
    DO jcd = 1, i_nchdom

      icid    =  p_patch%child_id(jcd)
      CALL delete_comm_pattern_collection( &
        p_patch_all(icid)%comm_pat_coll_interpol_scal_grf)
      CALL delete_comm_pattern_collection( &
        p_patch_all(icid)%comm_pat_coll_interpol_scal_ubc)
      CALL delete_comm_pattern_collection( &
        p_patch_all(icid)%comm_pat_coll_interpol_vec_grf)
      CALL delete_comm_pattern_collection( &
        p_patch_all(icid)%comm_pat_coll_interpol_vec_ubc)

    ENDDO ! child domains

  ENDDO ! domain ID

END SUBROUTINE destruct_interpol_patterns

!-------------------------------------------------------------------------
!
!! Deallocation of fields needed for grid refinement for a single state.
!!
SUBROUTINE deallocate_grf_state( ptr_patch, ptr_grf )
!
TYPE(t_patch), INTENT(IN) :: ptr_patch

TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: ptr_grf

CHARACTER(LEN=*), PARAMETER :: routine = modname//':deallocate_grf_state'
TYPE(t_gridref_single_state), POINTER :: ptr_int    => NULL()

INTEGER :: jcd, i_nchdom
INTEGER :: ist

!-----------------------------------------------------------------------

  i_nchdom = MAX(1,ptr_patch%n_childdom)

  DO jcd = 1, i_nchdom

    ptr_int => ptr_grf%p_dom(jcd)

    ! index arrays

    DEALLOCATE (ptr_int%grf_vec_ind_1a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_ind_1a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_ind_1b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_ind_1b failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_ind_2a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_ind_2a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_ind_2b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_ind_2b failed')
    ENDIF

    ! block index arrays

    DEALLOCATE (ptr_int%grf_vec_blk_1a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_ind_1a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_blk_1b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_ind_1b failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_blk_2a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_ind_2a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_blk_2b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_ind_2b failed')
    ENDIF

    ! arrays holding the actual number of cells/edges available for reconstruction

    DEALLOCATE (ptr_int%grf_vec_stencil_1a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_stencil_1a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_stencil_1b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_stencil_1b failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_stencil_2a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_stencil_2a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_stencil_2b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_stencil_2b failed')
    ENDIF

    ! arrays holding the interpolation coefficients for reconstruction related to grid refinement

    DEALLOCATE (ptr_int%grf_vec_coeff_1a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_coeff_1a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_coeff_1b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_coeff_1b failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_coeff_2a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_coeff_2a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_coeff_2b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_vec_coeff_2b failed')
    ENDIF

    DEALLOCATE (ptr_int%grf_dist_pc2cc, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_dist_pc2cc failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_dist_pe2ce, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,                       &
      &             'deallocation for grf_dist_pe2ce failed')
    ENDIF

  ENDDO

  DEALLOCATE (ptr_grf%fbk_wgt_aw, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine,                       &
      &             'deallocation for fbk_wgt_aw failed')
  ENDIF
  DEALLOCATE (ptr_grf%fbk_wgt_bln, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine,                       &
      &             'deallocation for fbk_wgt_bln failed')
  ENDIF
  DEALLOCATE (ptr_grf%fbk_wgt_e, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine,                       &
      &             'deallocation for fbk_wgt_e failed')
  ENDIF

  DEALLOCATE (ptr_grf%p_dom, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (routine,                       &
      &             'deallocation for p_dom failed')
  ENDIF

END SUBROUTINE deallocate_grf_state

!-------------------------------------------------------------------------
!
!! Deallocation of fields needed for grid refinement.
!!
SUBROUTINE destruct_2d_gridref_state( p_patch, p_patch_local_parent, p_grf_state, &
  &                                   p_grf_state_local_parent )
  !
  TYPE(t_patch),         INTENT(IN)    :: p_patch(n_dom_start:)
  TYPE(t_patch),         INTENT(IN)    :: p_patch_local_parent(n_dom_start+1:)
  TYPE(t_gridref_state), INTENT(INOUT) :: p_grf_state(n_dom_start:)
  TYPE(t_gridref_state), INTENT(INOUT) :: p_grf_state_local_parent(n_dom_start+1:)

  CHARACTER(LEN=*), PARAMETER :: routine = modname//':destruct_2d_gridref_state'
  INTEGER :: jg

!-----------------------------------------------------------------------

  CALL message(routine, 'start to destruct 2d_gridref state')

  !
  ! deallocate gridref_state
  !
  DO jg = n_dom_start, n_dom
    CALL deallocate_grf_state(p_patch(jg), p_grf_state(jg))
  END DO

  DO jg = n_dom_start+1, n_dom
    CALL deallocate_grf_state(p_patch_local_parent(jg), p_grf_state_local_parent(jg))
  END DO

  CALL message (routine, 'destruction of 2d_gridref state finished')

END SUBROUTINE destruct_2d_gridref_state


!-------------------------------------------------------------------------
END MODULE mo_grf_intp_state
