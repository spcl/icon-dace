! Contains the implementation of interpolation and reconstruction
! routines used by the shallow water model, including the RBF
! reconstruction routines.
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

#ifdef __xlC__
@PROCESS smp=noopt
@PROCESS noopt
#endif
#ifdef __PGI
!pgi$g opt=1
#endif

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_intp_rbf_coeffs
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
USE mo_exception,           ONLY: message, message_text, finish
USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int
USE mo_model_domain,        ONLY: t_patch, t_tangent_vectors
USE mo_math_types,          ONLY: t_cartesian_coordinates
USE mo_math_utilities,      ONLY: gvec2cvec, arc_length_v
USE mo_math_utility_solvers, ONLY: solve_chol_v, choldec_v
USE mo_parallel_config,     ONLY: nproma
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_vec_dim_e, rbf_vec_dim_v,     &
  &                               rbf_c2grad_dim, rbf_vec_kern_c, rbf_vec_kern_e,  &
  &                               rbf_vec_kern_v, rbf_vec_scale_c, rbf_vec_scale_e,&
  &                               rbf_vec_scale_v
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_idx
USE mo_grid_geometry_info,  ONLY: sphere_geometry
USE mo_physical_constants,  ONLY: earth_radius
USE mo_math_constants,      ONLY: deg2rad  ! = pi/180
USE mo_mpi,                 ONLY: get_my_mpi_work_id  ! only for debugging
IMPLICIT NONE

!> level of output verbosity
INTEGER, PARAMETER  :: dbg_level = 0

PRIVATE

PUBLIC :: rbf_vec_index_cell, rbf_c2grad_index, rbf_vec_compute_coeff_cell,     &
          & rbf_compute_coeff_c2grad, rbf_vec_index_vertex, rbf_vec_index_edge, &
          & rbf_vec_compute_coeff_vertex, rbf_vec_compute_coeff_edge,           &
          & gen_index_list_radius

!> module name
CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_rbf_coeffs'

CONTAINS

#include "intp_functions.inc"

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Part 1: routines for RBF vector reconstruction and related functions
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!! This routine initializes the indexes used to define the stencil.
!!
!! This routine initializes the indexes used to define the stencil
!! of the vector RBF interpolator. The stencil is cell based and includes
!! a variable number of edges (rbf_vec_dim) around each cell, in order
!! to reconstruct a vector at the cell center.
!!
SUBROUTINE rbf_vec_index_cell( ptr_patch, ptr_int )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), INTENT(inout) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER  :: nblks_c
INTEGER  :: jc, jb           ! loop indices
INTEGER  :: ilc, ibc         ! line and block index of neighbour cell

INTEGER  :: i_startblk       ! start block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index
REAL(wp) :: z_stencil(UBOUND(ptr_int%rbf_vec_stencil_c,1),UBOUND(ptr_int%rbf_vec_stencil_c,2))

!--------------------------------------------------------------------

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c
  !
  ! The stencil consists of 9 edges, the cell edges are taken
  ! and the edges of the nearest neighbours are taken
  !
  ! This stencil is not available when refin_ctrl=1
  i_startblk = ptr_patch%cells%start_blk(2,1)

  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, 2)
!$NEC ivdep
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      !
      ! get the global line and block indices of the edges of each neighbor
      ! cell and store them in the rbf_vec_idx_c and rbf_vec_blk_c components
      !
      ! Sorry for the clumsy code, but it was not possible to get it vectorized
      ! on the NEC with nested do loops and UNROLL directives
      !
      ! First neighbor cell
      ilc = ptr_patch%cells%neighbor_idx(jc,jb,1)
      ibc = ptr_patch%cells%neighbor_blk(jc,jb,1)
      ptr_int%rbf_vec_idx_c(1,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,1)
      ptr_int%rbf_vec_blk_c(1,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,1)
      ptr_int%rbf_vec_idx_c(2,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,2)
      ptr_int%rbf_vec_blk_c(2,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,2)
      ptr_int%rbf_vec_idx_c(3,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,3)
      ptr_int%rbf_vec_blk_c(3,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,3)

      ! Second neighbor cell
      ilc = ptr_patch%cells%neighbor_idx(jc,jb,2)
      ibc = ptr_patch%cells%neighbor_blk(jc,jb,2)
      ptr_int%rbf_vec_idx_c(4,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,1)
      ptr_int%rbf_vec_blk_c(4,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,1)
      ptr_int%rbf_vec_idx_c(5,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,2)
      ptr_int%rbf_vec_blk_c(5,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,2)
      ptr_int%rbf_vec_idx_c(6,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,3)
      ptr_int%rbf_vec_blk_c(6,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,3)

      ! Third neighbor cell
      ilc = ptr_patch%cells%neighbor_idx(jc,jb,3)
      ibc = ptr_patch%cells%neighbor_blk(jc,jb,3)
      ptr_int%rbf_vec_idx_c(7,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,1)
      ptr_int%rbf_vec_blk_c(7,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,1)
      ptr_int%rbf_vec_idx_c(8,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,2)
      ptr_int%rbf_vec_blk_c(8,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,2)
      ptr_int%rbf_vec_idx_c(9,jc,jb) = ptr_patch%cells%edge_idx(ilc,ibc,3)
      ptr_int%rbf_vec_blk_c(9,jc,jb) = ptr_patch%cells%edge_blk(ilc,ibc,3)

      ! take care of cells at patch boundaries, then the value of 
      ! ptr_int%rbf_vec_stencil_c might be smaller than rbf_vec_dim_c:
!$NEC unroll(9)
      ptr_int%rbf_vec_stencil_c(jc,jb) = COUNT(ptr_int%rbf_vec_idx_c(1:9,jc,jb) /= 0)

    END DO

  END DO
  
  DO jb = 1, rbf_vec_dim_c
    CALL sync_idx(SYNC_C, SYNC_E, ptr_patch, ptr_int%rbf_vec_idx_c(jb,:,:), &
      &           ptr_int%rbf_vec_blk_c(jb,:,:))
  ENDDO
  z_stencil(:,:) = REAL(ptr_int%rbf_vec_stencil_c(:,:),wp)
  CALL sync_patch_array(SYNC_C,ptr_patch,z_stencil)
  ptr_int%rbf_vec_stencil_c(:,:) = NINT(z_stencil(:,:))

END SUBROUTINE rbf_vec_index_cell

!! This routine creates the index list needed for RBF reconstruction
!! of gradients of scalar values at cell centers, using the values at the
!! surrounding cell centers.
!!
SUBROUTINE rbf_c2grad_index( ptr_patch, ptr_int )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER  :: jc, jb           ! loop indices
INTEGER  :: ilc, ibc         ! line and block index of neighbour cell
INTEGER  :: ic1, ic2         ! index counters

INTEGER  :: i_startblk       ! start block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index

! Pointers to neighbor indices/blocks
INTEGER, DIMENSION(:,:,:), POINTER :: inidx, inblk
INTEGER :: rl_start, rl_end, i_nchdom, i_endblk

!--------------------------------------------------------------------

  inidx => ptr_patch%cells%neighbor_idx 
  inblk => ptr_patch%cells%neighbor_blk

  ! values for the blocking
  rl_start = 2
  rl_end = min_rlcell_int
  
  i_nchdom   = MAX(1,ptr_patch%n_childdom)
  i_startblk = ptr_patch%cells%start_blk(rl_start,1)
  i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)    

  !
  ! The stencil has a total of 10 cells, namely the local cell,
  ! its neighbors, and the neighbors of the neighbors (excluding the local cell)
  !
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)
!$NEC ivdep
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! Local point
      ptr_int%rbf_c2grad_idx(1,jc,jb) = jc
      ptr_int%rbf_c2grad_blk(1,jc,jb) = jb

      ! Direct neighbors
      ptr_int%rbf_c2grad_idx(2,jc,jb) = inidx(jc,jb,1)
      ptr_int%rbf_c2grad_blk(2,jc,jb) = inblk(jc,jb,1)
      ptr_int%rbf_c2grad_idx(5,jc,jb) = inidx(jc,jb,2)
      ptr_int%rbf_c2grad_blk(5,jc,jb) = inblk(jc,jb,2)
      ptr_int%rbf_c2grad_idx(8,jc,jb) = inidx(jc,jb,3)
      ptr_int%rbf_c2grad_blk(8,jc,jb) = inblk(jc,jb,3)

      ! Neighbors of neighbor 1 (excluding local cell)
      ilc = inidx(jc,jb,1)
      ibc = inblk(jc,jb,1)
      IF (inidx(ilc,ibc,1) == jc .AND. inblk(ilc,ibc,1) == jb) THEN
        ic1 = 2
        ic2 = 3
      ELSE IF (inidx(ilc,ibc,2) == jc .AND. inblk(ilc,ibc,2) == jb) THEN
        ic1 = 1
        ic2 = 3
      ELSE
        ic1 = 1
        ic2 = 2
      ENDIF
      ptr_int%rbf_c2grad_idx(3,jc,jb) = inidx(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_blk(3,jc,jb) = inblk(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_idx(4,jc,jb) = inidx(ilc,ibc,ic2)
      ptr_int%rbf_c2grad_blk(4,jc,jb) = inblk(ilc,ibc,ic2)

      ! Neighbors of neighbor 2 (excluding local cell)
      ilc = inidx(jc,jb,2)
      ibc = inblk(jc,jb,2)
      IF (inidx(ilc,ibc,1) == jc .AND. inblk(ilc,ibc,1) == jb) THEN
        ic1 = 2
        ic2 = 3
      ELSE IF (inidx(ilc,ibc,2) == jc .AND. inblk(ilc,ibc,2) == jb) THEN
        ic1 = 1
        ic2 = 3
      ELSE
        ic1 = 1
        ic2 = 2
      ENDIF
      ptr_int%rbf_c2grad_idx(6,jc,jb) = inidx(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_blk(6,jc,jb) = inblk(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_idx(7,jc,jb) = inidx(ilc,ibc,ic2)
      ptr_int%rbf_c2grad_blk(7,jc,jb) = inblk(ilc,ibc,ic2)

      ! Neighbors of neighbor 3 (excluding local cell)
      ilc = inidx(jc,jb,3)
      ibc = inblk(jc,jb,3)
      IF (inidx(ilc,ibc,1) == jc .AND. inblk(ilc,ibc,1) == jb) THEN
        ic1 = 2
        ic2 = 3
      ELSE IF (inidx(ilc,ibc,2) == jc .AND. inblk(ilc,ibc,2) == jb) THEN
        ic1 = 1
        ic2 = 3
      ELSE
        ic1 = 1
        ic2 = 2
      ENDIF
      ptr_int%rbf_c2grad_idx(9,jc,jb) = inidx(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_blk(9,jc,jb) = inblk(ilc,ibc,ic1)
      ptr_int%rbf_c2grad_idx(10,jc,jb) = inidx(ilc,ibc,ic2)
      ptr_int%rbf_c2grad_blk(10,jc,jb) = inblk(ilc,ibc,ic2)
    END DO
  END DO

  DO jb = 1, rbf_c2grad_dim
    CALL sync_idx(SYNC_C, SYNC_C, ptr_patch, ptr_int%rbf_c2grad_idx(jb,:,:), &
      &           ptr_int%rbf_c2grad_blk(jb,:,:))
  ENDDO

END SUBROUTINE rbf_c2grad_index


  !! Description: 
  !!   For a given search_radius (in m)
  !!   find all triangles with indices (jc, jb), whose circumcenter lies
  !!   in a distance smaller than the search_radius to the circumcenter of
  !!   the center triangle (ic, ib).
  !! Output:
  !!   ptr_int contains the derived type variable 'cell_environ'; here, the following entries are set:
  !!      cell_environ%max_nmbr_nghbr_cells
  !!      cell_environ%nmbr_nghbr_cells (ic, ib)
  !!      cell_environ%idx       (ic, ib, l)
  !!      cell_environ%blk       (ic, ib, l)
  !!      cell_environ%area_norm (ic, ib, l)
  !!      cell_environ%radius          ! only for later control checks
  !!      cell_environ%max_nmbr_iter   ! only for later control checks
  !!   The purpose of such a set of indices usually is averaging for postprocessing 
  !!   (currently for SDI and LPI)
  !!
  !! Method:
  !!   iteratively inspect the neighbors of each cell vortices.
  !!   The maximum iteration depth is limited by max_nmbr_iter.
  !!   This subroutine is a generalization of the subroutine rbf_c2grad_index.
  !!
  SUBROUTINE gen_index_list_radius( ptr_int, ptr_patch, jg, search_radius, max_nmbr_iter )

    IMPLICIT NONE

    TYPE(t_int_state),     INTENT(INOUT) :: ptr_int
    TYPE(t_patch), TARGET, INTENT(INOUT) :: ptr_patch   ! is modified by sync_idx
    INTEGER,               INTENT(IN)    :: jg
    REAL (wp),             INTENT(IN)    :: search_radius
    INTEGER,               INTENT(IN)    :: max_nmbr_iter  ! max. number of iterations.
            ! =1: leads to a reproducible setup since two triangle rows are exchanged in any case.
            ! >1: possible, but does not longer deliver reproducible results.

    INTEGER :: i_rlstart, i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: max_nmbr_nghbr_cells2
    INTEGER :: il, l, k, i, kv
    INTEGER :: ic, ib
    INTEGER :: jc, jb
    INTEGER :: jc2, jb2
    INTEGER :: idx_v, blk_v
    INTEGER :: idx_start_list
    INTEGER :: iter
    LOGICAL :: flag_cell_found
    REAL(wp) :: dist2, radius2
    REAL(wp) :: dx, dy
    REAL(wp) :: cos_phi
    REAL(wp) :: sum
  
    INTEGER, ALLOCATABLE :: idx_new(:)
    INTEGER, ALLOCATABLE :: blk_new(:)

    INTEGER :: ist
    LOGICAL :: flag_out

    ! without halo or boundary  points:
    i_rlstart = 2     ! see subr. rbf_c2grad_index
    !i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    ptr_int%cell_environ%radius        = search_radius ! just for later validation purposes
    ptr_int%cell_environ%max_nmbr_iter = max_nmbr_iter ! "

    ptr_int%cell_environ%max_nmbr_nghbr_cells = 0

    radius2 = search_radius * search_radius

    ! Allocations

    max_nmbr_nghbr_cells2 = ptr_int%cell_environ%nmbr_nghbr_cells_alloc * 2    ! factor 2 is for 'security reasons'

    ALLOCATE( idx_new( ptr_patch%geometry_info%cell_type * max_nmbr_nghbr_cells2), STAT=ist )
    IF (ist /= SUCCESS)  CALL finish (modname//':gen_index_list_radius', 'allocation for idx_new failed')

    ALLOCATE( blk_new( ptr_patch%geometry_info%cell_type * max_nmbr_nghbr_cells2), STAT=ist )
    IF (ist /= SUCCESS)  CALL finish (modname//':gen_index_list_radius', 'allocation for blk_new failed')
       ! cell_type=3 or 6, i.e. should work for triangle and hexagon grid

!$OMP PARALLEL
!$OMP DO PRIVATE(ib,i_startidx,i_endidx,ic,idx_start_list,cos_phi,iter,il,l,jc,jb,kv,  &
!$OMP            idx_v,blk_v,k,jc2,jb2,dx,dy,dist2,idx_new,blk_new,i,flag_cell_found,sum) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = i_startblk, i_endblk

      CALL get_indices_c (ptr_patch, ib, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO ic = i_startidx, i_endidx

        ! (ic,ib) is the center cell

        ! this center cell is stored, too:
        ! (so to speak, it defines the 0-th 'layer')
        ptr_int%cell_environ%nmbr_nghbr_cells( ic, ib)    = 1
        ptr_int%cell_environ%idx             ( ic, ib, 1) = ic
        ptr_int%cell_environ%blk             ( ic, ib, 1) = ib

        idx_start_list = 1

        cos_phi = COS( ptr_patch%cells%center(ic,ib)%lat )

        ITER_LOOP: DO iter=1, max_nmbr_iter
          ! search iteratively the next 'layer' of cells around the former 'layer'

          ! determine next direct neighbours:
          il=0

          DO l= idx_start_list, ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib)

            jc = ptr_int%cell_environ%idx( ic, ib, l)
            jb = ptr_int%cell_environ%blk( ic, ib, l)

            IF ( ptr_patch%cells%decomp_info%owner_mask(jc,jb) ) THEN
              ! = this is a cell in the interior and therefore possesses neigbouring cells

              DO kv=1,3
                idx_v = ptr_patch%cells%vertex_idx(jc,jb,kv)
                blk_v = ptr_patch%cells%vertex_blk(jc,jb,kv)

                DO k=1, ptr_patch%verts%num_edges(idx_v, blk_v)
                  jc2 = ptr_patch%verts%cell_idx( idx_v, blk_v, k)
                  jb2 = ptr_patch%verts%cell_blk( idx_v, blk_v, k)

                  ! dist2 = distance (squared!) between the cell centers of (jc2,jb2) and (ic,ib)
                  !    (the following is only a quasi-cartesian *approximation* and only applicable over 
                  !     short distances and for not too coarse grids)
                  dx = ( ptr_patch%cells%center(jc2,jb2)%lon - ptr_patch%cells%center(ic,ib)%lon )   &
                    &        * deg2rad * earth_radius * cos_phi 
                  dy = ( ptr_patch%cells%center(jc2,jb2)%lat - ptr_patch%cells%center(ic,ib)%lat )   &
                    &        * deg2rad * earth_radius
                  dist2 = dx*dx + dy*dy 

                  IF ( dist2 <= radius2 ) THEN
                    il=il+1
                    IF ( il > max_nmbr_nghbr_cells2 ) THEN
                      CALL finish( modname//':gen_index_list_radius',&
                        &          "max_nmbr_nghbr_cells2 is chosen too small" )
                    END IF
                    idx_new(il) = jc2
                    blk_new(il) = jb2
                  END IF

                END DO
              END DO

            END IF

          END DO

          IF ( il == 0 ) THEN
            ! I didn't find new cells within the search_radius ==> leave the iteration process
            EXIT ITER_LOOP

          ELSE

            ! for the next iteration:
            idx_start_list = ptr_int%cell_environ%nmbr_nghbr_cells( ic, ib) + 1

            ! store only new cells:
            DO i=1, il

              flag_cell_found = .FALSE.

              ! Is this i-th cell already contained in the list 'cell_environ'?
              TEST_LOOP: DO l=1, ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib)
                ! l=1, ...==> this search inspects every cell found so far.
                ! A more clever search wouldn't inspect cells from the most inner layers...

                IF (      ( idx_new(i) == ptr_int%cell_environ%idx( ic, ib, l) )    &
                  & .AND. ( blk_new(i) == ptr_int%cell_environ%blk( ic, ib, l) ) ) THEN
                  ! I know this cell already ==> do not store again
                  flag_cell_found = .TRUE.

                  EXIT  TEST_LOOP
                END IF
              END DO  TEST_LOOP

              IF ( .NOT. flag_cell_found ) THEN
                ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib) =      &
                ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib) + 1

                l = ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib)

                IF ( l > ptr_int%cell_environ%nmbr_nghbr_cells_alloc ) THEN
                   CALL finish( modname//':gen_index_list_radius',&
                     &          "nmbr_nghbr_cells_alloc is chosen too small" )
                END IF

                ptr_int%cell_environ%idx( ic, ib, l ) = idx_new(i)
                ptr_int%cell_environ%blk( ic, ib, l ) = blk_new(i)

              END IF

            END DO

          END IF

        END DO  ITER_LOOP

        ptr_int%cell_environ%max_nmbr_nghbr_cells = MAX( ptr_int%cell_environ%max_nmbr_nghbr_cells,  &
          &                                              ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib) )

        ! Set area_norm
        sum = 0.0_wp
        DO l=1, ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib)
          jc2 = ptr_int%cell_environ%idx( ic, ib, l )
          jb2 = ptr_int%cell_environ%blk( ic, ib, l ) 
          ptr_int%cell_environ%area_norm( ic, ib, l) = ptr_patch%cells%area( jc2, jb2 )
          sum = sum + ptr_int%cell_environ%area_norm( ic, ib, l)
        END DO

        DO l=1, ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib)
          ptr_int%cell_environ%area_norm( ic, ib, l) = ptr_int%cell_environ%area_norm( ic, ib, l) / sum
        END DO

        ! set unused elements:
        DO l=ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib)+1, ptr_int%cell_environ%nmbr_nghbr_cells_alloc
          ptr_int%cell_environ%area_norm( ic, ib, l) = 0.0_wp
          ! just to have a defined index contained:
          ptr_int%cell_environ%idx( ic, ib, l ) = ptr_int%cell_environ%idx( ic, ib, 1 )
          ptr_int%cell_environ%blk( ic, ib, l ) = ptr_int%cell_environ%blk( ic, ib, 1 )
        END DO

      END DO  ! ic
    END DO  ! ib
!$OMP END PARALLEL


    ! control output
    IF ( .FALSE. ) THEN
      DO ib = i_startblk, i_endblk

        CALL get_indices_c (ptr_patch, ib, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend)

        DO ic = i_startidx, i_endidx

          !IF ( get_my_mpi_work_id() == 0 ) THEN
          !WRITE(*,'(A,6I6)') "GGG ", get_my_mpi_work_id(), jg, ic, ib, &
          !  & ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib), ptr_int%cell_environ%max_nmbr_nghbr_cells
          !END IF

          flag_out = .FALSE.

          IF ( ptr_int%cell_environ%max_nmbr_nghbr_cells /= 13 ) THEN
            flag_out = .TRUE.
          END IF
          DO l=1, ptr_int%cell_environ%max_nmbr_nghbr_cells
            IF ( ( ptr_int%cell_environ%area_norm( ic, ib, l) < 0.05 ) .OR.   &
              &  ( ptr_int%cell_environ%area_norm( ic, ib, l) > 0.09 ) ) THEN
              flag_out = .TRUE.
            END IF
          END DO

          IF ( flag_out ) THEN
            WRITE(*,'(A,4(A,I5),A,I5,I5)') "[gen_index_list_radius]: ", "ID=", get_my_mpi_work_id(),   &
              &     " jg=" , jg, " ic=", ic, " ib=", ib, &
              &     " nmbr_nghbr=", ptr_int%cell_environ%nmbr_nghbr_cells(ic,ib), ptr_int%cell_environ%max_nmbr_nghbr_cells
            DO l=1, ptr_int%cell_environ%max_nmbr_nghbr_cells
              jc = ptr_int%cell_environ%idx( ic, ib, l)
              jb = ptr_int%cell_environ%blk( ic, ib, l)
              WRITE(*,'(A,I4,A,I4,I6,I6,3F13.6)') "  ID=", get_my_mpi_work_id(), "  l=", l, jc, jb,         &
                ptr_patch%cells%center(jc,jb)%lon, ptr_patch%cells%center(jc,jb)%lat,   &
                ptr_int%cell_environ%area_norm( ic, ib, l)
            END DO
          END IF

        END DO  ! ic
      END DO  ! ib
    END IF


    DEALLOCATE( idx_new, STAT=ist )
    IF (ist /= SUCCESS)  CALL finish (modname//':gen_index_list_radius', 'deallocation for idx_new failed')
    DEALLOCATE( blk_new, STAT=ist )
    IF (ist /= SUCCESS)  CALL finish (modname//':gen_index_list_radius', 'deallocation for blk_new failed')

    DO l=1, ptr_int%cell_environ%nmbr_nghbr_cells_alloc
      CALL sync_idx(SYNC_C, SYNC_C, ptr_patch, ptr_int%cell_environ%idx(:,:,l),  &
        &                                      ptr_int%cell_environ%blk(:,:,l) )
    END DO


  END SUBROUTINE gen_index_list_radius


!! This routine initializes the indexes used to define the stencil.
!!
!! This routine initializes the indexes used to define the stencil
!! of the vector RBF interpolator e.g. for the Simpson rule. The stencil
!! is vertex based and includes a variable number of edges
!! (rbf_vec_dim_v) around each vertex, in order
!! to reconstruct a vector at the vertices of a triangle.
!!
SUBROUTINE rbf_vec_index_vertex( ptr_patch, ptr_int )

!
TYPE(t_patch), INTENT(inout) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER,  ALLOCATABLE :: ile(:), ibe(:) ! edge indices for the stencil

INTEGER :: nblks_v
INTEGER :: istencil                 ! number of edges for the stencil
INTEGER :: jv, jb                   ! loop indices
INTEGER :: jje                      ! additional loop indices
INTEGER :: iie, ipent
INTEGER  :: i_startblk       ! start block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index
LOGICAL :: ll_pent                  ! if .TRUE. vertex is center of a pentagon
INTEGER :: ist                      ! status variable
INTEGER :: jg

REAL(wp) :: z_stencil(UBOUND(ptr_int%rbf_vec_stencil_v,1),UBOUND(ptr_int%rbf_vec_stencil_v,2))

!--------------------------------------------------------------------

  ! values for the blocking
  nblks_v  = ptr_patch%nblks_v

  ! init pentagon counter for security check
  ipent = 0

  !
  ! The stencil consists of 6 edges, the edges sharing the vertex are taken.
  ! (for the central vertex of a pentagon we use only 5 edges)
  !
  ! This stencil is not available when refin_ctrl=1
  i_startblk = ptr_patch%verts%start_blk(2,1)

!$OMP PARALLEL PRIVATE(ile,ibe,ist)
  ALLOCATE( ile(rbf_vec_dim_v), ibe(rbf_vec_dim_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (modname//':rbf_vec_index_vertex',  &
      &             'allocation for ile, ibe failed')
  ENDIF
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jje,istencil,ll_pent,iie) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, 2)

    DO jv = i_startidx, i_endidx

      IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

      istencil = 6
      ll_pent = .FALSE.

      ! init counter
      iie = 0

      ! check if vertex is central vertex of a pentagon
      IF ( ptr_patch%verts%num_edges(jv,jb) == 5 ) THEN
        istencil = 5
        ll_pent = .TRUE.
!$OMP ATOMIC
        ipent = ipent+1
      END IF

      DO jje = 1, ptr_patch%verts%num_edges(jv,jb)
        !
        ! get global indices of the edges around the vertex
        !
        iie = iie+1
        ile(iie) = ptr_patch%verts%edge_idx(jv,jb,jje)
        ibe(iie) = ptr_patch%verts%edge_blk(jv,jb,jje)
      END DO

      IF (iie /= istencil) THEN
        PRINT *, "ERROR ==>  iie = ", iie,  " /= istencil = ", istencil
        CALL finish (modname//':rbf_vec_index_vertex',  &
          &             'wrong number of stencil points')
      ENDIF

      ! save number of edges for the stencil in rbf_vec_stencil_v
      ptr_int%rbf_vec_stencil_v(jv,jb) = istencil

      DO jje = 1, istencil
        !
        ! store line and block indices of edges for the stencil of vertex jv
        ! in rbf_vec_idx_v, rbf_vec_blk_v
        !
        ptr_int%rbf_vec_idx_v(jje,jv,jb) = ile(jje)
        ptr_int%rbf_vec_blk_v(jje,jv,jb) = ibe(jje)
      END DO

      ! for pentagons fill index arrays with dummy values
      IF ( ll_pent ) THEN
        DO jje = istencil+1, rbf_vec_dim_v
          ptr_int%rbf_vec_idx_v(jje,jv,jb) = ile(jje-istencil)
          ptr_int%rbf_vec_blk_v(jje,jv,jb) = ibe(jje-istencil)
        END DO
      END IF

    END DO ! end vertex loop

  END DO ! end block loop
!$OMP END DO NOWAIT

  DEALLOCATE( ile, ibe, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (modname//':rbf_vec_index_vertex',  &
      &             'deallocation for ile, ibe failed')
  ENDIF
!$OMP END PARALLEL

  DO jb = 1, rbf_vec_dim_v
    CALL sync_idx(SYNC_V, SYNC_E, ptr_patch, ptr_int%rbf_vec_idx_v(jb,:,:), &
      &           ptr_int%rbf_vec_blk_v(jb,:,:))
  ENDDO
  z_stencil(:,:) = REAL(ptr_int%rbf_vec_stencil_v(:,:),wp)
  CALL sync_patch_array(SYNC_V,ptr_patch,z_stencil)
  ptr_int%rbf_vec_stencil_v(:,:) = NINT(z_stencil(:,:))

END SUBROUTINE rbf_vec_index_vertex

!-------------------------------------------------------------------------
!
!! This routine initializes the indexes used to define the stencil.
!!
!! This routine initializes the indexes used to define the stencil
!! of the vector RBF interpolator. The stencil is edge based and
!! includes a variable number of edges (rbf_vec_dim_e) around
!! each edge, in order to reconstruct a vector at the edge midpoints.
!!
SUBROUTINE rbf_vec_index_edge(ptr_patch, ptr_int)
!
TYPE(t_patch), INTENT(inout) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER :: nblks_e
INTEGER :: istencil                 ! number of edges for the stencil
INTEGER :: je, jb                   ! loop indices

INTEGER  :: i_startblk       ! start block
INTEGER  :: i_startidx       ! start index
INTEGER  :: i_endidx         ! end index

REAL(wp) :: z_stencil(UBOUND(ptr_int%rbf_vec_stencil_e,1),UBOUND(ptr_int%rbf_vec_stencil_e,2))

!--------------------------------------------------------------------

  ! values for the blocking
  nblks_e  = ptr_patch%nblks_e

  !
  ! The stencil consists of 4 edges, the edges of two neighboring
  ! triangles, excluding the local point, are taken.
  !
  ! This stencil is not available when refin_ctrl=1
  i_startblk = ptr_patch%edges%start_blk(2,1)

  istencil = 4

  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, 2)
!$NEC ivdep
    DO je = i_startidx, i_endidx

      ! There is not much work to do because the required stencil points
      ! are already stored in quad_idx/quad_blk
      ptr_int%rbf_vec_idx_e(1,je,jb) = ptr_patch%edges%quad_idx(je,jb,1)
      ptr_int%rbf_vec_blk_e(1,je,jb) = ptr_patch%edges%quad_blk(je,jb,1)
      ptr_int%rbf_vec_idx_e(2,je,jb) = ptr_patch%edges%quad_idx(je,jb,2)
      ptr_int%rbf_vec_blk_e(2,je,jb) = ptr_patch%edges%quad_blk(je,jb,2)
      ptr_int%rbf_vec_idx_e(3,je,jb) = ptr_patch%edges%quad_idx(je,jb,3)
      ptr_int%rbf_vec_blk_e(3,je,jb) = ptr_patch%edges%quad_blk(je,jb,3)
      ptr_int%rbf_vec_idx_e(4,je,jb) = ptr_patch%edges%quad_idx(je,jb,4)
      ptr_int%rbf_vec_blk_e(4,je,jb) = ptr_patch%edges%quad_blk(je,jb,4)

      ! save number of edges for the stencil in rbf_vec_stencil_e
      ptr_int%rbf_vec_stencil_e(je,jb) = istencil

    END DO ! end edge loop

  END DO ! end block loop

  DO jb = 1, rbf_vec_dim_e
    CALL sync_idx(SYNC_E, SYNC_E, ptr_patch, ptr_int%rbf_vec_idx_e(jb,:,:), &
      &           ptr_int%rbf_vec_blk_e(jb,:,:))
  ENDDO
  ! Not really necessary, only for the case that rbf_vec_stencil_e should be changed:
  z_stencil(:,:) = REAL(ptr_int%rbf_vec_stencil_e(:,:),wp)
  CALL sync_patch_array(SYNC_E,ptr_patch,z_stencil)
  ptr_int%rbf_vec_stencil_e(:,:) = NINT(z_stencil(:,:))

END SUBROUTINE rbf_vec_index_edge

!-------------------------------------------------------------------------
!
!! This routine computes the coefficients needed for vector RBF interpolation,.
!!
!! This routine computes the coefficients needed for vector RBF interpolation,
!! which are then stored in the arrays <i>rbf_vec_coeff</i>.
!! This computation involves the inversion of the interpolation
!! matrix, which is performed by a Cholesky decomposition.
!! The Cholesky decomposition is currently implemented by a home made routine
!! which can be substituted by a call to a numerical library, if available.
!!
SUBROUTINE rbf_vec_compute_coeff_cell( ptr_patch, ptr_int )
!

!
TYPE(t_patch), INTENT(inout) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int

REAL(wp) :: cc_e1(3), cc_e2(3), cc_c(nproma,3), cc_aux(3)  ! coordinates of edge midpoints

REAL(wp) :: z_lon, z_lat          ! longitude and latitude

! 3d  normal velocity vectors at edge midpoints
REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2, z_nx3

REAL(wp) :: z_norm                ! norm of velocity vectors

REAL(wp) :: z_dist                ! distance between data points

REAL(wp) :: z_nxprod              ! scalar product of normal
                                  ! velocity vectors

REAL(wp),ALLOCATABLE :: z_rbfmat(:,:,:)    ! RBF interpolation matrix

REAL(wp),ALLOCATABLE :: z_diag(:,:)        ! diagonal of cholesky
                                           ! decomposition matrix

REAL(wp),ALLOCATABLE :: z_rbfval(:,:)      ! RBF function value

REAL(wp),ALLOCATABLE :: z_rhs1(:,:), &     ! right hand side of linear
                      & z_rhs2(:,:)        ! interpolation system in 2d

INTEGER :: nblks_c
INTEGER :: jc, jb          ! integer over cells, blocks and levels
INTEGER :: je1, je2        ! integer over edges
INTEGER :: ile1, ibe1, ile2, ibe2  ! edge indices
INTEGER :: ist             ! return value of array allocation
INTEGER :: i_startblk      ! start block
INTEGER :: i_startidx      ! start index
INTEGER :: i_endidx        ! end index
INTEGER :: i_rcstartlev    ! refinement control start level
INTEGER :: istencil(nproma) ! actual number of stencil points
INTEGER :: jg


REAL(wp) ::  checksum_u,checksum_v ! to check if sum of interpolation coefficients is correct

!--------------------------------------------------------------------

  CALL message(modname//':rbf_vec_compute_coeff_cell', '')

  jg = ptr_patch%id

  i_rcstartlev = 2

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist)
  ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_c,rbf_vec_dim_c),  &
            z_diag(nproma,rbf_vec_dim_c),                  &
            z_rbfval(nproma,rbf_vec_dim_c),                &
            z_rhs1(nproma,rbf_vec_dim_c),                  &
            z_rhs2(nproma,rbf_vec_dim_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (modname//':rbf_vec_compute_coeff_cell',  &
      &             'allocation for working arrays failed')
  ENDIF

!$OMP DO PRIVATE (jb,jc,i_startidx,i_endidx,je1,je2,istencil,      &
!$OMP             ist,ile1,ibe1,cc_e1,z_lon,z_lat,z_norm,      &
!$OMP             z_nx1,ile2,ibe2,cc_e2,cc_c,cc_aux,z_nx2,z_nxprod,z_dist, &
!$OMP             z_nx3,checksum_u,checksum_v) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each cell, build the vector RBF interpolation matrix
    !
    DO je1 = 1, rbf_vec_dim_c

      DO je2 = 1, je1

        IF (ptr_patch%geometry_info%geometry_type == sphere_geometry) THEN ! use vectorizable version
!$NEC ivdep
          DO jc = i_startidx, i_endidx

            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
              ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
              istencil(jc) = 0
              CYCLE
            ENDIF

            ! Get actual number of stencil points
            istencil(jc) = ptr_int%rbf_vec_stencil_c(jc,jb)
            !
            IF ( (je1 > istencil(jc)) .OR. (je2 > istencil(jc)) ) CYCLE
            !
            ! line and block indices for each edge je1 and je2 of RBF stencil
            !
            ile1 = ptr_int%rbf_vec_idx_c(je1,jc,jb)
            ibe1 = ptr_int%rbf_vec_blk_c(je1,jc,jb)
            ile2 = ptr_int%rbf_vec_idx_c(je2,jc,jb)
            ibe2 = ptr_int%rbf_vec_blk_c(je2,jc,jb)
            !
            ! get Cartesian coordinates and orientation vectors
            !
            cc_e1(:) = ptr_patch%edges%cartesian_center(ile1,ibe1)%x(:)
            cc_e2(:) = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
            !
            z_nx1(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)
            z_nx2(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)
            !
            ! compute dot product of normal vectors and distance between edge midpoints
            !
            z_nxprod = DOT_PRODUCT(z_nx1(jc,:),z_nx2(jc,:))
            z_dist   = arc_length_v(cc_e1,cc_e2)
            !
            ! set up interpolation matrix
            !
            IF      (rbf_vec_kern_c == 1) THEN
              z_rbfmat(jc,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_c == 3) THEN
              z_rbfmat(jc,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
            ENDIF
  
            IF (je1 > je2) z_rbfmat(jc,je2,je1) = z_rbfmat(jc,je1,je2)
          END DO
        ELSE ! use generic, unvectorizable version of arc_length_v
          DO jc = i_startidx, i_endidx

            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
              ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
              istencil(jc) = 0
              CYCLE
            ENDIF

            ! Get actual number of stencil points
            istencil(jc) = ptr_int%rbf_vec_stencil_c(jc,jb)
            !
            IF ( (je1 > istencil(jc)) .OR. (je2 > istencil(jc)) ) CYCLE
            !
            ! line and block indices for each edge je1 and je2 of RBF stencil
            !
            ile1 = ptr_int%rbf_vec_idx_c(je1,jc,jb)
            ibe1 = ptr_int%rbf_vec_blk_c(je1,jc,jb)
            ile2 = ptr_int%rbf_vec_idx_c(je2,jc,jb)
            ibe2 = ptr_int%rbf_vec_blk_c(je2,jc,jb)
            !
            ! get Cartesian coordinates and orientation vectors
            !
            cc_e1(:) = ptr_patch%edges%cartesian_center(ile1,ibe1)%x(:)
            cc_e2(:) = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
            !
            z_nx1(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)
            z_nx2(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)
            !
            ! compute dot product of normal vectors and distance between edge midpoints
            !
            z_nxprod = DOT_PRODUCT(z_nx1(jc,:),z_nx2(jc,:))
            z_dist   = arc_length_v(cc_e1,cc_e2,ptr_patch%geometry_info)
            !
            ! set up interpolation matrix
            !
            IF      (rbf_vec_kern_c == 1) THEN
              z_rbfmat(jc,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_c == 3) THEN
              z_rbfmat(jc,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
            ENDIF
  
            IF (je1 > je2) z_rbfmat(jc,je2,je1) = z_rbfmat(jc,je1,je2)
          END DO
        ENDIF
      END DO

    END DO

    ! apply Cholesky decomposition to matrix
    !
#ifdef __SX__
    CALL choldec_v(i_startidx,i_endidx,istencil,rbf_vec_dim_c,z_rbfmat,z_diag)
#else
    CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
#endif
!$NEC ivdep
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
      !
      ! Solve immediately for coefficients
      !
      !Get cartesian coordinate of the cell center
      cc_c(jc,1:3) = ptr_patch%cells%cartesian_center(jc,jb)%x(1:3)

      z_lon  = ptr_patch%cells%center(jc,jb)%lon
      z_lat  = ptr_patch%cells%center(jc,jb)%lat

      ! Zonal wind component
      CALL gvec2cvec(1._wp,0._wp,z_lon,z_lat,z_nx1(jc,1),z_nx1(jc,2),z_nx1(jc,3),ptr_patch%geometry_info)

      z_norm = SQRT( DOT_PRODUCT(z_nx1(jc,:),z_nx1(jc,:)) )
      z_nx1(jc,:)  = 1._wp/z_norm * z_nx1(jc,:)

      ! Meridional wind component
      CALL gvec2cvec(0._wp,1._wp,z_lon,z_lat,z_nx2(jc,1),z_nx2(jc,2),z_nx2(jc,3),ptr_patch%geometry_info)

      z_norm = SQRT( DOT_PRODUCT(z_nx2(jc,:),z_nx2(jc,:)) )
      z_nx2(jc,:)  = 1._wp/z_norm * z_nx2(jc,:)

    END DO

    !
    ! set up right hand side for interpolation system
    !
    DO je2 = 1, rbf_vec_dim_c
      IF (ptr_patch%geometry_info%geometry_type == sphere_geometry) THEN ! use vectorizable version
!$NEC ivdep
        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

          IF (je2 > istencil(jc)) CYCLE
          !
          ! get indices and coordinates of edge midpoints and compute distance
          ! to cell center
          !
          ile2   = ptr_int%rbf_vec_idx_c(je2,jc,jb)
          ibe2   = ptr_int%rbf_vec_blk_c(je2,jc,jb)
          !
          cc_e2(:)  = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
          cc_aux(1:3) = cc_c(jc,1:3)
          z_dist = arc_length_v(cc_aux, cc_e2)

          !
          ! get Cartesian orientation vector
          z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

          IF (rbf_vec_kern_c == 1) THEN
            z_rbfval(jc,je2) = gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_c == 3) THEN
            z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
          ENDIF
          !
          ! compute projection on target vector orientation
          !
          z_rhs1(jc,je2) = z_rbfval(jc,je2) * DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
          z_rhs2(jc,je2) = z_rbfval(jc,je2) * DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))

        END DO
      ELSE ! use generic, unvectorizable version of arc_length_v
        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

          IF (je2 > istencil(jc)) CYCLE
          !
          ! get indices and coordinates of edge midpoints and compute distance
          ! to cell center
          !
          ile2   = ptr_int%rbf_vec_idx_c(je2,jc,jb)
          ibe2   = ptr_int%rbf_vec_blk_c(je2,jc,jb)
          !
          cc_e2(:)  = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)

          z_dist = arc_length_v(cc_c(jc,:), cc_e2, ptr_patch%geometry_info)

          !
          ! get Cartesian orientation vector
          z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile2,ibe2)%x(:)

          IF (rbf_vec_kern_c == 1) THEN
            z_rbfval(jc,je2) = gaussi(z_dist,rbf_vec_scale_c(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_c == 3) THEN
            z_rbfval(jc,je2) = inv_multiq(z_dist,rbf_vec_scale_c(MAX(jg,1)))
          ENDIF
          !
          ! compute projection on target vector orientation
          !
          z_rhs1(jc,je2) = z_rbfval(jc,je2) * DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
          z_rhs2(jc,je2) = z_rbfval(jc,je2) * DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))

        END DO
      ENDIF
    END DO

    !
    ! compute vector coefficients
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_c, z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_int%rbf_vec_coeff_c(:,1,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_int%rbf_vec_coeff_c(:,1,:,jb))
#endif
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_c, z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_int%rbf_vec_coeff_c(:,2,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_int%rbf_vec_coeff_c(:,2,:,jb))
#endif

    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! Ensure that sum of interpolation coefficients is correct

      checksum_u = 0._wp
      checksum_v = 0._wp

      DO je1 = 1, istencil(jc)
        ile1   = ptr_int%rbf_vec_idx_c(je1,jc,jb)
        ibe1   = ptr_int%rbf_vec_blk_c(je1,jc,jb)

        ! get Cartesian orientation vector
        z_nx3(jc,:) = ptr_patch%edges%primal_cart_normal(ile1,ibe1)%x(:)

        checksum_u = checksum_u + ptr_int%rbf_vec_coeff_c(je1,1,jc,jb)* &
          DOT_PRODUCT(z_nx1(jc,:),z_nx3(jc,:))
        checksum_v = checksum_v + ptr_int%rbf_vec_coeff_c(je1,2,jc,jb)* &
          DOT_PRODUCT(z_nx2(jc,:),z_nx3(jc,:))
      ENDDO

      DO je1 = 1, istencil(jc)
        ptr_int%rbf_vec_coeff_c(je1,1,jc,jb) = &
          ptr_int%rbf_vec_coeff_c(je1,1,jc,jb) / checksum_u
        ptr_int%rbf_vec_coeff_c(je1,2,jc,jb) = &
          ptr_int%rbf_vec_coeff_c(je1,2,jc,jb) / checksum_v
      ENDDO

    END DO

  END DO
!$OMP END DO NOWAIT

  DEALLOCATE( z_rbfmat, z_diag, z_rbfval, z_rhs1, z_rhs2,  &
    STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (modname//':rbf_vec_compute_coeff_cell',  &
      &             'deallocation for working arrays failed')
  ENDIF
!$OMP END PARALLEL

  DO jb = 1, rbf_vec_dim_c
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_int%rbf_vec_coeff_c(jb,1,:,:))
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_int%rbf_vec_coeff_c(jb,2,:,:))
  ENDDO

! Optional debug output for RBF coefficients
#ifdef DEBUG_COEFF
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

      DO jc = i_startidx, i_endidx
      istencil(jc) = ptr_int%rbf_vec_stencil_c(jc,jb)
      WRITE(510+jg,'(2i5,25f12.6)') jb,jc,ptr_int%rbf_vec_coeff_c(1:istencil(jc),1,jc,jb)
      WRITE(510+jg,'(2i5,25f12.6)') jb,jc,ptr_int%rbf_vec_coeff_c(1:istencil(jc),2,jc,jb)
    END DO
  END DO
  CLOSE (510+jg)
#endif

END SUBROUTINE rbf_vec_compute_coeff_cell


!! This routine computes the coefficients needed for reconstructing the gradient
!! of a cell-based variable at the cell center. The operations performed here
!! combine taking the centered-difference gradient (like in grad_fd_norm) at
!! the 9 edges entering into the usual 9-point RBF vector reconstruction at
!! cell centers, followed by applying the RBF reconstruction.
!!
SUBROUTINE rbf_compute_coeff_c2grad (ptr_patch, ptr_int)

!
!  patch on which computation is performed
!
TYPE(t_patch), INTENT(inout) :: ptr_patch

! interpolation state
!
TYPE(t_int_state), INTENT(inout) :: ptr_int

INTEGER :: jc, jb, je, jcc
INTEGER :: i_rcstartlev
INTEGER :: i_startblk, nblks_c, i_startidx, i_endidx

INTEGER :: ile, ibe, ilc1, ibc1, ilc2, ibc2, ilcc, ibcc

REAL(wp), DIMENSION(nproma,rbf_c2grad_dim,2) :: aux_coeff

!--------------------------------------------------------------------

  i_rcstartlev = 2

  ! Values for the blocking
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)
  nblks_c    = ptr_patch%nblks_c

  ! loop through all patch cells (and blocks)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jcc,jc,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilcc,ibcc,aux_coeff) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

    aux_coeff(:,:,:) = 0._wp

    DO je = 1, rbf_vec_dim_c
      DO jcc = 1, rbf_c2grad_dim
!$NEC ivdep
        DO jc = i_startidx, i_endidx

          IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

          ile = ptr_int%rbf_vec_idx_c(je,jc,jb)
          ibe = ptr_int%rbf_vec_blk_c(je,jc,jb)

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          ilcc = ptr_int%rbf_c2grad_idx(jcc,jc,jb)
          ibcc = ptr_int%rbf_c2grad_blk(jcc,jc,jb)

          IF (ilcc == ilc1 .AND. ibcc == ibc1) THEN

            aux_coeff(jc,jcc,1) = aux_coeff(jc,jcc,1) - &
              ptr_int%rbf_vec_coeff_c(je,1,jc,jb) /    &
              ptr_patch%edges%dual_edge_length(ile,ibe)

            aux_coeff(jc,jcc,2) = aux_coeff(jc,jcc,2) - &
              ptr_int%rbf_vec_coeff_c(je,2,jc,jb) /    &
              ptr_patch%edges%dual_edge_length(ile,ibe)

          ELSE IF (ilcc == ilc2 .AND. ibcc == ibc2) THEN

            aux_coeff(jc,jcc,1) = aux_coeff(jc,jcc,1) + &
              ptr_int%rbf_vec_coeff_c(je,1,jc,jb) /    &
              ptr_patch%edges%dual_edge_length(ile,ibe)

            aux_coeff(jc,jcc,2) = aux_coeff(jc,jcc,2) + &
              ptr_int%rbf_vec_coeff_c(je,2,jc,jb) /    &
              ptr_patch%edges%dual_edge_length(ile,ibe)

          ENDIF

        ENDDO !cell loop
      ENDDO
    ENDDO

!CDIR NOLOOPCHG
    DO jcc = 1, rbf_c2grad_dim
      DO jc = i_startidx, i_endidx
        ptr_int%rbf_c2grad_coeff(jcc,1,jc,jb) = aux_coeff(jc,jcc,1)
        ptr_int%rbf_c2grad_coeff(jcc,2,jc,jb) = aux_coeff(jc,jcc,2)
      ENDDO
    ENDDO

  END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  DO jcc = 1, rbf_c2grad_dim
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_int%rbf_c2grad_coeff(jcc,1,:,:))
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_int%rbf_c2grad_coeff(jcc,2,:,:))
  ENDDO

END SUBROUTINE rbf_compute_coeff_c2grad

!! This routine computes the coefficients needed for vector RBF interpolation to.
!!
!! This routine computes the coefficients needed for vector RBF interpolation to
!! the vertices, which are then stored in the arrays <i>rbf_vec_coeff_v</i>.
!! This computation involves the inversion of the interpolation matrix, which
!! is performed by a Cholesky decomposition.
!! The Cholesky decomposition is currently implemented by a home made routine
!! which can be substituted by a call to a numerical library, if available.
!!
SUBROUTINE rbf_vec_compute_coeff_vertex( ptr_patch, ptr_int )
!

!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int

REAL(wp) :: cc_e1(3), cc_e2(3), cc_v(nproma,3), cc_aux(3) ! coordinates of edge midpoints

REAL(wp)           :: z_lon, z_lat          ! longitude and latitude

! 3d  normal velocity vectors at edge midpoints
REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2, z_nx3

REAL(wp)           :: z_norm                ! norm of velocity vectors

REAL(wp)           :: z_dist                ! distance between data points

REAL(wp)           :: z_nxprod              ! scalar product of normal
                                            ! velocity vectors

REAL(wp),ALLOCATABLE :: z_rbfmat(:,:,:)    ! RBF interpolation matrix

REAL(wp),ALLOCATABLE :: z_diag(:,:)        ! diagonal of cholesky
                                           ! decomposition matrix

REAL(wp),ALLOCATABLE :: z_rbfval(:,:)      ! RBF function value

REAL(wp),ALLOCATABLE :: z_rhs1(:,:), &     ! right hand side of linear
                      & z_rhs2(:,:)        ! interpolation system in 2d

INTEGER :: nblks_v
INTEGER :: istencil(nproma)          ! number of edges for the stencil
INTEGER :: jb                        ! index of current block
INTEGER :: jv                        ! index of current vertex
INTEGER :: je1, je2                  ! edge indices (loop)
INTEGER :: ile1, ibe1, ile2, ibe2    ! edge indices
INTEGER :: ist                       ! status variable
INTEGER :: i_startblk                ! start block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: i_rcstartlev              ! refinement control start level
INTEGER :: i_outunit                 ! base unit for optional debug output
INTEGER :: jg

REAL(wp) ::  checksum_u,checksum_v   ! to check if sum of interpolation coefficients is correct

TYPE(t_cartesian_coordinates), DIMENSION(:,:),   POINTER :: ptr_orient ! pointer to orientation vectors
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff  ! pointer to output coefficients

!--------------------------------------------------------------------

  CALL message(modname//':rbf_vec_compute_coeff_vertex', '')

  jg = ptr_patch%id

  i_rcstartlev = 2

  ! values for the blocking
  nblks_v  = ptr_patch%nblks_v

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%verts%start_blk(i_rcstartlev,1)

  ! orientation vectors of input components
  ptr_orient => ptr_patch%edges%primal_cart_normal

  ! RBF coefficients (output)
  ptr_coeff  => ptr_int%rbf_vec_coeff_v

  ! base unit for optional debug output
  i_outunit  = 520


!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval,z_rhs1,z_rhs2, ist)
  ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_v,rbf_vec_dim_v),  &
            z_diag(nproma,rbf_vec_dim_v),                  &
            z_rbfval(nproma,rbf_vec_dim_v),                &
            z_rhs1(nproma,rbf_vec_dim_v),                  &
            z_rhs2(nproma,rbf_vec_dim_v), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (modname//':rbf_compute_coeff_vertex',  &
      &             'allocation for working arrays failed')
  ENDIF

!$OMP DO PRIVATE (jb,jv,i_startidx,i_endidx,je1,je2,istencil,ist,ile1,ibe1, &
!$OMP             cc_e1,z_lon,z_lat,z_norm,z_nx1,ile2,ibe2,cc_e2,cc_v,      &
!$OMP             cc_aux,z_nx2,z_nxprod,z_dist,z_nx3,checksum_u,&
!$OMP checksum_v) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each vertex, build the vector RBF interpolation matrix
    !
    DO je1 = 1, rbf_vec_dim_v

      DO je2 = 1, je1
        IF (ptr_patch%geometry_info%geometry_type == sphere_geometry) THEN ! use vectorizable version
!$NEC ivdep
          DO jv = i_startidx, i_endidx

            IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) THEN
              ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
              istencil(jv) = 0
              CYCLE
            ENDIF

            ! Get actual number of stencil points
            istencil(jv) = ptr_int%rbf_vec_stencil_v(jv,jb)
            !
            IF ( (je1 > istencil(jv)) .OR. (je2 > istencil(jv)) ) CYCLE
            !
            ! line and block indices for each edge je1 and je2 of RBF stencil
            !
            ile1 = ptr_int%rbf_vec_idx_v(je1,jv,jb)
            ibe1 = ptr_int%rbf_vec_blk_v(je1,jv,jb)
            ile2 = ptr_int%rbf_vec_idx_v(je2,jv,jb)
            ibe2 = ptr_int%rbf_vec_blk_v(je2,jv,jb)
            !
            ! get Cartesian coordinates and orientation vectors
            !
            cc_e1(:) = ptr_patch%edges%cartesian_center(ile1,ibe1)%x(:)
            cc_e2(:) = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
            !
            z_nx1(jv,:) = ptr_orient(ile1,ibe1)%x(:)
            z_nx2(jv,:) = ptr_orient(ile2,ibe2)%x(:)
            !
            ! compute dot product of normal vectors and distance between edge midpoints
            !
            z_nxprod = DOT_PRODUCT(z_nx1(jv,:),z_nx2(jv,:))
            z_dist   = arc_length_v(cc_e1,cc_e2)
            !
            ! set up interpolation matrix
            !
            ! up to now only Gaussian or Inverse multiquadric kernel
            ! (without polynomial correction) are implemented
            IF      (rbf_vec_kern_v == 1) THEN
              z_rbfmat(jv,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_v(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_v == 3) THEN
              z_rbfmat(jv,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_v(MAX(jg,1)))
            ENDIF

            IF (je1 > je2) z_rbfmat(jv,je2,je1) = z_rbfmat(jv,je1,je2)

          END DO
        ELSE ! use generic, unvectorizable version of arc_length_v
          DO jv = i_startidx, i_endidx

            IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) THEN
              ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
              istencil(jv) = 0
              CYCLE
            ENDIF

            ! Get actual number of stencil points
            istencil(jv) = ptr_int%rbf_vec_stencil_v(jv,jb)
            !
            IF ( (je1 > istencil(jv)) .OR. (je2 > istencil(jv)) ) CYCLE
            !
            ! line and block indices for each edge je1 and je2 of RBF stencil
            !
            ile1 = ptr_int%rbf_vec_idx_v(je1,jv,jb)
            ibe1 = ptr_int%rbf_vec_blk_v(je1,jv,jb)
            ile2 = ptr_int%rbf_vec_idx_v(je2,jv,jb)
            ibe2 = ptr_int%rbf_vec_blk_v(je2,jv,jb)
            !
            ! get Cartesian coordinates and orientation vectors
            !
            cc_e1(:) = ptr_patch%edges%cartesian_center(ile1,ibe1)%x(:)
            cc_e2(:) = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
            !
            z_nx1(jv,:) = ptr_orient(ile1,ibe1)%x(:)
            z_nx2(jv,:) = ptr_orient(ile2,ibe2)%x(:)
            !
            ! compute dot product of normal vectors and distance between edge midpoints
            !
            z_nxprod = DOT_PRODUCT(z_nx1(jv,:),z_nx2(jv,:))
            z_dist   = arc_length_v(cc_e1,cc_e2,ptr_patch%geometry_info)
            !
            ! set up interpolation matrix
            !
            ! up to now only Gaussian or Inverse multiquadric kernel
            ! (without polynomial correction) are implemented
            IF      (rbf_vec_kern_v == 1) THEN
              z_rbfmat(jv,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_v(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_v == 3) THEN
              z_rbfmat(jv,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_v(MAX(jg,1)))
            ENDIF

            IF (je1 > je2) z_rbfmat(jv,je2,je1) = z_rbfmat(jv,je1,je2)

          END DO
        ENDIF
      END DO

    END DO

    ! apply Cholesky decomposition to matrix
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL choldec_v(i_startidx,i_endidx,istencil,rbf_vec_dim_v,z_rbfmat,z_diag)
#else
    CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
#endif
!$NEC ivdep
    DO jv = i_startidx, i_endidx

      IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

      !
      ! Solve immediately for coefficients
      !
      ! convert coordinates of vertex to cartesian vector
      !
      cc_v(jv,1:3) = ptr_patch%verts%cartesian(jv,jb)%x(1:3)

      z_lon = ptr_patch%verts%vertex(jv,jb)%lon
      z_lat = ptr_patch%verts%vertex(jv,jb)%lat

      ! Zonal wind component
      CALL gvec2cvec(1._wp,0._wp,z_lon,z_lat,z_nx1(jv,1),z_nx1(jv,2),z_nx1(jv,3),ptr_patch%geometry_info)

      z_norm = SQRT( DOT_PRODUCT(z_nx1(jv,:),z_nx1(jv,:)) )
      z_nx1(jv,:)  = 1._wp/z_norm * z_nx1(jv,:)

      ! Meridional wind component
      CALL gvec2cvec(0._wp,1._wp,z_lon,z_lat,z_nx2(jv,1),z_nx2(jv,2),z_nx2(jv,3),ptr_patch%geometry_info)

      z_norm = SQRT( DOT_PRODUCT(z_nx2(jv,:),z_nx2(jv,:)) )
      z_nx2(jv,:)  = 1._wp/z_norm * z_nx2(jv,:)

    END DO

    !
    ! set up right hand side for interpolation system
    !
    DO je2 = 1, rbf_vec_dim_v
      IF (ptr_patch%geometry_info%geometry_type == sphere_geometry) THEN ! use vectorizable version
!$NEC ivdep
        DO jv = i_startidx, i_endidx

          IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

          IF (je2 > istencil(jv)) CYCLE

          !
          ! get indices and coordinates of edge midpoints and compute distance
          ! to vertex
          !
          ile2   = ptr_int%rbf_vec_idx_v(je2,jv,jb)
          ibe2   = ptr_int%rbf_vec_blk_v(je2,jv,jb)
          !
          cc_e2(:)  = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
          cc_aux(1:3) = cc_v(jv,1:3)
          z_dist = arc_length_v(cc_aux, cc_e2)
          !
          ! get Cartesian orientation vector
          z_nx3(jv,:) = ptr_orient(ile2,ibe2)%x(:)

          ! up to now only Gaussian or Inverse multiquadric kernel
          ! (without polynomial correction) are implemented
          IF      (rbf_vec_kern_v == 1) THEN
            z_rbfval(jv,je2) = gaussi(z_dist,rbf_vec_scale_v(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_v == 3) THEN
            z_rbfval(jv,je2) = inv_multiq(z_dist,rbf_vec_scale_v(MAX(jg,1)))
          ENDIF
          !
          ! compute projection on target vector orientation
          !
          z_rhs1(jv,je2) = z_rbfval(jv,je2) * DOT_PRODUCT(z_nx1(jv,:),z_nx3(jv,:))
          z_rhs2(jv,je2) = z_rbfval(jv,je2) * DOT_PRODUCT(z_nx2(jv,:),z_nx3(jv,:))

        END DO
      ELSE  ! use generic, unvectorizable version of arc_length_v
        DO jv = i_startidx, i_endidx

          IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

          IF (je2 > istencil(jv)) CYCLE

          !
          ! get indices and coordinates of edge midpoints and compute distance
          ! to vertex
          !
          ile2   = ptr_int%rbf_vec_idx_v(je2,jv,jb)
          ibe2   = ptr_int%rbf_vec_blk_v(je2,jv,jb)
          !
          cc_e2(:)  = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
          z_dist = arc_length_v(cc_v(jv,:), cc_e2, ptr_patch%geometry_info)
          !
          ! get Cartesian orientation vector
          z_nx3(jv,:) = ptr_orient(ile2,ibe2)%x(:)

          ! up to now only Gaussian or Inverse multiquadric kernel
          ! (without polynomial correction) are implemented
          IF      (rbf_vec_kern_v == 1) THEN
            z_rbfval(jv,je2) = gaussi(z_dist,rbf_vec_scale_v(MAX(jg,1)))
          ELSE IF (rbf_vec_kern_v == 3) THEN
            z_rbfval(jv,je2) = inv_multiq(z_dist,rbf_vec_scale_v(MAX(jg,1)))
          ENDIF
          !
          ! compute projection on target vector orientation
          !
          z_rhs1(jv,je2) = z_rbfval(jv,je2) * DOT_PRODUCT(z_nx1(jv,:),z_nx3(jv,:))
          z_rhs2(jv,je2) = z_rbfval(jv,je2) * DOT_PRODUCT(z_nx2(jv,:),z_nx3(jv,:))

        END DO
      ENDIF
    END DO

    !
    ! compute vector coefficients
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_v, z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_coeff(:,1,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs1, ptr_coeff(:,1,:,jb))
#endif
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_v, z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_coeff(:,2,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat,  &
      &               z_diag, z_rhs2, ptr_coeff(:,2,:,jb))
#endif

    DO jv = i_startidx, i_endidx

      IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

      ! Ensure that sum of interpolation coefficients is correct

      checksum_u = 0._wp
      checksum_v = 0._wp

      DO je1 = 1, istencil(jv)
        ile1   = ptr_int%rbf_vec_idx_v(je1,jv,jb)
        ibe1   = ptr_int%rbf_vec_blk_v(je1,jv,jb)

        z_nx3(jv,:) = ptr_orient(ile1,ibe1)%x(:)

        checksum_u = checksum_u + ptr_coeff(je1,1,jv,jb)* &
          DOT_PRODUCT(z_nx1(jv,:),z_nx3(jv,:))
        checksum_v = checksum_v + ptr_coeff(je1,2,jv,jb)* &
          DOT_PRODUCT(z_nx2(jv,:),z_nx3(jv,:))
      ENDDO

      DO je1 = 1, istencil(jv)
        ptr_coeff(je1,1,jv,jb) = ptr_coeff(je1,1,jv,jb) / checksum_u
        ptr_coeff(je1,2,jv,jb) = ptr_coeff(je1,2,jv,jb) / checksum_v
      ENDDO

    END DO

  END DO
!$OMP END DO NOWAIT

  DEALLOCATE( z_rbfmat, z_diag, z_rbfval, z_rhs1, z_rhs2,  &
    STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (modname//':rbf_vec_compute_coeff_vertex',  &
      &             'deallocation for working arrays failed')
  ENDIF
!$OMP END PARALLEL

  DO jb = 1, rbf_vec_dim_v
    CALL sync_patch_array(SYNC_V, ptr_patch, ptr_int%rbf_vec_coeff_v(jb,1,:,:))
    CALL sync_patch_array(SYNC_V, ptr_patch, ptr_int%rbf_vec_coeff_v(jb,2,:,:))
  ENDDO

! Optional debug output for RBF coefficients
#ifdef DEBUG_COEFF
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, i_rcstartlev)

    DO jv = i_startidx, i_endidx

      istencil = ptr_int%rbf_vec_stencil_v(jv,jb)
      WRITE(i_outunit+jg,'(2i5,25f12.6)') jb,jv,ptr_coeff(1:istencil(jv),1,jv,jb)
      WRITE(i_outunit+jg,'(2i5,25f12.6)') jb,jv,ptr_coeff(1:istencil(jv),2,jv,jb)
    END DO
  END DO
  CLOSE (i_outunit+jg)
#endif

END SUBROUTINE rbf_vec_compute_coeff_vertex

!-------------------------------------------------------------------------
!
!! This routine computes the coefficients needed for vector RBF interpolation to.
!!
!! This routine computes the coefficients needed for vector RBF interpolation to
!! the edges, which are then stored in the array <i>rbf_vec_coeff_e</i>.
!! This computation involves the inversion of the interpolation matrix, which
!! is performed by a Cholesky decomposition.
!! The Cholesky decomposition is currently implemented by a home made routine
!! which can be substituted by a call to a numerical library, if available.
!!
SUBROUTINE rbf_vec_compute_coeff_edge( ptr_patch, ptr_int )
!

!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int

REAL(wp) :: cc_e1(3), cc_e2(3), cc_e(nproma,3), cc_aux(3) ! coordinates of edge midpoints

REAL(wp)           :: z_lon, z_lat          ! longitude and latitude
REAL(wp)           :: z_nu, z_nv            ! zonal and meridional component
                                            ! of normal velocity vectors
                                            ! at edge midpoints

! 3d  normal velocity vectors at edge midpoints
REAL(wp), DIMENSION (nproma,3) :: z_nx1, z_nx2

REAL(wp)           :: z_norm                ! norm of velocity vectors

REAL(wp)           :: z_dist                ! distance between data points

REAL(wp)           :: z_nxprod              ! scalar product of normal
                                            ! velocity vectors

REAL(wp), ALLOCATABLE :: z_rbfmat(:,:,:)    ! RBF interpolation matrix

REAL(wp), ALLOCATABLE :: z_diag(:,:)        ! diagonal of cholesky
                                            ! decomposition matrix

REAL(wp), ALLOCATABLE :: z_rbfval(:,:)      ! RBF function value


INTEGER :: nblks_e
INTEGER :: istencil(nproma)          ! number of edges for the stencil
INTEGER :: jb                        ! index of current block
INTEGER :: je                        ! index of current edge
INTEGER :: je1, je2                  ! edge indices (loop)
INTEGER :: ile1, ibe1, ile2, ibe2    ! edge indices
INTEGER :: ist                       ! status variable
INTEGER :: i_startblk                ! start block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: i_rcstartlev              ! refinement control start level
INTEGER :: i_outunit                 ! base unit for optional debug output
INTEGER :: jg

REAL(wp) ::  checksum_vt             ! to check if sum of interpolation coefficients is correct

TYPE(t_cartesian_coordinates), DIMENSION(:,:), POINTER :: ptr_orient ! pointer to orientation vectors
REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff  ! pointer to output coefficients

! pointer to orientation of output vectors
TYPE(t_tangent_vectors), DIMENSION(:,:), POINTER :: ptr_orient_out

!--------------------------------------------------------------------

  CALL message(modname//':rbf_vec_compute_coeff_edge', '')

  jg = ptr_patch%id


  i_rcstartlev = 2

  ! values for the blocking
  nblks_e  = ptr_patch%nblks_e

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%edges%start_blk(i_rcstartlev,1)

  ! orientation vectors of input components
  ptr_orient => ptr_patch%edges%primal_cart_normal

  ! orientation vectors of output components
  ptr_orient_out => ptr_patch%edges%dual_normal

  ! RBF coefficients (output)
  ptr_coeff  => ptr_int%rbf_vec_coeff_e

  ! base unit for optional debug output
  i_outunit  = 530

!$OMP PARALLEL PRIVATE (z_rbfmat,z_diag,z_rbfval, ist)
  ALLOCATE( z_rbfmat(nproma,rbf_vec_dim_e,rbf_vec_dim_e),  &
            z_diag(nproma,rbf_vec_dim_e),                  &
            z_rbfval(nproma,rbf_vec_dim_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (modname//':rbf_vec_compute_coeff_edge_t',  &
      &             'allocation for working arrays failed')
  ENDIF

!$OMP DO PRIVATE (jb,je,i_startidx,i_endidx,je1,je2,istencil,        &
!$OMP    ist,ile1,ibe1,cc_e1,z_nu,z_nv,z_lon,z_lat,z_norm,z_nx1,     &
!$OMP    ile2,ibe2,cc_e2,cc_e,cc_aux,z_nx2,z_nxprod,z_dist,&
!$OMP checksum_vt) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each edge, build the vector RBF interpolation matrix
    !
    DO je1 = 1, rbf_vec_dim_e

      DO je2 = 1, je1
        IF (ptr_patch%geometry_info%geometry_type == sphere_geometry) THEN ! use vectorizable version
!$NEC ivdep
          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) THEN
              ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
              istencil(je) = 0
              CYCLE
            ENDIF

            istencil(je) = ptr_int%rbf_vec_stencil_e(je,jb)
            !
            IF ( (je1 > istencil(je)) .OR. (je2 > istencil(je)) ) CYCLE
            !
            ! line and block indices for each edge je1 and je2 of RBF stencil
            !
            ile1 = ptr_int%rbf_vec_idx_e(je1,je,jb)
            ibe1 = ptr_int%rbf_vec_blk_e(je1,je,jb)
            ile2 = ptr_int%rbf_vec_idx_e(je2,je,jb)
            ibe2 = ptr_int%rbf_vec_blk_e(je2,je,jb)
            !
            ! get Cartesian coordinates and orientation vectors
            !
            cc_e1(:) = ptr_patch%edges%cartesian_center(ile1,ibe1)%x(:)
            cc_e2(:) = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
            !
            z_nx1(je,:) = ptr_orient(ile1,ibe1)%x(:)
            z_nx2(je,:) = ptr_orient(ile2,ibe2)%x(:)
            !
            ! compute dot product of normal vectors and distance between edge midpoints
            !
            z_nxprod = DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))
            z_dist   = arc_length_v(cc_e1,cc_e2)

            ! set up interpolation matrix
            !
            ! up to now only Gaussian or Inverse multiquadric kernel
            ! (without polynomial correction) are implemented
            IF      (rbf_vec_kern_e == 1) THEN
              z_rbfmat(je,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_e(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_e == 3) THEN
              z_rbfmat(je,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_e(MAX(jg,1)))
            ENDIF

            IF (je1 > je2) z_rbfmat(je,je2,je1) = z_rbfmat(je,je1,je2)

          END DO
        ELSE ! use generic, unvectorizable version of arc_length_v
          DO je = i_startidx, i_endidx

            IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) THEN
              ! Avoid the matrix decomposition for boundary cells since the matrix might get singular
              istencil(je) = 0
              CYCLE
            ENDIF

            istencil(je) = ptr_int%rbf_vec_stencil_e(je,jb)
            !
            IF ( (je1 > istencil(je)) .OR. (je2 > istencil(je)) ) CYCLE
            !
            ! line and block indices for each edge je1 and je2 of RBF stencil
            !
            ile1 = ptr_int%rbf_vec_idx_e(je1,je,jb)
            ibe1 = ptr_int%rbf_vec_blk_e(je1,je,jb)
            ile2 = ptr_int%rbf_vec_idx_e(je2,je,jb)
            ibe2 = ptr_int%rbf_vec_blk_e(je2,je,jb)
            !
            ! get Cartesian coordinates and orientation vectors
            !
            cc_e1(:) = ptr_patch%edges%cartesian_center(ile1,ibe1)%x(:)
            cc_e2(:) = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
            !
            z_nx1(je,:) = ptr_orient(ile1,ibe1)%x(:)
            z_nx2(je,:) = ptr_orient(ile2,ibe2)%x(:)
            !
            ! compute dot product of normal vectors and distance between edge midpoints
            !
            z_nxprod = DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))
            z_dist   = arc_length_v(cc_e1,cc_e2,ptr_patch%geometry_info)

            ! set up interpolation matrix
            !
            ! up to now only Gaussian or Inverse multiquadric kernel
            ! (without polynomial correction) are implemented
            IF      (rbf_vec_kern_e == 1) THEN
              z_rbfmat(je,je1,je2) = z_nxprod * gaussi(z_dist,rbf_vec_scale_e(MAX(jg,1)))
            ELSE IF (rbf_vec_kern_e == 3) THEN
              z_rbfmat(je,je1,je2) = z_nxprod * inv_multiq(z_dist,rbf_vec_scale_e(MAX(jg,1)))
            ENDIF

            IF (je1 > je2) z_rbfmat(je,je2,je1) = z_rbfmat(je,je1,je2)

          END DO
        ENDIF
      END DO

    END DO


    ! apply Cholesky decomposition to matrix
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL choldec_v(i_startidx,i_endidx,istencil,rbf_vec_dim_e,z_rbfmat,z_diag)
#else
    CALL choldec_v(i_startidx,i_endidx,istencil,              z_rbfmat,z_diag)
#endif
!$NEC ivdep
    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

      !
      ! Solve immediately for coefficients
      !
      ! convert coordinates of edge midpoint to cartesian vector
      !
      cc_e(je,1:3) = ptr_patch%edges%cartesian_center(je,jb)%x(1:3)

      z_nu  = ptr_orient_out(je,jb)%v1
      z_nv  = ptr_orient_out(je,jb)%v2
      z_lon = ptr_patch%edges%center(je,jb)%lon
      z_lat = ptr_patch%edges%center(je,jb)%lat

      CALL gvec2cvec(z_nu,z_nv,z_lon,z_lat,z_nx1(je,1),z_nx1(je,2),z_nx1(je,3),ptr_patch%geometry_info)

      z_norm = SQRT( DOT_PRODUCT(z_nx1(je,:),z_nx1(je,:)) )
      z_nx1(je,:)  = 1._wp/z_norm * z_nx1(je,:)

    END DO

    !
    ! set up right hand side for interpolation system
    !
    DO je2 = 1, rbf_vec_dim_e
      IF (ptr_patch%geometry_info%geometry_type == sphere_geometry) THEN ! use vectorizable version
!$NEC ivdep
        DO je = i_startidx, i_endidx

          IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

          IF (je2 > istencil(je)) CYCLE
          !
          ! get indices and coordinates of edge midpoints and compute distance
          ! to the edge where the vector is reconstructed
          !
          ile2   = ptr_int%rbf_vec_idx_e(je2,je,jb)
          ibe2   = ptr_int%rbf_vec_blk_e(je2,je,jb)

          cc_e2(:)  = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
          cc_aux(1:3) = cc_e(je,1:3)
          z_dist = arc_length_v(cc_aux, cc_e2)
          !
          ! get Cartesian orientation vector
          z_nx2(je,:) = ptr_orient(ile2,ibe2)%x(:)

          z_nxprod = DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))

          ! up to now only Gaussian or Inverse multiquadric kernel
          ! (without polynomial correction) are implemented
          IF      (rbf_vec_kern_e == 1) THEN
            z_rbfval(je,je2) = gaussi(z_dist,rbf_vec_scale_e(MAX(jg,1))) * z_nxprod
          ELSE IF (rbf_vec_kern_e == 3) THEN
            z_rbfval(je,je2) = inv_multiq(z_dist,rbf_vec_scale_e(MAX(jg,1))) * z_nxprod
          ENDIF

        END DO
      ELSE ! use generic, unvectorizable version of arc_length_v
        DO je = i_startidx, i_endidx

          IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

          IF (je2 > istencil(je)) CYCLE
          !
          ! get indices and coordinates of edge midpoints and compute distance
          ! to the edge where the vector is reconstructed
          !
          ile2   = ptr_int%rbf_vec_idx_e(je2,je,jb)
          ibe2   = ptr_int%rbf_vec_blk_e(je2,je,jb)

          cc_e2(:)  = ptr_patch%edges%cartesian_center(ile2,ibe2)%x(:)
          z_dist = arc_length_v(cc_e(je,:), cc_e2, ptr_patch%geometry_info)
          !
          ! get Cartesian orientation vector
          z_nx2(je,:) = ptr_orient(ile2,ibe2)%x(:)

          z_nxprod = DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))

          ! up to now only Gaussian or Inverse multiquadric kernel
          ! (without polynomial correction) are implemented
          IF      (rbf_vec_kern_e == 1) THEN
            z_rbfval(je,je2) = gaussi(z_dist,rbf_vec_scale_e(MAX(jg,1))) * z_nxprod
          ELSE IF (rbf_vec_kern_e == 3) THEN
            z_rbfval(je,je2) = inv_multiq(z_dist,rbf_vec_scale_e(MAX(jg,1))) * z_nxprod
          ENDIF

        END DO
      ENDIF
    END DO

    !
    ! compute vector coefficients
    !
!CDIR NOIEXPAND
#ifdef __SX__
    CALL solve_chol_v(i_startidx, i_endidx, istencil, rbf_vec_dim_e, z_rbfmat, &
                    z_diag, z_rbfval, ptr_coeff(:,:,jb))
#else
    CALL solve_chol_v(i_startidx, i_endidx, istencil,                z_rbfmat, &
                    z_diag, z_rbfval, ptr_coeff(:,:,jb))
#endif

    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

      ! Ensure that sum of interpolation coefficients is correct

      checksum_vt = 0._wp

      DO je1 = 1, istencil(je)
        ile1   = ptr_int%rbf_vec_idx_e(je1,je,jb)
        ibe1   = ptr_int%rbf_vec_blk_e(je1,je,jb)

        z_nx2(je,:) = ptr_orient(ile1,ibe1)%x(:)

        checksum_vt = checksum_vt+ptr_coeff(je1,je,jb)*&
          DOT_PRODUCT(z_nx1(je,:),z_nx2(je,:))
      ENDDO

      DO je1 = 1, istencil(je)
        ptr_coeff(je1,je,jb) = ptr_coeff(je1,je,jb) / checksum_vt
      ENDDO

    END DO

  END DO
!$OMP END DO NOWAIT

  DEALLOCATE( z_rbfmat, z_diag, z_rbfval, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish (modname//':rbf_vec_compute_coeff_edge',  &
      &             'deallocation for working arrays failed')
  ENDIF
!$OMP END PARALLEL

  DO jb = 1, rbf_vec_dim_e
    CALL sync_patch_array(SYNC_E, ptr_patch, ptr_int%rbf_vec_coeff_e(jb,:,:))
  ENDDO

! Optional debug output for RBF coefficients
#ifdef DEBUG_COEFF
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, i_rcstartlev)

    DO je = i_startidx, i_endidx

      istencil(je) = ptr_int%rbf_vec_stencil_e(je,jb)
      WRITE(i_outunit+jg,'(2i5,25f12.6)') jb,je,ptr_coeff(1:istencil(je),je,jb)
    END DO
  END DO
  CLOSE (i_outunit+jg)
#endif

END SUBROUTINE rbf_vec_compute_coeff_edge


END MODULE mo_intp_rbf_coeffs
