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

! Module for efficient promixity queries / nearest neighbor search
! with "geometric near-neighbor access trees" (GNATs), i.e.
! a tree structure based on generalized hyperplanes.
!
! This algorithm is divided into two parts: In the first phase we
! build a search tree containing the cell centers, which is then in a
! second phase traversed while searching for cells in the vicinity of
! a given search point.
!
! Both the tree construction and the point query can be performed
! with in multiple threads.
!
! See
! Brin, Sergey: "Near Neighbor Search in Large Metric Spaces"
! VLDB '95 : Proceedings of the 21st International Conference on
!            Very Large Data Bases,
! Zurich Switzerland, Sept. 11--15, 1995, pp. 574-584,
! Morgan Kaufmann Publishers, 1995.
!
! Basic usage example:
!
!    ! data structure
!    type(t_gnat) :: gnat
!
!    ! build fast search tree for proximity queries in unstructured triangular grid
!    ! (at the beginning of the program)
!    CALL gnat_init_grid(gnat, p_patch, .TRUE., 1, nblks_c)
!
!    ! perform a search of a single point
!    min_node_idx(1:2) = 0
!    vmin_dist         = MAX_RANGE
!    v                 = (/ plam * pi_180, pphi * pi_180 /)        ! search point
!    r                 = gnat_std_radius(gnat)                     ! search radius
!    CALL gnat_query(gnat, v, r, vmin_dist, min_node_idx)
!    jc = min_node_idx(1)
!    jb = min_node_idx(2)
!
!    ! finish (at the end of the program)
!    CALL gnat_destroy(grid%gnat)
!
! @todo On different platforms: test for beneficiary effects of NOT using REAL(wp)
!       for distance measurements.
!
! @note Load balancing issues:
!
!       When searching for the cells containing the points of a lon-lat
!       grid the meridian convergence affects this second phase only: The
!       MPI process which covers the pole region must handle a larger
!       number of search operations.

MODULE mo_gnat_gridsearch

#ifdef _OPENMP
  USE OMP_LIB
#endif

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_math_constants,      ONLY: pi_180
  USE mo_math_types,          ONLY: t_geographical_coordinates
  USE mo_model_domain,        ONLY: t_grid_cells, t_grid_vertices, t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_mpi,                 ONLY: get_my_mpi_work_id, p_pe_work, &
    &                               p_comm_work, my_process_is_mpi_test, &
    &                               p_bcast, p_max, p_comm_work_2_test
  USE mo_communication,       ONLY: idx_1d
  USE mo_icon_comm_lib,       ONLY: t_mpi_mintype, mpi_reduce_mindistance_pts
#ifndef NOMPI
  USE mpi
#endif

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_gnat_gridsearch'

  ! ---------------------------------------------------------
  ! constants
  ! ---------------------------------------------------------

  !> level of output verbosity
  INTEGER, PARAMETER  :: dbg_level    = 0

  ! GNAT degree
  !
  ! greater values lead to a time-consuming tree construction process,
  ! while the query time decreases.
  !
  ! On the SX9 architecture, distance computations
  ! are vectorized. Thus it is beneficial to choose
  ! a shallow tree with a large degree.
#ifdef __SX__
  INTEGER, PARAMETER :: gnat_k        = 30
#else
  INTEGER, PARAMETER :: gnat_k        = 10
#endif

  !> REAL kind value for GNAT data structure
  INTEGER, PARAMETER  :: gk           = wp

  !> maximum distance (unit sphere!)
  REAL(gk), PARAMETER :: MAX_RANGE    = 99._gk

  !> marks a free pointer
  INTEGER,  PARAMETER :: UNASSOCIATED = -1

  !> value for invalid result (algorithm failure)
  INTEGER, PARAMETER  :: INVALID_NODE = -1
  INTEGER, PARAMETER  :: UNUSED_NODE  = -2
  INTEGER, PARAMETER  :: SKIP_NODE    = -3

  !> dimension of coordinate system
  INTEGER, PARAMETER  :: icoord_dim   =  2


  ! ---------------------------------------------------------
  ! declaration of derived data types
  ! ---------------------------------------------------------

  !> coordinate type
  ! note: we do not use mo_base_geometry::t_geographical_coordinates
  !       here, since icoord_dim appears as a dimension in many lists
  TYPE  t_coord
    REAL(gk)                    :: p (icoord_dim)
    REAL(gk)                    :: sin_p, cos_p  ! precomputed sin/cos of latitude
    INTEGER                     :: idx(3)        ! dimensions: (idx,blk,glb_idx)
    ! (where global index is required for "equal distance decisions")
  END TYPE t_coord

  !> tree node
  TYPE t_gnat_node
    INTEGER                     :: isplit_pts
    TYPE (t_coord)              :: p(gnat_k)
    INTEGER                     :: child(gnat_k) ! index pointer
    REAL(gk)                    :: drange(gnat_k, gnat_k, 2)
  END TYPE t_gnat_node

  !> tree structure based on generalized hyperplanes
  TYPE t_gnat_tree
    !> GNAT-based tree (index of root node)
    INTEGER :: gnat_tree = UNASSOCIATED
    !> the linear array where the tree nodes are actually stored:
    TYPE (t_gnat_node), ALLOCATABLE :: node_storage(:)
    !> size of array "node_storage"
    INTEGER :: max_num_nodes      = 0
    !> expected total number of nodes (may be altered during init)
    INTEGER :: expected_num_nodes = 100000
    !> current usage of array "node_storage"
    INTEGER :: num_nodes          = 0
  END TYPE t_gnat_tree


  ! ---------------------------------------------------------
  ! public interface definition
  ! ---------------------------------------------------------

  PRIVATE

  ! functions and subroutines
  !
  ! tree initialization and destruction
  PUBLIC :: gnat_init_grid
  PUBLIC :: gnat_destroy
  ! querying single points / lists of points
  PUBLIC :: gnat_recursive_query            !< nearest cell circumcenter, recursive implementation
  PUBLIC :: gnat_recursive_proximity_query  !< test if any cell circumcenter within radius
  PUBLIC :: gnat_query                      !< nearest cell circumcenter, non-recursive implementation
  PUBLIC :: gnat_query_nnb                  !< nearest neighbor search for a list of points
  PUBLIC :: flag_ll_points                  !< divide&conquer algorithm for regular grids
  ! identifying GNAT subtrees
  PUBLIC :: gnat_find_subtrees
  PUBLIC :: gnat_find_subtrees2
  PUBLIC :: gnat_flag_subtree
  ! utility functions
  PUBLIC :: gnat_std_radius
  PUBLIC :: dist_p
  ! distributed point query with MPI
  PUBLIC :: gnat_query_containing_triangles
  PUBLIC :: gnat_merge_distributed_queries
  !
  ! constants and derived data types
  !
  PUBLIC :: gk
  PUBLIC :: gnat_k
  PUBLIC :: icoord_dim
  PUBLIC :: t_gnat_node
  PUBLIC :: t_gnat_tree
  PUBLIC :: t_coord
  PUBLIC :: MAX_RANGE
  PUBLIC :: UNASSOCIATED, INVALID_NODE, SKIP_NODE

CONTAINS

  ! ---------------------------------------------------------
  ! subroutines for tree initialization and destruction
  ! ---------------------------------------------------------


  !> Read ICON triangle center points into GNAT data structure.
  ! (Possible) improvement for the future:
  ! - implementation of min/max range boundaries if not
  !   all grid points are of interest.
  SUBROUTINE gnat_init_grid(gnat, p_patch, opt_ldegree, opt_startblk, opt_endblk)
    TYPE (t_gnat_tree),         INTENT(INOUT) :: gnat
    TYPE(t_patch),              INTENT(IN)    :: p_patch
    LOGICAL, INTENT(IN), OPTIONAL :: opt_ldegree
    INTEGER, INTENT(IN), OPTIONAL :: opt_startblk, opt_endblk
    ! local variables:
    CHARACTER(*), PARAMETER :: routine = modname//"::gnat_init_grid"
    INTEGER                      :: nproc

#ifdef _OPENMP
    nproc = OMP_GET_MAX_THREADS()
#else
    nproc = 1
#endif

    ! consistency check
    IF (gnat%gnat_tree /= UNASSOCIATED) THEN
      CALL message(routine, "Discarding existing data in GNAT structure!")
      CALL gnat_destroy(gnat)
    END IF

    ! tell the GNAT how many nodes are expected
    ! ("5" is for safety reasons)
    gnat%expected_num_nodes = 5*p_patch%n_patch_cells/gnat_k + 1

    ! build GNAT data structure based on clon, clat:
    CALL gnat_insert_mt(gnat, p_patch, nproc, opt_ldegree, opt_startblk, opt_endblk)

  END SUBROUTINE gnat_init_grid


  !>  Destroy GNAT data structure.
  RECURSIVE SUBROUTINE gnat_destroy(gnat)
    TYPE (t_gnat_tree), INTENT(INOUT) :: gnat
    INTEGER :: errstat

    IF (gnat%max_num_nodes > 0) THEN
      DEALLOCATE(gnat%node_storage, STAT=errstat)
      IF (errstat /= 0)  &
          CALL finish ('mo_gnat_gridsearch:gnat_destroy', &
          &            'Error in DEALLOCATE operation!')
      gnat%num_nodes     = 0
      gnat%max_num_nodes = 0
      gnat%gnat_tree     = UNASSOCIATED
    END IF
  END SUBROUTINE gnat_destroy


  ! ---------------------------------------------------------
  ! subroutines for querying single points / lists of points
  ! ---------------------------------------------------------


  !>  Identify points within given radius.
  RECURSIVE SUBROUTINE gnat_recursive_query(gnat, tree_idx, v, r, min_dist, min_node_idx)
    TYPE (t_gnat_tree), INTENT(IN), TARGET :: gnat
    INTEGER,            INTENT(IN)    :: tree_idx ! index pointer
    REAL(gk),           INTENT(IN)    :: v(icoord_dim)
    REAL(gk),           INTENT(IN)    :: r
    REAL(gk),           INTENT(INOUT) :: min_dist
    INTEGER,            INTENT(INOUT) :: min_node_idx(3)
    ! local variables
    LOGICAL                           :: pflag(gnat_k)
    INTEGER                           :: ip, isplit_pts
    REAL(gk)                          :: dist_vp
    TYPE(t_gnat_node),  POINTER       :: tree
#ifdef __SX__
    REAL(gk)                          :: pdist(gnat_k)
#endif

    tree => gnat%node_storage(tree_idx)
    isplit_pts = tree%isplit_pts
    IF (isplit_pts == 0) RETURN

    ! include all split points in set P
    ! note that values for indices > tree%isplit_pts are undefined!
    pflag(1:isplit_pts) = .TRUE.

    ! note: On the SX it might be cheaper to compute all distances (vectorized)
    ! outside of the following loop (though many won't be needed then)
#ifdef __SX__
!CDIR IEXPAND
    CALL dist_vect(tree%p, v, isplit_pts, pdist)
#endif

    ! remove some elements from P
    P : DO ip=1,isplit_pts
      IF (pflag(ip)) THEN

#ifdef __SX__
        dist_vp = pdist(ip)
#else
        dist_vp = dist(tree%p(ip), v)
#endif

        IF (dist_vp <= r) THEN
          IF (dist_vp < min_dist) THEN
            min_dist = dist_vp
            min_node_idx = tree%p(ip)%idx
          ELSE IF (dist_vp == min_dist) THEN
            IF (min_node_idx(3) > tree%p(ip)%idx(3)) THEN
              min_dist = dist_vp
              min_node_idx = tree%p(ip)%idx
            END IF
          END IF
        END IF

        WHERE (pflag(1:isplit_pts))
          pflag(1:isplit_pts) = (tree%drange(1:isplit_pts,ip,1) <= (dist_vp+r)) .AND. &
            &                   (tree%drange(1:isplit_pts,ip,2) >= (dist_vp-r))
        END WHERE

      END IF
    END DO P

    ! traverse subtrees for remaining ip in P
    DO ip=1,isplit_pts
      IF (pflag(ip) .AND. (tree%child(ip) /= UNASSOCIATED)) THEN
        CALL gnat_recursive_query(gnat, tree%child(ip), v, r, min_dist, min_node_idx(:))
      END IF
    END DO
  END SUBROUTINE gnat_recursive_query


  !>  @return .TRUE. if **any** point lies inside a given radius.
  !
  RECURSIVE SUBROUTINE gnat_recursive_proximity_query(gnat, tree_idx, v, r, is_pt_inside)
    TYPE (t_gnat_tree), INTENT(IN), TARGET :: gnat
    INTEGER,            INTENT(IN)    :: tree_idx ! index pointer
    REAL(gk),           INTENT(IN)    :: v(icoord_dim)
    REAL(gk),           INTENT(IN)    :: r
    LOGICAL,            INTENT(INOUT) :: is_pt_inside
    ! local variables
    LOGICAL                           :: pflag(gnat_k)
    INTEGER                           :: ip, isplit_pts
    REAL(gk)                          :: dist_vp
    TYPE(t_gnat_node), POINTER        :: tree
#ifdef __SX__
    REAL(gk)                          :: pdist(gnat_k)
#endif

    tree => gnat%node_storage(tree_idx)
    isplit_pts = tree%isplit_pts
    IF (isplit_pts == 0) RETURN

    ! include all split points in set P
    ! note that values for indices > tree%isplit_pts are undefined!
    pflag(1:isplit_pts) = .TRUE.

    ! note: On the SX it might be cheaper to compute all distances (vectorized)
    ! outside of the following loop (though many won't be needed then)
#ifdef __SX__
!CDIR IEXPAND
    CALL dist_vect(tree%p, v, isplit_pts, pdist)
#endif

    ! remove some elements from P
    P : DO ip=1,isplit_pts
      IF (pflag(ip)) THEN

#ifdef __SX__
        dist_vp = pdist(ip)
#else
        dist_vp = dist(tree%p(ip), v)
#endif

        IF (dist_vp <= r) THEN
          is_pt_inside = .TRUE.
          RETURN
        END IF

        WHERE (pflag(1:isplit_pts))
          pflag(1:isplit_pts) = (tree%drange(1:isplit_pts,ip,1) <= (dist_vp+r)) .AND. &
            &                   (tree%drange(1:isplit_pts,ip,2) >= (dist_vp-r))
        END WHERE

      END IF
    END DO P

    ! traverse subtrees for remaining ip in P
    DO ip=1,isplit_pts
      IF (pflag(ip) .AND. (tree%child(ip) /= UNASSOCIATED)) THEN
        CALL gnat_recursive_proximity_query(gnat, tree%child(ip), v, r, is_pt_inside)
      END IF
      IF (is_pt_inside) RETURN
    END DO
  END SUBROUTINE gnat_recursive_proximity_query


  !> Identify nearest point within given radius, non-recursive, partly
  !  vectorizable version of the algorithm.
  !
  SUBROUTINE gnat_query(gnat, v, r, min_dist, min_node_idx)
    TYPE (t_gnat_tree), INTENT(IN)    :: gnat
    REAL(gk),           INTENT(IN)    :: v(icoord_dim)
    REAL(gk),           INTENT(IN)    :: r
    REAL(gk),           INTENT(INOUT) :: min_dist
    INTEGER,            INTENT(INOUT) :: min_node_idx(3)
    ! local variables
    INTEGER, PARAMETER :: NLIST = 100 ! max size of traversal list
    LOGICAL                           :: pflag(gnat_k,NLIST)
    INTEGER                           :: traversal_list(NLIST),        &
      &                                  ntrv_list(gnat_k*NLIST),      &
      &                                  ip, isplit_pts, ntraversal0,  &
      &                                  i, ti, ntraversal, ihead
    TYPE(t_coord)                     :: p_v, p_ip
#ifdef __SX__
    REAL(gk)                          :: dist_v(gnat_k)
#endif
    REAL(gk)                          :: dist_vp, r0

    ! initialization
    ihead = gnat%gnat_tree
    p_v%p(:)  = v(:)
    p_v%sin_p = SIN(v(2))
    p_v%cos_p = COS(v(2))

    ntraversal        = 1
    traversal_list(1) = ihead
    r0                = MIN(r, min_dist)
    DO
      pflag(:,:) = .TRUE.

      ! remove some elements from P
      DO i=1,ntraversal
        ti         = traversal_list(i)
        isplit_pts = gnat%node_storage(ti)%isplit_pts

#ifdef __SX__
        CALL dist_vect(gnat%node_storage(ti)%p, v, isplit_pts, dist_v)
#endif

        P : DO ip=1,isplit_pts
          IF (pflag(ip,i)) THEN
#ifdef __SX__
            dist_vp = dist_v(ip)
#else
            p_ip    = gnat%node_storage(ti)%p(ip)
            dist_vp = dist_t_coord(p_ip,p_v)
#endif

            IF (dist_vp <= r0) THEN
#ifdef __SX__
              p_ip = gnat%node_storage(ti)%p(ip)
#endif
              IF (dist_vp == r0) THEN
                IF (min_node_idx(3) > p_ip%idx(3)) THEN
                  min_node_idx(:) = p_ip%idx(:)
                  r0              = dist_vp
                END IF
              ELSE
                min_node_idx(:) = p_ip%idx(:)
                r0              = dist_vp
              END IF
            END IF

            WHERE (pflag(1:isplit_pts,i))
              pflag(1:isplit_pts,i) = (gnat%node_storage(ti)%drange(1:isplit_pts,ip,1) <= (dist_vp+r0)) .AND. &
                &                     (gnat%node_storage(ti)%drange(1:isplit_pts,ip,2) >= (dist_vp-r0))
            END WHERE
          END IF
        END DO P
      END DO

      ! store subtrees for remaining ip in P (traversed later):
      ntraversal0 = ntraversal
      ntraversal  = 0
      DO ip=1,gnat_k
        DO i=1,ntraversal0
          IF (pflag(ip,i) .AND. (gnat%node_storage(traversal_list(i))%child(ip) /= UNASSOCIATED)) THEN
            ntraversal = ntraversal + 1
            ntrv_list(ntraversal) = gnat%node_storage(traversal_list(i))%child(ip)
          END IF
        END DO
      END DO

      IF (ntraversal == 0) THEN
        ! abort, if no more nodes to be traversed
        min_dist = r0
        EXIT
      ELSE
        ! copy new traversal list:
        traversal_list(1:ntraversal) = ntrv_list(1:ntraversal)
      END IF
    END DO
  END SUBROUTINE gnat_query


  ! Performs a nearest-neighbor query for a given list of points and
  ! returns the indices and block indices of the mesh triangles that
  ! contain these points.
  !
  SUBROUTINE gnat_query_nnb(gnat, v, iv_nproma, iv_nblks,   &
    &                       iv_npromz, radius, tri_idx, min_dist)

    TYPE (t_gnat_tree), INTENT(IN) :: gnat
    INTEGER,  INTENT(IN)    :: iv_nproma, iv_nblks, iv_npromz      ! list size
    REAL(wp), INTENT(IN)    :: radius
    REAL(gk), INTENT(IN)    :: v(iv_nproma, iv_nblks, icoord_dim)  ! list of search points
    INTEGER,  INTENT(INOUT) :: tri_idx(2,iv_nproma, iv_nblks)      ! containing triangle (idx,block)
    REAL(gk), INTENT(OUT)   :: min_dist(iv_nproma, iv_nblks)       ! minimal distance
    ! local parameters
    CHARACTER(*), PARAMETER :: routine = modname//"::gnat_query_nnb"

    INTEGER                 :: min_node_idx(3,iv_nproma, iv_nblks)  ! corresponding GNAT nodes
    INTEGER                 :: jb, jc, end_idx

    ! set default value ("failure notice")
    min_dist(:,:)  = MAX_RANGE

    ! query list of nearest neighbors
    ! TODO[FP] : For some test cases it might be reasonable to enable
    ! radius adaptation
    min_node_idx(1,:,:) = tri_idx(1,:,:)
    CALL gnat_query_list_mt(gnat, v, iv_nproma, iv_nblks, iv_npromz,     &
      &                     min_dist, min_node_idx, .FALSE., opt_rr=radius)
    WHERE (min_node_idx(1,:,:) == SKIP_NODE)
      min_node_idx(1,:,:) = INVALID_NODE
      min_dist(:,:)       = MAX_RANGE
    END WHERE

    ! loop over blocks
!$OMP PARALLEL DO private(jb, end_idx, jc)
    DO jb=1,iv_nblks

      ! set end index in current block:
      end_idx = MERGE(iv_nproma, iv_npromz, jb /= iv_nblks)

      ! loop over indices
      DO jc=1,end_idx
        tri_idx(1:2,jc,jb) = min_node_idx(1:2,jc,jb)
      END DO
    END DO
!$OMP END PARALLEL DO
  END SUBROUTINE gnat_query_nnb


  !> Query a sequence of points, exploiting previous search results.
  !
  SUBROUTINE gnat_query_list(gnat, v, iv_nproma, iv_nblks, iv_npromz, iv, istart, &
    &                        rr, min_dist, min_node_idx)

    TYPE (t_gnat_tree), INTENT(IN) :: gnat
    INTEGER,       INTENT(IN)    :: iv_nproma, iv_nblks, iv_npromz     ! list size
    INTEGER,       INTENT(IN)    :: iv                   ! local chunk size (blocks)
    INTEGER,       INTENT(IN)    :: istart               ! start idx
    REAL(gk),      INTENT(IN)    :: v(iv_nproma, iv_nblks, icoord_dim)    ! search point
    REAL(gk),      INTENT(IN)    :: rr                   ! initial search radius
    REAL(gk),      INTENT(INOUT) :: min_dist(iv_nproma, iv_nblks)         ! minimal distance
    INTEGER,       INTENT(INOUT) :: min_node_idx(3,iv_nproma, iv_nblks)   ! corresponding node
    ! local parameters
    INTEGER                      :: jb, jc, jc1, end_idx      ! block, index loop counter
    REAL(gk)                     :: r2, r, min_dist_old
    REAL(gk)                     :: p_old(icoord_dim), p_new(icoord_dim)
    REAL(gk)                     :: vmin_dist
    INTEGER                      :: vmin_node_idx(3), indices(iv_nproma), &
      &                             all_indices(iv_nproma), zero(iv_nproma), tree

    tree = gnat%gnat_tree

    ! Description:

    ! Often, proximity search for a sequence of points means that
    ! subsequent points are located next to each other. In this case we can
    ! shorten the tree traversal (and reduce the number of distance
    ! computations) by the following strategies:
    !    the search result from a previous run gives us a true upper
    !    bound for the new search radius that may be smaller than the
    !    previous search radius:
    !    radius_new = min(radius_old, min_dist_old + |v_old, v_new|)

    r                              = rr
    min_dist(:,istart:istart+iv-1) = MAX_RANGE
    min_dist_old                   = 0._gk
    p_old(:)                       = 0._gk
    all_indices(:)                 = (/ (jc, jc=1,iv_nproma) /)
    zero(:)                        = 0

    DO jb=istart, istart+iv-1
      ! set end index in current block:
      end_idx = iv_nproma
      IF (jb == iv_nblks) THEN
        end_idx = iv_npromz
        min_node_idx(1, (end_idx+1):, jb) = UNUSED_NODE
      END IF

      ! Build an index list: Skip search points that have been externally
      ! flagged with index "SKIP_NODE":
      indices(1:end_idx) = PACK(all_indices(1:end_idx), (min_node_idx(1,1:end_idx,jb) /= SKIP_NODE), zero(1:end_idx))

      DO jc1=1,end_idx
        jc = indices(jc1)
        IF (jc == 0)  EXIT

        p_new(1:icoord_dim) = v(jc,jb,1:icoord_dim)

        ! first query, block "istart", index "1": no previous result
        ! available, otherwise use the previous search result
        IF ((jb > istart) .OR. (jc>1)) THEN

          ! distance between old and new search point:
          r2 = dist_p(p_old(:), p_new(:))

          ! compute a true upper bound for search radius
          ! ("1.05" is just for safety)
          r = MIN(rr, 1.05_gk*(r2 + min_dist_old))

        END IF

        ! traverse tree
        vmin_dist        = MAX_RANGE
        vmin_node_idx(1) = INVALID_NODE

        CALL gnat_recursive_query(gnat, tree, p_new(:), r, vmin_dist, vmin_node_idx(:))
        min_dist(jc,jb)         = vmin_dist
        min_node_idx(:,jc,jb)   = vmin_node_idx(:)

        p_old(1:icoord_dim)     = p_new(1:icoord_dim)
        min_dist_old            = vmin_dist

      END DO
    END DO
  END SUBROUTINE gnat_query_list


  !> Multi-threaded query for a sequence of points
  !  @note:
  !  - the search radius parameter is optional
  !  - search radius is adapted in case of algorithm failure
  !
  SUBROUTINE gnat_query_list_mt(gnat, v, iv_nproma, iv_nblks, iv_npromz,       &
    &                           min_dist, min_node_idx, ladapt_radius, opt_rr)

    ! max. trial steps to adapt search radius, r_new = opt_rr*10**iadapt
    INTEGER, PARAMETER :: nadapt = 15

    TYPE (t_gnat_tree), INTENT(IN) :: gnat
    INTEGER,            INTENT(IN)    :: iv_nproma, iv_nblks, iv_npromz     ! list size
    REAL(gk),           INTENT(IN)    :: v(iv_nproma, iv_nblks, icoord_dim) ! list of search points
    REAL(gk),           INTENT(INOUT) :: min_dist(iv_nproma, iv_nblks)      ! minimal distance
    INTEGER,            INTENT(INOUT) :: min_node_idx(3,iv_nproma, iv_nblks)! corresponding node
    LOGICAL,            INTENT(IN)    :: ladapt_radius
    REAL(gk), OPTIONAL, INTENT(IN)    :: opt_rr                                 ! opt. search radius
    ! local parameters
    CHARACTER(*), PARAMETER :: routine = modname//"::gnat_query_list_mt"
    LOGICAL,      PARAMETER :: l_chunk = .FALSE. !< distributes bigger chunks to each thread

    INTEGER                           :: istart, ichunksize, iend,     &
    &                                    nproc, iproc, iadapt, jb, jc, &
    &                                    end_idx, nqueries, tree
    REAL(gk)                          :: r                    ! search radius
    REAL(gk)                          :: vmin_dist
    INTEGER                           :: vmin_node_idx(3)

    tree = gnat%gnat_tree

    IF (PRESENT(opt_rr)) THEN
      r = opt_rr
    ELSE
      ! determine a meaningful search radius
      r = gnat_std_radius(gnat) ! search radius ~ (query time)**(-1)
    END IF

#ifdef _OPENMP
    nproc = OMP_GET_MAX_THREADS()
#else
    nproc = 1
#endif

    IF (l_chunk) THEN
!$OMP PARALLEL private(istart, iend, ichunksize, iproc)
#ifdef _OPENMP
      iproc      = OMP_GET_THREAD_NUM() + 1
#else
      iproc      = 1
#endif

      ! proc "iproc" will query the entry blocks (istart,...,iend)
      ichunksize = (iv_nblks+nproc-1)/nproc
      istart     = (iproc-1)*ichunksize + 1
      iend       = MIN(iv_nblks, istart+ichunksize-1)
      IF (iend >= istart) THEN
        ichunksize = iend - istart + 1
        CALL gnat_query_list(gnat, v, iv_nproma, iv_nblks, iv_npromz, ichunksize, &
          &                  istart, r, min_dist, min_node_idx)
      END IF
!$OMP END PARALLEL
    ELSE
!$OMP PARALLEL
!$OMP DO SCHEDULE(DYNAMIC)
      DO istart=1,iv_nblks
        CALL gnat_query_list(gnat, v, iv_nproma, iv_nblks, iv_npromz, 1, &
          &                  istart, r, min_dist, min_node_idx)
      END DO
!$OMP END DO
!$OMP END PARALLEL
    END IF


    ! for all failed queries (if there are any) repeat query with a
    ! larger search radius:
    IF (.NOT. ladapt_radius) RETURN

    iadapt = 0
    DO
      IF ((iadapt > nadapt) .OR. &
        & (COUNT(min_node_idx(1,:,:) == INVALID_NODE) == 0)) EXIT
      ! adapt search radius
      iadapt = iadapt + 1
      r = r*2._gk
      IF (r > MAX_RANGE) EXIT
      IF (dbg_level > 1) CALL message(routine, "adapting radius")

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
!$OMP             private(end_idx,jb,jc,vmin_dist, vmin_node_idx)
      DO jb=1,iv_nblks

        ! set end index in current block:
        end_idx = iv_nproma
        IF (jb == iv_nblks) end_idx = iv_npromz

        DO jc=1,end_idx
          IF (min_node_idx(1,jc,jb) == INVALID_NODE) THEN
            vmin_dist        = MAX_RANGE
            vmin_node_idx(1) = INVALID_NODE
            CALL gnat_recursive_query(gnat, tree, v(jc,jb,:), r,     &
              &                       vmin_dist, vmin_node_idx(:))
            min_dist(jc,jb)         = vmin_dist
            min_node_idx(:,jc,jb)   = vmin_node_idx(:)
          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
    END DO

    IF (iadapt > nadapt) THEN
      WRITE(message_text,*) "failed to adapt search radius! Last radius=",r
      CALL message(routine, TRIM(message_text))
      ! CALL finish (routine, TRIM(message_text))
    END IF

    IF (dbg_level > 10) THEN
      IF (iadapt <= nadapt) THEN
        WRITE(message_text,*) "no. of adaptations: ", iadapt
        CALL message(routine, TRIM(message_text))
      END IF
      nqueries = (iv_nblks-1)*iv_nproma + iv_npromz
      WRITE(message_text,*) "  no. of queries: ", nqueries
      CALL message(routine, TRIM(message_text))
    END IF
  END SUBROUTINE gnat_query_list_mt


  !> Utility function: Divide-and-Conquer algorithm to determine
  !  points in lon-lat grid which are "far off" the current patch of
  !  the triangular grid.
  !
  RECURSIVE SUBROUTINE flag_ll_points(gnat, rotated_pts, s_lon, e_lon, s_lat, e_lat, &
    &                                 pts_flags, recursion_depth, max_recursion)
    TYPE (t_gnat_tree), INTENT(IN)    :: gnat
    REAL(wp),           INTENT(IN)    :: rotated_pts(:,:,:)  ! dim (lon,lat,1:2)
    INTEGER,            INTENT(IN)    :: s_lon, e_lon, s_lat, e_lat
    INTEGER,            INTENT(INOUT) :: pts_flags(:,:)
    INTEGER,            INTENT(IN)    :: recursion_depth, max_recursion
    ! local variables
    INTEGER  :: m_lon, m_lat, c_lon(4), c_lat(4), v_lon(4,4), v_lat(4,4), &
      &         s_lon2(4), s_lat2(4), e_lon2(4), e_lat2(4), icirc, ivert, &
      &         i1, i2
    REAL(gk) :: p1(4,2), p2(2), radius(4)
    LOGICAL  :: is_pt_inside(4)

    IF (recursion_depth == max_recursion) RETURN

    ! determine mid longitude and latitude in current range, thus
    ! dividing the regular grid into 4 subgrids:
    m_lon = (e_lon + s_lon)/2
    m_lat = (e_lat + s_lat)/2

    ! nomenclature
    !
    ! *-------*-------*   - e_lat
    ! |       |       |
    ! |   3   |   4   |
    ! |       |       |
    ! *-------*-------*   - m_lat
    ! |       |       |
    ! |   1   |   2   |
    ! |       |       |
    ! *-------*-------*   - s_lat
    !
    ! |       |       |
    !  s_lon   m_lon   e_lon

    s_lon2(1:4) = (/ s_lon, m_lon, s_lon, m_lon /)
    e_lon2(1:4) = (/ m_lon, e_lon, m_lon, e_lon /)
    s_lat2(1:4) = (/ s_lat, s_lat, m_lat, m_lat /)
    e_lat2(1:4) = (/ m_lat, m_lat, e_lat, e_lat /)
    ! for each of the four subblocks, determine centre point and
    ! vertex indices:
    DO icirc=1,4
      c_lon(icirc) = (s_lon2(icirc) + e_lon2(icirc)) / 2
      c_lat(icirc) = (s_lat2(icirc) + e_lat2(icirc)) / 2
      ! vertex numbering counter-clockwise
      v_lon(icirc,1:4) = (/ s_lon2(icirc), e_lon2(icirc), e_lon2(icirc), s_lon2(icirc) /)
      v_lat(icirc,1:4) = (/ s_lat2(icirc), s_lat2(icirc), e_lat2(icirc), e_lat2(icirc) /)
    END DO
    ! determine radii of circles overlapping each subblock:
    radius(:) = 0._gk
    DO icirc=1,4
      p1(icirc,1:2) = rotated_pts(c_lon(icirc),c_lat(icirc),1:2)
      DO ivert=1,4
        p2(1:2) = rotated_pts(v_lon(icirc,ivert),v_lat(icirc,ivert),1:2)
        radius(icirc) = MAX(radius(icirc), dist_p(p1(icirc,1:2), p2))
      END DO
    END DO
    radius(1:4) = 1.25_gk * radius(1:4) ! safety margin
    ! test, if there exists a triangle center inside the circles:
    DO icirc=1,4
      is_pt_inside(icirc) = .FALSE.
      CALL gnat_recursive_proximity_query(gnat, gnat%gnat_tree, p1(icirc,1:2), &
        &                                 radius(icirc), is_pt_inside(icirc))
    END DO
    ! for each of the four subblocks:
    DO icirc=1,4
      IF (.NOT. is_pt_inside(icirc)) THEN
        ! if patch of triangular grid does not intersect circle:
        ! mark corresponding points in regular grid with "SKIP_NODE"
        DO i2=s_lat2(icirc),e_lat2(icirc)
          DO i1=s_lon2(icirc),e_lon2(icirc)
            pts_flags(i1,i2) = SKIP_NODE
          END DO
        END DO
      ELSE
        ! if patch of triangular grid does intersect circle:
        ! recursion and further subdivision:
        CALL flag_ll_points(gnat, rotated_pts, s_lon2(icirc), e_lon2(icirc), s_lat2(icirc), e_lat2(icirc), &
          &                 pts_flags, (recursion_depth+1), max_recursion)
      END IF
    END DO
  END SUBROUTINE flag_ll_points


  ! ---------------------------------------------------------
  ! subroutines for identifying GNAT subtrees
  ! ---------------------------------------------------------


  !> Identify GNAT subtrees within given radius.
  !
  RECURSIVE SUBROUTINE gnat_find_subtrees(gnat, tree_idx, v, r, flag_field, flag_field_c, dist_cnt, &
    &                                     opt_ilevel, opt_imax_level)
    TYPE (t_gnat_tree), INTENT(IN), TARGET :: gnat
    INTEGER,            INTENT(IN)    :: tree_idx           !< index pointer
    REAL(gk),           INTENT(IN)    :: v(icoord_dim)      !< search point
    REAL(gk),           INTENT(IN)    :: r                  !< search radius
    INTEGER,            INTENT(INOUT) :: flag_field(:,:)    !< field (nnodes) where flags are set
    INTEGER,            INTENT(INOUT) :: flag_field_c(:,:)  !< field (nnodes) where flags are set
    INTEGER,            INTENT(INOUT) :: dist_cnt           !< no. of distance computations
    INTEGER,            INTENT(IN), OPTIONAL  :: opt_ilevel, opt_imax_level
    ! local variables
    LOGICAL                           :: pflag(gnat_k)
    INTEGER                           :: ip, isplit_pts, jc, jb, ilevel, imax_level
    REAL(gk)                          :: dist_vp, distv(gnat_k)
    TYPE(t_gnat_node), POINTER        :: tree

    ilevel     = 0
    imax_level = 1
    IF (PRESENT(opt_imax_level)) THEN
      ilevel     = opt_ilevel
      imax_level = opt_imax_level
    END IF

    tree => gnat%node_storage(tree_idx)
    isplit_pts = tree%isplit_pts
    IF (isplit_pts == 0)  RETURN

    ! include all split points in set P
    ! note that values for indices > tree%isplit_pts are undefined!
    pflag(1:isplit_pts) = .TRUE.

    CALL dist_vect(tree%p, v, isplit_pts, distv)
    dist_cnt = dist_cnt + isplit_pts
    ! remove some elements from P
    P : DO ip=1,isplit_pts
      dist_vp = distv(ip)
      jc = tree%p(ip)%idx(1)
      jb = tree%p(ip)%idx(2)
      IF (dist_vp > r) THEN
        IF (flag_field_c(jc,jb) == -1)  flag_field_c(jc,jb) = 1
      ELSE
        flag_field_c(jc,jb) = 0
      END IF

      WHERE (pflag(1:isplit_pts))
        pflag(1:isplit_pts) = (tree%drange(1:isplit_pts,ip,1) <= (dist_vp+r)) .AND. &
          &                   (tree%drange(1:isplit_pts,ip,2) >= (dist_vp-r))
      END WHERE
    END DO P

    DO ip=1,isplit_pts
      jc = tree%child(ip)
      IF (flag_field(tree_idx,ip) == -1)  flag_field(tree_idx,ip) = 0
      IF ((jc /= UNASSOCIATED) .AND. (pflag(ip))) THEN
        ! increase subtree counter (flag_field)
        flag_field(tree_idx,ip) = flag_field(tree_idx,ip) + 1
        ! traverse subtree
        IF (ilevel < imax_level) &
          &   CALL gnat_find_subtrees(gnat, jc, v, r, flag_field, flag_field_c, dist_cnt, &
          &                           ilevel+1, opt_imax_level)
      END IF
    END DO
  END SUBROUTINE gnat_find_subtrees


  !> Identify GNAT subtrees within given radius.
  !
  RECURSIVE SUBROUTINE gnat_find_subtrees2(gnat, tree_idx, v, r, flag_field, nflag_field, &
    &                                      flag_field_c, nflag_field_c,                   &
    &                                      dist_cnt, opt_ilevel, opt_imax_level)
    TYPE (t_gnat_tree), INTENT(IN), TARGET :: gnat
    INTEGER,            INTENT(IN)    :: tree_idx           !< index pointer
    REAL(gk),           INTENT(IN)    :: v(icoord_dim)      !< search point
    REAL(gk),           INTENT(IN)    :: r                  !< search radius
    INTEGER,            INTENT(INOUT) :: flag_field(:,:)    !< field (nnodes) where flags are set
    INTEGER,            INTENT(INOUT) :: nflag_field
    INTEGER,            INTENT(INOUT) :: flag_field_c(:,:)  !< field (nnodes) where flags are set
    INTEGER,            INTENT(INOUT) :: nflag_field_c
    INTEGER,            INTENT(INOUT) :: dist_cnt           !< no. of distance computations
    INTEGER,            INTENT(IN), OPTIONAL  :: opt_ilevel, opt_imax_level
    ! local variables
    LOGICAL                           :: pflag(gnat_k)
    INTEGER                           :: ip, isplit_pts, jc, jb, ilevel, imax_level
    REAL(gk)                          :: dist_vp, distv(gnat_k)
    TYPE(t_gnat_node), POINTER        :: tree

    ilevel     = 0
    imax_level = 1
    IF (PRESENT(opt_imax_level)) THEN
      ilevel     = opt_ilevel
      imax_level = opt_imax_level
    END IF

    tree => gnat%node_storage(tree_idx)
    isplit_pts = tree%isplit_pts
    IF (isplit_pts == 0)  RETURN

    ! include all split points in set P
    ! note that values for indices > tree%isplit_pts are undefined!
    pflag(1:isplit_pts) = .TRUE.

    CALL dist_vect(tree%p, v, isplit_pts, distv)
    dist_cnt = dist_cnt + isplit_pts
    ! remove some elements from P
    P : DO ip=1,isplit_pts
      dist_vp = distv(ip)
      jc = tree%p(ip)%idx(1)
      jb = tree%p(ip)%idx(2)
      IF (dist_vp <= r) THEN
        nflag_field_c = nflag_field_c + 1
        flag_field_c(:, nflag_field_c) = (/ jc, jb /)
      END IF

      WHERE (pflag(1:isplit_pts))
        pflag(1:isplit_pts) = (tree%drange(1:isplit_pts,ip,1) <= (dist_vp+r)) .AND. &
          &                   (tree%drange(1:isplit_pts,ip,2) >= (dist_vp-r))
      END WHERE
    END DO P

    DO ip=1,isplit_pts
      jc = tree%child(ip)
      IF (jc /= UNASSOCIATED) THEN
        IF (pflag(ip)) THEN
          ! traverse subtree
          IF (ilevel < imax_level) &
            &   CALL gnat_find_subtrees2(gnat, jc, v, r, flag_field, nflag_field, flag_field_c, nflag_field_c, &
            &                            dist_cnt, ilevel+1, opt_imax_level)
        ELSE
          ! increase subtree counter (flag_field)
          nflag_field = nflag_field + 1
          flag_field(1:2, nflag_field) = (/ tree_idx,ip /)
        END IF
      END IF
    END DO
  END SUBROUTINE gnat_find_subtrees2


  !> Mark all triangles in a subtree of the GNAT data structure
  !
  RECURSIVE SUBROUTINE gnat_flag_subtree(gnat, root_node, root_idx, flag_field_c)
    TYPE (t_gnat_tree), INTENT(INOUT), TARGET :: gnat
    INTEGER,            INTENT(IN)    :: root_node, root_idx  !< root node, split point idx of current subtree
    INTEGER,            INTENT(INOUT) :: flag_field_c(:,:)    !< field (nproma,nblks_c) where flags are set

    ! local variables
    INTEGER                     :: ip, isplit_pts, idx(2), child
    TYPE(t_gnat_node), POINTER  :: tree

    ! set flag:
    tree => gnat%node_storage(root_node)
    child = tree%child(root_idx)

    ! traverse subtree
    IF (child /= UNASSOCIATED) THEN
      idx(1:2) = tree%p(root_idx)%idx(1:2)
      flag_field_c(idx(1), idx(2)) = 1
      isplit_pts = gnat%node_storage(child)%isplit_pts
      DO ip=1,isplit_pts
        idx(1:2) = gnat%node_storage(child)%p(ip)%idx(1:2)
        flag_field_c(idx(1), idx(2)) = 1
        CALL gnat_flag_subtree(gnat, child, ip, flag_field_c)
      END DO
    END IF
  END SUBROUTINE gnat_flag_subtree


  ! ---------------------------------------------------------
  ! utility subroutines and functions
  ! ---------------------------------------------------------


  !> Determine a meaningful value for the search radius.
  !
  FUNCTION gnat_std_radius(gnat)
    REAL(gk) :: gnat_std_radius
    TYPE (t_gnat_tree), TARGET, INTENT(IN) :: gnat
    ! local parameters
    TYPE (t_gnat_node), POINTER            :: node
    INTEGER                                :: i

    ! move down tree an collect range values
    node => gnat%node_storage(gnat%gnat_tree)
    i = node%isplit_pts
    gnat_std_radius = MINVAL(node%drange(2:i, 1, 1))
    DO
      i = node%isplit_pts
      gnat_std_radius = MIN(gnat_std_radius, MINVAL(node%drange( 2:i, 1, 1)))

      IF (node%child(1) == UNASSOCIATED) EXIT
      node => gnat%node_storage(node%child(1))
    END DO
    ! increase radius for safety
    gnat_std_radius = 2._gk*gnat_std_radius
  END FUNCTION gnat_std_radius


  !> Utility function.
  !  @return index to a free node, allocates if necessary.
  !
  FUNCTION get_free_node(gnat)
    INTEGER :: get_free_node
    TYPE (t_gnat_tree), INTENT(INOUT) :: gnat

    gnat%num_nodes = gnat%num_nodes + 1
    get_free_node  = gnat%num_nodes
  END FUNCTION get_free_node


  !> Utility function: allocates space for tree nodes.
  !
  SUBROUTINE resize_node_storage(gnat)
    TYPE (t_gnat_tree), INTENT(INOUT) :: gnat

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::resize_node_storage"
    INTEGER :: new_max_num_nodes, errstat, i
    TYPE (t_gnat_node), TARGET, ALLOCATABLE :: tmp(:)

    ! triangle copy
    new_max_num_nodes = gnat%max_num_nodes + gnat%expected_num_nodes

    IF (ALLOCATED(gnat%node_storage)) THEN
      IF (gnat%max_num_nodes > SIZE(gnat%node_storage)) &
        & CALL finish (routine, 'Unexpected array size!')
      ALLOCATE(tmp(gnat%max_num_nodes), STAT=errstat)
      IF (errstat /= 0) CALL finish (routine, 'Error in ALLOCATE operation!')
      tmp(1:gnat%max_num_nodes) = gnat%node_storage(1:gnat%max_num_nodes)
      DEALLOCATE(gnat%node_storage, STAT=errstat)
      IF (errstat /= 0) CALL finish (routine, 'Error in DEALLOCATE operation!')
    END IF

    IF (dbg_level > 1) THEN
      WRITE(message_text,*) "allocating node storage of size=", new_max_num_nodes
      CALL message(routine, TRIM(message_text))
    END IF
    ALLOCATE(gnat%node_storage(new_max_num_nodes), STAT=errstat)
    IF (errstat /= 0)  CALL finish (routine, 'Error in ALLOCATE operation!')
    FORALL (i=(gnat%max_num_nodes+1):new_max_num_nodes)
      gnat%node_storage(i)%child(:)   = UNASSOCIATED
      gnat%node_storage(i)%isplit_pts = 0
    END FORALL
    IF (ALLOCATED(tmp)) THEN
      gnat%node_storage(1:gnat%max_num_nodes) = tmp(1:gnat%max_num_nodes)
      DEALLOCATE(tmp, STAT=errstat)
      IF (errstat /= 0)  &
        CALL finish (routine, 'Error in DEALLOCATE operation!')
    END IF
    gnat%max_num_nodes = new_max_num_nodes
  END SUBROUTINE resize_node_storage


  ! ---------------------------------------------------------
  ! geometric kernel functions, e.g. for distance measurement
  ! ---------------------------------------------------------


  !> distance computation on the sphere, uses precomputed sine, cosine
  !
  PURE FUNCTION dist(p1, p2)
    REAL(gk)                    :: dist
    TYPE (t_coord), INTENT(IN)  :: p1
    REAL(gk),       INTENT(IN)  :: p2(icoord_dim)
    ! local variables
    REAL(gk) :: val

    ! spherical distance:
    val = p1%sin_p*SIN(p2(2)) + p1%cos_p*COS(p2(2))*COS(p1%p(1)-p2(1))
    dist = ACOS( MIN(1._gk, MAX(-1._gk, val)) )
  END FUNCTION dist


  !> distance computation on the sphere, uses precomputed sine, cosine
  !
  PURE FUNCTION dist_t_coord(p1, p2)
    REAL(gk)                    :: dist_t_coord
    TYPE (t_coord), INTENT(IN)  :: p1, p2
    ! local variables
    REAL(gk) :: val

    ! spherical distance:
    val = p1%sin_p*p2%sin_p + p1%cos_p*p2%cos_p*COS(p1%p(1)-p2%p(1))
    dist_t_coord = ACOS( MIN(1._gk, MAX(-1._gk, val)) )
  END FUNCTION dist_t_coord


  !> distance computation on the sphere; no precomputed values
  !> required.
  !
  PURE FUNCTION dist_p(p1, p2)
    REAL(gk)              :: dist_p
    REAL(gk), INTENT(IN)  :: p1(icoord_dim), p2(icoord_dim)
    ! local variables
    REAL(gk) :: val

    ! spherical distance:
    val = SIN(p1(2))*SIN(p2(2)) + COS(p1(2))*COS(p2(2))*COS(p1(1)-p2(1))
    dist_p = ACOS( MIN(1._gk, MAX(-1._gk, val)) )
  END FUNCTION dist_p


  !> vector variant of distance computation on the sphere,
  !> uses precomputed sine, cosine.
  !
  PURE SUBROUTINE dist_vect(p1, p2, n, pdist)
    INTEGER       , INTENT(IN)     :: n
    TYPE (t_coord), INTENT(IN)     :: p1(n)
    REAL(gk)      , INTENT(IN)     :: p2(icoord_dim)
    REAL(gk)      , INTENT(INOUT)  :: pdist(n)
    INTEGER                        :: i
    REAL(gk)                       :: sin_p2, cos_p2
    ! local variables
    REAL(gk) :: val(n)

    sin_p2 = SIN(p2(2))
    cos_p2 = COS(p2(2))

    ! spherical distance:
    FORALL (i=1:n)
      val(i)   = (p1(i)%sin_p*sin_p2 +  &
        &         p1(i)%cos_p*cos_p2*COS(p1(i)%p(1)-p2(1)))
      pdist(i) = ACOS( MIN(1._gk, MAX(-1._gk, val(i))) )
    END FORALL
  END SUBROUTINE dist_vect


  ! ---------------------------------------------------------
  ! internal subroutines for tree construction
  ! ---------------------------------------------------------


  !> Dynamic insertion of a point into GNAT-based tree.
  !  Non-recursive version, suitable for parallel processing
  !  note: sequence of insertion points should be more or less random
  !  in order to provide split points which are fairly far apart.
  !
  SUBROUTINE gnat_insert(gnat, tree, p, idx, free_node, lcomplete)
    TYPE (t_gnat_tree), INTENT(INOUT), TARGET :: gnat
    INTEGER,            INTENT(INOUT) :: tree ! index pointer
    REAL(gk),           INTENT(IN)    :: p(icoord_dim)
    INTEGER,            INTENT(IN)    :: idx(3)
    INTEGER,            INTENT(INOUT) :: free_node
    LOGICAL,            INTENT(OUT)   :: lcomplete

    INTEGER                           :: i, j, imin(1)
    REAL(gk)                          :: pdist(gnat_k)
    TYPE (t_gnat_node), POINTER       :: node

    node => gnat%node_storage(tree)

    !-- case 1: tree is an empty/incomplete node
    ! compute distances
    CALL dist_vect(node%p, p, node%isplit_pts, pdist)

    ! add point p as split point
    IF (node%isplit_pts < gnat_k) THEN
      i = node%isplit_pts + 1
      node%isplit_pts = i
      node%p(i)%p(:)  = p(:)
      node%p(i)%sin_p = SIN(p(2))
      node%p(i)%cos_p = COS(p(2))
      node%p(i)%idx   = idx
      node%drange(i,i,1:2) = 0._gk

      FORALL (j=1:i-1)
        node%drange(i,j,1) = MIN( node%drange(i,j,1), pdist(j) )
        node%drange(i,j,2) = MAX( node%drange(i,j,2), pdist(j) )
        node%drange(j,i,1:2) = pdist(j)
      END FORALL
      lcomplete = .TRUE.
      RETURN
    END IF

    !-- case 2: traverse into subtree for closest split point
    imin  = MINLOC(pdist)
    FORALL (j=1:gnat_k)
      node%drange(imin,j,1) = MIN( node%drange(imin,j,1), pdist(j) )
      node%drange(imin,j,2) = MAX( node%drange(imin,j,2), pdist(j) )
    END FORALL

    ! instead of recursive traversal, like
    ! CALL gnat_insert(tree%child(imin(1))%ptr, p, idx)
    ! we return with status lcomplete == .FALSE.
    tree = node%child(imin(1))
    IF (tree == UNASSOCIATED) THEN
      ! use the free node provided by the calling
      ! function to store the new node child:
      node%child(imin(1)) = free_node
      free_node = UNASSOCIATED
      tree = node%child(imin(1))
      node => gnat%node_storage(tree)
      node%isplit_pts    = 0
      node%drange(:,:,1) = MAX_RANGE
      node%drange(:,:,2) = 0._gk
      node%child(:) = UNASSOCIATED
    END IF

    lcomplete = .FALSE.
  END SUBROUTINE gnat_insert


  !> Multi-threaded insertion of points into GNAT data structure
  !
  ! @note The construction of the GNAT data structure is aimed at
  !       interpolation purposes. Thus we insert only cells with
  !       "ref_ctrl" flags >= 2.
  !
  SUBROUTINE gnat_insert_mt(gnat, p_patch, nproc,  &
    &                       opt_ldegree, opt_startblk, opt_endblk)

    TYPE (t_gnat_tree), INTENT(INOUT), TARGET :: gnat
    TYPE(t_patch),      INTENT(IN)            :: p_patch
    INTEGER,            INTENT(IN)            :: nproc
    LOGICAL, INTENT(IN), OPTIONAL :: opt_ldegree
    INTEGER, INTENT(IN), OPTIONAL :: opt_startblk, opt_endblk

    ! Local parameters
    CHARACTER(*), PARAMETER :: routine = modname//"::gnat_insert_mt"

    ! reordering parameter (for permuting points before insertion):
    INTEGER                    :: nstrides
    INTEGER                    :: root ! index pointer
    INTEGER                    :: idx(3,nproc), work(nproc), depth(nproc)
    INTEGER                    :: node_proc(nproc) ! index pointer
    INTEGER                    :: free_node(nproc) ! short list of available nodes
    REAL(gk)                   :: p(icoord_dim,nproc)
    INTEGER                    :: j, new_idx(3), jb, jc
    INTEGER                    :: iproc
    LOGICAL                    :: lcomplete
    TYPE (t_geographical_coordinates) :: pcoord
    TYPE (t_gnat_node), POINTER:: node
    INTEGER                    :: i_startblk, i_endblk, &
      &                           i_startidx, i_endidx, &
      &                           rl_start, rl_end, i_nchdom
    LOGICAL                    :: l_loop_end, ldegree
    INTEGER, ALLOCATABLE       :: cell_indices(:,:), permutation(:) ! 1D array of cell points
    INTEGER                    :: ntotal, errstat, icount, tree_depth

    IF (dbg_level > 10) THEN
        WRITE(message_text,*) "Running with ", nproc, "thread(s)."
        CALL message(routine, TRIM(message_text))
    END IF

    ldegree = .FALSE.
    IF (PRESENT(opt_ldegree)) ldegree = opt_ldegree

    ! values for the blocking
    rl_start = 2
    rl_end = min_rlcell_int

    IF (PRESENT(opt_startblk) .AND. PRESENT(opt_endblk)) THEN
      i_startblk = opt_startblk
      i_endblk   = opt_endblk
    ELSE
      i_nchdom   = MAX(1,p_patch%n_childdom)
      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    END IF

    ! create a 1D list of cell indices:
    ALLOCATE(cell_indices(3,p_patch%n_patch_cells), &
      &      permutation(p_patch%n_patch_cells),STAT=errstat)
    IF (errstat /= 0)  CALL finish (routine, 'Error in ALLOCATE operation!')

    ntotal = 0
    IF (PRESENT(opt_startblk) .AND. PRESENT(opt_endblk)) THEN
      DO jb=i_startblk,i_endblk
        i_startidx = 1
        i_endidx   = UBOUND(p_patch%cells%center,1) ! = nproma
        IF (jb == i_endblk) i_endidx = p_patch%npromz_c
        DO jc=i_startidx,i_endidx
          ntotal = ntotal + 1
          cell_indices(1,ntotal) = jc
          cell_indices(2,ntotal) = jb
          cell_indices(3,ntotal) = p_patch%cells%decomp_info%glb_index(idx_1d(jc,jb))
        END DO
      END DO
    ELSE
      DO jb=i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb,  &
          &                i_startblk, i_endblk, &
          &                i_startidx, i_endidx, &
          &                rl_start, rl_end)
        DO jc=i_startidx,i_endidx
          IF(.NOT. p_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
          ntotal = ntotal + 1
          cell_indices(1,ntotal) = jc
          cell_indices(2,ntotal) = jb
          cell_indices(3,ntotal) = p_patch%cells%decomp_info%glb_index(idx_1d(jc,jb))
        END DO
      END DO
    END IF
    icount = 0

    ! set a (somewhat empirical) value to improve
    ! the distribution of split points:
    nstrides = NINT(SQRT(REAL(ntotal,wp)))

    DO jb=1,nstrides
      DO jc=jb,ntotal,nstrides
        icount = icount + 1
        permutation(icount) = jc
      END DO
    END DO

    ! create a root node
    CALL resize_node_storage(gnat)
    root = get_free_node(gnat)
    node => gnat%node_storage(root)
    node%isplit_pts    = 0
    node%drange(:,:,1) = MAX_RANGE
    node%drange(:,:,2) = 0._gk
    node%child(:) = UNASSOCIATED

    idx (:,:)    = -1  ! indices of nodes processed in parallel
    work(:)      =  0  ! counter: no. of pts. inserted by procs

    depth(:)     =  0
    tree_depth   =  0  ! maximum tree depth
    free_node(:) = UNASSOCIATED
    icount       =  0

    IF (p_patch%n_patch_cells > 0) &   ! needed to skip processes not owning any cells
      new_idx(:)   = cell_indices(:,permutation(1))

    ! Short description of parallel processing:

    ! 0 mark all processes as "idle" (idx == -1)
    ! 1 insert one point into the "pipeline" of points processed in parallel.
    ! 2 move all currently processed points one level down the tree
    !   this happens in parallel, since we know that all active
    !   nodes are on a different level of the GNAT structure.
    ! 3 if one point has reached its final position in the
    !   GNAT structure, mark the corresponding process as "idle" again
    ! 4 repeat steps 1-3 until all points have been inserted

    ! disadvantage of this simple pipelining: for a moderate number
    ! of nodes in the GNAT, the tree is rather flat, thus processes
    ! may be often idle since all levels are currently occupied.
    ! In other words: This algorithm does not scale well beyond ~8 threads!

    l_loop_end = .FALSE.

!$OMP PARALLEL num_threads(nproc)                    &
!$OMP          private(lcomplete,iproc)

#ifdef _OPENMP
      iproc = OMP_GET_THREAD_NUM() + 1
#else
      iproc = 1
#endif

    POINTS : DO

!$OMP BARRIER
!$OMP MASTER
      !-- block executed by first thread:
      !-- resize storage for tree nodes if necessary
      IF (gnat%num_nodes >= gnat%max_num_nodes-2*nproc) THEN
        CALL resize_node_storage(gnat)
      END IF

      !-- build a short list of available node
      !   indices for NEW tree nodes
      DO j=1,nproc
        IF (free_node(j) == UNASSOCIATED) THEN
          free_node(j) = get_free_node(gnat)
        END IF
      END DO

      !-- assign new work to one of the PEs
      ! find first non-associated PE in node_proc
      FINDFIRST : DO j = 1, nproc
        IF (idx(1,j) == -1) EXIT FINDFIRST
      END DO FINDFIRST

      IF (j<=nproc) THEN
        ! insert next point for this proc
        icount = icount + 1
        IF (icount <= ntotal) THEN
          new_idx(:)   = cell_indices(:,permutation(icount))
          pcoord       = p_patch%cells%center(new_idx(1), new_idx(2))
          p(1,j)       = REAL( pcoord%lon, gk)
          p(2,j)       = REAL( pcoord%lat, gk)
          IF (ldegree) p(:,j) = p(:,j)*pi_180
          idx(:,j)     = new_idx(:)
          work(j)      = work(j) + 1
          tree_depth   = MAX(tree_depth, depth(j))
          depth(j)     = 0
          node_proc(j) = root
        END IF
      END IF
      IF ((icount > ntotal) .AND. (ALL(idx(1,:) == -1))) THEN
        !EXIT POINTS
        l_loop_end = .TRUE.
      END IF
!$OMP END MASTER
!$OMP BARRIER

      IF (l_loop_end) EXIT POINTS

      ! step all active nodes one step down the tree
      IF (idx(1,iproc) /= -1) THEN
        CALL gnat_insert(gnat, node_proc(iproc), p(:,iproc), idx(:,iproc), &
          &              free_node(iproc), lcomplete)
        ! "lcomplete" : proc is free for new point
        IF (lcomplete) THEN
          idx(1,iproc) = -1
        ELSE
          depth(iproc) = depth(iproc) + 1
        END IF
      END IF
    END DO POINTS
!$OMP END PARALLEL

    gnat%gnat_tree = root
    IF (dbg_level > 10) THEN
      WRITE(0,*) "distribution of insertions over threads: ", work(1:nproc), &
        &        ", maximum tree depth:", tree_depth
    END IF

    DEALLOCATE(cell_indices, permutation, STAT=errstat)
    IF (errstat /= 0)  CALL finish (routine, 'Error in DEALLOCATE operation!')
  END SUBROUTINE gnat_insert_mt


  ! ---------------------------------------------------------
  ! subroutines for distributed point query with MPI
  ! ---------------------------------------------------------

  ! Perform a nearest-neighbor query for a given list of points and
  ! return the indices and block indices of the mesh triangles that
  ! contain these points.
  SUBROUTINE gnat_query_containing_triangles(gnat, p_patch, v, iv_nproma, iv_nblks,       &
    &                                        iv_npromz, grid_sphere_radius, l_p_test_run, &
    &                                        tri_idx, min_dist)

    TYPE (t_gnat_tree),    INTENT(IN)    :: gnat
    TYPE(t_patch), INTENT(IN) :: p_patch
    INTEGER,  INTENT(IN)    :: iv_nproma, iv_nblks, iv_npromz      ! list size
    REAL(wp), INTENT(IN)    :: grid_sphere_radius
    LOGICAL,  INTENT(IN)    :: l_p_test_run
    REAL(gk), INTENT(IN)    :: v(iv_nproma, iv_nblks, icoord_dim)  ! list of search points
    INTEGER,  INTENT(INOUT) :: tri_idx(2,iv_nproma, iv_nblks)      ! containing triangle (idx,block)
    REAL(gk), INTENT(OUT)   :: min_dist(iv_nproma, iv_nblks)       ! minimal distance
    ! local parameters
    INTEGER                 :: i_nv
    CHARACTER(*), PARAMETER :: routine = modname//"::gnat_query_containing_triangles"

    REAL(gk)                :: radius
    INTEGER                 :: bcast_source, i_startblk, &
      &                        i_startidx, i_endidx, &
      &                        rl_start, rl_end, i_nchdom
    LOGICAL :: l_my_process_is_mpi_test

    ! set default value ("failure notice")
    min_dist(:,:)  = MAX_RANGE

    IF (p_patch%n_patch_cells == 0) THEN
      ! ensure that processes without cells are involved in global communication 
      radius = 0._gk
      radius = p_max(radius, comm=p_comm_work)
      RETURN
    END IF

    i_nv  =  p_patch%geometry_info%cell_type

    IF (i_nv /= 3)  CALL finish(routine, "Wrong number of cell vertices!")

    ! make a sensible guess for search radius (just taking a "randomly
    ! chosen" edge length)
    rl_start = 2
    rl_end   = min_rlcell_int
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    CALL get_indices_e(p_patch, i_startblk, i_startblk, i_startblk, &
      &                i_startidx, i_endidx, rl_start, rl_end)
    radius = 3._gk * REAL(p_patch%edges%primal_edge_length(i_startidx,i_startblk) &
      & /grid_sphere_radius, gk)
    ! for MPI-independent behaviour: determine global max. of search radii
    radius = p_max(radius, comm=p_comm_work)
#ifndef NOMPI
    IF (l_p_test_run) THEN
      IF (my_process_is_mpi_test()) THEN
        bcast_source=0
      ELSE
        IF (p_pe_work == 0) THEN
          bcast_source=mpi_root
        ELSE
          bcast_source=mpi_proc_null
        END IF
      END IF
      CALL p_bcast(radius, bcast_source, comm=p_comm_work_2_test)
    END IF
#endif

    ! query list of nearest neighbors
    CALL gnat_query_nnb(gnat, v, iv_nproma, iv_nblks,   &
      &                 iv_npromz, radius, tri_idx, min_dist)

  END SUBROUTINE gnat_query_containing_triangles


  !> Merge the results of distributed proximity queries.
  !
  !  In this routine, the master process merges the results of several
  !  distributed nearest neighbor searches, i.e. we reduce multiple
  !  lists of pairs (distance, index) of the same length to a single
  !  list. All information except the one with minimal distance is
  !  discarded.
  !
  SUBROUTINE gnat_merge_distributed_queries(p_patch, total_dim, iv_nproma, iv_nblks, min_dist, &
    &                                       tri_idx, lonlat_points, global_idx,                &
    &                                       ithis_local_pts)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    INTEGER,               INTENT(IN)    :: total_dim
    INTEGER,               INTENT(IN)    :: iv_nproma, iv_nblks
    REAL(gk),              INTENT(IN)    :: min_dist(:,:)
    INTEGER,               INTENT(INOUT) :: tri_idx(:,:,:)
    REAL(gk),              INTENT(INOUT) :: lonlat_points(:,:,:)
    INTEGER,               INTENT(INOUT) :: global_idx(:)
    INTEGER,               INTENT(OUT)   :: ithis_local_pts !< no. of points on this PE
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::gnat_merge_distributed_queries"
    INTEGER  :: i, j, jc, jb
    INTEGER  :: new_tri_idx(2,total_dim)
    REAL(gk) :: new_lonlat_points(total_dim,2)
    INTEGER  :: iowner, my_id, gidx
    INTEGER  :: array_shape(2)
    TYPE(t_mpi_mintype) :: in(total_dim)

    my_id = get_my_mpi_work_id()

    ! 1. Perform an MPI_ALLREDUCE operation with the min_dist vector
    !    to find out which process is in charge of which in_point

    in(:)%rdist = RESHAPE(min_dist(:,:), (/ total_dim /) )
    in(:)%owner = my_id
    IF (p_patch%n_patch_cells == 0) in(:)%owner = -1;
    in(:)%glb_index = -1

    !NEC_RP: prevent processes without cells from accessing non initialized arrays
    IF (p_patch%n_patch_cells > 0) THEN
      DO j=1,total_dim
        ! convert global index into local idx/block pair:
        jb = (j-1)/iv_nproma + 1
        jc = j - (jb-1)*iv_nproma

        IF (tri_idx(1,jc,jb) /= INVALID_NODE) THEN
          gidx = idx_1d(tri_idx(1,jc,jb), tri_idx(2,jc,jb))
          in(j)%glb_index = p_patch%cells%decomp_info%glb_index(gidx)
        END IF
      END DO
    ENDIF

    ! call user-defined parallel reduction operation
    CALL mpi_reduce_mindistance_pts(in, total_dim, p_comm_work)

    ! 2. If we are a working PE, reduce the list of in_points and the
    !    tri_idx array to points which are actually located on this
    !    portion of the domain.

    ! store list of global indices owned by this PE:
    j=0
    DO i=1,total_dim
      iowner = in(i)%owner
      IF (iowner == my_id) THEN
        jb = (i-1)/iv_nproma + 1
        jc = i - (jb-1)*iv_nproma

        IF (tri_idx(1,jc,jb) /= INVALID_NODE) THEN
          j = j + 1
          global_idx(j)            = i
          new_tri_idx(1:2,j)       = tri_idx(1:2,jc,jb)
          new_lonlat_points(j,1:2) = lonlat_points(jc,jb,1:2)
        END IF
      END IF

    END DO

    ! now, j points are left for proc "get_my_mpi_work_id()"
    ithis_local_pts = j

    array_shape(:) = (/ iv_nproma, iv_nblks /)
    tri_idx(1,:,:) = RESHAPE( new_tri_idx(1,:), array_shape(:), (/ 0 /) )
    tri_idx(2,:,:) = RESHAPE( new_tri_idx(2,:), array_shape(:), (/ 0 /) )
    lonlat_points(:,:,1) = RESHAPE( new_lonlat_points(:,1), array_shape(:), (/ 0._gk /) )
    lonlat_points(:,:,2) = RESHAPE( new_lonlat_points(:,2), array_shape(:), (/ 0._gk /) )
  END SUBROUTINE gnat_merge_distributed_queries

END MODULE mo_gnat_gridsearch
