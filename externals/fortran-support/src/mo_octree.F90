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

!! Utility module: Puts geometric objects identified by indices into a
!! densely stored octree, based on their bounding boxes.

MODULE mo_octree

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64
  USE mo_exception, ONLY: finish
  USE mo_util_sort, ONLY: quicksort
#ifdef __SX__
  USE mo_util_sort, ONLY: radixsort
#endif
  IMPLICIT NONE

  ! subroutines and functions
  PUBLIC :: octree_init
  PUBLIC :: octree_finalize
  PUBLIC :: octree_count_point
  PUBLIC :: octree_query_point
  ! data types
  PUBLIC :: t_range_octree
  ! constants
  PUBLIC :: OCTREE_DEPTH

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_octree'
  INTEGER, PARAMETER :: SUCCESS = 0

  INTEGER, PARAMETER :: OCTREE_DEPTH = 7 !< octree depth, where root box is level 1
  INTEGER, PARAMETER :: OCTREE_SIZE = (8**OCTREE_DEPTH - 1)/7 !< no. of boxes in octree

  TYPE t_range_octree
    INTEGER, ALLOCATABLE :: object_idx(:) !< 1D list of objects inserted into the octree
    INTEGER, ALLOCATABLE :: object_data(:) !< 1D list of object data
    REAL(wp), ALLOCATABLE :: obj_min(:, :), obj_max(:, :) !< bounding boxes for inserted objects
    INTEGER, ALLOCATABLE :: box(:, :) !< pairs (box,startindex) for accessing object_idx
    REAL(wp)              :: brange0(2, 3) !< box range (min/max, dim=1,2,3)
    REAL(wp)              :: delta !< resolution of this octree (smallest box length)
  END TYPE t_range_octree

CONTAINS

  !> Initialize data structure.
  !
  SUBROUTINE octree_init(octree, brange, pmin, pmax, opt_index)
    TYPE(t_range_octree), INTENT(INOUT) :: octree !< octree data structure
    REAL(wp), INTENT(IN)    :: brange(2, 3) !< box range (min/max, dim=1,2,3)
    REAL(wp), INTENT(IN)    :: pmin(:, :), pmax(:, :) !< dim=(1,...,nobjects, x/y/z)  : range (corners)
    INTEGER, OPTIONAL, INTENT(IN)    :: opt_index(:) !< optional index array
    octree%brange0(:, :) = brange
    octree%delta = MINVAL(brange(2, :) - brange(1, :))/(2._wp**(OCTREE_DEPTH - 1))
    CALL octree_insert(octree, pmin, pmax, opt_index)
  END SUBROUTINE octree_init

  !> Deallocate data structure.
  !
  SUBROUTINE octree_finalize(octree)
    TYPE(t_range_octree), INTENT(INOUT) :: octree !< octree data structure
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::octree_finalize'
    INTEGER :: ierrstat

    DEALLOCATE (octree%object_idx, octree%object_data, octree%obj_min, &
      &        octree%obj_max, octree%box, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed')
  END SUBROUTINE octree_finalize

  !> Insert a list of objects into the octree.
  !
  SUBROUTINE octree_insert(octree, pmin, pmax, opt_index)
    TYPE(t_range_octree), INTENT(INOUT) :: octree !< octree data structure
    REAL(wp), INTENT(IN)    :: pmin(:, :), pmax(:, :) !< dim=(1,...,nobjects, x/y/z)  : range (corners)
    INTEGER, OPTIONAL, INTENT(IN)    :: opt_index(:) !< optional index array
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::octree_insert'
    INTEGER :: nobjects, nboxes, i, ierrstat
    INTEGER, ALLOCATABLE :: box(:)

    nobjects = SIZE(pmin, 1)
    IF ((SIZE(pmin, 1) /= SIZE(pmax, 1)) .OR. &
      & (SIZE(pmin, 2) /= 3) .OR. (SIZE(pmax, 2) /= 3)) THEN
      CALL finish(routine, "Inconsistent array sizes!")
    END IF
    IF (ALLOCATED(octree%object_idx)) THEN
      CALL finish(routine, "Internal error!")
    END IF
    ! allocate working arrays:
    ALLOCATE (box(nobjects), octree%object_idx(nobjects), octree%object_data(nobjects), &
      &      octree%obj_min(nobjects, 3), octree%obj_max(nobjects, 3), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed')
    DO i = 1, nobjects
      box(i) = octree_query_range(octree, pmin(i, :), pmax(i, :))
      octree%object_idx(i) = i
    END DO
    IF (PRESENT(opt_index)) THEN
      DO i = 1, nobjects
        octree%object_data(i) = opt_index(i)
      END DO
    ELSE
      DO i = 1, nobjects
        octree%object_data(i) = i
      END DO
    END IF
    octree%obj_min(1:nobjects, :) = pmin(1:nobjects, :)
    octree%obj_max(1:nobjects, :) = pmax(1:nobjects, :)
#ifdef __SX__
    CALL radixsort(box, perm=octree%object_idx)
#else
    CALL quicksort(box, permutation=octree%object_idx)
#endif
    ! count no. of occupied boxes
    nboxes = 1
    DO i = 2, nobjects
      IF (box(i) /= box(i - 1)) nboxes = nboxes + 1
    END DO
    ALLOCATE (octree%box(nboxes + 1, 2), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed')
    !NEC_RP: prevent processes not owning any grid points from accessing not allocated data
    IF (nobjects > 0) &
      octree%box(1, 1:2) = (/box(1), 1/)
    octree%box(nboxes + 1, 1:2) = (/OCTREE_SIZE + 1, nobjects + 1/)
    nboxes = 1
    DO i = 2, nobjects
      IF (box(i) /= box(i - 1)) THEN
        nboxes = nboxes + 1
        octree%box(nboxes, 1:2) = (/box(i), i/)
      END IF
    END DO
    ! clean up
    DEALLOCATE (box, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed')
  END SUBROUTINE octree_insert

  !> Insert a single object into the octree data structure.
  !
  !> Range query in octree data structure.
  INTEGER FUNCTION octree_query_range(octree, pmin, pmax)
    TYPE(t_range_octree), INTENT(INOUT) :: octree !< octree data structure
    REAL(wp), INTENT(IN)    :: pmin(3), pmax(3) !< range (corners)
    ! local variables
    INTEGER  :: ibox, l1, i, c
    LOGICAL  :: valid
    REAL(wp) :: brange(2, 3), mid(3) ! box range (min/max, dim=1,2,3)

    octree_query_range = -1
    ibox = 1 ! root box
    brange(:, :) = octree%brange0(:, :)
    ! move down through the octree until the object does not fit into
    ! any of the child boxes:
    INSERT_LOOP: DO l1 = 1, OCTREE_DEPTH
      IF (l1 == OCTREE_DEPTH) octree_query_range = ibox
      ! find containing child box and move down the tree:
      mid = (/((brange(1, i) + brange(2, i))/2._wp, i=1, 3)/)
      c = 0
      valid = .TRUE.
      DO i = 1, 3
        IF (pmin(i) >= mid(i)) THEN
          brange(1, i) = mid(i)
          c = c + ISHFT(1, i - 1)
        ELSE IF (pmax(i) <= mid(i)) THEN
          brange(2, i) = mid(i)
        ELSE
          valid = .FALSE.
        END IF
      END DO
      IF (.NOT. valid) THEN
        octree_query_range = ibox
        EXIT INSERT_LOOP
      ELSE
        ibox = 8*ibox - 6 + c
      END IF
    END DO INSERT_LOOP
  END FUNCTION octree_query_range

  !> searches in octree%box for start index corresponding to "ibox"
  INTEGER FUNCTION octree_get_box_offset(octree, ibox)
    TYPE(t_range_octree), INTENT(IN) :: octree !< octree data structure
    INTEGER, INTENT(IN) :: ibox !< octree box index
    ! local variables
    INTEGER :: res, s0, s1, mid

    ! octree%box is sorted by the first dimension
    !
    ! therefore we can perform a binary search
    s0 = 1
    s1 = SIZE(octree%box, 1)
    res = -1
    LOOP2: DO
      IF (s0 > s1) EXIT LOOP2
      mid = (s1 + s0)/2
      IF (octree%box(mid, 1) > ibox) THEN
        s1 = mid - 1
      ELSE IF (octree%box(mid, 1) < ibox) THEN
        s0 = mid + 1
      ELSE IF (octree%box(mid, 1) == ibox) THEN
        res = mid
        EXIT LOOP2
      END IF
    END DO LOOP2
    octree_get_box_offset = res
  END FUNCTION octree_get_box_offset

  !> Point query in octree data structure.
  !
  !  @return no. of objects contained in traversed octree boxes.
  INTEGER FUNCTION octree_count_point(octree, p)
    TYPE(t_range_octree), INTENT(IN) :: octree !< octree data structure
    REAL(wp), INTENT(IN) :: p(3) !< point to insert
    ! local variables
    INTEGER  :: ibox, l1, i, c, nobjects, ioffset
    REAL(wp) :: brange(2, 3), mid !< box range (min/max, dim=1,2,3)
    REAL(wp) :: pmin(3), pmax(3)

    nobjects = 0
    ibox = 1 ! root box
    brange(:, :) = octree%brange0(:, :)
    LOOP: DO l1 = 1, OCTREE_DEPTH
      ioffset = octree_get_box_offset(octree, ibox)
      !NEC_RP: prevent processes from accessing not allocated data
      IF ((ioffset /= -1) .AND. (SIZE(octree%object_idx) > 0)) THEN
        DO i = octree%box(ioffset, 2), (octree%box(ioffset + 1, 2) - 1)
          pmin = octree%obj_min(octree%object_idx(i), :)
          pmax = octree%obj_max(octree%object_idx(i), :)
          IF ((p(1) >= pmin(1)) .AND. (p(1) <= pmax(1)) .AND. &
            & (p(2) >= pmin(2)) .AND. (p(2) <= pmax(2)) .AND. &
            & (p(3) >= pmin(3)) .AND. (p(3) <= pmax(3))) THEN
            nobjects = nobjects + 1
          END IF
        END DO
      END IF
      ! find containing child box and move down the tree:
      c = 0
      DO i = 1, 3
        mid = (brange(1, i) + brange(2, i))*0.5_wp
        IF (p(i) > mid) THEN
          brange(1, i) = mid
          c = c + ISHFT(1, i - 1)
        ELSE
          brange(2, i) = mid
        END IF
      END DO
      ibox = 8*ibox - 6 + c
    END DO LOOP
    octree_count_point = nobjects
  END FUNCTION octree_count_point

  !> Point query in octree data structure.
  !
  !  @return list of objects contained in traversed octree boxes.
  SUBROUTINE octree_query_point(octree, p, obj_list)
    TYPE(t_range_octree), INTENT(IN)    :: octree !< octree data structure
    REAL(wp), INTENT(IN)    :: p(3) !< point to insert
    INTEGER, INTENT(INOUT) :: obj_list(:) !< result: list of objects in traversed boxes.
    ! local variables
    INTEGER  :: ibox, l1, i, c, nobjects, ioffset
    REAL(wp) :: brange(2, 3), mid !< box range (min/max, dim=1,2,3)
    REAL(wp) :: pmin(3), pmax(3)

    nobjects = 0
    ibox = 1 ! root box
    brange(:, :) = octree%brange0(:, :)
    LOOP: DO l1 = 1, OCTREE_DEPTH
      ioffset = octree_get_box_offset(octree, ibox)
      !NEC_RP: prevent processes from accessing not allocated data
      IF ((ioffset /= -1) .AND. (SIZE(octree%object_idx) > 0)) THEN
        DO i = octree%box(ioffset, 2), (octree%box(ioffset + 1, 2) - 1)
          pmin = octree%obj_min(octree%object_idx(i), :)
          pmax = octree%obj_max(octree%object_idx(i), :)
          IF (((p(1) >= pmin(1)) .AND. (p(1) <= pmax(1)) .AND. &
            &  (p(2) >= pmin(2)) .AND. (p(2) <= pmax(2)) .AND. &
            &  (p(3) >= pmin(3)) .AND. (p(3) <= pmax(3)))) THEN
            nobjects = nobjects + 1
            obj_list(nobjects) = octree%object_data(octree%object_idx(i))
          END IF
        END DO
      END IF
      ! find containing child box and move down the tree:
      c = 0
      DO i = 1, 3
        mid = (brange(1, i) + brange(2, i))*0.5_wp
        IF (p(i) > mid) THEN
          brange(1, i) = mid
          c = c + ISHFT(1, i - 1)
        ELSE
          brange(2, i) = mid
        END IF
      END DO
      ibox = 8*ibox - 6 + c
    END DO LOOP
  END SUBROUTINE octree_query_point

END MODULE mo_octree
