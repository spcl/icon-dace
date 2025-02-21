! Basic derived types for blocked and non-blocked 1D index lists.
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

MODULE mo_idx_list

  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_exception,           ONLY: finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_fortran_tools,       ONLY: DO_DEALLOCATE
  USE mo_communication,       ONLY: idx_1d, blk_no, idx_no
  USE mo_sync,                ONLY: global_sum_array



  USE mo_util_sort,           ONLY: quicksort


  IMPLICIT NONE

  PRIVATE

  ! subroutines
  PUBLIC :: copy_list_blocked
  PUBLIC :: compare_sets

  ! types
  PUBLIC :: t_idx_list_blocked
  PUBLIC :: t_idx_list1D

  ! standard (non-blocked) index list
  TYPE t_idx_list1D
    INTEGER, ALLOCATABLE :: idx(:)     ! 1D indices
    INTEGER              :: ncount     ! length
    LOGICAL              :: lopenacc   ! list is copied to GPU if true

  CONTAINS
    !
    ! initialize
    PROCEDURE  :: construct => idx_list1D__construct
    !
    ! finalize
    PROCEDURE  :: finalize => idx_list1D__finalize
    !
    ! get blocked list
    PROCEDURE  :: get_blocked_list => idx_list1D__get_blocked_list
  END TYPE t_idx_list1D


  ! blocked index list
  TYPE t_idx_list_blocked
    INTEGER, ALLOCATABLE :: idx(:,:)  ! cell indices for given block
    INTEGER, ALLOCATABLE :: ncount(:) ! length for given block
    LOGICAL              :: lopenacc  ! list is copied to GPU if true

  CONTAINS
    !
    ! initialize
    PROCEDURE  :: construct => idx_list_blocked__construct
    !
    ! finalize
    PROCEDURE  :: finalize => idx_list_blocked__finalize
    !
    ! get 1D list
    PROCEDURE  :: get_list1D => idx_list_blocked__get_list1D
    !
    ! get global sum (total number of of points in index list)
    PROCEDURE  :: get_sum_global => idx_list_blocked__get_sum_global
  END TYPE t_idx_list_blocked




CONTAINS


  !-------------------------------------------------------------------------
  !! Allocates all components of the object of type t_idx_list1D
  !!
  SUBROUTINE idx_list1D__construct(obj, size, lopenacc)
    CLASS(t_idx_list1D) :: obj
    INTEGER, INTENT(IN) :: size
    LOGICAL, INTENT(IN), OPTIONAL :: lopenacc
    !
    ! local
    INTEGER :: ist

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_idx_list: idx_list1D__construct'
  !----------------------------------------------------------

    IF (PRESENT(lopenacc)) THEN
      obj%lopenacc = .TRUE.
    ELSE
      obj%lopenacc = .FALSE.
    ENDIF

    ALLOCATE(obj%idx(size), STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ( TRIM(routine), 'allocation for idx failed' )
    ENDIF

    ! initialize
    obj%idx(:) = -1
    obj%ncount = 0

    !$ACC ENTER DATA COPYIN(obj%idx, obj%ncount) IF(obj%lopenacc)

  END SUBROUTINE idx_list1D__construct



  !-------------------------------------------------------------------------
  !! Deallocates all components of the object of type t_idx_list1D
  !!
  SUBROUTINE idx_list1D__finalize(obj)
    CLASS(t_idx_list1D) :: obj

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_idx_list: idx_list1D__finalize'
  !----------------------------------------------------------

    CALL DO_DEALLOCATE(obj%idx)
    obj%ncount = 0

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(obj%idx, obj%ncount) IF(obj%lopenacc)

  END SUBROUTINE idx_list1D__finalize




  !-------------------------------------------------------------------------
  !! Convert non-blocked 1D index list into blocked list
  !!
  SUBROUTINE idx_list1D__get_blocked_list(obj, list1D_blocked)
    CLASS(t_idx_list1D) :: obj
    TYPE(t_idx_list_blocked), INTENT(INOUT) :: list1D_blocked
    !
    ! local
    INTEGER :: jc, jb, ic

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_idx_list: get_blocked_list'
  !----------------------------------------------------------

    DO ic = 1, obj%ncount
      jb = blk_no(obj%idx(ic))
      jc = idx_no(obj%idx(ic))
      list1D_blocked%ncount(jb) = list1D_blocked%ncount(jb)+1
      list1D_blocked%idx(list1D_blocked%ncount(jb),jb) = jc
    ENDDO

  END SUBROUTINE idx_list1D__get_blocked_list


  !-------------------------------------------------------------------------
  !! Allocates all components of the object of type t_idx_list_blocked
  !!
  SUBROUTINE idx_list_blocked__construct(obj, nproma, nblks, lopenacc)
    CLASS(t_idx_list_blocked) :: obj
    INTEGER, INTENT(IN) :: nproma
    INTEGER, INTENT(IN) :: nblks
    LOGICAL, INTENT(IN), OPTIONAL :: lopenacc
    !
    ! local
    INTEGER :: ist

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_idx_list: construct'
  !----------------------------------------------------------

    IF (PRESENT(lopenacc)) THEN
      obj%lopenacc = .TRUE.
    ELSE
      obj%lopenacc = .FALSE.
    ENDIF

    ALLOCATE(obj%idx(nproma,nblks), &
      &      obj%ncount(nblks), STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish ( TRIM(routine), 'allocation for idx and ncount failed' )
    ENDIF

    ! initialize
    obj%idx(:,:)  = -1
    obj%ncount(:) = 0

    !$ACC ENTER DATA COPYIN(obj%idx, obj%ncount) IF(obj%lopenacc)

  END SUBROUTINE idx_list_blocked__construct



  !-------------------------------------------------------------------------
  !! Deallocates all components of the object of type t_idx_list_blocked
  !!
  SUBROUTINE idx_list_blocked__finalize(obj)
    CLASS(t_idx_list_blocked) :: obj
    !
    ! local
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_idx_list: finalize'
  !----------------------------------------------------------

    CALL DO_DEALLOCATE(obj%idx)
    CALL DO_DEALLOCATE(obj%ncount)

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(obj%idx, obj%ncount) IF(obj%lopenacc)

  END SUBROUTINE idx_list_blocked__finalize




  !-------------------------------------------------------------------------
  !! Get non-blocked 1D list from blocked list of type t_idx_list_blocked
  !!
  SUBROUTINE idx_list_blocked__get_list1D(obj, list1D)
    CLASS(t_idx_list_blocked)          :: obj
    TYPE(t_idx_list1D), INTENT(INOUT)  :: list1D    ! 1D (non-blocked) index list
    !
    ! local
    INTEGER :: jb, jc, ic
    INTEGER :: cnt       ! counter

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_idx_list: get_list1D'
  !----------------------------------------------------------
    cnt = 0
    DO jb=1,SIZE(obj%ncount)
      !     
      IF (obj%ncount(jb) == 0) CYCLE
      !
      DO ic=1,obj%ncount(jb)
        cnt = cnt + 1
        jc = obj%idx(ic,jb) 
        list1D%idx(cnt) = idx_1d(jc,jb)
      ENDDO 
    ENDDO
    list1D%ncount = cnt
  END SUBROUTINE idx_list_blocked__get_list1D



  !-------------------------------------------------------------------------
  !! Get global sum of points in index list
  !!
  FUNCTION idx_list_blocked__get_sum_global(obj, opt_startblk, opt_endblk) RESULT(ncount_global)
    CLASS(t_idx_list_blocked)     :: obj
    INTEGER, OPTIONAL, INTENT(IN) :: opt_startblk, opt_endblk ! start and end block for summation
    INTEGER                       :: ncount_global            ! total number of list entries
    !
    ! local
    INTEGER :: i_startblk, i_endblk

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_idx_list: get_sum_global'
  !----------------------------------------------------------

    IF (PRESENT(opt_startblk)) THEN
      i_startblk = opt_startblk
    ELSE
      i_startblk = 1
    ENDIF

    IF (PRESENT(opt_endblk)) THEN
      i_endblk = opt_endblk
    ELSE
      i_endblk = SIZE(obj%ncount)
    ENDIF

    ncount_global = SUM(obj%ncount(i_startblk:i_endblk))
    ncount_global = global_sum_array(ncount_global)

  END FUNCTION idx_list_blocked__get_sum_global


  !-------------------------------------------------------------------------
  !! Copies traditional blocked index list (consisting of 2 integer arrays) 
  !! to blocked list of type t_idx_list_blocked
  !!
  SUBROUTINE copy_list_blocked(source_list, source_count, target_list)

    INTEGER,                  INTENT(IN   ) :: source_list(:,:)
    INTEGER,                  INTENT(IN   ) :: source_count(:)
    TYPE(t_idx_list_blocked), INTENT(INOUT) :: target_list
    !
    ! local
    INTEGER :: jb, ic
  !----------------------------------------------------------

    DO jb=1,SIZE(source_count)
      target_list%ncount(jb) = source_count(jb)

      IF (source_count(jb)==0) CYCLE

      DO ic=1,source_count(jb)
        target_list%idx(ic,jb) = source_list(ic,jb)
      ENDDO
    ENDDO
  END SUBROUTINE copy_list_blocked



  !-------------------------------------------------------------------------
  !! Compares 2 sorted lists of type t_idx_list1D and groups  
  !! the elements into the following three sublists 
  !! list_intersect: elements which exist in both lists
  !! list1_only    : elements which exist in list 1 only
  !! list2_only    : elements which exist in list 2 only
  !!
  !! Output: sublists of type t_idx_list1D (i.e. non-blocked)
  !!
  SUBROUTINE compare_sorted_1Dlists(list1, list2, list_intersect, &
    &                               list1_only, list2_only)

    TYPE(t_idx_list1D), TARGET, INTENT(IN   ) :: list1, list2    ! lists to be compared
    TYPE(t_idx_list1D),         INTENT(INOUT) :: list_intersect  ! elements found in both lists
    TYPE(t_idx_list1D),         INTENT(INOUT) :: list1_only      ! elements found in list 1 only
    TYPE(t_idx_list1D),         INTENT(INOUT) :: list2_only      ! elements found in list 2 only

    !
    ! local
    INTEGER :: ic
    INTEGER :: cnt_l1, cnt_l2, cnt_lboth, cnt_1only, cnt_2only
    LOGICAL :: fall_off_list1, fall_off_list2                   ! .TRUE.: we fell off the list
    INTEGER, POINTER :: element_l1, element_l2                  ! pointer to list elements

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_idx_list: compare_sorted_1Dlists'
  !----------------------------------------------------------

    ! init
    element_l1 => NULL()
    element_l2 => NULL()

    cnt_l1    = 1
    cnt_l2    = 1
    cnt_lboth = 0
    cnt_1only = 0
    cnt_2only = 0

    fall_off_list1 = cnt_l1 > list1%ncount
    fall_off_list2 = cnt_l2 > list2%ncount


    DO WHILE ((.NOT.fall_off_list1) .AND. (.NOT.fall_off_list2))
      element_l1 => list1%idx(cnt_l1)
      element_l2 => list2%idx(cnt_l2)

      IF (element_l1 == element_l2) THEN
        ! element is found in both lists
        cnt_lboth = cnt_lboth + 1
        list_intersect%idx(cnt_lboth) = element_l1
        cnt_l1 = cnt_l1 + 1
        cnt_l2 = cnt_l2 + 1
      ELSE IF (element_l1 < element_l2) THEN
        cnt_1only = cnt_1only + 1
        list1_only%idx(cnt_1only) = element_l1
        cnt_l1 = cnt_l1 + 1
      ELSE  ! element_l1 > element_l2
        cnt_2only = cnt_2only + 1
        list2_only%idx(cnt_2only) = element_l2
        cnt_l2 = cnt_l2 + 1
      ENDIF
      ! check whether we reached the end of one of the lists
      fall_off_list1 = cnt_l1 > list1%ncount
      fall_off_list2 = cnt_l2 > list2%ncount
    END DO


    ! we reached the end of at least one list. 
    ! Check if there are exist remaining elements in the lists 
    ! and group them.
    IF (fall_off_list1 .AND. fall_off_list2) THEN
      ! nothing to do
    ELSE IF (fall_off_list1) THEN
      ! we fell off list 1, but not list 2
      ! => remaining elements in list2 appear only in list2
      ! => add them to list2_only
      DO ic = cnt_l2,list2%ncount
        cnt_2only = cnt_2only + 1
        list2_only%idx(cnt_2only) = list2%idx(ic)
      ENDDO
    ELSE IF (fall_off_list2) THEN
      ! we fell off list 2, but not list 1
      ! => remaining elements in list1 appear only in list1
      ! => add them to list1_only
      DO ic = cnt_l1,list1%ncount
        cnt_1only = cnt_1only + 1
        list1_only%idx(cnt_1only) = list1%idx(ic)
      ENDDO
    ELSE
      CALL finish (routine, 'Exited while-loop even though end of lists was not reached')
    ENDIF

    list_intersect%ncount = cnt_lboth
    list1_only%ncount     = cnt_1only
    list2_only%ncount     = cnt_2only

  END SUBROUTINE compare_sorted_1Dlists


  !-------------------------------------------------------------------------
  !! Wrapper routine for compare_sorted_1Dlists, which reads and writes 
  !! blocked lists and performs sorting.
  !!
  !! Compares 2 blocked lists of type t_idx_list_blocked and groups  
  !! the elements into the following three blocked sublists 
  !! list_intersect: elements which exist in both lists
  !! list1_only    : elements which exist in list 1 only
  !! list2_only    : elements which exist in list 2 only
  !!
  !! The basic concept is as follows:
  !! 1) convert blocked lists into unblocked 1D lists
  !! 2) sort lists using the quicksort/radixsort algorithm
  !! 3) compare the two lists and group the elements into 3 sublists
  !! 4) transform the sublists back into blocked lists
  !!
  SUBROUTINE compare_sets(p_patch, list1, list2, list_intersect, &
    &                     list1_only, list2_only)

    TYPE(t_patch),                    INTENT(IN   ) :: p_patch
    TYPE(t_idx_list_blocked), TARGET, INTENT(IN   ) :: list1, list2    ! lists to be compared
    TYPE(t_idx_list_blocked),         INTENT(INOUT) :: list_intersect  ! elements found in both lists
    TYPE(t_idx_list_blocked),         INTENT(INOUT) :: list1_only      ! elements found in list 1 only
    TYPE(t_idx_list_blocked),         INTENT(INOUT) :: list2_only      ! elements found in list 2 only

    ! local
    TYPE(t_idx_list1D) :: list1_1d, list2_1d
    TYPE(t_idx_list1D) :: list_intersect_1d, list1_only_1d, list2_only_1d
  !----------------------------------------------------------

    ! transform blocked lists into unblocked 1D lists
    ! 
    CALL list1_1d%construct(p_patch%n_patch_cells)
    CALL list2_1d%construct(p_patch%n_patch_cells)
    !
    CALL list1%get_list1D(list1_1d)
    CALL list2%get_list1D(list2_1d)


    ! sort




    CALL quicksort(a=list1_1d%idx, l_in=1, r_in=list1_1d%ncount)
    CALL quicksort(a=list2_1d%idx, l_in=1, r_in=list2_1d%ncount)



    ! compare sorted 1D lists
    ! elements are grouped into three 1D sublists
    CALL list_intersect_1d%construct(p_patch%n_patch_cells)
    CALL list1_only_1d%construct    (p_patch%n_patch_cells)
    CALL list2_only_1d%construct    (p_patch%n_patch_cells)
    !
    CALL compare_sorted_1Dlists( list1          = list1_1d,          & !in
             &                   list2          = list2_1d,          & !in
             &                   list_intersect = list_intersect_1d, & !inout
             &                   list1_only     = list1_only_1d,     & !inout
             &                   list2_only     = list2_only_1d      ) !inout

    ! transform non-blocked sublists into blocked sublists
    !
    CALL list_intersect_1d%get_blocked_list(list_intersect)
    CALL list1_only_1d%get_blocked_list(list1_only)
    CALL list2_only_1d%get_blocked_list(list2_only)


    ! cleanup
    CALL list1_1d%finalize()
    CALL list2_1d%finalize()
    !
    CALL list_intersect_1d%finalize()
    CALL list1_only_1d%finalize()
    CALL list2_only_1d%finalize()

  END SUBROUTINE compare_sets


END MODULE mo_idx_list

