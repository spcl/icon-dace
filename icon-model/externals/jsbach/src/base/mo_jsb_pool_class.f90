!> Contains class definitions for pools
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_pool_class
#ifndef __NO_JSBACH__

  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message
  USE mo_util,          ONLY: int2string
  USE mo_jsb_var_class, ONLY: t_jsb_var, t_jsb_var_p, t_jsb_var_real2d, t_jsb_var_real3d
#ifdef __QUINCY_STANDALONE__
#else
  USE mo_jsb_subset,    ONLY: t_subset, ON_DOMAIN, ON_CHUNK
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_jsb_pool, t_jsb_pool_p
  PUBLIC :: ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15  ! element_id (element enumerator) for pools

  TYPE, ABSTRACT :: t_jsb_pool
    CHARACTER(LEN=:),  ALLOCATABLE :: shortname                                ! No spaces!
    CHARACTER(LEN=:),  ALLOCATABLE :: name                                     ! No spaces!
    CHARACTER(LEN=:),  ALLOCATABLE :: path                                     !
    CHARACTER(LEN=:),  ALLOCATABLE :: element_unit                             ! unit of pool element variables (used for output stream)
    INTEGER                        :: id                            = -1
      !< a pool id that - if set /= -1 - should be unique (from mo_lnd_bgcm_store_class)
    INTEGER                        :: owner_model_id                = 0        ! ID of the jsbach model this pool is running at
    LOGICAL                        :: contains_elements             = .FALSE.  ! contains element variables
    LOGICAL                        :: contains_bgcm                 = .FALSE.  ! contains bgc_material (t_jsb_pool)
    LOGICAL                        :: contains_elem_1d_ptr_to_chunk = .FALSE.  ! t_jsb_pool_elem2d -> set TRUE in the pool init routines in the PROC_memory_class
    LOGICAL                        :: contains_elem_2d_ptr_to_chunk = .FALSE.  ! t_jsb_pool_elem3d -> set TRUE in the pool init routines in the PROC_memory_class
    CLASS(t_jsb_pool_p),   POINTER :: pool_list(:) => NULL()                   ! children
    CLASS(t_jsb_pool_p),   POINTER :: leaf_list(:) => NULL()                   ! leaves
    TYPE(t_jsb_var_p), ALLOCATABLE :: element_list(:)                          ! elements of this pool (elements used in PROC are defined in 'all_elements(:)' in PROC_memory_class)
    TYPE(t_jsb_var_p), ALLOCATABLE :: var_list(:)                              ! other associated variables - currently not used
  CONTAINS
    PROCEDURE(init_iface), DEFERRED :: Init
    PROCEDURE(get_el_name_iface_id), DEFERRED :: Get_element_name_by_id
    GENERIC                         :: Get_element_name => Get_element_name_by_id
    PROCEDURE                       :: Add_element           => t_jsb_pool_add_element
    PROCEDURE                       :: Add_var               => t_jsb_pool_add_var
    PROCEDURE                       :: Has_elements          => t_jsb_pool_has_elements
    PROCEDURE                       :: Add_pool              => t_jsb_pool_add_pool
    PROCEDURE                       ::                          t_jsb_pool_find_element_position_by_id
    PROCEDURE                       ::                          t_jsb_pool_find_element_position_by_name
    PROCEDURE                       :: Set_paths             => t_jsb_pool_set_paths
    GENERIC                         :: Find_element_position => t_jsb_pool_find_element_position_by_id,   &
      &                                                         t_jsb_pool_find_element_position_by_name
#ifdef __QUINCY_STANDALONE__
#else
    PROCEDURE                       :: Find_pool             => t_jsb_pool_find_pool
    PROCEDURE                       :: Associate_var_pointers => t_jsb_pool_associate_var_pointers
    PROCEDURE                       :: Sum_leaves            => t_jsb_pool_sum_leaves
    PROCEDURE                       :: Print                 => t_jsb_pool_print
    PROCEDURE                       :: Get_subset            => t_jsb_pool_get_subset         ! get chunk and block and type information for the pool and thread
#endif
    PROCEDURE                       :: force_finalization    => t_jsb_pool_force_finalization
  END TYPE t_jsb_pool

  TYPE :: t_jsb_pool_p
    CLASS(t_jsb_pool), POINTER :: p
  END TYPE t_jsb_pool_p

  ABSTRACT INTERFACE
    SUBROUTINE init_iface(this, name, shortname, elements_index_map, l_elements, element_list, element_unit)
      IMPORT t_jsb_pool
      CLASS(t_jsb_pool),          INTENT(inout) :: this
      CHARACTER(LEN=*),           INTENT(in)    :: name
      CHARACTER(LEN=*),           INTENT(in)    :: shortname
      INTEGER,                    INTENT(in)    :: elements_index_map(:)
      LOGICAL,          OPTIONAL, INTENT(in)    :: l_elements      ! Add elements yes or no?
      INTEGER,          OPTIONAL, INTENT(in)    :: element_list(:) ! Select elements to add
      CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: element_unit
    END SUBROUTINE init_iface
    FUNCTION get_el_name_iface_id(this, id) RESULT(name)
      IMPORT t_jsb_pool
      CLASS(t_jsb_pool),          INTENT(in) :: this
      INTEGER,           INTENT(in) :: id
      CHARACTER(LEN=:), ALLOCATABLE :: name
    END FUNCTION get_el_name_iface_id
  END INTERFACE

  !>element_id (element enumerator) for bgc_material
  !>
  !! used in the PROC_memory_class for initialization of pools
  !!
  !! names of element variables used in JSBACH are:
  !! (sequence accordingly to the below enumerator)
  !!    carbon, nitrogen, phosphorus, carbon13, carbon14, nitrogen15
  !!
  !! but see element ID ENUM in mo_lnd_bgcm_class
  !!
  ENUM, BIND(C)
    ENUMERATOR :: ELEM_C = 1, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_pool_class'

CONTAINS

  SUBROUTINE t_jsb_pool_add_element(this, name, shortname, id, dim)

    CLASS(t_jsb_pool), INTENT(inout) :: this
    CHARACTER(LEN=*),  INTENT(in)    :: name
    CHARACTER(LEN=*),  INTENT(in)    :: shortname
    INTEGER,           INTENT(in)    :: id
    INTEGER,           INTENT(in)    :: dim

    CLASS(t_jsb_var),  POINTER     :: element
    TYPE(t_jsb_var_p), ALLOCATABLE :: element_list(:)
    INTEGER :: n_elements, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_add_element'

    IF (dim==2) THEN
      ALLOCATE(t_jsb_var_real2d :: element)
    ELSE IF (dim==3) THEN
      ALLOCATE(t_jsb_var_real3d :: element)
    END IF
    element%element_name      = name
    element%element_shortname = shortname
    element%element_id        = id

    IF (ALLOCATED(this%element_list)) THEN
      n_elements = SIZE(this%element_list)
    ELSE
      n_elements = 0
    END IF

    ALLOCATE(element_list(n_elements+1))
    DO i=1,n_elements
      element_list(i)%p => this%element_list(i)%p
    END DO
    element_list(n_elements+1)%p => element

    CALL MOVE_ALLOC(element_list, this%element_list)

  END SUBROUTINE t_jsb_pool_add_element

  FUNCTION t_jsb_pool_has_elements(this) RESULT(has_elements)

    CLASS(t_jsb_pool), INTENT(in) :: this
    LOGICAL                       :: has_elements

    has_elements = ALLOCATED(this%element_list)

  END FUNCTION t_jsb_pool_has_elements

  FUNCTION t_jsb_pool_find_element_position_by_id(this, element_id) RESULT(pos)

    CLASS(t_jsb_pool), INTENT(in) :: this
    INTEGER,           INTENT(in) :: element_id
    INTEGER                       :: pos

    INTEGER :: i

    pos = 0

    IF (.NOT. this%Has_elements()) RETURN

    i = 0
    DO
      i = i + 1
      IF (.NOT. ASSOCIATED(this%element_list(i)%p)) THEN   ! End of list
        EXIT
      ELSE IF (this%element_list(i)%p%element_id == element_id) THEN
        pos = i
        EXIT
      END IF
    END DO

  END FUNCTION t_jsb_pool_find_element_position_by_id

  FUNCTION t_jsb_pool_find_element_position_by_name(this, element_name) RESULT(pos)

    USE mo_util_string, ONLY: tolower

    CLASS(t_jsb_pool), INTENT(in) :: this
    CHARACTER(LEN=*),  INTENT(in) :: element_name
    INTEGER                       :: pos

    INTEGER :: i

    pos = 0

    IF (.NOT. ALLOCATED(this%element_list)) RETURN

    i = 0
    DO
      i = i + 1
      IF (tolower(this%element_list(i)%p%element_name) == tolower(element_name)) THEN
        pos = i
        EXIT
      END IF
    END DO

  END FUNCTION t_jsb_pool_find_element_position_by_name

  SUBROUTINE t_jsb_pool_add_var(this, var, dim)

    CLASS(t_jsb_pool),        INTENT(inout) :: this
    CLASS(t_jsb_var), TARGET, INTENT(in)    :: var
    INTEGER,                  INTENT(in)    :: dim

    TYPE(t_jsb_var_p), ALLOCATABLE :: var_list(:)
    INTEGER :: n_vars, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_add_var'

    IF (ALLOCATED(this%var_list)) THEN
      n_vars = SIZE(this%var_list)
    ELSE
      n_vars = 0
    END IF

    ALLOCATE(var_list(n_vars+1))
    DO i=1,n_vars
      var_list(i)%p => this%var_list(i)%p
    END DO
    var_list(n_vars+1)%p => var

    CALL MOVE_ALLOC(var_list, this%var_list)

  END SUBROUTINE t_jsb_pool_add_var

  SUBROUTINE t_jsb_pool_add_pool(this, pool, id)

    CLASS(t_jsb_pool),         INTENT(inout) :: this
    CLASS(t_jsb_pool), TARGET, INTENT(inout) :: pool
    INTEGER, OPTIONAL,         INTENT(in)    :: id  !< unique id of this bgcm

    TYPE(t_jsb_pool_p), POINTER :: pool_list(:), leaf_list(:)
    INTEGER :: n_pools, n_leaves, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_add_pool'

    ! if bgc materials should be used within the bgcm store, they need to have a unique id for quick association
    IF( PRESENT(id)) pool%id = id

    IF (ASSOCIATED(this%pool_list)) THEN
      n_pools = SIZE(this%pool_list)
    ELSE
      n_pools = 0
    END IF

    ALLOCATE(pool_list(n_pools+1))
    DO i=1,n_pools
      pool_list(i)%p => this%pool_list(i)%p
    END DO
    pool_list(n_pools+1)%p => pool

    IF (ASSOCIATED(this%pool_list)) DEALLOCATE(this%pool_list)
    this%pool_list => pool_list

    IF (ASSOCIATED(this%leaf_list)) THEN
      n_leaves = SIZE(this%leaf_list)
    ELSE
      n_leaves = 0
    END IF

    IF (ASSOCIATED(pool%leaf_list)) THEN
      ALLOCATE(leaf_list(n_leaves + SIZE(pool%leaf_list)))
    ELSE
      ALLOCATE(leaf_list(n_leaves + 1))
    END IF

    DO i=1,n_leaves
      leaf_list(i)%p => this%leaf_list(i)%p
    END DO

    IF (ASSOCIATED(pool%leaf_list)) THEN
      DO i=1,SIZE(pool%leaf_list)
        leaf_list(n_leaves+i)%p => pool%leaf_list(i)%p
      END DO
    ELSE
      leaf_list(n_leaves+1)%p => pool
    END IF

    IF (ASSOCIATED(this%leaf_list)) DEALLOCATE(this%leaf_list)
    this%leaf_list => leaf_list

  END SUBROUTINE t_jsb_pool_add_pool

  RECURSIVE SUBROUTINE t_jsb_pool_set_paths(this, path)

    CLASS(t_jsb_pool), INTENT(inout) :: this
    CHARACTER(LEN=*),     OPTIONAL   :: path

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_set_paths'

    IF (PRESENT(path)) THEN
      IF (path /= '') THEN
        this%path = path // ':' // this%shortname
      ELSE
        this%path = this%shortname
      END IF
    ELSE
      this%path = this%shortname
    END IF
    ! print*,'DDD ', this%shortname, path, this%path

    IF (ASSOCIATED(this%pool_list)) THEN
      DO i=1,SIZE(this%pool_list)
        CALL this%pool_list(i)%p%Set_paths(this%path)
      END DO
    END IF

  END SUBROUTINE t_jsb_pool_set_paths

  RECURSIVE FUNCTION t_jsb_pool_find_pool(this, path) RESULT(pool)

    CLASS(t_jsb_pool), TARGET, INTENT(in) :: this
    CHARACTER(LEN=*),          INTENT(in) :: path
    CLASS(t_jsb_pool), POINTER            :: pool

    INTEGER :: i, ipos
    CHARACTER(LEN=:), ALLOCATABLE :: path_left, path_right

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_find_pool'

    IF (.NOT. ALLOCATED(this%shortname)) &
      & CALL finish(routine, 'Pool must have a shortname')

    IF (path == '') THEN
      pool => this
      RETURN
    END IF

    ipos = INDEX(path, ':')
    IF (ipos > 1) THEN
      path_left = path(1:ipos-1)
      path_right = path(ipos+1:LEN(path))
    ELSE IF (path /= '') THEN
      path_left = path
      path_right = ''
    ELSE
      path_left = ''
      path_right = ''
    END IF

    IF (ASSOCIATED(this%pool_list)) THEN
      DO i=1,SIZE(this%pool_list)
        IF (path_left == this%pool_list(i)%p%shortname) THEN
          pool => this%pool_list(i)%p%Find_pool(path_right)
          EXIT
        END IF
      END DO
    ELSE IF (path_right /= '') THEN
      CALL finish(routine, 'Pool not found '//this%shortname//' '//path)
    END IF

    IF (.NOT. ASSOCIATED(pool)) &
      CALL finish(routine, 'Pool not found '//this%shortname//' '//path)

  END FUNCTION t_jsb_pool_find_pool

  RECURSIVE SUBROUTINE t_jsb_pool_print(this)

    CLASS(t_jsb_pool), INTENT(in) :: this

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_print'

    IF (.NOT. ALLOCATED(this%path)) &
      & CALL finish(routine, 'Pool must have a unique identifier ... call %Set_paths on top pool container!')

    CALL message('', this%path//' ('//this%name//')')

    IF (ALLOCATED(this%element_list)) THEN
      DO i=1,SIZE(this%element_list)
        CALL message('  element', this%element_list(i)%p%element_name           // &
          &            ', id=' // int2string(this%element_list(i)%p%element_id) // &
          &            ', unit=' // this%element_unit                              &
          & )
      END DO
    END IF

    IF (ASSOCIATED(this%pool_list)) THEN
      DO i=1,SIZE(this%pool_list)
        CALL this%pool_list(i)%p%Print()
      END DO
    END IF

  END SUBROUTINE t_jsb_pool_print

  !================================================================================================================================
  ! get subset
  ! chunk and block information per model & thread from mo_jsb_subset:jsbach_subsets(model_id)%sub(thread)
  !================================================================================================================================
  FUNCTION t_jsb_pool_get_subset(this) RESULT(subset)

    USE mo_jsb_subset,   ONLY: jsbach_subsets
    USE mo_jsb_parallel, ONLY: Get_omp_thread

    CLASS(t_jsb_pool), INTENT(in) :: this
    TYPE(t_subset)                :: subset

    subset = jsbach_subsets(this%owner_model_id)%sub(Get_omp_thread())

  END FUNCTION t_jsb_pool_get_subset

  RECURSIVE SUBROUTINE t_jsb_pool_associate_var_pointers(this, ic_start, ic_end, iblk_start, iblk_end)

    USE mo_jsb_var_class,  ONLY: t_jsb_var_real2d, t_jsb_var_real3d
    USE mo_jsb_subset,     ONLY: ON_DOMAIN, ON_CHUNK

    CLASS(t_jsb_pool),  INTENT(inout) :: this
    INTEGER, OPTIONAL,  INTENT(in)    :: ic_start
    INTEGER, OPTIONAL,  INTENT(in)    :: ic_end
    INTEGER, OPTIONAL,  INTENT(in)    :: iblk_start
    INTEGER, OPTIONAL,  INTENT(in)    :: iblk_end

    INTEGER :: i, itype

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_associate_var_pointers'

    IF (iblk_start == iblk_end) THEN
      itype = ON_CHUNK
    ELSE
      itype = ON_DOMAIN
    END IF

    ! loop over elements in element_list(:)
    !   associate pointer %ptrm1
    !   set subset_type
    !   if applicable: associate REAL element vectors with %ptrm1 & %ptr2d/3d
    IF (ALLOCATED(this%element_list)) THEN
      DO i=1,SIZE(this%element_list)
        ASSOCIATE (element => this%element_list(i)%p)
          CALL element%Associate_pointers(ic_start, ic_end, iblk_start, iblk_end)
          element%subset_type = itype
        END ASSOCIATE
      END DO
    END IF

    IF (ALLOCATED(this%var_list)) THEN
      DO i=1,SIZE(this%var_list)
        ASSOCIATE (var => this%var_list(i)%p)
          CALL var%Associate_pointers(ic_start, ic_end, iblk_start, iblk_end)
          var%subset_type = itype
        END ASSOCIATE
      END DO
    END IF

    IF (ASSOCIATED(this%pool_list)) THEN
      DO i=1,SIZE(this%pool_list)
        CALL this%pool_list(i)%p%Associate_var_pointers(ic_start, ic_end, iblk_start, iblk_end)
      END DO
    END IF

  END SUBROUTINE t_jsb_pool_associate_var_pointers

  SUBROUTINE t_jsb_pool_force_finalization(this)

    CLASS(t_jsb_pool), INTENT(inout) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_force_finalization'

  END SUBROUTINE t_jsb_pool_force_finalization

  SUBROUTINE t_jsb_pool_sum_leaves(this, element_list)

    CLASS(t_jsb_pool), INTENT(inout) :: this
    INTEGER,           INTENT(in)    :: element_list(:)

    INTEGER :: i, i_leaf, i_element

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_pool_sum_leaves'

    IF (.NOT. ASSOCIATED(this%leaf_list)) CALL finish(routine, 'Pool has no leaves')

    DO i=1,SIZE(element_list)
      i_element = this%Find_element_position(element_list(i))
      IF (i_element < 1) CALL finish(routine, 'Element '//this%Get_element_name(element_list(i))//' not found')

      SELECT TYPE (element1 => this%element_list(i_element)%p)
      CLASS IS (t_jsb_var_real2d)
        element1 = 0._wp
        DO i_leaf=1,SIZE(this%leaf_list)
          SELECT TYPE (element2 => this%leaf_list(i_leaf)%p%element_list(i_element)%p)
          CLASS IS (t_jsb_var_real2d)
            element1 = element1 + element2
          END SELECT
        END DO
      CLASS IS (t_jsb_var_real3d)
        element1 = 0._wp
        DO i_leaf=1,SIZE(this%leaf_list)
          SELECT TYPE (element2 => this%leaf_list(i_leaf)%p%element_list(i_element)%p)
          CLASS IS (t_jsb_var_real3d)
            element1 = element1 + element2
          END SELECT
        END DO
      END SELECT
    END DO

  END SUBROUTINE t_jsb_pool_sum_leaves

#endif
END MODULE mo_jsb_pool_class
