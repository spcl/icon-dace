!> Definition of Hierarchical State Machine.
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
MODULE mo_hsm_class

  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_exception,          ONLY: message, message_text, finish, warning
  ! USE mo_util,      ONLY: int2string
#ifdef _OPENMP
  USE omp_lib,      ONLY: omp_get_max_threads, omp_get_thread_num
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_Message, t_Action, t_State, t_Hsm

  INTEGER, PARAMETER :: MAX_STATE_NESTING = 10
  INTEGER, PARAMETER :: MAX_ACTIONS = 5

  ENUM, BIND(C)
    ENUMERATOR :: &
      & HSM_full = 1, &
      & HSM_simple
  END ENUM

  TYPE t_Action
    CHARACTER(len=SHORT_NAME_LEN) :: name = ""
  CONTAINS
    PROCEDURE :: Is_equal        => Is_equal_action
    GENERIC   :: OPERATOR(==)    => Is_equal
    PROCEDURE :: Is_not_equal    => Is_not_equal_action
    GENERIC   :: OPERATOR(/=)    => Is_not_equal
  END TYPE t_Action

  INTERFACE t_Action
    PROCEDURE Construct_action
  END INTERFACE t_Action

  TYPE t_Message
    CHARACTER(len=50)      :: name = ''
    TYPE(t_Action)         :: action
  END TYPE t_Message

  INTERFACE t_Message
    PROCEDURE Construct_message
  END INTERFACE t_Message
  !>
  !! Type for States of the State Machine
  !!
  !! This derived type defines the state class
  !!
  TYPE, ABSTRACT :: t_State
    CHARACTER(len=SHORT_NAME_LEN) :: name
    CHARACTER(len=50)             :: description
!!$    INTEGER                     :: id
    INTEGER                       :: level = 0                      ! Level in hierarchy (top=1)
    INTEGER                       :: path(MAX_STATE_NESTING) = -1
    INTEGER, PRIVATE              :: no_of_children = -1
    !
    LOGICAL,        ALLOCATABLE   :: visited(:)
    CLASS(t_State), POINTER       :: parent       => NULL()
    CLASS(t_State), POINTER       :: first_child  => NULL()
    CLASS(t_State), POINTER       :: last_child   => NULL()
    !
    CLASS(t_State), POINTER       :: next_sibling => NULL()
    CLASS(t_State), POINTER       :: prev_sibling => NULL()
  CONTAINS
    PROCEDURE :: Init_state
    PROCEDURE :: Get_parent         => Get_parent_state
    PROCEDURE :: Get_parent_name    => Get_parent_name_state
    PROCEDURE :: Set_no_of_children => Set_no_of_children_state
    PROCEDURE :: Get_no_of_children => Get_no_of_children_state
    PROCEDURE :: Get_pos_of_child   => Get_pos_of_child_state
    PROCEDURE :: Has_children       => Has_children_state
    PROCEDURE :: Has_siblings       => Has_siblings_state
    PROCEDURE :: Is_first_child     => Is_first_child_state
    PROCEDURE :: Is_last_child      => Is_last_child_state
    PROCEDURE :: Is_last_leaf       => Is_last_leaf_state
    PROCEDURE :: Reset_all          => Reset_all_state
    PROCEDURE :: Print_state
    PROCEDURE :: Is_equal           => Is_equal_state
    GENERIC   :: OPERATOR(==)       => Is_equal
    PROCEDURE :: Is_not_equal       => Is_not_equal_state
    GENERIC   :: OPERATOR(/=)       => Is_not_equal
    PROCEDURE :: Process_message    => Process_message_state
    PROCEDURE(Handler_interface), DEFERRED :: Handler             ! State's handler function
  END TYPE t_State

  ! Derived-type constructor
!!$  INTERFACE t_State
!!$    PROCEDURE Construct_state
!!$  END INTERFACE t_State

  TYPE t_State_p
    CLASS(t_State), POINTER :: p
  END TYPE t_State_p
  !>
  !! Type for Hierarchical State Machine
  !!
  !! This derived type defines the structure of a basic Hierarchical State Machine
  !!
  TYPE, ABSTRACT :: t_Hsm
    CHARACTER(len=30)            :: name
    CLASS(t_State), POINTER      :: top     => NULL()    !< Topmost state
    !
    TYPE(t_State_p), ALLOCATABLE :: current(:)     !< Current state
    TYPE(t_State_p), ALLOCATABLE :: next(:)        !< Next state (non_null if transition taken)
    TYPE(t_State_p), ALLOCATABLE :: source(:)      !< Source state during last transition
    TYPE(t_Action),  ALLOCATABLE :: actions(:)         !< List of possible actions (events)
    INTEGER,         ALLOCATABLE :: levels_to_lca(:)   !< Temporary storage
    INTEGER                      :: mode               !< HSM_full or HSM_simple
  CONTAINS
    PROCEDURE :: Init_hsm
    PROCEDURE :: Set_mode          => Set_mode_hsm
    PROCEDURE :: Get_current       => Get_current_hsm
    PROCEDURE :: Goto              => Goto_state_hsm
    PROCEDURE :: Goto_next         => Goto_next_state_hsm
    PROCEDURE :: Get_next          => Get_next_state_hsm
    PROCEDURE :: Get_state         => Get_state_hsm
    PROCEDURE :: Start_state       => Start_state_hsm
    PROCEDURE :: Reset_state       => Reset_state_hsm
    PROCEDURE :: Reset_all         => Reset_all_hsm
    PROCEDURE :: Register_action   => Register_action_hsm
    PROCEDURE :: Get_action        => Get_action_hsm
    PROCEDURE :: Print_actions     => Print_actions_hsm
    PROCEDURE :: Process_message   => Process_message_hsm        ! State machine engine
    PROCEDURE :: Process_entry_to_next
    PROCEDURE :: Get_levels_to_LCA => Get_levels_to_LCA_hsm
    PROCEDURE :: Process_exit_to_LCA
    PROCEDURE :: Start             => Start_hsm                  ! Start state machine
  END TYPE t_Hsm

  ! Derived-type constructor
!!$  INTERFACE t_Hsm
!!$    PROCEDURE Construct_hsm
!!$  END INTERFACE t_Hsm

!!$  ABSTRACT INTERFACE
!!$    FUNCTION Handler_interface(this, hsm, msg_in) RESULT(msg_out)
!!$      IMPORT :: t_State, t_Hsm, t_Message
!!$      CLASS(t_State),   INTENT(inout) :: this
!!$      CLASS(t_Hsm),     POINTER       :: hsm
!!$      CLASS(t_Message), POINTER       :: msg_in
!!$      CLASS(t_Message), POINTER       :: msg_out
!!$    END FUNCTION Handler_interface
!!$  END INTERFACE
  ABSTRACT INTERFACE
    FUNCTION Handler_interface(this, msg_in) RESULT(msg_out)
      IMPORT :: t_State, t_Hsm, t_Message
      CLASS(t_State),   INTENT(inout) :: this
      CLASS(t_Message), INTENT(in)    :: msg_in
      CLASS(t_Message), POINTER       :: msg_out
    END FUNCTION Handler_interface
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hsm_class'

CONTAINS

  ! ===============================================================================
  ! Methods for state class
  ! ===============================================================================
!!$  FUNCTION Construct_state(name, description, parent, handler) RESULT(return_ptr)
!!$
!!$    CHARACTER(len=*), INTENT(in)      :: name
!!$    CHARACTER(len=*), INTENT(in)      :: description
!!$    CLASS(t_State), POINTER, OPTIONAL :: parent
!!$    PROCEDURE(Handler_interface), POINTER :: handler
!!$    PROCEDURE(), POINTER :: handler
!!$    CLASS(t_State), POINTER           :: return_ptr
!!$
!!$    INTEGER :: depth
!!$
!!$    ALLOCATE(return_ptr)
!!$
!!$    return_ptr%name = TRIM(name)
!!$    return_ptr%description = TRIM(description)
!!$    return_ptr%Handler => handler
!!$
!!$    IF (.NOT. PRESENT(parent)) THEN
!!$      return_ptr%path(:) = -1
!!$      return_ptr%path(1) = 1
!!$      CALL return_ptr%Print()
!!$      RETURN
!!$    END IF
!!$
!!$    return_ptr%path(:) = parent%path(:)
!!$    depth = MINLOC(parent%path, DIM=1)  ! Index of first -1 in path of parent
!!$
!!$    IF (.NOT. ASSOCIATED(parent%first_child)) THEN
!!$      parent%first_child => return_ptr
!!$      return_ptr%path(depth) = 1
!!$    END IF
!!$
!!$    IF (ASSOCIATED(parent%last_child)) THEN
!!$      parent%last_child%next_sibling => return_ptr
!!$      return_ptr%prev_sibling        => parent%last_child
!!$      return_ptr%path(depth) = return_ptr%prev_sibling%path(depth) + 1
!!$    END IF
!!$
!!$    parent%last_child  => return_ptr
!!$
!!$    return_ptr%parent => parent
!!$
!!$    CALL return_ptr%Print()
!!$
!!$  END FUNCTION Construct_state

!!$  SUBROUTINE Init_state(this, name, description, parent, handler)
  SUBROUTINE Init_state(this, name, description, parent)
    CLASS(t_State),   INTENT(inout), TARGET :: this
    CHARACTER(len=*), INTENT(in)            :: name
    CHARACTER(len=*), INTENT(in)            :: description
    CLASS(t_State), POINTER, OPTIONAL       :: parent
!!$    PROCEDURE(Handler_interface), POINTER   :: handler

    INTEGER :: depth
    INTEGER :: i, nthreads

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif
    ALLOCATE(this%visited(nthreads))
    DO i=1,nthreads
      this%visited(i) = .FALSE.
    END DO


    this%name = TRIM(name)
    this%description = TRIM(description)
!!$    this%Handler => handler

    IF (.NOT. PRESENT(parent)) THEN
      this%path(:) = -1
      this%path(1) = 1
      this%level = 1
      !CALL this%Print_state()
      RETURN
    END IF

    SELECT TYPE(parent)
    CLASS IS (t_State)

      this%path(:) = parent%path(:)
      depth = parent%level + 1
      ! depth = MINLOC(parent%path, DIM=1)  ! Index of first -1 in path of parent
      this%level = depth

      IF (.NOT. ASSOCIATED(parent%first_child)) THEN
        parent%first_child => this
        this%path(depth) = 1
      END IF

      IF (ASSOCIATED(parent%last_child)) THEN
        parent%last_child%next_sibling => this
        this%prev_sibling        => parent%last_child
        this%path(depth) = this%prev_sibling%path(depth) + 1
      END IF

      parent%last_child  => this

      this%parent => parent

    END SELECT

    ! CALL this%Print()

  END SUBROUTINE Init_state

  FUNCTION Get_parent_state(this) RESULT(return_ptr)

    CLASS(t_State), INTENT(in) :: this
    CLASS(*), POINTER    :: return_ptr

    CLASS(*), POINTER :: p

!!$    return_ptr => this%parent
    IF (ASSOCIATED(this%parent)) THEN
      p => this%parent
      SELECT TYPE(p)
      CLASS IS (t_State)
        return_ptr => p
      END SELECT
    ELSE
      return_ptr => NULL()
    END IF

  END FUNCTION Get_parent_state

  FUNCTION Get_parent_name_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in)    :: this
    CHARACTER(len=SHORT_NAME_LEN) :: return_value

    CLASS(*), POINTER :: p

    p => this%parent
    IF (ASSOCIATED(p)) THEN
      SELECT TYPE(p)
      CLASS IS (t_State)
        return_value = p%name
      END SELECT
    ELSE
      return_value = ''
    END IF

  END FUNCTION Get_parent_name_state

  FUNCTION Has_children_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                   :: return_value

    return_value = ASSOCIATED(this%first_child)

  END FUNCTION Has_children_state

  FUNCTION Has_siblings_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                    :: return_value

    return_value = ASSOCIATED(this%next_sibling) .OR. ASSOCIATED(this%prev_sibling)

  END FUNCTION Has_siblings_state

  FUNCTION Is_last_child_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                   :: return_value

    return_value = .NOT. ASSOCIATED(this%next_sibling)

  END FUNCTION Is_last_child_state

  ! TODO This only works if the tree is walked through downwards!
  FUNCTION Is_last_leaf_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                   :: return_value

    IF (ASSOCIATED(this%parent)) THEN
      return_value = .NOT. this%Has_children() .AND. this%Is_last_child() .AND. this%parent%Is_last_child()
    ELSE
      return_value = .NOT. this%Has_children() .AND. this%Is_last_child()
    END IF

  END FUNCTION Is_last_leaf_state

  FUNCTION Is_first_child_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                   :: return_value

    return_value = .NOT. ASSOCIATED(this%prev_sibling)

  END FUNCTION Is_first_child_state

  SUBROUTINE Set_no_of_children_state(this)

    CLASS(t_State), INTENT(inout) :: this

    CLASS(t_State), POINTER :: current

    IF (this%no_of_children >= 0) THEN
      RETURN
    END IF

    IF (this%Has_children()) THEN
      this%no_of_children = 1
      current => this%first_child
      DO WHILE (.NOT. current%Is_last_child())
        current => current%next_sibling
        this%no_of_children = this%no_of_children + 1
      END DO
    ELSE
      this%no_of_children = 0
    END IF

  END SUBROUTINE Set_no_of_children_state

  FUNCTION Get_no_of_children_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(inout) :: this
    INTEGER                       :: return_value

    IF (this%no_of_children < 0) THEN
      CALL this%Set_no_of_children()
    END IF

    return_value = this%no_of_children

  END FUNCTION Get_no_of_children_state

  FUNCTION Get_pos_of_child_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    INTEGER                    :: return_value

    INTEGER :: depth

    IF (ALL(this%path > 0)) THEN
      depth = SIZE(this%path)
    ELSE
      depth = MINLOC(this%path, DIM=1) - 1
    END IF

    return_value = this%path(depth)

  END FUNCTION

  ! Resets %visited(:) to .FALSE. for state and all states below
  ! So, calling this on the top state resets all states
  RECURSIVE SUBROUTINE Reset_all_state(this)

    CLASS(t_State), INTENT(inout) :: this

    INTEGER :: nthreads

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    this%visited(1:nthreads) = .FALSE.

    IF (this%Has_children()) THEN
      CALL this%first_child%Reset_all()
    END IF
    IF (ASSOCIATED(this%next_sibling)) THEN
      CALL this%next_sibling%Reset_all()
    END IF

  END SUBROUTINE Reset_all_state

  SUBROUTINE Print_state(this)

    CLASS(t_State), INTENT(in) :: this

    CHARACTER(LEN=3*MAX_STATE_NESTING) :: c_path

    INTEGER :: i, pos1, pos2, pos_sum

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    c_path = ''
    pos1 = 0
    pos2 = 0
    pos_sum = 0
    DO i=1,MAX_STATE_NESTING
      IF (this%path(i) /= -1) THEN
        pos1 = pos2 + 1
        IF (this%path(i) < 10) THEN
          pos2 = pos1 + 1
          WRITE(c_path(pos1:pos2),'(I1,A1)') this%path(i), '-'
        ELSE IF (this%path(i) < 100) THEN
          pos2 = pos1 + 2
          WRITE(c_path(pos1:pos2),'(I2,A1)') this%path(i), '-'
        END IF
        pos_sum = pos_sum + pos2 - pos1 + 1
      ELSE
        c_path(pos_sum:3*MAX_STATE_NESTING) = ' '  ! Delete last dash
        EXIT
      END IF
    END DO

    CALL message('State', this%name)
    CALL message('  Description', this%description)
    CALL message('  Tree path', TRIM(c_path))
    IF (this%visited(no_omp_thread)) THEN
      CALL message('  Visited', '.TRUE.')
    ELSE
      CALL message('  Visited', '.FALSE.')
    END IF

  END SUBROUTINE Print_state

  ! Only works like this for two states in the same state machine!
  FUNCTION Is_equal_state(this, that) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = ALL(this%path == that%path)

  END FUNCTION Is_equal_state

  FUNCTION Is_not_equal_state(this, that) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = ANY(this%path /= that%path)

  END FUNCTION Is_not_equal_state

  FUNCTION Process_message_state(this, msg) RESULT(return_ptr)

    CLASS(t_State),   INTENT(inout) :: this
    CLASS(t_Message), INTENT(in),   POINTER :: msg
    CLASS(t_Message), POINTER       :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = modname//':Process_message_state'

    INTEGER :: no_omp_thread

    IF (.NOT. ASSOCIATED(msg)) RETURN

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

!!$    this%msg => msg
!!$    this%hsm => hsm
    !CALL this%Print_state()

    return_ptr => this%Handler(msg)

    SELECT CASE (TRIM(msg%action%name))
    CASE ('START')
      this%visited(no_omp_thread) = .TRUE.
    CASE ('ENTER')
      this%visited(no_omp_thread) = .TRUE.
    CASE ('EXIT')
      this%visited(no_omp_thread) = .FALSE.
    END SELECT

  END FUNCTION Process_message_state

  ! ===============================================================================
  ! Methods for hsm class
  ! ===============================================================================
  SUBROUTINE Init_hsm(self, name, debug)

    CLASS(t_Hsm),      INTENT(inout) :: self
    CHARACTER(len=*),  INTENT(in)    :: name
    LOGICAL, OPTIONAL, INTENT(in)    :: debug

    !TYPE(t_Action), POINTER :: action
    INTEGER :: i, nthreads
    LOGICAL :: l_debug

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    l_debug = .FALSE.
    IF (PRESENT(debug)) l_debug = debug

    CALL self%Set_mode('full')

    ALLOCATE(self%current(nthreads))
    ALLOCATE(self%next   (nthreads))
    ALLOCATE(self%source (nthreads))
    ALLOCATE(self%levels_to_lca(nthreads))
    DO i=1,nthreads
      self%current(i)%p => NULL()
      self%next   (i)%p => NULL()
      self%source (i)%p => NULL()
      self%levels_to_lca(i) = -1
    END DO

    ALLOCATE(self%actions(MAX_ACTIONS))
    !NULLIFY(action)
    DO i=1,MAX_ACTIONS
      self%actions(i) = t_Action("")
      !action => t_Action("")
      !self%actions(i) = t_Action("")
      !DEALLOCATE(action)
    END DO
    CALL self%Register_action(t_Action("ENTER"), l_debug)
    CALL self%Register_action(t_Action("EXIT"),  l_debug)
    CALL self%Register_action(t_Action("START"), l_debug)

    self%name = TRIM(name)

  END SUBROUTINE Init_hsm

  SUBROUTINE Set_mode_hsm(self, mode)

    CLASS(t_Hsm),     INTENT(inout) :: self
    CHARACTER(len=*), INTENT(in), OPTIONAL :: mode

    SELECT CASE (TRIM(mode))
    CASE ('full', 'FULL')
      self%mode = HSM_full
    CASE ('simple', 'SIMPLE')
      self%mode = HSM_simple
   END SELECT

  END SUBROUTINE Set_mode_hsm

  FUNCTION Get_current_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm),   INTENT(in) :: this
    CLASS(t_State), POINTER    :: return_ptr

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    return_ptr => this%current(no_omp_thread)%p

  END FUNCTION Get_current_hsm

  SUBROUTINE Reset_state_hsm(this)

    CLASS(t_hsm), INTENT(inout) :: this

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    this%current(no_omp_thread)%p => this%top
    this%top%visited(no_omp_thread) = .FALSE.

  END SUBROUTINE Reset_state_hsm

  ! Resets %visited(:) to .FALSE. and sets hsm%current to top state for all tiles and all OMP threasds
  SUBROUTINE Reset_all_hsm(this)

    CLASS(t_hsm), INTENT(inout) :: this

    INTEGER :: nthreads, i

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    CALL this%top%Reset_all()

    DO i=1,nthreads
      this%current(i)%p => this%top
    END DO

  END SUBROUTINE Reset_all_hsm

  ! Get next state in pre-order traversal
  FUNCTION Get_next_state_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm), INTENT(inout) :: this
    CLASS(t_State), POINTER     :: return_ptr

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    IF (this%current(no_omp_thread)%p%visited(no_omp_thread)) THEN   ! State has been visited before on traversal
      this%current(no_omp_thread)%p%visited(no_omp_thread) = .FALSE.
      IF (this%current(no_omp_thread)%p == this%top) THEN
        return_ptr => NULL()
        RETURN
      ELSE
        IF (ASSOCIATED(this%current(no_omp_thread)%p%next_sibling)) THEN
          return_ptr => this%current(no_omp_thread)%p%next_sibling
        ELSE
          return_ptr => this%current(no_omp_thread)%p%parent
        END IF
      END IF
    ELSE                             ! First visit of this state during traversal
      IF (this%current(no_omp_thread)%p%Has_children()) THEN
        this%current(no_omp_thread)%p%visited(no_omp_thread) = .TRUE.
        return_ptr => this%current(no_omp_thread)%p%first_child
      ELSE IF (ASSOCIATED(this%current(no_omp_thread)%p%next_sibling)) THEN
        return_ptr => this%current(no_omp_thread)%p%next_sibling
      ELSE
        return_ptr => this%current(no_omp_thread)%p%parent
      END IF
    END IF

  END FUNCTION Get_next_state_hsm

  ! Set Hsm to next state in pre-order traversal
  FUNCTION Goto_next_state_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm), INTENT(inout) :: this
    CLASS(t_State), POINTER     :: return_ptr

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    !this%current%visited = .NOT. this%current%visited
    this%current(no_omp_thread)%p => this%Get_next()
    return_ptr => this%current(no_omp_thread)%p

  END FUNCTION Goto_next_state_hsm

  FUNCTION Get_state_hsm(this, path) RESULT(return_ptr)

    CLASS(t_hsm), INTENT(in) :: this
    INTEGER,      INTENT(in) :: path(:)
    CLASS(t_State), POINTER  :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_state_hsm'

    INTEGER :: depth, i, j

    IF (SIZE(path) > MAX_STATE_NESTING) &
      & CALL finish(TRIM(routine), 'Path too long')

    return_ptr => this%top

    IF (ALL(path > 0)) THEN
      depth = SIZE(path)
    ELSE
      depth = MINLOC(path, DIM=1) - 1
    END IF
    DO i=2,depth
      return_ptr => return_ptr%first_child
      DO j=2,path(i)
        return_ptr => return_ptr%next_sibling
      END DO
    END DO

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Path to invalid state')

  END FUNCTION Get_state_hsm

  ! State transition
  SUBROUTINE Goto_state_hsm(this, that, debug)

    CLASS(t_hsm), INTENT(inout) :: this
    CLASS(t_State), INTENT(in), TARGET :: that
    LOGICAL, OPTIONAL, INTENT(in) :: debug

    LOGICAL :: l_debug

    CHARACTER(len=*), PARAMETER :: routine = modname//':Goto_state_hsm'

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    l_debug = .FALSE.
    IF (PRESENT(debug)) l_debug = debug

    IF (l_debug) CALL message(TRIM(routine), 'Taking transition from '//TRIM(this%current(no_omp_thread)%p%name)//&
      &                                                         ' to '//TRIM(that%name))

    IF (this%mode == HSM_full) THEN

      IF (ASSOCIATED(this%next(no_omp_thread)%p)) &
        & CALL finish(TRIM(routine), 'Transition already in progress')

      this%levels_to_lca(no_omp_thread) = this%Get_levels_to_LCA(that)

    END IF

    this%next(no_omp_thread)%p => that

  END SUBROUTINE Goto_state_hsm
  !>
  !! Get number of levels from state machines source state to the least common
  !! ancestor with some other target state.
  !!
  FUNCTION Get_levels_to_LCA_hsm(this, that) RESULT(return_value)

    CLASS(t_hsm),   INTENT(in)         :: this
    CLASS(t_State), INTENT(in), TARGET :: that
    INTEGER                            :: return_value

    CLASS(t_State), POINTER :: sstate, tstate

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_levels_to_LCA_hsm'

    INTEGER :: no_omp_thread

    return_value = 0

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    sstate => this%source(no_omp_thread)%p
    IF (sstate == that) THEN
      return_value = 1
      RETURN
    END IF

    DO WHILE (ASSOCIATED(sstate))
      tstate => that
      DO WHILE (ASSOCIATED(tstate))
        IF (sstate == tstate) THEN
!!$          CALL message(TRIM(routine), 'Levels to LCA from '//TRIM(this%source%name)//' to '//&
!!$            &                                                TRIM(that%name)//': '//int2string(return_value))
          RETURN
        END IF
        tstate => tstate%parent
      END DO
      return_value = return_value + 1
      sstate => sstate%parent
    END DO

  END FUNCTION Get_levels_to_LCA_hsm

  ! Take start transition (no states need to be exited)
  SUBROUTINE Start_state_hsm(this, that)

    CLASS(t_Hsm),   INTENT(inout)      :: this
    CLASS(t_State), INTENT(in), TARGET :: that

    CHARACTER(len=*), PARAMETER :: routine = modname//':Start_state_hsm'

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

!!$    CALL message(TRIM(routine), 'Taking start transition from '//TRIM(this%current%name)//&
!!$      &                                                  ' to '//TRIM(that%name))

    IF (ASSOCIATED(this%next(no_omp_thread)%p)) &
      & CALL warning(TRIM(routine), 'Transition to '//TRIM(this%next(no_omp_thread)%p%name)//' already in progress')

    this%next(no_omp_thread)%p => that

  END SUBROUTINE Start_state_hsm

  SUBROUTINE Start_hsm(this)

    CLASS(t_Hsm), INTENT(inout) :: this

    CLASS(t_Message), POINTER :: start_msg, entry_msg
    CLASS(t_Message), POINTER :: tmp_msg

    CHARACTER(len=*), PARAMETER :: routine = modname//':Start_hsm'

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    IF (this%mode == HSM_simple) THEN
      this%current(no_omp_thread)%p => this%top
      RETURN
    END IF

!!$    CALL message(TRIM(routine), 'Starting hsm ...')

    !ALLOCATE(t_Message::start_msg)
    !ALLOCATE(t_Message::entry_msg)
    !start_msg%action = this%Get_action('START')
    !entry_msg%action = this%Get_action('ENTER')
    !entry_msg%action = t_action('ENTER')

    start_msg => t_Message('', this%Get_action('START'))
    entry_msg => t_Message('', this%Get_action('ENTER'))

    this%current(no_omp_thread)%p => this%top
    !this%next    => NULL()
    tmp_msg => this%current(no_omp_thread)%p%Process_message(entry_msg)
    IF(ASSOCIATED(tmp_msg)) DEALLOCATE(tmp_msg)

    tmp_msg => this%current(no_omp_thread)%p%Process_message(start_msg)
    IF (ASSOCIATED(tmp_msg)) DEALLOCATE(tmp_msg)
    DO WHILE (ASSOCIATED(this%next(no_omp_thread)%p))
      CALL this%Process_entry_to_next()
      this%current(no_omp_thread)%p => this%next(no_omp_thread)%p
      NULLIFY(this%next(no_omp_thread)%p)
      tmp_msg => this%current(no_omp_thread)%p%Process_message(start_msg)
      IF (ASSOCIATED(tmp_msg)) DEALLOCATE(tmp_msg)
    END DO

    !DEALLOCATE(start_msg%action, entry_msg%action)
    DEALLOCATE(start_msg, entry_msg)

!!$    CALL message(TRIM(routine), 'Hsm started.')

  END SUBROUTINE Start_hsm
  !
  ! State machine engine
  !
  SUBROUTINE Process_message_hsm(this, msg, debug)

    CLASS(t_Hsm),     TARGET, INTENT(inout) :: this
!!$    CLASS(t_Message), POINTER       :: msg
    CLASS(t_Message), POINTER, INTENT(in) :: msg
    LOGICAL, OPTIONAL, INTENT(in) :: debug

    CLASS(t_State),   POINTER :: state
    CLASS(t_Message), POINTER :: new_msg, tmp_msg
    CLASS(t_Message),  POINTER :: start_msg
    !CLASS(t_Hsm),     POINTER :: hsm
    LOGICAL :: l_debug

    CHARACTER(len=*), PARAMETER :: routine = modname//':Process_message_hsm'

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    l_debug = .FALSE.
    IF (PRESENT(debug)) l_debug = debug

!!$    CALL message(TRIM(routine), 'Processing message '//TRIM(msg%name)//' on hsm ...')

    IF (.NOT. ASSOCIATED(msg)) RETURN

    SELECT CASE (this%mode)
    CASE (HSM_simple)

      ! Transitions are only taken to parent, next sibling or first child
      ALLOCATE(new_msg, source=msg)
      state => this%current(no_omp_thread)%p
      IF (l_debug) CALL message(TRIM(routine), 'Starting loop, message action is '//TRIM(msg%action%name))
      DO WHILE (ASSOCIATED(state))
        ! Process message on current state
        IF (l_debug) CALL message(TRIM(routine), 'State is '//TRIM(state%name))
        tmp_msg => state%Process_message(new_msg)
        IF (ASSOCIATED(new_msg)) THEN
          DEALLOCATE(new_msg)
        END IF
        new_msg => tmp_msg    ! new_msg not null: message not processed, or message processed and new message_text
                              ! new_msg is null : message processed and no new message
        !
        ! Message not processed or new message and transition: goto next state with message
        ! Message not processed or new message and no transition: goto parent with message
        ! Message processed and transition: goto next state and exit loop
        ! Message processed and no transition: exit loop
        IF (ASSOCIATED(this%next(no_omp_thread)%p)) THEN  ! Take transition to next
          ! CALL this%Process_entry_to_next()
          IF (this%next(no_omp_thread)%p%level == this%current(no_omp_thread)%p%level) THEN
            ! Goto to a sibling: exit current state
            this%current(no_omp_thread)%p%visited = .FALSE.   ! Exit current on way up
          END IF
          this%current(no_omp_thread)%p => this%next(no_omp_thread)%p
          NULLIFY(this%next(no_omp_thread)%p)
          this%current(no_omp_thread)%p%visited = .TRUE.  ! Enter new state
          IF (.NOT. ASSOCIATED(new_msg)) THEN ! Message processed
            EXIT
          ELSE
            state => this%current(no_omp_thread)%p
          END IF
        ELSE                                ! No state transition taken
          IF (ASSOCIATED(new_msg)) THEN
            state%visited = .FALSE.         ! Exit current state
            state => state%parent           ! Goto parent
          ELSE
            EXIT
          END IF
        END IF
      END DO


    CASE (HSM_full)

      !hsm => this
      ALLOCATE(t_Message::start_msg)
      start_msg%action = this%Get_action('START')

  !!$    CALL this%Process_exit_to_LCA()

      ALLOCATE(new_msg, source=msg)
      state => this%current(no_omp_thread)%p
      IF (l_debug) CALL message(TRIM(routine), 'Starting loop, message action is '//TRIM(msg%action%name))
      DO WHILE (ASSOCIATED(state))
        IF (l_debug) CALL message(TRIM(routine), 'State is '//TRIM(state%name))
        this%source(no_omp_thread)%p => state  ! Set source always to the previous state. This is to find back again later.
        tmp_msg => state%Process_message(new_msg)
        IF (ASSOCIATED(new_msg)) THEN
        !  IF (ASSOCIATED(new_msg%action)) DEALLOCATE(new_msg%action)
          DEALLOCATE(new_msg)
        END IF
        new_msg => tmp_msg ! Set new message to the "result" of last message
        CALL this%Process_exit_to_LCA()
        IF (ASSOCIATED(this%next(no_omp_thread)%p)) THEN      ! State transition taken (i.e. already exited to LCA)
          ! Execute entry and start messages to target
          CALL this%Process_entry_to_next()
          this%current(no_omp_thread)%p => this%next(no_omp_thread)%p
          NULLIFY(this%next(no_omp_thread)%p)
          tmp_msg => this%current(no_omp_thread)%p%Process_message(start_msg)
          IF (ASSOCIATED(tmp_msg)) THEN
            !IF (ASSOCIATED(tmp_msg%action)) DEALLOCATE(tmp_msg%action)
            DEALLOCATE(tmp_msg)
          END IF
          DO WHILE (ASSOCIATED(this%next(no_omp_thread)%p))
            CALL this%Process_entry_to_next()
            this%current(no_omp_thread)%p => this%next(no_omp_thread)%p
            NULLIFY(this%next(no_omp_thread)%p)
            tmp_msg => this%current(no_omp_thread)%p%Process_message(start_msg)
            IF (ASSOCIATED(tmp_msg)) THEN
              !IF (ASSOCIATED(tmp_msg%action)) DEALLOCATE(tmp_msg%action)
              DEALLOCATE(tmp_msg)
            END IF
          END DO
          IF (.NOT. ASSOCIATED(new_msg)) THEN ! Message processed
            EXIT
          ELSE
            state => this%current(no_omp_thread)%p
          END IF
        ELSE                                ! No state transition taken
          IF (ASSOCIATED(new_msg)) THEN
            state => state%parent
          ELSE
            EXIT
          END IF
        END IF
      END DO

      !DEALLOCATE(start_msg%action)
      DEALLOCATE(start_msg)

    END SELECT

    IF (ASSOCIATED(new_msg)) THEN
      !IF (ASSOCIATED(new_msg%action)) DEALLOCATE(new_msg%action)
      DEALLOCATE(new_msg)
    END IF

    IF (l_debug) THEN
      IF (ASSOCIATED(new_msg)) THEN
        CALL message(TRIM(routine), 'Exited loop, new message action is '//TRIM(new_msg%action%name))
      ELSE
        CALL message(TRIM(routine), 'Exited loop, new message action is NULL')
      END IF
    END IF

!!$    CALL message(TRIM(routine), 'Processing finished.')

  END SUBROUTINE Process_message_hsm

  ! Process entry messages from LCA to target
  SUBROUTINE Process_entry_to_next(this)

    CLASS(t_Hsm), INTENT(inout) :: this

    CLASS(t_Message), POINTER :: entry_msg
    CLASS(t_Message), POINTER :: tmp_msg
    CLASS(t_State), POINTER :: state
    TYPE(t_State_p) :: entrypath(MAX_STATE_NESTING)
    INTEGER         :: i, j

    CHARACTER(len=*), PARAMETER :: routine = modname//':Process_entry_to_next'

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

!!$    CALL message(TRIM(routine), 'Processing entry messages to next ...')

    IF (.NOT. ASSOCIATED(this%current(no_omp_thread)%p)) THEN
      RETURN
    END IF

    SELECT CASE (this%mode)
    CASE (HSM_simple)
      this%next(no_omp_thread)%p%visited = .TRUE.
    CASE (HSM_full)
      ALLOCATE(t_Message::entry_msg)
      entry_msg%action = this%Get_action('ENTER')

      i = 0
      state => this%next(no_omp_thread)%p
      DO WHILE (state /= this%current(no_omp_thread)%p)
        i = i + 1
        entrypath(i)%p => state
        state => state%parent
        IF (.NOT. ASSOCIATED(state)) EXIT
      END DO

      DO j=i,1,-1
        state => entrypath(j)%p
        tmp_msg => state%Process_message(entry_msg)
        IF (ASSOCIATED(tmp_msg)) DEALLOCATE(tmp_msg)
      END DO

      !DEALLOCATE(entry_msg%action)
      DEALLOCATE(entry_msg)
    END SELECT

!!$    CALL message(TRIM(routine), 'Processing of entry messages finished.')

  END SUBROUTINE Process_entry_to_next

  ! Process exit messages from current state up to LCA with target state
  ! (set during a previous state transition)
!  SUBROUTINE Process_exit_to_lca(this, that)
  SUBROUTINE Process_exit_to_lca(this)

    CLASS(t_Hsm),   INTENT(inout) :: this
!    CLASS(t_State), INTENT(in)    :: that

    CLASS(t_Message), POINTER :: exit_msg
    CLASS(t_Message), POINTER :: tmp_msg
    CLASS(t_State),   POINTER :: state
    INTEGER                   :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Process_exit_to_lca'

    INTEGER :: no_omp_thread

#ifdef _OPENMP
    no_omp_thread = omp_get_thread_num() + 1
#else
    no_omp_thread = 1
#endif

    IF (this%levels_to_lca(no_omp_thread) == -1) RETURN

!!$    CALL message(TRIM(routine), 'Processing exit messages to lca ...'//int2string(this%levels_to_lca))

    ALLOCATE(t_Message::exit_msg)
    exit_msg%action = this%Get_action('EXIT')

    state => this%current(no_omp_thread)%p

    IF (ASSOCIATED(this%source(no_omp_thread)%p)) THEN
      DO WHILE (state /= this%source(no_omp_thread)%p)
        tmp_msg => state%Process_message(exit_msg)
        IF (ASSOCIATED(tmp_msg)) DEALLOCATE(tmp_msg)
        state => state%parent
        IF (.NOT. ASSOCIATED(state)) EXIT
      END DO
    END IF

    state => this%source(no_omp_thread)%p
!    DO i=this%Get_levels_to_LCA(that),1,-1
    DO i=this%levels_to_lca(no_omp_thread),1,-1
      tmp_msg => state%Process_message(exit_msg)
      IF (ASSOCIATED(tmp_msg)) DEALLOCATE(tmp_msg)
      state => state%parent
      IF (.NOT. ASSOCIATED(state)) EXIT
    END DO

    this%current(no_omp_thread)%p => state

    this%levels_to_lca(no_omp_thread) = -1

    !DEALLOCATE(exit_msg%action)
    DEALLOCATE(exit_msg)

!!$    CALL message(TRIM(routine), 'Processing of exit messages finished.')

  END SUBROUTINE Process_exit_to_lca

  SUBROUTINE Register_action_hsm(this, action, debug)

    CLASS(t_Hsm),   INTENT(inout) :: this
    TYPE(t_Action), INTENT(in)    :: action
    LOGICAL, OPTIONAL, INTENT(in) :: debug

    INTEGER :: i
    LOGICAL :: l_debug

    CHARACTER(len=*), PARAMETER :: routine = modname//':Register_action_hsm'

    l_debug = .FALSE.
    IF (PRESENT(debug)) l_debug = debug

    i = 1
    DO WHILE (TRIM(this%actions(i)%name) /= "")
      i = i + 1
    END DO
    IF (i > MAX_ACTIONS) THEN
      CALL finish(TRIM(routine), 'Maximum number of actions exceeded!')
    ELSE
      this%actions(i) = action
      IF (l_debug) CALL message(TRIM(routine), 'Registering state machine action '//TRIM(this%actions(i)%name))
    END IF

  END SUBROUTINE Register_action_hsm

  FUNCTION Get_action_hsm(this, name) RESULT(return_value)

    CLASS(t_Hsm),     TARGET, INTENT(in) :: this
    CHARACTER(len=*),         INTENT(in) :: name
    TYPE(t_Action)                       :: return_value

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_action_hsm'

    !return_ptr => NULL()

    DO i=1,MAX_ACTIONS
      IF (TRIM(this%actions(i)%name) == TRIM(name)) THEN
        !ALLOCATE(return_ptr, source=this%actions(i))
        return_value = this%actions(i)
        EXIT
      END IF
    END DO

    !IF (.NOT. ASSOCIATED(return_ptr)) &
    IF (return_value%name == '') &
      & CALL finish(TRIM(routine), 'Action '//TRIM(name)//' not found')

  END FUNCTION Get_action_hsm

  SUBROUTINE Print_actions_hsm(this)

    CLASS(t_Hsm), INTENT(in) :: this

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Print_actions_hsm'

    message_text = "State machine actions:"
    DO i=1,MAX_ACTIONS
      IF (TRIM(this%actions(i)%name) /= "") &
        & message_text = TRIM(message_text)//" "//TRIM(this%actions(i)%name)
    END DO
    CALL message(TRIM(routine), TRIM(message_text)//'.')

  END SUBROUTINE Print_actions_hsm

  ! ===============================================================================
  ! Methods for t_Message class
  ! ===============================================================================
  FUNCTION Construct_message(name, action) RESULT(return_ptr)

    CHARACTER(len=*), INTENT(in) :: name
    TYPE(t_Action),   INTENT(in) :: action
    TYPE(t_Message),  POINTER    :: return_ptr

    ALLOCATE(return_ptr)
    !ALLOCATE(t_Action::return_ptr)
    return_ptr%name = TRIM(name)
    return_ptr%action%name = TRIM(action%name)

  END FUNCTION Construct_message

  ! ===============================================================================
  ! Methods for t_Action class
  ! ===============================================================================
  FUNCTION Construct_action(name) RESULT(return_ptr)

    CHARACTER(len=*), INTENT(in) :: name
    TYPE(t_Action), POINTER :: return_ptr

    ALLOCATE(return_ptr)
    !ALLOCATE(t_Action::return_ptr)
    return_ptr%name = TRIM(name)

  END FUNCTION Construct_action

  FUNCTION Is_equal_action(this, that) RESULT(return_value)

    CLASS(t_Action), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = TRIM(this%name) == TRIM(that%name)

  END FUNCTION Is_equal_action

  FUNCTION Is_not_equal_action(this, that) RESULT(return_value)

    CLASS(t_Action), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = TRIM(this%name) /= TRIM(that%name)

  END FUNCTION Is_not_equal_action

END MODULE mo_hsm_class
