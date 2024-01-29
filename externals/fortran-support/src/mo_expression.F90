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

! ---------------------------------------------------------------------
! Evaluation of basic arithmetic expressions specified by
! character-strings (infix expressions) based on a Finite State
! Machine (FSM) and Dijkstra's shunting yard algorithm.
! It is possible to include mathematical functions, operators, and
! constants, see the LaTeX documentation for this module in the
! appendix of the namelist documentaion. Besides, Fortran variables
! can be linked to the expression and used in the evaluation. The
! implementation supports scalar input variables as well as 2D and 3D
! fields, where it is implicitly assumed that 2D fields are embedded
! in 3D fields as "3D(:,level,:) = 2D(:,:)".
!
! Basic usage example:
!
!  formula = expression("sin(45*pi/180.) * 10 + 5")
!  CALL formula%evaluate(val)
!
!
! Note that the implementation of the expression parser is contained
! in the C module "util_arithmetic_expr".
! ---------------------------------------------------------------------

MODULE mo_expression
  USE, INTRINSIC :: ISO_C_BINDING
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: real64
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: real_kind = real64

  INTEGER, PARAMETER :: MAX_BUF_LEN = 1024 ! max length of expression operator list
  INTEGER, PARAMETER :: RESULT_STACKSIZE = 1024 ! max size of result stack
  INTEGER, PARAMETER :: MAX_NAME_LEN = 32 ! max length of variable name

  PUBLIC :: expression
  PUBLIC :: parse_expression_string

  ! ----------------------------------------------------------------------
  ! ISO-C bindings
  ! ----------------------------------------------------------------------

  ENUM, BIND(c)
    ENUMERATOR :: VARIABLE, VALUE, OPERATOR, FUNCTION
  END ENUM

  ENUM, BIND(c)
    ENUMERATOR :: fEXP = 0, fLOG = 1, fSIN = 2, fCOS = 3, &
      &           fMIN = 4, fMAX = 5, fIF = 6, fSQRT = 7, &
      &           fERF = 8
  END ENUM

  TYPE, BIND(c) :: t_item
    INTEGER(kind=c_int)               :: etype
    CHARACTER(kind=c_char)            :: op
    REAL(kind=c_double)               :: val
    INTEGER(kind=c_int)               :: fct
    CHARACTER(kind=c_char, len=1)     :: field(32)
  END TYPE t_item

  TYPE, BIND(c) :: t_list
    INTEGER(kind=c_int)               :: isize
    TYPE(t_item)                      :: list(MAX_BUF_LEN)
  END TYPE t_list

  INTERFACE
    FUNCTION my_do_parse_infix(in_parse_line, queue) RESULT(ierr) BIND(c, name='do_parse_infix')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPORT t_list
      CHARACTER(c_char), DIMENSION(*) :: in_parse_line
      TYPE(t_list)                    :: queue
      INTEGER(c_int)                  :: ierr
    END FUNCTION my_do_parse_infix
  END INTERFACE

  ! ----------------------------------------------------------------------
  ! Derived data types (Fortran)
  ! ----------------------------------------------------------------------

  ! single field on result stack: scalar, 2D or 3D results
  TYPE t_result_item
    REAL(real_kind), POINTER :: ptr_0D => NULL()
    REAL(real_kind), POINTER :: ptr_2D(:, :) => NULL()
    REAL(real_kind), POINTER :: ptr_3D(:, :, :) => NULL()
  END TYPE t_result_item

  ! stack for result computation, base class
  TYPE t_result_stack
    INTEGER                :: ridx = 0
    TYPE(t_result_item)    :: rstack(RESULT_STACKSIZE)
  END TYPE t_result_stack

  ! stack for result computation, scalar results
  TYPE, EXTENDS(t_result_stack) :: t_result_stack_0D
  CONTAINS
    PROCEDURE, PUBLIC :: pop => result_stack_pop_0D
    PROCEDURE, PUBLIC :: push => result_stack_push_0D
  END TYPE t_result_stack_0D

  ! stack for result computation, 2D results
  TYPE, EXTENDS(t_result_stack) ::  t_result_stack_2D
    INTEGER                :: dim1, dim2
  CONTAINS
    PROCEDURE, PUBLIC :: pop => result_stack_pop_2D
    PROCEDURE, PUBLIC :: push => result_stack_push_2D
  END TYPE t_result_stack_2D

  ! user-defined constructor
  INTERFACE t_result_stack_2D
    MODULE PROCEDURE t_result_stack_2D_init
  END INTERFACE t_result_stack_2D

  ! stack for result computation, 3D results
  TYPE, EXTENDS(t_result_stack) ::  t_result_stack_3D
    INTEGER                :: dim1, dim2, dim3
  CONTAINS
    PROCEDURE, PUBLIC :: pop => result_stack_pop_3D
    PROCEDURE, PUBLIC :: push => result_stack_push_3D
  END TYPE t_result_stack_3D

  ! user-defined constructor
  INTERFACE t_result_stack_3D
    MODULE PROCEDURE t_result_stack_3D_init
  END INTERFACE t_result_stack_3D

  ! table of error codes
  ENUM, BIND(C)
    ENUMERATOR :: ERR_NONE,               &
      &           ERR_VARIABLE_NOT_FOUND, &
      &           ERR_EXPRESSION_EMPTY,   &
      &           ERR_PARSE_ERROR,        &
      &           ERR_DIMENSION_MISMATCH, &
      &           ERR_INTERNAL
  END ENUM

  ! base class for single expression operating on result stack
  TYPE, ABSTRACT :: t_stack_op
    INTEGER :: err_no = ERR_NONE
  CONTAINS
    PROCEDURE(stack_op_eval_0D), PUBLIC, DEFERRED :: eval_0D
    PROCEDURE(stack_op_eval_2D), PUBLIC, DEFERRED :: eval_2D
    PROCEDURE(stack_op_eval_3D), PUBLIC, DEFERRED :: eval_3D
    GENERIC, PUBLIC :: eval => eval_0D, eval_2D, eval_3D
  END TYPE t_stack_op

  TYPE t_expr_ptr
    CLASS(t_stack_op), POINTER :: p => NULL()
  END TYPE t_expr_ptr

  ! pointer (link) to an associated variable
  TYPE t_var_ptr
    CHARACTER(len=MAX_NAME_LEN) :: name
    TYPE(t_result_item)         :: p
  END TYPE t_var_ptr

  ! list of expressions executed first-to-last
  TYPE expression
    INTEGER                                  :: isize, nvars
    TYPE(t_expr_ptr)                         :: list(MAX_BUF_LEN)
    TYPE(t_var_ptr), ALLOCATABLE :: var(:)

    INTEGER                                  :: err_no = ERR_NONE
    CHARACTER(len=MAX_NAME_LEN), ALLOCATABLE :: used_vars(:)
  CONTAINS
    PROCEDURE, PUBLIC :: finalize => expr_list_finalize
    PROCEDURE, PUBLIC :: evaluate_0D => expr_list_evaluate_0D
    PROCEDURE, PUBLIC :: evaluate_2D => expr_list_evaluate_2D
    PROCEDURE, PUBLIC :: evaluate_3D => expr_list_evaluate_3D
    PROCEDURE, PUBLIC :: link_0D => expr_list_link_0D
    PROCEDURE, PUBLIC :: link_2D => expr_list_link_2D
    PROCEDURE, PUBLIC :: link_3D => expr_list_link_3D
    PROCEDURE, PUBLIC :: resize_varlist => expr_list_resize_varlist
    PROCEDURE, PUBLIC :: get_result_shape => expr_list_get_result_shape
    GENERIC, PUBLIC   :: link => link_0D, link_2D, link_3D
    GENERIC, PUBLIC   :: evaluate => evaluate_0D, evaluate_2D, evaluate_3D
    ! FINAL :: expr_list_finalize  ! gfortran < 4.9 does not support finalizers :(
  END TYPE expression

  ! user-defined constructor
  INTERFACE expression
    MODULE PROCEDURE expression_init
  END INTERFACE expression

  ! interface definition for expressions operating on result stack
  ABSTRACT INTERFACE
    SUBROUTINE stack_op_eval_0D(this, result_stack)
      IMPORT t_stack_op, t_result_stack_0D
      CLASS(t_stack_op), INTENT(INOUT)        :: this
      CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    END SUBROUTINE stack_op_eval_0D
  END INTERFACE

  ! interface definition for expressions operating on result stack
  ABSTRACT INTERFACE
    SUBROUTINE stack_op_eval_2D(this, result_stack)
      IMPORT t_stack_op, t_result_stack_2D
      CLASS(t_stack_op), INTENT(INOUT)        :: this
      CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    END SUBROUTINE stack_op_eval_2D
  END INTERFACE

  ! interface definition for expressions operating on result stack
  ABSTRACT INTERFACE
    SUBROUTINE stack_op_eval_3D(this, result_stack)
      IMPORT t_stack_op, t_result_stack_3D
      CLASS(t_stack_op), INTENT(INOUT)        :: this
      CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    END SUBROUTINE stack_op_eval_3D
  END INTERFACE

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_value ! add value to result stack
    REAL(real_kind) :: val
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_value_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_value_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_value_eval_3D
  END TYPE t_op_value

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_variable ! add variable to result stack
    TYPE(t_var_ptr) :: var
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_variable_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_variable_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_variable_eval_3D
  END TYPE t_op_variable

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_plus ! addition "+"
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_plus_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_plus_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_plus_eval_3D
  END TYPE t_op_plus

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_minus ! subtraction "-"
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_minus_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_minus_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_minus_eval_3D
  END TYPE t_op_minus

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_mul ! multiplication "*"
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_mul_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_mul_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_mul_eval_3D
  END TYPE t_op_mul

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_div ! division "/"
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_div_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_div_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_div_eval_3D
  END TYPE t_op_div

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_pow ! exponentiation "^"
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_pow_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_pow_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_pow_eval_3D
  END TYPE t_op_pow

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_gt ! greater-than ">"
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_gt_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_gt_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_gt_eval_3D
  END TYPE t_op_gt

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_lt ! less-than "<"
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_lt_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_lt_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_lt_eval_3D
  END TYPE t_op_lt

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_exp ! EXP()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_exp_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_exp_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_exp_eval_3D
  END TYPE t_op_exp

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_log ! LOG()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_log_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_log_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_log_eval_3D
  END TYPE t_op_log

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_sin ! SIN()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_sin_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_sin_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_sin_eval_3D
  END TYPE t_op_sin

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_cos ! COS()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_cos_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_cos_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_cos_eval_3D
  END TYPE t_op_cos

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_min ! MIN()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_min_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_min_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_min_eval_3D
  END TYPE t_op_min

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_max ! MAX()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_max_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_max_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_max_eval_3D
  END TYPE t_op_max

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_if ! IF()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_if_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_if_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_if_eval_3D
  END TYPE t_op_if

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_sqrt ! SQRT()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_sqrt_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_sqrt_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_sqrt_eval_3D
  END TYPE t_op_sqrt

  TYPE, PUBLIC, EXTENDS(t_stack_op) :: t_op_erf ! ERF()
  CONTAINS
    PROCEDURE, PUBLIC :: eval_0D => stack_op_erf_eval_0D
    PROCEDURE, PUBLIC :: eval_2D => stack_op_erf_eval_2D
    PROCEDURE, PUBLIC :: eval_3D => stack_op_erf_eval_3D
  END TYPE t_op_erf

CONTAINS

  FUNCTION t_result_stack_2D_init(dim1, dim2) RESULT(rstack)
    INTEGER, INTENT(IN) :: dim1, dim2
    TYPE(t_result_stack_2D) :: rstack
    rstack%dim1 = dim1
    rstack%dim2 = dim2
  END FUNCTION t_result_stack_2D_init

  FUNCTION t_result_stack_3D_init(dim1, dim2, dim3) RESULT(rstack)
    INTEGER, INTENT(IN) :: dim1, dim2, dim3
    TYPE(t_result_stack_3D) :: rstack
    rstack%dim1 = dim1
    rstack%dim2 = dim2
    rstack%dim3 = dim3
  END FUNCTION t_result_stack_3D_init

  SUBROUTINE result_stack_pop_0D(this, res)
    CLASS(t_result_stack_0D), INTENT(INOUT) :: this
    REAL(real_kind), POINTER :: res
    res => this%rstack(this%ridx)%ptr_0D
    NULLIFY (this%rstack(this%ridx)%ptr_0D)
    this%ridx = this%ridx - 1
  END SUBROUTINE result_stack_pop_0D

  SUBROUTINE result_stack_pop_2D(this, res)
    CLASS(t_result_stack_2D), INTENT(INOUT) :: this
    REAL(real_kind), POINTER :: res(:, :)
    res => this%rstack(this%ridx)%ptr_2D
    NULLIFY (this%rstack(this%ridx)%ptr_2D)
    this%ridx = this%ridx - 1
  END SUBROUTINE result_stack_pop_2D

  SUBROUTINE result_stack_pop_3D(this, res)
    CLASS(t_result_stack_3D), INTENT(INOUT) :: this
    REAL(real_kind), POINTER :: res(:, :, :)
    res => this%rstack(this%ridx)%ptr_3D
    NULLIFY (this%rstack(this%ridx)%ptr_3D)
    this%ridx = this%ridx - 1
  END SUBROUTINE result_stack_pop_3D

  SUBROUTINE result_stack_push_0D(this, val)
    CLASS(t_result_stack_0D), INTENT(INOUT) :: this
    REAL(real_kind), POINTER :: val
    this%ridx = this%ridx + 1
    this%rstack(this%ridx)%ptr_0D => val
  END SUBROUTINE result_stack_push_0D

  SUBROUTINE result_stack_push_2D(this, val)
    CLASS(t_result_stack_2D), INTENT(INOUT) :: this
    REAL(real_kind), POINTER                :: val(:, :)
    this%ridx = this%ridx + 1
    this%rstack(this%ridx)%ptr_2D => val
  END SUBROUTINE result_stack_push_2D

  SUBROUTINE result_stack_push_3D(this, val)
    CLASS(t_result_stack_3D), INTENT(INOUT) :: this
    REAL(real_kind), POINTER                :: val(:, :, :)
    this%ridx = this%ridx + 1
    this%rstack(this%ridx)%ptr_3D => val
  END SUBROUTINE result_stack_push_3D

  ! parse arithmetic expression; build a list of expression operators
  FUNCTION expression_init(string) RESULT(expr_list)
    CHARACTER(LEN=*), INTENT(IN)    :: string
    TYPE(expression) :: expr_list
    CALL parse_expression_string(expr_list, string)
  END FUNCTION expression_init

  SUBROUTINE parse_expression_string(expr_list, string)
    CHARACTER(LEN=*), INTENT(IN)    :: string
    TYPE(expression), INTENT(out)   :: expr_list

    ! local variables
    TYPE(t_list)                      :: cqueue
    INTEGER                           :: i, j, ierr, nvars_used
    CLASS(t_op_value), POINTER        :: op_value
    CLASS(t_op_variable), POINTER     :: op_variable
    CHARACTER(len=MAX_NAME_LEN)       :: used_vars(MAX_BUF_LEN)

    expr_list%err_no = ERR_NONE
    ! call C expression parser
    ierr = my_do_parse_infix(TRIM(ADJUSTL(string))//c_null_char, cqueue)
    ! translate C postfix queue into sequence of procedure pointers
    nvars_used = 0
    expr_list%nvars = 0
    expr_list%isize = cqueue%isize
    IF (ierr /= 0) THEN
      CALL expr_list_finalize(expr_list)
      expr_list%err_no = ERR_PARSE_ERROR
      RETURN
    END IF
    DO i = 1, cqueue%isize
      SELECT CASE (cqueue%list(i)%etype)
      CASE (VALUE)
        ALLOCATE (op_value)
        op_value%val = REAL(cqueue%list(i)%val)
        expr_list%list(i)%p => op_value
      CASE (VARIABLE)
        ALLOCATE (op_variable)
        op_variable%var%name = TRIM(convert_c_string(MAX_NAME_LEN, cqueue%list(i)%field))
        expr_list%list(i)%p => op_variable
        ! insert variable into "used_variables" list
        DO j = 1, nvars_used
          IF (TRIM(used_vars(j)) == TRIM(op_variable%var%name)) EXIT
        END DO
        IF (j >= nvars_used) THEN
          nvars_used = nvars_used + 1
          used_vars(nvars_used) = TRIM(op_variable%var%name)
        END IF
      CASE (OPERATOR)
        SELECT CASE (cqueue%list(i)%op)
        CASE ('+')
          ALLOCATE (t_op_plus::expr_list%list(i)%p)
        CASE ('-')
          ALLOCATE (t_op_minus::expr_list%list(i)%p)
        CASE ('*')
          ALLOCATE (t_op_mul::expr_list%list(i)%p)
        CASE ('/')
          ALLOCATE (t_op_div::expr_list%list(i)%p)
        CASE ('^')
          ALLOCATE (t_op_pow::expr_list%list(i)%p)
        CASE ('>')
          ALLOCATE (t_op_gt::expr_list%list(i)%p)
        CASE ('<')
          ALLOCATE (t_op_lt::expr_list%list(i)%p)
        END SELECT
      CASE (FUNCTION)
        SELECT CASE (cqueue%list(i)%fct)
        CASE (fEXP)
          ALLOCATE (t_op_exp::expr_list%list(i)%p)
        CASE (fLOG)
          ALLOCATE (t_op_log::expr_list%list(i)%p)
        CASE (fSIN)
          ALLOCATE (t_op_sin::expr_list%list(i)%p)
        CASE (fCOS)
          ALLOCATE (t_op_cos::expr_list%list(i)%p)
        CASE (fMIN)
          ALLOCATE (t_op_min::expr_list%list(i)%p)
        CASE (fMAX)
          ALLOCATE (t_op_max::expr_list%list(i)%p)
        CASE (fIF)
          ALLOCATE (t_op_if::expr_list%list(i)%p)
        CASE (fSQRT)
          ALLOCATE (t_op_sqrt::expr_list%list(i)%p)
        CASE (fERF)
          ALLOCATE (t_op_erf::expr_list%list(i)%p)
        END SELECT
      END SELECT
    END DO
    ALLOCATE (expr_list%used_vars(nvars_used))
    expr_list%used_vars(:) = used_vars(1:nvars_used)
  END SUBROUTINE parse_expression_string

  ! deallocate list of expression operators
  SUBROUTINE expr_list_finalize(this)
    CLASS(expression) :: this
    ! local variables
    INTEGER :: i
    DO i = 1, this%isize
      IF (ASSOCIATED(this%list(i)%p)) DEALLOCATE (this%list(i)%p)
    END DO
    this%isize = 0
    IF (ALLOCATED(this%used_vars)) DEALLOCATE (this%used_vars)
  END SUBROUTINE expr_list_finalize

  ! Find the dimensions of expression result. Besides, test if the
  ! variable has not yet been linked to this expression list.
  SUBROUTINE expr_list_get_result_shape(this, var_shape, dim1, dim2, dim3)
    CLASS(expression), INTENT(INOUT) :: this
    TYPE(t_var_ptr), INTENT(OUT)   :: var_shape ! input variable that has the result shape
    INTEGER, INTENT(OUT)   :: dim1, dim2, dim3 ! variable dimensions
    ! local variables
    INTEGER                    :: i, j
    LOGICAL                    :: found
    CLASS(t_stack_op), POINTER :: stack_op

    ! first step: determine the dimensions involved
    dim1 = 1
    dim2 = 1
    dim3 = 1
    DO i = 1, this%isize
      stack_op => this%list(i)%p
      ! variable expressions: test if the variable has not yet been
      ! linked to this expression list
      SELECT TYPE (stack_op)
      TYPE IS (t_op_variable)
        found = .FALSE.
        VARLOOP: DO j = 1, this%nvars
          IF (TRIM(this%var(j)%name) == TRIM(stack_op%var%name)) THEN
            stack_op%var = this%var(j)
            IF (ASSOCIATED(this%var(j)%p%ptr_0D)) THEN
              IF ((dim1 == 1) .AND. (dim2 == 1) .AND. (dim3 == 1)) THEN
                var_shape = this%var(j)
              END IF
            END IF
            IF (ASSOCIATED(this%var(j)%p%ptr_2D)) THEN
              dim1 = MAX(dim1, SIZE(this%var(j)%p%ptr_2d, 1))
              dim2 = MAX(dim2, SIZE(this%var(j)%p%ptr_2d, 2))
              IF ((dim1 == SIZE(this%var(j)%p%ptr_2d, 1)) .AND.  &
                & (dim2 == SIZE(this%var(j)%p%ptr_2d, 2)) .AND.  &
                & (dim3 == 1)) THEN
                var_shape = this%var(j)
              END IF
            ELSE IF (ASSOCIATED(this%var(j)%p%ptr_3D)) THEN
              dim1 = MAX(dim1, SIZE(this%var(j)%p%ptr_3d, 1))
              dim2 = MAX(dim2, SIZE(this%var(j)%p%ptr_3d, 2))
              dim3 = MAX(dim3, SIZE(this%var(j)%p%ptr_3d, 3))
              IF ((dim1 == SIZE(this%var(j)%p%ptr_3d, 1)) .AND.  &
                & (dim2 == SIZE(this%var(j)%p%ptr_3d, 2)) .AND.  &
                & (dim3 == SIZE(this%var(j)%p%ptr_3d, 3))) THEN
                var_shape = this%var(j)
              END IF
            END IF
            found = .TRUE.
            EXIT VARLOOP
          END IF
        END DO VARLOOP
        IF (.NOT. found) THEN ! error
          this%err_no = ERR_VARIABLE_NOT_FOUND
          RETURN
        END IF
      END SELECT
    END DO
  END SUBROUTINE expr_list_get_result_shape

  ! evaluate an expression, i.e. evaluate its list of operators.
  SUBROUTINE expr_list_evaluate_0D(this, evaluate)
    CLASS(expression), INTENT(INOUT) :: this
    REAL(real_kind), POINTER :: evaluate
    ! local variables
    INTEGER                    :: i, dim1, dim2, dim3
    TYPE(t_result_stack_0D)    :: result_stack
    TYPE(t_var_ptr)            :: var_shape
    REAL(real_kind), POINTER   :: ptr

    IF (this%isize <= 0) THEN
      this%err_no = ERR_EXPRESSION_EMPTY
      RETURN
    ELSE
      this%err_no = ERR_NONE
    END IF
    ! first step: determine the dimensions involved
    CALL this%get_result_shape(var_shape, dim1, dim2, dim3)

    IF ((dim1 /= 1) .OR. (dim2 /= 1) .OR. (dim3 /= 1)) THEN
      this%err_no = ERR_DIMENSION_MISMATCH
      RETURN
    END IF
    ! then evaluate the list of operators:
    DO i = 1, this%isize
      CALL this%list(i)%p%eval(result_stack)
    END DO
    IF (.NOT. ASSOCIATED(evaluate)) THEN
      IF (result_stack%ridx == 1) THEN
        evaluate => result_stack%rstack(1)%ptr_0D
      ELSE
        evaluate => NULL()
      END IF
    ELSE
      IF (result_stack%ridx == 1) THEN
        evaluate = result_stack%rstack(1)%ptr_0D
        CALL result_stack%pop(ptr)
        DEALLOCATE (ptr)
      END IF
    END IF
  END SUBROUTINE expr_list_evaluate_0D

  SUBROUTINE expr_list_evaluate_2D(this, evaluate)
    CLASS(expression), INTENT(INOUT) :: this
    REAL(real_kind), POINTER :: evaluate(:, :)
    ! local variables
    INTEGER                    :: i, dim1, dim2, dim3, size1, size2
    TYPE(t_result_stack_2D)    :: result_stack
    TYPE(t_var_ptr)            :: var_shape
    REAL(real_kind), DIMENSION(:, :), POINTER :: ptr

    IF (this%isize <= 0) THEN
      this%err_no = ERR_EXPRESSION_EMPTY
      RETURN
    ELSE
      this%err_no = ERR_NONE
    END IF
    ! first step: determine the dimensions involved
    CALL this%get_result_shape(var_shape, dim1, dim2, dim3)
    IF (dim3 /= 1) THEN
      this%err_no = ERR_DIMENSION_MISMATCH
      RETURN
    END IF
    result_stack = t_result_stack_2D(dim1, dim2)
    ! then evaluate the list of operators:
    DO i = 1, this%isize
      CALL this%list(i)%p%eval(result_stack)
    END DO

    IF (.NOT. ASSOCIATED(evaluate)) THEN
      IF (result_stack%ridx == 1) THEN
        evaluate => result_stack%rstack(1)%ptr_2D
      ELSE
        evaluate => NULL()
      END IF
    ELSE
      IF (result_stack%ridx == 1) THEN
        size1 = SIZE(evaluate, 1)
        size2 = SIZE(evaluate, 2)
        IF ((dim1 == 1) .AND. (dim2 == 1)) THEN
          evaluate(:, :) = result_stack%rstack(1)%ptr_2D(1, 1)
        ELSE
          IF ((dim1 == size1) .AND. (dim2 == size2)) THEN
            evaluate = result_stack%rstack(1)%ptr_2D
          ELSE
            this%err_no = ERR_DIMENSION_MISMATCH
            evaluate(:, :) = 0._real_kind
          END IF
        END IF
      END IF
      CALL result_stack%pop(ptr)
      DEALLOCATE (ptr)
    END IF
  END SUBROUTINE expr_list_evaluate_2D

  SUBROUTINE expr_list_evaluate_3D(this, evaluate)
    CLASS(expression), INTENT(INOUT) :: this
    REAL(real_kind), POINTER :: evaluate(:, :, :)
    ! local variables
    INTEGER                    :: i, dim1, dim2, dim3, size1, size2, size3
    TYPE(t_result_stack_3D)    :: result_stack
    TYPE(t_var_ptr)            :: var_shape
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: ptr

    IF (this%isize <= 0) THEN
      this%err_no = ERR_EXPRESSION_EMPTY
      RETURN
    ELSE
      this%err_no = ERR_NONE
    END IF
    ! first step: determine the dimensions involved
    CALL this%get_result_shape(var_shape, dim1, dim2, dim3)
    result_stack = t_result_stack_3D(dim1, dim2, dim3)
    ! then evaluate the list of operators:
    DO i = 1, this%isize
      CALL this%list(i)%p%eval(result_stack)
    END DO

    IF (.NOT. ASSOCIATED(evaluate)) THEN
      IF (result_stack%ridx == 1) THEN
        evaluate => result_stack%rstack(1)%ptr_3D
      ELSE
        evaluate => NULL()
      END IF
    ELSE
      IF (result_stack%ridx == 1) THEN
        size1 = SIZE(evaluate, 1)
        size2 = SIZE(evaluate, 2)
        size3 = SIZE(evaluate, 3)
        IF ((dim1 == 1) .AND. (dim2 == 1)) THEN
          evaluate(:, :, :) = result_stack%rstack(1)%ptr_3D(1, 1, 1)
        ELSE IF ((dim3 == 1) .AND. (dim1 == size1) .AND. (dim2 == size3)) THEN
          DO i = 1, size2
            evaluate(:, i, :) = result_stack%rstack(1)%ptr_3D(:, :, 1)
          END DO
        ELSE IF ((dim1 == size1) .AND. (dim2 == size2) .AND. (dim3 == size3)) THEN
          evaluate = result_stack%rstack(1)%ptr_3D
        ELSE
          this%err_no = ERR_DIMENSION_MISMATCH
          evaluate(:, :, :) = 0._real_kind
        END IF
      END IF
      CALL result_stack%pop(ptr)
      DEALLOCATE (ptr)
    END IF
    IF (result_stack%ridx == 1) evaluate => result_stack%rstack(1)%ptr_3D
  END SUBROUTINE expr_list_evaluate_3D

  ! resize variable list to size "new_size".
  ! new elements are not initialized
  SUBROUTINE expr_list_resize_varlist(this, new_size)
    CLASS(expression), INTENT(INOUT) :: this
    INTEGER, INTENT(IN)    :: new_size
    ! local variables
    TYPE(t_var_ptr), ALLOCATABLE :: tmp(:)

    IF (.NOT. ALLOCATED(this%var)) THEN
      ALLOCATE (this%var(new_size))
    ELSE
      ALLOCATE (tmp(new_size))
      tmp(1:SIZE(this%var)) = this%var
      CALL MOVE_ALLOC(tmp, this%var)
    END IF
  END SUBROUTINE expr_list_resize_varlist

  ! link a variable name string to a given variable
  SUBROUTINE expr_list_link_0D(this, string, variable)
    CLASS(expression), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN)    :: string
    REAL(real_kind), TARGET          :: variable
    ! local variables
    INTEGER :: i

    ! overwrite existing links with the same name
    DO i = 1, this%nvars
      IF (TRIM(this%var(i)%name) == TRIM(string)) EXIT
    END DO
    IF (i == (this%nvars + 1)) THEN
      IF (.NOT. ALLOCATED(this%var)) THEN
        CALL this%resize_varlist(1)
      ELSE
        CALL this%resize_varlist(2*SIZE(this%var))
      END IF
      this%nvars = this%nvars + 1
    END IF
    this%var(i)%name = TRIM(string)
    this%var(i)%p%ptr_0D => variable
  END SUBROUTINE expr_list_link_0D

  ! link a variable name string to a given variable
  SUBROUTINE expr_list_link_2D(this, string, variable)
    CLASS(expression), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN)    :: string
    REAL(real_kind), TARGET          :: variable(:, :)
    ! local variables
    INTEGER :: i

    ! overwrite existing links with the same name
    DO i = 1, this%nvars
      IF (TRIM(this%var(i)%name) == TRIM(string)) EXIT
    END DO
    IF (i == (this%nvars + 1)) THEN
      IF (.NOT. ALLOCATED(this%var)) THEN
        CALL this%resize_varlist(1)
      ELSE
        CALL this%resize_varlist(2*SIZE(this%var))
      END IF
      this%nvars = this%nvars + 1
    END IF
    this%var(i)%name = TRIM(string)
    this%var(i)%p%ptr_2D => variable
  END SUBROUTINE expr_list_link_2D

  ! link a variable name string to a given variable
  SUBROUTINE expr_list_link_3D(this, string, variable)
    CLASS(expression), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN)    :: string
    REAL(real_kind), TARGET          :: variable(:, :, :)
    ! local variables
    INTEGER :: i

    ! overwrite existing links with the same name
    DO i = 1, this%nvars
      IF (TRIM(this%var(i)%name) == TRIM(string)) EXIT
    END DO
    IF (i == (this%nvars + 1)) THEN
      IF (.NOT. ALLOCATED(this%var)) THEN
        CALL this%resize_varlist(1)
      ELSE
        CALL this%resize_varlist(2*SIZE(this%var))
      END IF
      this%nvars = this%nvars + 1
    END IF
    this%var(i)%name = TRIM(string)
    this%var(i)%p%ptr_3D => variable
  END SUBROUTINE expr_list_link_3D

  SUBROUTINE stack_op_value_eval_0D(this, result_stack)
    CLASS(t_op_value), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: ptr
    ALLOCATE (ptr)
    ptr = this%val
    CALL result_stack%push(ptr)
  END SUBROUTINE stack_op_value_eval_0D

  SUBROUTINE stack_op_value_eval_2D(this, result_stack)
    CLASS(t_op_value), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: ptr
    ALLOCATE (ptr(result_stack%dim1, result_stack%dim2))
    ptr(:, :) = this%val
    CALL result_stack%push(ptr)
  END SUBROUTINE stack_op_value_eval_2D

  SUBROUTINE stack_op_value_eval_3D(this, result_stack)
    CLASS(t_op_value), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: ptr
    ALLOCATE (ptr(result_stack%dim1, result_stack%dim2, result_stack%dim3))
    ptr(:, :, :) = this%val
    CALL result_stack%push(ptr)
  END SUBROUTINE stack_op_value_eval_3D

  SUBROUTINE stack_op_variable_eval_0D(this, result_stack)
    CLASS(t_op_variable), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: ptr
    IF (ASSOCIATED(this%var%p%ptr_0D)) THEN
      ALLOCATE (ptr)
      ptr = this%var%p%ptr_0D
      CALL result_stack%push(ptr)
    ELSE
      this%err_no = ERR_INTERNAL
    END IF
  END SUBROUTINE stack_op_variable_eval_0D

  SUBROUTINE stack_op_variable_eval_2D(this, result_stack)
    CLASS(t_op_variable), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: ptr
    INTEGER :: pdim1, pdim2
    ALLOCATE (ptr(result_stack%dim1, result_stack%dim2))
    IF (ASSOCIATED(this%var%p%ptr_0D)) THEN
      ptr = this%var%p%ptr_0D
    ELSE
      pdim1 = SIZE(this%var%p%ptr_2D, 1)
      pdim2 = SIZE(this%var%p%ptr_2D, 2)
      IF ((pdim1 == 1) .AND. (pdim2 == 1)) THEN
        ptr = this%var%p%ptr_2D(1, 1)
      ELSE IF ((pdim1 == result_stack%dim1) .AND. (pdim2 == result_stack%dim2)) THEN
        ptr = this%var%p%ptr_2D
      ELSE
        DEALLOCATE (ptr)
        this%err_no = ERR_DIMENSION_MISMATCH
        RETURN
      END IF
    END IF
    CALL result_stack%push(ptr)
  END SUBROUTINE stack_op_variable_eval_2D

  SUBROUTINE stack_op_variable_eval_3D(this, result_stack)
    CLASS(t_op_variable), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: ptr
    INTEGER :: pdim1, pdim2, pdim3, i
    ALLOCATE (ptr(result_stack%dim1, result_stack%dim2, result_stack%dim3))
    IF (ASSOCIATED(this%var%p%ptr_0D)) THEN
      ptr = this%var%p%ptr_0D
    ELSE IF (ASSOCIATED(this%var%p%ptr_2D)) THEN
      pdim1 = SIZE(this%var%p%ptr_2D, 1)
      pdim2 = SIZE(this%var%p%ptr_2D, 2)
      IF ((pdim1 == 1) .AND. (pdim2 == 1)) THEN
        ptr = this%var%p%ptr_2D(1, 1)
      ELSE IF ((pdim1 == result_stack%dim1) .AND. (pdim2 == result_stack%dim3)) THEN
        DO i = 1, result_stack%dim2
          ptr(:, i, :) = this%var%p%ptr_2D(:, :)
        END DO
      ELSE
        DEALLOCATE (ptr)
        this%err_no = ERR_DIMENSION_MISMATCH
        RETURN
      END IF
    ELSE IF (ASSOCIATED(this%var%p%ptr_3D)) THEN
      pdim1 = SIZE(this%var%p%ptr_3D, 1)
      pdim2 = SIZE(this%var%p%ptr_3D, 2)
      pdim3 = SIZE(this%var%p%ptr_3D, 3)
      IF ((pdim1 == 1) .AND. (pdim2 == 1) .AND. (pdim3 == 1)) THEN
        ptr = this%var%p%ptr_3D(1, 1, 1)
      ELSE IF ((pdim1 == result_stack%dim1) .AND. &
        &      (pdim2 == result_stack%dim2) .AND. &
        &      (pdim3 == result_stack%dim3)) THEN
        ptr(:, :, :) = this%var%p%ptr_3D(:, :, :)
      ELSE
        DEALLOCATE (ptr)
        this%err_no = ERR_DIMENSION_MISMATCH
        RETURN
      END IF
    END IF
    CALL result_stack%push(ptr)
  END SUBROUTINE stack_op_variable_eval_3D

  SUBROUTINE stack_op_plus_eval_0D(this, result_stack)
    CLASS(t_op_plus), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1 + arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_plus_eval_0D

  SUBROUTINE stack_op_plus_eval_2D(this, result_stack)
    CLASS(t_op_plus), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1 + arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_plus_eval_2D

  SUBROUTINE stack_op_plus_eval_3D(this, result_stack)
    CLASS(t_op_plus), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1 + arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_plus_eval_3D

  SUBROUTINE stack_op_minus_eval_0D(this, result_stack)
    CLASS(t_op_minus), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1 - arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_minus_eval_0D

  SUBROUTINE stack_op_minus_eval_2D(this, result_stack)
    CLASS(t_op_minus), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1 - arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_minus_eval_2D

  SUBROUTINE stack_op_minus_eval_3D(this, result_stack)
    CLASS(t_op_minus), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1 - arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_minus_eval_3D

  SUBROUTINE stack_op_mul_eval_0D(this, result_stack)
    CLASS(t_op_mul), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1*arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_mul_eval_0D

  SUBROUTINE stack_op_mul_eval_2D(this, result_stack)
    CLASS(t_op_mul), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1*arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_mul_eval_2D

  SUBROUTINE stack_op_mul_eval_3D(this, result_stack)
    CLASS(t_op_mul), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1*arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_mul_eval_3D

  SUBROUTINE stack_op_div_eval_0D(this, result_stack)
    CLASS(t_op_div), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1/arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_div_eval_0D

  SUBROUTINE stack_op_div_eval_2D(this, result_stack)
    CLASS(t_op_div), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1/arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_div_eval_2D

  SUBROUTINE stack_op_div_eval_3D(this, result_stack)
    CLASS(t_op_div), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1/arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_div_eval_3D

  SUBROUTINE stack_op_pow_eval_0D(this, result_stack)
    CLASS(t_op_pow), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1**arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_pow_eval_0D

  SUBROUTINE stack_op_pow_eval_2D(this, result_stack)
    CLASS(t_op_pow), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1**arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_pow_eval_2D

  SUBROUTINE stack_op_pow_eval_3D(this, result_stack)
    CLASS(t_op_pow), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = arg1**arg2
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_pow_eval_3D

  SUBROUTINE stack_op_gt_eval_0D(this, result_stack)
    CLASS(t_op_gt), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    IF (arg1 > arg2) THEN
      arg1 = 1._real_kind
    ELSE
      arg1 = 0._real_kind
    END IF
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_gt_eval_0D

  SUBROUTINE stack_op_gt_eval_2D(this, result_stack)
    CLASS(t_op_gt), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    WHERE (arg1(:, :) > arg2(:, :))
      arg1 = 1._real_kind
    ELSEWHERE
      arg1 = 0._real_kind
    END WHERE
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_gt_eval_2D

  SUBROUTINE stack_op_gt_eval_3D(this, result_stack)
    CLASS(t_op_gt), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    WHERE (arg1(:, :, :) > arg2(:, :, :))
      arg1 = 1._real_kind
    ELSEWHERE
      arg1 = 0._real_kind
    END WHERE
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_gt_eval_3D

  SUBROUTINE stack_op_lt_eval_0D(this, result_stack)
    CLASS(t_op_lt), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    IF (arg1 < arg2) THEN
      arg1 = 1._real_kind
    ELSE
      arg1 = 0._real_kind
    END IF
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_lt_eval_0D

  SUBROUTINE stack_op_lt_eval_2D(this, result_stack)
    CLASS(t_op_lt), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    WHERE (arg1(:, :) < arg2)
      arg1(:, :) = 1._real_kind
    ELSEWHERE
      arg1(:, :) = 0._real_kind
    END WHERE
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_lt_eval_2D

  SUBROUTINE stack_op_lt_eval_3D(this, result_stack)
    CLASS(t_op_lt), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    WHERE (arg1(:, :, :) < arg2(:, :, :))
      arg1 = 1._real_kind
    ELSEWHERE
      arg1 = 0._real_kind
    END WHERE
    CALL result_stack%push(arg1)
    DEALLOCATE (arg2)
  END SUBROUTINE stack_op_lt_eval_3D

  SUBROUTINE stack_op_exp_eval_0D(this, result_stack)
    CLASS(t_op_exp), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = EXP(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_exp_eval_0D

  SUBROUTINE stack_op_exp_eval_2D(this, result_stack)
    CLASS(t_op_exp), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = EXP(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_exp_eval_2D

  SUBROUTINE stack_op_exp_eval_3D(this, result_stack)
    CLASS(t_op_exp), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = EXP(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_exp_eval_3D

  SUBROUTINE stack_op_log_eval_0D(this, result_stack)
    CLASS(t_op_log), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = LOG(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_log_eval_0D

  SUBROUTINE stack_op_log_eval_2D(this, result_stack)
    CLASS(t_op_log), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = LOG(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_log_eval_2D

  SUBROUTINE stack_op_log_eval_3D(this, result_stack)
    CLASS(t_op_log), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = LOG(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_log_eval_3D

  SUBROUTINE stack_op_sin_eval_0D(this, result_stack)
    CLASS(t_op_sin), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = SIN(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_sin_eval_0D

  SUBROUTINE stack_op_sin_eval_2D(this, result_stack)
    CLASS(t_op_sin), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = SIN(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_sin_eval_2D

  SUBROUTINE stack_op_sin_eval_3D(this, result_stack)
    CLASS(t_op_sin), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = SIN(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_sin_eval_3D

  SUBROUTINE stack_op_cos_eval_0D(this, result_stack)
    CLASS(t_op_cos), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = COS(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_cos_eval_0D

  SUBROUTINE stack_op_cos_eval_2D(this, result_stack)
    CLASS(t_op_cos), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = COS(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_cos_eval_2D

  SUBROUTINE stack_op_cos_eval_3D(this, result_stack)
    CLASS(t_op_cos), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = COS(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_cos_eval_3D

  SUBROUTINE stack_op_min_eval_0D(this, result_stack)
    CLASS(t_op_min), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = MIN(arg1, arg2)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_min_eval_0D

  SUBROUTINE stack_op_min_eval_2D(this, result_stack)
    CLASS(t_op_min), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = MIN(arg1, arg2)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_min_eval_2D

  SUBROUTINE stack_op_min_eval_3D(this, result_stack)
    CLASS(t_op_min), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = MIN(arg1, arg2)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_min_eval_3D

  SUBROUTINE stack_op_max_eval_0D(this, result_stack)
    CLASS(t_op_max), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = MAX(arg1, arg2)
    DEALLOCATE (arg2)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_max_eval_0D

  SUBROUTINE stack_op_max_eval_2D(this, result_stack)
    CLASS(t_op_max), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = MAX(arg1, arg2)
    DEALLOCATE (arg2)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_max_eval_2D

  SUBROUTINE stack_op_max_eval_3D(this, result_stack)
    CLASS(t_op_max), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    arg1 = MAX(arg1, arg2)
    DEALLOCATE (arg2)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_max_eval_3D

  SUBROUTINE stack_op_if_eval_0D(this, result_stack)
    CLASS(t_op_if), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1, arg2, arg3
    CALL result_stack%pop(arg3)
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    IF (arg1 > 0.) THEN
      CALL result_stack%push(arg2)
      DEALLOCATE (arg3)
    ELSE
      CALL result_stack%push(arg3)
      DEALLOCATE (arg2)
    END IF
    DEALLOCATE (arg1)
  END SUBROUTINE stack_op_if_eval_0D

  SUBROUTINE stack_op_if_eval_2D(this, result_stack)
    CLASS(t_op_if), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1, arg2, arg3
    CALL result_stack%pop(arg3)
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    WHERE (arg1(:, :) <= 0._real_kind)
      arg2(:, :) = arg3(:, :)
    END WHERE
    CALL result_stack%push(arg2)
    DEALLOCATE (arg1, arg3)
  END SUBROUTINE stack_op_if_eval_2D

  SUBROUTINE stack_op_if_eval_3D(this, result_stack)
    CLASS(t_op_if), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1, arg2, arg3
    CALL result_stack%pop(arg3)
    CALL result_stack%pop(arg2)
    CALL result_stack%pop(arg1)
    WHERE (arg1(:, :, :) <= 0._real_kind)
      arg2(:, :, :) = arg3(:, :, :)
    END WHERE
    CALL result_stack%push(arg2)
    DEALLOCATE (arg1, arg3)
  END SUBROUTINE stack_op_if_eval_3D

  SUBROUTINE stack_op_sqrt_eval_0D(this, result_stack)
    CLASS(t_op_sqrt), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: p
    CALL result_stack%pop(p)
    p = SQRT(p)
    CALL result_stack%push(p)
  END SUBROUTINE stack_op_sqrt_eval_0D

  SUBROUTINE stack_op_sqrt_eval_2D(this, result_stack)
    CLASS(t_op_sqrt), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: p(:, :)
    CALL result_stack%pop(p)
    p = SQRT(p)
    CALL result_stack%push(p)
  END SUBROUTINE stack_op_sqrt_eval_2D

  SUBROUTINE stack_op_sqrt_eval_3D(this, result_stack)
    CLASS(t_op_sqrt), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: p(:, :, :)
    CALL result_stack%pop(p)
    p = SQRT(p)
    CALL result_stack%push(p)
  END SUBROUTINE stack_op_sqrt_eval_3D

  SUBROUTINE stack_op_erf_eval_0D(this, result_stack)
    CLASS(t_op_erf), INTENT(INOUT) :: this
    CLASS(t_result_stack_0D), INTENT(INOUT) :: result_stack
    REAL(real_kind), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = ERF(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_erf_eval_0D

  SUBROUTINE stack_op_erf_eval_2D(this, result_stack)
    CLASS(t_op_erf), INTENT(INOUT) :: this
    CLASS(t_result_stack_2D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = ERF(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_erf_eval_2D

  SUBROUTINE stack_op_erf_eval_3D(this, result_stack)
    CLASS(t_op_erf), INTENT(INOUT) :: this
    CLASS(t_result_stack_3D), INTENT(INOUT) :: result_stack
    REAL(real_kind), DIMENSION(:, :, :), POINTER :: arg1
    CALL result_stack%pop(arg1)
    arg1 = ERF(arg1)
    CALL result_stack%push(arg1)
  END SUBROUTINE stack_op_erf_eval_3D

  !> Convert c-style character array into Fortran string.
  !
  FUNCTION convert_c_string(length, c_str) RESULT(regular_string)
    USE ISO_C_BINDING, ONLY: c_char, c_null_char
    INTEGER, INTENT(IN) :: length
    CHARACTER(kind=c_char, len=1), DIMENSION(length), INTENT(IN) :: c_str
    CHARACTER(len=length) :: regular_string
    ! local variables
    INTEGER :: i

    regular_string = " "
    loop_string: DO i = 1, length
      IF (c_str(i) == c_null_char) THEN
        EXIT loop_string
      ELSE
        regular_string(i:i) = c_str(i)
      END IF
    END DO loop_string
  END FUNCTION convert_c_string

END MODULE mo_expression
