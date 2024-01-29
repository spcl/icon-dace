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


#include "omp_definitions.inc"

MODULE mo_variable

  USE mo_kind, ONLY : wp, sp
  USE mo_exception, ONLY: finish
  USE mo_fortran_tools, ONLY: init
  USE mo_math_types, ONLY : t_geographical_coordinates

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  TYPE t_variable
    ! Meta-data
    CHARACTER(LEN=:), ALLOCATABLE :: name
    CHARACTER(LEN=:), ALLOCATABLE :: units
    CHARACTER(LEN=:), ALLOCATABLE :: type_id   ! "bool", "int", "real", "single", "geocoord"
    INTEGER :: dim = 0
    INTEGER :: dims(5) = [ 1, 1, 1, 1, 1 ]
    INTEGER :: starts(5) = [ 1, 1, 1, 1, 1 ]
    LOGICAL :: l_opt     ! Whether optional or not
    LOGICAL :: bound = .false.

    LOGICAL, POINTER             :: l0d            => NULL()
    LOGICAL, POINTER, PRIVATE    :: l0d_alloc(:)   => NULL()
    LOGICAL, POINTER             :: l1d(:)         => NULL()

    INTEGER, POINTER             :: i0d            => NULL()
    INTEGER, POINTER, PRIVATE    :: i0d_alloc(:)   => NULL()
    INTEGER, POINTER             :: i1d(:)         => NULL()
    INTEGER, POINTER             :: i2d(:,:)       => NULL()
    INTEGER, POINTER             :: i3d(:,:,:)     => NULL()
    INTEGER, POINTER             :: i4d(:,:,:,:)   => NULL()
    INTEGER, POINTER             :: i5d(:,:,:,:,:) => NULL()

    REAL(wp), POINTER            :: r0d            => NULL()
    REAL(wp), POINTER, PRIVATE   :: r0d_alloc(:)   => NULL()
    REAL(wp), POINTER            :: r1d(:)         => NULL()
    REAL(wp), POINTER            :: r2d(:,:)       => NULL()
    REAL(wp), POINTER            :: r3d(:,:,:)     => NULL()
    REAL(wp), POINTER            :: r4d(:,:,:,:)   => NULL()
    REAL(wp), POINTER            :: r5d(:,:,:,:,:) => NULL()

    REAL(sp), POINTER            :: s2d(:,:)       => NULL()
    REAL(sp), POINTER            :: s3d(:,:,:)     => NULL()

    TYPE(t_geographical_coordinates ), POINTER :: gc2d(:,:) => NULL()
    TYPE(t_geographical_coordinates ), POINTER :: gc4d(:,:,:,:) => NULL()

 END TYPE t_variable

  INTERFACE t_variable
    MODULE PROCEDURE t_variable_constructor
  END INTERFACE t_variable
  
  INTERFACE bind_variable
    MODULE PROCEDURE bind_variable_l0d
    MODULE PROCEDURE bind_variable_l1d
    MODULE PROCEDURE bind_variable_i0d
    MODULE PROCEDURE bind_variable_i1d
    MODULE PROCEDURE bind_variable_i2d
    MODULE PROCEDURE bind_variable_i3d
    MODULE PROCEDURE bind_variable_i4d
    MODULE PROCEDURE bind_variable_i5d
    MODULE PROCEDURE bind_variable_r0d
    MODULE PROCEDURE bind_variable_r1d
    MODULE PROCEDURE bind_variable_r2d
    MODULE PROCEDURE bind_variable_r3d
    MODULE PROCEDURE bind_variable_r4d
    MODULE PROCEDURE bind_variable_r5d
    MODULE PROCEDURE bind_variable_s2d
    MODULE PROCEDURE bind_variable_s3d
    MODULE PROCEDURE bind_variable_gc2d
    MODULE PROCEDURE bind_variable_gc4d
  END INTERFACE bind_variable

  PUBLIC t_variable, bind_variable, allocate_variable, deallocate_variable, unbind_variable

  CHARACTER(len=*), PARAMETER :: modname = 'mo_variable'

CONTAINS

  ! FUNCTION t_variable_constructor(name, dim, d, units, type_id, l_opt, s) RESULT(tv)
  !   TYPE(t_variable) :: tv
  !   CHARACTER(LEN=*), INTENT(IN) :: name
  !   CHARACTER(LEN=*), INTENT(IN) :: units
  !   CHARACTER(LEN=*), INTENT(IN) :: type_id
  !   INTEGER, INTENT(IN) :: dim     ! Number of dimensions
  !   INTEGER, INTENT(IN) :: d(:)
  !   LOGICAL, INTENT(IN) :: l_opt
  !   INTEGER, OPTIONAL, INTENT(IN) :: s(:)
  !   tv%name  = name
  !   tv%units = units
  !   tv%dim   = dim
  !   tv%dims(:dim) = d(:dim)
  !   if (PRESENT(s)) THEN
  !     tv%starts(:dim) = s(:dim)
  !   ENDIF
  !   tv%l_opt = l_opt
  !   tv%type_id = type_id
  ! END FUNCTION t_variable_constructor

  FUNCTION t_variable_constructor(name, d, units, type_id, l_opt, s) RESULT(tv)
    TYPE(t_variable) :: tv
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(IN) :: units
    CHARACTER(LEN=*), INTENT(IN) :: type_id
    INTEGER, INTENT(IN) :: d(:)        ! Shape
    LOGICAL, OPTIONAL, INTENT(IN) :: l_opt
    INTEGER, OPTIONAL, INTENT(IN) :: s(:)

    INTEGER :: dim     ! Number of dimensions

    dim = SIZE(d)
    IF (dim == 1 .AND. d(1) == 0) dim = 0

    tv%name  = name
    tv%units = units
    tv%dim   = dim
    tv%dims(:dim) = d(:dim)
    if (PRESENT(s)) THEN
      tv%starts(:dim) = s(:dim)
      !$ACC ENTER DATA COPYIN(tv%starts)
    ENDIF
    IF (PRESENT(l_opt)) THEN
      tv%l_opt = l_opt
    ELSE
      tv%l_opt = .FALSE.
    END IF
    tv%type_id = type_id

  END FUNCTION t_variable_constructor

  ! Doesn't actually bind to v but only copies!
  SUBROUTINE bind_variable_l0d( tv, v )
    CLASS(t_variable), POINTER :: tv
    LOGICAL, TARGET :: v
    IF ( ASSOCIATED( tv%l0d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      ALLOCATE(tv%l0d_alloc(1))
      tv%l0d_alloc(1) = v
      tv%l0d => tv%l0d_alloc(1)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 0 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not a scalar"
    ENDIF
  END SUBROUTINE bind_variable_l0d

  ! SUBROUTINE bind_variable_l0d( tv, v )
  !   CLASS(t_variable), POINTER :: tv
  !   LOGICAL, TARGET :: v
  !   IF ( ASSOCIATED( tv%l0d ) ) THEN
  !     PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
  !   ELSE
  !     tv%l0d => v
  !     tv%bound = .true.
  !   ENDIF
  !   IF ( tv%dim /= 0 ) THEN
  !     PRINT *, "ERROR: ", TRIM(tv%name), " is not a scalar"
  !   ENDIF
  ! END SUBROUTINE bind_variable_l0d

  SUBROUTINE bind_variable_l1d( tv, v )
    CLASS(t_variable), POINTER :: tv
    LOGICAL, POINTER :: v(:)
    IF ( ASSOCIATED( tv%l1d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%l1d => v
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 1 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not one dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_l1d

  ! Doesn't actually bind to v but only copies!
  SUBROUTINE bind_variable_i0d( tv, v )
    CLASS(t_variable), POINTER :: tv
    INTEGER, TARGET :: v
    IF ( ASSOCIATED( tv%i0d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      ALLOCATE(tv%i0d_alloc(1))
      tv%i0d_alloc(1) = v
      !$ACC ENTER DATA COPYIN(tv%i0d_alloc)
      tv%i0d => tv%i0d_alloc(1)
      __acc_attach(tv%i0d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 0 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not a scalar"
    ENDIF
  END SUBROUTINE bind_variable_i0d

  SUBROUTINE bind_variable_i1d( tv, v )
    CLASS(t_variable), POINTER :: tv
    INTEGER, POINTER :: v(:)
    IF ( ASSOCIATED( tv%i1d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%i1d => v
      __acc_attach(tv%i1d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 1 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not one dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_i1d

  SUBROUTINE bind_variable_i2d( tv, v )
    CLASS(t_variable), POINTER :: tv
    INTEGER, POINTER :: v(:,:)
    IF ( ASSOCIATED( tv%i2d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%i2d => v
      __acc_attach(tv%i2d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 2 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not two dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_i2d

  SUBROUTINE bind_variable_i3d( tv, v )
    CLASS(t_variable), POINTER :: tv
    INTEGER, POINTER :: v(:,:,:)
    IF ( ASSOCIATED( tv%i3d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%i3d => v
      __acc_attach(tv%i3d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 3 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not three dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_i3d

  SUBROUTINE bind_variable_i4d( tv, v )
    CLASS(t_variable), POINTER :: tv
    INTEGER, POINTER :: v(:,:,:,:)
    IF ( ASSOCIATED( tv%i4d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%i4d => v
      __acc_attach(tv%i4d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 4 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not four dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_i4d

  SUBROUTINE bind_variable_i5d( tv, v )
    CLASS(t_variable), POINTER :: tv
    INTEGER, POINTER :: v(:,:,:,:,:)
    IF ( ASSOCIATED( tv%i5d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%i5d => v
      __acc_attach(tv%i5d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 5 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not five dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_i5d

  ! Doesn't actually bind to v but only copies!
  SUBROUTINE bind_variable_r0d( tv, v )
    CLASS(t_variable), POINTER :: tv
    REAL(wp), TARGET :: v
    IF ( ASSOCIATED( tv%r0d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      ALLOCATE(tv%r0d_alloc(1))
      tv%r0d_alloc(1) = v
      !$ACC ENTER DATA COPYIN(tv%r0d_alloc)
      tv%r0d => tv%r0d_alloc(1)
      __acc_attach(tv%r0d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 0 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not a scalar"
    ENDIF
  END SUBROUTINE bind_variable_r0d

  SUBROUTINE bind_variable_r1d( tv, v )
    CLASS(t_variable), POINTER :: tv
    REAL(wp), POINTER :: v(:)
    IF ( ASSOCIATED( tv%r1d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%r1d => v
      __acc_attach(tv%r1d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 1 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not one dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_r1d

  SUBROUTINE bind_variable_r2d( tv, v )
    CLASS(t_variable), POINTER :: tv
    REAL(wp), POINTER :: v(:,:)
    IF ( ASSOCIATED( tv%r2d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%r2d => v
      __acc_attach(tv%r2d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 2 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not two dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_r2d

  SUBROUTINE bind_variable_r3d( tv, v )
    CLASS(t_variable), POINTER :: tv
    REAL(wp), POINTER :: v(:,:,:)
    IF ( ASSOCIATED( tv%r3d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%r3d => v
      __acc_attach(tv%r3d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 3 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not three dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_r3d

  SUBROUTINE bind_variable_r4d( tv, v )
    CLASS(t_variable), POINTER :: tv
    REAL(wp), POINTER :: v(:,:,:,:)
    IF ( ASSOCIATED( tv%r4d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%r4d => v
      __acc_attach(tv%r4d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 4 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not four dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_r4d

  SUBROUTINE bind_variable_r5d( tv, v )
    CLASS(t_variable), POINTER :: tv
    REAL(wp), POINTER :: v(:,:,:,:,:)
    IF ( ASSOCIATED( tv%r5d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%r5d => v
      __acc_attach(tv%r5d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 5 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not five dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_r5d

  SUBROUTINE bind_variable_s2d( tv, v )
    CLASS(t_variable), POINTER :: tv
    REAL(sp), POINTER :: v(:,:)
    IF ( ASSOCIATED( tv%r2d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%s2d => v
      __acc_attach(tv%s2d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 2 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not two dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_s2d

  SUBROUTINE bind_variable_s3d( tv, v )
    CLASS(t_variable), POINTER :: tv
    REAL(sp), POINTER :: v(:,:,:)
    IF ( ASSOCIATED( tv%s3d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%s3d => v
      __acc_attach(tv%s3d)
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 3 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not three dimensional"
    ENDIF
    IF (.not. ASSOCIATED(v) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " array not associated"
      STOP
    ENDIF
    IF ( ANY(SHAPE(v) /= tv%dims(1:tv%dim)) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " must have size ", tv%dims(1:tv%dim), " but has size ", SHAPE(v)
    ENDIF
  END SUBROUTINE bind_variable_s3d

  SUBROUTINE bind_variable_gc2d( tv, v )
    TYPE(t_variable), INTENT(INOUT) :: tv
    TYPE(t_geographical_coordinates), POINTER :: v(:,:)
    IF ( ASSOCIATED( tv%gc2d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%gc2d => v
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 2 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not two dimensional"
    ENDIF
  END SUBROUTINE bind_variable_gc2d

  SUBROUTINE bind_variable_gc4d( tv, v )
    TYPE(t_variable), INTENT(INOUT) :: tv
    TYPE(t_geographical_coordinates), POINTER :: v(:,:,:,:)
    IF ( ASSOCIATED( tv%gc4d ) ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is already associated"
    ELSE
      tv%gc4d => v
      tv%bound = .true.
    ENDIF
    IF ( tv%dim /= 4 ) THEN
      PRINT *, "ERROR: ", TRIM(tv%name), " is not four dimensional"
    ENDIF
  END SUBROUTINE bind_variable_gc4d
  
  SUBROUTINE allocate_variable(tv)
    CLASS(t_variable), INTENT(INOUT) :: tv

    CHARACTER(len=*), PARAMETER :: routine = modname//':allocate_variable'

    IF (tv%bound) THEN
      ! variable has been bound to already allocated memory
      RETURN
    END IF

    SELECT CASE( tv%dim )
    CASE (0)
      IF (tv%type_id == "bool") THEN 
        ALLOCATE(tv%l0d_alloc(1)) ; tv%l0d => tv%l0d_alloc(1) ; tv%l0d = .FALSE.
        !$ACC ENTER DATA CREATE(tv%l0d_alloc)
        __acc_attach(tv%l0d)
      ENDIF
      IF (tv%type_id == "int") THEN
        ALLOCATE(tv%i0d_alloc(1)) ; tv%i0d => tv%i0d_alloc(1) ; tv%i0d = 0
        !$ACC ENTER DATA CREATE(tv%i0d_alloc)
        __acc_attach(tv%i0d)
      ENDIF
      IF (tv%type_id == "real") THEN 
        ALLOCATE(tv%r0d_alloc(1)) ; tv%r0d => tv%r0d_alloc(1) ; tv%r0d = 0._wp
        !$ACC ENTER DATA COPYIN(tv%r0d_alloc)
        __acc_attach(tv%r0d)
      ENDIF
    CASE (1)
      IF (tv%type_id == "int") THEN
        CALL finish(routine, 'not implemented yet for i1d')
        ! ALLOCATE(tv%i1d(tv%starts(1):tv%starts(1)-1+tv%dims(1)))
        ! !$ACC ENTER DATA CREATE(tv%i1d)
        ! !ICON_OMP PARALLEL
        ! CALL init(tv%i1d)
        ! !ICON_OMP END PARALLEL
        ! !$ACC UPDATE DEVICE(tv%i1d)
      END IF
      IF (tv%type_id == "real") THEN
        ALLOCATE(tv%r1d(tv%starts(1):tv%starts(1)-1+tv%dims(1)))
        !$ACC ENTER DATA CREATE(tv%r1d)
        !ICON_OMP PARALLEL
        CALL init(tv%r1d)
        !ICON_OMP END PARALLEL
      END IF
      IF (tv%type_id == "bool") THEN
        ALLOCATE(tv%l1d(tv%starts(1):tv%starts(1)-1+tv%dims(1)))
        !$ACC ENTER DATA CREATE(tv%l1d)
        tv%l1d = .FALSE.
        !$ACC UPDATE DEVICE(tv%l1d)
      END IF

    CASE (2)
      IF (tv%type_id == "int") THEN
        ALLOCATE(tv%i2d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2)))
        !$ACC ENTER DATA CREATE(tv%i2d)
        !ICON_OMP PARALLEL
        CALL init(tv%i2d)
        !ICON_OMP END PARALLEL
      END IF
      IF (tv%type_id == "real") THEN
        ALLOCATE(tv%r2d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2)))
        !$ACC ENTER DATA CREATE(tv%r2d)
        !ICON_OMP PARALLEL
        CALL init(tv%r2d)
        !ICON_OMP END PARALLEL
      END IF
      IF (tv%type_id == "geocoord") THEN
        ALLOCATE(tv%gc2d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2)))
      END IF

    CASE (3)
      IF (tv%type_id == "int") THEN
        ALLOCATE(tv%i3d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2),tv%starts(3):tv%starts(3)-1+tv%dims(3)))
        !$ACC ENTER DATA CREATE(tv%i3d)
        !ICON_OMP PARALLEL
        CALL init(tv%i3d)
        !ICON_OMP END PARALLEL
      END IF
      IF (tv%type_id == "real") THEN
        ALLOCATE(tv%r3d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2),tv%starts(3):tv%starts(3)-1+tv%dims(3)))
        !$ACC ENTER DATA CREATE(tv%r3d)
        !ICON_OMP PARALLEL
        CALL init(tv%r3d)
        !ICON_OMP END PARALLEL
      END IF

    CASE (4)
      IF (tv%type_id == "int") THEN
        ALLOCATE(tv%i4d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2), &
          &      tv%starts(3):tv%starts(3)-1+tv%dims(3),tv%starts(4):tv%starts(4)-1+tv%dims(4)))
        !$ACC ENTER DATA CREATE(tv%i4d)
        !ICON_OMP PARALLEL
        CALL init(tv%i4d)
        !ICON_OMP END PARALLEL
      END IF
      IF (tv%type_id == "real") THEN
        ALLOCATE(tv%r4d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2), &
          &      tv%starts(3):tv%starts(3)-1+tv%dims(3),tv%starts(4):tv%starts(4)-1+tv%dims(4)))
        !$ACC ENTER DATA CREATE(tv%r4d)
        !ICON_OMP PARALLEL
        CALL init(tv%r4d)
        !ICON_OMP END PARALLEL
      END IF
      IF (tv%type_id == "geocoord") THEN
        ALLOCATE(tv%gc4d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2), &
          &      tv%starts(3):tv%starts(3)-1+tv%dims(3),tv%starts(4):tv%starts(4)-1+tv%dims(4)))
      END IF

    CASE (5)
      IF (tv%type_id == "int") THEN
        ALLOCATE(tv%i5d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2), &
          &      tv%starts(3):tv%starts(3)-1+tv%dims(3),tv%starts(4):tv%starts(4)-1+tv%dims(4),        &
          &      tv%starts(5):tv%starts(5)-1+tv%dims(5)))
        !$ACC ENTER DATA CREATE(tv%i5d)
        !ICON_OMP PARALLEL
        CALL init(tv%i5d, 0)
        !ICON_OMP END PARALLEL
      END IF
      IF (tv%type_id == "real") THEN
        ALLOCATE(tv%r5d(tv%starts(1):tv%starts(1)-1+tv%dims(1),tv%starts(2):tv%starts(2)-1+tv%dims(2), &
          &      tv%starts(3):tv%starts(3)-1+tv%dims(3),tv%starts(4):tv%starts(4)-1+tv%dims(4),        &
          &      tv%starts(5):tv%starts(5)-1+tv%dims(5)))
        !$ACC ENTER DATA CREATE(tv%r5d)
        !ICON_OMP PARALLEL
        CALL init(tv%r5d, 0._wp)
        !ICON_OMP END PARALLEL
      END IF
    CASE DEFAULT
      PRINT *, "allocation: dimension ", tv%dim, " not currently supported"     
    END SELECT

    tv%bound = .true.

  END SUBROUTINE allocate_variable

  SUBROUTINE deallocate_variable(tv)
    TYPE(t_variable), INTENT(INOUT) :: tv

    SELECT CASE( tv%dim )
    CASE (0)
      IF (tv%type_id == "bool") THEN 
        DEALLOCATE(tv%l0d_alloc) ; NULLIFY( tv%l0d )
      ENDIF
      IF (tv%type_id == "int") THEN 
        DEALLOCATE(tv%i0d_alloc) ; NULLIFY( tv%i0d )
      ENDIF
      IF (tv%type_id == "real") THEN
        DEALLOCATE(tv%r0d_alloc) ; NULLIFY( tv%r0d )
      ENDIF
    CASE (1)
      IF (tv%type_id == "int") DEALLOCATE(tv%i1d)
      IF (tv%type_id == "real")  DEALLOCATE(tv%r1d)
      IF (tv%type_id == "bool")  DEALLOCATE(tv%l1d)

    CASE (2)
      IF (tv%type_id == "int") DEALLOCATE(tv%i2d)
      IF (tv%type_id == "real") DEALLOCATE(tv%r2d)
      IF (tv%type_id == "geocoord") DEALLOCATE(tv%gc2d)

    CASE (3)
      IF (tv%type_id == "int") DEALLOCATE(tv%i3d)
      IF (tv%type_id == "real") DEALLOCATE(tv%r3d)

    CASE (4)
      IF (tv%type_id == "int") DEALLOCATE(tv%i4d)
      IF (tv%type_id == "real") DEALLOCATE(tv%r4d)
      IF (tv%type_id == "geocoord") DEALLOCATE(tv%gc4d)

    CASE (5)
      IF (tv%type_id == "int") DEALLOCATE(tv%i5d)
      IF (tv%type_id == "real") DEALLOCATE(tv%r5d)
    CASE DEFAULT
      PRINT *, "deallocation: dimension ", tv%dim, " not currently supported"     
    END SELECT
    tv%bound = .false.
  END SUBROUTINE deallocate_variable

  SUBROUTINE unbind_variable(tv)
    ! TYPE(t_variable), INTENT(INOUT) :: tv
    TYPE(t_variable) :: tv

    IF (.NOT. tv%bound) print*, "unbind_variable: "//tv%name//" not bound!"

    SELECT CASE( tv%dim )
    CASE (0)
      IF (tv%type_id == "bool") NULLIFY(tv%l0d)
      IF (tv%type_id == "int") NULLIFY(tv%i0d)
      IF (tv%type_id == "real") NULLIFY(tv%r0d)

    CASE (1)
      IF (tv%type_id == "int") NULLIFY(tv%i1d)
      IF (tv%type_id == "real")  NULLIFY(tv%r1d)
      IF (tv%type_id == "bool")  NULLIFY(tv%l1d)

    CASE (2)
      IF (tv%type_id == "int") NULLIFY(tv%i2d)
      IF (tv%type_id == "real") NULLIFY(tv%r2d)
      IF (tv%type_id == "geocoord") NULLIFY(tv%gc2d)

    CASE (3)
      IF (tv%type_id == "int") NULLIFY(tv%i3d)
      IF (tv%type_id == "real") NULLIFY(tv%r3d)

    CASE (4)
      IF (tv%type_id == "int") NULLIFY(tv%i4d)
      IF (tv%type_id == "real") NULLIFY(tv%r4d)

    CASE (5)
      IF (tv%type_id == "int") NULLIFY(tv%i5d)
      IF (tv%type_id == "real") NULLIFY(tv%r5d)
    CASE DEFAULT
      PRINT *, "unbind: dimension ", tv%dim, " not currently supported"     
    END SELECT
    tv%bound = .false.
  END SUBROUTINE unbind_variable

END MODULE mo_variable
