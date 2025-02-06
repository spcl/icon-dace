!> Contains basic definitions for variable class
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
MODULE mo_jsb_var_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp, dp
  USE mo_exception,         ONLY: finish, message
  USE mo_jsb_varlist_iface, ONLY: VARNAME_LEN
  USE mo_jsb_subset,        ONLY: t_subset, ON_DOMAIN, ON_CHUNK

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var, t_jsb_var_p, REAL1D, REAL2D, REAL3D

  ENUM, BIND(C)
    ENUMERATOR :: REAL1D=1, REAL2D, REAL3D
  END ENUM

  !================================================================================================================================
  ! t_jsb_var
  !================================================================================================================================
  TYPE, ABSTRACT :: t_jsb_var
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr1d(:) => NULL()
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr2d(:,:) => NULL()
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr3d(:,:,:) => NULL()
    CHARACTER(len=VARNAME_LEN)    :: name                      ! Name of variable
    CHARACTER(len=VARNAME_LEN)    :: full_name                 ! Full name of variable in varlist, including prefix and suffix
    CHARACTER(len=:), ALLOCATABLE :: unit
    INTEGER                       :: type                      ! REAL1D / REAL2D / REAL3D
    INTEGER                       :: subset_type               ! ON_DOMAIN / ON_CHUNK
    CHARACTER(len=:), ALLOCATABLE :: element_name              !
    CHARACTER(len=:), ALLOCATABLE :: element_shortname         ! Short name of element
    INTEGER                       :: element_id                !
    INTEGER                       :: owner_model_id            !
    INTEGER                       :: owner_proc_id = -1        !
    INTEGER, ALLOCATABLE          :: owner_tile_path(:)        !
    INTEGER, ALLOCATABLE          :: child_idx(:)
      !< Position-numbers of the same var (as this var) on the children of the current tile (i.e. within mem%vars(:) of that child).
    REAL(dp)                      :: missval
    LOGICAL                       :: l_aggregate_all = .FALSE. ! When aggregating, consider fractions of all child tiles even
                                                               ! of those where variable doesn't exist
    LOGICAL                       :: is_conserved_quan = .FALSE. ! Whether variable carries a to-be-conserved quantity
    INTEGER                       :: cons_quan_type_id = -1    ! id identifying the type of the conserved quantity
    LOGICAL                       :: is_in_output  = .FALSE.   ! Whether variable is requested to be written to output
    LOGICAL                       :: is_in_restart = .FALSE.   ! Whether variable is requested to be written to restart file
  CONTAINS
    PROCEDURE :: Get_subset         => t_jsb_var_get_subset         ! get chunk and block information for the var and thread
    PROCEDURE :: force_finalization => t_jsb_var_force_finalization ! finalization after operator overload avoiding memory leaks
    PROCEDURE(associate_pointers_iface), DEFERRED :: Associate_pointers
  END TYPE t_jsb_var

  ABSTRACT INTERFACE
    SUBROUTINE Associate_pointers_iface(this,ic_start,ic_end,iblk_start,iblk_end)
      IMPORT t_jsb_var
      CLASS(t_jsb_var), INTENT(inout), TARGET :: this
      INTEGER,          INTENT(in)    :: ic_start
      INTEGER,          INTENT(in)    :: ic_end
      INTEGER,          INTENT(in)    :: iblk_start
      INTEGER,          INTENT(in)    :: iblk_end
    END SUBROUTINE
  END INTERFACE

!#pragma GCC diagnostic push
!#pragma GCC diagnostic ignored "-Wsurprising"

  !================================================================================================================================
  ! t_jsb_var_real1d
  !================================================================================================================================
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real1d
    REAL(wp), POINTER &
#ifdef __ICON__
      :: ptr(:) => NULL()
#else
      ! 1d stream elements are not supported with ECHAM
      :: ptr(:,:) => NULL()
#endif
  CONTAINS
    PROCEDURE :: Associate_pointers => t_jsb_var_real1d_associate_pointers
  END TYPE t_jsb_var_real1d

  !================================================================================================================================
  ! t_jsb_var_real2d
  !================================================================================================================================
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
      :: ptr(:,:) => NULL(), &    ! associated with ptr2d(:,:)
      &  ptrm1(:) => NULL()       ! points to ptr(ics:ice,iblk) updated each timestep with Associate_pointers
  CONTAINS
    PROCEDURE               :: Associate_pointers => t_jsb_var_real2d_associate_pointers
    !! operator overload
    ! +
    PROCEDURE                 :: Add_real2d_and_scalar            !< var_real2d = var_real2d + scalar     | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(var2)     :: Add_scalar_and_real2d            !< var_real2d = scalar     + var_real2d | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_real2d_and_field1d           !< var_real2d = var_real2d + field1d    | ON_CHUNK
    PROCEDURE, PASS(var2)     :: Add_field1d_and_real2d           !< var_real2d = field1d    + var_real2d | ON_CHUNK
    PROCEDURE                 :: Add_real2d_and_field2d           !< var_real2d = var_real2d + field2d    | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Add_real2d_and_real2d            !< var_real2d = var_real2d + var_real2d | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(+) => Add_real2d_and_scalar,  Add_scalar_and_real2d, &
      &                                         Add_real2d_and_field1d, Add_field1d_and_real2d, &
      &                                         Add_real2d_and_field2d, &
      &                                         Add_real2d_and_real2d
    ! -
    PROCEDURE                 :: Subtract_real2d_minus_scalar     !< var_real2d = var_real2d - scalar     | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_real2d_minus_field1d    !< var_real2d = var_real2d - field1d    | ON_CHUNK
    PROCEDURE                 :: Subtract_real2d_minus_field2d    !< var_real2d = var_real2d - field2d    | ON_CHUNK & ON_DOMAIN ! field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Subtract_real2d_minus_real2d     !< var_real2d = var_real2d - var_real2d | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(-) => Subtract_real2d_minus_scalar, &
      &                                         Subtract_real2d_minus_field1d, &
      &                                         Subtract_real2d_minus_field2d, &
      &                                         Subtract_real2d_minus_real2d
    ! *
    PROCEDURE                 :: Multiply_real2d_with_scalar      !< var_real2d = var_real2d * scalar     | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(var2)     :: Multiply_scalar_with_real2d      !< var_real2d = scalar     * var_real2d | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_real2d_with_field1d     !< var_real2d = var_real2d * field1d    | ON_CHUNK
    PROCEDURE, PASS(var2)     :: Multiply_field1d_with_real2d     !< var_real2d = field1d    * var_real2d | ON_CHUNK
    PROCEDURE                 :: Multiply_real2d_with_field2d     !< var_real2d = var_real2d * field2d    | ON_CHUNK & ON_DOMAIN ! field2d(:,1) & field2d(:,:)
    PROCEDURE, PASS(var2)     :: Multiply_field2d_with_real2d     !< var_real2d = field2d    * var_real2d | ON_CHUNK & ON_DOMAIN ! field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Multiply_real2d_with_real2d      !< var_real2d = var_real2d * var_real2d | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(*) => Multiply_real2d_with_scalar,  Multiply_scalar_with_real2d, &
      &                                         Multiply_real2d_with_field1d, Multiply_field1d_with_real2d, &
      &                                         Multiply_real2d_with_field2d, Multiply_field2d_with_real2d, &
      &                                         Multiply_real2d_with_real2d
    ! /
    PROCEDURE                 :: Div_real2d_by_scalar             !< var_real2d = var_real2d / scalar     | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_real2d_by_field1d            !< var_real2d = var_real2d / field1d    | ON_CHUNK
    PROCEDURE                 :: Div_real2d_by_field2d            !< var_real2d = var_real2d / field2d    | ON_CHUNK & ON_DOMAIN ! field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Div_real2d_by_real2d             !< var_real2d = var_real2d / var_real2d | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(/) => Div_real2d_by_scalar, &
      &                                         Div_real2d_by_field1d, &
      &                                         Div_real2d_by_field2d, &
      &                                         Div_real2d_by_real2d
    ! =
    PROCEDURE                 :: Assign_scalar_to_var_real2d      !< var_real2d = scalar       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_field1d_to_var_real2d     !< var_real2d = field1d(:)   | ON_CHUNK
    PROCEDURE                 :: Assign_field2d_to_var_real2d     !< var_real2d = field2d(:,:) | ON_CHUNK & ON_DOMAIN ! field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Assign_var_real2d_to_var_real2d  !< var_real2d = var_real2d   | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: ASSIGNMENT(=) => Assign_scalar_to_var_real2d, &
      &                                           Assign_field1d_to_var_real2d, &
      &                                           Assign_field2d_to_var_real2d, &
      &                                           Assign_var_real2d_to_var_real2d
    ! finalize
    FINAL                     :: finalize_jsb_var_real2d
  END TYPE t_jsb_var_real2d

  !================================================================================================================================
  ! t_jsb_var_real3d
  !================================================================================================================================
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d
    REAL(wp), POINTER :: ptr(:,:,:) => NULL()     ! associated with ptr3d(:,:,:)
    REAL(wp), POINTER :: ptrm1(:,:) => NULL()     ! points to ptr(ics:ice,:,iblk) updated each timestep with Associate_pointers
  CONTAINS
    PROCEDURE :: Associate_pointers => t_jsb_var_real3d_associate_pointers
    !! operator overload
    ! +
    PROCEDURE                 :: Add_real3d_and_scalar            !< var_real3d = var_real3d + scalar     | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(var2)     :: Add_scalar_and_real3d            !< var_real3d = scalar     + var_real3d | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_real3d_and_field2d           !< var_real3d = var_real3d + field2d    | ON_CHUNK             ! field2d(nc,vgrid)
    PROCEDURE, PASS(var2)     :: Add_field2d_and_real3d           !< var_real3d = field2d    + var_real3d | ON_CHUNK             ! field2d(nc,vgrid)
    PROCEDURE                 :: Add_real3d_and_field3d           !< var_real3d = var_real3d + field3d    | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_real3d_and_real3d            !< var_real3d = var_real3d + var_real3d | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(+) => Add_real3d_and_scalar,  Add_scalar_and_real3d, &
                                                Add_real3d_and_field2d, Add_field2d_and_real3d, &
                                                Add_real3d_and_field3d, &
                                                Add_real3d_and_real3d
    ! -
    PROCEDURE                 :: Subtract_real3d_minus_scalar     !< var_real3d = var_real3d - scalar     | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_real3d_minus_field2d    !< var_real3d = var_real3d - field2d    | ON_CHUNK             ! field2d(nc,vgrid)
    PROCEDURE                 :: Subtract_real3d_minus_field3d    !< var_real3d = var_real3d - field3d    | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_real3d_minus_real3d     !< var_real3d = var_real3d - var_real3d | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(-) => Subtract_real3d_minus_scalar, &
                                                Subtract_real3d_minus_field2d, &
                                                Subtract_real3d_minus_field3d, &
                                                Subtract_real3d_minus_real3d
    ! *
    PROCEDURE                 :: Multiply_real3d_with_scalar      !< var_real3d = var_real3d * scalar     | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(var2)     :: Multiply_scalar_with_real3d      !< var_real3d = scalar     * var_real3d | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_real3d_with_field2d     !< var_real3d = var_real3d * field2d    | ON_CHUNK             ! field2d(nc,vgrid)
    PROCEDURE, PASS(var2)     :: Multiply_field2d_with_real3d     !< var_real3d = field2d    * var_real3d | ON_CHUNK             ! field2d(nc,vgrid)
    PROCEDURE                 :: Multiply_real3d_with_field3d     !< var_real3d = var_real3d * field3d    | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(var2)     :: Multiply_field3d_with_real3d     !< var_real3d = field3d    * var_real3d | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_real3d_with_real3d      !< var_real3d = var_real3d * var_real3d | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(*) => Multiply_real3d_with_scalar,  Multiply_scalar_with_real3d, &
                                                Multiply_real3d_with_field2d, Multiply_field2d_with_real3d, &
                                                Multiply_real3d_with_field3d, Multiply_field3d_with_real3d, &
                                                Multiply_real3d_with_real3d
    ! /
    PROCEDURE                 :: Div_real3d_by_scalar             !< var_real3d = var_real3d / scalar     | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_real3d_by_field2d            !< var_real3d = var_real3d / field2d    | ON_CHUNK             ! field2d(nc,vgrid)
    PROCEDURE                 :: Div_real3d_by_field3d            !< var_real3d = var_real3d / field3d    | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_real3d_by_real3d             !< var_real3d = var_real3d / var_real3d | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(/) => Div_real3d_by_scalar, &
                                                Div_real3d_by_field2d, &
                                                Div_real3d_by_field3d, &
                                                Div_real3d_by_real3d
    ! =
    PROCEDURE                 :: Assign_scalar_to_var_real3d      !< var_real3d = scalar         | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_field2d_to_var_real3d     !< var_real3d = field2d(:,:)   | ON_CHUNK field2d(:,:)   | ON_DOMAIN / finish /
    PROCEDURE                 :: Assign_field3d_to_var_real3d     !< var_real3d = field3d(:,:,:) | ON_CHUNK field3d(:,:,1) | ON_DOMAIN field3d(:,:,:)
    PROCEDURE                 :: Assign_var_real3d_to_var_real3d  !< var_real3d = var_real3d     | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: ASSIGNMENT(=) => Assign_scalar_to_var_real3d, &
      &                                           Assign_field2d_to_var_real3d, &
      &                                           Assign_field3d_to_var_real3d, &
      &                                           Assign_var_real3d_to_var_real3d
    ! finalize
    FINAL     :: finalize_jsb_var_real3d
  END TYPE t_jsb_var_real3d

!#pragma GCC diagnostic pop

  !================================================================================================================================
  ! t_jsb_var_p
  !================================================================================================================================
  TYPE t_jsb_var_p
    CLASS(t_jsb_var), POINTER :: p => NULL()
  END TYPE t_jsb_var_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_var_class'

CONTAINS

  !================================================================================================================================
  ! get subset
  ! chunk and block information per model & thread from mo_jsb_subset:jsbach_subsets(model_id)%sub(thread)
  !================================================================================================================================
  FUNCTION t_jsb_var_get_subset(this) RESULT(subset)

    USE mo_jsb_subset,   ONLY: jsbach_subsets
    USE mo_jsb_parallel, ONLY: Get_omp_thread

    CLASS(t_jsb_var), INTENT(in) :: this
    TYPE(t_subset)               :: subset

    ! print*,'FFF ',TRIM(this%full_name), this%owner_model_id, Get_omp_thread()
    subset = jsbach_subsets(this%owner_model_id)%sub(Get_omp_thread())
    ! print*,'GGG ',subset%type,subset%ics,subset%ice,subset%iblk

  END FUNCTION t_jsb_var_get_subset

  !================================================================================================================================
  ! associate "dimension-reduced" pointers ptrm1
  ! var_real2d & var_real3d
  ! called for each timestep in the mo_jsb_interface "CALL model%Associate_var_pointers(ics, ice, iblk, iblk)"
  !================================================================================================================================
  SUBROUTINE t_jsb_var_real1d_associate_pointers(this,ic_start,ic_end,iblk_start,iblk_end)

    CLASS(t_jsb_var_real1d), INTENT(inout), TARGET :: this
    INTEGER,                 INTENT(in)    :: ic_start
    INTEGER,                 INTENT(in)    :: ic_end
    INTEGER,                 INTENT(in)    :: iblk_start
    INTEGER,                 INTENT(in)    :: iblk_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_var_real1d_associate_pointers'

    ! Do nothing for 1d variables

  END SUBROUTINE t_jsb_var_real1d_associate_pointers

  SUBROUTINE t_jsb_var_real2d_associate_pointers(this,ic_start,ic_end,iblk_start,iblk_end)

    CLASS(t_jsb_var_real2d), INTENT(inout), TARGET :: this
    INTEGER,                 INTENT(in)    :: ic_start
    INTEGER,                 INTENT(in)    :: ic_end
    INTEGER,                 INTENT(in)    :: iblk_start
    INTEGER,                 INTENT(in)    :: iblk_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_var_real2d_associate_pointers'

    IF (ASSOCIATED(this%ptrm1)) NULLIFY(this%ptrm1)
    IF (iblk_start == iblk_end) THEN
      this%ptrm1 => this%ptr2d(ic_start:ic_end, iblk_start)
    END IF

  END SUBROUTINE t_jsb_var_real2d_associate_pointers

  SUBROUTINE t_jsb_var_real3d_associate_pointers(this,ic_start,ic_end,iblk_start,iblk_end)

    CLASS(t_jsb_var_real3d), INTENT(inout), TARGET :: this
    INTEGER,                 INTENT(in)    :: ic_start
    INTEGER,                 INTENT(in)    :: ic_end
    INTEGER,                 INTENT(in)    :: iblk_start
    INTEGER,                 INTENT(in)    :: iblk_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_var_real3d_associate_pointers'

    IF (ASSOCIATED(this%ptrm1)) NULLIFY(this%ptrm1)
    IF (iblk_start == iblk_end) THEN
      this%ptrm1 => this%ptr3d(ic_start:ic_end, :, iblk_start)
    END IF

  END SUBROUTINE t_jsb_var_real3d_associate_pointers

  !================================================================================================================================
  ! force finalization
  !================================================================================================================================
  SUBROUTINE t_jsb_var_force_finalization(this)

    CLASS(t_jsb_var), INTENT(inout) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_var_force_finalization'

  END SUBROUTINE t_jsb_var_force_finalization

  !================================================================================================================================
  ! finalize after operator overload
  ! var_real2d & var_real3d
  !================================================================================================================================
  SUBROUTINE finalize_jsb_var_real2d(this)

    TYPE(t_jsb_var_real2d) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_jsb_var_real2d'

    !CALL message(routine, 'Finalizing '//TRIM(this%name))

    IF (ALLOCATED(this%owner_tile_path)) DEALLOCATE(this%owner_tile_path)
    IF (ALLOCATED(this%child_idx)) DEALLOCATE(this%child_idx)

    IF (ASSOCIATED(this%ptr))   NULLIFY(this%ptr)
    IF (ASSOCIATED(this%ptr2d)) DEALLOCATE(this%ptr2d)

  END SUBROUTINE finalize_jsb_var_real2d

  SUBROUTINE finalize_jsb_var_real3d(this)

    TYPE(t_jsb_var_real3d) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_jsb_var_real3d'

    !CALL message(routine, 'Finalizing '//TRIM(this%name))

    IF (ALLOCATED(this%owner_tile_path)) DEALLOCATE(this%owner_tile_path)
    IF (ALLOCATED(this%child_idx)) DEALLOCATE(this%child_idx)

    IF (ASSOCIATED(this%ptr))   NULLIFY(this%ptr)
    IF (ASSOCIATED(this%ptr3d)) DEALLOCATE(this%ptr3d)

  END SUBROUTINE finalize_jsb_var_real3d

  !================================================================================================================================
  ! create var_real2d/3d
  ! used as temporary variable for operator overload functions ?
  !================================================================================================================================
!!$  FUNCTION Create_jsb_var(i_shape, i_type) RESULT(return_ptr)
!!$
!!$    INTEGER, INTENT(in) :: i_shape(2)
!!$    INTEGER, INTENT(in), OPTIONAL :: i_type
!!$    TYPE(t_jsb_var), POINTER     :: return_ptr
!!$
!!$    ALLOCATE(return_ptr)
!!$
!!$    IF (SIZE(i_shape) == 2) THEN
!!$      return_value%type = REAL2D
!!$      ELSE IF (SIZE(i_shape) == 3) THEN
!!$        return_value%type = REAL3D
!!$      END IF
!!$    END IF
!!$
!!$    SELECT CASE (return_value%type)
!!$    CASE (REAL2D)
!!$      ALLOCATE(return_value%real2d_ptr(i_shape(1), i_shape(2)))
!!$    CASE (REAL3D)
!!$      ALLOCATE(return_value%real3d_ptr(i_shape(1), i_shape(2), i_shape(3)))
!!$    END SELECT
!!$
!!$  END FUNCTION Create_jsb_var

  !================================================================================================================================
  ! assign
  !================================================================================================================================
  !-----------------------------------------------------------------------------------------------------
  ! 2D
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d = scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_scalar_to_var_real2d(lhs, rhs)

    CLASS(t_jsb_var_real2d), INTENT(inout) :: lhs
    REAL(wp),                INTENT(in)    :: rhs

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_var_real2d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      lhs%ptr2d(ics:ice,iblk) = rhs
    ELSE
      lhs%ptr2d(:,:) = rhs
    END IF

  END SUBROUTINE Assign_scalar_to_var_real2d
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d = field1d(:)
  !   ON_CHUNK   field1d(:)
  !   ON_DOMAIN  / finish /
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field1d_to_var_real2d(lhs, rhs)

    CLASS(t_jsb_var_real2d), INTENT(inout) :: lhs
    REAL(wp),                INTENT(in)    :: rhs(:)

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field1d_to_var_real2d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      lhs%ptr2d(ics:ice,iblk) = rhs(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot assign 1d-field to jsb_var real2d on domain')
    END IF

  END SUBROUTINE Assign_field1d_to_var_real2d
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d = field2d(:,:)
  !   ON_CHUNK   field2d(:,1)
  !   ON_DOMAIN  field2d(:,:)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field2d_to_var_real2d(lhs, rhs)

    CLASS(t_jsb_var_real2d), INTENT(inout) :: lhs
    REAL(wp),                INTENT(in)    :: rhs(:,:)

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_var_real2d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(rhs,2) > 1) CALL finish(routine, 'Dimension mismatch, field with nblks > 1, assignment to jsb_var real2d on chunk')
      lhs%ptr2d(ics:ice,iblk) = rhs(:,1)
    ELSE
      lhs%ptr2d(:,:) = rhs(:,:)
    END IF

  END SUBROUTINE Assign_field2d_to_var_real2d
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d = var_real2d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_var_real2d_to_var_real2d(lhs, rhs)

    CLASS(t_jsb_var_real2d), INTENT(inout) :: lhs
    CLASS(t_jsb_var_real2d), INTENT(in)    :: rhs

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_var_real2d_to_var_real2d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      lhs%ptr2d(ics:ice,iblk) = rhs%ptr2d(ics:ice,iblk)
    ELSE
      lhs%ptr2d(:,:) = rhs%ptr2d(:,:)
    END IF

  END SUBROUTINE Assign_var_real2d_to_var_real2d

  !-----------------------------------------------------------------------------------------------------
  ! 3D
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d = scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_scalar_to_var_real3d(lhs, rhs)

    CLASS(t_jsb_var_real3d), INTENT(inout) :: lhs
    REAL(wp),                INTENT(in)    :: rhs

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_var_real3d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      lhs%ptr3d(ics:ice,:,iblk) = rhs
    ELSE
      lhs%ptr3d(:,:,:) = rhs
    END IF

  END SUBROUTINE Assign_scalar_to_var_real3d
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d = field2d(:,:) | field2d(chunk,vgrid)
  !   ON_CHUNK   field2d(:,:)
  !   ON_DOMAIN  / finish /
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field2d_to_var_real3d(lhs, rhs)

    CLASS(t_jsb_var_real3d), INTENT(inout) :: lhs
    REAL(wp),                INTENT(in)    :: rhs(:,:)

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_var_real3d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(lhs%ptr3d,2) /= SIZE(rhs,2)) CALL finish(routine, 'Dimension mismatch, vgrid size differs.')
      lhs%ptr3d(ics:ice,:,iblk) = rhs(:,:)
    ELSE
      CALL finish(routine, 'Cannot assign 2d-field to jsb_var real3d on domain')
    END IF

  END SUBROUTINE Assign_field2d_to_var_real3d
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d = field3d(:,:,:) | field2d(chunk,vgrid,nblks)
  !   ON_CHUNK   field3d(:,:,1)
  !   ON_DOMAIN  field3d(:,:,:)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field3d_to_var_real3d(lhs, rhs)

    CLASS(t_jsb_var_real3d), INTENT(inout) :: lhs
    REAL(wp),                INTENT(in)    :: rhs(:,:,:)

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field3d_to_var_real3d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (SIZE(lhs%ptr3d,2) /= SIZE(rhs,2)) CALL finish(routine, 'Dimension mismatch, vgrid size differs.')
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(rhs,3) > 1) CALL finish(routine, 'Dimension mismatch, field with nblks > 1, assignment to jsb_var real3d on chunk')
      lhs%ptr3d(ics:ice,:,iblk) = rhs(:,:,1)
    ELSE
      lhs%ptr3d(:,:,:) = rhs(:,:,:)
    END IF

  END SUBROUTINE Assign_field3d_to_var_real3d
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d = var_real3d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_var_real3d_to_var_real3d(lhs, rhs)

    CLASS(t_jsb_var_real3d), INTENT(inout) :: lhs
    CLASS(t_jsb_var_real3d), INTENT(in)    :: rhs

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_var_real3d_to_var_real3d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (SIZE(lhs%ptr3d,2) /= SIZE(rhs%ptr3d,2)) CALL finish(routine, 'Dimension mismatch, vgrid size differs.')
    IF (sub%type == ON_CHUNK) THEN
      lhs%ptr3d(ics:ice,:,iblk) = rhs%ptr3d(ics:ice,:,iblk)
    ELSE
      lhs%ptr3d(:,:,:) = rhs%ptr3d(:,:,:)
    END IF

  END SUBROUTINE Assign_var_real3d_to_var_real3d

  !================================================================================================================================
  ! add
  !================================================================================================================================
  !-----------------------------------------------------------------------------------------------------
  ! 2D
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d + scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real2d_and_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real2d_and_scalar'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) + scalar
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) + scalar
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Add_real2d_and_scalar
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = scalar + var_real2d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_real2d(scalar, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: scalar
    CLASS(t_jsb_var_real2d), INTENT(in) :: var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_real2d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr2d, source=var2%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = scalar + var2%ptr2d(ics:ice,iblk)
    ELSE
      var3%ptr2d(:,:)          = scalar + var2%ptr2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Add_scalar_and_real2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,iblk) = var_real2d(:,iblk) + field1d
  !   ON_CHUNK   field1d(:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Add_real2d_and_field1d(var1, field1d) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field1d(:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real2d_and_field1d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) + field1d(:)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 1d and 2d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Add_real2d_and_field1d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,iblk) = field1d + var_real2d(:,iblk)
  !   ON_CHUNK   field1d(:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Add_field1d_and_real2d(field1d, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: field1d(:)
    CLASS(t_jsb_var_real2d), INTENT(in) :: var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field1d_and_real2d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr2d, source=var2%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = field1d(:) + var2%ptr2d(ics:ice,iblk)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 1d and 2d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Add_field1d_and_real2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d + field2d
  !   ON_CHUNK   field2d(:,1)
  !   ON_DOMAIN  field2d(:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Add_real2d_and_field2d(var1, field2d) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field2d(:,:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real2d_and_field2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field with nblks > 1, op-ov jsb_var real2d on chunk')
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) + field2d(:,1)
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) + field2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Add_real2d_and_field2d
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d + var_real2d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real2d_and_real2d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real2d_and_real2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) + var2%ptr2d(ics:ice,iblk)
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) + var2%ptr2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Add_real2d_and_real2d

  !-----------------------------------------------------------------------------------------------------
  ! 3D
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d + scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real3d_and_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real3d_and_scalar'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) + scalar
    ELSE
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) + scalar
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Add_real3d_and_scalar
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = scalar + var_real3d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_real3d(scalar, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: scalar
    CLASS(t_jsb_var_real3d), INTENT(in) :: var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_real3d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr3d, source=var2%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) =  scalar + var2%ptr3d(ics:ice,:,iblk)
    ELSE
      var3%ptr3d(:,:,:)          =  scalar + var2%ptr3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Add_scalar_and_real3d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,iblk) = var_real3d(:,:,iblk) + field2d | ON_CHUNK
  !   ON_CHUNK   field2d(nc,vgrid)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Add_real3d_and_field2d(var1, field2d) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field2d(:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real3d_and_field2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field2d,2) /= SIZE(var1%ptr,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch of nc or vgrid between field2d and var')
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) + field2d(:,:)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 2d and 3d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Add_real3d_and_field2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,iblk) = field2d + var_real3d(:,:,iblk) | ON_CHUNK
  !   ON_CHUNK   field2d(nc,vgrid)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Add_field2d_and_real3d(field2d, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: field2d(:,:)
    CLASS(t_jsb_var_real3d), INTENT(in) :: var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field2d_and_real3d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr3d, source=var2%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,1) /= SIZE(var2%ptr,1) .OR. SIZE(field2d,2) /= SIZE(var2%ptr,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch of nc or vgrid between field2d and var')
      var3%ptr3d(ics:ice,:,iblk) = field2d(:,:) + var2%ptr3d(ics:ice,:,iblk)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 2d and 3d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Add_field2d_and_real3d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d + field3d | ON_CHUNK & ON_DOMAIN
  !   ON_CHUNK   field3d(nc,:,1)
  !   ON_DOMAIN  field3d(: ,:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Add_real3d_and_field3d(var1, field3d) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field3d(:,:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real3d_and_field3d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var1%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var1%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) + field3d(:,:,1)
    ELSE
      IF (SIZE(field3d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var1%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var1%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) + field3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Add_real3d_and_field3d
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d + var_real3d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real3d_and_real3d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real3d_and_real3d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) + var2%ptr3d(ics:ice,:,iblk)
    ELSE
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) + var2%ptr3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Add_real3d_and_real3d

  !================================================================================================================================
  ! subtract
  !================================================================================================================================
  !-----------------------------------------------------------------------------------------------------
  ! 2D
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d - scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real2d_minus_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_real2d_minus_scalar'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) - scalar
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) - scalar
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Subtract_real2d_minus_scalar
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,iblk) = var_real2d(:,iblk) - field1d
  !   ON_CHUNK   field1d(:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real2d_minus_field1d(var1, field1d) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field1d(:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_real2d_minus_field1d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) - field1d(:)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 1d and 2d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Subtract_real2d_minus_field1d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d - field2d
  !   ON_CHUNK   field2d(:,1)
  !   ON_DOMAIN  field2d(:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real2d_minus_field2d(var1, field2d) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field2d(:,:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_real2d_minus_field2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field with nblks > 1, op-ov jsb_var real2d on chunk')
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) - field2d(:,1)
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) - field2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Subtract_real2d_minus_field2d
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d - var_real2d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real2d_minus_real2d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_real2d_minus_real2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) - var2%ptr2d(ics:ice,iblk)
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) - var2%ptr2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Subtract_real2d_minus_real2d

  !-----------------------------------------------------------------------------------------------------
  ! 3D
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d - scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real3d_minus_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_real3d_minus_scalar'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) - scalar
    ELSE
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) - scalar
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Subtract_real3d_minus_scalar
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,iblk) = var_real3d(:,:,iblk) - field2d | ON_CHUNK
  !   ON_CHUNK   field2d(nc,vgrid)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real3d_minus_field2d(var1, field2d) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field2d(:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_real3d_minus_field2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field2d,2) /= SIZE(var1%ptr,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch of nc or vgrid between field2d and var')
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) - field2d(:,:)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 2d and 3d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Subtract_real3d_minus_field2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d - field3d | ON_CHUNK & ON_DOMAIN
  !   ON_CHUNK   field3d(nc,:,1)
  !   ON_DOMAIN  field3d(: ,:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real3d_minus_field3d(var1, field3d) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field3d(:,:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_real3d_minus_field3d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var1%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var1%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) - field3d(:,:,1)
    ELSE
      IF (SIZE(field3d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var1%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var1%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) - field3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Subtract_real3d_minus_field3d
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d - var_real3d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real3d_minus_real3d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_real3d_minus_real3d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) - var2%ptr3d(ics:ice,:,iblk)
    ELSE
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) - var2%ptr3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Subtract_real3d_minus_real3d

  !================================================================================================================================
  ! multiply
  !================================================================================================================================
  !-----------------------------------------------------------------------------------------------------
  ! 2D
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d * scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real2d_with_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real2d_with_scalar'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) * scalar
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) * scalar
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Multiply_real2d_with_scalar
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = scalar * var_real2d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_with_real2d(scalar, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: scalar
    CLASS(t_jsb_var_real2d), INTENT(in) :: var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_real2d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr2d, source=var2%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) =  scalar * var2%ptr2d(ics:ice,iblk)
    ELSE
      var3%ptr2d(:,:)          =  scalar * var2%ptr2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Multiply_scalar_with_real2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,iblk) = var_real2d(:,iblk) * field1d
  !   ON_CHUNK   field1d(:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real2d_with_field1d(var1, field1d) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field1d(:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real2d_with_field1d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) * field1d(:)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 1d and 2d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Multiply_real2d_with_field1d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,iblk) = field1d * var_real2d(:,iblk)
  !   ON_CHUNK   field1d(:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field1d_with_real2d(field1d, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: field1d(:)
    CLASS(t_jsb_var_real2d), INTENT(in) :: var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field1d_with_real2d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr2d, source=var2%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = field1d(:) * var2%ptr2d(ics:ice,iblk)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 1d and 2d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Multiply_field1d_with_real2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d * field2d
  !   ON_CHUNK   field2d(:,1)
  !   ON_DOMAIN  field2d(:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real2d_with_field2d(var1, field2d) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field2d(:,:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real2d_with_field2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field with nblks > 1, op-ov jsb_var real2d on chunk')
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) * field2d(:,1)
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) * field2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Multiply_real2d_with_field2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = field2d * var_real2d
  !   ON_CHUNK   field2d(:,1)
  !   ON_DOMAIN  field2d(:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field2d_with_real2d(field2d, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: field2d(:,:)
    CLASS(t_jsb_var_real2d), INTENT(in) :: var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_real2d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr2d, source=var2%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field with nblks > 1, op-ov jsb_var real2d on chunk')
      var3%ptr2d(ics:ice,iblk) = field2d(:,1) * var2%ptr2d(ics:ice,iblk)
    ELSE
      var3%ptr2d(:,:)          = field2d(:,:) * var2%ptr2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Multiply_field2d_with_real2d
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d * var_real2d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real2d_with_real2d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real2d_with_real2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) * var2%ptr2d(ics:ice,iblk)
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) * var2%ptr2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Multiply_real2d_with_real2d

  !-----------------------------------------------------------------------------------------------------
  ! 3D
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d * scalar | ON_CHUNK & ON_DOMAIN
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real3d_with_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real3d_with_scalar'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) * scalar
    ELSE
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) * scalar
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Multiply_real3d_with_scalar
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = scalar * var_real3d | ON_CHUNK & ON_DOMAIN
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_with_real3d(scalar, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: scalar
    CLASS(t_jsb_var_real3d), INTENT(in) :: var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_real3d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr3d, source=var2%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = scalar * var2%ptr3d(ics:ice,:,iblk)
    ELSE
      var3%ptr3d(:,:,:)          = scalar * var2%ptr3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Multiply_scalar_with_real3d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,iblk) = var_real3d(:,:,iblk) * field2d | ON_CHUNK
  !   ON_CHUNK   field2d(nc,vgrid)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real3d_with_field2d(var1, field2d) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field2d(:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real3d_with_field2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field2d,2) /= SIZE(var1%ptr,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch of nc or vgrid between field2d and var')
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) * field2d(:,:)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 2d and 3d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Multiply_real3d_with_field2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,iblk) = field2d * var_real3d(:,:,iblk) | ON_CHUNK
  !   ON_CHUNK   field2d(nc,vgrid)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field2d_with_real3d(field2d, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: field2d(:,:)
    CLASS(t_jsb_var_real3d), INTENT(in) :: var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_real3d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr3d, source=var2%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,1) /= SIZE(var2%ptr,1) .OR. SIZE(field2d,2) /= SIZE(var2%ptr,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch of nc or vgrid between field2d and var')
      var3%ptr3d(ics:ice,:,iblk) = field2d(:,:) * var2%ptr3d(ics:ice,:,iblk)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 2d and 3d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Multiply_field2d_with_real3d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d * field3d | ON_CHUNK & ON_DOMAIN
  !   ON_CHUNK   field3d(nc,:,1)
  !   ON_DOMAIN  field3d(: ,:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real3d_with_field3d(var1, field3d) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field3d(:,:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real3d_with_field3d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var1%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var1%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) * field3d(:,:,1)
    ELSE
      IF (SIZE(field3d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var1%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var1%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) * field3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Multiply_real3d_with_field3d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = field3d * var_real3d | ON_CHUNK & ON_DOMAIN
  !   ON_CHUNK   field3d(nc,:,1)
  !   ON_DOMAIN  field3d(: ,:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field3d_with_real3d(field3d, var2) RESULT(var3)

    REAL(wp),                INTENT(in) :: field3d(:,:,:)
    CLASS(t_jsb_var_real3d), INTENT(in) :: var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field3d_with_real3d'

    sub = var2%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var2%name
    var3%full_name        = var2%full_name
    var3%owner_proc_id    = var2%owner_proc_id
    var3%owner_model_id   = var2%owner_model_id
    var3%owner_tile_path  = var2%owner_tile_path
    var3%missval          = var2%missval
    ALLOCATE(var3%ptr3d, source=var2%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,1) /= SIZE(var2%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var2%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var2%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(ics:ice,:,iblk) = field3d(:,:,1) * var2%ptr3d(ics:ice,:,iblk)
    ELSE
      IF (SIZE(field3d,1) /= SIZE(var2%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var2%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var2%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(:,:,:)          = field3d(:,:,:) * var2%ptr3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Multiply_field3d_with_real3d
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d * var_real3d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real3d_with_real3d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real3d_with_real3d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) * var2%ptr3d(ics:ice,:,iblk)
    ELSE
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) * var2%ptr3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Multiply_real3d_with_real3d

  !================================================================================================================================
  ! divide
  !================================================================================================================================
  !-----------------------------------------------------------------------------------------------------
  ! 2D
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d / scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Div_real2d_by_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_real2d_by_scalar'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) / scalar
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) / scalar
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Div_real2d_by_scalar
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,iblk) = var_real2d(:,iblk) / field1d
  !   ON_CHUNK   field1d(:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Div_real2d_by_field1d(var1, field1d) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field1d(:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_real2d_by_field1d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) / field1d(:)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 1d and 2d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Div_real2d_by_field1d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d / field2d
  !   ON_CHUNK   field2d(:,1)
  !   ON_DOMAIN  field2d(:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Div_real2d_by_field2d(var1, field2d) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field2d(:,:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_real2d_by_field2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field with nblks > 1, op-ov jsb_var real2d on chunk')
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) / field2d(:,1)
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) / field2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Div_real2d_by_field2d
  !-----------------------------------------------------------------------------------------------------
  ! var_real2d(:,:) = var_real2d / var_real2d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Div_real2d_by_real2d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_real2d_by_real2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr2d, source=var1%ptr2d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr2d(ics:ice,iblk) = var1%ptr2d(ics:ice,iblk) / var2%ptr2d(ics:ice,iblk)
    ELSE
      var3%ptr2d(:,:)          = var1%ptr2d(:,:) / var2%ptr2d(:,:)
    END IF

    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Div_real2d_by_real2d

  !-----------------------------------------------------------------------------------------------------
  ! 3D
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d / scalar | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Div_real3d_by_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_real3d_by_scalar'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) / scalar
    ELSE
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) / scalar
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Div_real3d_by_scalar
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,iblk) = var_real3d(:,:,iblk) / field2d | ON_CHUNK
  !   ON_CHUNK   field2d(nc,vgrid)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Div_real3d_by_field2d(var1, field2d) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field2d(:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_real3d_by_field2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field2d,2) /= SIZE(var1%ptr,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch of nc or vgrid between field2d and var')
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) / field2d(:,:)
    ELSE
      CALL finish(routine, 'Ambiguous dimension-matching 2d and 3d ON_DOMAIN')
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Div_real3d_by_field2d
  ! -----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d / field3d | ON_CHUNK & ON_DOMAIN
  !   ON_CHUNK   field3d(nc,:,1)
  !   ON_DOMAIN  field3d(: ,:,:)
  ! -----------------------------------------------------------------------------------------------------
  FUNCTION Div_real3d_by_field3d(var1, field3d) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field3d(:,:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_real3d_by_field3d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var1%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var1%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) / field3d(:,:,1)
    ELSE
      IF (SIZE(field3d,1) /= SIZE(var1%ptr,1) .OR. SIZE(field3d,2) /= SIZE(var1%ptr,2) .OR. SIZE(field3d,3) /= SIZE(var1%ptr,3)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) / field3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Div_real3d_by_field3d
  !-----------------------------------------------------------------------------------------------------
  ! var_real3d(:,:,:) = var_real3d / var_real3d | ON_CHUNK & ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Div_real3d_by_real3d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_real3d_by_real3d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    var3%name             = var1%name
    var3%full_name        = var1%full_name
    var3%owner_proc_id    = var1%owner_proc_id
    var3%owner_model_id   = var1%owner_model_id
    var3%owner_tile_path  = var1%owner_tile_path
    var3%missval          = var1%missval
    ALLOCATE(var3%ptr3d, source=var1%ptr3d)

    IF (sub%type == ON_CHUNK) THEN
      var3%ptr3d(ics:ice,:,iblk) = var1%ptr3d(ics:ice,:,iblk) / var2%ptr3d(ics:ice,:,iblk)
    ELSE
      var3%ptr3d(:,:,:)          = var1%ptr3d(:,:,:) / var2%ptr3d(:,:,:)
    END IF

    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Div_real3d_by_real3d

#endif
END MODULE mo_jsb_var_class
