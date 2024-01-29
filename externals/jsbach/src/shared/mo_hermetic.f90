!> Memory management module for automatic destruction of temporary objects to avoid memory leaks
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
!> After Rouson, Xia and Xu (2011). Scientific Software Design - The Object-Oriented Way
!> See also Stewart, G. W. (2007). Memory leaks in derived types revisited.
!>          ACM SIGPLAN Fortran Forum, 22(3), 25–27. https://doi.org/10.1145/962180.962183
!>
MODULE mo_hermetic

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_hermetic ! Expose type and type-bound procedures

  TYPE, ABSTRACT :: t_hermetic
    PRIVATE
    INTEGER, POINTER :: temporary => NULL() ! Null marks non-temporary data
  CONTAINS
    PROCEDURE                            :: SetTemp   ! Mark object as temporary
    PROCEDURE                            :: GuardTemp ! Increment the depth count
    PROCEDURE                            :: CleanTemp ! Decrement depth count / free memory if 1
    PROCEDURE(final_interface), DEFERRED :: force_finalization
  END TYPE t_hermetic

  ABSTRACT INTERFACE
    SUBROUTINE final_interface(this)
      IMPORT :: t_hermetic
      CLASS(t_hermetic), INTENT(inout) :: this
    END SUBROUTINE
  END INTERFACE

CONTAINS

  SUBROUTINE SetTemp(this)

    CLASS(t_hermetic), INTENT(inout) :: this

    IF (.NOT. ASSOCIATED(this%temporary)) ALLOCATE(this%temporary)
    this%temporary = 1

  END SUBROUTINE SetTemp

  SUBROUTINE GuardTemp(this)

    CLASS(t_hermetic), INTENT(inout) :: this

    IF (ASSOCIATED(this%temporary)) this%temporary = this%temporary + 1

  END SUBROUTINE GuardTemp

  SUBROUTINE CleanTemp(this)

    CLASS(t_hermetic), INTENT(inout) :: this

    IF (ASSOCIATED(this%temporary)) THEN
      IF (this%temporary > 1) this%temporary = this%temporary - 1
      IF (this%temporary == 1) THEN
        CALL this%force_finalization()
        DEALLOCATE(this%temporary)
      END IF
    END IF

  END SUBROUTINE CleanTemp

END MODULE mo_hermetic
