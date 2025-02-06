!> fage (forest age) utilities
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
!>#### Contains utilities for the fage (forest age) proc
!>
MODULE mo_fage_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_mean_age

  INTERFACE get_mean_age
    MODULE PROCEDURE get_mean_age_1d
    MODULE PROCEDURE get_mean_age_2d
    MODULE PROCEDURE get_mean_age_3d
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_util'

CONTAINS


  ! ====================================================================================================== !
  !
  !> Calculate the mean age for one cell
  !
  FUNCTION get_mean_age_1d(fract_per_age) RESULT(mean_age)
    !-----------------------------------------------------------------------
    REAL(wp), INTENT(in) :: &
      & fract_per_age (:)     !< fraction for each age (nidx,z)
    REAL(wp) :: mean_age      !< mean age
    !-----------------------------------------------------------------------
    INTEGER  :: i, n
    REAL(wp) :: fract_sum
    CHARACTER(len=*), PARAMETER :: routine = modname//':get_mean_age_1d'
    !-----------------------------------------------------------------------

    fract_sum = SUM(fract_per_age)
    n = SIZE(fract_per_age)

    IF (fract_sum >= 0._wp) THEN
      mean_age = 0._wp
      DO i = 1,n
          mean_age  = mean_age + i * (fract_per_age (i) /  fract_sum)
      END DO
    END IF
  END FUNCTION get_mean_age_1d

  ! ====================================================================================================== !
  !
  !> Calculate the mean age for a block
  !
  FUNCTION get_mean_age_2d(fract_per_age) RESULT(mean_age)
    !-----------------------------------------------------------------------
    REAL(wp), INTENT(in) :: &
      & fract_per_age (:,:)     !< fraction for each age (nidx,z)
    REAL(wp) :: mean_age(SIZE(fract_per_age,1))  !< mean age
    !-----------------------------------------------------------------------
    INTEGER  :: i, n
    REAL(wp) :: fract_sum(SIZE(fract_per_age,1))
    CHARACTER(len=*), PARAMETER :: routine = modname//':get_mean_age_2d'
    !-----------------------------------------------------------------------

    n  = SIZE(fract_per_age,2)       ! max tracked age

    fract_sum(:) = 0._wp
    DO i=1,n
      fract_sum(:) = fract_sum(:) + fract_per_age (:,i)
    END DO

    mean_age(:)  = 0._wp
    DO i=1,n
      WHERE(fract_sum(:) > 0._wp)
        mean_age(:)  = mean_age (:) + i * (fract_per_age (:,i) /  fract_sum(:))
      END WHERE
    END DO

  END FUNCTION get_mean_age_2d

  ! ====================================================================================================== !
  !
  !> Calculate the mean age on domain
  !
  FUNCTION get_mean_age_3d(fract_per_age) RESULT(mean_age)
    !-----------------------------------------------------------------------
    REAL(wp), INTENT(in) :: &
      & fract_per_age (:,:,:)     !< fraction for each age (nidx,z,nbkls)
    REAL(wp) :: mean_age(SIZE(fract_per_age,1),SIZE(fract_per_age,3))  !< mean age
    !-----------------------------------------------------------------------
    INTEGER :: i,n
    REAL(wp) :: fract_sum(SIZE(fract_per_age,1),SIZE(fract_per_age,3))
    CHARACTER(len=*), PARAMETER :: routine = modname//':get_mean_age_2d'
    !-----------------------------------------------------------------------

    n  = SIZE(fract_per_age,2)       ! max tracked age

    fract_sum(:,:) = 0._wp
    DO i=1,n
      fract_sum(:,:) = fract_sum(:,:) + fract_per_age (:,i,:)
    END DO

    mean_age (:,:) = 0._wp
    DO i=1,n
      WHERE(fract_sum(:,:) > 0._wp)
        mean_age (:,:) = mean_age (:,:) + i * (fract_per_age (:,i,:) / fract_sum(:,:))
      END WHERE
    END DO

  END FUNCTION get_mean_age_3d


#endif
END MODULE mo_fage_util
