!> Contains basic structure for sub-sets
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
MODULE mo_jsb_subset

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_subsets
  PUBLIC :: t_subset, ON_DOMAIN, ON_CHUNK

  ENUM, BIND(C)
    ENUMERATOR :: ON_DOMAIN=1, ON_CHUNK
  END ENUM

  TYPE t_subset
    INTEGER :: type = 0   ! is actually never set to ON_DOMAIN but only set to ON_CHUNK with interface_full() -> model%Set_subset() after PROC_init were called
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: nb
    INTEGER :: nc
  END TYPE t_subset

  TYPE t_model_subsets
    TYPE(t_subset), ALLOCATABLE :: sub (:)
  END TYPE t_model_subsets

  TYPE(t_model_subsets), ALLOCATABLE, TARGET, SAVE :: jsbach_subsets(:)

END MODULE mo_jsb_subset
