!> QUINCY phenology constants
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### declare and define phenology constants
!>
MODULE mo_q_pheno_constants
#ifndef __NO_QUINCY__

  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PUBLIC

  !> photosynthetic pathways and growth forms in the QUINCY model
  !>
  !> TODO re-write as ENUM
  !>
  INTEGER, PARAMETER ::                     &
      ievergreen              = 1         , &
      isummergreen            = 2         , &
      iraingreen              = 3         , &
      iperennial              = 4

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_pheno_constants'

#endif
END MODULE mo_q_pheno_constants
