!> constants for the tcq process (test conserved quantities)
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
MODULE mo_tcq_constants
#ifndef __NO_JSBACH__

  USE mo_jsb_varlist,            ONLY: VARNAME_LEN

  IMPLICIT NONE
  PUBLIC

  ! Definitions for lcc

  !> list of potential active variables for the tcq process: For these variables the tcq interface provides active relocation
  CHARACTER(len=VARNAME_LEN) :: tcq_potential_active_vars(1) = [character(len=VARNAME_LEN) :: 'a_veg_c_tcq' ]

  !> list of required passive vars if vars are active (because these are used to relocate the active variables)
  CHARACTER(len=VARNAME_LEN) :: tcq_required_passive_vars(1) = [character(len=VARNAME_LEN) :: 'a_dead_c_tcq' ]

#endif
END MODULE mo_tcq_constants
