# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/parallel_infrastructure/mo_communication_yaxt.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/parallel_infrastructure/mo_communication_yaxt.f90"
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

! This module provides the yaxt based communication routines
! for parallel runs

!----------------------------

# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/include/icon_definitions.inc" 1
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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
!

















# 17 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/parallel_infrastructure/mo_communication_yaxt.f90" 2

# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/include/crayftn_ptr_fail.inc" 1
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

! Cray ftn compilers 8.4 and 8.6 are known to misidentify argument INTENT
! of pointer components, i.e. will disallow changes to an array pointed to
! by a pointer component of a TYPE
! therefore this case needs to be handled specially
!
# 18 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/parallel_infrastructure/mo_communication_yaxt.f90" 2
!----------------------------
MODULE mo_communication_yaxt
# 2985 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/parallel_infrastructure/mo_communication_yaxt.f90"
! HAVE_YAXT
END MODULE mo_communication_yaxt
