# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/shared/mo_openacc.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/shared/mo_openacc.f90"
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

! This module provides fortran interface for OpenACC routines,
! because the OpenACC standard only defines C-Interfaces.

! (GZ, 2013-08-30): So far, the Cray compiler is the only one for which an OpenMP parallelization
! of copying data into / back from the MPI-buffer seems to give a benefit. Further compilers may
! be added here once the OpenMP implementation is sufficiently efficient

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

















# 21 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/shared/mo_openacc.f90" 2
!----------------------------
MODULE mo_openacc
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
# 404 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/shared/mo_openacc.f90"

END MODULE mo_openacc
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
