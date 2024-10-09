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

! Module containing types for the reading / on-the-fly generation and the writing of the vertical grid.

MODULE mo_util_vgrid_types

  USE mo_kind, ONLY: wp
  USE mo_util_uuid_types, ONLY: t_uuid
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_vgrid_buffer
  PUBLIC :: vgrid_buffer

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_vgrid_types'


  TYPE t_vgrid_buffer
    ! Geom. height at vertical interface of cells and vertices (nproma,nlevp1,nblks_c/nblks_v): 
    !
    ! Note: This array is deallocated after being used in
    !       mo_vertical_grid::set_nh_metrics
    !
    REAL(wp), POINTER, CONTIGUOUS :: z_ifc(:,:,:)

    ! UUID of vertical grid
    TYPE(t_uuid) :: uuid
  END TYPE t_vgrid_buffer

  ! module variable: temporary buffer for coordinate arrays
  TYPE (t_vgrid_buffer), ALLOCATABLE :: vgrid_buffer(:)


END MODULE mo_util_vgrid_types
