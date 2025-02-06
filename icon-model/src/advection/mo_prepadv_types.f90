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

! Type definition for preparation of transport with optional reduced
! calling frequency.

MODULE mo_prepadv_types

  USE mo_kind,                 ONLY: wp
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: t_prepare_adv
  PUBLIC :: t_step_adv


  ! for preparation of transport with optional reduced calling frequency
  !
  TYPE :: t_prepare_adv

    REAL(wp), POINTER, CONTIGUOUS ::                    &
      mass_flx_me(:,:,:) ,  & !< mass flux at full level edges   [kg/m^2/s]
                              !< averaged over dynamics substeps
      mass_flx_ic(:,:,:) ,  & !< mass flux at half level centers [kg/m^2/s]
                              !< averaged over dynamics substeps
      vol_flx_ic(:,:,:) ,   & !< volume flux at half level centers [m/s]
                              !< averaged over dynamics substeps
      vn_traj    (:,:,:) ,  & !< horizontal velocity at edges for computation of backward trajectories [m/s]
                              !< averaged over dynamics substeps
      q_int      (:,:,:) ,  & !< Storage field for vertical nesting: 
                              !< q at child nest interface level [kg/kg]
      q_ubc      (:,:,:)      !< Storage field for vertical nesting:
                              !< q at (nest) upper boundary      [kg/kg]

    TYPE(t_ptr_2d3d),ALLOCATABLE ::   &
      q_int_ptr(:),         & !< pointer array: one pointer for each tracer
      q_ubc_ptr(:)            !< pointer array: one pointer for each tracer

  END TYPE t_prepare_adv


  ! counter for Marchuk splitting
  TYPE :: t_step_adv
    ! Determines sequence of operations for Marchuk-splitting (for transport)
    INTEGER :: marchuk_order
  END TYPE t_step_adv

END MODULE mo_prepadv_types

