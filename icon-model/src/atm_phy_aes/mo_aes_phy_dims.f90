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

! Dimensions for the AES physics package.

MODULE mo_aes_phy_dims

  USE mo_impl_constants,  ONLY: max_dom
  USE mo_parallel_config, ONLY: nproma
  USE mo_run_config,      ONLY: num_lev, ntracer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: aes_phy_dims, init_aes_phy_dims

  !>
  !! Dimensions for physics, for multiple domains/grids.
  !!
  TYPE t_aes_phy_dims
     !
     INTEGER :: nproma  !< size of cells dimension
     INTEGER :: nlev    !< size of levels dimension
     INTEGER :: ntracer !< size of tracers dimension
     !
  END type t_aes_phy_dims

  TYPE(t_aes_phy_dims)  , TARGET :: aes_phy_dims   (max_dom)

CONTAINS

  !----

  !>
  !! Initialize the dimensions
  !!
  SUBROUTINE init_aes_phy_dims
    !
    aes_phy_dims(:)%nproma  = nproma
    aes_phy_dims(:)%ntracer = ntracer
    aes_phy_dims(:)%nlev    = num_lev(:)
    !
  END SUBROUTINE init_aes_phy_dims

  !----

END MODULE mo_aes_phy_dims
