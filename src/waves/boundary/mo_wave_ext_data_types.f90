! Data type defintion for wave external data state
!
! Defines the data type for storing wave-specific external parameter
! fields such as bathymetry.
!
!
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

MODULE mo_wave_ext_data_types

  USE mo_kind,               ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: t_external_wave


  ! wave-specific external data type
  !
  TYPE :: t_external_wave

    ! ocean topography <=> bathymetric height used in the ocean
    ! cell centers and edges only
    !
    REAL(wp), POINTER ::   &  !<  bathymetric height at cell centers [m]
      &  bathymetry_c(:,:)    !  index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &  !< topographic height at cell edges    [m]
      &  bathymetry_e(:,:)    !  index1=1,nproma, index2=1,nblks_e

    REAL(wp), POINTER ::   &  !< topographic height at cell vertices [m]
      &  bathymetry_v(:,:)    !  index1=1,nproma, index2=1,nblks_v

    ! *** Land-Sea-Mask ***
    INTEGER, POINTER  ::   &  !< land-sea-mask for cell centers          [ ]
      &  lsm_wave_ctr_c(:,:)       !  index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::    &  !< land-sea-mask for cell edges
      &  lsm_wave_ctr_e(:,:)       !  index1=1,nproma, index2=1,nblks_e

  END TYPE t_external_wave

END MODULE mo_wave_ext_data_types
