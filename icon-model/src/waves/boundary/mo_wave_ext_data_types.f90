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

! Data type defintion for wave external data state
!
! Defines the data type for storing wave-specific external parameter
! fields such as bathymetry.

MODULE mo_wave_ext_data_types

  USE mo_kind,               ONLY: wp
  USE mo_fortran_tools,      ONLY: t_ptr_2d3d
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: t_external_wave


  ! wave-specific external data type
  !
  TYPE :: t_external_wave
    ! ocean topography <=> bathymetric height used in the ocean
    ! cell centers and edges only
    !
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  bathymetry_c(:,:),          &  !< topographic height at cell centers [m]
                                        !  index1=1,nproma, index2=1,nblks_c
      &  bathymetry_e(:,:)              !< topographic height at cell edges    [m]
                                        !  index1=1,nproma, index2=1,nblks_e

    REAL(wp), POINTER, CONTIGUOUS :: &
      &  geo_depth_grad_c(:,:,:)       !< bathymetry geographical gradient [m/m]
                                       !  index1=2, index2=1,nproma, index3=1,nblks_c

    REAL(wp), POINTER, CONTIGUOUS :: & !<  water depth at cell center [m]
      &  depth_c(:,:)                    !   index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER, CONTIGUOUS :: & !<  water depth at cell edges [m]
      &  depth_e(:,:)                    !   index1=1,nproma, index2=1,nblks_e

    TYPE(t_ptr_2d3d), ALLOCATABLE :: grad_ptr(:)

  END TYPE t_external_wave

END MODULE mo_wave_ext_data_types
