!
! mo_art_OH_chem_types
! This module provides types for the simplified OH chemistry
! (lart_chemtracer == .TRUE. and <c_solve>OH</c_solve> for one of the tracers)
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

MODULE mo_art_OH_chem_types
! ICON
  USE mo_kind,                          ONLY: wp
! ART

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_art_OH_chem_meta
  
  TYPE t_art_OH_chem_meta
    INTEGER :: num_elements        !< number of elements in the list
    REAL(wp), ALLOCATABLE ::  &
     &  ozone_Nconc(:,:,:)         !< ozone number concentration (cm-3) (either GEMS climatology 
                                   !  or linearized ozone chemistry)
    LOGICAL ::   &
     &  is_init = .FALSE.          !< flag if structure is initialized
    REAL(wp), POINTER ::          &
     &  OH_Nconc(:,:,:) => NULL(),& !< diagnostic number concentration of OH in cm-3
     &  CH4_star(:,:,:) => NULL(),& !< pointers to the star values of CH4 and CO
     &  CO_star(:,:,:) => NULL()

    INTEGER, ALLOCATABLE :: &
      &  level_CH4_gt_1ppm(:,:)   !< level index for which CH4 is greater than 1 ppmv
    REAL(wp), ALLOCATABLE :: &
      &  k_CH4_OH(:,:,:),    &    !< reaction rate constants of CH4 and CO with OH (cm3/s)
      &  k_CO_OH(:,:,:)
  END TYPE t_art_OH_chem_meta

END MODULE mo_art_OH_chem_types
