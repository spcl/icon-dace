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

! Arrays used to store ice variables and organize coupling

MODULE mo_ice_fem_types

  USE mo_kind,    ONLY: wp

  IMPLICIT NONE

  PUBLIC
  SAVE

  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: u_ice, v_ice, m_ice, a_ice  
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: rhs_u, rhs_v, m_snow
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: rhs_m, rhs_a, rhs_mis
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: u_w, v_w
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: elevation
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: mass_matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: lmass_matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: sigma11, sigma12, sigma22  
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_x         
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_y
  ! auxiliary arrays required for fct
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: m_icel, a_icel, m_snowl
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: dm_ice, da_ice, dm_snow
  REAL(wp), ALLOCATABLE, DIMENSION(:,:)       :: icefluxes
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: icepplus, icepminus

! There are, in principle, 3 matrices which share the same structure stored in icestiff.
! We use icestiff%values for both implicit advection schemes and VP algorithms.
! Since they are re-assembled at each time step, there is no problem with sharing it.
! Additionally, there is mass matrix with the same structure, which is used in FCT advection scheme.
! So we use mass_matrix of the size of icestiff%values only to keep  the values.

  TYPE sparse_matrix
     INTEGER :: nza
     INTEGER :: dim
     REAL(wp),          POINTER, DIMENSION(:)   :: values
     INTEGER(KIND=4),   POINTER, DIMENSION(:)   :: colind
     INTEGER(KIND=4),   POINTER, DIMENSION(:)   :: rowptr
  END TYPE sparse_matrix

  TYPE(sparse_matrix) :: icestiff

END MODULE mo_ice_fem_types
