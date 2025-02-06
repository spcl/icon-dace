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

!>
!!  Contains the definition of basic structures and geometry parameters
!!  These are included in the grid/patch info
!!

MODULE mo_lib_grid_geometry_info

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64
  USE mo_math_types, ONLY: t_cartesian_coordinates

  IMPLICIT NONE

  PRIVATE

  ! public parameters
  PUBLIC :: sphere_geometry, planar_torus_geometry, planar_channel_geometry, planar_geometry
  PUBLIC :: triangular_cell, hexagonal_cell
  PUBLIC :: cut_off_grid, refined_bisection_grid, dualy_refined_grid
  PUBLIC :: undefined

  ! public structures
  PUBLIC :: t_grid_geometry_info

  ! -----------------------------
  ! types of grid geometries
  INTEGER, PARAMETER ::  sphere_geometry = 1
  INTEGER, PARAMETER ::  planar_torus_geometry = 2
  INTEGER, PARAMETER ::  planar_channel_geometry = 3
  INTEGER, PARAMETER ::  planar_geometry = 4

  ! -----------------------------
  ! types of grids
  INTEGER, PARAMETER ::  triangular_cell = 3
  INTEGER, PARAMETER ::  hexagonal_cell = 6

  ! -----------------------------
  ! types of grid creation
  INTEGER, PARAMETER ::  cut_off_grid = 1
  INTEGER, PARAMETER ::  refined_bisection_grid = 2
  INTEGER, PARAMETER ::  dualy_refined_grid = 3

  !--------------------------------------------------------------
  INTEGER, PARAMETER ::  undefined = -1
  !--------------------------------------------------------------
  !> Holds the grid geometry parameters
  TYPE t_grid_geometry_info
    INTEGER :: cell_type ! triangular_cell, etc
    INTEGER :: geometry_type ! sphere_geometry, etc
    !> The creation process of the grid (cut_off, refined, etc, see parameters).
    INTEGER :: grid_creation_process
    !> The grid optimization process
    INTEGER :: grid_optimization_process

    TYPE(t_cartesian_coordinates) :: center
    REAL(wp) :: mean_edge_length ! (meters)
    REAL(wp) :: mean_dual_edge_length ! (meters)
    REAL(wp) :: mean_cell_area ! (meters^2)
    REAL(wp) :: mean_dual_cell_area ! (meters^2)
    REAL(wp) :: domain_length ! (meters)
    REAL(wp) :: domain_height ! (meters)

    !> Sphere/Cylinder parameters
    REAL(wp) :: sphere_radius

    !> derived info, the following can be calculated from the previous
    REAL(wp) :: mean_characteristic_length ! the sqrt(mean_cell_area)

  END TYPE t_grid_geometry_info

END MODULE mo_lib_grid_geometry_info
!----------------------------------------------------------------------------

