! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
! Get the definition of the 'YAC_MPI_FINT_FC_KIND':
#include "config.h"
#endif

module yac_utils
  use, intrinsic :: iso_c_binding, only : &
    c_int, c_long, c_long_long, c_short, c_char

  implicit none

  public

  !------------------------------------------------
  ! Constants for Fortran-C interoperability
  !------------------------------------------------

  integer, parameter :: YAC_MPI_FINT_KIND = YAC_MPI_FINT_FC_KIND

  !---------------
  ! C free routine
  !---------------

  interface

    subroutine yac_free_c ( ptr ) bind ( c, NAME='free' )

       use, intrinsic :: iso_c_binding

       type(c_ptr), value :: ptr

     end subroutine yac_free_c

  end interface

  !-------------------------------------------------------------
  ! Duplicates stencils
  ! * can be used in case duplicated points have been masked out
  !-------------------------------------------------------------

  interface

    subroutine yac_duplicate_stencils_c( &
      weights, tgt_grid, tgt_orig_global_id, tgt_duplicated_idx, &
      nbr_duplicated, location) &
      bind ( c, name='yac_duplicate_stencils_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: weights
      type(c_ptr), value            :: tgt_grid
      type(c_ptr), value            :: tgt_orig_global_id ! yac_int *
      type(c_ptr), value            :: tgt_duplicated_idx ! size_t *
      integer(kind=c_size_t), value :: nbr_duplicated
      integer(kind=c_int), value    :: location

    end subroutine yac_duplicate_stencils_c
  end interface

  !------------------------
  ! Cubed sphere generation
  !------------------------

  interface

    subroutine yac_generate_cubed_sphere_grid_information_c( &
      n, num_cells, num_vertices, x_vertices, y_vertices, z_vertices, &
      vertices_of_cell, face_id) &
      bind ( c, name='yac_generate_cubed_sphere_grid_information' )

      use, intrinsic :: iso_c_binding

      integer(kind=c_int), value :: n
      integer(kind=c_int)        :: num_cells
      integer(kind=c_int)        :: num_vertices
      type(c_ptr)                :: x_vertices ! double **
      type(c_ptr)                :: y_vertices ! double **
      type(c_ptr)                :: z_vertices ! double **
      type(c_ptr)                :: vertices_of_cell ! unsigned **
      type(c_ptr)                :: face_id ! unsigned **
    end subroutine yac_generate_cubed_sphere_grid_information_c

    subroutine yac_generate_part_cube_grid_information_c( &
      n, nbr_vertices, nbr_cells, num_vertices_per_cell, cell_to_vertex, &
      x_vertices, y_vertices, x_cells, y_cells, global_cell_id, &
      cell_core_mask, global_corner_id, corner_core_mask, rank, size) &
      bind ( c, name='yac_generate_part_cube_grid_information' )

      use, intrinsic :: iso_c_binding

      integer(kind=c_int), value :: n
      integer(kind=c_int)        :: nbr_vertices
      integer(kind=c_int)        :: nbr_cells
      type(c_ptr)                :: num_vertices_per_cell ! unsigned **
      type(c_ptr)                :: cell_to_vertex ! unsigned **
      type(c_ptr)                :: x_vertices ! double **
      type(c_ptr)                :: y_vertices ! double **
      type(c_ptr)                :: x_cells ! double **
      type(c_ptr)                :: y_cells ! double **
      type(c_ptr)                :: global_cell_id ! int **
      type(c_ptr)                :: cell_core_mask ! int **
      type(c_ptr)                :: global_corner_id ! int **
      type(c_ptr)                :: corner_core_mask ! int **
      integer(kind=c_int), value :: rank
      integer(kind=c_int), value :: size
    end subroutine yac_generate_part_cube_grid_information_c

    function yac_generate_cubed_sphere_basic_grid_c(name, n) &
      bind ( c, name='yac_generate_cubed_sphere_basic_grid' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: name(*)
      integer(kind=c_size_t), value :: n

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_generate_cubed_sphere_basic_grid_c
    end function yac_generate_cubed_sphere_basic_grid_c

    subroutine yac_write_cubed_sphere_grid_c(n, filename) &
      bind ( c, name='yac_write_cubed_sphere_grid' )

      use, intrinsic :: iso_c_binding

      integer(kind=c_int), value :: n
      character(kind=c_char) :: filename(*)
    end subroutine yac_write_cubed_sphere_grid_c

  end interface

  !----------------------------
  ! 2D decomposition generation
  !----------------------------

  interface

    subroutine yac_generate_reg2d_decomp_c( &
      num_points, total_num_procs, num_procs) &
      bind ( c, name='yac_generate_reg2d_decomp' )

      use, intrinsic :: iso_c_binding

      integer(kind=c_int)        :: num_points(2)
      integer(kind=c_int), value :: total_num_procs
      integer(kind=c_int)        :: num_procs(2)
    end subroutine yac_generate_reg2d_decomp_c

  end interface

  !-------------------------
  ! Writing grid to VTK file
  !-------------------------

  interface

    subroutine yac_write_basic_grid_data_to_file_c(grid, name) &
      bind ( c, name='yac_write_basic_grid_data_to_file' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value     :: grid ! struct yac_basic_grid_data
      character(kind=c_char) :: name(*)
    end subroutine yac_write_basic_grid_data_to_file_c

    subroutine yac_write_basic_grid_to_file_c(grid, name) &
      bind ( c, name='yac_write_basic_grid_to_file' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value     :: grid ! struct yac_basic_grid
      character(kind=c_char) :: name(*)
    end subroutine yac_write_basic_grid_to_file_c

  end interface

  !-------------------
  ! Reading FESOM grid
  !-------------------

  interface

    function yac_read_fesom_basic_grid_c(filename, gridname) &
      bind ( c, name='yac_read_fesom_basic_grid' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: filename(*)
      character(kind=c_char) :: gridname(*)

      type(c_ptr) :: yac_read_fesom_basic_grid_c ! struct yac_basic_grid *
    end function yac_read_fesom_basic_grid_c

    subroutine yac_read_fesom_grid_information_c( &
      filename, nbr_vertices, nbr_cells, num_vertices_per_cell, &
      cell_to_vertex, x_vertices, y_vertices, x_cells, y_cells) &
      bind ( c, name='yac_read_fesom_grid_information' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: filename(*)
      integer(kind=c_int)    :: nbr_vertices
      integer(kind=c_int)    :: nbr_cells
      type(c_ptr)            :: num_vertices_per_cell ! int **
      type(c_ptr)            :: cell_to_vertex ! int **
      type(c_ptr)            :: x_vertices ! double **
      type(c_ptr)            :: y_vertices ! double **
      type(c_ptr)            :: x_cells ! double **
      type(c_ptr)            :: y_cells ! double **
    end subroutine yac_read_fesom_grid_information_c

  end interface

  !------------------
  ! Reading ICON grid
  !------------------

  interface

    function yac_read_icon_basic_grid_c(filename, gridname) &
      bind ( c, name='yac_read_icon_basic_grid' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: filename(*)
      character(kind=c_char) :: gridname(*)

      type(c_ptr) :: yac_read_icon_basic_grid_c ! struct yac_basic_grid *
    end function yac_read_icon_basic_grid_c

    function yac_read_icon_basic_grid_parallel_c(filename, gridname, comm) &
      bind ( c, name='yac_read_icon_basic_grid_parallel_f2c' )

      use, intrinsic :: iso_c_binding
      import :: YAC_MPI_FINT_KIND

      character(kind=c_char)   :: filename(*)
      character(kind=c_char)   :: gridname(*)
      integer(kind=YAC_MPI_FINT_KIND) :: comm

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_read_icon_basic_grid_parallel_c
    end function yac_read_icon_basic_grid_parallel_c

    subroutine yac_read_icon_grid_information_c( &
      filename, num_vertices, num_cells, num_vertices_per_cell, &
      cell_to_vertex, x_vertices, y_vertices, x_cells, y_cells, &
      cell_mask) &
      bind ( c, name='yac_read_icon_grid_information' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: filename(*)
      integer(kind=c_int)    :: num_vertices
      integer(kind=c_int)    :: num_cells
      type(c_ptr)            :: num_vertices_per_cell ! int **
      type(c_ptr)            :: cell_to_vertex ! int **
      type(c_ptr)            :: x_vertices ! double **
      type(c_ptr)            :: y_vertices ! double **
      type(c_ptr)            :: x_cells ! double **
      type(c_ptr)            :: y_cells ! double **
      type(c_ptr)            :: cell_mask ! int **
    end subroutine yac_read_icon_grid_information_c

    subroutine yac_read_part_icon_grid_information_c( &
      filename, num_vertices, num_cells, num_vertices_per_cell, &
      cell_to_vertex, x_vertices, y_vertices, x_cells, y_cells, &
      global_cell_id, cell_mask, cell_core_mask, global_corner_id, &
      corner_core_mask, rank, size) &
      bind ( c, name='yac_read_part_icon_grid_information' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char)     :: filename(*)
      integer(kind=c_int)        :: num_vertices
      integer(kind=c_int)        :: num_cells
      type(c_ptr)                :: num_vertices_per_cell ! int **
      type(c_ptr)                :: cell_to_vertex ! int **
      type(c_ptr)                :: x_vertices ! double **
      type(c_ptr)                :: y_vertices ! double **
      type(c_ptr)                :: x_cells ! double **
      type(c_ptr)                :: y_cells ! double **
      type(c_ptr)                :: global_cell_id ! int **
      type(c_ptr)                :: cell_mask ! int **
      type(c_ptr)                :: cell_core_mask ! int **
      type(c_ptr)                :: global_corner_id ! int **
      type(c_ptr)                :: corner_core_mask ! int **
      integer(kind=c_int), value :: rank
      integer(kind=c_int), value :: size
    end subroutine yac_read_part_icon_grid_information_c

    subroutine yac_read_icon_grid_information_parallel_c( &
      filename, comm, num_vertices, num_cells, num_vertices_per_cell, &
      cell_to_vertex, global_cell_id, cell_owner, global_vertex_ids, &
      vertex_owner, x_vertices, y_vertices, x_cells, y_cells, &
      cell_mask) &
      bind ( c, name='yac_read_icon_grid_information_parallel_f2c' )

      use, intrinsic :: iso_c_binding
      import :: YAC_MPI_FINT_KIND

      character(kind=c_char)                 :: filename(*)
      integer(kind=YAC_MPI_FINT_KIND), value :: comm
      integer(kind=c_int)                    :: num_vertices
      integer(kind=c_int)                    :: num_cells
      type(c_ptr)                            :: num_vertices_per_cell ! int **
      type(c_ptr)                            :: cell_to_vertex ! int **
      type(c_ptr)                            :: global_cell_id ! int **
      type(c_ptr)                            :: cell_owner ! int **
      type(c_ptr)                            :: global_vertex_ids ! int **
      type(c_ptr)                            :: vertex_owner ! int **
      type(c_ptr)                            :: x_vertices ! double **
      type(c_ptr)                            :: y_vertices ! double **
      type(c_ptr)                            :: x_cells ! double **
      type(c_ptr)                            :: y_cells ! double **
      type(c_ptr)                            :: cell_mask ! int **
    end subroutine yac_read_icon_grid_information_parallel_c

    subroutine yac_delete_icon_grid_data_c( &
      cell_mask, global_cell_id, cell_core_mask, num_vertices_per_cell, &
      global_corner_id, corner_core_mask, cell_to_vertex, &
      x_cells, y_cells, x_vertices, y_vertices) &
      bind ( c, name='yac_delete_icon_grid_data' )

      use, intrinsic :: iso_c_binding

      type(c_ptr)                :: cell_mask ! int **
      type(c_ptr)                :: global_cell_id ! int **
      type(c_ptr)                :: cell_core_mask ! int **
      type(c_ptr)                :: num_vertices_per_cell ! int **
      type(c_ptr)                :: global_corner_id ! int **
      type(c_ptr)                :: corner_core_mask ! int **
      type(c_ptr)                :: cell_to_vertex ! int **
      type(c_ptr)                :: x_cells ! double **
      type(c_ptr)                :: y_cells ! double **
      type(c_ptr)                :: x_vertices ! double **
      type(c_ptr)                :: y_vertices ! double **
    end subroutine yac_delete_icon_grid_data_c

  end interface

  !-------------------
  ! Reading MPIOM grid
  !-------------------

  interface

    function yac_read_mpiom_basic_grid_c(filename, gridname) &
      bind ( c, name='yac_read_mpiom_basic_grid' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: filename(*)
      character(kind=c_char) :: gridname(*)

      type(c_ptr) :: yac_read_mpiom_basic_grid_c ! struct yac_basic_grid *
    end function yac_read_mpiom_basic_grid_c

    subroutine yac_read_mpiom_grid_information_c( &
      filename, num_vertices, num_cells, num_vertices_per_cell, &
      cell_to_vertex, x_vertices, y_vertices, x_cells, y_cells, &
      cell_mask) &
      bind ( c, name='yac_read_mpiom_grid_information' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: filename(*)
      integer(kind=c_int)    :: num_vertices
      integer(kind=c_int)    :: num_cells
      type(c_ptr)            :: num_vertices_per_cell ! int **
      type(c_ptr)            :: cell_to_vertex ! int **
      type(c_ptr)            :: x_vertices ! double **
      type(c_ptr)            :: y_vertices ! double **
      type(c_ptr)            :: x_cells ! double **
      type(c_ptr)            :: y_cells ! double **
      type(c_ptr)            :: cell_mask ! int **
    end subroutine yac_read_mpiom_grid_information_c

    subroutine yac_read_part_mpiom_grid_information_c( &
      filename, num_vertices, num_cells, num_vertices_per_cell, &
      cell_to_vertex, x_vertices, y_vertices, x_cells, y_cells, &
      global_cell_id, cell_mask, cell_core_mask, global_corner_id, &
      corner_core_mask, rank, size) &
      bind ( c, name='yac_read_part_mpiom_grid_information' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char)     :: filename(*)
      integer(kind=c_int)        :: num_vertices
      integer(kind=c_int)        :: num_cells
      type(c_ptr)                :: num_vertices_per_cell ! int **
      type(c_ptr)                :: cell_to_vertex ! int **
      type(c_ptr)                :: x_vertices ! double **
      type(c_ptr)                :: y_vertices ! double **
      type(c_ptr)                :: x_cells ! double **
      type(c_ptr)                :: y_cells ! double **
      type(c_ptr)                :: global_cell_id ! int **
      type(c_ptr)                :: cell_mask ! int **
      type(c_ptr)                :: cell_core_mask ! int **
      type(c_ptr)                :: global_corner_id ! int **
      type(c_ptr)                :: corner_core_mask ! int **
      integer(kind=c_int), value :: rank
      integer(kind=c_int), value :: size
    end subroutine yac_read_part_mpiom_grid_information_c

  end interface

  !-------------------
  ! Reading MPIOM grid
  !-------------------

  interface

    function yac_read_scrip_basic_grid_c( &
      grid_filename, mask_filename, grid_name, valid_mask_value, name, &
      use_ll_edges, cell_coord_idx, duplicated_cell_idx, &
      orig_cell_global_id, nbr_duplicated_cells) &
      bind ( c, name='yac_read_scrip_basic_grid' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char)     :: grid_filename(*)
      character(kind=c_char)     :: mask_filename(*)
      character(kind=c_char)     :: grid_name(*)
      integer(kind=c_int), value :: valid_mask_value
      character(kind=c_char)     :: name(*)
      integer(kind=c_int), value :: use_ll_edges

      integer(kind=c_size_t)     :: cell_coord_idx
      type(c_ptr)                :: duplicated_cell_idx ! size_t **
      type(c_ptr)                :: orig_cell_global_id ! yac_int **
      integer(kind=c_size_t)     :: nbr_duplicated_cells

      type(c_ptr) :: yac_read_scrip_basic_grid_c ! struct yac_basic_grid *
    end function yac_read_scrip_basic_grid_c

    function yac_read_scrip_basic_grid_parallel_c( &
      grid_filename, mask_filename, comm, grid_name, valid_mask_value, name, &
      use_ll_edges, cell_coord_idx, duplicated_cell_idx, &
      orig_cell_global_id, nbr_duplicated_cells) &
      bind ( c, name='yac_read_scrip_basic_grid_parallel_f2c' )

      use, intrinsic :: iso_c_binding
      import :: YAC_MPI_FINT_KIND

      character(kind=c_char)                 :: grid_filename(*)
      character(kind=c_char)                 :: mask_filename(*)
      integer(kind=YAC_MPI_FINT_KIND), value :: comm
      character(kind=c_char)                 :: grid_name(*)
      integer(kind=c_int), value             :: valid_mask_value
      character(kind=c_char)                 :: name(*)
      integer(kind=c_int), value             :: use_ll_edges

      integer(kind=c_size_t)                 :: cell_coord_idx
      type(c_ptr)                            :: duplicated_cell_idx ! size_t **
      type(c_ptr)                            :: orig_cell_global_id ! yac_int **
      integer(kind=c_size_t)                 :: nbr_duplicated_cells

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_read_scrip_basic_grid_parallel_c
    end function yac_read_scrip_basic_grid_parallel_c

    subroutine yac_read_scrip_grid_information_c( &
      grid_filename, mask_filename, grid_name, valid_mask_value, &
      num_vertices, num_cells, num_vertices_per_cell, cell_to_vertex, &
      x_vertices, y_vertices, x_cells, y_cells, cell_mask) &
      bind ( c, name='yac_read_scrip_grid_information' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char)     :: grid_filename(*)
      character(kind=c_char)     :: mask_filename(*)
      character(kind=c_char)     :: grid_name(*)
      integer(kind=c_int), value :: valid_mask_value
      integer(kind=c_int)        :: num_vertices
      integer(kind=c_int)        :: num_cells
      type(c_ptr)                :: num_vertices_per_cell ! int **
      type(c_ptr)                :: cell_to_vertex ! int **
      type(c_ptr)                :: x_vertices ! double **
      type(c_ptr)                :: y_vertices ! double **
      type(c_ptr)                :: x_cells ! double **
      type(c_ptr)                :: y_cells ! double **
      type(c_ptr)                :: cell_mask ! int **
    end subroutine yac_read_scrip_grid_information_c

  end interface

  !-----------------------
  ! Generating test fields
  !-----------------------

  interface

    function yac_test_func_c(lon, lat) &
      bind ( c, name='yac_test_func' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: lon
      real(kind=c_double), value :: lat

      real(kind=c_double) :: yac_test_func_c
    end function yac_test_func_c

    function yac_test_func_deg_c(lon, lat) &
      bind ( c, name='yac_test_func_deg' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: lon
      real(kind=c_double), value :: lat

      real(kind=c_double) :: yac_test_func_deg_c
    end function yac_test_func_deg_c

    function yac_test_ana_fcos_c(lon, lat) &
      bind ( c, name='yac_test_ana_fcos' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: lon
      real(kind=c_double), value :: lat

      real(kind=c_double) :: yac_test_ana_fcos_c
    end function yac_test_ana_fcos_c

    function yac_test_ana_fcossin_c(lon, lat) &
      bind ( c, name='yac_test_ana_fcossin' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: lon
      real(kind=c_double), value :: lat

      real(kind=c_double) :: yac_test_ana_fcossin_c
    end function yac_test_ana_fcossin_c

    function yac_test_one_c(lon, lat) &
      bind ( c, name='yac_test_one' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: lon
      real(kind=c_double), value :: lat

      real(kind=c_double) :: yac_test_one_c
    end function yac_test_one_c

    function yac_test_gulfstream_c(lon, lat) &
      bind ( c, name='yac_test_gulfstream' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: lon
      real(kind=c_double), value :: lat

      real(kind=c_double) :: yac_test_gulfstream_c
    end function yac_test_gulfstream_c

    function yac_test_harmonic_c(lon, lat) &
      bind ( c, name='yac_test_harmonic' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: lon
      real(kind=c_double), value :: lat

      real(kind=c_double) :: yac_test_harmonic_c
    end function yac_test_harmonic_c

    function yac_test_vortex_c(lon, lat) &
      bind ( c, name='yac_test_vortex' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: lon
      real(kind=c_double), value :: lat

      real(kind=c_double) :: yac_test_vortex_c
    end function yac_test_vortex_c

  end interface

  !------------------
  ! Writing VTK files
  !------------------

  interface

    function yac_vtk_open_c(vtk_file, filename, title) &
      bind ( c, name='yac_vtk_open' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value  :: vtk_file
      character(kind=c_char) :: filename(*)
      character(kind=c_char) :: title(*)

      type(c_ptr) :: yac_vtk_open_c ! YAC_VTK_FILE *
    end function yac_vtk_open_c

    subroutine yac_vtk_write_point_data_c(vtk_file, point_data, num_points) &
      bind ( c, name='yac_vtk_write_point_data' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value  :: vtk_file
      real(kind=c_double) :: point_data(*)
      integer(kind=c_int) :: num_points
    end subroutine yac_vtk_write_point_data_c

    subroutine yac_vtk_write_cell_data_c( &
      vtk_file, cell_corners, num_points_per_cell, num_cells) &
      bind ( c, name='yac_vtk_write_cell_data' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      integer(kind=c_int)        :: cell_corners(*)
      integer(kind=c_int)        :: num_points_per_cell(*)
      integer(kind=c_int), value :: num_cells
    end subroutine yac_vtk_write_cell_data_c

    subroutine yac_vtk_write_cell_scalars_uint_c( &
      vtk_file, scalars, num_cells, name) &
      bind ( c, name='yac_vtk_write_cell_scalars_uint' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      integer(kind=c_int)        :: scalars(*)
      integer(kind=c_int), value :: num_cells
      character(kind=c_char)     :: name(*)
    end subroutine yac_vtk_write_cell_scalars_uint_c

    subroutine yac_vtk_write_point_scalars_uint_c( &
      vtk_file, scalars, num_points, name) &
      bind ( c, name='yac_vtk_write_point_scalars_uint' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      integer(kind=c_int)        :: scalars(*)
      integer(kind=c_int), value :: num_points
      character(kind=c_char)     :: name(*)
    end subroutine yac_vtk_write_point_scalars_uint_c

    subroutine yac_vtk_write_cell_scalars_int_c( &
      vtk_file, scalars, num_cells, name) &
      bind ( c, name='yac_vtk_write_cell_scalars_int' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      integer(kind=c_int)        :: scalars(*)
      integer(kind=c_int), value :: num_cells
      character(kind=c_char)     :: name(*)
    end subroutine yac_vtk_write_cell_scalars_int_c

    subroutine yac_vtk_write_point_scalars_int_c( &
      vtk_file, scalars, num_points, name) &
      bind ( c, name='yac_vtk_write_point_scalars_int' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      integer(kind=c_int)        :: scalars(*)
      integer(kind=c_int), value :: num_points
      character(kind=c_char)     :: name(*)
    end subroutine yac_vtk_write_point_scalars_int_c

    subroutine yac_vtk_write_cell_scalars_float_c( &
      vtk_file, scalars, num_cells, name) &
      bind ( c, name='yac_vtk_write_cell_scalars_float' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      real(kind=c_float)         :: scalars(*)
      integer(kind=c_int), value :: num_cells
      character(kind=c_char)     :: name(*)
    end subroutine yac_vtk_write_cell_scalars_float_c

    subroutine yac_vtk_write_point_scalars_float_c( &
      vtk_file, scalars, num_points, name) &
      bind ( c, name='yac_vtk_write_point_scalars_float' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      real(kind=c_float)         :: scalars(*)
      integer(kind=c_int), value :: num_points
      character(kind=c_char)     :: name(*)
    end subroutine yac_vtk_write_point_scalars_float_c

    subroutine yac_vtk_write_cell_scalars_double_c( &
      vtk_file, scalars, num_cells, name) &
      bind ( c, name='yac_vtk_write_cell_scalars_double' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      real(kind=c_double)        :: scalars(*)
      integer(kind=c_int), value :: num_cells
      character(kind=c_char)     :: name(*)
    end subroutine yac_vtk_write_cell_scalars_double_c

    subroutine yac_vtk_write_point_scalars_double_c( &
      vtk_file, scalars, num_points, name) &
      bind ( c, name='yac_vtk_write_point_scalars_double' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: vtk_file
      real(kind=c_double)        :: scalars(*)
      integer(kind=c_int), value :: num_points
      character(kind=c_char)     :: name(*)
    end subroutine yac_vtk_write_point_scalars_double_c

    subroutine yac_vtk_close_c(vtk_file) bind ( c, name='yac_vtk_close' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: vtk_file
    end subroutine yac_vtk_close_c

  end interface

end module yac_utils
