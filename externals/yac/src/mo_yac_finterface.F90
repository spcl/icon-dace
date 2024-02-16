!>
!! @file mo_yac_finterface.F90
!! @brief Implementation of the YAC Fortran user interface description
!!
!! @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
!!                                 Rene Redler <rene.redler@mpimet.mpg.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Moritz Hanke <hanke@dkrz.de>
!!         Rene Redler <rene.redler@mpimet.mpg.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Moritz Hanke <hanke@dkrz.de>
!             Rene Redler <rene.redler@mpimet.mpg.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
!
! This file is part of YAC.
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!

!> \example test_dummy_coupling_dble.F90
!! This example simulates a whole model setup with three components (ocean,
!! atmosphere, io). It uses one process for each component.

!> \example test_dummy_coupling_real.F90
!! This example simulates a whole model setup with three components (ocean,
!! atmosphere, io). It uses one process for each component.

!> \example test_dummy_coupling3_dble.F90
!! Fortran version of \ref test_dummy_coupling3_c.c

!> \example test_dummy_coupling3_real.F90
!! Fortran version of \ref test_dummy_coupling3_c.c

!> \example test_dummy_coupling5_dble.F90
!! Fortran version of \ref test_dummy_coupling5_c.c

!> \example test_dummy_coupling5_real.F90
!! Fortran version of \ref test_dummy_coupling5_c.c

!> \example test_dummy_coupling6_dble.F90
!! Fortran version of \ref test_dummy_coupling6_c.c

!> \example test_dummy_coupling6_real.F90
!! Fortran version of \ref test_dummy_coupling6_c.c

!> \example test_dummy_coupling7_dble.F90
!! This test checks the frac mask feature.

!> \example test_dummy_coupling7_real.F90
!! This test checks the frac mask feature.

!> \example test_restart_dble.F90
!! Fortran version of \ref test_restart_c.c

!> \example test_dynamic_config.F90
!! This example tests the usage of the interface for the
!! dynamic configuration.

#ifdef HAVE_CONFIG_H
! Get the definition of the 'YAC_MPI_FINT_FC_KIND' macro.
#include "config.h"
#endif

module mo_yac_finterface
  use, intrinsic :: iso_c_binding, only : c_int, c_long, &
                                        & c_long_long, c_short, c_char

  !----------------------------------------------------------------------
  !>
  !!   Constants
  !!
  !----------------------------------------------------------------------

  public

  integer, parameter :: YAC_MAX_CHARLEN = 132
  integer, parameter :: YAC_MPI_FINT_KIND = YAC_MPI_FINT_FC_KIND


  !----------------------------------------------------------------------
  !>
  !!   Point types
  !!
  !----------------------------------------------------------------------

  enum, bind(c)
   enumerator :: YAC_LOCATION_CELL   = 0
   enumerator :: YAC_LOCATION_CORNER = 1
   enumerator :: YAC_LOCATION_EDGE   = 2
  end enum

  !----------------------------------------------------------------------
  !>
  !!   Exchange types
  !!
  !----------------------------------------------------------------------

  enum, bind(c)
   enumerator :: YAC_EXCHANGE_TYPE_NONE   = 0
   enumerator :: YAC_EXCHANGE_TYPE_SOURCE = 1
   enumerator :: YAC_EXCHANGE_TYPE_TARGET = 2
  end enum

  !----------------------------------------------------------------------
  !>
  !!   Action types (see also event.h)
  !!
  !----------------------------------------------------------------------

  enum, bind(c)
   enumerator :: YAC_ACTION_NONE               = 0
   enumerator :: YAC_ACTION_REDUCTION          = 1
   enumerator :: YAC_ACTION_COUPLING           = 2
   enumerator :: YAC_ACTION_GET_FOR_RESTART    = 4
   enumerator :: YAC_ACTION_PUT_FOR_RESTART    = 5
   enumerator :: YAC_ACTION_OUT_OF_BOUND       = 8
  end enum

  !----------------------------------------------------------------------
  !>
  !!   Time units
  !!
  !----------------------------------------------------------------------
  enum, bind(c)
   enumerator :: YAC_TIME_UNIT_MILLISECOND = 0
   enumerator :: YAC_TIME_UNIT_SECOND = 1
   enumerator :: YAC_TIME_UNIT_MINUTE = 2
   enumerator :: YAC_TIME_UNIT_HOUR = 3
   enumerator :: YAC_TIME_UNIT_DAY = 4
   enumerator :: YAC_TIME_UNIT_MONTH = 5
   enumerator :: YAC_TIME_UNIT_YEAR = 6
   enumerator :: YAC_TIME_UNIT_ISO_FORMAT = 7
  end enum

  !----------------------------------------------------------------------
  !>
  !!   Reduction types (see also event.h)
  !!
  !----------------------------------------------------------------------

  enum, bind(c)
   enumerator :: YAC_REDUCTION_TIME_NONE       = 0
   enumerator :: YAC_REDUCTION_TIME_ACCUMULATE = 1
   enumerator :: YAC_REDUCTION_TIME_AVERAGE    = 2
   enumerator :: YAC_REDUCTION_TIME_MINIMUM    = 3
   enumerator :: YAC_REDUCTION_TIME_MAXIMUM    = 4
  end enum

  !----------------------------------------------------------------------
  !>
  !!   Calendar types for use in yac_fdef_calendar
  !!
  !----------------------------------------------------------------------
  enum, bind(c)
   enumerator :: YAC_CALENDAR_NOT_SET = 0
   enumerator :: YAC_PROLEPTIC_GREGORIAN = 1
   enumerator :: YAC_YEAR_OF_365_DAYS = 2
   enumerator :: YAC_YEAR_OF_360_DAYS = 3
  end enum

  !----------------------------------------------------------------------
  !>
  !!   Parameters for interpolations
  !!
  !----------------------------------------------------------------------

  enum, bind(c)
     enumerator :: YAC_AVG_ARITHMETIC = 0
     enumerator :: YAC_AVG_DIST       = 1
     enumerator :: YAC_AVG_BARY       = 2
  end enum

  enum, bind(c)
     enumerator :: YAC_NNN_AVG    = 0
     enumerator :: YAC_NNN_DIST   = 1
     enumerator :: YAC_NNN_GAUSS  = 2
     enumerator :: YAC_NNN_RBF    = 3
  end enum

  enum, bind(c)
     enumerator :: YAC_CONSERV_DESTAREA = 0
     enumerator :: YAC_CONSERV_FRACAREA = 1
  end enum

  enum, bind(c)
     enumerator :: YAC_SPMAP_AVG  = 0
     enumerator :: YAC_SPMAP_DIST = 1
  end enum

  !----------------------------------------------------------------------
  !>
  !!   Flag paramters for emitting of coupling configurations
  !!
  !----------------------------------------------------------------------

  integer :: YAC_YAML_EMITTER_DEFAULT_F
  integer :: YAC_YAML_EMITTER_JSON_F

  !----------------------------------------------------------------------
  !>
  !!   pointer data types
  !!
  !----------------------------------------------------------------------

  type yac_real_ptr
    real, pointer :: p(:)
  end type yac_real_ptr

  type yac_dble_ptr
    double precision, pointer :: p(:)
 end type yac_dble_ptr

 type :: yac_string
     character(len=:), allocatable :: string
  end type  yac_string

  !----------------------------------------------------------------------
  !>
  !!   Fortran2C interface for yac collective routines
  !!
  !----------------------------------------------------------------------


  !> \example test_mpi_handshake.F90
  !! Fortran version of \ref test_mpi_handshake_c.c

  interface yac_fmpi_handshake
    subroutine yac_fmpi_handshake ( comm, group_names, group_comms )
      import :: YAC_MAX_CHARLEN
      integer, intent(in) :: comm      !< [IN]  MPI communicator
      character(len=YAC_MAX_CHARLEN), intent(in) :: &
        group_names(:)                 !< [IN]  group names
      integer, intent(out) :: &
        group_comms(SIZE(group_names)) !< [OUT] group communicators
    end subroutine yac_fmpi_handshake
  end interface yac_fmpi_handshake

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the coupler initialisation
  !!
  !----------------------------------------------------------------------

  !> \example test_init_final.F90
  !! This contains an example on how to use yac_finit and yac_ffinalize.

  !> \example test_init_comm_final.F90
  !! This contains an example on how to use yac_finit_comm and yac_ffinalize.

  interface yac_finit_comm

     subroutine yac_finit_comm ( comm )

       integer, intent(in) :: comm !< [IN] MPI communicator

     end subroutine yac_finit_comm

     subroutine yac_finit_comm_instance ( comm, yac_instance_id )

       integer, intent(in)  :: comm            !< [IN] MPI communicator
       integer, intent(out) :: yac_instance_id !< [OUT] id of the YAC instance initialised by this call

     end subroutine yac_finit_comm_instance

  end interface yac_finit_comm

  interface yac_finit

     subroutine yac_finit ( )

     end subroutine yac_finit

     subroutine yac_finit_instance ( yac_instance_id)

       integer, intent(out) :: yac_instance_id !< [OUT] id of the YAC instance initialised by this call

     end subroutine yac_finit_instance

  end interface yac_finit

  interface yac_finit_dummy
     subroutine yac_finit_dummy ()
     end subroutine yac_finit_dummy

     subroutine yac_finit_comm_dummy ( comm )

       integer, intent(in) :: comm !< [IN] MPI communicator

     end subroutine yac_finit_comm_dummy
  end interface yac_finit_dummy


  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the reading of configuration files
  !!
  !----------------------------------------------------------------------

  interface yac_fread_config_yaml
     subroutine yac_fread_config_yaml (yaml_filename)
       character(len=*), intent(in) :: yaml_filename !< [IN] YAML configuration file
     end subroutine yac_fread_config_yaml

     subroutine yac_fread_config_yaml_instance (yac_instance_id, yaml_filename)
       integer, intent(in) :: yac_instance_id        !< [IN] YAC instance id
       character(len=*), intent(in) :: yaml_filename !< [IN] YAML configuration file
     end subroutine yac_fread_config_yaml_instance
  end interface yac_fread_config_yaml

  interface yac_fread_config_json
     subroutine yac_fread_config_json (json_filename)
       character(len=*), intent(in) :: json_filename !< [IN] JSON configuration file
     end subroutine yac_fread_config_json

     subroutine yac_fread_config_json_instance (yac_instance_id, json_filename)
       integer, intent(in) :: yac_instance_id        !< [IN] YAC instance id
       character(len=*), intent(in) :: json_filename !< [IN] JSON configuration file
     end subroutine yac_fread_config_json_instance
  end interface yac_fread_config_json

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the coupler cleanup before restart
  !!
  !----------------------------------------------------------------------

  interface yac_fcleanup

     subroutine yac_fcleanup ()
     end subroutine yac_fcleanup

     subroutine yac_fcleanup_instance ( yac_instance_id )

      integer, intent(in) :: yac_instance_id !< [IN] YAC instance identifier

     end subroutine yac_fcleanup_instance

  end interface yac_fcleanup

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the coupler termination
  !!
  !----------------------------------------------------------------------

  interface yac_ffinalize

     subroutine yac_ffinalize ()
     end subroutine yac_ffinalize

     subroutine yac_ffinalize_instance ( yac_instance_id )

      integer, intent(in) :: yac_instance_id !< [IN] YAC instance identifier

     end subroutine

  end interface yac_ffinalize

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for getting the yac version
  !!
  !----------------------------------------------------------------------

  !> \example test_version.F90
  !! Test for correct setting of YAC version.

  interface
     function yac_fget_version () result (version_string)
     character(len=:), ALLOCATABLE :: version_string
     end function yac_fget_version
  end interface

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the component definition
  !!
  !----------------------------------------------------------------------

  INTERFACE yac_fpredef_comp

     SUBROUTINE yac_fpredef_comp ( comp_name, comp_id )

       character(len=*), intent(in) :: comp_name !< [IN]   comp name
       integer, intent(out)         :: comp_id   !< [OUT]  component identifier

     END SUBROUTINE yac_fpredef_comp

     SUBROUTINE yac_fpredef_comp_instance ( yac_instance_id, comp_name, comp_id )

       integer, intent(in)          :: yac_instance_id !< [IN] YAC instance identifier
       character(len=*), intent(in) :: comp_name       !< [IN]   comp name
       integer, intent(out)         :: comp_id         !< [OUT]  component identifier

     END SUBROUTINE yac_fpredef_comp_instance

  END INTERFACE yac_fpredef_comp

  !> \example test_def_comps.F90
  !! This contains an example on how to use yac_fdef_comp.

  interface yac_fdef_comp

     subroutine yac_fdef_comp ( comp_name, comp_id )

       character(len=*), intent(in) :: comp_name !< [IN]   comp name
       integer, intent(out)         :: comp_id   !< [OUT]  component identifier

     end subroutine yac_fdef_comp

     subroutine yac_fdef_comp_instance ( yac_instance_id, comp_name, comp_id )

       integer, intent(in)          :: yac_instance_id !< [IN] YAC instance identifier
       character(len=*), intent(in) :: comp_name       !< [IN]   comp name
       integer, intent(out)         :: comp_id         !< [OUT]  component identifier

     end subroutine yac_fdef_comp_instance

  end interface yac_fdef_comp

  interface yac_fdef_comps

     subroutine yac_fdef_comps ( comp_names, num_comps, comp_ids )

        use, intrinsic :: iso_c_binding, only : c_char
        import :: YAC_MAX_CHARLEN

        integer, intent(in)   :: num_comps
        character(kind=c_char, len=YAC_MAX_CHARLEN), intent(in) :: &
                                 comp_names(num_comps)
        integer, intent(out)  :: comp_ids(num_comps)

     end subroutine yac_fdef_comps

     subroutine yac_fdef_comps_instance ( yac_instance_id, &
                                          comp_names,      &
                                          num_comps,       &
                                          comp_ids )

        use, intrinsic :: iso_c_binding, only : c_char
        import :: YAC_MAX_CHARLEN

        integer, intent(in)   :: yac_instance_id
        integer, intent(in)   :: num_comps
        character(kind=c_char, len=YAC_MAX_CHARLEN), intent(in) :: &
                                 comp_names(num_comps)
        integer, intent(out)  :: comp_ids(num_comps)

     end subroutine yac_fdef_comps_instance

  end interface yac_fdef_comps

  !> \example test_def_datetime.F90
  !! This shows how to use yac_fdef_datetime.

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the definition of time parameters
  !!
  !----------------------------------------------------------------------

  interface yac_fdef_datetime

     subroutine yac_fdef_datetime ( start_datetime, end_datetime )

       character(len=*), intent(in), optional :: start_datetime !< [IN] start datetime of job
       character(len=*), intent(in), optional :: end_datetime   !< [IN] end datetime of job

     end subroutine yac_fdef_datetime

     subroutine yac_fdef_datetime_instance ( yac_instance_id, &
                                             start_datetime,  &
                                             end_datetime )

       integer, intent(in)                    :: yac_instance_id !< [IN] YAC instance identifier
       character(len=*), intent(in), optional :: start_datetime  !< [IN] start datetime of job
       character(len=*), intent(in), optional :: end_datetime    !< [IN] end datetime of job

     end subroutine yac_fdef_datetime_instance

  end interface yac_fdef_datetime

  interface yac_fdef_calendar

     subroutine yac_fdef_calendar( calendar )
       integer, intent(in)     :: calendar
     end subroutine yac_fdef_calendar

  end interface yac_fdef_calendar

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for getting the start- and end datetime
  !!
  !----------------------------------------------------------------------

  interface yac_fget_start_datetime

     function yac_fget_start_datetime () result (start_datetime_string)
     character(len=:), ALLOCATABLE :: start_datetime_string
     end function yac_fget_start_datetime

     function yac_fget_start_datetime_instance (yac_instance_id) &
       result (start_datetime_string)
     integer, intent(in)            :: yac_instance_id !< [IN] YAC instance identifier
     character(len=:), ALLOCATABLE :: start_datetime_string
     end function yac_fget_start_datetime_instance

  end interface yac_fget_start_datetime

  interface yac_fget_end_datetime

     function yac_fget_end_datetime () result (end_datetime_string)
     character(len=:), ALLOCATABLE :: end_datetime_string
     end function yac_fget_end_datetime

     function yac_fget_end_datetime_instance (yac_instance_id) &
       result (end_datetime_string)
     integer, intent(in)            :: yac_instance_id !< [IN] YAC instance identifier
     character(len=:), ALLOCATABLE :: end_datetime_string
     end function yac_fget_end_datetime_instance

  end interface yac_fget_end_datetime

  !> \example test_def_grid.F90
  !! This shows how to use yac_fdef_grid

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the definition of grids
  !!
  !----------------------------------------------------------------------

  interface yac_fdef_grid

     !> Definition of a non-uniform unstructured grid (cells have
     !! varying numbers of vertices)
     !! @param[in]  grid_name             grid name
     !! @param[in]  nbr_vertices          number of vertices
     !! @param[in]  nbr_cells             number of cells
     !! @param[in]  nbr_connections       total size of cell_to_vertex
     !! @param[in]  nbr_vertices_per_cell number of vertices for each cell
     !! @param[in]  x_vertices            longitudes of vertices
     !! @param[in]  y_vertices            latitudes of vertices
     !! @param[in]  cell_to_vertex        connectivity between vertices and cells\n
     !!                                   (the vertex indices per cell have to be in
     !!                                   clockwise or counterclockwise ordering)
     !! @param[out] grid_id               grid identifier
     !! @param[in]  use_ll_edges          use lonlat edges
     !!
     !! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
     !!         and determine whether they are on circles of longitudes
     !!         (same x coordinate) or latitudes (same y coordinate). An edge that
     !!         does not fullfill this condition will cause an error.
     subroutine yac_fdef_grid_nonuniform_real ( grid_name,             &
                                                nbr_vertices,          &
                                                nbr_cells,             &
                                                nbr_connections,       &
                                                nbr_vertices_per_cell, &
                                                x_vertices,            &
                                                y_vertices,            &
                                                cell_to_vertex,        &
                                                grid_id,               &
                                                use_ll_edges)

       character(len=*), intent(in)  :: grid_name
       integer, intent(in)           :: nbr_vertices
       integer, intent(in)           :: nbr_cells
       integer, intent(in)           :: nbr_connections
       integer, intent(in)           :: nbr_vertices_per_cell(nbr_cells)
       real, intent(in)              :: x_vertices(nbr_vertices)
       real, intent(in)              :: y_vertices(nbr_vertices)
       integer, intent(in)           :: cell_to_vertex(nbr_connections)
       integer, intent(out)          :: grid_id
       logical, optional, intent(in) :: use_ll_edges

     end subroutine yac_fdef_grid_nonuniform_real

     !> Definition of a non-uniform unstructured grid (cells have
     !! varying numbers of vertices)
     !! @param[in]  grid_name             grid name
     !! @param[in]  nbr_vertices          number of vertices
     !! @param[in]  nbr_cells             number of cells
     !! @param[in]  nbr_connections       total size of cell_to_vertex
     !! @param[in]  nbr_vertices_per_cell number of vertices for each cell
     !! @param[in]  x_vertices            longitudes of vertices
     !! @param[in]  y_vertices            latitudes of vertices
     !! @param[in]  cell_to_vertex        connectivity between vertices and cells\n
     !!                                   (the vertex indices per cell have to be in
     !!                                   clockwise or counterclockwise ordering)
     !! @param[out] grid_id               grid identifier
     !! @param[in]  use_ll_edges          use lonlat edges
     !!
     !! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
     !!         and determine whether they are on circles of longitudes
     !!         (same x coordinate) or latitudes (same y coordinate). An edge that
     !!         does not fullfill this condition will cause an error.
     subroutine yac_fdef_grid_nonuniform_dble ( grid_name,             &
                                                nbr_vertices,          &
                                                nbr_cells,             &
                                                nbr_connections,       &
                                                nbr_vertices_per_cell, &
                                                x_vertices,            &
                                                y_vertices,            &
                                                cell_to_vertex,        &
                                                grid_id,               &
                                                use_ll_edges)

       character(len=*), intent(in)  :: grid_name
       integer, intent(in)           :: nbr_vertices
       integer, intent(in)           :: nbr_cells
       integer, intent(in)           :: nbr_connections
       integer, intent(in)           :: nbr_vertices_per_cell(nbr_cells)
       double precision, intent(in)  :: x_vertices(nbr_vertices)
       double precision, intent(in)  :: y_vertices(nbr_vertices)
       integer, intent(in)           :: cell_to_vertex(nbr_connections)
       integer, intent(out)          :: grid_id
       logical, optional, intent(in) :: use_ll_edges

     end subroutine yac_fdef_grid_nonuniform_dble

     !> Definition of a uniform unstructured grid (all cells have the
     !! number of vertices)
     !! @param[in]  grid_name             grid name
     !! @param[in]  nbr_vertices          number of vertices
     !! @param[in]  nbr_cells             number of cells
     !! @param[in]  nbr_vertices_per_cell number of vertices for each cell
     !! @param[in]  x_vertices            longitudes of vertices
     !! @param[in]  y_vertices            latitudes of vertices
     !! @param[in]  cell_to_vertex        connectivity between vertices and cells\n
     !!                                   (the vertex indices per cell have to be
     !!                                   in clockwise or counterclockwise ordering)
     !! @param[out] grid_id               grid identifier
     !! @param[in]  use_ll_edges          use lonlat edges
     !!
     !! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
     !!         and determine whether they are on circles of longitudes
     !!         (same x coordinate) or latitudes (same y coordinate). An edge that
     !!         does not fullfill this condition will cause an error.
     subroutine yac_fdef_grid_unstruct_real ( grid_name,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              nbr_vertices_per_cell, &
                                              x_vertices,            &
                                              y_vertices,            &
                                              cell_to_vertex,        &
                                              grid_id,               &
                                              use_ll_edges)

       character(len=*), intent(in)  :: grid_name
       integer, intent(in)           :: nbr_vertices
       integer, intent(in)           :: nbr_cells
       integer, intent(in)           :: nbr_vertices_per_cell
       real, intent(in)              :: x_vertices(nbr_vertices)
       real, intent(in)              :: y_vertices(nbr_vertices)
       integer, intent(in)           :: cell_to_vertex( &
                                          nbr_vertices_per_cell,nbr_cells)
       integer, intent(out)          :: grid_id
       logical, optional, intent(in) :: use_ll_edges

     end subroutine yac_fdef_grid_unstruct_real

     !> Definition of a uniform unstructured grid (all cells have the
     !! number of vertices)
     !! @param[in]  grid_name             grid name
     !! @param[in]  nbr_vertices          number of vertices
     !! @param[in]  nbr_cells             number of cells
     !! @param[in]  nbr_vertices_per_cell number of vertices for each cell
     !! @param[in]  x_vertices            longitudes of vertices
     !! @param[in]  y_vertices            latitudes of vertices
     !! @param[in]  cell_to_vertex        connectivity between vertices and cells\n
     !!                                   (the vertex indices per cell have to be
     !!                                   in clockwise or counterclockwise ordering)
     !! @param[out] grid_id               grid identifier
     !! @param[in]  use_ll_edges          use lonlat edges
     !!
     !! @remark If (use_ll_edges == .TRUE.) YAC will check all edges of the grid
     !!         and determine whether they are on circles of longitudes
     !!         (same x coordinate) or latitudes (same y coordinate). An edge that
     !!         does not fullfill this condition will cause an error.
     subroutine yac_fdef_grid_unstruct_dble ( grid_name,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              nbr_vertices_per_cell, &
                                              x_vertices,            &
                                              y_vertices,            &
                                              cell_to_vertex,        &
                                              grid_id,               &
                                              use_ll_edges)

       character(len=*), intent(in)  :: grid_name
       integer, intent(in)           :: nbr_vertices
       integer, intent(in)           :: nbr_cells
       integer, intent(in)           :: nbr_vertices_per_cell
       double precision, intent(in)  :: x_vertices(nbr_vertices)
       double precision, intent(in)  :: y_vertices(nbr_vertices)
       integer, intent(in)           :: cell_to_vertex( &
                                          nbr_vertices_per_cell,nbr_cells)
       integer, intent(out)          :: grid_id
       logical, optional, intent(in) :: use_ll_edges

     end subroutine yac_fdef_grid_unstruct_dble

     !> Definition of a 2d curvilinear grid
     !! @param[in]  grid_name    grid name
     !! @param[in]  nbr_vertices number of cells in each dimension
     !! @param[in]  cyclic       cyclic behavior of cells in each dimension
     !! @param[in]  x_vertices   longitudes of vertices
     !! @param[in]  y_vertices   latitudes of vertices
     !! @param[out] grid_id      grid identifier
     subroutine yac_fdef_grid_curve2d_real ( grid_name,             &
                                             nbr_vertices,          &
                                             cyclic,                &
                                             x_vertices,            &
                                             y_vertices,            &
                                             grid_id)

       character(len=*), intent(in) :: grid_name
       integer, intent(in)          :: nbr_vertices(2)
       integer, intent(in)          :: cyclic(2)
       real, intent(in)             :: &
         x_vertices(nbr_vertices(1),nbr_vertices(2))
       real, intent(in)             :: &
         y_vertices(nbr_vertices(1),nbr_vertices(2))
       integer, intent(out)         :: grid_id

     end subroutine yac_fdef_grid_curve2d_real

     !> Definition of a 2d curvilinear grid
     !! @param[in]  grid_name    grid name
     !! @param[in]  nbr_vertices number of cells in each dimension
     !! @param[in]  cyclic       cyclic behavior of cells in each dimension
     !! @param[in]  x_vertices   longitudes of vertices
     !! @param[in]  y_vertices   latitudes of vertices
     !! @param[out] grid_id      grid identifier
     subroutine yac_fdef_grid_curve2d_dble ( grid_name,             &
                                             nbr_vertices,          &
                                             cyclic,                &
                                             x_vertices,            &
                                             y_vertices,            &
                                             grid_id)

       character(len=*), intent(in) :: grid_name
       integer, intent(in)          :: nbr_vertices(2)
       integer, intent(in)          :: cyclic(2)
       double precision, intent(in) :: &
         x_vertices(nbr_vertices(1),nbr_vertices(2))
       double precision, intent(in) :: &
         y_vertices(nbr_vertices(1),nbr_vertices(2))
       integer, intent(out)         :: grid_id

     end subroutine yac_fdef_grid_curve2d_dble

     !> Definition of a 2d regular grid
     !! @param[in]  grid_name    grid name
     !! @param[in]  nbr_vertices number of cells in each dimension
     !! @param[in]  cyclic       cyclic behavior of cells in each dimension
     !! @param[in]  x_vertices   longitudes of vertices
     !! @param[in]  y_vertices   latitudes of vertices
     !! @param[out] grid_id      grid identifier
     subroutine yac_fdef_grid_reg2d_real ( grid_name,             &
                                           nbr_vertices,          &
                                           cyclic,                &
                                           x_vertices,            &
                                           y_vertices,            &
                                           grid_id)

       character(len=*), intent(in) :: grid_name
       integer, intent(in)          :: nbr_vertices(2)
       integer, intent(in)          :: cyclic(2)
       real, intent(in)             :: x_vertices(nbr_vertices(1))
       real, intent(in)             :: y_vertices(nbr_vertices(2))
       integer, intent(out)         :: grid_id

     end subroutine yac_fdef_grid_reg2d_real

     !> Definition of a 2d regular grid
     !! @param[in]  grid_name    grid name
     !! @param[in]  nbr_vertices number of cells in each dimension
     !! @param[in]  cyclic       cyclic behavior of cells in each dimension
     !! @param[in]  x_vertices   longitudes of vertices
     !! @param[in]  y_vertices   latitudes of vertices
     !! @param[out] grid_id      grid identifier
     subroutine yac_fdef_grid_reg2d_dble ( grid_name,             &
                                           nbr_vertices,          &
                                           cyclic,                &
                                           x_vertices,            &
                                           y_vertices,            &
                                           grid_id)

       character(len=*), intent(in) :: grid_name
       integer, intent(in)          :: nbr_vertices(2)
       integer, intent(in)          :: cyclic(2)
       double precision, intent(in) :: x_vertices(nbr_vertices(1))
       double precision, intent(in) :: y_vertices(nbr_vertices(2))
       integer, intent(out)         :: grid_id

     end subroutine yac_fdef_grid_reg2d_dble

  end interface yac_fdef_grid

  !> \example test_def_points.F90
  !! This shows how to use yac_fdef_points.

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the definition of points
  !!
  !----------------------------------------------------------------------

  interface yac_fdef_points

     subroutine yac_fdef_points_reg2d_real ( grid_id,       &
                                             nbr_points,    &
                                             location,      &
                                             x_points_real, &
                                             y_points_real, &
                                             point_id )

       integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
       integer, intent(in)  :: nbr_points(2)                !< [IN]  number of points
       integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex
       real, intent(in)     :: x_points_real(nbr_points(1)) !< [IN]  longitudes of points
       real, intent(in)     :: y_points_real(nbr_points(2)) !< [IN]  latitudes of points
       integer, intent(out) :: point_id                     !< [OUT] point identifier

     end subroutine yac_fdef_points_reg2d_real

     subroutine yac_fdef_points_reg2d_dble ( grid_id,       &
                                             nbr_points,    &
                                             location,      &
                                             x_points,      &
                                             y_points,      &
                                             point_id )

       integer, intent(in)  :: grid_id                         !< [IN]  grid identifier
       integer, intent(in)  :: nbr_points(2)                   !< [IN]  number of points
       integer, intent(in)  :: location                        !< [IN]  location, one of center/edge/vertex
       double precision, intent(in) :: x_points(nbr_points(1)) !< [IN]  longitudes of points
       double precision, intent(in) :: y_points(nbr_points(2)) !< [IN]  latitudes of points
       integer, intent(out) :: point_id                        !< [OUT] point identifier

     end subroutine yac_fdef_points_reg2d_dble

     subroutine yac_fdef_points_curve2d_real ( grid_id,       &
                                               nbr_points,    &
                                               location,      &
                                               x_points_real, &
                                               y_points_real, &
                                               point_id )

       integer, intent(in)  :: grid_id              !< [IN]  grid identifier
       integer, intent(in)  :: nbr_points(2)        !< [IN]  number of points
       integer, intent(in)  :: location             !< [IN]  location, one of center/edge/vertex
       real, intent(in)     :: &
         x_points_real(nbr_points(1),nbr_points(2)) !< [IN]  longitudes of points
       real, intent(in)     :: &
         y_points_real(nbr_points(1),nbr_points(2)) !< [IN]  latitudes of points
       integer, intent(out) :: point_id             !< [OUT] point identifier

     end subroutine yac_fdef_points_curve2d_real

     subroutine yac_fdef_points_curve2d_dble ( grid_id,       &
                                               nbr_points,    &
                                               location,      &
                                               x_points,      &
                                               y_points,      &
                                               point_id )

       integer, intent(in)  :: grid_id         !< [IN]  grid identifier
       integer, intent(in)  :: nbr_points(2)   !< [IN]  number of points
       integer, intent(in)  :: location        !< [IN]  location, one of center/edge/vertex
       double precision, intent(in) :: &
         x_points(nbr_points(1),nbr_points(2)) !< [IN]  longitudes of points
       double precision, intent(in) :: &
         y_points(nbr_points(1),nbr_points(2)) !< [IN]  latitudes of points
       integer, intent(out) :: point_id        !< [OUT] point identifier

     end subroutine yac_fdef_points_curve2d_dble

     subroutine yac_fdef_points_unstruct_real ( grid_id,       &
                                                nbr_points,    &
                                                location,      &
                                                x_points_real, &
                                                y_points_real, &
                                                point_id )

       integer, intent(in)  :: grid_id                   !< [IN]  grid identifier
       integer, intent(in)  :: nbr_points                !< [IN]  number of points
       integer, intent(in)  :: location                  !< [IN]  location, one of center/edge/vertex
       real, intent(in)     :: x_points_real(nbr_points) !< [IN]  longitudes of points
       real, intent(in)     :: y_points_real(nbr_points) !< [IN]  latitudes of points
       integer, intent(out) :: point_id                  !< [OUT] point identifier

     end subroutine yac_fdef_points_unstruct_real

     subroutine yac_fdef_points_unstruct_dble ( grid_id,       &
                                                nbr_points,    &
                                                location,      &
                                                x_points,      &
                                                y_points,      &
                                                point_id )

       integer, intent(in)  :: grid_id                      !< [IN]  grid identifier
       integer, intent(in)  :: nbr_points                   !< [IN]  number of points
       integer, intent(in)  :: location                     !< [IN]  location, one of center/edge/vertex
       double precision, intent(in) :: x_points(nbr_points) !< [IN]  longitudes of points
       double precision, intent(in) :: y_points(nbr_points) !< [IN]  latitudes of points
       integer, intent(out) :: point_id                     !< [OUT] point identifier

     end subroutine yac_fdef_points_unstruct_dble

  end interface yac_fdef_points

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the setting of grid global ids
  !!
  !----------------------------------------------------------------------

  interface

     subroutine yac_fset_global_index ( global_index, &
                                        location,     &
                                        grid_id )

       integer, intent(in)  :: global_index(*)
       integer, intent(in)  :: location
       integer, intent(in)  :: grid_id

     end subroutine yac_fset_global_index

  end interface

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the setting of a grid core masks
  !!
  !----------------------------------------------------------------------

  interface yac_fset_core_mask

     subroutine yac_fset_core_imask ( is_core,  &
                                      location, &
                                      grid_id )

       integer, intent(in)  :: is_core(*)
       integer, intent(in)  :: location
       integer, intent(in)  :: grid_id

     end subroutine yac_fset_core_imask

     subroutine yac_fset_core_lmask ( is_core,  &
                                      location, &
                                      grid_id )

       logical, intent(in)  :: is_core(*) !< [IN] mask array
       integer, intent(in)  :: location   !< [IN] location
       integer, intent(in)  :: grid_id    !< [IN] grid identifier

     end subroutine yac_fset_core_lmask

  end interface yac_fset_core_mask

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the setting of default pointset masks
  !!
  !----------------------------------------------------------------------

  interface yac_fset_mask

     subroutine yac_fset_imask ( is_valid, points_id )

       integer, intent(in)  :: is_valid(*)  !< [IN]  mask array
       integer, intent(in)  :: points_id    !< [IN]  associated points id

     end subroutine yac_fset_imask

     subroutine yac_fset_lmask ( is_valid, points_id )

       logical, intent(in)  :: is_valid(*)  !< [IN]  mask array
       integer, intent(in)  :: points_id    !< [IN]  associated points id

     end subroutine yac_fset_lmask

  end interface

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the definition of masks
  !!
  !----------------------------------------------------------------------

  !> \example test_def_mask.F90
  !! This contains an example on how to use yac_fdef_mask.

  interface yac_fdef_mask

     subroutine yac_fdef_imask ( grid_id,    &
                                 nbr_points, &
                                 location,   &
                                 is_valid,   &
                                 mask_id )

       integer, intent(in)  :: grid_id     !< [IN] grid identifier
       integer, intent(in)  :: nbr_points  !< [IN] number of points
       integer, intent(in)  :: location    !< [IN] location, one of center/edge/vertex
       integer, intent(in)  :: is_valid(*) !< [IN] logical mask
                                           !< false, point is masked out
                                           !< true, point is valid
       integer, intent(out) :: mask_id     !< [OUT] mask identifier

     end subroutine yac_fdef_imask

     subroutine yac_fdef_lmask ( grid_id,    &
                                 nbr_points, &
                                 location,   &
                                 is_valid,   &
                                 mask_id )

       integer, intent(in)  :: grid_id     !< [IN] grid identifier
       integer, intent(in)  :: nbr_points  !< [IN] number of points
       integer, intent(in)  :: location    !< [IN] location, one of center/edge/vertex
       logical, intent(in)  :: is_valid(*) !< [IN] logical mask
                                           !< false, point is masked out
                                           !< true, point is valid
       integer, intent(out) :: mask_id     !< [OUT] mask identifier

     end subroutine yac_fdef_lmask

  end interface

  interface yac_fdef_mask_named

     subroutine yac_fdef_imask_named ( grid_id,    &
                                       nbr_points, &
                                       location,   &
                                       is_valid,   &
                                       name,       &
                                       mask_id )

       integer, intent(in)  :: grid_id      !< [IN] grid identifier
       integer, intent(in)  :: nbr_points   !< [IN] number of points
       integer, intent(in)  :: location     !< [IN] location, one of center/edge/vertex
       integer, intent(in)  :: is_valid(*)  !< [IN] logical mask
                                            !< false, point is masked out
                                            !< true, point is valid
       character(len=*), intent(in) :: name !< [IN] name of the mask
       integer, intent(out) :: mask_id      !< [OUT] mask identifier

     end subroutine yac_fdef_imask_named

     subroutine yac_fdef_lmask_named ( grid_id,    &
                                       nbr_points, &
                                       location,   &
                                       is_valid,   &
                                       name,       &
                                       mask_id )

       integer, intent(in)  :: grid_id      !< [IN] grid identifier
       integer, intent(in)  :: nbr_points   !< [IN] number of points
       integer, intent(in)  :: location     !< [IN] location, one of center/edge/vertex
       logical, intent(in)  :: is_valid(*)  !< [IN] logical mask
                                            !< false, point is masked out
                                            !< true, point is valid
       character(len=*), intent(in) :: name !< [IN] name of the mask
       integer, intent(out) :: mask_id      !< [OUT] mask identifier

     end subroutine yac_fdef_lmask_named

  end interface

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the definition of coupling fields using
  !!   default masks
  !!
  !----------------------------------------------------------------------

  interface

     subroutine yac_fdef_field ( field_name,     &
                                 component_id,   &
                                 point_ids,      &
                                 num_pointsets,  &
                                 collection_size,&
                                 timestep,       &
                                 time_unit,      &
                                 field_id )

       character(len=*), intent (in) :: field_name      !< [IN] short name of the field
       integer, intent (in)          :: component_id    !< [IN] component identifier
       integer, intent (in)          :: point_ids(*)    !< [IN] point identifier
       integer, intent (in)          :: num_pointsets   !< [IN] number of pointsets per grid
       integer, intent (in)          :: collection_size !< [IN] collection size
       character(len=*), intent (in) :: timestep        !< [IN] timestep
       integer, intent (in)          :: time_unit       !< [IN] unit of timestep
       integer, intent (out)         :: field_id        !< [OUT] returned field handle

     end subroutine yac_fdef_field

  end interface

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for the definition of coupling fields using
  !!   explicit masks
  !!
  !----------------------------------------------------------------------

  interface

     subroutine yac_fdef_field_mask ( field_name,     &
                                      component_id,   &
                                      point_ids,      &
                                      mask_ids,       &
                                      num_pointsets,  &
                                      collection_size,&
                                      timestep,       &
                                      time_unit,      &
                                      field_id )

       character(len=*), intent (in) :: field_name      !< [IN] short name of the field
       integer, intent (in)          :: component_id    !< [IN] component identifier
       integer, intent (in)          :: point_ids(*)    !< [IN] point identifier
       integer, intent (in)          :: mask_ids(*)     !< [IN] mask identifier
       integer, intent (in)          :: num_pointsets   !< [IN] number of pointsets per grid
       integer, intent (in)          :: collection_size !< [IN] collection size
       character(len=*), intent (in) :: timestep        !< [IN] timestep
       integer, intent (in)          :: time_unit       !< [IN] unit of timestep
       integer, intent (out)         :: field_id        !< [OUT] returned field handle

     end subroutine yac_fdef_field_mask

  end interface

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for checking the dimensions of a field
  !!
  !----------------------------------------------------------------------

  interface

    subroutine yac_fcheck_field_dimensions( field_id,          &
                                            collection_size,   &
                                            num_interp_fields, &
                                            interp_field_sizes )

      integer, intent (in) :: field_id                              !<[IN] field handle
      integer, intent (in) :: collection_size                       !<[IN] collection size
      integer, intent (in) :: num_interp_fields                     !<[IN] number of interpolation fields
                                                                    !!     (number of pointsets)
      integer, intent (in) :: interp_field_sizes(num_interp_fields) !<[IN] data size of each
                                                                    !!     interpolation field

    end subroutine yac_fcheck_field_dimensions

  end interface

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for sending coupling fields
  !!
  !----------------------------------------------------------------------

  interface yac_fput

     subroutine yac_fput_real ( field_id,         &
                                nbr_hor_points,   &
                                nbr_pointsets,    &
                                collection_size,  &
                                send_field,       &
                                info,             &
                                ierror )

       integer, intent (in)          :: field_id        !< [IN] field identifier
       integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)          :: nbr_pointsets   !< [IN] number of point sets
       integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
       real, intent (in)             :: send_field(nbr_hor_points, &
                                                   nbr_pointsets,  &
                                                   collection_size)
                                                        !< [IN] send field
       integer, intent (out)         :: info            !< [OUT] returned info argument
       integer, intent (out)         :: ierror          !< [OUT] returned error handler

     end subroutine yac_fput_real

     subroutine yac_fput_real_ptr ( field_id,         &
                                    nbr_pointsets,    &
                                    collection_size,  &
                                    send_field,       &
                                    info,             &
                                    ierror )

       import :: yac_real_ptr

       integer, intent (in)     :: field_id        !< [IN] field identifier
       integer, intent (in)     :: nbr_pointsets   !< [IN] number of point sets
       integer, intent (in)     :: collection_size !< [IN] number of vertical level or bundles
       type(yac_real_ptr), intent (in) :: send_field(nbr_pointsets, collection_size)
                                                   !< [IN] send field
       integer, intent (out)    :: info            !< [OUT] returned info
       integer, intent (out)    :: ierror          !< [OUT] returned error

     end subroutine yac_fput_real_ptr

     subroutine yac_fput_single_pointset_real ( field_id,         &
                                                nbr_hor_points,   &
                                                collection_size,  &
                                                send_field,       &
                                                info,             &
                                                ierror )

       integer, intent (in)          :: field_id        !< [IN] field identifier
       integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
       real, intent (in)             :: send_field(nbr_hor_points, collection_size)
                                                        !< [IN] send field
       integer, intent (out)         :: info            !< [OUT] returned info argument
       integer, intent (out)         :: ierror          !< [OUT] returned error handler

     end subroutine yac_fput_single_pointset_real

     subroutine yac_fput_dble ( field_id,         &
                                nbr_hor_points,   &
                                nbr_pointsets,    &
                                collection_size,  &
                                send_field,       &
                                info,             &
                                ierror )

       integer, intent (in)          :: field_id        !< [IN] field identifier
       integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)          :: nbr_pointsets   !< [IN] number of point sets
       integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
       double precision, intent (in) :: send_field(nbr_hor_points, nbr_pointsets, collection_size)
                                                        !< [IN] send field
       integer, intent (out)         :: info            !< [OUT] returned info argument
       integer, intent (out)         :: ierror          !< [OUT] returned error handler

     end subroutine yac_fput_dble

     subroutine yac_fput_dble_ptr ( field_id,         &
                                    nbr_pointsets,    &
                                    collection_size,  &
                                    send_field,       &
                                    info,             &
                                    ierror )

         import :: yac_dble_ptr

         integer, intent (in)     :: field_id        !< [IN] field identifier
         integer, intent (in)     :: nbr_pointsets   !< [IN] number of point sets
         integer, intent (in)     :: collection_size !< [IN] number of vertical level or bundles
         type(yac_dble_ptr), intent (in) :: send_field(nbr_pointsets, collection_size)
                                                     !< [IN] send field
         integer, intent (out)    :: info            !< [OUT] returned info
         integer, intent (out)    :: ierror          !< [OUT] returned error

     end subroutine yac_fput_dble_ptr

     subroutine yac_fput_single_pointset_dble ( field_id,         &
                                                nbr_hor_points,   &
                                                collection_size,  &
                                                send_field,       &
                                                info,             &
                                                ierror )

       integer, intent (in)          :: field_id        !< [IN] field identifier
       integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
       double precision, intent (in) :: send_field(nbr_hor_points, collection_size)
                                                        !< [IN] send field
       integer, intent (out)         :: info            !< [OUT] returned info argument
       integer, intent (out)         :: ierror          !< [OUT] returned error handler

     end subroutine yac_fput_single_pointset_dble

     subroutine yac_fput_frac_real ( field_id,         &
                                     nbr_hor_points,   &
                                     nbr_pointsets,    &
                                     collection_size,  &
                                     send_field,       &
                                     send_frac_mask,   &
                                     info,             &
                                     ierror )

       integer, intent (in)          :: field_id        !< [IN] field identifier
       integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)          :: nbr_pointsets   !< [IN] number of point sets
       integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
       real, intent (in)             :: send_field(nbr_hor_points, &
                                                   nbr_pointsets,  &
                                                   collection_size)
                                                        !< [IN] send field
       real, intent (in)             :: send_frac_mask(nbr_hor_points, &
                                                       nbr_pointsets,  &
                                                       collection_size)
                                                        !< [IN] fractional mask
       integer, intent (out)         :: info            !< [OUT] returned info argument
       integer, intent (out)         :: ierror          !< [OUT] returned error handler

     end subroutine yac_fput_frac_real

     subroutine yac_fput_frac_real_ptr ( field_id,         &
                                         nbr_pointsets,    &
                                         collection_size,  &
                                         send_field,       &
                                         send_frac_mask,   &
                                         info,             &
                                         ierror )

       import :: yac_real_ptr

       integer, intent (in)     :: field_id        !< [IN] field identifier
       integer, intent (in)     :: nbr_pointsets   !< [IN] number of point sets
       integer, intent (in)     :: collection_size !< [IN] number of vertical level or bundles
       type(yac_real_ptr), intent (in) :: send_field(nbr_pointsets, collection_size)
                                                   !< [IN] send field
       type(yac_real_ptr), intent (in) :: send_frac_mask(nbr_pointsets, collection_size)
                                                   !< [IN] fractional mask
       integer, intent (out)    :: info            !< [OUT] returned info
       integer, intent (out)    :: ierror          !< [OUT] returned error

     end subroutine yac_fput_frac_real_ptr

     subroutine yac_fput_frac_single_pointset_real ( field_id,         &
                                                     nbr_hor_points,   &
                                                     collection_size,  &
                                                     send_field,       &
                                                     send_frac_mask,   &
                                                     info,             &
                                                     ierror )

       integer, intent (in)          :: field_id        !< [IN] field identifier
       integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
       real, intent (in)             :: send_field(nbr_hor_points, collection_size)
                                                        !< [IN] send field
       real, intent (in)             :: send_frac_mask(nbr_hor_points, collection_size)
                                                        !< [IN] fractional mask
       integer, intent (out)         :: info            !< [OUT] returned info argument
       integer, intent (out)         :: ierror          !< [OUT] returned error handler

     end subroutine yac_fput_frac_single_pointset_real

     subroutine yac_fput_frac_dble ( field_id,         &
                                     nbr_hor_points,   &
                                     nbr_pointsets,    &
                                     collection_size,  &
                                     send_field,       &
                                     send_frac_mask,   &
                                     info,             &
                                     ierror )

       integer, intent (in)          :: field_id        !< [IN] field identifier
       integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)          :: nbr_pointsets   !< [IN] number of point sets
       integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
       double precision, intent (in) :: send_field(nbr_hor_points, nbr_pointsets, collection_size)
                                                        !< [IN] send field
       double precision, intent (in) :: send_frac_mask(nbr_hor_points, nbr_pointsets, collection_size)
                                                        !< [IN] fractional mask
       integer, intent (out)         :: info            !< [OUT] returned info argument
       integer, intent (out)         :: ierror          !< [OUT] returned error handler

     end subroutine yac_fput_frac_dble

     subroutine yac_fput_frac_dble_ptr ( field_id,         &
                                         nbr_pointsets,    &
                                         collection_size,  &
                                         send_field,       &
                                         send_frac_mask,   &
                                         info,             &
                                         ierror )

         import :: yac_dble_ptr

         integer, intent (in)     :: field_id        !< [IN] field identifier
         integer, intent (in)     :: nbr_pointsets   !< [IN] number of point sets
         integer, intent (in)     :: collection_size !< [IN] number of vertical level or bundles
         type(yac_dble_ptr), intent (in) :: send_field(nbr_pointsets, collection_size)
                                                     !< [IN] send field
         type(yac_dble_ptr), intent (in) :: send_frac_mask(nbr_pointsets, collection_size)
                                                     !< [IN] fractional mask
         integer, intent (out)    :: info            !< [OUT] returned info
         integer, intent (out)    :: ierror          !< [OUT] returned error

     end subroutine yac_fput_frac_dble_ptr

     subroutine yac_fput_frac_single_pointset_dble ( field_id,         &
                                                     nbr_hor_points,   &
                                                     collection_size,  &
                                                     send_field,       &
                                                     send_frac_mask,   &
                                                     info,             &
                                                     ierror )

       integer, intent (in)          :: field_id        !< [IN] field identifier
       integer, intent (in)          :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)          :: collection_size !< [IN] number of vertical level or bundles
       double precision, intent (in) :: send_field(nbr_hor_points, collection_size)
                                                        !< [IN] send field
       double precision, intent (in) :: send_frac_mask(nbr_hor_points, collection_size)
                                                        !< [IN] fractional mask
       integer, intent (out)         :: info            !< [OUT] returned info argument
       integer, intent (out)         :: ierror          !< [OUT] returned error handler

     end subroutine yac_fput_frac_single_pointset_dble

  end interface yac_fput

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for receiving coupling fields
  !!
  !----------------------------------------------------------------------

  interface yac_fget

     subroutine yac_fget_real ( field_id,         &
                                nbr_hor_points,   &
                                collection_size,  &
                                recv_field,       &
                                info,             &
                                ierror )

       integer, intent (in)  :: field_id        !< [IN] field identifier
       integer, intent (in)  :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
       real, intent (inout)  :: recv_field(nbr_hor_points, collection_size) !< [INOUT] returned field
       integer, intent (out) :: info            !< [OUT] returned info
       integer, intent (out) :: ierror          !< [OUT] returned error handler

     end subroutine yac_fget_real

     subroutine yac_fget_real_ptr ( field_id,         &
                                    collection_size,  &
                                    recv_field,       &
                                    info,             &
                                    ierror )

       import :: yac_real_ptr

       integer, intent (in)  :: field_id        !< [IN] field identifier
       integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
       type(yac_real_ptr)      :: recv_field(collection_size) !< [OUT] returned field
       integer, intent (out) :: info            !< [OUT] returned info
       integer, intent (out) :: ierror          !< [OUT] returned error handler

     end subroutine yac_fget_real_ptr

     subroutine yac_fget_dble ( field_id,         &
                                nbr_hor_points,   &
                                collection_size,  &
                                recv_field,       &
                                info,             &
                                ierror )

       integer, intent (in)  :: field_id        !< [IN] field identifier
       integer, intent (in)  :: nbr_hor_points  !< [IN] number of horizontal points
       integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
       double precision, intent (inout) :: &
          recv_field(nbr_hor_points, collection_size) !< [INOUT] returned field
       integer, intent (out) :: info            !< [OUT] returned info
       integer, intent (out) :: ierror          !< [OUT] returned error handler

     end subroutine yac_fget_dble



     subroutine yac_fget_dble_ptr ( field_id,         &
                                    collection_size,  &
                                    recv_field,       &
                                    info,             &
                                    ierror )

       import :: yac_dble_ptr

       integer, intent (in)  :: field_id        !< [IN] field identifier
       integer, intent (in)  :: collection_size !< [IN] number of vertical level or bundles
       type(yac_dble_ptr)    :: recv_field(collection_size) !< [OUT] returned field
       integer, intent (out) :: info            !< [OUT] returned info
       integer, intent (out) :: ierror          !< [OUT] returned error handler

     end subroutine yac_fget_dble_ptr

  end interface yac_fget

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for exchanging coupling fields
  !!
  !----------------------------------------------------------------------

  interface yac_fexchange

     subroutine yac_fexchange_real ( send_field_id,       &
                                     recv_field_id,       &
                                     send_nbr_hor_points, &
                                     send_nbr_pointsets,  &
                                     recv_nbr_hor_points, &
                                     collection_size,     &
                                     send_field,          &
                                     recv_field,          &
                                     send_info,           &
                                     recv_info,           &
                                     ierror )

       integer, intent (in)  :: send_field_id       !< [IN] field identifier
       integer, intent (in)  :: recv_field_id       !< [IN] field identifier
       integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
       integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
       integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
       integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
       real, intent (in)     :: send_field(send_nbr_hor_points, &
                                           send_nbr_pointsets,  &
                                           collection_size)
                                                    !< [IN] send field
       real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                           collection_size)
                                                    !< [INOUT] returned recv field
       integer, intent (out) :: send_info           !< [OUT] returned send info
       integer, intent (out) :: recv_info           !< [OUT] returned recv info
       integer, intent (out) :: ierror              !< [OUT] returned error

     end subroutine yac_fexchange_real

     subroutine yac_fexchange_real_ptr ( send_field_id,      &
                                         recv_field_id,      &
                                         send_nbr_pointsets, &
                                         collection_size,    &
                                         send_field,         &
                                         recv_field,         &
                                         send_info,          &
                                         recv_info,          &
                                         ierror )

       import :: yac_real_ptr

       integer, intent (in)  :: send_field_id       !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
       integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
       type(yac_real_ptr), intent (in) ::                        &
                                send_field(send_nbr_pointsets, &
                                           collection_size)
                                                    !< [IN] send field
       type(yac_real_ptr)      :: recv_field(collection_size)
                                                    !< [INOUT] returned recv field
       integer, intent (out) :: send_info           !< [OUT] returned send info
       integer, intent (out) :: recv_info           !< [OUT] returned recv info
       integer, intent (out) :: ierror              !< [OUT] returned error

     end subroutine yac_fexchange_real_ptr

     subroutine yac_fexchange_single_pointset_real ( send_field_id,       &
                                                     recv_field_id,       &
                                                     send_nbr_hor_points, &
                                                     recv_nbr_hor_points, &
                                                     collection_size,     &
                                                     send_field,          &
                                                     recv_field,          &
                                                     send_info,           &
                                                     recv_info,           &
                                                     ierror )

       integer, intent (in)  :: send_field_id        !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id        !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_hor_points  !< [IN] number of horizontal send points
       integer, intent (in)  :: recv_nbr_hor_points  !< [IN] number of horizontal recv points
       integer, intent (in)  :: collection_size      !< [IN] number of vertical level or bundles
       real, intent (in)     :: send_field(send_nbr_hor_points, &
                                           collection_size)
                                                     !< [IN] send field
       real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                           collection_size)
                                                     !< [INOUT] returned recv field
       integer, intent (out) :: send_info            !< [OUT] returned send info
       integer, intent (out) :: recv_info            !< [OUT] returned recv info
       integer, intent (out) :: ierror               !< [OUT] returned error

     end subroutine yac_fexchange_single_pointset_real

     subroutine yac_fexchange_dble ( send_field_id,       &
                                     recv_field_id,       &
                                     send_nbr_hor_points, &
                                     send_nbr_pointsets,  &
                                     recv_nbr_hor_points, &
                                     collection_size,     &
                                     send_field,          &
                                     recv_field,          &
                                     send_info,           &
                                     recv_info,           &
                                     ierror )

       integer, intent (in)  :: send_field_id       !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
       integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
       integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
       integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
       double precision, intent (in) :: &
                                send_field(send_nbr_hor_points, &
                                           send_nbr_pointsets,  &
                                           collection_size)
                                                    !< [IN] send field
       double precision, intent (inout):: &
                                recv_field(recv_nbr_hor_points, &
                                           collection_size)
                                                    !< [INOUT] returned recv field
       integer, intent (out) :: send_info           !< [OUT] returned send info
       integer, intent (out) :: recv_info           !< [OUT] returned recv info
       integer, intent (out) :: ierror              !< [OUT] returned error

     end subroutine yac_fexchange_dble

     subroutine yac_fexchange_dble_ptr ( send_field_id,      &
                                         recv_field_id,      &
                                         send_nbr_pointsets, &
                                         collection_size,    &
                                         send_field,         &
                                         recv_field,         &
                                         send_info,          &
                                         recv_info,          &
                                         ierror )

       import :: yac_dble_ptr

       integer, intent (in)  :: send_field_id      !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id      !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_pointsets !< [IN] number of send point sets
       integer, intent (in)  :: collection_size    !< [IN] number of vertical level or bundles
       type(yac_dble_ptr), intent (in) ::                        &
                                send_field(send_nbr_pointsets, &
                                           collection_size)
                                                   !< [IN] send field
       type(yac_dble_ptr)      :: recv_field(collection_size)
                                                   !< [OUT] returned recv field
       integer, intent (out) :: send_info          !< [OUT] returned send info
       integer, intent (out) :: recv_info          !< [OUT] returned recv info
       integer, intent (out) :: ierror             !< [OUT] returned error

     end subroutine yac_fexchange_dble_ptr

     subroutine yac_fexchange_single_pointset_dble ( send_field_id,       &
                                                     recv_field_id,       &
                                                     send_nbr_hor_points, &
                                                     recv_nbr_hor_points, &
                                                     collection_size,     &
                                                     send_field,          &
                                                     recv_field,          &
                                                     send_info,           &
                                                     recv_info,           &
                                                     ierror )

       integer, intent (in)  :: send_field_id       !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
       integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
       integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
       double precision, intent (in) ::                         &
                                send_field(send_nbr_hor_points, &
                                           collection_size)
                                                    !< [IN] send field
       double precision, intent (inout)::                       &
                                recv_field(recv_nbr_hor_points, &
                                           collection_size)
                                                    !< [INOUT] returned recv field
       integer, intent (out) :: send_info           !< [OUT] returned send info
       integer, intent (out) :: recv_info           !< [OUT] returned recv info
       integer, intent (out) :: ierror              !< [OUT] returned error

     end subroutine yac_fexchange_single_pointset_dble

     subroutine yac_fexchange_frac_real ( send_field_id,       &
                                          recv_field_id,       &
                                          send_nbr_hor_points, &
                                          send_nbr_pointsets,  &
                                          recv_nbr_hor_points, &
                                          collection_size,     &
                                          send_field,          &
                                          send_frac_mask,      &
                                          recv_field,          &
                                          send_info,           &
                                          recv_info,           &
                                          ierror )

       integer, intent (in)  :: send_field_id       !< [IN] field identifier
       integer, intent (in)  :: recv_field_id       !< [IN] field identifier
       integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
       integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
       integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
       integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
       real, intent (in)     :: send_field(send_nbr_hor_points, &
                                           send_nbr_pointsets,  &
                                           collection_size)
                                                    !< [IN] send field
       real, intent (in)     :: send_frac_mask(send_nbr_hor_points, &
                                               send_nbr_pointsets,  &
                                               collection_size)
                                                    !< [IN] fractional mask
       real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                           collection_size)
                                                    !< [INOUT] returned recv field
       integer, intent (out) :: send_info           !< [OUT] returned send info
       integer, intent (out) :: recv_info           !< [OUT] returned recv info
       integer, intent (out) :: ierror              !< [OUT] returned error

     end subroutine yac_fexchange_frac_real

     subroutine yac_fexchange_frac_real_ptr ( send_field_id,      &
                                              recv_field_id,      &
                                              send_nbr_pointsets, &
                                              collection_size,    &
                                              send_field,         &
                                              send_frac_mask,     &
                                              recv_field,         &
                                              send_info,          &
                                              recv_info,          &
                                              ierror )

       import :: yac_real_ptr

       integer, intent (in)  :: send_field_id       !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
       integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
       type(yac_real_ptr), intent (in) ::                      &
                                send_field(send_nbr_pointsets, &
                                           collection_size)
                                                    !< [IN] send field
       type(yac_real_ptr), intent (in) ::                          &
                                send_frac_mask(send_nbr_pointsets, &
                                               collection_size)
                                                    !< [IN] fractional mask
       type(yac_real_ptr)    :: recv_field(collection_size)
                                                    !< [INOUT] returned recv field
       integer, intent (out) :: send_info           !< [OUT] returned send info
       integer, intent (out) :: recv_info           !< [OUT] returned recv info
       integer, intent (out) :: ierror              !< [OUT] returned error

     end subroutine yac_fexchange_frac_real_ptr

     subroutine yac_fexchange_frac_single_pointset_real ( send_field_id,       &
                                                          recv_field_id,       &
                                                          send_nbr_hor_points, &
                                                          recv_nbr_hor_points, &
                                                          collection_size,     &
                                                          send_field,          &
                                                          send_frac_mask,      &
                                                          recv_field,          &
                                                          send_info,           &
                                                          recv_info,           &
                                                          ierror )

       integer, intent (in)  :: send_field_id        !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id        !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_hor_points  !< [IN] number of horizontal send points
       integer, intent (in)  :: recv_nbr_hor_points  !< [IN] number of horizontal recv points
       integer, intent (in)  :: collection_size      !< [IN] number of vertical level or bundles
       real, intent (in)     :: send_field(send_nbr_hor_points, &
                                           collection_size)
                                                     !< [IN] send field
       real, intent (in)     :: send_frac_mask(send_nbr_hor_points, &
                                               collection_size)
                                                     !< [IN] fractional mask
       real, intent (inout)  :: recv_field(recv_nbr_hor_points, &
                                           collection_size)
                                                     !< [INOUT] returned recv field
       integer, intent (out) :: send_info            !< [OUT] returned send info
       integer, intent (out) :: recv_info            !< [OUT] returned recv info
       integer, intent (out) :: ierror               !< [OUT] returned error

     end subroutine yac_fexchange_frac_single_pointset_real

     subroutine yac_fexchange_frac_dble ( send_field_id,       &
                                          recv_field_id,       &
                                          send_nbr_hor_points, &
                                          send_nbr_pointsets,  &
                                          recv_nbr_hor_points, &
                                          collection_size,     &
                                          send_field,          &
                                          send_frac_mask,      &
                                          recv_field,          &
                                          send_info,           &
                                          recv_info,           &
                                          ierror )

       integer, intent (in)  :: send_field_id       !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
       integer, intent (in)  :: send_nbr_pointsets  !< [IN] number of send point sets
       integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
       integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
       double precision, intent (in) :: &
                                send_field(send_nbr_hor_points, &
                                           send_nbr_pointsets,  &
                                           collection_size)
                                                    !< [IN] send field
       double precision, intent (in) :: &
                                send_frac_mask(send_nbr_hor_points, &
                                               send_nbr_pointsets,  &
                                               collection_size)
                                                    !< [IN] fractional mask
       double precision, intent (inout):: &
                                recv_field(recv_nbr_hor_points, &
                                           collection_size)
                                                    !< [INOUT] returned recv field
       integer, intent (out) :: send_info           !< [OUT] returned send info
       integer, intent (out) :: recv_info           !< [OUT] returned recv info
       integer, intent (out) :: ierror              !< [OUT] returned error

     end subroutine yac_fexchange_frac_dble

     subroutine yac_fexchange_frac_dble_ptr ( send_field_id,      &
                                              recv_field_id,      &
                                              send_nbr_pointsets, &
                                              collection_size,    &
                                              send_field,         &
                                              send_frac_mask,     &
                                              recv_field,         &
                                              send_info,          &
                                              recv_info,          &
                                              ierror )

       import :: yac_dble_ptr

       integer, intent (in)  :: send_field_id      !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id      !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_pointsets !< [IN] number of send point sets
       integer, intent (in)  :: collection_size    !< [IN] number of vertical level or bundles
       type(yac_dble_ptr), intent (in) ::                      &
                                send_field(send_nbr_pointsets, &
                                           collection_size)
                                                   !< [IN] send field
       type(yac_dble_ptr), intent (in) ::                          &
                                send_frac_mask(send_nbr_pointsets, &
                                               collection_size)
                                                   !< [IN] fractional mask
       type(yac_dble_ptr)      :: recv_field(collection_size)
                                                   !< [OUT] returned recv field
       integer, intent (out) :: send_info          !< [OUT] returned send info
       integer, intent (out) :: recv_info          !< [OUT] returned recv info
       integer, intent (out) :: ierror             !< [OUT] returned error

     end subroutine yac_fexchange_frac_dble_ptr

     subroutine yac_fexchange_frac_single_pointset_dble ( send_field_id,       &
                                                          recv_field_id,       &
                                                          send_nbr_hor_points, &
                                                          recv_nbr_hor_points, &
                                                          collection_size,     &
                                                          send_field,          &
                                                          send_frac_mask,      &
                                                          recv_field,          &
                                                          send_info,           &
                                                          recv_info,           &
                                                          ierror )

       integer, intent (in)  :: send_field_id       !< [IN] send field identifier
       integer, intent (in)  :: recv_field_id       !< [IN] recv field identifier
       integer, intent (in)  :: send_nbr_hor_points !< [IN] number of horizontal send points
       integer, intent (in)  :: recv_nbr_hor_points !< [IN] number of horizontal recv points
       integer, intent (in)  :: collection_size     !< [IN] number of vertical level or bundles
       double precision, intent (in) ::                         &
                                send_field(send_nbr_hor_points, &
                                           collection_size)
                                                    !< [IN] send field
       double precision, intent (in) ::                             &
                                send_frac_mask(send_nbr_hor_points, &
                                               collection_size)
                                                    !< [IN] fractional mask
       double precision, intent (inout)::                       &
                                recv_field(recv_nbr_hor_points, &
                                           collection_size)
                                                    !< [INOUT] returned recv field
       integer, intent (out) :: send_info           !< [OUT] returned send info
       integer, intent (out) :: recv_info           !< [OUT] returned recv info
       integer, intent (out) :: ierror              !< [OUT] returned error

     end subroutine yac_fexchange_frac_single_pointset_dble

  end interface yac_fexchange

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for testing fields for active communicaitons
  !!
  !----------------------------------------------------------------------

  interface yac_ftest

     subroutine yac_ftest_i ( field_id, flag )

       integer, intent (in)  :: field_id !< [IN]  field identifier
       integer, intent (out) :: flag     !< [OUT] "0" if there is an uncompleted asynchronous
                                         !!       communication associated with the field,
                                         !!       "1" otherwise

     end subroutine yac_ftest_i

     subroutine yac_ftest_l ( field_id, flag )

       integer, intent (in)  :: field_id !< [IN]  field identifier
       logical, intent (out) :: flag     !< [OUT] .FLAGS. if there is an uncompleted asynchronous
                                         !!       communication associated with the field,
                                         !!       .TRUE. otherwise

     end subroutine yac_ftest_l

  end interface yac_ftest

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for getting back a local MPI communicator
  !!
  !----------------------------------------------------------------------

  interface

     subroutine yac_fget_comp_comm ( comp_id, comp_comm )

       integer, intent (in)  :: comp_id   !< [IN]  component ID
       integer, intent (out) :: comp_comm !< [OUT] component MPI communicator

     end subroutine yac_fget_comp_comm

  end interface

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for getting back a MPI communicator for
  !!   communication between components
  !!
  !----------------------------------------------------------------------

  interface yac_fget_comps_comm

     subroutine yac_fget_comps_comm ( comp_names, &
                                      num_comps,  &
                                      comps_comm)

       use, intrinsic :: iso_c_binding, only : c_char
       import :: YAC_MAX_CHARLEN

       integer, intent(in)   :: num_comps              !< [IN]  number of components
       character(kind=c_char, len=YAC_MAX_CHARLEN), intent(in) :: &
                                comp_names(num_comps)  !< [IN]  name of components
       integer, intent (out) :: comps_comm             !< [OUT] MPI communicator

     end subroutine yac_fget_comps_comm

     subroutine yac_fget_comps_comm_instance ( yac_instance_id, &
                                               comp_names,      &
                                               num_comps,       &
                                               comps_comm)

       use, intrinsic :: iso_c_binding, only : c_char
       import :: YAC_MAX_CHARLEN

       integer, intent(in)          :: yac_instance_id !< [IN]  YAC instance identifier
       integer, intent(in)   :: num_comps              !< [IN]  number of components
       character(kind=c_char, len=YAC_MAX_CHARLEN), intent(in) :: &
                                comp_names(num_comps)  !< [IN]  name of components
       integer, intent (out) :: comps_comm             !< [OUT] MPI communicator

     end subroutine yac_fget_comps_comm_instance

  end interface yac_fget_comps_comm


  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for invoking the end of the definition phase
  !!
  !----------------------------------------------------------------------

  interface yac_fsync_def

     subroutine yac_fsync_def ( )

     end subroutine yac_fsync_def

     subroutine yac_fsync_def_instance ( yac_instance_id )

       integer, intent(in)  :: yac_instance_id !< [IN]  YAC instance identifier

     end subroutine yac_fsync_def_instance

  end interface yac_fsync_def

 !----------------------------------------------------------------------
 !>
 !!   Fortran interfaces for the definition of an interpolation stack
 !!
 !----------------------------------------------------------------------

 interface

    subroutine yac_fget_interp_stack_config(interp_stack_config_id)
      integer, intent(out) :: interp_stack_config_id
    end subroutine yac_fget_interp_stack_config

    subroutine yac_ffree_interp_stack_config(interp_stack_config_id)
      integer, intent(in) :: interp_stack_config_id
    end subroutine yac_ffree_interp_stack_config

    subroutine yac_fadd_interp_stack_config_average(interp_stack_config_id, &
         reduction_type, partial_coverage)
      integer, intent(in)     :: interp_stack_config_id
      integer, intent(in)     :: reduction_type
      integer, intent(in)     :: partial_coverage
    end subroutine yac_fadd_interp_stack_config_average

    subroutine yac_fadd_interp_stack_config_nnn(interp_stack_config_id, &
         type, n, scale)
      integer, intent(in)          :: interp_stack_config_id
      integer, intent(in)          :: type
      integer, intent(in)          :: n
      double precision, intent(in) :: scale
    end subroutine yac_fadd_interp_stack_config_nnn

    subroutine yac_fadd_interp_stack_config_conservative(interp_stack_config_id, &
         order, enforced_conserv, partial_coverage, normalization)
      integer, intent(in)     :: interp_stack_config_id
      integer, intent(in)     :: order
      integer, intent(in)     :: enforced_conserv
      integer, intent(in)     :: partial_coverage
      integer, intent(in)     :: normalization
    end subroutine yac_fadd_interp_stack_config_conservative

    subroutine yac_fadd_interp_stack_config_spmap(interp_stack_config_id, &
         spread_distance, max_search_distance, weight_type)
      integer, intent(in)          :: interp_stack_config_id
      double precision, intent(in) :: spread_distance
      double precision, intent(in) :: max_search_distance
      integer, intent(in)          :: weight_type
    end subroutine yac_fadd_interp_stack_config_spmap

    subroutine yac_fadd_interp_stack_config_hcsbb(interp_stack_config_id)
      integer, intent(in)     :: interp_stack_config_id
    end subroutine yac_fadd_interp_stack_config_hcsbb

    subroutine yac_fadd_interp_stack_config_user_file(interp_stack_config_id, &
         filename, src_grid_name, tgt_grid_name)
      integer, intent(in)     :: interp_stack_config_id
      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: src_grid_name
      character (len=*), intent(in) :: tgt_grid_name
    end subroutine yac_fadd_interp_stack_config_user_file

    subroutine yac_fadd_interp_stack_config_fixed(interp_stack_config_id, &
         val)
      integer, intent(in)          :: interp_stack_config_id
      double precision, intent(in) :: val
    end subroutine yac_fadd_interp_stack_config_fixed

    subroutine yac_fadd_interp_stack_config_creep(interp_stack_config_id, &
         creep_distance)
      integer, intent(in) :: interp_stack_config_id
      integer, intent(in) :: creep_distance
    end subroutine yac_fadd_interp_stack_config_creep

 end interface

 !----------------------------------------------------------------------
 !>
 !!   Fortran interface for definition of a couple
 !!
 !----------------------------------------------------------------------

 interface yac_fdef_couple

    subroutine yac_fdef_couple (                     &
       src_comp_name, src_grid_name, src_field_name, &
       tgt_comp_name, tgt_grid_name, tgt_field_name, &
       coupling_timestep, time_unit, time_reduction, &
       interp_stack_config_id, src_lag, tgt_lag,     &
       weight_file, mapping_side, scale_factor,      &
       scale_summand, src_mask_names, tgt_mask_name )

      import :: yac_string
      character ( len=* ), intent(in)           :: src_comp_name
      character ( len=* ), intent(in)           :: src_grid_name
      character ( len=* ), intent(in)           :: src_field_name
      character ( len=* ), intent(in)           :: tgt_comp_name
      character ( len=* ), intent(in)           :: tgt_grid_name
      character ( len=* ), intent(in)           :: tgt_field_name
      character ( len=* ), intent(in)           :: coupling_timestep
      integer, intent(in)                       :: time_unit
      integer, intent(in)                       :: time_reduction
      integer, intent(in)                       :: interp_stack_config_id
      integer, intent(in), optional             :: src_lag
      integer, intent(in), optional             :: tgt_lag
      character ( len=* ), intent(in), optional :: weight_file
      integer, intent(in), optional             :: mapping_side
      double precision, intent(in), optional    :: scale_factor
      double precision, intent(in), optional    :: scale_summand
      type(yac_string), intent(in), optional    :: src_mask_names(:)
      character ( len=* ), intent(in), optional :: tgt_mask_name
    end subroutine yac_fdef_couple

    subroutine yac_fdef_couple_instance ( instance_id, &
         src_comp_name, src_grid_name, src_field_name, &
         tgt_comp_name, tgt_grid_name, tgt_field_name, &
         coupling_timestep, time_unit, time_reduction, &
         interp_stack_config_id, src_lag, tgt_lag,     &
         weight_file, mapping_side, scale_factor,      &
         scale_summand, src_mask_names, tgt_mask_name )

      import :: yac_string
      integer, intent(in)                       :: instance_id
      character ( len=* ), intent(in)           :: src_comp_name
      character ( len=* ), intent(in)           :: src_grid_name
      character ( len=* ), intent(in)           :: src_field_name
      character ( len=* ), intent(in)           :: tgt_comp_name
      character ( len=* ), intent(in)           :: tgt_grid_name
      character ( len=* ), intent(in)           :: tgt_field_name
      character ( len=* ), intent(in)           :: coupling_timestep
      integer, intent(in)                       :: time_unit
      integer, intent(in)                       :: time_reduction
      integer, intent(in)                       :: interp_stack_config_id
      integer, intent(in), optional             :: src_lag
      integer, intent(in), optional             :: tgt_lag
      character ( len=* ), intent(in), optional :: weight_file
      integer, intent(in), optional             :: mapping_side
      double precision, intent(in), optional    :: scale_factor
      double precision, intent(in), optional    :: scale_summand
      type(yac_string), intent(in), optional    :: src_mask_names(:)
      character ( len=* ), intent(in), optional :: tgt_mask_name
    end subroutine yac_fdef_couple_instance

 end interface yac_fdef_couple

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for invoking the end of the definition phase
  !!
  !----------------------------------------------------------------------

  interface yac_fenddef

    subroutine yac_fenddef ( )

    end subroutine yac_fenddef

    subroutine yac_fenddef_instance ( yac_instance_id )

      integer, intent(in)  :: yac_instance_id !< [IN] YAC instance identifier

    end subroutine yac_fenddef_instance

    subroutine yac_fenddef_and_emit_config(emit_flags, config)

      integer, intent (in)           :: emit_flags      !< [IN] flags for emitting the config
      character (len=:), ALLOCATABLE :: config          !< [IN,OUT] configuration string

    end subroutine yac_fenddef_and_emit_config

    subroutine yac_fenddef_and_emit_config_instance( &
      yac_instance_id, emit_flags, config)

      integer, intent (in)           :: yac_instance_id !< [IN] YAC instance identifier
      integer, intent (in)           :: emit_flags      !< [IN] flags for emitting the config
      character (len=:), ALLOCATABLE :: config          !< [IN,OUT] configuration string

    end subroutine yac_fenddef_and_emit_config_instance

  end interface yac_fenddef

  !----------------------------------------------------------------------
  !>
  !!   Fortran interface for invoking query functions
  !!
  !----------------------------------------------------------------------

  interface yac_fget_field_id

     function yac_fget_field_id ( comp_name, grid_name, field_name ) result(field_id)

       character(len=*), intent(in)  :: comp_name
       character(len=*), intent(in)  :: grid_name
       character(len=*), intent(in)  :: field_name
       integer :: field_id

     end function yac_fget_field_id

     function yac_fget_field_id_instance ( yac_id,    &
                                           comp_name, &
                                           grid_name, &
                                           field_name ) result(field_id)

       integer, intent(in) :: yac_id
       character(len=*), intent(in)  :: comp_name
       character(len=*), intent(in)  :: grid_name
       character(len=*), intent(in)  :: field_name
       integer :: field_id

     end function yac_fget_field_id_instance
  end interface yac_fget_field_id

  ! ---------------------------------------------------------------------

  interface yac_fget_component_name

     function yac_fget_component_name_from_field_id ( field_id ) &
          result( comp_name )

       integer, intent (in)           :: field_id
       character(len=:), ALLOCATABLE :: comp_name

     end function yac_fget_component_name_from_field_id

  end interface yac_fget_component_name

  interface yac_fget_grid_name

     function yac_fget_grid_name_from_field_id ( field_id ) &
          result( grid_name )

       integer, intent (in)           :: field_id
       character(len=:), ALLOCATABLE :: grid_name

     end function yac_fget_grid_name_from_field_id

  end interface yac_fget_grid_name

  interface yac_fget_field_name

     function yac_fget_field_name_from_field_id ( field_id ) &
          result( field_name )

       integer, intent (in)           :: field_id
       character(len=:), ALLOCATABLE :: field_name

     end function yac_fget_field_name_from_field_id

  end interface yac_fget_field_name

  ! ---------------------------------------------------------------------

  interface yac_fget_comp_names

     function yac_fget_comp_names( ) result(comp_names)
       import :: yac_string
       type(yac_string), allocatable :: comp_names(:)
     end function  yac_fget_comp_names

     function yac_fget_comp_names_instance( yac_instance_id) &
         result( comp_names )
       import :: yac_string
       integer, intent(in) :: yac_instance_id
       type(yac_string), allocatable :: comp_names(:)
     end function yac_fget_comp_names_instance

  end interface yac_fget_comp_names

  interface yac_fget_grid_names

     function yac_fget_grid_names( ) result( grid_names )
       import :: yac_string
       type(yac_string), allocatable :: grid_names(:)
     end function yac_fget_grid_names

     function yac_fget_grid_names_instance( yac_instance_id ) &
         result ( grid_names )
       import :: yac_string
       integer, intent(in) :: yac_instance_id
       type(yac_string), allocatable :: grid_names(:)
     end function yac_fget_grid_names_instance

  end interface yac_fget_grid_names

  interface yac_fget_comp_grid_names

     function yac_fget_comp_grid_names( comp_name ) result( grid_names )
       import :: yac_string
       character(len=*), intent(in)  :: comp_name
       type(yac_string), allocatable :: grid_names(:)
     end function yac_fget_comp_grid_names

     function yac_fget_comp_grid_names_instance( yac_instance_id, comp_name ) &
         result ( grid_names )
       import :: yac_string
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in)  :: comp_name
       type(yac_string), allocatable :: grid_names(:)
     end function yac_fget_comp_grid_names_instance

  end interface yac_fget_comp_grid_names

  interface yac_fget_field_names

     function yac_fget_field_names( comp_name, grid_name ) result ( field_names )
       import :: yac_string
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       type(yac_string), allocatable :: field_names(:)
     end function yac_fget_field_names

     function yac_fget_field_names_instance( yac_instance_id, &
                                             comp_name,       &
                                             grid_name )      &
          result( field_names )
       import :: yac_string
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       type(yac_string), allocatable :: field_names(:)
     end function yac_fget_field_names_instance

  end interface yac_fget_field_names

  ! ---------------------------------------------------------------------

  interface yac_fget_field_role

     function yac_fget_role_from_field_id ( field_id )

       integer, intent (in)           :: field_id
       integer                        :: yac_fget_role_from_field_id

     end function yac_fget_role_from_field_id

     function yac_fget_field_role ( comp_name, grid_name, field_name )
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       integer :: yac_fget_field_role
     end function yac_fget_field_role

     function yac_fget_field_role_instance ( yac_instance_id, &
                                             comp_name,       &
                                             grid_name,       &
                                             field_name )
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       integer :: yac_fget_field_role_instance
     end function yac_fget_field_role_instance

  end interface yac_fget_field_role

  ! ---------------------------------------------------------------------

  interface yac_fget_field_timestep


     function yac_fget_timestep_from_field_id ( field_id ) result( timestep )
       integer, intent (in)           :: field_id   !< [IN]  field name
       character(len=:), ALLOCATABLE :: timestep   !< [OUT] timestep in iso format

     end function yac_fget_timestep_from_field_id

     function yac_fget_field_timestep ( comp_name, grid_name, field_name ) &
          result( timestep )

       character(len=*), intent(in)   :: comp_name
       character(len=*), intent(in)   :: grid_name
       character(len=*), intent(in)   :: field_name
       character(len=:), ALLOCATABLE :: timestep

     end function  yac_fget_field_timestep

     function yac_fget_field_timestep_instance ( yac_instance_id, &
                                                 comp_name,       &
                                                 grid_name,       &
                                                field_name )      &
          result( timestep )

       integer, intent(in)            :: yac_instance_id
       character(len=*), intent(in)   :: comp_name
       character(len=*), intent(in)   :: grid_name
       character(len=*), intent(in)   :: field_name
       character(len=:), ALLOCATABLE :: timestep

     end function  yac_fget_field_timestep_instance

  end interface yac_fget_field_timestep

  ! ---------------------------------------------------------------------

  interface yac_fget_field_datetime

     function yac_fget_field_datetime ( field_id ) &
          result ( datetime )
       integer, intent(in) :: field_id
       character(len=:), allocatable :: datetime
     end function yac_fget_field_datetime

  end interface yac_fget_field_datetime

  ! ---------------------------------------------------------------------

  interface yac_fenable_field_frac_mask

     subroutine yac_fenable_field_frac_mask ( &
       comp_name, grid_name, field_name, frac_mask_fallback_value )
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       double precision, intent(in) :: frac_mask_fallback_value
     end subroutine yac_fenable_field_frac_mask

     subroutine yac_fenable_field_frac_mask_instance ( &
       yac_instance_id, comp_name, grid_name, field_name, &
       frac_mask_fallback_value )
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       double precision, intent(in) :: frac_mask_fallback_value
     end subroutine yac_fenable_field_frac_mask_instance

  end interface yac_fenable_field_frac_mask

  ! ---------------------------------------------------------------------

  interface yac_fget_field_frac_mask_fallback_value

     function yac_fget_field_frac_mask_fallback_value ( &
       comp_name, grid_name, field_name )
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       double precision :: yac_fget_field_frac_mask_fallback_value
     end function yac_fget_field_frac_mask_fallback_value

     function yac_fget_field_frac_mask_fallback_value_instance ( &
       yac_instance_id, comp_name, grid_name, field_name )
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       double precision :: yac_fget_field_frac_mask_fallback_value_instance
     end function yac_fget_field_frac_mask_fallback_value_instance

  end interface yac_fget_field_frac_mask_fallback_value

  ! ---------------------------------------------------------------------

  interface yac_fget_field_collection_size

     function yac_fget_collection_size_from_field_id ( field_id ) &
          RESULT(collection_size)

       integer, intent (in)  :: field_id        !< [IN]  field identifier
       integer :: collection_size !< [OUT] collection size

     end function  yac_fget_collection_size_from_field_id

     function yac_fget_field_collection_size ( comp_name, grid_name, field_name )
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       integer :: yac_fget_field_collection_size
     end function yac_fget_field_collection_size

     function yac_fget_field_collection_size_instance ( yac_instance_id, &
                                                        comp_name,       &
                                                        grid_name,       &
                                                        field_name )
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       integer :: yac_fget_field_collection_size_instance
     end function yac_fget_field_collection_size_instance

  end interface yac_fget_field_collection_size

  ! ---------------------------------------------------------------------

  interface yac_fdef_component_metadata
     subroutine yac_fdef_component_metadata(comp_name, metadata)
       character(len=*), intent(in)    :: comp_name
       character(len=*), intent(in)    :: metadata
     end subroutine yac_fdef_component_metadata

     subroutine yac_fdef_component_metadata_instance(yac_instance_id, comp_name, &
          metadata)
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in)    :: comp_name
       character(len=*), intent(in)    :: metadata
     end subroutine yac_fdef_component_metadata_instance
  end interface yac_fdef_component_metadata

  interface yac_fcomponent_has_metadata
    function yac_fcomponent_has_metadata(comp_name) result( has_metadata )
      character(len=*), intent(in) :: comp_name
      logical :: has_metadata
    end function yac_fcomponent_has_metadata

    function yac_fcomponent_has_metadata_instance(yac_instance_id, comp_name) &
      result( has_metadata )
      integer, intent(in) :: yac_instance_id
      character(len=*), intent(in) :: comp_name
      logical :: has_metadata
    end function yac_fcomponent_has_metadata_instance
  end interface yac_fcomponent_has_metadata

  interface yac_fget_component_metadata
     function yac_fget_component_metadata(comp_name) result( metadata )
       character(len=*), intent(in)  :: comp_name
       character(len=:), allocatable :: metadata
     end function yac_fget_component_metadata

     function yac_fget_component_metadata_instance(yac_instance_id, comp_name) &
          result( metadata )
       integer, intent(in)           :: yac_instance_id
       character(len=*), intent(in)  :: comp_name
       character(len=:), allocatable :: metadata
     end function yac_fget_component_metadata_instance
  end interface yac_fget_component_metadata

  interface yac_fdef_grid_metadata
     subroutine yac_fdef_grid_metadata(grid_name, metadata)
       character(len=*), intent(in)    :: grid_name
       character(len=*), intent(in)    :: metadata
     end subroutine yac_fdef_grid_metadata

     subroutine yac_fdef_grid_metadata_instance(yac_instance_id, grid_name, metadata)
       integer, intent(in)             :: yac_instance_id
       character(len=*), intent(in)    :: grid_name
       character(len=*), intent(in)    :: metadata
     end subroutine yac_fdef_grid_metadata_instance
  end interface yac_fdef_grid_metadata

  interface yac_fgrid_has_metadata
    function yac_fgrid_has_metadata(grid_name) result( has_metadata )
      character(len=*), intent(in) :: grid_name
      logical :: has_metadata
    end function yac_fgrid_has_metadata

    function yac_fgrid_has_metadata_instance(yac_instance_id, grid_name) &
      result( has_metadata )
      integer, intent(in) :: yac_instance_id
      character(len=*), intent(in) :: grid_name
      logical :: has_metadata
    end function yac_fgrid_has_metadata_instance
  end interface yac_fgrid_has_metadata

  interface yac_fget_grid_metadata
     function yac_fget_grid_metadata(grid_name) result( metadata )
       character(len=*), intent(in) :: grid_name
       character(len=:), allocatable :: metadata
     end function yac_fget_grid_metadata

     function yac_fget_grid_metadata_instance(yac_instance_id, grid_name) &
          result( metadata )
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in) :: grid_name
       character(len=:), allocatable :: metadata
     end function yac_fget_grid_metadata_instance
  end interface yac_fget_grid_metadata

  interface yac_fdef_field_metadata
     subroutine yac_fdef_field_metadata(comp_name, grid_name, field_name, metadata)
       character(len=*), intent(in)    :: comp_name
       character(len=*), intent(in)    :: grid_name
       character(len=*), intent(in)    :: field_name
       character(len=*), intent(in)    :: metadata
     end subroutine yac_fdef_field_metadata

     subroutine yac_fdef_field_metadata_instance(yac_instance_id, comp_name, &
          grid_name, field_name, metadata)
       integer, intent(in)             :: yac_instance_id
       character(len=*), intent(in)    :: comp_name
       character(len=*), intent(in)    :: grid_name
       character(len=*), intent(in)    :: field_name
       character(len=*), intent(in)    :: metadata
     end subroutine yac_fdef_field_metadata_instance
  end interface yac_fdef_field_metadata

  interface yac_ffield_has_metadata
    function yac_ffield_has_metadata(comp_name, grid_name, field_name) &
        result( has_metadata )
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       logical :: has_metadata
    end function yac_ffield_has_metadata

    function yac_ffield_has_metadata_instance( &
        yac_instance_id, comp_name, grid_name, field_name) &
      result( has_metadata )
      integer, intent(in) :: yac_instance_id
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       logical :: has_metadata
    end function yac_ffield_has_metadata_instance
  end interface yac_ffield_has_metadata

  interface yac_fget_field_metadata
     function yac_fget_field_metadata(comp_name, grid_name, field_name) &
          result( metadata )
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       character(len=:), allocatable :: metadata
     end function yac_fget_field_metadata

     function yac_fget_field_metadata_instance(yac_instance_id, comp_name, &
          grid_name, field_name) &
          result( metadata )
       integer, intent(in) :: yac_instance_id
       character(len=*), intent(in) :: comp_name
       character(len=*), intent(in) :: grid_name
       character(len=*), intent(in) :: field_name
       character(len=:), allocatable :: metadata
     end function yac_fget_field_metadata_instance
  end interface yac_fget_field_metadata

  ! ---------------------------------------------------------------------

  interface

    subroutine yac_fget_action ( field_id, action )

      integer, intent (in)  :: field_id !< [IN]  field identifier
      integer, intent (out) :: action   !< [OUT] action for the current timestep\n
                                        !!       (\ref YAC_ACTION_NONE,
                                        !!        \ref YAC_ACTION_COUPLING,
                                        !!        \ref YAC_ACTION_GET_FOR_RESTART,
                                        !!        \ref YAC_ACTION_PUT_FOR_RESTART, or
                                        !!        \ref YAC_ACTION_OUT_OF_BOUND)

    end subroutine yac_fget_action

  end interface

  ! ---------------------------------------------------------------------

  interface

    subroutine yac_fupdate ( field_id )

      integer, intent (in)  :: field_id !< [IN]  field identifier

    end subroutine yac_fupdate

 end interface

! ---------------------------------------------------------------------

  interface yac_redirstdout

     subroutine yac_redirstdout_f ( filestem, parallel, my_pe, npes )

       character(len=*), intent (in) :: filestem
       integer, intent(in) :: parallel
       integer, intent(in) :: my_pe
       integer, intent(in) :: npes

     end subroutine yac_redirstdout_f

  end interface yac_redirstdout

  ! --- ISO_C interface -------------------------------------------------

  !> \example test_abort.F90
  !! This contains an example on how to use yac_abort_message.

  interface

     subroutine yac_abort_message ( text, file, line ) &
           bind ( c, name='yac_abort_message' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: text       !< [IN] error message
       character ( kind=c_char ), dimension(*) :: file       !< [IN] error message
       integer ( kind=c_int ), value           :: line       !< [IN] line number

     end subroutine yac_abort_message

  end interface

end module mo_yac_finterface
